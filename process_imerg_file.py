import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from glob import glob
import pandas as pd
import xesmf as xe

# Years to process
years = [2013,2014, 2015, 2016]

# Define date range for plotting
date_low_plot = np.datetime64("2013-04-02T00:00:00Z")
date_high_plot = np.datetime64("2013-04-03T00:00:00Z")

var = "precipitation"

# Define domain
latmin, latmax = 42, 51.5
lonmin, lonmax = 5, 16
res_grid = 0.04

# Define target grid for regridding
target_grid = {
    "lon": np.arange(lonmin, lonmax + res_grid, res_grid),
    "lat": np.arange(latmin, latmax + res_grid, res_grid),
}

# Process each year
for year in years:
    print(f"📅 Processing Year: {year}")

    # Define input/output directories
    data_dir = f"/data/sat/msg/IMERG/{year}"
    fig_output_dir = f"/work/dcorradi/IMERG/fig/{year}"
    nc_output_dir = f"/work/dcorradi/IMERG/data/{year}"
    os.makedirs(fig_output_dir, exist_ok=True)
    os.makedirs(nc_output_dir, exist_ok=True)

    # Get NetCDF files
    nc_files = sorted(glob(os.path.join(data_dir, "*.nc4")))
    print(f"📂 Found {len(nc_files)} NetCDF files in {data_dir}")

    daily_datasets = {}  # Dictionary to store datasets by date

    # Loop through each file
    for file in nc_files:
        file_path = os.path.join(data_dir, file)

        # Open dataset
        ds = xr.open_dataset(file_path)
        print(f"📂 Processing File: {file}")

        # Check if precipitation variable exists
        if var in ds.variables:
            ds_prec = ds[var]
        else:
            print(f"⚠️ Variable {var} not found in {file}")
            ds.close()
            continue

        # Select first time step
        if "time" in ds_prec.dims:
            precip = ds_prec.isel(time=0)
            #convert time to pandas datetime
            time_np = np.datetime64(ds.time.values[0])
            time_pd = pd.to_datetime(str(time_np))
            print(f"🕒 Time: {time_pd}")
            exit()

        # Extract date (YYYY-MM-DD) for daily merging
        month = time_pd.month
        day = time_pd.day
        date_str = f"{year}-{month:02d}-{day:02d}"

        # Create a regridding target dataset
        ds_target = xr.Dataset(
            {"lat": (["lat"], target_grid["lat"]), "lon": (["lon"], target_grid["lon"])}
        )

        # Create regridder
        regridder = xe.Regridder(ds, ds_target, "linear", periodic=True)

        # Apply regridding
        precip_regridded = regridder(precip)

        # Save the dataset for merging
        if date_str not in daily_datasets:
            daily_datasets[date_str] = []
        daily_datasets[date_str].append(precip_regridded)

        if time_pd > date_low_plot and time_pd < date_high_plot:
        
            # Plot side-by-side comparison
            fig, axes = plt.subplots(ncols=2, figsize=(12, 5), subplot_kw={"projection": ccrs.PlateCarree()})

            # Define plotting function
            def plot_precip(ax, data, title, cmap="tab20c_r"):
                lon_grid, lat_grid = np.meshgrid(data.lon.values, data.lat.values, indexing="ij")
                im = ax.pcolormesh(lon_grid, lat_grid, data.values, transform=ccrs.PlateCarree(), cmap=cmap, vmin=0, vmax=10)
                ax.coastlines()
                ax.add_feature(cfeature.BORDERS, linestyle=":")
                ax.add_feature(cfeature.LAND, facecolor="lightgray")
                ax.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())
                cbar = plt.colorbar(im, ax=ax, orientation="horizontal", pad=0.05, shrink=0.7, extend="max")
                cbar.set_label("Precipitation (mm/hr)")
                ax.set_title(title)

            # Plot original data
            plot_precip(axes[0], precip, "Original Grid")

            # Plot regridded data
            plot_precip(axes[1], precip_regridded, f"Regridded ({res_grid}°)")

            # Save figure
            output_file = os.path.join(fig_output_dir, f"IMERG_{date_str}.png")
            plt.savefig(output_file, bbox_inches="tight")
            plt.close()
            print(f"📊 Saved plot: {output_file}")

        # Close dataset
        ds.close()

    # Merge and save daily datasets
    for date_str, datasets in daily_datasets.items():
        merged_ds = xr.concat(datasets, dim="time")  # Merge all time steps for the same day
        merged_ds = merged_ds.mean(dim="time")  # Compute daily mean

        # Define output directory and filename
        year_month = date_str[:7]  # YYYY-MM format
        month_output_dir = os.path.join(nc_output_dir, year_month)
        os.makedirs(month_output_dir, exist_ok=True)
        output_nc_file = os.path.join(month_output_dir, f"IMERG_daily_{date_str}.nc")

        # Save with compression
        comp = dict(zlib=True, complevel=5)  # Compression level 5 (balance between speed and size)
        encoding = {var: comp}

        merged_ds.to_netcdf(output_nc_file, encoding=encoding)
        print(f"💾 Saved daily dataset: {output_nc_file}")
