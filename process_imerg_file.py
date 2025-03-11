import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from glob import glob
import pandas as pd
from scipy.interpolate import griddata

def regrid_with_scipy(ds, target_grid, var_name):
    """
    Perform regridding using SciPy's griddata.
    """
    # Extract original lat/lon and precipitation data
    orig_lats = ds.lat.values.squeeze()
    orig_lons = ds.lon.values.squeeze()
    orig_precip = ds.values.squeeze()
    time = ds.time.values[0]  # Preserve time coordinate
    
    # Create meshgrid for original data
    lon_orig_mesh, lat_orig_mesh = np.meshgrid(orig_lons, orig_lats, indexing="ij")
    points = np.array([lon_orig_mesh.ravel(), lat_orig_mesh.ravel()]).T
    values = orig_precip.ravel()
    
    # Create target grid mesh
    lon_target_mesh, lat_target_mesh = np.meshgrid(target_grid["lon"], target_grid["lat"], indexing="ij")
    target_points = np.array([lon_target_mesh.ravel(), lat_target_mesh.ravel()]).T
    
    # Perform interpolation
    precip_regridded = griddata(points, values, target_points, method="nearest")
    precip_regridded = precip_regridded.reshape(lon_target_mesh.shape)
    
    # Convert to xarray DataArray (without time initially)
    da_regridded = xr.DataArray(
        precip_regridded.astype('float32'),
        dims=("lon", "lat"),
        coords={
            "lon": target_grid["lon"].astype('float32'),
            "lat": target_grid["lat"].astype('float32')
           
        },
        name=var_name
    )
    
    # Expand to include time dimension
    da_regridded = da_regridded.expand_dims(time=[time])

    return da_regridded

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

# Years to process
years = [2013, 2014, 2015, 2016]

# Define date range for plotting
plot = False
date_low_plot = np.datetime64("2013-04-01T00:00:00Z")
date_high_plot = np.datetime64("2013-04-02T00:00:00Z")

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
    nc_output_dir = f"/work/dcorradi/IMERG/data/"
    os.makedirs(fig_output_dir, exist_ok=True)
    os.makedirs(nc_output_dir, exist_ok=True)

    # Get NetCDF files
    nc_files = sorted(glob(os.path.join(data_dir, "*.nc4")))
    print(f"📂 Found {len(nc_files)} NetCDF files in {data_dir}")

    current_date = None
    daily_datasets = []

    for file in nc_files:
        file_path = os.path.join(data_dir, file)
        ds = xr.open_dataset(file_path)
        #print(ds)
        
        print(f"📂 Processing File: {file}")

        if var in ds.variables:
            ds_prec = ds[var]
            #print(ds_prec)
        else:
            print(f"⚠️ Variable {var} not found in {file}")
            ds.close()
            continue

        if "time" in ds_prec.dims:
            #precip = ds_prec.isel(time=0)
            time_np = np.datetime64(ds.time.values[0])
            time_pd = pd.to_datetime(str(time_np))
            print(f"🕒 Time: {time_pd}")

        date_str = time_pd.strftime("%Y-%m-%d")

        if current_date is None:
            current_date = date_str

        # If the day changes, merge and save the previous day's dataset
        if date_str != current_date and daily_datasets:
            merged_ds = xr.concat(daily_datasets, dim="time")
            merged_ds.attrs = daily_datasets[0].attrs
            
            month_output_dir = os.path.join(nc_output_dir, current_date[:4], current_date[5:7])
            os.makedirs(month_output_dir, exist_ok=True)
            output_nc_file = os.path.join(month_output_dir, f"IMERG_daily_{current_date}.nc")

            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp}
            
            merged_ds.to_netcdf(output_nc_file, encoding=encoding)
            print(f"💾 Saved daily dataset: {output_nc_file}")

            # Reset for new day
            daily_datasets = []
            current_date = date_str

        # Process current file
        precip_regridded = regrid_with_scipy(ds_prec, target_grid, var)
        precip_regridded = precip_regridded.where(precip_regridded != -9999.9)
        #print(precip_regridded.time)
       
        daily_datasets.append(precip_regridded)

        if plot and (date_low_plot <= time_pd <= date_high_plot):
            fig, axes = plt.subplots(ncols=2, figsize=(10, 5), subplot_kw={"projection": ccrs.PlateCarree()})
            plot_precip(axes[0], ds_prec, "Original Grid (0.1°)")
            plot_precip(axes[1], precip_regridded, f"Regridded ({res_grid}°)")
            plt.suptitle(f"{time_pd.strftime('%Y-%m-%d %H:%M')} UTC")
            plt.subplots_adjust(hspace=0.5)
            output_file = os.path.join(fig_output_dir, f"IMERG_{time_pd}.png")
            plt.savefig(output_file, bbox_inches="tight")
            plt.close()
            print(f"📊 Saved plot: {output_file}")

        ds.close()

    # Save the last day's dataset
    if daily_datasets:
        merged_ds = xr.concat(daily_datasets, dim="time")
        merged_ds.attrs = daily_datasets[0].attrs
        
        month_output_dir = os.path.join(nc_output_dir, current_date[:4], current_date[5:7])
        os.makedirs(month_output_dir, exist_ok=True)
        output_nc_file = os.path.join(month_output_dir, f"IMERG_daily_{current_date}.nc")

        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp}
        
        merged_ds.to_netcdf(output_nc_file, encoding=encoding)
        print(f"💾 Saved daily dataset: {output_nc_file}")

#2896554 nohup
