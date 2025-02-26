import os
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from glob import glob
import numpy as np

years = [2014,2015,2016]

date_low_plot = "2014-07-14T00:00:00Z"
date_high_plot = "2014-07-15T00:00:00Z"

var = "precipitation"

# Define domain
latmin = 42
latmax = 51.5
lonmin = 5
lonmax = 16

for year in years:
    print(f"📅 Year: {year}")

    # Define the directory containing the IMERG NetCDF files
    data_dir = f"/data/sat/msg/IMERG/{year}"

    output_dir = f"/work/dcorradi/IMERG/fig/{year}"
    os.makedirs(output_dir, exist_ok=True)

    # Get a list of all NetCDF files in the directory
    nc_files = sorted(glob(os.path.join(data_dir, "*.nc4")))
    print(f"📂 Found {len(nc_files)} NetCDF files in {data_dir}")

    # Loop through each file
    for file in nc_files:
        file_path = os.path.join(data_dir, file)
        
        # Open NetCDF file using xarray
        ds = xr.open_dataset(file_path)
        
        # Print file metadata and variables
        print(f"📂 File: {file}")

        if var in ds.variables:
            ds_prec = ds[var]
        else:
            raise ValueError(f"Variable {var} not found in {file}")
        #print(ds_prec)        

        # Select the first time step if data has a time dimension
        if "time" in ds_prec.dims:
            precip = ds_prec.isel(time=0)
            time = ds.time.values[0]
            print(f"🕒 Time: {time}")
            # print(precip.lat.values.shape)
            # print(precip.lon.values.shape)
            # print(precip.values.shape)
        
        if time > np.datetime64(date_low_plot) and time < np.datetime64(date_high_plot):

            # Plot the data
            fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})
            ax.set_global()

            lon_grid, lat_grid = np.meshgrid(precip.lon.values, precip.lat.values, indexing="ij")
            # print(lat_grid.shape)
            # print(lon_grid.shape)
            
            # Plot precipitation
            im = ax.pcolormesh(lon_grid, lat_grid, precip.values, transform=ccrs.PlateCarree(), cmap="tab20c_r", vmin=0, vmax=10)
            
            # Add features
            ax.coastlines()
            ax.add_feature(cfeature.BORDERS, linestyle=":")
            ax.add_feature(cfeature.LAND, facecolor="lightgray")

            # Set extent
            ax.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())

            # Colorbar
            cbar = plt.colorbar(im, ax=ax, orientation="horizontal", pad=0.05, shrink=0.7, extend="max")
            cbar.set_label("Precipitation (mm/hr)")

            # Title
            plt.title(f"{str(ds.time.values[0]).split('.')[0]}")
            
            # Show plot
            plt.savefig(os.path.join(output_dir, f"IMERG_merged-microwave-infrared-gauge-adjusted_{str(ds.time.values[0]).split('.')[0]}.png"), bbox_inches="tight")
            plt.close()
            print(f"📊 Saved plot to {output_dir}")

        # Close dataset
        ds.close()
        
