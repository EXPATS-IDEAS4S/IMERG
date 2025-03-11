import os
import xarray as xr
from glob import glob

def check_merged_nc(years):
    base_dir = "/work/dcorradi/IMERG/data"
    
    for year in years:
        year_dir = os.path.join(base_dir, str(year))
        if not os.path.exists(year_dir):
            print(f"⚠️ Year directory not found: {year_dir}")
            continue
        
        months = sorted(next(os.walk(year_dir))[1])  # List all month directories
        for month in months:
            month_dir = os.path.join(year_dir, month)
            nc_files = sorted(glob(os.path.join(month_dir, "IMERG_daily_*.nc")))
            
            if not nc_files:
                print(f"⚠️ No NetCDF files found in {month_dir}")
                continue
            
            for nc_file in nc_files:
                try:
                    ds = xr.open_dataset(nc_file)
                    print(f"✅ Successfully opened: {nc_file}")
                    print(ds)
                    #print(ds.precipitation.values)
                    ds.close()
                except Exception as e:
                    print(f"❌ Error opening {nc_file}: {e}")

# Example usage
check_merged_nc([2013])