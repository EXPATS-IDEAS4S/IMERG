import os
import re
from datetime import datetime, timedelta

def generate_expected_timestamps(year):
    """
    Generate all expected timestamps from April 1st to September 30th of the given year,
    every 30 minutes.
    """
    start_date = datetime(year, 4, 1, 0, 0)
    end_date = datetime(year, 9, 30, 23, 30)
    delta = timedelta(minutes=30)
    
    timestamps = []
    current_time = start_date
    while current_time <= end_date:
        timestamps.append(current_time.strftime('%Y%m%d-S%H%M%S'))
        current_time += delta
    
    return set(timestamps)

def extract_timestamps_from_filenames(files):
    """
    Extract timestamps from filenames using regex.
    """
    timestamp_pattern = re.compile(r'3B-HHR\.MS\.MRG\.3IMERG\.(\d{8}-S\d{6})')
    found_timestamps = set()
    
    for file in files:
        match = timestamp_pattern.search(file)
        if match:
            found_timestamps.add(match.group(1))
    
    return found_timestamps

def check_missing_files(base_path):
    """
    Check for missing timestamps in the /data/sat/msg/IMERG/YYYY directories.
    """
    years = [d for d in os.listdir(base_path) if d.isdigit()]
    print(years)
    missing_files_report = "missing_files.txt"
    
    with open(missing_files_report, 'w') as report:
        for year in years:
            year_path = os.path.join(base_path, year)
            if not os.path.isdir(year_path):
                continue
            
            expected_timestamps = generate_expected_timestamps(int(year))
            found_timestamps = extract_timestamps_from_filenames(os.listdir(year_path))
            
            missing_timestamps = expected_timestamps - found_timestamps
            if missing_timestamps:
                report.write(f"Missing timestamps for {year}:")
                for timestamp in sorted(missing_timestamps):
                    report.write(f"{timestamp}\n")
                report.write("\n")

if __name__ == "__main__":
    base_directory = "/data/sat/msg/IMERG"
    check_missing_files(base_directory)
    print("Check completed. Missing timestamps (if any) are recorded in missing_files.txt.")
