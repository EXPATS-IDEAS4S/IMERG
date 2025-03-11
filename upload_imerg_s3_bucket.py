import os
import boto3
from glob import glob
import logging
from botocore.exceptions import ClientError
import xarray	

from s3_credential import S3_BUCKET_NAME, S3_ACCESS_KEY, S3_SECRET_ACCESS_KEY, S3_ENDPOINT_URL


def upload_file(s3_client, file_name, bucket, object_name=None):
    """Upload a file to an S3 bucket"""

    if object_name is None:
        object_name = os.path.basename(file_name)

    try:
        with open(file_name, "rb") as f:
            #file_size = os.path.getsize(file_name)  # Get file size
            #print(f"Uploading {file_name} -> S3://{bucket}/{object_name} ({file_size} bytes)")

            # Upload file with explicit ContentLength
            s3_client.upload_fileobj(f, bucket, object_name)
            
    except ClientError as e:
        logging.error(f"S3 Upload Error: {e}")
        return False
    return True



# Initialize the S3 client
s3 = boto3.client(
    's3',
    endpoint_url=S3_ENDPOINT_URL,
    aws_access_key_id=S3_ACCESS_KEY,
    aws_secret_access_key=S3_SECRET_ACCESS_KEY
)

# # List the objects in our bucket
# response = s3.list_objects(Bucket=S3_BUCKET_NAME)
# for item in response['Contents']:
#     print(item['Key'])

# exit()

#Directory with the data to uplad
years = [2013, 2014, 2015, 2016]
months = range(4,10)
path_dir = f"/work/dcorradi/IMERG/data"

for year in years:
    for month in months:
        data_filepattern = f"{path_dir}/{year:04d}/{month:02d}/*.nc"
        file_list = sorted(glob(data_filepattern))

        for file in file_list:
            if not os.path.exists(file):
                print(f"Error: File {file} does not exist.")
                continue  # Skip this file

            file_size = os.path.getsize(file)
            if file_size == 0:
                print(f"Error: File {file} is empty (0 bytes).")
                continue  # Skip empty files

            print(f"Uploading {file} ({file_size} bytes) to S3...")
            upload_file(s3, file, S3_BUCKET_NAME, os.path.basename(file))

         


# List the objects in our bucket
response = s3.list_objects(Bucket=S3_BUCKET_NAME)
for item in response['Contents']:
    print(item['Key'])