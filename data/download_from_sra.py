#!/usr/bin/env python 
import subprocess
import os

def fetch_data(sra_id):
    # Execute the command (prefetching) and capture its output
    prefetch_cmd = ["prefetch", sra_id]
    print(f"\nRunning: {' '.join(prefetch_cmd)}")
    prefetch = subprocess.Popen(prefetch_cmd, stdout=subprocess.PIPE)
    output, _ = prefetch.communicate()

    # Convert the output bytes to a string
    output = output.decode('utf-8')
    print(output)

    # Split the output into words
    words = output.split()

    # Find the word that starts with "SRR"
    srr_id = None
    for word in words:
        if word.startswith("'SRR"):
            srr_id = word.split("'")[1]
            break

    # fetch for data
    if srr_id:
        # Now download the fastq and capture its output
        fetch_cmd = ["fastq-dump", "--outdir", "raw", "--gzip", "--skip-technical",  "--readids", "--read-filter", "pass", "--dumpbase", "--split-3", "--clip", f"{srr_id}/{srr_id}.sra"]
        print(f"\nRunning: {' '.join(fetch_cmd)}")
        fetch = subprocess.Popen(fetch_cmd, stdout=subprocess.PIPE)
        output, _ = fetch.communicate()

        # Convert the output bytes to a string
        output = output.decode('utf-8')
        print(output)
    else:
        print("No word starting with 'SRR' found.")


def subsample(sra_id):
    subsample_cmd = ["rasusa", "-n", "20m", "-i", f"raw/{sra_id}_pass_1.fastq.gz", f"raw/{sra_id}_pass_2.fastq.gz", "-o", f"subsampled/{sra_id}_1.fastq.gz", f"subsampled/{sra_id}_2.fastq.gz"]
    print(f"\nRunning: {' '.join(subsample_cmd)}")
    subsample = subprocess.Popen(subsample_cmd, stdout=subprocess.PIPE)
    output, _ = subsample.communicate()

    # Convert the output bytes to a string
    output = output.decode('utf-8')
    print(output)

    os.remove(f"raw/{sra_id}_pass_1.fastq.gz")
    os.remove(f"raw/{sra_id}_pass_2.fastq.gz")


downsample=False
with open("list_ids.txt", 'r') as f:
    for line in f:
        sra_id = line.strip()
        if sra_id == "subsampled":
            downsample = True
            continue
        print(f"Fetching {sra_id}...")
        fetch_data(sra_id)
        os.remove(f"{sra_id}/{sra_id}.sra")
        os.rmdir(sra_id)
        if downsample:
            subsample(sra_id)  

