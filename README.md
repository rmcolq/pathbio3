# pathbio3

To run locally, first clone the repository using 
```
git clone https://github.com/rmcolq/pathbio3.git
```

This may take up to 5 minutes because the subsampled read files are still quite large.

Next install a conda environment to run in. It is faster to use mamba to do this, but conda should also work if you don't have mamba installed.
```
cd pathbio3
mamba env create -f environment.yml
conda activate pathbio3
```

Finally open up the notebooks in your browser with jupyter using 
```
jupyter lab
```

Note: this is a work in progress and is likely to change a lot...
