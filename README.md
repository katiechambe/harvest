Tools to collect and analyze TNG galaxy/subhalo pairs. 

# Instructions for download and set-up

Step 1:
-------
Clone this repository on your local machine:
```git clone https://github.com/katiechambe/harvest.git```
or, if you have ssh keys set up on your local machine:
```git clone git@github.com:katiechambe/harvest.git```
This will download the necessary Python code and a small subset of data 
used in the example Jupyter Notebooks 1 & 2.

Step 2:
-------
Add the utility scripts to your python path. In your .bashrc file (or .zshrc, if using zshell),
append the following:
``` export PYTHONPATH=<path-to-harvest-utils>/harvest:$PYTHONPATH```
For example, this line reads "export PYTHONPATH=:/xdisk/gbesla/katiechambe/harvest/harvest:$PYTHONPATH" 
in my personal bashrc file. *Make sure to add the path to the harvest/harvest directory where the harvesting-tools directory lives. 


You can also append this path to your .bashrc file with 
```echo "export PYTHONPATH=<path-to-harvest-utils>/harvest:$PYTHONPATH" >> ~/.bashrc```

Step 3:
-------
Before you run any of the code, you MUST set up a conda environment that 
downloads all necessary script dependencies prior to running. 
Dependencies are specified in SETUP/environ.yml, which is a file that can be used 
to create a conda environment that mimics my own. 

More details about setting up your conda environment can be found in SETUP.md. 
After completing the steps therein, make sure you activate the new conda environment:
```conda activate harvest```

Step 4:
-------
Start a Jupyter Notebook session in the harvest/ directory and any of it's parent directories, and get started with the first example notebook in ```example-notebooks/``` 

# Directories
```
------------
- harvest
    |
    |- data : Location of data to create subsamples, or for analysis
    |   | 
    |   |- pairs
    |   |- orbits
    |
    |- harvest
    |- notebooks
    |- SETUP
    |- _dev
```
