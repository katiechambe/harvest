# Instructions for environment setup

In order to make sure that the correct modules are installed and loaded into Python for the analysis tools, you should create a new conda environment so that the necessary dependencies here do not affect other ongoing projects you may have. 

- First, make sure you've installed conda. 
Instructions can be found here: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
* Note: there are two options, Anaconda and Miniconda. Anaconda uses >3GB of storage space due to the additional GUI and packages that make this install a little bit more "user friendly." Miniconda uses <500MB of disk space, and primarily only downloads code that is necesssary to make the conda ecosystem run. The remainder of the setup instructions here can be followed with either install, however I would recommend miniconda since it's less bulky! 

- Second, create a new environment for projects that use the code in this repository. 
I have included a environ.yml file that specifies all of the dependencies necessary. You can create a new environment with all these dependencies installed by running 
    ```conda env create -f environ.yml```
The name of the environment is specified in the first line of the yml file and can be changed, but default name is harvest.
* More information about managing environments can be found on the conda website at https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html 

- Third, before opening or running any scripts from this project, you MUST activate the harvest environment by:
    ```conda activate harvest```
This step must be done EACH TIME before using the code herein. 