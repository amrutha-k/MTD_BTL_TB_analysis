# Lab5015Analysis
collection of programs for the analysis of Lab5015 measurements



### Login to your favourite machine
```sh
ssh username@hostname
```



### Fresh installation of the analysis package
```sh
export MYNAME=putYourNameHere  #use your name for development
mkdir $MYNAME
cd $MYNAME
git clone --recursive https://github.com/Lab5015/Lab5015Analysis
cd Lab5015Analysis
source scripts/setup.sh
make
make exe
```


### Package structure
This package is structured as follows
- Utilities are organized in the `interface` and `src` folders:
    - `interface`: contains the class or functions headers
    - `src`: contains the class or functions implementation
- `main`: contains the main analysis code that make use of the Utilities
- `cfg`: contains the config files which are used to pass parameters to the executables

After the compilation, each step of the analysis can be executed from the main folder `Lab5015Analysis` with a command like:
`./bin/executable cfg/configfile`



### Before running the analysis
The output filename and output location of each analysis step are defined in the cfg files. Before running the analysis make sure you have checked and if needed updated the relevant output paths in order not to overwrite the work of others.

Other than that, every time you login remember to source the setup script:
```
cd Lab5015Analysis
source scripts/setup.sh
```


### Run the analysis
The analysis of the collected data is structured in three steps
1. `moduleCharacterization_step1.cpp`:
   This step performs the first loop over the events and fills base histograms such as energy and ToT for each bar, each theshold and each over-voltage. A skim of the events according to the selections defined in the config file is also performed. The output is in the form of TTrees and histograms.

1. `moduleCharacterization_step2.cpp`:
   This step performs several loops over the events:
    1. loop over the histos filled in step1 and define the energy ranges
    1. loop over the skimmed TTrees and fill raw distribution of energy and energy ratio according to the predefined bins
    1. loop over the skimmed TTrees and compute time walk corrections
    1. loop over the skimmed TTrees and apply the time walk corrections
    
   The output of this step are histograms.

1. `moduleCharacterizationSummaryPlots.py`:
   This step takes the output of step2 as an input and displays summary plots in a website. Loop over the events doesn't belong here.
   example:
   ```sh
   python moduleCharacterizationSummaryPlots.py -m 2 -i run2071 -o /var/www/html/TOFHIR2X/ModuleCharacterization/run2071
   ```


An additional code which is useful to plot the pulse shape for a given channel is `drawPulseShape.exe`. The pulse shape is computed with respect to the trigger.
example:
```sh
./bin/drawPulseShape.exe cfg/drawPulseShape_HPK_2E14_52deg_T-40C_Vov1.50.cfg`
```


### Submit the analysis in parallel on pcfatis
A script under the `scripts` folder allows the submission of multiple jobs (e.g. a set of runs corresponding to an overvoltage/threshold scan) in parallel on `pcfatis`. one job per run. The script is used as follows:
```
python scripts/submit_moduleCharacterization.py --label myTask_stepX -b /path/of/Lab5015Analysis -e bin/moduleCharacterization_stepX.exe -r run1-run2,run3 -c cfg/moduleCharacterization.cfg --submit
```
The config file is used as a template for each job, and the run number and file names in it are modified by the script as required. The command can be executed without the `--submit` flag first to inspect the jobs (in the `scripts/jobs` folder) before actually running them.



### Visualize the results
All the results are available for inspection on a website hosted on pcfatis. The [link](http://pcfatis.mib.infn.it) is accessible from the INFN network or via tunnel.



