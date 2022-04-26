# MTD BTL TB analysis
Based on Lab5015/Lab5015Analysis. Codes adapted for spatial uniformity studies of sensor modules in MTD BTL TB at CERN (Oct2021)



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
```
./bin/executable cfg/configfile
eg: ./bin/moduleCharacterization_step1.cpp cfg/moduleCharacterization.cfg
```


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

Note: moduleCharacterization_step1_new.cpp and moduleCharacterization_step2_new.cpp are 1 and 2 adapted to include MCP information. Use config file cfg/moduleCharacterization_HPK_MCP.cfg

### Results and Documentation

1. [BTL TB analysis working meeting, 3 December 2021](https://indico.cern.ch/event/1102712/contributions/4638951/attachments/2358593/4025820/BTL%20OCT%20TB%20analysis%20meeting.pdf).
2. [BTL TB analysis working meeting, 18 February 2022](https://indico.cern.ch/event/1126254/contributions/4727432/attachments/2393938/4092895/BTL%20OCT%20TB%20analysis%20meeting%20-%20status%20update.pdf).
3. [BTL TB analysis working meeting, 18 March 2022](https://indico.cern.ch/event/1140648/contributions/4786478/attachments/2410506/4124651/BTL_OCT_TB%20analysis%20meeting_status%20update%20_18March.pdf).
 



