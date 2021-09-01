# TremRoutine 1.0

Repository of automatic processing routine for volcanic tremor signals based on polarization analysis.

The code here was made for the master thesis named "Detection and Location of Tremor signals: A case study from East Java, Indonesia".

## Table of Contents
* [General Info](#general-info)
* [Installation](#installation)
* [References](#references)

## General Info

![Architecture](./code_architecture.png)

The root directory or folder DO NOT make reference to the root folder of any UNIX-core type device. Therefore, "root directory" means a main folder of the user's preference on where to store the cloned repository (this one).

Due to confidentiality of ERC-funded Lusi Lab project, the initial dataset cannot be shared (waveform files and response files) as well as the results. Therefore, the "data" and "results" folders shown in the architecture scheme of the program are not included. However, any used dataset must be organized in the specified hierarchical order of the folders.

The content inside "results" will be automatically generated. Files generated for statistics under "stats" folder and the database of the events in the "events" folder will be as .JSON files. Each one in "events" will be a list of dictionaries, where each dictionary corresponds to one event:

   - Stations and components involved in the detection of that specific events.
   - Waveform file paths involved in the detection.
   - Arrival times per station (tart-end of the detection).
   - Start-End times for waveform extraction of the event in each station involved (for post-processing cases).

The folder "pyFunctions" contains most of the Python files that has the methods for manipulation files around the architecture, between the third-party codes and the other folders, detection and association steps, polarization analysis, plotting functions, and statistics of the detections achieved.

All images that can be produced will be stored in the folder called "figures" if running the plotting functions inside "pyFunctions/plotting.py" and deciding to save them. Currently no folder named "figures" is presented by the same confidentiality agreement, as the files also present results. A folder called with the same name can just be created from scratch by the user.

The "metadata" folder contains files related to the position of stations and volcanoes. Also other priory information can be saved here.

The "programs" folder is designed to contain third-party codes and external programs such as NonLinLoc (Lomax, 2004) that can be used for the location part. Normally this code is designed to work with NonLinLoc installed/compiled within a folder called "NonLinLoc" inside of the "programs" folder. Internal executions for travel-time calculations and probailistic model must be made before with the commands of NonLinLoc locally on the folder of the code.

The other Python files outside of the folders are made for managing functions that are mainly in "pyFunctions" folder. This includes:

  - _"locate.py"_: code that manages the detection and association of detected signals of possible events. It then generates the database of the events in the "results" and "stats" folder. Control parameters are defined inside the file. It has the option for parallelizing the process.
  - _"make_prestats.py"_: code that generates additional information as the polarization analysis, P and Rayleigh wave classification, Maximum Amplitude, and others. Once runned, it will update the files inside the database of the events with the new calculations. It has the option for parallelizing the process.
  - _"NNLarange.py"_: code for managing files during the Location inversion with NonLinLoc, such as copying or deleting files, updating the database with the locations, errors, etc. Location inversions will be done inside "NonLinLoc", under "programs" folder.
  - _"playground.py"_: not a relevat file for processing, but just a space for playing with the code, calling other functions such as plotting, printing statistics, exploring detection methods and benchmark tests.
  - _"clone_superpc.sh"_: shell script for updating files to the Sigma2 HPC servers.
  - _"loc_inv.sh"_: shell script for that once run, it automatically starts the location inversion for all the database by batches. This involves the creation of .obs files for NonLinLoc, the copy of those files into NonLinLoc respective folder, running the inversion, copying the output and writing it to updated files of the database, and cleaning of the NonLinLoc .obs files.
  - _"routine.sh"_: shell script for running the entire routine, from start (detection and association) to end (location of events in database). Not recommended to run as it will be wise to run the codes separately so the user can supervise each process and have quality control in all of them.


## Installation

## References

  - Lomax, A. (2004). Probabilistic, non-linear, global-search earthquake location in 3d media.AnthonyLomax Scientific Software, Mouans–Sartoux, France.
