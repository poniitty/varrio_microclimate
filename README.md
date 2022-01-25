# varrio_microclimate
This repository contains microclimate data collected from Värriö area, Finland and the code to combine the data and correct for some problems in the data


The ./data folder contain all raw logger data from TOMST TMS4 and LogTag HAXO-8 loggers.

### The R scripts in ./scripts folder execute the following tasks:

**initial_file_sorting.R** resolves study site names and organises the datafiles  
**read_haxo_varrio.R** Reads in and combines HAXO logger data and masks periods when logger under snow  
**read_tomst_varrio.R** Reads in and combines TOMST logger data. Screens the data for suspicious peaks   
**correct_haxo_varrio.R** Corrects HAXO data for nonmatching time stamps  
**microclimate_derivatives.R** Calculates daily and monthly statistics from the data  
