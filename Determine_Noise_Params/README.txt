This folder is used to determine the process noise covariances for the position and velocity states. 
Running the simKnownFexInput.m file is the place to start. This function will generate a structure called
runStruct.mat that you can place in the data folder if you want to update this information, this structure contains
information on the noise in each run, wave parameters as well as scaled factors. This structure also contains information
on the steepness of the wave and the ratio of characteristic buoy dimension and wavelength which can be used for 
gaining intuition on the wave situation.

