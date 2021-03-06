# stimulation-motion

Latency analysis code used in

Cerebral organoids at the air-liquid interface generate diverse nerve tracts with functional output
Stefano L Giandomenico, Susanna B Mierau, George M Gibbons, Lea MD Wenger, Laura Masullo, Timothy Sit, Magdalena Sutcliffe, Jerome Boulanger, Marco Tripodi, Emmanuel Derivery, Ole Paulsen, Andras Lakatos, Madeline Lancaster
bioRxiv 353151; doi: https://doi.org/10.1101/353151 

The sample is stimulated using electrical pulses at a given frequency and the twitching is monitored using a microscope. The stimulation and the fire TTL signals are recorded using the picoscope https://www.picotech.com/  in parallel of the image acquisition. Motion can then be quantified using the Quantify_Motion.ijm imagej macro and latency analyzed using the latencyanalysis.m script.

The imageJ macro collect the timestamps from the nikon nd2 files, the timestamps are then synchronized with the start of the TTL fire signal acquired with the picoscope so that a global time reference is defined. Motion is measure as the average over the field of view of the absolute value of the frame difference. A baseline is then extracted and an animation is optionally generated allowing a visual impection of the result. Finally, the timestamps and measured motion are saved into a result table that can be saved as filename'-motion.csv'.

The latencyanalysis.m script will look for these files '-motion.csv' and match the corresponding picoscope .mat files located in a picodata/ folder. After synchronisation of the fire TTL and the timestamps, the latency is then defined as the time after a TTL pulse that exhibit some motion amplitude above 2 x S.D of the motion signal.

## step by step analysis procedure

1. open the acquired image sequence (eg filename.nd2) in imagej (as a virtual stack if the file is big) and run the macro Quantify_Motion.ijm (drap n drop the file and press run).
2. check if the default information displayed (frame interval) is correct and tick 'animation' if you want to record  an animation with the motion ampltiude overlayed with the image sequence.
2. save the result table as filename-motion.csv
3. export the picoscope file using picoscope software (https://www.picotech.com/downloads) to matlab .mat files in a folder picodata
4. run the latencynaalysis.m script in matlab or octave and select the folder where the -motion.csv files are located.

