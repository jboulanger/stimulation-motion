/*
 *	 Quantify motion using a difference between the 
 *	 image and a running average (3d blur)
 *	 
 *	 The signal is decomposed into 3 components:
 *	 	motion, baseline and stimulus and saved in the results table
 *	 
 *	 1. Load the image sequence using import>avi
 *	 Use virtual stacks and convert to gray scale
 *	 
 *	 2. Select a (small) region of interest and duplicate to crop the image
 *	 
 *	 3. Set the Time interval in Image>Properties
 *	 
 *	 4. Run the macro on the image
 *	 
 *	 Jerome Boulanger (2018) for Stefano Giandomenico
 *	 
 *	 18-05-17 signal decomposition
 *	 
 */
 
macro "Frame Differences" {
	saveSettings();
	dT = Stack.getFrameInterval();	
	if (dT == 0) {
		dT = 1;
	}
	Stack.getUnits(xu, yu, zu, tu, vu);
	Dialog.create("Frame Difference");
	Dialog.addNumber("start frame",1);
	Dialog.addNumber("end frame",nSlices);
	Dialog.addNumber("frame interval",dT);
	Dialog.addString("frame interval unit",tu);
	Dialog.addCheckbox("animation",true);
	Dialog.addNumber("scale",1);	
	Dialog.addNumber("plot max ampltiude ", 1.5);	
	Dialog.show();
	t0 = Dialog.getNumber();
	t1 = Dialog.getNumber();
	dT = Dialog.getNumber();
	tu = Dialog.getString();
	anim = Dialog.getCheckbox();
	binning = Dialog.getNumber();
	ymax = Dialog.getNumber();
	
	setBatchMode(true);	
	N = t1-t0-1;	
	y = newArray(N); // y-axis: motion amplitude over time
	x = newArray(N); // x-axis: time	
	z = newArray(N); // residual motion
	b = newArray(N); // base line
	p = newArray(N); // peaks
	computeImageDifferences(x,y,t0,t1,dT);
	analyzeSignal(y,z,b,p,5);
	if (anim) {
		createAnimation(x,z,b,p,t0,t1,binning,ymax);
	}	
	toResults(x,y,z,b,p,t0,t1);	
	setBatchMode(false);
	restoreSettings;
}

macro "Frame Differences Batch" {
	folder = getDirectory("Input folder");
	list = getFileList(folder);	
	processFolder(folder, list, ".*.nd2")
}

function processFolder(_folder, _list, _pattern) {
	for (i = 0; i < _list.length; i++) {
		if (matches(_list[i],_pattern)) {
			path = _folder + _list[i];
			processFile(path);
		}
	}
}

function processFile(_path) {		
	run("Bio-Formats Importer", "open=["+_path+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack");
	id = getImageID;
	setBatchMode(true);	
	dT = Stack.getFrameInterval();	
	Stack.getUnits(xu, yu, zu, tu, vu);
	t0 = 1;
	t1 = nSlices;
	N = t1-t0-1;	
	y = newArray(N); // y-axis: motion amplitude over time
	x = newArray(N); // x-axis: time	
	z = newArray(N); // residual motion
	b = newArray(N); // base line
	p = newArray(N); // peaks
	computeImageDifferences(x,y,t0,t1,dT);
	analyzeSignal(y,z,b,p,5);	
	toResults(x,y,z,b,p,t0,t1);	
	setBatchMode(false);
	saveAs("Results", replace(_path,".nd2","-motion.csv"));	
	selectWindow("Results");
	run("Close");
	selectImage(id);
	close();
}

// Analyze signal y by detecting jumps above a running median
function analyzeSignal(y,z,b,p,w) {	
	Array.getStatistics(y, ymin, ymax,ymean,ystd);
	// arbitrary large threshold to detect very high motion due to illumination change
	T = ymean + 12 * ystd;
	// compute a baseline a 1d closing	
	b0 = dilate1dsig(erode1dsig(y,w),w);	
	// split the signal in 3 components
	for (i = 0; i < y.length; i++) {	
		b[i] = b0[i]; // the baseline
		if (y[i] >  T) { // if large change
			p[i] = y[i] - b[i]; // peak
			z[i] = 0;			
		} else {
			p[i] = 0; 
			z[i] = y[i] - b[i]; // residual motion
		}
	}
}

// 1d signal erosion of y by a struturing elemets of length 2*s+1
function erode1dsig(_x,_s) {
	z = newArray(x.length);
	for (i = 0; i < _x.length; i++) {
		j0 = maxOf(i-_s,0);
		j1 = minOf(i+_s,_x.length-1);
		block = Array.slice(_x,j0,j1);
		Array.getStatistics(block, bmin);
		z[i] = bmin;
	}	
	return z;
}

// 1d signal dilation of y by a struturing elemets of length 2*s+1
function dilate1dsig(_x,_s) {
	z = newArray(x.length);
	for (i = 0; i < _x.length; i++) {
		j0 = maxOf(i-_s,0);
		j1 = minOf(i+_s,_x.length-1);
		block = Array.slice(_x,j0,j1);
		Array.getStatistics(block, bmin, bmax);
		z[i] = bmax;
	}	
	return z;
}

// median of an array
function median(_A) {
	B = Array.copy(_A);
	Array.sort(B);
	return B[round(B.length-1)/2];
}

// Median of Absolute Deviation esimation of the standard deviation
// unused and might be incorrect
function MAD(_A) {
	B = Array.copy(_A);
	Array.sort(B);
	m = B[round(B.length-1)/2];		
	for (i = 0; i < B.length; i++) {
		B[i] = abs(_A[i] - m);
	}
	Array.sort(B);
	m = B[round(B.length-1)/2];	
	return 1.4825 * m;
}

// compute image differences and store the results as a 1d signal in (x,y) 
// from frame from t0 to t1 with dT is the time sampling used to compute x.
function computeImageDifferences(x,y,t0,t1,dT) {
	if (endsWith(getTitle(),".nd2")) {
		hwts = getND2TimeStamps();		
		usehwts = true;
	} else {
		usehwts = false;
	}
	id = getImageID;
	setSlice(t0);
	run("Duplicate...", "use title=A");	
	run("32-bit");
	//run("Square Root", "stack");
	for (t = t0+1; t < t1; t++) {		
		selectImage(id);
		setSlice(t+1);
		run("Duplicate...", "use title=B");
		run("32-bit");
		//run("Square Root", "stack");
		imageCalculator("Difference", "A","B");
		selectWindow("A");
		List.setMeasurements();
		if (usehwts) {
			x[t-t0-1] = hwts[t];
		} else {
			x[t-t0-1] = (t-t0)*dT;
		}
		y[t-t0-1] = List.getValue("Mean");		
		selectWindow("A"); close();
		selectWindow("B"); rename("A");
	}
	selectWindow("A"); close();	
	selectImage(id);
}

// store the results into the table results
function toResults(x,y,z,b,p,t0,t1){
	for (t = t0; t < t1-1; t++) {
		setResult("Time",t-t0,x[t-t0]);
		setResult("Signal",t-t0,y[t-t0]);
		setResult("Motion",t-t0,z[t-t0]);
		setResult("Baseline",t-t0,b[t-t0]);
		setResult("Stimulus",t-t0,p[t-t0]);
	}
	updateResults();
}

// create an animation with the graph blended on the image sequence
function createAnimation(x,y,b,p,t0,t1,binning,ymax) {
	id = getImageID;
	name0 = getTitle;
 	// size of the original image
	w0 = getWidth;
	h0 = getHeight;
	ytmp = newArray(y.length);	
	for (i = 0; i < y.length; i++) {
		ytmp[i] = y[i] + b[i];	
	}
	if (ymax < 0) {
		Array.getStatistics(ytmp, ymin, ymax);
	} else {
		ymin = 0;		
	}
	Array.getStatistics(b, bmin, bmax);
	Array.getStatistics(x, xmin, xmax);
	ptmp = newArray(p.length);	
	for (i = 0; i < p.length; i++) {
		if (p[i] > 0) {
			ptmp[i] = ymax;
		} else {
			ptmp[i] = ymin;
		}		
	}	
	yline=newArray(ymin,ymax);
	// size of the destination image
	w3=round(w0/binning);
	h3=round(h0/binning);
	setPasteMode("Blend");
	// width and height added to the graph by the axis
	waxis=83;
	haxis=55;	
	for (t = t0; t < t1-1; t++) {
		// create the plot
		xline=newArray(x[t-t0],x[t-t0]);		
		Plot.create("Motion", "time ["+tu+"]", "Motion");		
		Plot.setColor("#ec7063");		
		Plot.setLineWidth(5);
		Plot.add("line", xline, yline);
		Plot.setColor("#aa00ff");
		Plot.setLineWidth(1);
		Plot.add("line", x, ytmp);
		Plot.setColor("#5dade2");
		Plot.add("line", x, b);
		Plot.setColor("#00ffff");
		Plot.add("line", x, ptmp);		
		Plot.setLimits(xmin,xmax,ymin,ymax);
		Plot.setFrameSize(w3-waxis, h3/4-haxis);		
		Plot.show();
		Plot.makeHighResolution("Motion_HiRes-1",2.0);		
		run("Scale...", "x="+w3/getWidth+" y="+0.25*h3/getHeight+" width="+w3+" height="+0.25*h3+" bicubic average create");
		// size of the plot to paste
		w1 = getWidth();
		h1 = getHeight();		
		run("Macro...", "code=v=255-v");
		run("Copy");
		
		// extract a frame from the original movie
		selectImage(id);		
		setSlice(t);
		run("Select None");		
		run("Duplicate...", "use");				
		id2=getImageID;			
		// scale it if needed
		if (abs(binning - 1) > 0.1) {
			run("Scale...", "x="+(1/binning)+" y="+(1/binning)+" width="+w3+" height="+h3+" interpolation=Bicubic average create");
		}
		
		if (t==t0) {		
			id3=getImageID;
			rename("Movie");
			run("RGB Color");			
			makeRectangle(0, getHeight-h1, w1, h1);
			run("Paste");
		} else {			
			id4=getImageID;
			rename("Frame");
			run("RGB Color");			
			makeRectangle(0, getHeight-h1, w1, h1);	
			run("Paste");
			run("Concatenate...", "  title=[Movie] image1=[Movie] image2=[Frame]");			
		}
	}	
	selectWindow("Movie");
	rename("Motion analysis of "+ name0);
}

// old macro (use too much memory)
macro "Quantifiy Motion" {
	setBatchMode(true);
	name0=getTitle;
	run("Duplicate...","duplicate title=tmp");	
	run("Gaussian Blur 3D...", "x=0 y=1 z=10");
	imageCalculator("Difference create 32-bit stack", name0, "tmp");	
	rename("difference");
	selectWindow("difference");	
	makeRectangle(10, 10, getHeight-11, getWidth-11);
	run("Plot Z-axis Profile");
	Plot.setXYLabels("Time [frame]", "Motion [a.u.]");
	rename("Motion vs Time +[" + name0+"]");
	setBatchMode(false);
}

// Split a string into an array of substrings and keep the 
function getND2TimeStamps() {	
	_src = getMetadata("Info");
	_separator = "\n";
	_pattern = ".*timestamp.*";
	i = 0;
	j = 0;	
	k = 0;	
	dst = newArray(nSlices);
	while(i < lengthOf(_src)) {
		j = indexOf(_src, _separator, i);
		if (j == -1) {
			j = lengthOf(_src);
		}		
		line = substring(_src,i,j);
		if (matches(line,_pattern)) {			
			v = split(line, "=");
			value = parseFloat(v[1]);			
			dst[k] = value;			
			k = k + 1;
		}
		i = j + 1;	
	}
	return dst;
}

