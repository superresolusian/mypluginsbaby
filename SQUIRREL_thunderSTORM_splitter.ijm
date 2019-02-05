// make dialog
makeDialog();

// grab values from dialog
start = Dialog.getNumber();
end = Dialog.getNumber();
increment = Dialog.getNumber();

topLeftX = Dialog.getNumber();
topLeftY = Dialog.getNumber();
originalWidth = Dialog.getNumber();
originalHeight = Dialog.getNumber();
magnification = Dialog.getNumber();

// set up output stack
width = originalWidth*magnification;
height = originalHeight*magnification;
nFrames = floor((end-start)/increment);
newImage("SR stack", "32-bit black", width, height, nFrames);

// populate stack
for(n=1; n<=nFrames; n++){
	thisEnd = start + (n*increment);
	run("Show results table", "action=filter formula=[frame>"+(start-1)+" & frame<"+thisEnd+"]");
	run("Visualization", "imleft="+topLeftX+" imtop="+topLeftY+" imwidth="+originalWidth+" imheight="+originalHeight+" renderer=[Averaged shifted histograms] magnification="+magnification+" colorize=false threed=false shifts=2");
	run("Select All");
	run("Copy");
	close();
	selectWindow("SR stack");
	setSlice(n);
	setMetadata("Label", "Frames "+start+"-"+(thisEnd-1));
	run("Paste");
}



function makeDialog(){
	Dialog.create("Split up particles table");

	Dialog.addNumber("Start frame", 1);
	Dialog.addNumber("End frame", 10000);
	Dialog.addNumber("Increment", 50);

	Dialog.addNumber("Top left x", 0);
	Dialog.addNumber("Top left y", 0);
	Dialog.addNumber("Original image width", 128);
	Dialog.addNumber("Original image height", 128);
	Dialog.addNumber("Magnification", 5);

	Dialog.show();
}
