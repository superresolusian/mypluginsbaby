package gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.measure.Calibration;
import ij.process.FloatProcessor;
import org.apache.commons.math3.random.RandomDataGenerator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Random;

public class DutyCycleSimulator_ extends _Dialog_{

    LinkedHashMap<String, double[]> fluorophoreDictionary;

    String[] windowTitles;
    String none = "**Use built in test-structure**";
    String defaultIm;
    String imageGT;

    String fluorophore, buffer;
    int experimentLength, frameDuration;
    double wavelength, photons, dutyCycle, survivalFraction, pixelSizeGT;

    int w, h;
    int xc, yc;
    int rad;
    int nParticles;

    Random random = new Random();
    RandomDataGenerator rnd = new RandomDataGenerator();

    @Override
    public boolean beforeSetupDialog(String arg) {
        String[] openWindowTitles = WindowManager.getImageTitles();

        windowTitles = new String[openWindowTitles.length+1];
        windowTitles[0] = none;

        for(int i=1; i<windowTitles.length; i++){
            windowTitles[i] = openWindowTitles[i-1];
        }

        if(windowTitles.length>1) defaultIm = windowTitles[1];
        else defaultIm = none;

        return true;
    }

    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Duty cycle simulator");
        gd.addChoice("Fluorophore", fluorophores, getPrefs("fluorophore", fluorophores[0]));
        gd.addChoice("Buffer", buffers, getPrefs("buffer", buffers[0]));
        gd.addNumericField("Experiment length (seconds)", getPrefs("experimentLength", 180), 0);
        gd.addNumericField("Frame duration (milliseconds)", getPrefs("frameDuration", 30), 0);
        gd.addChoice("Ground truth binary image to simulate", windowTitles, defaultIm);
        gd.addNumericField("GT pixel size in nm (if providing own image)", getPrefs("pixelSizeGT", 10), 0);

    }

    @Override
    public boolean loadSettings() {
        fluorophore = gd.getNextChoice();
        buffer = gd.getNextChoice();
        experimentLength = (int) gd.getNextNumber();
        frameDuration = (int) gd.getNextNumber();
        imageGT = gd.getNextChoice();
        pixelSizeGT = gd.getNextNumber();

        setPrefs("fluorophore", fluorophore);
        setPrefs("buffer", buffer);
        setPrefs("experimentLength", experimentLength);
        setPrefs("frameDuration", frameDuration);
        setPrefs("pixelSizeGT", pixelSizeGT);

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        populateFluorophoreDictionary();
        String key = fluorophore+buffer;

        experimentLength *= 1000; // convert to ms
        double[] properties = fluorophoreDictionary.get(key);
        wavelength = properties[0];
        photons = properties[1];
        dutyCycle = properties[2];
        survivalFraction = properties[3];

        int nFrames = experimentLength/frameDuration;

        double bleachGradient = (survivalFraction-1)/(400*1000);

        int[] xPos = new int[0];
        int[] yPos = new int[0];

        if(imageGT==none){
            w = 501;
            h = 501;
            xc = 250;
            yc = 250;
            rad = 50;
            nParticles = 50;
        }
        else{
            ImagePlus impGT = WindowManager.getImage(imageGT);
            FloatProcessor fpGT = impGT.getProcessor().convertToFloatProcessor();
            w = impGT.getWidth();
            h = impGT.getHeight();
            ArrayList<int[]> locs = getParticles(fpGT);
            xPos = locs.get(0);
            yPos = locs.get(1);
        }

        ArrayList<float[]> switchingLists = new ArrayList<float[]>(nParticles);

        double meanNSwitches = 0;

        for(int n=0; n<nParticles; n++){
            IJ.showStatus("getting trace for molecule "+(n+1)+"/"+nParticles+"...");
            ArrayList<Integer> thisStateMap = new ArrayList<Integer>();
            int currentState = 1;
            int timeInState = 0;

            for(int t=0; t<experimentLength; t++){
                // deal with on molecules
                if(currentState==1){
                    double diceRoll = random.nextDouble();

                    // see if molecule turned off
                    if(diceRoll>dutyCycle){
                        //log.msg("n "+(n+1)+":"+timeInState);
                        thisStateMap.add(timeInState);
                        currentState = 0;
                        timeInState = 0;
                        continue;
                    }

                    // see if molecule bleached
                    else{
                        // get probability of bleaching
                        double probBleach = t*bleachGradient + 1;
                        //log.msg("probBleach = "+probBleach);
                        diceRoll = random.nextDouble();

                        if(diceRoll>probBleach){
                            // molecule bleached
//                            log.msg("n "+(n+1)+":"+timeInState);
//                            log.msg("n "+(n+1)+" bleached");
                            thisStateMap.add(timeInState);
                            break;
                        }

                        // molecule did not bleach
                        timeInState++;
                    }
                }

                // deal with off molecules
                else{
                    double diceRoll = random.nextDouble();

                    // see if molecule turned on
                    if(diceRoll<dutyCycle){
                        //log.msg("n "+(n+1)+":"+(-timeInState));
                        thisStateMap.add(-timeInState);
                        currentState = 1;
                        timeInState = 0;
                    }
                    else{
                        timeInState++;
                    }
                }
            }

            int counter = 0;
            int f=0;
            float[] photonTrace = new float[nFrames];
            meanNSwitches += thisStateMap.size();

            // convert map into list of photons at each frame
            while(f<nFrames){
                int thisState = thisStateMap.get(counter);
                //log.msg("f="+f+", state="+thisState);

                // off molecules
                if(thisState<0){
                    int framesOff = -thisState/frameDuration;
                    framesOff = Math.max(framesOff, 1);
                    f += framesOff;
                }

                // on molecules
                else{
                    int framesOn = thisState/frameDuration;
                    framesOn = Math.max(framesOn, 1);
                    float photonsPerFrame = (float) (photons/framesOn);
                    for(int f_=f; f_<f+framesOn; f_++){
                        photonTrace[f_] = photonsPerFrame;
                    }
                    f += framesOn;
                }
                counter++;
                if(counter==thisStateMap.size()) break;
            }

            switchingLists.add(photonTrace);
        }


        if(imageGT==none) {
            // grab molecule coordinates
            double thetaInc = 2 * Math.PI / nParticles;
            xPos = new int[nParticles];
            yPos = new int[nParticles];

            for (int n = 0; n < nParticles; n++) {
                xPos[n] = (int) (xc + rad * Math.cos(thetaInc * n));
                yPos[n] = (int) (yc + rad * Math.sin(thetaInc * n));
            }
        }


        double fwhm = wavelength/2.8;
        double sigma = (fwhm/2.35482)/10;

        int w_, h_;

        if(imageGT==none) {
            w_ = 50;
            h_ = 50;
        }
        else{
            w_ = (int) ((w*pixelSizeGT)/100);
            h_ = (int) ((h*pixelSizeGT)/100);
        }

        ImageStack ims = new ImageStack(w_, h_, nFrames);

        for(int f=0; f<nFrames; f++){
            IJ.showStatus("rendering frame "+(f+1)+"/"+nFrames);
            FloatProcessor fp = new FloatProcessor(w, h);

            for(int n=0; n<nParticles; n++){
                int x = xPos[n];
                int y = yPos[n];
                float v = switchingLists.get(n)[f];
                // poisson noise
                if(v>0){
                    v = (float) rnd.nextPoisson(v);
                }
                fp.setf(x, y, v);
            }

            fp.blurGaussian(sigma);
            // assume QE = 1 => 1 photon = 1 electron
            // EMCCD gain = 100
            fp.multiply(100);
            // downsample
            fp = fp.resize(w_, h_).convertToFloatProcessor();
            // add read noise std 50
            float[] pixels = (float[]) fp.getPixels();
            for(int i=0; i<pixels.length; i++) pixels[i] = (float) Math.max(0, pixels[i] + random.nextGaussian()*50);
            fp.setPixels(pixels);
            // analog to digital conversion
            fp.add(100); //base of 100
            // EM noise sqrt2
            fp.multiply(1.4);

            ims.setProcessor(fp, f+1);
        }

        //IJ.log("mean switching cycles = "+meanNSwitches/nParticles);

        ImagePlus impBlinking = new ImagePlus(fluorophore+" in "+buffer, ims);
        Calibration calibration = new Calibration();
        calibration.setUnit("micron");
        calibration.pixelWidth = 0.1;
        calibration.pixelHeight = 0.1;
        calibration.setTimeUnit("ms");
        calibration.frameInterval = frameDuration;
        impBlinking.setCalibration(calibration);

        impBlinking.show();


    }

    ArrayList<int[]> getParticles(FloatProcessor fp){
        ArrayList<Integer> xLocs = new ArrayList<Integer>();
        ArrayList<Integer> yLocs = new ArrayList<Integer>();
        for(int y=0; y<h; y++){
            for(int x=0; x<w; x++){
                if(fp.getf(x,y)>0){
                    xLocs.add(x);
                    yLocs.add(y);
                }
            }
        }

        nParticles = xLocs.size();
        int[] xLocs_ = new int[nParticles];
        int[] yLocs_ = new int[nParticles];

        for(int i=0; i<nParticles; i++){
            xLocs_[i] = xLocs.get(i);
            yLocs_[i] = yLocs.get(i);
        }

        ArrayList<int[]> output = new ArrayList<int[]>(2);
        output.add(xLocs_);
        output.add(yLocs_);
        return output;
    }

    public void populateFluorophoreDictionary(){
        fluorophoreDictionary = new LinkedHashMap<String, double[]>();
        fluorophoreDictionary.put(fluorophores[0]+buffers[0], new double[]{519, 1193, 0.00055, 0.94});
        fluorophoreDictionary.put(fluorophores[0]+buffers[1], new double[]{519,  427,  0.0017,    1});
        fluorophoreDictionary.put(fluorophores[1]+buffers[0], new double[]{603, 2826, 0.00058, 0.58});
        fluorophoreDictionary.put(fluorophores[1]+buffers[1], new double[]{603, 1686,  0.0027, 0.99});
        fluorophoreDictionary.put(fluorophores[2]+buffers[0], new double[]{665, 3823,  0.0005, 0.83});
        fluorophoreDictionary.put(fluorophores[2]+buffers[1], new double[]{665, 5202,  0.0012, 0.73});
        fluorophoreDictionary.put(fluorophores[3]+buffers[0], new double[]{775,  437, 0.00006, 0.36});
        fluorophoreDictionary.put(fluorophores[3]+buffers[1], new double[]{775,  703,  0.0001, 0.68});
        fluorophoreDictionary.put(fluorophores[4]+buffers[0], new double[]{810,  591, 0.00049, 0.54});
        fluorophoreDictionary.put(fluorophores[4]+buffers[1], new double[]{810,  740,  0.0014, 0.62});
        fluorophoreDictionary.put(fluorophores[5]+buffers[0], new double[]{523, 1341, 0.00065, 0.98});
        fluorophoreDictionary.put(fluorophores[5]+buffers[1], new double[]{523, 1110,  0.0022, 0.99});
        fluorophoreDictionary.put(fluorophores[6]+buffers[0], new double[]{538, 1231,  0.0015, 0.92});
        fluorophoreDictionary.put(fluorophores[6]+buffers[1], new double[]{538,  868, 0.00062, 0.86});
        fluorophoreDictionary.put(fluorophores[7]+buffers[0], new double[]{592,19714, 0.00058, 0.17});
        fluorophoreDictionary.put(fluorophores[7]+buffers[1], new double[]{592,13294, 0.00037, 0.55});
        fluorophoreDictionary.put(fluorophores[8]+buffers[0], new double[]{669, 1526,  0.0021, 0.46});
        fluorophoreDictionary.put(fluorophores[8]+buffers[1], new double[]{669,  944,  0.0016, 0.84});
        fluorophoreDictionary.put(fluorophores[9]+buffers[0], new double[]{669, 3254,  0.0012, 0.24});
        fluorophoreDictionary.put(fluorophores[9]+buffers[1], new double[]{669, 4433,  0.0035, 0.65});
        fluorophoreDictionary.put(fluorophores[10]+buffers[0], new double[]{684, 1105,  0.0006, 0.65});
        fluorophoreDictionary.put(fluorophores[10]+buffers[1], new double[]{684,  657,  0.0011, 0.78});
        fluorophoreDictionary.put(fluorophores[11]+buffers[0], new double[]{700, 1656,  0.0019, 0.65});
        fluorophoreDictionary.put(fluorophores[11]+buffers[1], new double[]{700,  987,  0.0024, 0.91});
        fluorophoreDictionary.put(fluorophores[12]+buffers[0], new double[]{764,  779, 0.00047, 0.31});
        fluorophoreDictionary.put(fluorophores[12]+buffers[1], new double[]{764,  463,  0.0014, 0.96});
        fluorophoreDictionary.put(fluorophores[13]+buffers[0], new double[]{506, 6241, 0.00012, 0.12});
        fluorophoreDictionary.put(fluorophores[13]+buffers[1], new double[]{506, 4583, 0.00045, 0.19});
        fluorophoreDictionary.put(fluorophores[14]+buffers[0], new double[]{570, 1365,  0.0003,    1});
        fluorophoreDictionary.put(fluorophores[14]+buffers[1], new double[]{570, 2057,  0.0004, 0.89});
        fluorophoreDictionary.put(fluorophores[15]+buffers[0], new double[]{570,11022,  0.0001, 0.17});
        fluorophoreDictionary.put(fluorophores[15]+buffers[1], new double[]{570, 8158,  0.0003, 0.55});
        fluorophoreDictionary.put(fluorophores[16]+buffers[0], new double[]{596, 4968,  0.0017, 0.89});
        fluorophoreDictionary.put(fluorophores[16]+buffers[1], new double[]{596, 8028,  0.0005, 0.61});
        fluorophoreDictionary.put(fluorophores[17]+buffers[0], new double[]{670, 4254,  0.0004, 0.75});
        fluorophoreDictionary.put(fluorophores[17]+buffers[1], new double[]{670, 5873,  0.0007, 0.83});
        fluorophoreDictionary.put(fluorophores[18]+buffers[0], new double[]{694, 5831,  0.0069, 0.87});
        fluorophoreDictionary.put(fluorophores[18]+buffers[1], new double[]{694, 6337,  0.0073, 0.85});
        fluorophoreDictionary.put(fluorophores[19]+buffers[0], new double[]{776,  852,  0.0003, 0.48});
        fluorophoreDictionary.put(fluorophores[19]+buffers[1], new double[]{776,  997,  0.0004, 0.49});
        fluorophoreDictionary.put(fluorophores[20]+buffers[0], new double[]{778,  712,  0.0006, 0.55});
        fluorophoreDictionary.put(fluorophores[20]+buffers[1], new double[]{778,  749,  0.0002, 0.58});
        fluorophoreDictionary.put(fluorophores[21]+buffers[0], new double[]{675, 3653,  0.0011, 0.79});
        fluorophoreDictionary.put(fluorophores[21]+buffers[1], new double[]{675, 3014,  0.0018, 0.64});
        fluorophoreDictionary.put(fluorophores[22]+buffers[0], new double[]{518, 1493, 0.00032, 0.51});
        fluorophoreDictionary.put(fluorophores[22]+buffers[1], new double[]{518,  776, 0.00034, 0.83});
        fluorophoreDictionary.put(fluorophores[23]+buffers[0], new double[]{518,  639, 0.00041, 0.75});
        fluorophoreDictionary.put(fluorophores[23]+buffers[1], new double[]{518, 1086, 0.00031, 0.90});
        fluorophoreDictionary.put(fluorophores[24]+buffers[0], new double[]{794, 2753,  0.0018, 0.60});
        fluorophoreDictionary.put(fluorophores[24]+buffers[1], new double[]{794, 2540,   0.038,    1});
        fluorophoreDictionary.put(fluorophores[25]+buffers[0], new double[]{575, 4884,  0.0017, 0.85});
        fluorophoreDictionary.put(fluorophores[25]+buffers[1], new double[]{575, 2025,  0.0049, 0.99});
    }

    String[] fluorophores =
            {"Alexa Fluor 488", //0
                    "Alexa Fluor 568", //1
                    "Alexa Fluor 647", //2
                    "Alexa Fluor 750", //3
                    "Alexa Fluor 790", //4
                    "Atto 488",
                    "Atto 520",
                    "Atto 565",
                    "Atto 647",
                    "Atto 647N",
                    "Atto 655",
                    "Atto 680",
                    "Atto 740",
                    "Cy2",
                    "Cy3B",
                    "Cy3",
                    "Cy3.5",
                    "Cy5",
                    "Cy5.5",
                    "Cy7",
                    "DyLight 750",
                    "Dyomics 654",
                    "Fluorescein",
                    "FITC",
                    "IRDye 800 CW",
                    "TAMRA"
            };

    String[] buffers = new String[] {"MEA", "BME"};
}