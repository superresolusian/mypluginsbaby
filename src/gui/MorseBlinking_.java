package gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import org.apache.commons.math3.random.RandomDataGenerator;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Random;

import static java.lang.Math.max;

/**
 * Created by sculley on 18/12/2018.
 */
public class MorseBlinking_ implements PlugIn {

    NonBlockingGenericDialog gd;

    static final double[] A = new double[]{100, 0, 100, 100, 100};
    static final double[] B = new double[]{100, 100, 100, 0, 100, 0, 100, 0, 100};
    static final double[] C = new double[]{100, 100, 100, 0, 100, 0, 100, 100, 100, 0, 100};
    static final double[] D = new double[]{100, 100, 100, 0, 100, 0, 100};
    static final double[] E = new double[]{100};
    static final double[] F = new double[]{100, 0, 100, 0, 100, 100, 100, 0, 100};
    static final double[] G = new double[]{100, 100, 100, 0, 100, 100, 100, 0, 100};
    static final double[] H = new double[]{100, 0, 100, 0, 100, 0, 100};
    static final double[] I = new double[]{100, 0, 100};
    static final double[] J = new double[]{100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100, 100, 100};
    static final double[] K = new double[]{100, 100, 100, 0, 100, 0, 100, 100, 100};
    static final double[] L = new double[]{100, 0, 100, 100, 100, 0, 100, 0, 100};
    static final double[] M = new double[]{100, 100, 100, 0, 100, 100, 100};
    static final double[] N = new double[]{100, 100, 100, 0, 100};
    static final double[] O = new double[]{100, 100, 100, 0, 100, 100, 100, 0, 100, 100, 100};
    static final double[] P = new double[]{100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100};
    static final double[] Q = new double[]{100, 100, 100, 0, 100, 100, 100, 0, 100, 0, 100, 100, 100};
    static final double[] R = new double[]{100, 0, 100, 100, 100, 0, 100};
    static final double[] S = new double[]{100, 0, 100, 0, 100};
    static final double[] T = new double[]{100, 100, 100};
    static final double[] U = new double[]{100, 0, 100, 0, 100, 100, 100};
    static final double[] V = new double[]{100, 0, 100, 0, 100, 0, 100, 100, 100};
    static final double[] W = new double[]{100, 0, 100, 100, 100, 0, 100, 100, 100};
    static final double[] X = new double[]{100, 100, 100, 0, 100, 0, 100, 0, 100, 100, 100};
    static final double[] Y = new double[]{100, 100, 100, 0, 100, 0, 100, 100, 100, 0, 100, 100, 100};
    static final double[] Z = new double[]{100, 100, 100, 0, 100, 100, 100, 0, 100, 0, 100};
    static final double[] _1 = new double[]{100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100, 100, 100};
    static final double[] _2 = new double[]{100, 0, 100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100, 100, 100};
    static final double[] _3 = new double[]{100, 0, 100, 0, 100, 0, 100, 100, 100, 0, 100, 100, 100};
    static final double[] _4 = new double[]{100, 0, 100, 0, 100, 0, 100, 0, 100, 100, 100};
    static final double[] _5 = new double[]{100, 0, 100, 0, 100, 0, 100, 0, 100, 0};
    static final double[] _6 = new double[]{100, 100, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0};
    static final double[] _7 = new double[]{100, 100, 100, 0, 100, 100, 100, 0, 100, 0, 100, 0, 100};
    static final double[] _8 = new double[]{100, 100, 100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100, 0, 100, 0};
    static final double[] _9 = new double[]{100, 100, 100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100};
    static final double[] _0 = new double[]{100, 100, 100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100, 100, 100};

    String message;
    LinkedHashMap<String, double[]> dict;
    ArrayList<String> swearDictionary;
    String[] titles;

    Random rand = new Random();
    RandomDataGenerator rdg = new RandomDataGenerator();

    public void setupDialog(){
        gd = new NonBlockingGenericDialog("Inspector Morse");
        gd.addTextAreas("Message to blink", null, 10, 50);
    }

    public boolean loadSettings() {
        message = gd.getNextText();
        return true;
    }

    public void execute() {

        populateSwearDictionary();

        for(int i=0; i<swearDictionary.size(); i++){
            if(message.contains(swearDictionary.get(i))){
                IJ.log("oooh you dirty bastard");
                break;
            }
        }

        populateDictionary();
        String[] lines = message.split("\r\n|\r|\n");
        int nLines = lines.length;
        int[] lineLengths = new int[nLines];
        int longestLine = 0;
        for(int i=0; i<nLines; i++){
            lineLengths[i] = lines[i].length();
            longestLine = max(longestLine, lineLengths[i]);
        }
        int w = longestLine * 51;
        int h = nLines*51;

        // need to work out number of frames - lowest common multiple of the letter lengths
        int[] nMorseLengths = new int[message.length()];
        for(int i=0; i<message.length(); i++){
            String substring = message.substring(i, i+1);
            if(dict.containsKey(substring)){
                nMorseLengths[i] = dict.get(substring).length+3;
            }
            else{
                nMorseLengths[i] = 1;
            }
        }

        int nFrames = lcm(nMorseLengths);
        int nCharacters = message.length() - (nLines-1);

        //create delays
        int[] delays = new int[nCharacters];
        for(int i=0; i<nCharacters; i++){
            int delay = rand.nextInt(3);
            delays[i] = delay;
        }

        ImageStack ims = null;

        for(int n=1; n<=nFrames; n++) {
            IJ.showStatus("Inspecting the Morse... ("+n+"/"+nFrames+")");
            IJ.showProgress(n, nFrames);
            float[] pixels = new float[w * h];

            for (int m = 0; m < nCharacters; m++) {

                //find line this character belongs on
                int lineInd = 0;
                int thisLine = 0;
                for (int l = 0; l < nLines; l++) {
                    lineInd += lineLengths[l];
                    if (m < lineInd) {
                        thisLine = l;
                        break;
                    }
                }

                int thisY = 25 + (51 * thisLine);

                //find column this character belongs on
                int thisColumn = m % lineLengths[thisLine];
                int thisX = 25 + (51 * thisColumn);

                String thisString = lines[thisLine];
                String subString = thisString.substring(thisColumn, thisColumn + 1);

                if (!dict.containsKey(subString)) continue;

                double[] morse = dict.get(subString);

                int rawIndex = (n - 1) % (morse.length + 3);
                int rawIndexDelay = rawIndex - delays[m];

                if (rawIndexDelay < 0 || rawIndexDelay >= morse.length) continue;

                float val = (float) morse[rawIndexDelay];

                pixels[thisX + thisY * w] = val;

            }

            FloatProcessor fp = new FloatProcessor(w, h, pixels);
            fp.blurGaussian(6.75);

            float[] blurPixels = (float[]) fp.getPixels();
            for(int i=0; i<blurPixels.length; i++) blurPixels[i] = blurPixels[i]*500 + 10;

            for(int p=0; p<blurPixels.length; p++){
                float v = blurPixels[p];
                blurPixels[p] = Math.max(rdg.nextPoisson(v), 0);
            }
            fp.setPixels(blurPixels);

            FloatProcessor fpSmall = fp.resize(w/5, h/5).convertToFloatProcessor();
            if(ims==null){
                ims = new ImageStack(fpSmall.getWidth(), fpSmall.getHeight(), nFrames);
            }
            ims.setProcessor(fpSmall, n);
        }

        getListOfTitles();

        new ImagePlus(titles[rand.nextInt(titles.length)], ims).show();

    }

    public int lcm(int[] array){

        int nElements = array.length;
        int lcm_ = array[0];

        for(int i=1; i<nElements; i++){
            int b = array[i];
            lcm_ = (lcm_*b)/gcd(lcm_,b);
        }
        return lcm_;
    }

    public int gcd(int a, int b) {
        if (a < 1 || b < 1) {
            throw new IllegalArgumentException("a or b is less than 1");
        }
        int remainder = 0;
        do {
            remainder = a % b;
            a = b;
            b = remainder;
        } while (b != 0);
        return a;
    }

    public void populateDictionary() {
        dict = new LinkedHashMap<String, double[]>();
        dict.put("a", A);
        dict.put("A", A);
        dict.put("b", B);
        dict.put("B", B);
        dict.put("c", C);
        dict.put("C", C);
        dict.put("d", D);
        dict.put("D", D);
        dict.put("e", E);
        dict.put("E", E);
        dict.put("f", F);
        dict.put("F", F);
        dict.put("g", G);
        dict.put("G", G);
        dict.put("h", H);
        dict.put("i", I);
        dict.put("I", I);
        dict.put("j", J);
        dict.put("J", J);
        dict.put("k", K);
        dict.put("K", K);
        dict.put("l", L);
        dict.put("L", L);
        dict.put("m", M);
        dict.put("M", M);
        dict.put("n", N);
        dict.put("N", N);
        dict.put("o", O);
        dict.put("O", O);
        dict.put("p", P);
        dict.put("P", P);
        dict.put("q", Q);
        dict.put("Q", Q);
        dict.put("r", R);
        dict.put("R", R);
        dict.put("s", S);
        dict.put("S", S);
        dict.put("t", T);
        dict.put("T", T);
        dict.put("u", U);
        dict.put("U", U);
        dict.put("v", V);
        dict.put("V", V);
        dict.put("w", W);
        dict.put("W", W);
        dict.put("x", X);
        dict.put("X", X);
        dict.put("y", Y);
        dict.put("Y", Y);
        dict.put("z", Z);
        dict.put("Z", Z);
        dict.put("0", _0);
        dict.put("1", _1);
        dict.put("2", _2);
        dict.put("3", _3);
        dict.put("4", _4);
        dict.put("5", _5);
        dict.put("6", _6);
        dict.put("7", _7);
        dict.put("8", _8);
        dict.put("9", _9);
        return;
    }

    public void populateSwearDictionary(){
        swearDictionary = new ArrayList<String>();
        swearDictionary.add("shit");
        swearDictionary.add("Shit");
        swearDictionary.add("SHIT");
        swearDictionary.add("fuck");
        swearDictionary.add("Fuck");
        swearDictionary.add("FUCK");
        swearDictionary.add("cunt");
        swearDictionary.add("Cunt");
        swearDictionary.add("CUNT");
        swearDictionary.add("piss");
        swearDictionary.add("Piss");
        swearDictionary.add("PISS");
        swearDictionary.add("cock");
        swearDictionary.add("Cock");
        swearDictionary.add("COCK");
        swearDictionary.add("twat");
        swearDictionary.add("Twat");
        swearDictionary.add("TWAT");
        swearDictionary.add("arse");
        swearDictionary.add("Arse");
        swearDictionary.add("ARSE");
        swearDictionary.add("bollocks");
        swearDictionary.add("Bollocks");
        swearDictionary.add("BOLLOCKS");
    }

    public void getListOfTitles(){
        String title1 = "Full of reMorse";
        String title2 = "Stop Morse-ing around and get back to work";
        String title3 = "I fought the Morse and the Morse won";
        String title4 = "Wild Morses couldn't drag me away";
        String title5 = "Gimme gimme Morse, gimme Morse, gimme gimme Morse";
        String title6 = "Morse Morse Morse, how do you like it, how do you like it?";
        String title7 = "Morse than a feeling";
        String title8 = "Morse than a woman to me";
        String title9 = "Hit me baby one Morse time";
        String title10 = "Sexual interMorse";
        String title11 = "Morse of nature";
        String title12 = "Air Morse one";

        titles = new String[]{title1, title2, title3, title4, title5, title6, title7, title8, title9, title10, title11, title12};

    }

    @Override
    public void run(String s) {
        setupDialog();

        if (gd != null) {
            // Show dialog
            gd.showDialog();
            if (gd.wasCanceled()) {
                return;
            }
        }

        loadSettings();

        execute();
    }



}
