package gui;

import ij.*;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.FolderOpener;
import ij.plugin.PlugIn;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.Map;

/**
 * NanoJ-Core _BaseDialog_ (by Ricardo Henriques) without NanoJ dependencies
 *
 **/

public abstract class _Dialog_ implements PlugIn, DialogListener {

    public NonBlockingGenericDialog gd = null;

    public ImagePlus imp = null;

    protected boolean autoOpenImp = true;

    public String prefsHeader = null, impPath = null;
    protected Prefs prefs = new Prefs();
    protected String arg;

    public void run() {
        run("");
    }

    @Override
    public void run(String arg) {
        this.arg = arg;

        if (!beforeSetupDialog(arg)) return;
        if (autoOpenImp && imp == null) {
            setupImp();
            if (imp == null) return;
        }

        setupDialog();

        if (gd != null) {
            // Add listener to dialog
            gd.addDialogListener(this);
            // Show dialog
            gd.showDialog();
            if (gd.wasCanceled()) {
                return;
            }
        } else {
            if (!loadSettings()) return;
        }

        loadSettings();

        try {
            execute();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        prefs.savePreferences();
    }

    // ~~~~~BASIC DIALOG METHODS~~~~~
    abstract public boolean beforeSetupDialog(String arg);

    abstract public void setupDialog();

    abstract public boolean loadSettings();

    abstract public void execute() throws InterruptedException, IOException;

    protected void setupImp() {
        //get open image processor
        imp = WindowManager.getCurrentImage();

        if (imp == null) {
            openImp();
        }
        if (imp == null) return;
    }

    protected void openImp() {

        impPath = IJ.getFilePath("Choose data to load...");
        if (impPath == null || impPath.equals("")) return;

        else if (impPath.endsWith(".tif")) {
            if (impPath.contains("00")) {
                //open image sequence in directory selected and show
                imp = FolderOpener.open(new File(impPath).getParent());
            } else {
                imp = IJ.openImage(impPath);
            }
        } else {
            imp = IJ.openImage(impPath);
        }
        if (imp == null) return;
        imp.show();
    }

    public String getClassName() {
        return this.getClass().getName();
    }

    // ~~~~~PREFS HANDLING~~~~~

    public int getPrefs(String key, int defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return (int) prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public float getPrefs(String key, float defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return (float) prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public double getPrefs(String key, double defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return (double) prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public boolean getPrefs(String key, boolean defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public String getPrefs(String key, String defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public void setPrefs(String key, int value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, float value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, double value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, boolean value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, String value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void savePrefs() {
        prefs.savePreferences();
    }
}