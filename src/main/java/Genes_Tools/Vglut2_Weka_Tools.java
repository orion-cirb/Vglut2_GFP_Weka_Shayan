package Genes_Tools;


import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import ij.plugin.RGBStackMerge;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;



public class Vglut2_Weka_Tools {
        
    // min size for dots
    public double minDots = 1;
    public double minDotsGFP = 5;
    public int slicemin = 1;
    public int slicemax = 1;
    // max size for dots
    public double maxDots = Double.MAX_VALUE;
    public double maxDotsGFP = Double.MAX_VALUE;
    private Calibration cal = new Calibration();
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    public boolean doPreprocess;
    public boolean doWeka;
    public boolean weka3D;
     /**
     * check  installed modules
     * @return 
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }

        
   /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No Image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
      /**
     * Find images in folder
     */
    public String findWekaModel(String imagesFolder, boolean gfp) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No Model found in "+imagesFolder);
            return null;
        }
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals("model")) {
                if (gfp) {
                    if (f.endsWith("-GFP.model")) return(imagesFolder + File.separator + f);
                }
                else {
                    if (!f.endsWith("-GFP.model")) return(imagesFolder + File.separator + f);
                }
            }
        }
        return"";
    }
    
    public void setNSlice(int z) { slicemax = z; }
    
    /**
     * Dialog
     */
    public int[] dialog(String[] channels) {
        String[] chNames = {"GFP : ", "Vglut2 : "};
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 15, 0);
        gd.addImage(icon);
        
        gd.addMessage("Channels selection", Font.getFont("Monospace"), Color.blue);
        for (int n = 0; n < chNames.length; n++) {
            gd.addChoice(chNames[n], channels, channels[0]);
        }
        gd.addMessage("Preprocessing (crop+normalisation)", Font.getFont("Monospace"), Color.blue);
        gd.addCheckbox("Preprocess images", true);
        gd.addMessage("Slices selection", Font.getFont("Monospace"), Color.black);
        gd.addNumericField("Min slice : ", slicemin, 0);
        gd.addToSameRow();
        gd.addNumericField("Max slice : ", slicemax, 0);
        gd.addMessage("Weka segmentation", Font.getFont("Monospace"), Color.blue);
        gd.addCheckbox("Do segmentation + results", true);
        gd.addCheckbox("Weka in 3D", true);
        gd.addMessage("Vglut2 dots filter", Font.getFont("Monospace"), Color.black);
        gd.addNumericField("Min dots : ", minDots, 2, 6, "µm3");
        gd.addToSameRow();
        gd.addNumericField("Max dots : ", maxDots, 2, 6, "µm3");
        gd.addMessage("GFP filter", Font.getFont("Monospace"), Color.black);
        gd.addNumericField("Min vol : ", minDotsGFP, 2, 6, "µm3");
        gd.addToSameRow();
        gd.addNumericField("Max vol : ", maxDotsGFP, 2, 6, "µm3");
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("xy calibration (µm)", cal.pixelHeight,3);
        gd.addNumericField("z calibration (µm)", cal.pixelDepth,3);
        
        gd.showDialog();
        int[] chChoices = new int[channels.length];
        for (int n = 0; n < chNames.length; n++) {
            chChoices[n] = ArrayUtils.indexOf(channels, gd.getNextChoice());
        }
        if (gd.wasCanceled())
                chChoices = null;
        doPreprocess = gd.getNextBoolean();
        slicemin = (int) gd.getNextNumber();
        slicemax = (int) gd.getNextNumber();
        doWeka = gd.getNextBoolean();
        weka3D = gd.getNextBoolean();
        minDots = gd.getNextNumber();
        maxDots = gd.getNextNumber();
        minDotsGFP = gd.getNextNumber();
        maxDotsGFP = gd.getNextNumber();
        
        cal.pixelHeight = gd.getNextNumber();
        cal.pixelWidth = cal.pixelHeight;
        cal.pixelDepth = gd.getNextNumber();
        return(chChoices);
    } 
    
    /**
     * Find channels name
     * @param imageName
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs+1];
        channels[0] = "None";
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n+1] = Integer.toString(n);
                    else 
                        channels[n+1] = meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n+1] = Integer.toString(n);
                    else 
                        channels[n+1] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n+1] = Integer.toString(n);
                    else 
                        channels[n+1] = meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n+1] = Integer.toString(n);
                    else 
                        channels[n+1] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;    
            default :
                for (int n = 0; n < chs; n++)
                    channels[n+1] = Integer.toString(n);
        }
        return(channels);         
    }
    
    public double scaleInZ() { return cal.pixelDepth; }
     public double scaleArea() { return cal.pixelWidth*cal.pixelHeight; }
    
     public void setCalibration(ImagePlus imp) {
         imp.setCalibration(cal);
     }
     
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public Calibration findImageCalib(IMetadata meta) {
        // read image calibration
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
        return(cal);
    }
    
      public Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
      
      public double stdArray(double[] arr, double mean, int[] coloc, boolean withcoloc) {
          double res = 0;
          if (arr.length==0) return 0;
          int n=0;
          for (int i = 0; i < arr.length; i++) {
              if (!withcoloc || (coloc[i]==1)) {
                res += (arr[i]-mean)*(arr[i]-mean);
                n++;
              }
        }
          if (n==0) return 0;
        return Math.sqrt(res)/n;
      }
    
    /**
     *
     * @param img
     */
    public void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    public void drawPop(Objects3DPopulation pop, ImagePlus gfp) {
        IJ.run(gfp, "Select All", "");
        IJ.run(gfp, "Clear", "stack");
        ImageHandler imh = ImageHandler.wrap(gfp);
        pop.draw(imh, 255);
        gfp = imh.getImagePlus();
    }
    
    public void drawPopulation(Objects3DPopulation pop, int[] imsize, ImagePlus gfp, String name) {
                ImagePlus red =  IJ.createImage("Vglut2", "8-bit black", imsize[0], imsize[1], 1, imsize[2], 1);
                ImageHandler imh = ImageHandler.wrap(red);
                pop.draw(imh, 255);
                red = imh.getImagePlus();
                ImagePlus res = null;
                if (gfp!=null){
                    ImagePlus[] images = {red, gfp};
                    res = RGBStackMerge.mergeChannels(images, false); 
                } else {
                    res = red;
                }
                res.setCalibration(cal);
               IJ.saveAs(res, "Tiff", name);
               closeImages(res);
    }

    
}
