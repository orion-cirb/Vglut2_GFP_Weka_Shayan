package Genes_Tools;


import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.plugin.filter.RankFilters;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Color;
import java.awt.Font;
import java.awt.Point;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import org.apache.commons.io.FilenameUtils;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.lang.ArrayUtils;



/**
 *
 * @author phm
 */

public class Vglut2_GFP_Processing {
    
    public CLIJ2 clij2 = CLIJ2.getInstance();
    
    // min size for dots
    public double minDots = 0.4;
    // max size for dots
    public double maxDots = Double.MAX_VALUE;
    // min size for fibers
    public double minFiber = 5;
    // max size for fiber
    public double maxFiber = Double.MAX_VALUE;
    private boolean noise = false;
    private double minDOGDots = 4;
    private double maxDOGDots = 6;
    private double minDOGFibers = 4;
    private double maxDOGFibers = 6;
    private String[] thresholdMethod = AutoThresholder.getMethods();
    private AutoThresholder.Method FibersThreshold = AutoThresholder.Method.Triangle;
    private Calibration cal = new Calibration();
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
     /**
     * check  installed modules
     * @return 
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
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
     * Dialog
     */
    public int[] dialog(String[] channels) {
        String[] chNames = {"Fibers : ", "Vglut2 : "};
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 20, 0);
        gd.addImage(icon);
        gd.addMessage("Channels selection", Font.getFont("Monospace"), Color.blue);
        for (int n = 0; n < chNames.length; n++) {
            gd.addChoice(chNames[n], channels, channels[0]);
        }
        gd.addMessage("Dots filter", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min dots : ", minDots, 2, 6, "µm3");
        gd.addNumericField("Min DOG : ", minDOGDots);
        gd.addNumericField("Max DOG : ", maxDOGDots);
        gd.addMessage("Fibers filter", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min DOG : ", minDOGFibers);
        gd.addNumericField("Max DOG : ", maxDOGFibers);
        gd.addChoice("Fibers Threshold method :", thresholdMethod, FibersThreshold.toString());
        gd.addCheckbox("Noise filter", noise);
        gd.showDialog();
        int[] chChoices = new int[channels.length];
        for (int n = 0; n < chNames.length; n++) {
            chChoices[n] = ArrayUtils.indexOf(channels, gd.getNextChoice());
        }
        minDots = gd.getNextNumber();
        minDOGDots= gd.getNextNumber();
        maxDOGDots= gd.getNextNumber();
        minDOGFibers = gd.getNextNumber();
        maxDOGFibers = gd.getNextNumber();
        FibersThreshold = AutoThresholder.Method.valueOf(gd.getNextChoice());
        noise = gd.getNextBoolean();
        if (gd.wasCanceled())
                chChoices = null;
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
    
    
    /**
     *
     * @param img
     */
    public void closeImages(ImagePlus img) {
        img.flush();
        img.close();
    }

    
   /**
     * return objects population in an binary image
     * @param img
     * @return pop objects population
     */

    public  Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    } 
    
    
    
    
    /**
     * Difference of Gaussians 
     * Using CLIJ2
     * @param imgCL
     * @param sizeX1
     * @param sizeY1
     * @param sizeX2
     * @param sizeY2

     * @return imgGauss
     */ 
    public ClearCLBuffer DOG(ClearCLBuffer imgCL, double sizeX1, double sizeY1, double sizeZ1, double sizeX2, double sizeY2, double sizeZ2) {
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian3D(imgCL, imgCLDOG, sizeX1, sizeY1, sizeZ1, sizeX2, sizeY2, sizeZ2);
        clij2.release(imgCL);
        return(imgCLDOG);
    }
    
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param imgCL
     * @param thMed
     * @param fill 
     */
    public ClearCLBuffer threshold(ClearCLBuffer imgCL, String thMed) {
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        return(imgCLBin);
    }
    
    /*Median filter 
     * 
     * @param img
     * @param size
     */ 
    public void median_filter(ImagePlus img, double size) {
        RankFilters median = new RankFilters();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            median.rank(img.getProcessor(), size, RankFilters.MEDIAN);
            img.updateAndDraw();
        }
    }
    
    /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public void clearOutSide(ImagePlus img, Roi roi) {
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(roi);
            ip.setBackgroundValue(0);
            ip.setColor(0);
            ip.fillOutside(roi);
        }
        img.updateAndDraw();
    }
    
    
    /**
     * return objects population in an binary image
     * Using CLIJ2
     * @param imgCL
     * @return pop
     */

    private Objects3DPopulation getPopFromClearBuffer(ClearCLBuffer imgCL, Roi roi, double min, double max) {
        ClearCLBuffer output = clij2.create(imgCL);
        clij2.connectedComponentsLabelingBox(imgCL, output);
        ImagePlus imgLab  = clij2.pull(output);
        imgLab.setCalibration(cal);
        if (roi != null) {
            clearOutSide(imgLab, roi);
        }   
        ImageInt labels = new ImageLabeller().getLabels(ImageHandler.wrap(imgLab));
        labels.setCalibration(cal);
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        Objects3DPopulation popFilter = new Objects3DPopulation(pop.getObjectsWithinVolume(min, max, true));
        clij2.release(output);
        pop = null;
        return (popFilter);
    }   

    /**
     * Find pml population
     * @param imgGene
     * @return pmlPop
     */
    public Objects3DPopulation findVglut2Pop(ImagePlus imgGene, Roi roiCrop) {
        IJ.showStatus("Finding gene dots ...");
        if (noise) {
            median_filter(imgGene, 2);
        }
        ClearCLBuffer imgCLMed = clij2.push(imgGene);
        ClearCLBuffer imgCLDOG = DOG(imgCLMed, minDOGDots, minDOGDots, minDOGDots, maxDOGDots, maxDOGDots, maxDOGDots);
        clij2.release(imgCLMed);
        ClearCLBuffer imgCLBin = threshold(imgCLDOG, "Otsu"); 
        clij2.release(imgCLDOG);
        Objects3DPopulation pmlPop = getPopFromClearBuffer(imgCLBin, roiCrop, minDots, maxDots);
        clij2.release(imgCLBin);       
        return(pmlPop);
    }
    
    /**
     * Nucleus segmentation
     * @param imgNuc
     * @return 
     */
    public Objects3DPopulation findFibers(ImagePlus imgGFP, Roi roiCrop) {
        IJ.showStatus("Finding fibers ...");
        int med = 1;
        if (noise) {
            med = 2;
        }
        median_filter(imgGFP, med);
        ClearCLBuffer imgCL = clij2.push(imgGFP);
        ClearCLBuffer imgCLDOG = DOG(imgCL, minDOGFibers, minDOGFibers, minDOGFibers/2, maxDOGFibers, maxDOGFibers, maxDOGFibers/2);
        clij2.release(imgCL);
        ClearCLBuffer imgCLBin = threshold(imgCLDOG, FibersThreshold.toString()); 
        clij2.release(imgCLDOG);
        Objects3DPopulation fibersPop = getPopFromClearBuffer(imgCLBin, roiCrop, minFiber, maxFiber);
        clij2.release(imgCLBin);   
        return(fibersPop);
    }
    
   
    
    /**
     * Find dots population in fibers
     * @param cellsPop
     * @param dotsPop
     * @return dotsFiberPop
     */
    public Objects3DPopulation findDotsPop(Objects3DPopulation fibersPop, Objects3DPopulation dotsPop) {
        Objects3DPopulation dotsInFibers = new Objects3DPopulation();
        for (int i = 0; i < fibersPop.getNbObjects(); i++) {
            Object3D fiberObj = fibersPop.getObject(i);
            for (int c = 0; c < dotsPop.getNbObjects(); c++) {
                Object3D dotObj = dotsPop.getObject(c);
                if (dotObj.hasOneVoxelColoc(fiberObj))
                    dotsInFibers.addObject(dotObj);
            }
        } 
        return(dotsInFibers);
    }

    /***
     * Roi Volume
     */
    public double roiVolume(ImagePlus img, Roi roi) {
        Calibration cal = img.getCalibration();
        img.setRoi(roi);
        ImageStatistics stats = img.getStatistics(Measurements.AREA);
        double depth = img.getNSlices() / cal.pixelDepth;
        double vol = stats.area * depth;
        return(vol);
    }
    
    /**
     * Find distance
     */
    public double distance(Object3D obj, Point pt1, Point pt2) {
        Double dist = Double.MAX_VALUE;
        LinkedList<Voxel3D> contour = obj.getContours();
        double x2x1 = pt2.x - pt1.x;
        double y2y1 = pt2.y - pt1.y;
        double pt1pt2 = Math.sqrt(Math.pow(x2x1, 2) + Math.pow(y2y1, 2));
        // find min distance
        for (Voxel3D vox : contour) {
            double d = Math.abs(x2x1*(pt1.y - vox.y) - (pt1.x - vox.x)*y2y1) / pt1pt2;
            if (d <= dist)
                dist = d;
        }
        return(dist*cal.pixelWidth); 
    }
    
    
    /**
     * Find distance from vglut extremity to line
     */
    public ArrayList<Double> findDistances(Roi roiLine, Objects3DPopulation dotsPop) {
        ArrayList<Double> dist = new ArrayList<>();
        Point[] points = roiLine.getContainedPoints();
        for (int i = 0; i < dotsPop.getNbObjects(); i++) {
            Object3D vglutObj =  dotsPop.getObject(i);
            dist.add(distance(vglutObj, points[0], points[1]));
        }
        return(dist);
    }

    /**
     * Label object
     * @param popObj
     * @param img 
     */
    public void labelsObject (Objects3DPopulation popObj, ImagePlus img, int fontSize) {
        Font tagFont = new Font("SansSerif", Font.PLAIN, fontSize);
        String name;
        for (int n = 0; n < popObj.getNbObjects(); n++) {
            Object3D obj = popObj.getObject(n);
            name = Integer.toString(n+1);
            int[] box = obj.getBoundingBox();
            int z = (int)obj.getCenterZ();
            int x = box[0] - 1;
            int y = box[2] - 1;
            img.setSlice(z+1);
            ImageProcessor ip = img.getProcessor();
            ip.setFont(tagFont);
            ip.setColor(255);
            ip.drawString(name, x, y);
            img.updateAndDraw();
        }
    }
    

    
   
}
