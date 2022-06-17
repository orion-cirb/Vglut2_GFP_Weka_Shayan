/*
 * Find vglut2 dots in GFP fiber 
 * Author Philippe Mailly 
 */



import Genes_Tools.Vglut2_GFP_Processing;
import ij.*;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.frame.RoiManager;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import java.util.ArrayList;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;


public class Vglut2_GFP implements PlugIn {

    private final boolean canceled = false;
    public String outDirResults = "";
    public Calibration cal = new Calibration();    
    private String imageDir = "";

    private Genes_Tools.Vglut2_GFP_Processing genes = new Vglut2_GFP_Processing();
      

    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir += IJ.getDirectory("Choose images directory")+File.separator;
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            ArrayList<String> imageFiles = genes.findImages(imageDir, "nd");
            if (imageFiles == null) {
                return;
            }
            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find chanels, image calibration
            reader.setId(imageFiles.get(0));
            String[] channels = genes.findChannels(imageFiles.get(0), meta, reader);
            int[] channelIndex = genes.dialog(channels);
            if(channelIndex == null)
                return;
            cal = genes.findImageCalib(meta);

            // create output folder
            outDirResults = inDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // Write headers results for results files
            FileWriter fileResults = new FileWriter(outDirResults + "results.xls", false);
            BufferedWriter outPutResults = new BufferedWriter(fileResults);
            outPutResults.write("ImageName\t#Vglut2 Dot\tDots volume (µm3)\tDots distance to Purkinje (%)\n");
            outPutResults.flush();
            
            FileWriter globalFileResults = new FileWriter(outDirResults + "globalResults.xls", false);
            BufferedWriter outPutGlobalResults = new BufferedWriter(globalFileResults);
            outPutGlobalResults.write("ImageName\tFibers volume (µm3)\tMean Vglut2 Volume\tVglut2 density\tMean Dots distance to Purkinje (%)\tRoi Volume (µm3)\n");
            outPutGlobalResults.flush();
            
            
            int series = 0;
            // Read images
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                
                // find rois
                String roiFile = imageDir+rootName+".roi";
                if (!new File(roiFile).exists()) {
                    roiFile = imageDir+rootName+".zip";
                    if (!new File(roiFile).exists()) 
                        IJ.showStatus("No roi file found !");
                    return;
                }
                RoiManager rm = new RoiManager(false);
                rm.runCommand("Open", roiFile);
                Roi[] rois = rm.getRoisAsArray();
                Roi roiLineL = null;
                Roi roiCrop = null;
                Roi roiLineP = null;
                for (Roi roi : rois) {
                    if (roi.getType() == Roi.LINE)
                        if (roi.getName().equals("P") || roi.getName().equals("p") )
                            roiLineP = roi;
                        else
                            roiLineL = roi;
                    else
                        roiCrop = roi;
                }
                
                reader.setId(f);
                ImporterOptions options = new ImporterOptions();
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setId(f);
                options.setSplitChannels(true);
                
                Objects3DPopulation fibersPop = new Objects3DPopulation();
                ImagePlus imgGFP = null;
                if (channelIndex[0] != 0) {
                    options.setCBegin(series, channelIndex[0] - 1);
                    options.setCEnd(series, channelIndex[0] - 1);
                
                    // Open GFP Channel detect fiber
                    imgGFP = BF.openImagePlus(options)[0];
                    fibersPop = genes.findFibers(imgGFP, roiCrop);
                    System.out.println("Total fibers found = "+fibersPop.getNbObjects());
                }

                // Open vglut image
                options.setCBegin(series, channelIndex[1] - 1);
                options.setCEnd(series, channelIndex[1] - 1);
                
                ImagePlus imgVglut2 = BF.openImagePlus(options)[0];
                Objects3DPopulation vglut2Pop = genes.findVglut2Pop(imgVglut2, roiCrop);
                System.out.println("Total vglut2 dots = "+vglut2Pop.getNbObjects());
                
                Objects3DPopulation vglut2InFibers = new Objects3DPopulation();
                if (channelIndex[0] != 0) {
                    // Find vglut2 touching fibers
                    vglut2InFibers = genes.findDotsPop(fibersPop, vglut2Pop);
                    System.out.println("Total vglut2 dots in fibers = "+vglut2InFibers.getNbObjects());
                }
                
                
                // find roi volume
                double roiVol = genes.roiVolume(imgVglut2, roiCrop);
                ArrayList<Double> distP = new ArrayList<>();
                ArrayList<Double> distL = new ArrayList<>();
                
                // find distance dots to lineP
                if (roiLineL != null || roiLineP != null) {
                    if (channelIndex[0] != 0) {
                        distP = genes.findDistances(roiLineP, vglut2InFibers);
                        distL = genes.findDistances(roiLineL, vglut2InFibers);
                    }
                    else {
                        distP = genes.findDistances(roiLineP, vglut2Pop);
                        distL = genes.findDistances(roiLineL, vglut2Pop);
                    }
                }
                // Write parameters
                IJ.showStatus("Writing parameters ...");
                DescriptiveStatistics distPourStat = new DescriptiveStatistics();
                DescriptiveStatistics vGlut2VolStats = new DescriptiveStatistics();
                

                // Vglut2 density
                double vglutDensity = 0;
                       
                //fiberVol;
                double fiberVol = 0;
                double distPour = 0;
                // If fiber exists
                if (channelIndex[0] != 0) {
                    // Vglut2 density in fibers
                    vglutDensity = vglut2InFibers.getNbObjects();
                    // Fibers volume
                    for (int i = 0; i < fibersPop.getNbObjects(); i++)
                        fiberVol += fibersPop.getObject(i).getVolumeUnit();
                    for (int i = 0; i < vglut2InFibers.getNbObjects(); i++) {
                        Object3D dotObj = vglut2InFibers.getObject(i);
                        double vglut2Vol = dotObj.getVolumeUnit();
                        vGlut2VolStats.addValue(vglut2Vol);
                        if (roiLineL != null || roiLineP != null) {
                            distPour = (distP.get(i)/(distL.get(i)+distP.get(i)))*100;
                            distPourStat.addValue(distPour);
                        }
                        outPutResults.write(rootName+"\t"+(i+1)+"\t"+vglut2Vol+"\t"+distPour+"\n");
                        outPutResults.flush();
                    }
                    outPutGlobalResults.write(rootName+"\t"+fiberVol+"\t"+vGlut2VolStats.getMean()+"\t"+vglutDensity+"\t"+distPourStat.getMean()+"\t"
                            +roiVol+"\n");
                    outPutGlobalResults.flush();

                    // Save objects image
                    ImageHandler imhFiber = ImageHandler.wrap(imgGFP).createSameDimensions();
                    // draw fibers
                    fibersPop.draw(imhFiber, 255);
                    ImageHandler imhDots = ImageHandler.wrap(imgVglut2).createSameDimensions();
                    vglut2InFibers.draw(imhDots, 255);
                    ImagePlus[] imgColors = {imhDots.getImagePlus(), imhFiber.getImagePlus(), null, imgGFP};
                    ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                    imgObjects.setCalibration(cal);
                    FileSaver ImgObjectsFile = new FileSaver(imgObjects);
                    ImgObjectsFile.saveAsTiff(outDirResults+rootName+"_Objects.tif"); 
                    imhFiber.closeImagePlus();
                    imhDots.closeImagePlus();
                    genes.closeImages(imgGFP);
                    genes.closeImages(imgVglut2);
                }
                else {
                    vglutDensity = vglut2Pop.getNbObjects();
                    for (int i = 0; i < vglut2Pop.getNbObjects(); i++) {
                        Object3D dotObj = vglut2Pop.getObject(i);
                        double vglut2Vol = dotObj.getVolumeUnit();
                        vGlut2VolStats.addValue(vglut2Vol);
                         if (roiLineL != null || roiLineP != null) {
                             distPour = (distP.get(i)/(distL.get(i)+distP.get(i)))*100;
                             distPourStat.addValue(distPour);
                         }
                        outPutResults.write(rootName+"\t"+(i+1)+"\t"+vglut2Vol+"\t"+distPour+"\n");
                        outPutResults.flush();
                    }
                    outPutGlobalResults.write(rootName+"\t"+fiberVol+"\t"+vGlut2VolStats.getMean()+"\t"+vglutDensity+"\t"+distPourStat.getMean()+"\t"
                            +roiVol+"\n");
                    outPutGlobalResults.flush();
                    // Save objects image
                    ImageHandler imhDots = ImageHandler.wrap(imgVglut2).createSameDimensions();
                    vglut2Pop.draw(imhDots, 255);
                    ImagePlus[] imgColors = {imhDots.getImagePlus(), null, null, imgVglut2};
                    ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                    imgObjects.setCalibration(cal);
                    FileSaver ImgObjectsFile = new FileSaver(imgObjects);
                    ImgObjectsFile.saveAsTiff(outDirResults+rootName+"_Objects.tif"); 
                    imhDots.closeImagePlus();
                    genes.closeImages(imgVglut2);
                }
            }
            IJ.showStatus("Process done");
        } catch (DependencyException | ServiceException | FormatException | IOException ex) {
            Logger.getLogger(Vglut2_GFP.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
