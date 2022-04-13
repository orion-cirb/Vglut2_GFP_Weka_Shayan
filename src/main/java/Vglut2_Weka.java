/*
 * Find vglut2 dots in GFP, segmentation with Weka plugin
 * Author Orion cirb
 */



import Genes_Tools.Vglut2_Weka_Tools;
import ij.*;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.awt.Rectangle;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import java.util.ArrayList;
import loci.common.Region;
import loci.common.services.ServiceFactory;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import org.apache.commons.io.FilenameUtils;
import trainableSegmentation.WekaSegmentation;


public class Vglut2_Weka implements PlugIn {

    private final boolean canceled = false;
    private String imageDir = "";
    RoiManager rm;

    private Genes_Tools.Vglut2_Weka_Tools genes = new Vglut2_Weka_Tools();
    private ImageProcessorReader reader = new ImageProcessorReader();
             
    private int[] channelIndex;
    private String processDir;
    public String resultsDir = "";
    
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
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find chanels, image calibration
            reader.setId(imageFiles.get(0));
            String[] channels = genes.findChannels(imageFiles.get(0), meta, reader);
            genes.setNSlice(reader.getSizeZ());
            genes.findImageCalib(meta);
            channelIndex = genes.dialog(channels);
            if(channelIndex == null)
                return;
            

            // create output folder for preprocessed normalised file
            processDir = inDir + File.separator+ "Preprocessed"+ File.separator;
            File procDir = new File(processDir);
            if (!Files.exists(Paths.get(processDir))) {
                procDir.mkdir();
            }
            
            // create output folder
            resultsDir = inDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(resultsDir);
            if (!Files.exists(Paths.get(resultsDir))) {
                outDir.mkdir();
            }
            
            // Write headers results for results files
            FileWriter fileResults = new FileWriter(resultsDir + "results.xls", false);
            BufferedWriter outPutResults = new BufferedWriter(fileResults);
            outPutResults.write("ImageName\tRoi\tVglut2 Dot index\tDot volume (um3)\tDot with GFP\n");
            outPutResults.flush();
            
            FileWriter globalFileResults = new FileWriter(resultsDir + "globalResults.xls", false);
            BufferedWriter outPutGlobalResults = new BufferedWriter(globalFileResults);
            outPutGlobalResults.write("ImageName\tRoi\tRoi Volume (um3)\tGFP Volume (um3)\tNb total Vglut2 dots\tMean Vol Vglut2\tStd Vol Vglut2\tDensityVglut2 in Roi(nb/um3)\tDensityVglut2 in GFP");
            outPutGlobalResults.write("\t Nb Vglut2 dots with GFP\t Mean Vol Vglut2 with GFP\t Std Vol Vglut2-GFP\tDensity Vglut2-GFP in Roi\t Density Vglut2-GFP in GFP \n");
            outPutGlobalResults.flush();      
            rm = new RoiManager(false);
                
            // Preprocess images
            String[] fileList = new String[imageFiles.size()];
            int i = 0;
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                String roiName = "";
                if (genes.doPreprocess) roiName = preprocessFile(f, rootName);
                fileList[i] = rootName+"_"+roiName;
                i++;
            }
             
            // Normalised all images
            if (genes.doPreprocess) {
                IJ.showStatus("Normalisation starting...");
                QuantileBasedNormalization qbn = new QuantileBasedNormalization();
                qbn.run(processDir, fileList, "");
                if (channelIndex[0] !=0 ) qbn.run(processDir, fileList, "-GFP");
                IJ.showStatus("Normalisation done");
            }
            
            // do Weka on all files
            IJ.showStatus("Segmentation starting...");
            if (genes.doWeka) goWeka(fileList, outPutResults, outPutGlobalResults);
            IJ.showStatus("Finished");
            outPutResults.close();
            outPutGlobalResults.close();
           
           } catch (Exception ex) {
            Logger.getLogger(Vglut2_Weka.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public String preprocessFile(String f, String rootName) throws Exception {
        String roiName = ""; 
        try {
        // find rois
            String roiFile = imageDir+rootName+".roi";
            if (!new File(roiFile).exists()) {
                roiFile = imageDir+rootName+".zip";
                if (!new File(roiFile).exists()) {
                    IJ.showStatus("No roi file found !");
                    return("");
                }
            }
            rm.reset();
            rm.runCommand("Open", roiFile);
            for (Roi roi : rm.getRoisAsArray()) {
                roiName = roi.getName();
                reader.setId(f);
                ImporterOptions options = new ImporterOptions();
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setId(f);
                options.setSplitChannels(true);
                options.setCrop(true);
                Region reg = new Region(roi.getBounds().x, roi.getBounds().y, roi.getBounds().width, roi.getBounds().height);
                options.setCropRegion(0, reg);
                int Zbegin = (genes.slicemin > reader.getSizeZ()) ? 1 :  genes.slicemin;
                int Zend = (genes.slicemax > reader.getSizeZ()) ? reader.getSizeZ() :  genes.slicemax;
                options.setZBegin(0, Zbegin+1);
                options.setZEnd(0, Zend-1);
                options.doCrop();
                // Open and crop GFP if any
                if (channelIndex[0] != 0) {
                    options.setCBegin(0, channelIndex[0] - 1);
                    options.setCEnd(0, channelIndex[0] - 1);
                    
                    ImagePlus imgGFP = BF.openImagePlus(options)[0];
                    genes.setCalibration(imgGFP);
//                        imgGFP.setRoi(roi);
//                        ImagePlus croppedGFP = imgGFP.crop("1-"+imgGFP.getNSlices());
//                        genes.closeImages(imgGFP);            
//                        ImagePlus img = croppedGFP.crop(""+genes.slicemin+"-"+genes.slicemax);
//                        genes.closeImages(croppedGFP);
                    IJ.saveAs(imgGFP, "Tiff", processDir+rootName+"_"+roiName+"-GFP.tif");
                    genes.closeImages(imgGFP);
                }

                // Open vglut image
                options.setCBegin(0, channelIndex[1] - 1);
                options.setCEnd(0, channelIndex[1] - 1);

                ImagePlus imgVglut = BF.openImagePlus(options)[0];
                genes.setCalibration(imgVglut);
//                imgVglut.setRoi(roi);
//                ImagePlus cropped = imgVglut.crop("1-"+imgVglut.getNSlices());
                roi.setLocation(0, 0); // update roi to cropped image
                rm.setRoi(roi, 0);
                rm.runCommand("Save", processDir+rootName+"_"+roiName+".roi");
//                genes.closeImages(imgVglut);            
//                ImagePlus img = cropped.crop(""+genes.slicemin+"-"+genes.slicemax);
//                genes.closeImages(cropped);
                IJ.saveAs(imgVglut, "Tiff", processDir+rootName+"_"+roiName+".tif");
                genes.closeImages(imgVglut);
            }
        }
        catch (Exception e) { throw e; }
        return(roiName);
    }
    
    public void goWeka(String[] files, BufferedWriter outRes, BufferedWriter globRes){
        // do weka on each file
        String model = genes.findWekaModel(imageDir, false);
        String modelGFP = genes.findWekaModel(imageDir, true);
        
        for (String f: files) {
            File file = new File(processDir+f+"-normalized.tif");
            if (file.exists()){
                // load Roi
                rm.reset();
                rm.runCommand("Open", processDir+f+".roi");
                Roi roi = rm.getRoi(0);
                String roiName = rm.getName();
                // Segment GFP channel with Weka
                ImagePlus resGFP = null;
                double volGFP = 0;
                 if (channelIndex[0]!=0){
                    ImagePlus impGFP = IJ.openImage(processDir+f+"-GFP-normalized.tif");
                    WekaSegmentation wekaGFP = new WekaSegmentation(genes.weka3D);    
                    wekaGFP.setTrainingImage(impGFP);
                    wekaGFP.loadClassifier(modelGFP);
                    wekaGFP.applyClassifier(false);
                    resGFP = wekaGFP.getClassifiedImage();
                    genes.closeImages(impGFP);
                    IJ.setAutoThreshold(resGFP, "MaxEntropy dark stack");
                    IJ.run(resGFP, "Convert to Mask", "method=MaxEntropy background=Dark");
                    if (resGFP.isInvertedLut()) IJ.run(resGFP, "Invert LUT", "");
                    IJ.run(resGFP, "Options...", "iterations=1 count=1 black do=Nothing");
                    IJ.run(resGFP, "Close-", "stack");
                    IJ.run(resGFP, "Open", "stack");
                    wekaGFP = null;
                    
                    // keep only dots inside Roi
                    resGFP.setRoi(roi);
                    IJ.run(resGFP, "Clear Outside", "stack");
                    // get dots 3D objects
                    Objects3DPopulation gfpPop = genes.getPopFromImage(resGFP);  
                    // remove small objects
                    Objects3DPopulation gfpDots = new Objects3DPopulation(gfpPop.getObjectsWithinVolume(genes.minDotsGFP, genes.maxDotsGFP, true));
                   
                    // calculate GFP volume in ROI
                    for ( int o=0; o<gfpDots.getNbObjects(); o++) {
                        Object3D obj = gfpDots.getObject(o);
                        volGFP += obj.getVolumeUnit();
                    } 
                    genes.drawPop(gfpDots, resGFP);
                }
                
                ImagePlus imp = IJ.openImage(processDir+f+"-normalized.tif");
                WekaSegmentation weka = new WekaSegmentation(genes.weka3D);
                    
                weka.setTrainingImage(imp);
                weka.loadClassifier(model);
                weka.applyClassifier(false);
                ImagePlus res = weka.getClassifiedImage();
                genes.closeImages(imp);
                IJ.setAutoThreshold(res, "MaxEntropy dark stack");
                IJ.run(res, "Convert to Mask", "method=MaxEntropy background=Dark");
                if (res.isInvertedLut()) IJ.run(res, "Invert LUT", "");
                IJ.run(res, "Options...", "iterations=1 count=1 black do=Nothing");
                IJ.run(res, "Open", "stack");
                IJ.run(res, "Watershed", "stack");
                
                // keep only dots inside Roi
                res.setRoi(roi);
                IJ.run(res, "Clear Outside", "stack");
                res.setRoi(roi);
                double roiVol = (res.getStatistics().area) * (res.getNSlices()*genes.scaleInZ());
                if (volGFP == 0) volGFP = roiVol;
                
                // get dots 3D objects
                Objects3DPopulation alldots = genes.getPopFromImage(res);     
                int[] imgsize = {res.getWidth(), res.getHeight(), res.getNSlices()};
                genes.closeImages(res);
   
                try{
                    Objects3DPopulation dots = new Objects3DPopulation(alldots.getObjectsWithinVolume(genes.minDots, genes.maxDots, true));
                    alldots = null;
                    double meanVol = 0;
                    double meanVolGFP = 0;
                    int nGFP = 0;
                    double[] vols = new double[dots.getNbObjects()];
                    int[] colocs = new int[vols.length];
                    for ( int i=0; i < dots.getNbObjects(); i++ ) {
                        Object3D dot = dots.getObject(i);
                        double dotVol = dot.getVolumeUnit();
                        boolean coloc = false;
                        if ( (resGFP!=null) && (dot.getPixMeanValue(ImageHandler.wrap(resGFP)) > 1) ) coloc = true;
                        outRes.write(""+f+"\t"+i+"\t"+dotVol+"\t"+(coloc?"Yes":"No")+"\n");
                        outRes.flush();
                        if (coloc) {
                            meanVolGFP += dotVol;
                            nGFP++;
                        }
                        meanVol += dotVol;
                        vols[i] = dotVol;
                        colocs[i] = coloc?1:0;
                    }
                    meanVol /= vols.length;
                    if (nGFP>0) meanVolGFP /= nGFP;
                    double stdVol = genes.stdArray(vols, meanVol, colocs, false);
                    double stdVolGFP = genes.stdArray(vols, meanVolGFP, colocs, true);
                    globRes.write(""+f+"\t"+roiName+"\t"+roiVol+"\t"+volGFP+"\t"+vols.length+"\t"+meanVol+"\t"+stdVol+"\t"+(vols.length/roiVol)+"\t"+(vols.length)/volGFP+
                        "\t"+nGFP+"\t"+meanVolGFP+"\t"+stdVolGFP+"\t"+(nGFP/roiVol)+"\t"+(nGFP/volGFP)+"\n");
                    globRes.flush();
                
                    // draw objects
                    genes.drawPopulation(dots, imgsize, resGFP, resultsDir+f+"-results.tif");
                    if (resGFP !=null) genes.closeImages(resGFP);
                    
                } catch (Exception ex) {
                    Logger.getLogger(Vglut2_Weka.class.getName()).log(Level.SEVERE, null, ex);
                }
               
                
            }
        }
    }
}
