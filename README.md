# Vglut2_GFP_Weka_Shayan

**Developed for:** Shayan

**Team:** Selimi

**Date:** October 2022

**Software:** Fiji

### Images description

3D images taken with a x60 objective

2 channels:
  1. *CSU_488:* GFP (not mandatory)
  2. *CSU_488:* Vglut2

With each image should be provided a *.roi* or *.zip* file containing one or multiple ROI(s).

### Plugin description

* If required, normalize images with quantile-based normalization

* Detect Vglut2 dots with Weka + watershed

* If GFP channel provided, detect GFP signal with Weka + colocalize it with Vglut2 dots


### Dependencies

* **3DImageSuite** Fiji plugin

* **Weka pixel classifier** model named *classifier.model*

* If GFP channel provided, **Weka pixel classifier** model named *classifier-GFP.model*
 

### Version history

Version 1 released on October 17, 2022.
