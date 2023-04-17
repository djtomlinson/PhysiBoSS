# PhysiBoSS 2.0: a sustainable integration of stochastic Boolean and agent-based modelling frameworks

**Version:** 2.2.0

**Release date:** March 2022

## Overview: 
PhysiBoSS 2.0 is a redesign and reimplementation of PhysiBoSS ([doi:10.1093/bioinformatics/bty766](https://doi.org/10.1093/bioinformatics/bty766)). It has been conceived as an add-on that expands the PhysiCell ([doi:10.1371/journal.pcbi.1005991](https://dx.doi.org/10.1371/journal.pcbi.1005991)) agent-based functionalities with intracellular cell signalling using MaBoSS having a decoupled, maintainable and model-agnostic design. PhysiBoSS 2.0 reproduces simulations reported in the original PhysiBoSS publications and can be used with other Boolean models, for instance to predict drug synergy in a gastric adenocarcinoma cell line.

**Reference paper:** [Ponce-de-Leon, Miguel, Arnau Montagud, Vincent Noël, Gerard Pradas, Annika Meert, Emmanuel Barillot, Laurence Calzone, and Alfonso Valencia. "PhysiBoSS 2.0: a sustainable integration of stochastic Boolean and agent-based modelling frameworks." bioRxiv (2022): 2022-01.](https://www.biorxiv.org/content/10.1101/2022.01.06.468363v3)

**Reference paper doi:** [10.1101/2022.01.06.468363v3](https://www.biorxiv.org/content/10.1101/2022.01.06.468363v3)

### How to run a PhysiBoSS sample_project inside PhysiCell:
~~~bash
git clone https://github.com/PhysiBoSS/PhysiBoSS.git
cd PhysiBoSS
make physiboss-tnf-model
make
./spheroid_TNF_model
~~~

### Key makefile rules, from [PhysiCell repository](https://github.com/MathCancer/PhysiCell):

**`make`**: compiles the current project. If no 
                     project has been defined, it first 
                     populates the cancer heterogeneity 2D 
                     sample project and compiles it 
   
**`make project-name`**: populates the indicated sample project. 
                     Use "make" to compile it. 

   * **`project-name`** choices:
      * template 
      * biorobots-sample 
      * cancer-biorobots-sample 
      * cancer-immune-sample
      * celltypes3-sample 
      * heterogeneity-sample 
      * pred-prey-farmer 
      * virus-macrophage-sample 
      * worm-sample
      * ode-energy-sample 
      * physiboss-cell-lines-sample 
      * cancer-metabolism-sample
      * interaction-sample
      * mechano-sample

**`make list-projects`** : list all available sample projects 

**`make clean`**         : removes all .o files and the executable, so that the next "make" recompiles the entire project 

**`make data-cleanup`**  : clears out all simulation data 

**`make reset`**         : de-populates the sample project and returns to the original PhysiCell state. Use this when switching to a new PhysiCell sample project. 

**`make jpeg`**          : uses ImageMagick to convert the SVG files in the output directory to JPG (with appropriate sizing to make movies). Supply `OUTPUT=foldername` to select a different folder. 

**`make movie`**         : uses ffmpeg to convert the JPG files in the output directory an mp4 movie. Supply `OUTPUT=foldername` to select a different folder, or `FRAMERATE=framerate` to override the frame rate.

**`make upgrade`**       : fetch the latest release of PhysiCell and overwrite the core library and sample projects. 

