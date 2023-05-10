# Individualized neonatal functional connectivity

This repository contains the code for  

> Derman D., Pham D. D., Mejia A. F., & Ferradal S. L. (2023). Surfaced-based Bayesian modeling reveals individual patterns of functional connectivity in neonates. ~INSERT DOI LINK~ (preprint)

## Overview

The code is organized in different files as per the block diagram present in supplementary figure 1 below. 

![readme_figure.png](/readme_figure.png)

0. Minimally preprocessed MRI data is available from the [second release of the developing Human Connectome Project (dHCP)](http://www.developingconnectome.org/data-release/second-data-release/information-registration-and-download/) database. See details in [Fitzgibbon et al., 2020](https://doi.org/10.1016/j.neuroimage.2020.117303).
1. Frame censoring.
2. Volume-to-surface mapping was performed first to the individual surface in native space and then to the 40-week atlas. This was based on the [script](https://git.fmrib.ox.ac.uk/seanf/dhcp-neonatal-fmri-pipeline/-/blob/master/dhcp/func/hcp_surface.sh) by the DHCP group, in turn a modified version of the [HCP pipeline](https://github.com/Washington-University/HCPpipelines/blob/master/fMRISurface/scripts/RibbonVolumeToSurfaceMapping.sh). 
3. group ICA was performed using the [melodic](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MELODIC) tool within the FSL toolbox. 
4. Dual regression was calculated by the `dual_reg()` function within the [templateICAr library](https://github.com/mandymejia/templateICAr/) and checked with the [FSL Toolbox](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/DualRegression) implementation.
5. Template estimation.
6. Individual Bayesian inference and masks.

--------------------------------------------------------------------------
- [ ] Citations
- [ ] References
- [ ] Summary of the paper?
