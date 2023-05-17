#' Population priors
#' 
#' This script estimate the template of population priors on DHCP data (2nd release)
#' Should run in parallel by: rel2_estimate-template_term.py
#' 
#'
#' @author Diego Derman, FRG - IUB
#' v1: 2021-08-12
#' v2: 2022-02-01
#' v3: 2022-02-16
#' duration: 7 minutes for 35 subjects, using 48 threads which are leverage by WB.
#' 
#' input:
#' @param template_fname the filename root for the output template files.
#' @param ica_file the filename for the input group ICA cifti file.
#' @param normA boolean: Wheather the A matrix is normalized in between regressions, see Nickerson et al., 2017 appendix
#' @param var_method: "nn" or "ub", non-negative or unbiased variance estimation
#' @param scale: "global" or "local", smooth standard deviation estimates if local
#' @param bold_scale_fwhm:'If not empty, the FWHM for the local area of BOLD normalization by local variance. By default it is NULL and there is no normalization of input BOLD timeseries' 
#' @param smooth_fwhm: Smoothing kernel for the output maps of the kernel, in mm
#' @param cifti_names: filenames for each individual BOLD file composing the template
#' 
#' output:
#' Writes three filetypes with the same root filename, the R object, the cifti files without smoothing and the optional cifti files with smoothing.
#' 
#' @import dhcpy


# command line arguments ####
args = commandArgs(trailingOnly=TRUE)
cifti_fnames = args[-(1:8)]
smooth_fwhm = as.double(args[7])
bold_scale_fwhm = as.double(args[6])
scale = args[5]
var_method = args[4]
normA = as.logical(args[3])
ica_file = args[2]
template_fname = args[1]
verbose = TRUE

if (file.exists(paste0(template_fname, "_smooth-", smooth_fwhm , ".RDS"))){
  cat("Template exists, skipping...")
  quit(status=0)
}

# parse smoothing arguments
smoothing = FALSE
if (as.logical(smooth_fwhm)){
  smoothing = TRUE
} else {
  if (file.exists(paste0(template_fname, ".RDS"))){
  cat("Template exists, skipping...")
  quit(status=0)
  }
}
cat(paste0("Smoothing is ", smoothing, " kernel: ", smooth_fwhm, "\n"))

cat(paste0("sd scaling is ", scale, " kernel: ", bold_scale_fwhm, "\n"))


# load dependencies ####
# load utilities
devtools::load_all(path=paste0("./utils/dhpcr"), quiet = TRUE)
library("ciftiTools")
# set path for workbench tools
ciftiTools.setOption('wb_path', wb_path())
library("templateICAr")
check_tica('0.4.0')

# set up dhcp default surfaces
dhcpr_env = environment()
default_surfaces(dhcpr_env)

# estimate template ####
# calculate
if (verbose) print(paste0("scale=",scale, "\nscale_sm_FWHM=", bold_scale_fwhm))
template = estimate_template(cifti_fnames,
                             GICA = ica_file,
                             normA = normA,
                             varTol = 1e-100, # very important, this could generate NAs in template.
                             scale = scale,
                             scale_sm_FWHM = bold_scale_fwhm,
                             #scale_sm_surfL = dhcpr_env$surfL,
                             scale_sm_surfL = midthicknessL_path(),
                             scale_sm_surfR = midthicknessR_path())

# save results to output filename
saveRDS(template, file = paste0(template_fname, ".RDS"))
export_template(template, out_fname=paste0(template_fname, "_object"))

# smoothing ####

# debugging
template = readRDS(paste0(template_fname, ".RDS"))

if (smoothing){
  # Resample the data.
  template$template <- lapply(template$template, function(y){
    as.matrix(smooth_cifti(newdata_xifti(template$dat_struct, y), surf_FWHM = smooth_fwhm, surfL_fname=midthicknessL_path(), surfR_fname=midthicknessR_path()))
  })
  
  # Replace `NaN` values with NA values.
  template$template <- lapply(template$template, function(y){y[] <- ifelse(is.nan(y), NA, y)})
  
  
  # save smoothed ####  
  
  saveRDS(template, file = paste0(template_fname, "_smooth-", smooth_fwhm , ".RDS"))
}


