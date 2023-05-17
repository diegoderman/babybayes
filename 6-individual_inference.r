#' Bayesian inference
#' 
#' This script estimate the individual FC mean and standard error maps of a subject, based on a previous calculated population prior on DHCP data (2nd release)
#' It is written to be run in parallelized fashion
#' *warning*: Each thread will allocate ~100 GB of memory.#' 
#' 125 s per subject using 10 threads.
#'
#' @author Diego Derman, FRG - IUB
#' 
#' CLI input:
#' @param subid DHCP Subject ID.
#' @param sesid DHCP Session ID.
#' @param ticadir OUTDIR where all results (mean and var cifti) will be saved
#' @param template_name: name of the input template
#' @param surfdir: path of the surfaces directory. Assumes DHCP sub-ses canonical file structure.
#' @param cifti_fname: Suffix of the individual BOLD cifti file. Assumes dHCP sub-ses file structure. If smoothed include here.
#' @param kernel: Smoothing kernel for the output maps of the kernel, in mm
#' @param toggle_dr: toggle calculate DR instead of Bayesian ICA with templateICA method.
#' @param normA: toggle normalization of the A matrix in between regressions of DR.
#' @param scale: Whether to normalize BOLD input. options: GLOBAL, local, none.
#' @param bold_scale_fwhm: If not empty, the FWHM for the local area of BOLD normalization by local variance. By default it is NULL and there is no normalization of input BOLD timeseries.
#' @param smooth_fwhm: Smoothing kernel for the output maps of the kernel, in mm.
#' @param tis: kernel in mm FWHM of the ciftis that are input to template estimation.
#' @param verbose: verbosity toggle.

#' output:
#' 
#' @import dhcpy
#' 


# cli input arguments ####
args = commandArgs(trailingOnly=TRUE)
subid = args[1]
sesid = args[2]
ticadir = args[3] # all results will be saved with my path here */sub-*/subses
template_name = args[4]
templatedir = args[5]
surfdir = args[6]
cifti_fname = args[7]
kernel = args[8]
toggle_dr = as.logical(args[9])
normA = as.logical(args[10])
scale=args[11]
bold_scale_fwhm = as.double(args[12]) # fwhm of local var scaling
smooth_fwhm = as.double(args[13])
tis = as.double(args[14])
verbose = as.logical(args[15])

# dependencies #####

## setting environment variables for workbench **you may need to change this according to your scheduler**
Sys.setenv(OMP_NUM_THREADS = "2")

## library unixtools to change tempdir
library("unixtools")
## library stringi for random strings
library("stringi")
random_string = stri_rand_strings(1, 5)
# create and set a different temp directory 
dir.create(paste0("./temp/", random_string), showWarnings = FALSE)
set.tempdir(paste0("./temp/", random_string))

## library dhcpr and load all paths
devtools::load_all("./dhcpr", quiet = TRUE)
## library ciftitools and set workbench
library("ciftiTools")
wbpath = wb_path()
ciftiTools.setOption('wb_path', wbpath)
## library template ica
##  you cannot make request to github in a parallel script.
library("templateICAr")
check_tica('0.4.0')

# load default DHCP surfaces
dhcpr_env = environment()
default_surfaces(dhcpr_env)

# build file paths and DHCP specific parameters #####

# root name for output files.
fpart_tICA = paste0('sub-', subid, '_ses-', sesid, '_desc-', template_name, "_dris-", kernel)# replace with , '_bold', kernel
# directories for results
baby_ICA = dhcp_path() # **System dependent, only change once in dhcpr**
resultdir = paste0("./results")
## bold input filename
bold_fname = paste0(surfdir, "/sub-", subid, "/ses-", sesid, "/", cifti_fname) # already includes ,".dtseries.nii"
## template input filename assuming smoothed template
template_mean_fname = paste0(templatedir, "/", template_name, "_mean_smooth-2.dscalar.nii")
template_var_fname = paste0(templatedir, "/", template_name, "_var_smooth-2.dscalar.nii")
template_fname = paste0(templatedir, "/", template_name, ".RDS")
# by default cifti_templateica needs a prefix to which it needs to add _var and _mean
# I therefore manually modified the file names from above to comply:
invisible(file.copy(template_mean_fname, paste0(templatedir, "/", template_name, "_smooth-2_mean.dscalar.nii"), overwrite = TRUE))
invisible(file.copy(template_var_fname, paste0(templatedir, "/", template_name, "_smooth-2_var.dscalar.nii"), overwrite = TRUE))
if (smooth_fwhm != 0){
    template_fname = paste0(templatedir, "/", template_name, "_smooth-", smooth_fwhm)
} else {
    template_fname = paste0(templatedir, "/", template_name)
}
## create dir for results
subresultdir = paste0(ticadir, '/sub-', subid)
dir.create(subresultdir, showWarnings = FALSE)
## set midthickness surfaces
surfaceL_fname = midthicknessL_path()
surfaceR_fname = midthicknessR_path()

# Main Bayesian inference ####

if (file.exists(paste0(subresultdir, '/', fpart_tICA,'.RDS'))){
  # check if rds exists, otherwise go on 
  if (verbose) cat(paste0(subid, ": RDS found, loading RDS...\n"))
  error("Output RDS file exists, exiting.")

} else {

# otherwise, prepare for Bayesian inference proper 
  
# read bold file
BOLD <- read_cifti(bold_fname)

# get number of voxels in file
V <- nrow(do.call(rbind, BOLD$data))

# get surfaces
Lsurf = read_surf(surfaceL_fname)
Rsurf = read_surf(surfaceR_fname)
# for symmetrical atlas, number of vertices per hemisphere is half of total
hemi_v = V/2

# DUAL REGRESSION (option to check results)  #####
if(toggle_dr){

  ### Get files
  # get previously calculated template
  tmean = read_cifti(template_mean_fname)
  BOLD_mat = as.matrix(BOLD)
  tmean_mat <- rbind(tmean$data$cortex_left, tmean$data$cortex_right)

  # dual regression between subject BOLD and template mean
  DR_shv <- dual_reg(BOLD_mat, tmean_mat, scale=TRUE)

  # Write dual regression results to cifti
  tmp <- t(DR_shv$S)
  dr_cifti = as.xifti(cortexL = tmp[1:hemi_v,], cortexR = tmp[-(1:hemi_v),], cortexL_mwall = BOLD$meta$cortex$medial_wall_mask$left, cortexR_mwall = BOLD$meta$cortex$medial_wall_mask$right , surfL=Lsurf, surfR=Rsurf)
  #write out the xifti
  write_cifti(dr_cifti, paste0(fpart_tICA, "_dual-regression.dscalar.nii"))
  # Save as RDS file as well
  saveRDS(DR_shv, file=paste0(fpart_tICA,'_dual-regression.RDS'))
}

# RUN TEMPLATE ICA #####

if (verbose) cat(paste0(subid, ": loading template...\n"))

template = readRDS(paste0(template_fname, ".RDS"))

if (verbose) cat(paste0(subid, ": templateICA starting...\n"))

# main templateICAr function
tICA <- try(templateICA(
  BOLD=bold_fname,
  template=template,
  tvar_method="non-negative",
  spatial_model=FALSE,
  resamp_res=NULL,
  # Adding local scaling of bold input - new template mode
  scale=scale,
  scale_sm_FWHM = bold_scale_fwhm,
  scale_sm_surfL = midthicknessL_path(),
  scale_sm_surfR = midthicknessR_path(),
  #Q2=NULL,
  maxiter=300,
  epsilon=0.001,
  verbose=FALSE,
  kappa_init=NULL,
  varTol=1e-100,
  normA = normA
))

if("try-error" %in% class(tICA)) stop(paste0(subid, ": ERROR IN TEMPLATE ICA, STOP EXECUTION."))

if (verbose) cat(paste0(subid, ": templateICA completed.\n"))

if (verbose) cat(paste0(subid, ": saving RDS...\n"))
# Save: elapsed time 41s
saveRDS(tICA, file=paste0(subresultdir, '/', fpart_tICA,'.RDS'))
if (verbose) cat(paste0(subid, ": RDS saved\n"))

# SAVE CIFTI RESULTS #################################################

# check that file was not written previously
if (!file.exists(paste0(subresultdir, '/', fpart_tICA, "_mean.dscalar.nii"))){

  if (verbose) cat(paste0(subid, ": writing mean cifti...\n"))
  write_cifti(tICA$subjICmean, paste0(subresultdir, '/', fpart_tICA, "_mean.dscalar.nii"),  surfL_fname = surfaceL_fname, surfR_fname = surfaceR_fname)
  if (verbose) cat(paste0(subid, ": mean cifti saved.\n"))
  
  if (verbose) cat(paste0(subid, ": writing variance cifti...\n"))
  write_cifti(tICA$subjICse, paste0(subresultdir, '/', fpart_tICA, "_var.dscalar.nii"),  surfL_fname = surfaceL_fname, surfR_fname = surfaceR_fname)
  if (verbose) cat(paste0(subid, ": variance cifti completed.\n"))
}
}
# remove temp directory contents.
unlink(paste0("/N/project/baby_ICA/temp/", random_string), recursive = TRUE)
