#' utils.r
#'
#' Utility functions for the dhcp dataset
#'
#'

#' Make thresholded t-statistic maps.
#'
#' Performs one-sample t-test and writes tmap rds, cifti, and png for reporting
#' # SAVE ACTIVATIONS (BIN MASK) AS CIFTI NOT CURRENTLY WORKING (templateICAr 0.2.3)
#'
#' @import ciftiTools
#'
#' @param ses_pair DHCP subject session ids pair
#' @param template_name DHCP template name of the input dual regression estimates
#' @param save_png toggle to save tstat and binary images. Default: \code{FALSE}
#' @param save_rds path to save RDS file out of \code{templateICAr::activations()}
#' @param method_p multiple comparison p correction method from \code{templateICAr::activations()}
#' @param alpha threshold value for one-sample t-test
#' @param env an environment with the default surfaces, to save time.
#'  If null, will load dhcp symatlas surf.
#'
#'
#' @return Void, save images if set.
#'
#' @export

tmap = function(subid, sesid, template_name, save_png = FALSE, method_p = "fdr", alpha = 0.05, env = NULL, save_thresholded = TRUE){



  # get surfaces and mwall
  if (is.null(env)){
    env = environment()
    default_surfaces(env)
    default_mwall(env)
  }

  rds_fname = paste0("/N/project/baby_ICA/partial_products/individual_t-test/sub-",
                     subid, "/", template_fullname, "/sub-", subid, "_ses-", sesid, "_desc_",
                     template_fullname, "_", method_p, "-", alpha, ".RDS")

  if (file.exists(rds_fname)){

  cat(paste0("Activations RDS for ", subid, " exists. Skipping...\n"))
  quit(status=0)
  }


  # results filename
  sub_tica_rds = paste0(dhcp_path(), "/template_ICA/sub-", subid, "/sub-", subid, "_ses-", sesid, "_desc-", template_name, ".RDS")
  # load tica individual results (mean and var maps)
  sub_tica_results = readRDS(sub_tica_rds)
  # get activations
  activations = templateICAr::activations(sub_tica_results, method_p = method_p, alpha = alpha)
  # save rds
  dir.create(paste0("/N/project/baby_ICA/partial_products/individual_t-test/sub-",
                    subid, "/", template_fullname), recursive=TRUE, showWarnings = FALSE)
  saveRDS(activations, rds_fname)


  #create subdir
  dir.create(paste0(dhcp_path(), "/partial_products/individual_t-test/sub-", subid), showWarnings = FALSE)
  # create template dir
  dir.create(paste0(dhcp_path(), "/partial_products/individual_t-test/sub-", subid, "/", template_name), showWarnings = FALSE)


  # cifti save fname
  tstats_xifti_fname = paste0("/N/project/baby_ICA/partial_products/individual_t-test/sub-", subid,
                              "/", template_name, "/sub-", subid, "_ses-", sesid, "_desc_", template_name, "_tstat.dscalar.nii")
  #activations_xifti_fname = paste0("/N/project/baby_ICA/partial_products/individual_t-test/sub-", subid, "/", template_name, "sub-", subid, "_ses-", sesid, "_desc_", template_name, "_activations.dscalar.nii")
  activations_xifti_fname = paste0("/N/project/baby_ICA/partial_products/individual_t-test/sub-", subid,
                                   "/", template_name, "/sub-", subid, "_ses-", sesid, "_desc_", template_name, "_", method_p, "-", alpha, ".dscalar.nii")


  #tstats from matrix to xifti
  tstats_xifti = as.xifti(cortexL = activations$tstats[1:(V()/2),],
                          cortexR = activations$tstats[-(1:V()/2),],
                          cortexL_mwall = env$mwallL,
                          cortexR_mwall = env$mwallR, surfL = env$surfL, surfR = env$surfR)
  # write cifti surfaces

  if(!file.exists(tstats_xifti_fname)) write_cifti(tstats_xifti, cifti_fname = tstats_xifti_fname,
                                                   surfL_fname = midthicknessL_path(), surfR_fname = midthicknessR_path())
  write_cifti(activations$active, cifti_fname = activations_xifti_fname,
              surfL_fname = midthicknessL_path(), surfR_fname = midthicknessR_path())

  if (save_thresholded){

    tstats_xifti_thresholded = as.xifti(cortexL = activations$tstats[1:(V()/2),] * activations$active$data$cortex_left,
                                        cortexR = activations$tstats[-(1:V()/2),] * activations$active$data$cortex_right,
                                        cortexL_mwall = env$mwallL,
                                        cortexR_mwall = env$mwallR, surfL = env$surfL, surfR = env$surfR)

    tstat_thresholded_fname  = paste0(dhcp_path(), "/partial_products/individual_t-test/sub-",
                                      subid, "/", template_name, "/sub-", subid, "_ses-", sesid, "_desc_",
                                      template_name, "_tstat-", method_p, "-", alpha, ".dscalar.nii")

    write_cifti(tstats_xifti_thresholded, cifti_fname = tstat_thresholded_fname,
                surfL_fname = midthicknessL_path(), surfR_fname = midthicknessR_path())


  }

  #print(tstat_fname)
  # write images IMPORTANT NEEDS X SESSION

  if (save_png){

    # img save fname
    tstat_fname = paste0(dhcp_path(), "/partial_products/individual_t-test/sub-",
                         subid, "/", template_name, "/sub-", subid, "_ses-", sesid, "_desc_", template_name, "_tstat")
    activations_fname = paste0(dhcp_path(), "/partial_products/individual_t-test/sub-",
                               subid, "/", template_name, "/sub-", subid, "_ses-", sesid, "_desc_", template_name, "_", method_p, "-", alpha)



    if(!file.exists(paste0(tstat_fname, "_1.png"))) view_cifti_surface(tstats_xifti,surfL = env$surfL, surfR = env$surfR,
                                                                       idx = dhcpr::component_index.relevant(),
                                                                       fname = tstat_fname,
                                                                       title = paste0(template_name, "-t-stat sub- ", subid, "-", dhcpr::component_names()),
                                                                       zlim = c(-5, 5))#, color_mode = "diverging", colors = "PiYG", zlim = c(0,1,2))#, mwallL = env$mwallL, mwallR = env$mwallR)
    view_cifti_surface(activations$active, surfL = env$surfL, surfR = env$surfR,
                       idx = dhcpr::component_index.relevant(), fname = activations_fname,
                       title = paste0(template_name, "-", method_p, " < ", alpha, " sub- ", subid, "-", dhcpr::component_names()))#, zlim = c(-5, 5))#, color_mode = "diverging", colors = "PiYG", zlim = c(0,1,2))#, mwallL = env$mwallL, mwallR = env$mwallR)
  }


  return(activations)
}


#' individualSNR
#'
#' get individual snr as a matrix or cifti
#'
#'
#' @import ciftiTools
#'
#' @param bold_fname Filename of cifti timeseries for a single individual
#' @param output_fname Filename is null if matrix, cifti output.
#' @param use_symatlas To build output cifti
#' @param as_matrix Return matrix with SNR (V()x1) vector
#'
#'
#' @return matrix if set, void otherwise.
#'
#' @export




individualSNR = function (bold_fname, output_fname = NULL, use_symatlas = TRUE, as_matrix = TRUE){

  # define helper function rowVar
  rowVar <- function(x, ...) {
    var = rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
    #stopifnot(!is.na(var))
  }

  # parse input and catch errors ####
  if (missing(bold_fname)){
    stop("A cifti timeseries filename must be provided as input.")
  }


  if (use_symatlas){

    # get symatlas surfaces
    if (!exists("env")){
      env = environment()
      default_surfaces(env)
      # legacy surface asign
      Lsurf = env$surfL
      Rsurf = env$surfR
    }

    bold = read_cifti(bold_fname, surfL_fname = midthicknessL_path(), surfR_fname = midthicknessR_path())

    vertices = V() # Get number of vertices
    hemi_v = vertices / 2 # Get number of vertices in each hemishpere

  } else {

    bold = read_cifti(bold_fname) # Read in cifti file

  }

  bold_matrix = as.matrix(bold) # Take matrix of both hemispheres

  bold_timeavg = rowMeans(bold_matrix) # Average over time of the bold timeseries
  #check for NAs
  stopifnot(!any(is.na(bold_timeavg)))

  bold_timeVar = rowVar(bold_matrix) # Variance over time of the bold timeseries

  snr_matrix = bold_timeavg / bold_timeVar # Calculate snr according to Yeo et al. 2011
  # replace NAs for zeros
  snr_matrix[bold_timeVar == 0] = 0
  stopifnot(!any(is.na(snr_matrix)))

  # Get medialwall from input

  mwallL = env$mwallL #bold$meta$cortex$medial_wall_mask$left
  mwallR = env$mwallR #bold$meta$cortex$medial_wall_mask$right

  # Build SNR map cifti
  if (!is.null(output_fname)){
    write_cifti(snr_cifti, output_fname)
  }

  # Return matrix or cifti object depending on the parameter as_matrix

  if (as_matrix == FALSE){
    snr_cifti = as.xifti(cortexL = snr_matrix[1:hemi_v], cortexR = snr_matrix[-(1:hemi_v)],
                         cortexL_mwall = mwallL, cortexR_mwall = mwallR , surfL=Lsurf, surfR=Rsurf)
    return(snr_cifti)
  } else {
    return(snr_matrix)
  }
}
