#' Workbench Paths and Definitions
#'
#' Get workbench path depending on the host
#'
#' @description This function asks for the hosts and returns the appropriate path to wb_command binary
#'
#' @param par This is a test
#'
#' @return The path to workbench command
#'
#'
#' @keywords internal
#'
#' @export
#'

V = function(){
  return(57700)
}

wb_path = function(){
  home = Sys.getenv("HOME")

  ## getenv, select path

  if (substring(home, 1, 2) == "/N"){
    host = "carbonate"
  } else {
    host = "viper"
  }

  if (host == "viper"){
    baby_ICA = "~/baby_ICA"
    wbpath = "~/disco/Research/workbench/bin_linux64/wb_command"
  } else if (host == "carbonate"){
    baby_ICA = "/N/project/baby_ICA"
    wbpath = "~/dev-dhcp/workbench/bin_rh_linux64/wb_command"
  } else if (host == "job_manager"){
    baby_ICA = "/N/project/baby_ICA"
    wbpath = "~/dev-dhcp/workbench/bin_rh_linux64/wb_command"
  } else stop("Hostname not recognized.", call.=FALSE)

  return(wbpath)
}

#' check_tica
#'
#' check version of templateICAr library, if not desired, choose if stop
#'
#' @description This functions check version of templateICAr library. If compares false, it stops if stop_on_false = TRUE
#'
#' @import utils
#'
#' @param version char with version number
#' @param stop_on_false bool
#'
#' @return silent
#'
#'
#' @keywords internal
#'
#' @export
#'
check_tica = function(version = '0.2.2', stop_on_false = TRUE){
  if (packageVersion("templateICAr") != version){
    if (stop_on_false) stop(paste0("Error: templateICAr is version",  packageVersion("templateICAr"), "not ", version))
    else warning(paste0("Warning: templateICAr is version",  packageVersion("templateICAr"), "not ", version))
  }
}

#' DHCP path
#'
#' Get dhcp path depending on the host
#'
#' @description This function asks for the hosts and returns the appropriate path to wb_command binary
#'
#' @param par This is a test
#'
#' @return The path to workbench command
#'
#'
#' @keywords internal
#'
#' @export
#'
dhcp_path = function(){

  home = Sys.getenv("HOME")

  ## getenv, select path

  if (substring(home, 1, 2) == "/N"){
    host = "carbonate"
  } else {
    host = "viper"
  }

  if (host == "viper"){
    baby_ICA = "~/baby_ICA"
    wbpath = "~/disco/Research/workbench/bin_linux64/wb_command"
  } else if (host == "carbonate"){
    baby_ICA = "/N/project/baby_ICA"
    wbpath = "~/stICA/workbench/bin_linux64/wb_command"
  } else if (host == "job_manager"){
    baby_ICA = "~/dHCP_dataset"
    wbpath = "~/stICA/workbench/bin_linux64/wb_command"
  } else stop("Hostname not recognized.", call.=FALSE)

  return(baby_ICA)
}

midthicknessL_path = function(){
  return(paste0(dhcp_path(), "/surf_symatlas/week-40_hemi-left_space-dhcpSym_dens-32k_midthickness.surf.gii"))
}

inflatedL_path = function(){
  return(paste0(dhcp_path(), "/surf_symatlas/from-wm_inflated-L_dhcpSym_dens-32k_week-40.surf.gii"))
}

midthicknessR_path = function(){
  return(paste0(dhcp_path(), "/surf_symatlas/week-40_hemi-right_space-dhcpSym_dens-32k_midthickness.surf.gii"))
}

inflatedR_path = function(){
  return(paste0(dhcp_path(), "/surf_symatlas/from-wm_inflated-R_dhcpSym_dens-32k_week-40.surf.gii"))
}


#' Medial wall
#'
#' Get dhcp midthickness symatalas medialwall
#'
#' @description This function asks for the hosts and returns the appropriate path to wb_command binary
#'
#' @param par This is a test
#'
#' @return The path to workbench command
#'
#' @import ciftiTools
#'
#' @keywords internal
#'
#' @export
#'
mwallL = function(){
  requireNamespace("ciftiTools")
  mwall = read_cifti(paste0(dhcp_path(), "/groupICA/results/dhcp_s-24_c-20_k-4/melodic_IC.dscalar.nii"))$meta$cortex$medial_wall_mask$left
  return(mwall)
}

#' Medial wall
#'
#' Get dhcp midthickness symatalas medialwall
#'
#' @description This function asks for the hosts and returns the appropriate path to wb_command binary
#'
#' @param par This is a test
#'
#' @return The path to workbench command
#'
#' @import ciftiTools
#'
#' @keywords internal
#'
#' @export
#'
mwallR = function(){
  requireNamespace("ciftiTools")
  mwall = read_cifti(paste0(dhcp_path(), "/groupICA/results/dhcp_s-24_c-20_k-4/melodic_IC.dscalar.nii"))$meta$cortex$medial_wall_mask$right
  return(mwall)
}

#' Medial wall
#'
#' Get dhcp midthickness symatalas medialwall
#'
#' @description This function asks for the hosts and returns the appropriate path to wb_command binary
#'
#' @param par This is a test
#'
#' @return The path to workbench command
#'
#' @import ciftiTools
#'
#' @keywords internal
#'
#' @export
#'
default_mwall = function(env){
  env$mwallR = mwallR()
  env$mwallL = mwallL()
}

#' Midthickness surfaces
#'
#' Get dhcp midthickness symatalas surfaces
#'
#' @description This function asks for the hosts and returns the appropriate path to wb_command binary
#'
#' @param par This is a test
#'
#' @return The path to workbench command
#'
#' @import ciftiTools
#'
#' @keywords internal
#'
#' @export
#'
default_surfaces = function(env){
  env$surfR = read_surf(midthicknessR_path())
  env$surfL = read_surf(midthicknessL_path())
}

inflated_surfaces = function(env){
  env$surfR = read_surf(inflatedR_path())
  env$surfL = read_surf(inflatedL_path())
}

#' Definitions
#'
#' Get components names for the dhcp analysis (CHANGES WITH EACH run of gICA)
#'
#' @description This functions returns a vector with the names of the RSN components.
#' It should be modified according to the groupICA results of that specific run
#'
#' @param NULL
#'
#' @return vector with RSN names
#'
#'
#' @keywords internal
#'
#' @export
#'

#TODO: add ommit explanation and warning

component_names = function(ommit = TRUE, spaces = TRUE){

  if (ommit){
    if (spaces){
      components = c("Medial Motor",
                     "Lateral Motor",
                     "PCC",
                     "Posterior Parietal",
                     "Motor Association",
                     "Frontoparietal",
                     "Somatosensory",
                     "Visual Association R",
                     "Inferior Parietal",
                     "Auditory",
                     "Language",
                     "Visual",
                     "Visual Association L",
                     "Prefrontal")}
    else{
      components = c("Medial_Motor",
                     "Lateral_Motor",
                     "PCC",
                     "Posterior_Parietal",
                     "Motor_Association",
                     "Frontoparietal",
                     "Somatosensory",
                     "Visual_Association_R",
                     "Inferior_Parietal",
                     "Auditory",
                     "Language",
                     "Visual",
                     "Visual_Association_L",
                     "Prefrontal")}
  } else {
    if (spaces){
    components = c("Medial Motor",
                   "Lateral Motor",
                   "PCC",
                   "Posterior Parietal",
                   "Motor Association",
                   "Frontoparietal",
                   "Somatosensory",
                   "Visual Association R",
                   "Inferior Parietal",
                   "Auditory",
                   "Language",
                   "Visual",
                   "NA",
                   "Visual Association L",
                   "Prefrontal")
    } else {
      components = c("Medial_Motor",
                     "Lateral_Motor",
                     "PCC",
                     "Posterior_Parietal",
                     "Motor_Association",
                     "Frontoparietal",
                     "Somatosensory",
                     "Visual_Association_R",
                     "Inferior_Parietal",
                     "Auditory",
                     "Language",
                     "Visual",
                     "Nuisance",
                     "Visual_Association_L",
                     "Prefrontal")
    }
  }

  return(components)
}

#' relevant components index
#'
#' Get components indices for the dhcp analysis (CHANGES WITH EACH run of gICA)
#'
#' @description This functions returns a vector with the names of the RSN components.
#' It should be modified according to the groupICA results of that specific run
#'
#' @param NULL
#'
#' @return vector with RSN indices
#'
#'
#' @keywords internal
#'
#' @export
#'
component_index.relevant = function(){
  return(c(1:12, 14, 15))
}
