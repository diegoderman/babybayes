#' List of Subjects
#'
#' Get subid and sesid of a set of subjects, given a csv. TODO: Add filtering by conditions?
#'
#' @description This function takes a table of subjects in csv format and returns a list of (subid, sesid) tuples(lists)
#'
#' @param table: input table
#'
#' @return list of lists with the duple (subid, sesid)
#'
#'
#' @keywords internal
#'
#' @export
dhcp_list = function(table_fname){

  # read in csv file
  table = read.csv(table_fname)

  # number of individual sessions (=rows)
  n = nrow(table)

  # initialize list of pairs
  session_list = list()

  # format session id as character
  options("scipen"=100, "digits"=4)
  table$session_id = as.character(table$session_id)

  # extract fromated vectors from the table (for expliciting only, could be ommited)
  subid = as.character(table$participant_id)
  sesid = table$session_id

  # make a list of 'n' lists, each containing a subid and the corresponding sesid

  for (i in 1:n){
    session_list[[i]] = c(subid[i], sesid[i])
  }

  return(session_list)
}


