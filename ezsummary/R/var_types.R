#' Attach the variable type information with the dataset
#'
#' @description In order to analyze variables in the most appropriate way using this
#' \code{ezsummary} package, you'd better let the computer know what types of data
#' (quantitative or categorical) you are asking it to compute. This function will
#' attach a list of types you entered with the datasets so functions down the stream
#' line can read these information and analyze based on that. The information is stored
#' in the attributes of the dataset
#'
#' @param tbl A data.frame
#' @param types Character vector of length equal to the number of variables in the
#' dataset. Use "q" and "c" to denote quantitative and categorical variables.
#'
#' @export

var_types <- function(tbl, types){
  if(!is.data.frame(tbl))stop("Please supply a data.frame/data.table as the value of tbl")
  if(nchar(types) != length(attributes(tbl)$names))stop("The length of the string you entered doesn't match the number of variables")
  if(grepl("[^cq]", types)) stop('Unrecognizable character(s) detected!! Please review your input and use "q" and "c" to denote quantitative and categorical variables')
  attributes(tbl)$var_types <- unlist(strsplit(types, ""))
  return(tbl)
}

#' Automaticall assign var_types to the attributes
#'
#' @description If the user did not provide var_types, the function will preassume every variables to be quantitative
#' variables. If a variable's type is character but either the user or the automatic step said it's a quantitative
#' variable, the var_types attribute for that variable will be overwritten as categorical variable. At the same time,
#' a warning message will be printed on the screen.
#'
#' @param tbl The imported data.frame
#'
#' @export
auto_var_types <- function(tbl){
  col_class <- sapply(tbl, class)
  if(is.null(attributes(tbl)$var_types)){
    attributes(tbl)$var_types <- ifelse(col_class == "numeric", "q", "c")
    }
  if (length(attributes(tbl)$vars) != 0){
    attributes(tbl)$var_types[names(tbl) %in% attributes(tbl)$vars]<-"g"
  }
  return(tbl)
}
