#' @name factorFormula
#' @importFrom stringr str_extract_all
#' @importFrom stringr str_split_fixed
#' 
#' @title Convert Factor Levels in Formula to Numeric Values
#' @description When working in R, it is often more convenient to work in
#'   terms of the factor labels rather than the underlying numeric values.
#'   JAGS, however, requires that the numeric values be used.  
#'   \code{factorFormula} permits the user to define formulae to be passed
#'   to JAGS using R style coding, and having factor levels translated 
#'   to the underlying values as determined by the network structure.
#'   
#' @param form A formula object.
#' @param network A \code{HydeNetwork} object.
#' 
#' @details It is assumed that factor variables will be used in logical
#'   comparisons of the format \code{[variable_name] == '[factor_level]'} and
#'   only this pattern is recognized in the text search.  Single or 
#'   double quotes may be used around the level, and the spaces aroudn the
#'   \code{==} are optional.  
#'   
#'   While there are certainly limitations to this function that we have
#'   not yet found, we believe it covers the majority of cases in which
#'   it is useful. More complex cases that can't be handled by 
#'   \code{factorFormula} may require writing native JAGS code.
#'   
#' @author Jarrod Dalton and Benjamin Nutter
#' 
#' @examples 
#' \dontrun{
#' Net <- HydeNetwork(~ wells +
#'                     pe | wells +
#'                     d.dimer | pregnant*pe +
#'                     angio | pe +
#'                     treat | d.dimer*angio +
#'                     death | pe*treat,
#'                   data = PE)
#' factorFormula(form = payoff ~ (death == 'No') + (pe == 'Yes'),
#'               network = Net)
#' }
#' 
factorFormula <- function(form, network){
  form <- deparse(form)
  form <- trimws(form)
  form <- paste0(form, collapse = "")
  
  relabel <- extractFactors(form)
  
  relabel_mat <- isolateVariableFromLabel(relabel, network)
  
  new_label <- 
    mapply(getNumericLevel, 
           varname = relabel_mat[, 1], 
           label = relabel_mat[, 2],
           nodeType = relabel_mat[, 3],
           MoreArgs = list(network = network))
  
  if (any(vapply(new_label, length, numeric(1)) == 0))
  {
    noFactors <- unique(names(new_label[vapply(new_label, length, numeric(1)) == 0]))
    stop(paste0("The following nodes do not have factor levels defined ",
                "in the 'factorLevels' element of the HydeNetwork object: ", 
                paste0(noFactors, collapse = ", ")))
  }
  
  form <- rewriteFormula(relabel, new_label, form)
  
  as.formula(form)
}

extractFactors <- function(form){
  #* This pattern looks for any combination of numbers, letters, 
  #* periods or underscores, 
  #* followed by a space (or not), 
  #* followed by ==
  #* followed by a space (or not)
  #* followed by a quote (single or double)
  #* followed by any character string
  #* followed by a quote (single or double)
  #* It is intened to catch a [variable_name] == [factor_level]
  relabel <- 
    stringr::str_extract_all(string = form, 
        pattern = "[[:alpha:],[0-9],[.],[_]]+( |)[=][=]( |)('|\").*?('|\")")
  unlist(relabel)
}

isolateVariableFromLabel <- function(relabel, network){
  relabel_mat <- stringr::str_split_fixed(relabel, "[=][=]", 2)
  relabel_mat <- trimws(relabel_mat)
  relabel_mat <- gsub("('|\")", "", relabel_mat)
  cbind(relabel_mat,
        unlist(network$nodeType[relabel_mat[, 1]]))
}

getNumericLevel <- function(varname, label, nodeType, network){
  value <- which(network$factorLevels[[varname]] == label)
  if (nodeType == "dbern") value <- as.numeric(value) - 1
  as.character(value)
}

rewriteFormula <- function(relabel, new_label, form){
  for (i in seq_along(relabel)){
    new_label[i] <- sub("(?<=('|\")).*?(?=('|\"))", 
                        new_label[i], 
                        relabel[i], 
                        perl = TRUE)
    new_label[i] <- gsub("('|\")", "", new_label[i])
    form <- sub(relabel[i], new_label[i], form)
  }
  form
}

