#' Creates numeric variables
#'
#' R implementation of the SPSS \code{NUMERIC} function. Creates new numeric variables, which get appended at the end of the dataset. 
#' 
#' xpssNumeric creates new numeric variables, which get appended at the end of the dataset. The new variables are as long as the selected dataset. By default the new variables are blank and get filled with NA, otherwise every case for the selected variable get filled with the speficied value. 
#'
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @param varname atomic character or character vector with the name of the variables which should be created.
#' @param fill atomic numeric or atomic character values, which fill the variables. By default the value is NA for all new variables. It is possible to assign each new variable an own value to fill with.
#' @return Returns the input data extended by the new variables.
#' @author Andreas Wygrabek
#' @seealso \code{\link{xpssString}}
#' @examples \dontrun{
#' xpssNumeric(fromXPSS, varname = c("A"), fill = c(NA))
#' }
#' @export 


xpssNumeric <- function(x, varname = NULL, fill = NA){
  
  
  ####################################################################
  ####################################################################
  
  functiontype <- "DM"
  x <- applyMetaCheck(x)
  
  ####################################################################
  ####################################################################
  ####################################################################

  
    if(is.matrix(x)){x <- as.data.frame(x)}
    if (length(fill) == 1){
    attributes(x$porsche)
            for(i in 1:length(varname)){
                eval(parse(text = paste("x$",varname[i]," <- ", fill, sep = "")))
                eval(parse(text = paste("attr(x$",varname[i],",'varname') <- ",'varname[i]',sep="")))
                eval(parse(text = paste("attr(x$",varname[i],",'variable.label') <- ",'varname[i]',sep="")))
            } }else {
                    for(i in 1:length(varname)){
                        eval(parse(text = paste("x$",varname[i]," <- ", fill[i], sep = "")))
                        eval(parse(text = paste("attr(x$",varname[i],",'varname') <- ",'varname[i]',sep="")))
                        eval(parse(text = paste("attr(x$",varname[i],",'variable.label') <- ",'varname[i]',sep="")))
                    }
    }
  
  
  
  x <- applyAttributeDemerge(x)
  
    return(x)
}
