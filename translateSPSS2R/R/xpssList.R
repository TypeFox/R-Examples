#' Displays the content of variables.
#'
#' R implementation of the SPSS \code{LIST} Function
#'
#' LIST displays the content of selected variables. It is possible to display a sequence with the \code{cases} argument.\code{from} determine the begin of the sequence, \code{to} determine the end of the sequence. \code{by} determine the increment of the sequence.
#'
#' @usage xpssList(x, variables = colnames(x), cases = list(from = 1, to = nrow(x), by = 1)) 
#'
#' @param x a (non-empty) data.frame or input data of class \code{xpssFrame}. 
#' @param variables atomic character or character vector with the names of the variables.
#' @param cases list containing the arguments from, to, by. All parameters are atomic numerics. See Details for more.

#' @return A data.frame with case values for specified variables in the dataset. If cases and variables are not specified, List return the complete dataset. If cases are specified the output is a user-defined sequence.
#' @author Bastian Wiessner
#' @examples
#' data(fromXPSS)
#' 
#' xpssList(x=fromXPSS)
#' 
#' xpssList(x=fromXPSS, 
#'    variables = "V1")
#' 
#' xpssList(x=fromXPSS, 
#'    variables = c("V1","V2"))
#' 
#' 
#' xpssList(x=fromXPSS, 
#'    variables = span(fromXPSS,
#'                  from="V1",
#'                  to="V4"),
#'    cases =list(from=2,
#'                to=18,
#'                by=2))
#' 
#' @export

xpssList <- function(x,
                     variables = colnames(x),
                     cases = list(from = 1,
                                  to = nrow(x),
                                  by = 1)) 
{
 
  functiontype <- "SB"
  x <- applyMetaCheck(x)

  if(!is.numeric(cases$from) || !is.numeric(cases$to) || !is.numeric(cases$by))  {
    stop("the arguments from, to and by have to numeric")
  }
  if(cases$to > nrow(x)) {
    stop("argument to is bigger than the dataset")
  }
  
    pos <- seq(cases$from,cases$to,cases$by)
    erg <- data.frame(1:length(pos))
    if(length(variables) == 1){      
      erg <- as.data.frame(x[,variables][pos],stringsAsFactors=F)
      names(erg) <- names(x[which(colnames(x)%in%variables)])
      
    }else{
      for(i in 1:length(variables)){
        erg[[i]] <-   as.data.frame(x[,variables[i]][pos])
        names(erg[[i]]) <- colnames(x[i])
      }
    }
  return(erg) 
}



 
