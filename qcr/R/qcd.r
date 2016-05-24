#-----------------------------------------------------------------------------#
#                                                                             #
#                     QUALITY CONTROL DATA IN R                               #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores Sánchez                                       #
#              Student Master of Statistical Techniques                       #
#              University of The Coruña, SPAIN                                #
#              mflores@outlook.com                                            #
#                                                                             #
#-----------------------------------------------------------------------------#

#
#  Main function to create a 'qcd' object
#



##' Quality Control Data
##' 
##' Create an object of class 'qcd' to perform statistical quality control.
##' This object may then be used to plot Shewhart charts, Multivariate Control Charts,
##' and more.
##' 
##' 
##' @aliases qcd 
##' @param data a matrix or data-frame where it should contain data, index sample and,
##' optionally, covariate(s).
##' @param var.index a scalar with the column number corresponding the observed data for
##' the variable (the variable quality).  Alternativelly can be a string with the
##' name of the quality variable.
##' @param sample.index a scalar with the column number corresponding the index each
##' group (sample).
##' @param covar.index  optional. A scalar or numeric vector with the column number(s)
##' corresponding to the covariate(s). Alternativelly can be a character vector with
##' the names of the covariates.
##' @param covar.names  optional. A string or vector of strings with names for the
##' covariate columns.  Only valid if there is more than one column of data. By
##' default, takes the names from the original object.
##' @param data.name  a string specifying the name of the variable which appears on the
##' plots. If not provided is taken from the object given as data.
##' @param type.data  a string specifying el type de data.
##' @param sizes  optional. A value or a vector of values specifying the sample sizes
##' associated with each group. For continuous data the sample sizes are obtained counting the non-\code{NA} elements of
##' the sample.index vector. For attribute
##' variable the argument sizes is required.
##' @export

# data <- orangejuice[trial,]
# sizes <- sizes[trial]
# type.data <- "atributte"
# var.index=2
# sample.index=1

qcd <- function(data, var.index  =  1, sample.index  =  2,
                covar.index = NULL, covar.names = NULL,
                data.name = NULL,
                type.data  =  c("continuous", "atributte", "dependence"),
                sizes = NULL) 
  #.........................................................................
{
  
  if (!is.matrix(data) & !is.data.frame(data))
    stop("object must be a matrix or data.frame")
  
  type.data <- match.arg(type.data)
  
  
  if (type.data  ==  "atributte")
    if (is.null(sizes)) 
      stop("sample sizes must be given for a attribute variable")
  
  if (is.null(sizes)) sizes <- table(data[ ,sample.index])
  
  if (!is.null(covar.index)){
    result <- data[, c(var.index, sample.index, covar.index)]
  } else {
    result <- data[, c(var.index, sample.index)]  
  }
  
  result$sizes <- sizes
  
  
  if (!is.null(covar.index)) {
    if (!is.null(covar.names))
      names(result) <- c("x", "sample", covar.names,"sizes") 
    else 
      names(result)[1:2] <- c("x", "sample")
    names(result)[length(result)] <- c("sizes")
  } else 
    names(result) <- c("x", "sample","sizes")
  
  
  if (is.null(data.name))
    data.name <- deparse(substitute(data))
  
  attr(result, "data.name") <- data.name
  attr(result, "type.data") <- type.data
  
  oldClass(result) <- c("qcd", "data.frame")
  
  return(result)
} # qcd
#.........................................................................


#-------------------------------------------------------------------------
# state.control
#-------------------------------------------------------------------------
##' This function is used for the state of the process is under control. 
##' Were removed from the sample observations that violate the rules of a process under control
##' @aliases state.control
##' @param x  Object qcs (Quality Control Statitical)
##' @export
##' @examples
##' 
##' ##
##' ##  Continuous data 
##' ##
##'library(qcr)
##'data(pistonrings)
##'str(pistonrings)
##'pistonrings.qcd<-qcd(pistonrings)
##'
##'class(pistonrings.qcd)
##'
##'res.qcs <- qcs.xbar(pistonrings.qcd)
##'plot(res.qcs,title="Control Chart Xbar for pistonrings I")
##'summary(res.qcs)  
##'
##'res.qcd <- state.control(res.qcs)
##'res.qcs <- qcs.xbar(res.qcd)
##'plot(res.qcs,title="Control Chart Xbar for pistonrings II")
##'summary(res.qcs)  
##'
##'res.qcd <- state.control(res.qcs)
##'res.qcs <- qcs.xbar(res.qcd)
##'plot(res.qcs,title="Control Chart Xbar for pistonrings III")
##'summary(res.qcs)  
##'
state.control <- function(x) 
  #.........................................................................  
{
  if (length(x$violations[[1]])>0 || length(x$violations[[2]]>0)){
    ii<-which(is.na(match(x$qcd$sample,unlist(x$violations))))  
    result<-x$qcd[ii,]  
    result<-droplevels(result)
  } else {
    cat("The process is under control")
  }
  
  oldClass(result) <- c("qcd", "data.frame")
  
  invisible(result)
  
} #sate.control
#.........................................................................