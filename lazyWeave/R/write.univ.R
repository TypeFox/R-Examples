#' @name write.univ
#' @export write.univ
#' 
#' @title Write Univariate Table to a File
#' @description Write the output of \code{univ} to a file
#' 
#' @param x An object of type \code{univ}.
#' @param round Number of significant digits to be printed.
#' @param Factor Determines if the Factor (variable) description is printed.
#' @param Group Determines if the Group is printed
#' @param N Determines if the number of non missing values is printed
#' @param Missing Determines if the number of missing values is printed
#' @param Mean Determines if the mean is printed
#' @param SD Determines if the standard deviation is printed
#' @param LCL Determines if the lower confidence limit is printed
#' @param UCL Determines if the upper confidence limit is printed
#' @param Min Determines if the minimum value is printed
#' @param P25 Determines if the 25th percentile is printed
#' @param Median Determines if the median value is printed
#' @param P75 Determines if the 75th percentile is printed
#' @param Max Determines if the maximum value is printed
#' @param CV Determines if the coefficient of variation is printed
#' @param Pval Determines if the p-value is printed
#' @param pvalFormat Character string passed to \code{pvalString} and determines
#'   the pvalue style to be printed.
#' @param pvalArgs A list of additional arguments to be passed to \code{pvalString}
#' @param cat Logical. Determines if the output is returned as a character string
#'   or returned via the \code{cat} function (printed to console).  The default
#'   value is set by \code{options()$lazyWeave_cat}.  This argument allows for
#'   selective override of the default.
#' @param ... additional arguments to be passed to \code{lazy.matrix}
#' 
#' @examples
#' #output will be written to the working directory
#' getwd()
#' 
#' #write.univ function must be written to either a LaTeX
#' #or HTML file.  HTML format is through the lazyHTML package.
#' options(lazyReportFormat="html")
#' 
#' #Delivery dataset from CCFmisc library
#' data(Delivery)
#' 
#' #label the variables that will be used
#' Hmisc::label(Delivery$maternal.age) <- "Maternal Age"
#' Hmisc::label(Delivery$ga.weeks) <- "Gestation weeks"
#' Hmisc::label(Delivery$wt.gram) <- "Weight (g)"
#' 
#' #summaries of the continuous variables
#' #maternal.age, ga.weeks and wt.gram in the 
#' #Delivery dataset.
#' deliv.uni <- univ(Delivery,
#'                   vars=c("maternal.age", "ga.weeks", "wt.gram")
#' )
#' 
#' #summaries of continuous variables
#' #by level of delivery.type
#' delivBy.uni <- univ(Delivery,
#'                     vars=c("maternal.age", "ga.weeks", "wt.gram"),
#'                     byVar="delivery.type"
#' )
#' 
#' #to write univ based table to an HTML file enclose the
#' #write.univ() in the html_write function as below.
#' #see documentation for other options.
#' 
#' #To print byVariable group names in the table, 
#' #set the Group=T flag in the write.univ() function.
#' 
#' \dontrun{
#' lazy.write(
#'     lazy.file.start(),
#'     write.univ(deliv.uni),
#'     write.univ(delivBy.uni, Group=TRUE),
#'     lazy.file.end(),
#'     OutFile="ExampleFile.html"
#'   )
#'   
#'   unlink("ExampleFile.html")
#' }
#' 


write.univ <- function(x, round=1, 
    Factor=TRUE, Group=FALSE, N=TRUE, Missing=FALSE,
    Mean=TRUE, SD=TRUE, LCL=FALSE, UCL=FALSE, Min=TRUE,
    P25=TRUE, Median=TRUE, P75=TRUE, Max=TRUE, CV=FALSE, Pval=FALSE, 
    pvalFormat="default", pvalArgs=list(),
    cat=getOption("lazyWeave_cat"), ...){

  reportFormat <- getOption("lazyReportFormat")
  
  x <- x[, c(Factor, Group, N, Missing, Mean, SD, LCL, UCL,
             Min, P25, Median, P75, Max, CV, Pval), drop=FALSE]
  
  x$PVAL <- do.call("pvalString", 
                    c(list(p=x$PVAL, 
                           format=pvalFormat), 
                      pvalArgs))

  num <- sapply(x, is.numeric)
  x[, num] <- lapply(x[, num], round, round)
  x[is.na(x)] <- ""
  rownames(x) <- NULL

  if (cat) cat(lazy.matrix(x, ...))
  else return(lazy.matrix(x, ...))
}

