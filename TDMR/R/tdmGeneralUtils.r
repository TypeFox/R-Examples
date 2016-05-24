######################################################################################
######################################################################################
# Package Description for Roxygene:
#' Tuned Data Mining in R
#'
#' \tabular{ll}{
#' Package: \tab TDMR\cr
#' Type: \tab Package\cr
#' Version: \tab 1.3\cr
#' Date: \tab 30.08.2014\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' TDMR is a package for tuned data mining (predictive analytics, i.e. \bold{classification} and \bold{regression}). Its main features are: \cr
#' 1) A variety of tuners, with special emphasis on \link{SPOT} (a well-known R package for parameter tuning), but also CMA-ES 
#'    (package \code{\link[rCMA]{rCMA}}) and other tuning algorithms. \cr
#' 2) Tuning of preprocessing (feature generation) parameters and model building parameters simultaneously.  \cr
#' 3) Support for multiple tuning experiments (different settings, repetitions with different resamplings, ...).  \cr
#' 4) Easy parallelization of those experiments with the help of R packages \code{snow} and \code{\link{parallel}}.  \cr
#' 5) Extensibility: New tuning parameters, new feature preprocessing tools, model builders and even new tuners can be added easily.
#' 
#' The main entry point functions are \code{\link{tdmClassifyLoop}}, \code{\link{tdmRegressLoop}} and \code{\link{tdmBigLoop}}. 
#' See \code{\link{tdmOptsDefaultsSet}} and \code{\link{tdmDefaultsFill}} for an overview of adjustable TDMR-parameters.
#'                                                
#' @name TDMR-package
#' @aliases TDMR
#' @docType package
#' @title Tuned Data Mining in R
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de}), Patrick Koch
#' @references \url{http://gociop.de/research-projects/tuned-data-mining/}
#' @keywords package tuning data mining machine learning
#' @import SPOT 
#' @import twiddler 
#' @import testit
#' @import tcltk
#### @import adabag
#### no longer needed:  @import e1071     (because e1071 moved from "Depends" to "Suggests" in DESCRIPTION)


#End of Package Description
NA #NULL, ends description without hiding first function
######################################################################################
######################################################################################



######################################################################################
#
# GENERAL UTILITY FUNCTIONS
#
######################################################################################
######################################################################################

######################################################################################
#' Bind a column to a data frame. 
#'
#' Bind the column with name \code{response.predict} and 
#' contents \code{vec} as last column to data frame \code{d}
#' @param d data frame
#' @param response.predict name of new column 
#' @param vec the contents for the last column bound to data frame \code{d}
#' @return data frame \code{d} with column added
#' @export
######################################################################################
tdmBindResponse <- function(d,response.predict,vec)
{
    if (is.na(match(response.predict,names(d)))) {
      # bind column response.predict as last column to data frame d
      eval(parse(text=paste("d <- cbind(d, ",response.predict,"=vec)")))
    } else {
      # replace contents of existing column response.predict
      d[,response.predict] <- vec;
    }

    return(d)
}
#
# this older version does the same, but requires three copy-replacements of
# data frame d (instead of one):
tdmBindResponse_OLD <- function(d,response.predict,vec)
{
    # drop column response.predict if there, do nothing if not there
    d <- d[,setdiff(names(d),response.predict)]
    # bind column response.predict as last column to data frame d
    d <- cbind(d, prediction=vec)
    names(d)[ncol(d)] <- response.predict

    return(d)
}

######################################################################################
# printout functions for different verbosity levels
######################################################################################
#' Output objects to \code{cat} if \code{opts$VERBOSE>=1}.
#'
#' @param opts from which we need the element VERBOSE
#' @param ...  objects
#' @return None
#' @seealso   \code{\link{cat}}
#' @export
#' @keywords internal
cat1 <- function(opts, ...) {  if (opts$VERBOSE>=1) cat(...); }
######################################################################################
#' Output objects to \code{cat} if \code{opts$VERBOSE>=2}.
#'
#' @param opts from which we need the element VERBOSE
#' @param ...  objects
#' @return None
#' @seealso   \code{\link{cat}}
#' @export
#' @keywords internal
cat2 <- function(opts, ...) {  if (opts$VERBOSE>=2) cat(...); }

######################################################################################
#' Print objects using \code{print} if \code{opts$VERBOSE>=1}.
#'
#' @param opts from which we need the element VERBOSE
#' @param ...  objects
#' @return None
#' @seealso   \code{\link{print}}
#' @export
#' @keywords internal
print1 <- function(opts, ...) {  if (opts$VERBOSE>=1) print(...); }

######################################################################################
#' Print objects using \code{print} if \code{opts$VERBOSE>=2}.
#'
#' @param opts from which we need the element VERBOSE
#' @param ...  objects
#' @return None
#' @seealso   \code{\link{print}}
#' @export
#' @keywords internal
print2 <- function(opts, ...) {  if (opts$VERBOSE>=2) print(...); }


######################################################################################
#
# UTILITY FUNCTIONS FOR TDMenvir, TDMclassifier AND TDMregressor OBJECTS
#
######################################################################################
######################################################################################

######################################################################################
#' Make a prediction using the last model.
#'
#' Make a prediction with objects of class \code{\link{TDMenvir}}, \code{\link{TDMclassifier}},
#' \code{\link{TDMregressor}}. The prediction is based on the (last) model trained during 
#' \code{\link{unbiasedRun}}.
#'
#'   @param object  an object of class \code{\link{TDMenvir}}, \code{\link{TDMclassifier}},
#'      \code{\link{TDMregressor}} containing in element \code{lastModel} the relevant model.
#'   @param ... arguments passed on to the model's predict function. Usually the first argument of \code{...} should be 
#'        \code{newdata}, a data frame for which new predictions are desired.
#'   @return a vector with length \code{nrow(newdata)} containing the new predictions.
#' @method predict TDMenvir
#' @examples
#'    \dontrun{
#'        ## This example requires that demo04cpu.r is executed first (it will write demo04cpu.RData)
#'        path <- paste(find.package("TDMR"), "demo01cpu/",sep="/");
#'        tdm <- list(  filenameEnvT="demo04cpu.RData" );   # file with environment envT 
#'        load(paste(path,tdm$filenameEnvT,sep="/"));
#'                   
#'        # take only the first 15 records:
#'        newdata=read.csv2(file=paste(path,"data/cpu.csv", sep=""), dec=".")[1:15,];     
#'        z=predict(envT,newdata);
#'        print(z);
#'    }
#' @export
predict.TDMenvir <- function(object,...) {
  if (is.null(object$result)) stop("TDMenvir has no element 'result'");
  if (is.null(object$result$lastRes)) stop("TDMenvir has no element 'result$lastRes'");
  if (is.null(object$result$lastRes$lastModel)) stop("TDMenvir has no element 'result$lastRes$lastModel'. Consider to set tdm$U.saveModel=TRUE.");
  predict(object$result$lastRes$lastModel,...);
}
######################################################################################
#' @rdname predict.TDMenvir
#' @method predict TDMclassifier
#' @export
# @keywords internal
predict.TDMclassifier <- function(object,...) {
  predict(object$lastRes$lastModel,...);
}
######################################################################################
#' @rdname predict.TDMenvir
#' @method predict TDMregressor
#' @export
predict.TDMregressor <- function(object,...) {
  predict(object$lastRes$lastModel,...);
}


######################################################################################
#' Return the list 'opts'.
#'
#' Returns the list \code{opts} from objects of class \code{\link{TDMenvir}}, \code{\link{TDMclassifier}},
#'                  \code{\link{TDMregressor}}, \code{\link[=tdmClassify]{tdmClass}} or \code{\link[=tdmRegress]{tdmRegre}}.
#'   @param x  an object of class \code{\link{TDMclassifier}}, \code{\link[=tdmClassify]{tdmClass}}, 
#'                                \code{\link{TDMregressor}} or \code{\link[=tdmRegress]{tdmRegre}}.
#'   @param ... -- currently not used -- 
#'   @return the list \code{opts} with DM-specific settings contained in the specified object
#' @export
Opts <- function(x, ...)  UseMethod("Opts");

######################################################################################
#' @rdname Opts
#' @method Opts TDMenvir
#' @export
# @keywords internal
Opts.TDMenvir  <- function(x, ...) {
  if (is.null(x$result)) stop("TDMenvir has no element 'result'");
  if (is.null(x$result$lastRes)) stop("TDMenvir has no element 'result$lastRes'");
  x$result$lastRes$opts;
}
######################################################################################
#' @rdname Opts
#' @method Opts TDMclassifier
#' @export
# @keywords internal
Opts.TDMclassifier  <- function(x, ...) x$lastRes$opts;
#' @rdname Opts
#' @method Opts TDMregressor
#' @export
# @keywords internal
Opts.TDMregressor  <- function(x, ...) x$lastRes$opts;
#' @rdname Opts
#' @method Opts tdmClass
#' @export
# @keywords internal
Opts.tdmClass <- function(x,...) x$opts;
#' @rdname Opts
#' @method Opts tdmRegre
#' @export
# @keywords internal
Opts.tdmRegre <- function(x,...) x$opts;
#' @rdname Opts
#' @method Opts default
#' @export
# @keywords internal
Opts.default <- function(x, ...)  cat("This is Opts.default\n");

######################################################################################

#summary.TDMclassifier <- function(result,...) {
#  cat("This is the TDMclassifier summary\n");
#}


RGainTST <- function(result) {
  if (!inherits(result, "TDMclassifier"))
        stop("This function is permitted only for objects of class `TDMclassifier'");
  mean(result$R_vali);
}

RGainTRN <- function(result) {
  if (!inherits(result, "TDMclassifier"))
        stop("This function is permitted only for objects of class `TDMclassifier'");
  mean(result$R_train);
}

SRF <- function(result) {
  if (!inherits(result, "TDMclassifier"))
        stop("This function is permitted only for objects of class `TDMclassifier'");
  result$lastRes$SRF;
}

