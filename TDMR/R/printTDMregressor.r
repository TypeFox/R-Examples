######################################################################################
# print.TDMregressor     and print.tdmRegre
#
#'   Print an overview for a \code{\link{TDMregressor}} or \code{\link[=tdmRegress]{tdmRegre}} object.
#'
#'   @method print TDMregressor
#'   @param x  an object of class \code{\link[=tdmRegress]{tdmRegre}}, as returned from a prior call to \code{\link{tdmRegress}}, \cr
#'          or  an object of class \code{\link{TDMregressor}}, as returned from a prior call to \code{\link{tdmRegressLoop}}.
#'   @param ... e.g. 'type'    which information to print:
#'      \describe{
#'      \item{\code{"overview"}}{ (def.) RMAE on training/test set, number of records, see \code{\link{tdmRegressSummary}}}
#'      \item{\code{"..."}}{ ... other choices, TODO ...}
#'      \item{\code{"?"}}{ help on this method}
#'      }
#' @seealso   \code{\link{tdmRegress}}, \code{\link{tdmRegressSummary}}, \code{\link{TDMregressor}}
#' @export
######################################################################################
print.TDMregressor <- function(x,...) {
  internalPrintR <- function(result,type) {
    opts = result$lastRes$opts;
    opts$VERBOSE = 2;
    z <- switch(type
      , "overview"= { tdmRegressSummary(result,opts);
              cat("\nUse > print(result,type=\"?\") and > result$lastRes   for more info on TDMregressor object result."); 
              cat("\nUse > names(result) and > result[]   to see all names and all contents (may be long) of object result."); 
              1;  # a value for z
              }
      , "?"={cat("Help for print(<TDMregressor>,type=t). Possible values for 't' are:\n"
               ,"\"overview\": see tdmRegressSummary\n"
               ,"\"?\" : display this help message\n"
               ); 1;}   # a value for z
      , "invalid type"
      );
    if (z[1]=="invalid type") warning("Invalid type = ",type,". Allowed types are: overview, ?.");
    cat("\n");
  }

  vaarg <- list(...)
  #alternative: vavalues <- c(...)

  if (is.null(vaarg$type)) vaarg$type="overview";
  internalPrintR(x,vaarg$type);
}

######################################################################################
# print.tdmRegre
# Print an overview for a \code{\link[=tdmRegress]{tdmRegre}} object.
#'
#' @rdname print.TDMregressor
#' @method print tdmRegre
#' @export
######################################################################################
print.tdmRegre <- function(x,...) {
  internalPrintC <- function(res,type) {
    opts = res$opts;
    opts$VERBOSE = 2;
    z <- switch(type
      , "overview"= {show.tdmRegre(res);
              cat("\nUse < print(res,type=\"?\")   for more info on tdmRegre object res.\n");
              1;    # a value for z
              }
      , "?"={cat("Help for print(<tdmRegre>,type=t). Possible values for 't' are:\n"
               ,"\"overview\": [def.] info on model and datasets in tdmRegre object\n"
               ,"\"?\" : display this help message"
               ,"\n\nThe commands > res or > print(res) invoke the default print(res,type=\"overview\") for tdmRegre object res\n"
               ); 1;   # a value for z
            }
      , "invalid type"
      );
    if (z[1]=="invalid type") warning("Invalid type = ",type,". Allowed types are: overview, ?.");
  }

  vaarg <- list(...)
  if (is.null(vaarg$type)) vaarg$type="overview";
  internalPrintC(x,vaarg$type);
}

show.tdmRegre <- function(res) {
    opts = res$opts;
    cat("Model (for last response variable) of tdmRegre object:")
    print(res$lastModel);
    cat(sprintf("\nDatasets (# rows) of tdmRegre object:\n  d_train (%d), d_test (%d)\n",
                      nrow(res$d_train),nrow(res$d_test)));
}

