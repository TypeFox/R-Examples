######################################################################################
# print.TDMclassifier  and print.tdmClass
#
#'  Print an overview for a \code{\link{TDMclassifier}} or \code{\link[=tdmClassify]{tdmClass}} object.
#'
#'   @method print TDMclassifier
#'   @param x  an object of class \code{\link[=tdmClassify]{tdmClass}}, as returned from a prior call to \code{\link{tdmClassify}}, \cr
#'          or  an object of class \code{\link{TDMclassifier}}, as returned from a prior call to \code{\link{tdmClassifyLoop}}.
#'   @param ... e.g. 'type'    which information to print:
#'      \describe{
#'      \item{\code{"overview"}}{ (default) relative gain on training/test set, number of records, see \code{\link{tdmClassifySummary}}}
#'      \item{\code{"cm.train"}}{ confusion matrix on train set}
#'      \item{\code{"cm.vali"}}{ confusion matrix on test set}
#'      \item{\code{"?"}}{ help on this method}
#'      }
#'
#' @seealso   \code{\link{tdmClassify}}, \code{\link{tdmClassifySummary}}, \code{\link{TDMclassifier}}
#' @author Wolfgang Konen, FHK, 2011 - 2013
#' @export
print.TDMclassifier <- function(x,...) {
  internalPrintC <- function(result,type) {
    opts = result$lastRes$opts;
    opts$VERBOSE = 2;
    z <- switch(type
      , "overview"= {cat("Filename of task:",opts$filename);
              tdmClassifySummary(result,opts);
              cat("\nUse > print(result,type=\"?\") and > result$lastRes   for more info on TDMclassifier object result."); 
              cat("\nUse > names(result) and > result[]   to see all names and all contents (may be long) of object result."); 
              1;  # a value for z
              }
      , "cm.train"= show.cm.train(result)
      , "cm.vali"= show.cm.vali(result)
      , "?"={cat("Help for print(<TDMclassifier>,type=t). Possible values for 't' are:\n"
               ,"\"overview\": [def.] see help(tdmClassifySummary)\n"
               ,"\"cm.train\": confusion matrix on training data\n"
               ,"\"cm.vali\" : confusion matrix on validation data\n"
               ,"\"?\" : display this help message"
               ,"\nThe commands > result or > print(result) invoke the default print(result,type=\"overview\") for TDMclassifier object result"
               ); 1;}     # a value for z
      , "invalid type"
      );
    if (z[1]=="invalid type") warning("Invalid type = ",type,". Allowed types are: overview, cm.train, cm.vali, ?.");
    cat("\n");
  }

  vaarg <- list(...)
  #alternatively:
  # vavalues <- c(...)

  if (is.null(vaarg$type)) vaarg$type="overview";
  internalPrintC(x,vaarg$type);
}
######################################################################################
# print.tdmClass
# Print an overview for a \code{\link[=tdmClassify]{tdmClass}} object.
#'
#' @rdname print.TDMclassifier
#' @method print tdmClass
#' @export
print.tdmClass <- function(x,...) {
  internalPrintC <- function(res,type) {
    opts = res$opts;
    opts$VERBOSE = 2;
    z <- switch(type
      , "overview"= {show.tdmClass(res);
              cat("\nUse > print(res,type=\"?\")   for more info on tdmClass object res.");
              cat("\nUse > names(res) and > res[]   to see all names and all contents (may be long) of object res.\n"); 
              1;    # a value for z
              }
      , "cm.train"= show.cm.train(res)
      , "cm.vali"= show.cm.vali(res)
      , "?"={cat("Help for print(<tdmClass>,type=t). Possible values for 't' are:\n"
               ,"\"overview\": [def.] info on model and datasets in tdmClass object\n"
               ,"\"cm.train\": confusion matrix on training data\n"
               ,"\"cm.vali\" : confusion matrix on validation data\n"
               ,"\"?\" : display this help message"
               ,"\n\nThe commands > res or > print(res) invoke the default print(res,type=\"overview\") for tdmClass object res\n"
               ); 1;   # a value for z
            }
      , "invalid type"
      );
    if (z[1]=="invalid type") warning("Invalid type = ",type,". Allowed types are: overview, cm.train, cm.vali, ?.");
  }

  vaarg <- list(...)
  if (is.null(vaarg$type)) vaarg$type="overview";
  internalPrintC(x,vaarg$type);
}

show.tdmClass <- function(res) {
    opts = res$opts;
    cat("Model (for last response variable) of tdmClass object:")
    print(res$lastModel);
    cat("\nConfusion matrix on validation set\n");
    print(res$lastCmVali$mat);
    cat(sprintf("\nDatasets (# rows) of tdmClass object:\n  d_train (%d), d_test (%d), d_dis (%d)\n",
                      nrow(res$d_train),nrow(res$d_test),nrow(res$d_dis)));
}


# this helper fct works for both classes tdmClass and TDMClassifier in object result:
show.cm.train <- function(result) {
  cls <- class(result)[1];
  z <- switch(cls
    , "TDMclassifier"= {show.cm.train(result$lastRes); 1; }
    , "tdmClass"={
        # note that in here the object 'result' is of class ***tdmClass*** !
        opts = result$opts;            
        opts$VERBOSE = 2;
        cm.train <- result$lastCmTrain;
        cat1(opts,"Training cases (",nrow(result$d_train)
            ,ifelse(opts$TST.kind=="cv",", last fold )",")")
            ,sprintf("on last response variable '%s' :", tail(row.names(result$allEVAL),1))
            ,"\n");
        print1(opts,cm.train$mat)                     # confusion matrix on training set
        print1(opts,cm.train$gain.vector)
        cat1(opts,sprintf("total gain: %7.1f (is %7.3f%% of max. gain = %7.1f)\n",
                          cm.train$gain,cm.train$gain/cm.train$gainmax*100,cm.train$gainmax));
        cat1(opts,sprintf("Relative gain (rgain.type='%s') is %7.2f%%\n",opts$rgain.type,cm.train$rgain));
        1;    # a value for z
       }
    , "invalid class"
    );
  if (z[1]=="invalid class") warning("Invalid class = ",cls,". Allowed classes are: TDMclassifier, tdmClass.");
  result;
}

# this helper fct works for both classes tdmClass and TDMClassifier in object result:
show.cm.vali <- function(result) {
  cls <- class(result)[1];
  z <- switch(cls
    , "TDMclassifier"= {show.cm.vali(result$lastRes); 1; }
    , "tdmClass"={
        # note that in here the object 'result' is of class ***tdmClass*** !
        opts = Opts(result);
        opts$VERBOSE = 2;
        cm.vali <- result$lastCmVali;
        n.vali <- nrow(result$d_test);
        if (n.vali>0) {
          cat1(opts,"Vali cases (",n.vali
              ,ifelse(opts$TST.kind=="cv",", last fold )",")")
              ,sprintf("on last response variable '%s' :", tail(row.names(result$allEVAL),1))
              ,"\n");
          print1(opts,cm.vali$mat)                      # confusion matrix on test set
          print1(opts,cm.vali$gain.vector)
          cat1(opts,sprintf("total gain : %7.1f (is %7.3f%% of max. gain = %7.1f)\n",
                            cm.vali$gain,cm.vali$gain/cm.vali$gainmax*100,cm.vali$gainmax));
          cat1(opts,sprintf("Relative gain (rgain.type='%s') is %7.2f%%\n",opts$rgain.type,cm.vali$rgain));
        } else {
          cat1(opts,"There are no test cases in data set!\n");
        }
        1;    # a value for z
       }
    , "invalid class"
    );
  if (z[1]=="invalid class") warning("Invalid class = ",cls,". Allowed classes are: TDMclassifier, tdmClass.");
  result;
}

