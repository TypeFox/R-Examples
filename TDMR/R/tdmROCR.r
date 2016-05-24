######################################################################################
# tdmROCR
#
#'   Interactive plot of ROC, lift or other charts for a \code{\link{TDMclassifier}} object.
#'   See \code{\link{tdmROCR.TDMclassifier}} for details
#'
#'   @param x  return value from a prior call to \code{\link{tdmClassifyLoop}}, an object of class \code{\link{TDMclassifier}}.
#'   @param ... -- currently not used --
#' 
#' @seealso   \code{\link{tdmROCR.TDMclassifier}}   \code{\link{tdmROCRbase}}
#' @export
#' @keywords internal
tdmROCR <- function(x, ...)  UseMethod("tdmROCR");
tdmROCR.default <- function(x, ...)  cat("This is tdmROCR.default\n");

## @export
#tdmROCRbase <- function(x,dataset="validation",nRun=1,typ="ROC",...)  UseMethod("tdmROCRbase");
## @export
#tdmROCRbase.default <- function(x,dataset="validation",nRun=1,typ="ROC",...)  cat("This is tdmROCRbase.default\n");
# --- for unclear reasons it does not work to declare this as a class method with extra arguments, problems when building package
# --- ->> so we let tdmROCRbase as a 'normal' function /WK/


######################################################################################
# tdmROCR.TDMclassifier
#
#'   Interactive plot of ROC, lift or other charts for a \code{\link{TDMclassifier}} object.
#'
#'   Brings up a \code{\link[twiddler]{twiddle}} user interface, where the user may select a part of the dataset
#'   ("training" or "validation"), a run number (if \code{\link{Opts}}(x)$NRUN>1) 
#'   and a type-of-chart, see \code{\link{tdmROCRbase}} for details. Using \code{\link{tdmROCRbase}}, 
#'   the appropriate chart is plotted on the current graphics device.
#'
#'   @method tdmROCR TDMclassifier
#'   @param x  return value from a prior call to \code{\link{tdmClassifyLoop}}, an object of class \code{\link{TDMclassifier}}.
#'   @param ... -- currently not used --
#' 
#'   @return The area under the curve plotted most recently.
#'
#'   @note Side effect: Ror each chart, calculate and print the area between the curve and the bottom line (y=1.0 for \code{typ=="lift"}, y=0.0 else).
#'
#' @examples
#'    \dontrun{
#'      setwd(paste(find.package("TDMR"), "demo02sonar",sep="/"));
#'      source("main_sonar.r");
#'      result = main_sonar();
#'      tdmROCR(result);
#'    }
#' @seealso   \code{\link{tdmClassifyLoop}}   \code{\link{tdmROCRbase}}
#' @export
tdmROCR.TDMclassifier <- function(x,...) {
  #require(twiddler);   # Not needed anymore, because we have twiddler now in 'Depends'.
                        # We need to load library twiddler (via 'Depends' or via 'require'), 
                        # because it calls internally the function tcltk::tktoplevel, 
                        # which cannot be resolved otherwise.
  if (length(x$predProbList[[1]]$Val)==1 |
      length(x$predProbList[[1]]$Trn)==1) {
    warning("Object x of class TDMclassifier has no prediction score in 'predProbList' --> can not show a ROC curve");
  } else {
    cat1(Opts(x),"Showing ROC curves for TDMclassifier object\n");
    areaR = 0;
    tdmROCR_ <- function(dataset,nRun,typ) {
      areaR <<- tdmROCRbase(x,dataset,nRun,typ);
      cat(sprintf("Area under curve = %f\n",areaR));
      flush.console();
    }
    tne = length(x$predProbList);
    twiddleCmd <- paste("twiddle(tdmROCR_(dataset,nRun",sep="");
    if (tne==1)   twiddleCmd <- paste(twiddleCmd,"=1",sep="");
    twiddleCmd <- paste(twiddleCmd,",typ",sep="");
    twiddleCmd <- paste(twiddleCmd,"), eval=FALSE",sep="");
                  # eval=FALSE triggers two buttons "EVAL" and "CLOSE"  and inhibits auto-evaluation in twiddle
    twiddleCmd <- paste(twiddleCmd,", dataset=combo(\"validation\",\"training\")",sep="");
    twiddleCmd <- paste(twiddleCmd,", typ=combo(\"ROC\",\"lift\",\"precRec\")",sep="");
    if (tne>1) {
      twiddleCmd <- paste(twiddleCmd,", nRun=knob(c(1,",tne,"), res=1, label=\"nRun\")",sep="");
    }

    twiddleCmd <- paste(twiddleCmd,")",sep="");
    eval(parse(text=twiddleCmd));
    areaR;
  }
}

######################################################################################
# tdmROCRbase
#
#'   Single plot of ROC, lift or other chart for a \code{\link{TDMclassifier}} object.
#'
#  @method tdmROCR TDMclassifier
#'  @param x  return value from a prior call to \code{\link{tdmClassifyLoop}}, an object of class \code{\link{TDMclassifier}}.
#'  @param dataset ["validation"] which part of the data to use, either "training" or "validation"
#'  @param nRun    [1] if x contains multiple runs, which run to show  (1,...,\code{\link{Opts}}(x)$NRUN)
#'  @param typ     ["ROC"] which chart type, one out of ("ROC","lift","precRec") for 
#'                 (ROC, lift, precision-recall)-chart (see \code{\link[ROCR]{performance}} in package ROCR for more details):
#'    \itemize{
#'      \item "ROC":      receiver operating curve, TPR vs. FPR, with TPR=TP/(TP+FN)=TP/P and FPR=FP/(FP+TN)=FP/N (true and false positive rate).
#'      \item "lift":     lift chart, LIFT vs. RPP, with LIFT=TPR/RPR with random positive rate RPR=P/(P+N) and RPP=(TP+FP)/(P+N) (rate of pos. predictions).
#'      \item "precRec":  precision-recall-chart, PREC vs. RECALL, with PREC=TP/(TP+FP) and RECALL=TP/P (same as TPR).
#'     } 
#'  @param noPlot   [FALSE] if TRUE, suppress the plot, return only the area under curve
#'  @param ...      currently not used
#'
#'  @return The area between the curve and the bottom line y=0.0 in the case of \code{typ=="ROC" | typ=="precRec"} \cr
#'      or  the area between the curve and the bottom line y=1.0 in the case of \code{typ=="lift"}. \cr
#'      If object \code{x} does not contain a prediction score, a warning is issued and the return value is NULL.
#'
#'  @example  demo/demo06ROCR.r
#'  
#' @seealso   \code{\link{tdmClassifyLoop}}   \code{\link{tdmROCR.TDMclassifier}}
#' @export
tdmROCRbase  <- function(x,dataset="validation",nRun=1,typ="ROC",noPlot=FALSE,...) {
#tdmROCRbase.TDMclassifier  <- function(x,dataset="validation",nRun=1,typ="ROC",noPlot=FALSE,...) {
  if (!inherits(x, "TDMclassifier"))
        stop("This function is permitted only for objects of class `TDMclassifier'");
  if (length(x$predProbList[[1]]$Val)==1 |
      length(x$predProbList[[1]]$Trn)==1) {
    warning("Object x of class TDMclassifier has no prediction score in 'predProbList' --> can not show a ROC curve");
    NULL;
  } else {
    typList = c("ROC","lift","precRec");
    if (length(which(typ==typList))==0) stop(sprintf("Invalid value for param typ: %s",typ));
    if (nRun>length(x$predProbList) | nRun<1)  stop(sprintf("Invalid value for param nRun: %d",nRun));
    ppVal = switch(dataset
      , "validation" = x$predProbList[[nRun]]$Val
      , "training" = x$predProbList[[nRun]]$Trn
      , "invalidSwitch");
    if (!is.data.frame(ppVal)) stop(sprintf("Invalid value for param dataset: %s",dataset));
    titList = c("ROC","Lift","Precision/Recall");
    ymeas = c("tpr","lift","prec");
    xmeas = c("fpr","rpp","rec");

    perf <- tdmROCR_calc(ppVal,ymeas[typ==typList],xmeas[typ==typList]);
    areaR <- tdmROCR_area(perf,typ);
    if (!noPlot) { 
      ROCR::plot(perf,colorize=T,lwd=2,main=sprintf("%s on %s set",titList[typ==typList],dataset));
      if (.Devices[[dev.cur()]]=="windows") bringToTop();
    }
    areaR;
  }
}

tdmROCR_calc <- function(ppVal,ymeasure,xmeasure) {
    #require(ROCR);  # now via direct call 'ROCR::'
    #
    # TODO: extend for multiple response variables
    #    
    ll = ppVal[,2];         # the true class labels (for the first response variable)
    pp = ppVal[,4];         # the prediction score (for the first response variable)
    lo =levels(ll);
    # estimate which class label has the higher average score: order 'lo' in such a way that this label is the last one:
    if ( mean(pp[ll==lo[1]]) > mean(pp[ll==lo[2]]) ) lo = rev(lo);
    pred = ROCR::prediction(pp,ll,label.ordering=lo);
    perf = ROCR::performance(pred,ymeasure,xmeasure);
}
tdmROCR_area <- function(perf,typ="ROC") {
    #cat(sprintf("AUC = %f\n",ROCR::performance(pred,"auc")@y.values));   # area under curve (always ROC)
    #
    # area for ROC, lift of precision-recall chart:
    baseline = ifelse(typ=="lift",1.0,0.0);
    xv=perf@x.values[[1]];
    yv=perf@y.values[[1]];
    d=length(xv);
    dx=xv[2:d]-xv[1:(d-1)];
    areaR=sum(dx*(yv[2:d]-baseline));
}

#
# --- experimental, not yet ready for export
#
tdmROCR.tdmClass <- function(x,...) {
  #require(ROCR);  # now specific call with 'ROCR::'
  if (is.null(x$d_train$votes)) {
    warning("Object of class tdmClass has no prediction score d_train$votes --> can not show a ROC curve");
  } else {
    cat1(x$opts,"Showing ROC curve for training data (tdmClass)\n");
    ll = x$d_train$Class;         # the true class labels
    pp = x$d_train$votes;         # the prediction score
    lo =levels(ll);
    # estimate which class label has the higher average score: order 'lo' in such a way that this label is the last one:
    if ( mean(pp[ll==lo[1]]) > mean(pp[ll==lo[2]]) ) lo = rev(lo);
    predTr = ROCR::prediction(pp,ll,label.ordering=lo);
    perfTr = ROCR::performance(predTr,"tpr","fpr");
    ROCR::plot(perfTr,colorize=T,lwd=2,main="ROC Training Set")
  }
}

