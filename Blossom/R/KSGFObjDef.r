#############################################################
##### KSGF Object

setClass("KSGFObj",
  representation=list(inputData="vector",ObsDelta="numeric",ExpectDelta="numeric",VarDelta="numeric",
              StandDelta="numeric",DeltaSkew="numeric",
              P_value="numeric",DoArc="logical",ArcInterv="numeric",NumCases="numeric",
              Call="character"),
  prototype=list(inputData=0,ObsDelta=0,ExpectDelta=0,VarDelta=0,
              StandDelta=0,DeltaSkew=0,
              P_value=0,DoArc=FALSE,ArcInterv=0,NumCases=0)
 )
#######################################################
#### MRPP Summary

setMethod("summary",
    signature(object = "KSGFObj"),
    function (object, ...)
    {

     cv <- c(.001,.01,.05,.1,1)


    cat("\tKendall-Sherman Goodness of Fit Test \n")
    cat("\nCall: \n")
      cat(paste(strwrap(object@Call,width=70),collapse="\n\t"),"\n")

  cat("\n\nSpecification of Analysis:")

        cat("\n\tNumber of Observations    :", object@NumCases, "\n")
          cat("\tNumber of Intervals       :", object@NumCases-object@DoArc, "\n")
     if(object@DoArc==TRUE) 
          cat("\tArc distance used         :", object@ArcInterv, "\n")
  cat("\n\n Results:")
        cat("\n\tObserved Statistic T      : ",object@ObsDelta)
        cat("\n\tExpected Statistic T      : ",object@ExpectDelta)
        cat("\n\tVariance of Statistic T   : ",object@VarDelta)
        cat("\n\tStandardized statistic T  : ",object@StandDelta)
        cat("\n\tSkewness of statistic T   : ",object@DeltaSkew)
      cat("\n\n\tP-value of observed statistic,")
        cat("\n\t P(Expected >= Observed)  : ",object@P_value)
      switch(min(which(cv>=object@P_value,arr.ind=TRUE)),
  cat("***\n"),
  cat("**\n"),
  cat("*\n"),
  cat(".\n"),
  cat("\n")
)
     Signif <- symnum(cv, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
        cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
    invisible(object)
    }
)


setMethod("print", signature(x="KSGFObj"), function(x) {
summary(x)
})

setMethod("show", "KSGFObj", function(object) {
print(object)
})

if (!isGeneric("pvalue")) {
  setGeneric("pvalue", function(x){
  standardGeneric("pvalue")
 })
 }

setMethod("pvalue", "KSGFObj", function(x) {
 x@P_value
 })



