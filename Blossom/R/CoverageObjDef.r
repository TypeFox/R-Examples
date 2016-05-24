#############################################################
##### MRPP Object
setClass("CoverageObj",
  representation=list(NumPerms="numeric",DistExp="numeric",NumGrps="numeric",
              GpSizes="vector",inputData="vector",ObsDelta="numeric",VarDelta="numeric",ExpectDelta="numeric",
              DeltaSkew="numeric",Z_value="numeric",Skt="numeric",
              P_value="numeric",PZ="numeric",
              NumObs="numeric",PermVals="numeric",exact="logical",group.names="vector",
              Call="character"),
  prototype=list(NumPerms=0,DistExp=0,NumGrps=0,
              GpSizes=0,inputData=0,ObsDelta=0,DeltaVar=0,ExpectDelta=0,
              DeltaSkew=0,Z_value=0,Skt=0,
              P_value=0,PZ=0,
              NumObs=0,PermVals=0,exact=FALSE,
              group.names=0)
 )
#######################################################
#### MRPP Summary

setMethod("summary",
    signature(object = "CoverageObj"),
    function (object, ...)
    {

     cv <- c(.001,.01,.05,.1,1)

    if(object@exact) cat("\tExact Univariate G-Sample Empirical Coverage Test \n")
    else cat("\tUnivariate G-Sample Empirical Coverage Test \n")
    cat("\nCall: \n")
       cat(paste(strwrap(object@Call,width=70),collapse="\n\t"),"\n")
      cat("\n")
  cat("\n\nSpecification of Analysis:")

                   cat("\n\tNumber of Observations: ", sum(object@GpSizes), "\n")
                     cat("\tNumber of Groups      : ", object@NumGrps, "\n")
                     cat("\tDistance Exponent     : ", object@DistExp, "\n")
    if(!object@exact)cat("\tNumber of Permutations: ",object@NumPerms, "\n")

  cat("\n\nGroup Summary:\n")

        group.table <- cbind(object@group.names,signif(object@GpSizes,digits=getOption("digits")))

    colnames(group.table) <- list("Group Value","Group Size")
    rownames(group.table)=rep("",times=length(object@GpSizes))
    print.default(group.table,print.gap=2,quote=FALSE,row.names=FALSE)


  cat("\n\n Results:")
        cat("\n\tObserved coverage statistic             : ",object@ObsDelta)

    if(object@exact==FALSE){
        cat("\n\tMean of coverage statistic              : ",object@ExpectDelta)
        cat("\n\tEstimated variance of coverage statistic: ",object@VarDelta)
        cat("\n\tStandard deviation of the variance") 
        cat("\n\t of the coverage statistic              : ",object@DeltaSkew)

      cat("\n\n\tObserved standardized coverage statistic: ",object@Z_value)
        cat("\n\tSkewness of observed coverage statistic : ",object@Skt)
        cat("\n\tProbability (Pearson Type III)")
        cat("\n\tof a larger or equal coverage statistic : ",object@P_value)
        cat("\n\tProbability (Resampled)")
        cat("\n\tof a larger or equal coverage statistic : ",
        object@PZ)} 
      else {
       cat("\n\tProbability (Exact)")
       cat("\n\tof a larger or equal coverage statistic  : ",
        object@P_value)}
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


setMethod("print", signature(x="CoverageObj"), function(x) {
summary(x)
})

setMethod("show", "CoverageObj", function(object) {
print(object)
})

if (!isGeneric("pvalue")) {
  setGeneric("pvalue", function(x){
  standardGeneric("pvalue")
 })
 }

setMethod("pvalue", "CoverageObj", function(x) {
 x@P_value
 })

if (!isGeneric("ResampVals")) {
  setGeneric("ResampVals", function(x){
  standardGeneric("ResampVals")
 })
 }

setMethod("ResampVals", "CoverageObj", function(x) {
 x@PermVals
 })

