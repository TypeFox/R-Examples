#############################################################
##### MRSP Object

setClass("MRSPObj",
  representation=list(NumObs="numeric",NumVars="numeric",
    DistExp="numeric",NumPerm="numeric",DoResamp="logical",Exact="logical",Commens="logical",inputData="matrix",
    TestStat="numeric",ObsDelta="numeric",ExpectDelta="numeric",DeltaVar="numeric",
    DeltaSkew="numeric",RhoAgreement="numeric",P_value="numeric",
    Call="character",SaveTest="logical",PermVals="vector",CommAvgDist="vector"),
  prototype=list(NumObs=0,NumVars=0,
    DistExp=0,NumPerm=0,DoResamp=FALSE,Exact=FALSE,Commens=FALSE,inputData=matrix(0,0,0),
    TestStat=0,ObsDelta=0,ExpectDelta=0,DeltaVar=0,
    DeltaSkew=0,RhoAgreement=0,P_value=0,
    SaveTest=FALSE,PermVals=matrix(0,0,0)))

#######################################################
#### MRSP Summary

setMethod("summary",
    signature(object = "MRSPObj"),
    function (object, ...)
    {
     cv <- c(.001,.01,.05,.1,1)
      name.vect<-names(object@inputData)
      
    if(object@Exact==TRUE) cat("\tExact Multiresponse Sequence Procedure (EMRSP) \n")
      else cat("\tMultiresponse Sequence Procedure (MRSP) \n")
    cat("\nCall: \n")
      cat(paste(strwrap(object@Call,width=70),collapse="\n\t"),"\n")

  cat("\nSpecification of Analysis:")
   cat("\n\tNumber of Observations:", object@NumObs, "\n")
     cat("\tDistance Exponent     :", object@DistExp, "\n")


   if(object@Commens==TRUE & dim(object@inputData)[2]>=2){

     cat("\n\nVariable Commensuration Summary\n")
             group.size<-as.vector(table(object@inputData[,1]))
             group.value<-unique(object@inputData[,1])

              if(is.null(name.vect)){
                name.vect<-c(paste("Variable",seq(1:object@NumVars)))
                }
              
            align.table<-cbind(name.vect,format(object@CommAvgDist,digits=getOption("digits")))
            colnames(align.table) <- list("Variable Name","Average Euclidian Distance")
            rownames(align.table)<-rep("",times=dim(align.table)[1])

           print.default(align.table,print.gap=2,quote=FALSE,row.names=FALSE)
  } else if(object@Commens==FALSE) cat("\n\tVariables are not commensurated")


  cat("\n\n Results:")
    cat("\n\tDelta Observed: ",object@ObsDelta)

    if(object@Exact==FALSE & object@DoResamp==FALSE){
      cat("\n\tDelta Expected: ",object@ExpectDelta)
      cat("\n\tDelta Variance: ",object@DeltaVar)
      cat("\n\tDelta Skewness: ",object@DeltaSkew)
    cat("\n\n\tStandardized test statistic    : ", object@TestStat)
    cat("\n\n\tAgreement measure              : ", object@RhoAgreement)
      cat("\n\tProbability (Pearson Type III)") 
      cat("\n\tof a smaller or equal delta    : ",object@P_value)
    } else { if(object@DoResamp==TRUE){
      cat("\n\tProbability (Resample) ")
      cat("\n\tof a smaller or equal delta    : ",
            object@P_value)}
   else {
      cat("\n\tExact Probability") 
      cat("\n\tof a smaller or equal delta    : ",object@P_value)}
            }
      switch(min(which(cv>=object@P_value,arr.ind=TRUE)),
  cat("***\n"),
  cat("**\n"),
  cat("*\n"),
  cat(".\n"),
  cat("\n")
)
    #as.vector(seq(from=0,to=1,length=5))
     Signif <- symnum(cv, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
        cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
    invisible(object)
    }
)

#################################################################
#### MRSP Print

setMethod("print",
    signature(x = "MRSPObj"),
    function (x, ...)
    {
    cat("\nCall: \n")
       cat(paste(strwrap(x@Call,width=70),collapse="\n\t"),"\n")
     cv <- c(.001,.01,.05,.1,1)

      
  cat("\n\n Results:")
 cat("\n\tDelta Observed: ",x@ObsDelta)

    if(x@Exact==FALSE & x@DoResamp==FALSE){
        cat("\n\tDelta Expected                : ",x@ExpectDelta)
        cat("\n\tDelta Variance                : ",x@DeltaVar)
        cat("\n\tDelta Skewness                : ",x@DeltaSkew)

      cat("\n\n\tStandardized test statistic   : ", x@TestStat)
      cat("\n\n\tAgreement measure             : ",x@RhoAgreement)
        cat("\n\tProbability (Pearson Type III)")
        cat("\n\tof a smaller or equal delta   : ",
          x@P_value)
          } else { if(x@DoResamp==TRUE){
        cat("\n\tProbability (Resample)")
        cat("\n\tof a smaller or equal delta   : ",
            x@P_value)}
              else {
        cat("\n\tExact Probability")
        cat("\n\tof a smaller or equal delta   : ",
              x@P_value)}
            }
      switch(min(which(cv>=x@P_value,arr.ind=TRUE)),
  cat("***\n"),
  cat("**\n"),
  cat("*\n"),
  cat(".\n"),
  cat("\n")
)
    #as.vector(seq(from=0,to=1,length=5))
     Signif <- symnum(cv, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
        cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
    invisible(x)
    }
)
    

setMethod("show", "MRSPObj", function(object) {
print(object)
})

if (!isGeneric("pvalue")) {
  setGeneric("pvalue", function(x){
  standardGeneric("pvalue")
 })
 }

setMethod("pvalue", "MRSPObj", function(x) {
 x@P_value
 })

if (!isGeneric("ResampVals")) {
  setGeneric("ResampVals", function(x){
  standardGeneric("ResampVals")
 })
 }

setMethod("ResampVals", "MRSPObj", function(x) {
 x@PermVals
 })


