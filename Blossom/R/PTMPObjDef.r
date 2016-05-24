#############################################################
##### PTMP Object

setClass("PTMPObj",
  representation=list(NumPairs="numeric",expon="numeric",
    Data1="numeric",Data2="numeric",ExpectDelta="numeric",DeltaVar="numeric",
    DeltaSkew="numeric",StdStat="numeric",Rho="numeric",ObsDelta="numeric",
    P_value="numeric",Resample="logical",Exact="logical",NumPerms="numeric",
    Call="character",SaveTest="logical",PermVals="vector",NamesUsed="character",
    GroupNames="character",BlockNames="character",AlignVals="vector"),
  prototype=list(NumPairs=0,expon=0,
    PowerOfRanks=0,ExponOfRanks=0,Data1=0,Data2=0,ExpectDelta=0,DeltaVar=0,
    DeltaSkew=0,StdStat=0,Rho=0,ObsDelta=0,
    P_value=0,Resample=FALSE,Exact=FALSE,NumPerms=0,
    SaveTest=FALSE,PermVals=matrix(0,0,0)))

#######################################################
#### PTMP Summary

setMethod("summary",
    signature(object = "PTMPObj"),
    function (object, ...)
    {
     cv <- c(.001,.01,.05,.1,1)
     if(object@Exact==TRUE) cat("\tExact Permutation Tests for Matched Pairs (PTMP) \n")
        else  cat("\tPermutation Tests for Matched Pairs (PTMP) \n")
    cat("\nCall: \n")
      cat(paste(strwrap(object@Call,width=70),collapse="\n\t"),"\n")

       if(object@NamesUsed[1]!="") cat("\n   Grouping Variable: ", object@NamesUsed[1])
       if(object@NamesUsed[2]!="") cat("\n   Blocking Variable: ", object@NamesUsed[2])
       if(object@NamesUsed[3]!="") cat("\n   Response Variable: ", object@NamesUsed[3])
               NumNonZeroPairs<-length(object@Data1)
      
  cat("\n\nSpecification of Analysis:\n")
  cat("\n\tNumber of Observations        : ", 2*object@NumPairs, "\n")
    cat("\tNumber of Groups              : ", 2,"\n")
    cat("\tNumber of Pairs               : ", object@NumPairs, "\n")
    cat("\tNumber of Non-Zero Differences: ", NumNonZeroPairs, "\n")
    cat("\tDistance Exponent             : ", object@expon, "\n")

     cat("\n\nGroup Summary:\n")
                    
     group.table <- cbind(object@GroupNames,NumNonZeroPairs)

    colnames(group.table) <- list("Group Value","Group Size")
    rownames(group.table)=rep("",times=2)
    print.default(group.table,print.gap=2,quote=FALSE,row.names=FALSE)

  cat("\n\n Results:")

  cv <- c(.001,.01,.05,.1,1)
    cat("\n\tDelta Observed              : ",object@ObsDelta)
  if(object@Exact==FALSE & object@Resample==FALSE){
    cat("\n\tDelta Expected              : ",object@ExpectDelta)
    cat("\n\tDelta Variance              : ",object@DeltaVar)
    cat("\n\tDelta Skewness              : ",object@DeltaSkew)
  cat("\n\n\tAgreement measure           : ", object@Rho)
  }
    if(object@Exact==TRUE){
     cat("\n\tExact Probability")    
     cat("\n\tof a smaller or equal delta: ",
      object@P_value)
    }else if(object@Resample==FALSE){
    cat("\n\tProbability (Pearson Type III)")
    cat("\n\tof a smaller or equal delta : ",
      object@P_value)} else {
    cat("\n\tProbability (Resample)")
    cat("\n\tof a smaller or equal delta : ",
        object@P_value)}

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
#### PTMP Print

setMethod("print",
    signature(x = "PTMPObj"),
    function (x, ...)
    {
    cat("\nCall: \n")
      cat(deparse(x@Call,width.cutoff=40L),sep="\n")
    cat("\n\n Results:")
    cv <- c(.001,.01,.05,.1,1)
    cat("\n\tDelta Observed: ",x@ObsDelta)
    if(x@Exact==FALSE & x@Resample==FALSE){
    cat("\n\tDelta Expected: ",x@ExpectDelta)
    cat("\n\tDelta Variance: ",x@DeltaVar)
    cat("\n\tDelta Skewness: ",x@DeltaSkew)

     cat("\n\n\tAgreement measure          : ",x@Rho)
     }
    if(x@Exact==TRUE){
       cat("\n\tExact Probability")
       cat("\n\tof a smaller or equal delta: ",
      x@P_value)
    }else if(x@Resample==FALSE){
    cat("\n\tProbability (Pearson Type III)")
    cat("\n\tof a smaller or equal delta   : ",
      x@P_value)} else {
    cat("\n\t Probability (Resample)")
    cat("\n\t of a smaller or equal delta  : ",
        x@P_value)}

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

setMethod("show", "PTMPObj", function(object) {
print(object)
})

if (!isGeneric("pvalue")) {
  setGeneric("pvalue", function(x){
  standardGeneric("pvalue")
 })
 }

setMethod("pvalue", "PTMPObj", function(x) {
 x@P_value
 })

if (!isGeneric("ResampVals")) {
  setGeneric("ResampVals", function(x){
  standardGeneric("ResampVals")
 })
 }

setMethod("ResampVals", "PTMPObj", function(x) {
 x@PermVals
 })


