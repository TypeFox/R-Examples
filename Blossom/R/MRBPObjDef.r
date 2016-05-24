#############################################################
##### MRBP Object

setClass("MRBPObj",
  representation=list(inputData="data.frame",DistExp="numeric",NumVars="numeric",NumBlocks="numeric",
          NumGrps="numeric",NumPerm="numeric",NumObs="numeric",
          AlignVals="matrix",Align="logical",Exact="logical",Resample="logical",Commensurate="logical",
          ObsDelta="numeric", ExpectDelta="numeric",DeltaVar="numeric",
          DeltaSkew="numeric",AgreeVal="numeric",StdStat="numeric",P_value="numeric",
          Call="character",CommenAvgDist="vector",
          SaveTest="logical",PermVals="vector",group.names="character"),
  prototype=list(inputData=data.frame(),DistExp=0,NumVars=0,NumBlocks=0,
          NumGrps=0,NumPerm=0,NumObs=0,
          AlignVals=matrix(0,0,0),Align=FALSE,Exact=FALSE,Resample=FALSE,Commensurate=FALSE,
          ObsDelta=0, ExpectDelta=0,DeltaVar=0,
          DeltaSkew=0,AgreeVal=0,StdStat=0,P_value=0,
          Call="",CommenAvgDist="vector",
          SaveTest=FALSE,PermVals=matrix(0,0,0)))


#######################################################
#### MRBP Summary

setMethod("summary",
    signature(object = "MRBPObj"),
    function (object, ...)
    {
      group.size<-as.vector(table(object@inputData[,1]))
      group.value<-unique(object@inputData[,1])
      name.vect<-names(object@inputData)
      cv <- c(.001,.01,.05,.1,1)
      
     if(object@Exact==TRUE) cat("\tExact Multiresponse Randomized Block Procedure (EMRBP) \n")
        else  cat("\tMultiresponse Randomized Block Procedure (MRBP) \n")

    cat("\nCall: \n")
      cat(paste(strwrap(object@Call,width=70),collapse="\n\t"),"\n")

  cat("\nSpecification of Analysis:")

   cat("\n\tResponse Variables      : ", name.vect[3:length(name.vect)] )
   cat("\n\tNumber of Observations  : ", object@NumObs)
   cat("\n\tGrouping Variable       : ", name.vect[1], "\n")
     cat("\t   Number of Groups     : ", object@NumGrps)

   cat("\n\tBlocking Variable       : ", name.vect[2], "\n")
     cat("\t   Number of Blocks     : ", object@NumBlocks)
   cat("\n\tDistance Exponent       : ", object@DistExp)

    if(object@Align==TRUE) cat("\n\tData Aligned")
    if(object@Commensurate==TRUE) cat("\n\tHotelling's Commensuration Applied")

  cat("\n\nGroup Summary: \n")

        group.table <- cbind(group.value,group.size)


    colnames(group.table) <- list("Group Value","Group Size")
    rownames(group.table)=rep("",times=dim(group.table)[1])
    print.default(group.table,print.gap=2,quote=FALSE,row.names=FALSE)
  if(object@Align==TRUE){
      cat("\n\nBlock Alignment Summary:\n")

          block.table<-cbind(as.vector(rbind(unique(object@inputData[,2]),
              matrix(data="",ncol=object@NumBlocks,nrow=(length(name.vect)-3)))),
          rep(name.vect[3:length(name.vect)],times=object@NumBlocks),
          as.vector(t(format(object@AlignVals,digits=getOption("digits")))))

         colnames(block.table) <- list("Block Value","Variable Name","Alignment Value")
         rownames(block.table)=rep("",times=dim(block.table)[1])

        print.default(block.table,print.gap=2,quote=FALSE,row.names=FALSE)
   } else cat("\n\tData are not aligned within blocks")
   
 if(object@Commensurate==TRUE & dim(object@inputData)[2]>=2){
     cat("\n\nVariable Commensuration Summary\n")

            align.table<-cbind(name.vect[3:length(name.vect)],signif(object@CommenAvgDist,digits=getOption("digits")))
            colnames(align.table) <- list("Variable Name","Average Euclidean Distance")
            rownames(align.table)<-rep("",times=dim(align.table)[1])

           print.default(align.table,print.gap=2,quote=FALSE,row.names=FALSE)
  } else cat("\n\tVariables are not commensurated")
  
  cat("\n\n Results: ")
          cat("\n\tDelta Observed                : ",object@ObsDelta)
    if(object@Exact==FALSE & object@Resample==FALSE){
          cat("\n\tDelta Expected                : ",object@ExpectDelta)
          cat("\n\tDelta Variance                : ",object@DeltaVar)
          cat("\n\tDelta Skewness                : ",object@DeltaSkew)

          cat("\n\tAgreement measure among blocks: ",object@AgreeVal)
          cat("\n\tStandardized test statistic   : ",object@StdStat)
    }
      if(object@Exact==TRUE){
          cat("\n\tExact Probability ")
          cat("\n\tof a smaller or equal delta   : ",
      object@P_value)
    }else if(object@Resample==FALSE){
          cat("\n\tProbability (Pearson Type III)")
          cat("\n\tof a smaller or equal delta   : ",
      object@P_value)} else {
          cat("\n\tProbability (Resample)")
          cat("\n\t of a smaller or equal delta  : ",
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
#### MRBP Print

setMethod("print",
    signature(x = "MRBPObj"),
    function (x, ...)
    {
    cv <- c(.001,.01,.05,.1,1)
    cat("\nCall:\n")
      cat(x@Call,"\n")

  cat("\n\n Results:")
    cat("\n\tDelta Observed: ",x@ObsDelta)
    if(x@Exact==FALSE & x@Resample==FALSE){
          cat("\n\tDelta Expected: ",x@ExpectDelta)
          cat("\n\tDelta Variance: ",x@DeltaVar)
          cat("\n\tDelta Skewness: ",x@DeltaSkew)

          cat("\n\t Agreement measure among blocks: ",x@AgreeVal)
          cat("\n\t Standardized test statistic   : ",x@StdStat)
    }
      if(x@Exact==TRUE){
         cat("\n\tExact Probability")
         cat("\n\tof a smaller or equal delta     : ",
      x@P_value)
    }else if(x@Resample==FALSE){
         cat("\n\t Probability (Pearson Type III)")
         cat("\n\t of a smaller or equal delta    : ",
      x@P_value)} else {
         cat("\n\t Probability (Resample)")
         cat("\n\t of a smaller or equal delta    : ",
        x@P_value)}

      switch(min(which(cv>=x@P_value,arr.ind=TRUE)),
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
    invisible(x)
    }
)

setMethod("show", "MRBPObj", function(object) {
print(object)
})


if (!isGeneric("pvalue")) {
  setGeneric("pvalue", function(x){
  standardGeneric("pvalue")
 })
 }

setMethod("pvalue", "MRBPObj", function(x) {
 x@P_value
 })

if (!isGeneric("ResampVals")) {
  setGeneric("ResampVals", function(x){
  standardGeneric("ResampVals")
 })
 }

setMethod("ResampVals", "MRBPObj", function(x) {
 x@PermVals
 })


