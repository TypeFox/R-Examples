#############################################################
##### MRPP Object
setClass("MRPPObj",
  representation=list(NumObs="numeric",NumVars="numeric",
    DistExp="numeric",MaxDist="numeric",CForm="numeric",Hotelling="logical",Commens="logical",NumGrps="numeric",Interval="numeric",
    GpSizes="numeric",NumPerm="numeric",DoResamp="logical",inputData="data.frame",AvgDist="numeric",
    StandTestStat="numeric",ObsDelta="numeric",ExpectDelta="numeric",DeltaVar="numeric",
    DeltaSkew="numeric",P_value="numeric",YHot="matrix",
    d_ExcessVal="numeric",l_HasExcess="logical",da_GpVals="numeric",ia_GrpValTag="numeric",
    da_GroupV="numeric",Call="character",group.names="vector",SaveTest="logical",PermVals="vector",CommAvgDist="vector"),
  prototype=list(NumObs=0,NumVars=0,
    DistExp=0,MaxDist=0,CForm=0,Hotelling=FALSE,Commens=TRUE,NumGrps=0,Interval=0.0,
    GpSizes=0,NumPerm=0,DoResamp=FALSE,inputData=data.frame(),AvgDist=0,
    StandTestStat=0,ObsDelta=0,ExpectDelta=0,DeltaVar=0,
    DeltaSkew=0,P_value=0,YHot=matrix(0,0,0),
    d_ExcessVal=0,l_HasExcess=FALSE,da_GpVals=0,ia_GrpValTag=0,
    da_GroupV=0,SaveTest=FALSE,PermVals=matrix(0,0,0),CommAvgDist=0))

#######################################################
#### MRPP Summary

setMethod("summary",
    signature(object = "MRPPObj"),
    function (object, ...)
    {

     cv <- c(.001,.01,.05,.1,1)
     name.vect<-names(object@inputData)

    cat("\tMulti-Response Permutation Procedure (MRPP) \n")

    cat("\nCall: \n")
    cat(paste(strwrap(object@Call,width=70),collapse="\n\t"),"\n")
      if(name.vect[1]!="V1"){
      cat("\n  Grouping Variable : ", name.vect[1])
      cat("\n  Response Variables: ", name.vect[2:length(name.vect)])}

  cat("\n\nSpecification of Analysis:")

                          cat("\n\tNumber of Observations:", object@NumObs, "\n")
                            cat("\tNumber of Groups      :", object@NumGrps, "\n")
                            cat("\tDistance Exponent     :", object@DistExp, "\n")
        if(object@CForm==1) cat("\tWeighting Factor: n(I)/sum(n(I))=C(I)= ",object@CForm)
        if(object@CForm==2) cat("\tWeighting Factor: (n(I)-1)/sum(n(I)-1) = C(I) = ",object@CForm)
        if(object@CForm==3) cat("\tWeighting Factor: 1/sum(1) = C(I) = ",object@CForm)
        if(object@CForm==4) cat("\tWeighting Factor: (n(I)*(n(I)-1))/sum(n(I)*(n(I)-1)) = C(I) = ",object@CForm)
    if(object@MaxDist!=0) cat("\n\tDistance Cutoff       :",object@MaxDist)
  if(object@Interval!=0){ cat("\n\tARC distances used    :",object@Interval)
                          cat("\n\tIntervals in unit circle")
                         }
  cat("\n\nGroup Summary:\n")

     if(object@l_HasExcess==TRUE){
     group.names<-c(object@group.names[-c(which(object@da_GpVals==object@d_ExcessVal,arr.ind=TRUE))],
               paste(object@group.names[(which(object@da_GpVals==object@d_ExcessVal,arr.ind=TRUE))],"*",sep=""))
     AvgDist<-c(signif(object@AvgDist,digits=getOption("digits")),"")

     GpSizes<-c(object@GpSizes[-c(which(object@da_GpVals==object@d_ExcessVal,arr.ind=TRUE))],
               paste(object@GpSizes[(which(object@da_GpVals==object@d_ExcessVal,arr.ind=TRUE))],"*",sep=""))
     group.table <- cbind(group.names,GpSizes,AvgDist)
     }  else {
        group.table <- cbind(object@group.names,object@GpSizes,signif(object@AvgDist,digits=getOption("digits")))
       }
    colnames(group.table) <- list("Group Value","Group Size","Group Distance")
    rownames(group.table)=rep("",times=length(object@GpSizes))
    print.default(group.table,print.gap=2,quote=FALSE,row.names=FALSE)



  if(object@Hotelling==TRUE) {
    cat("\n  Hotelling's Commensuration Applied\n")

    cat("   Variance/Covariance Matrix:\n")
    colnames(object@YHot)=rep("",times=dim(object@YHot)[2])
    rownames(object@YHot)=paste(rep("     Variable",times=length(name.vect)-1),
        seq(from=1,to=length(name.vect)-1),rep(":",times=length(name.vect)-1),
        name.vect[2:length(name.vect)])
    print.default(format(object@YHot,digits=getOption("digit")),quote=FALSE)


    } else{
     if(object@Commens==TRUE & dim(object@inputData)[2]>2){

     cat("\n\nVariable Commensuration Summary\n")
             group.size<-as.vector(table(object@inputData[,1]))
             group.value<-unique(object@inputData[,1])

            align.table<-cbind(name.vect[2:length(name.vect)],signif(object@CommAvgDist,digits=getOption("digits")))
            colnames(align.table) <- list("Variable Name","Average Distance (Euclidian if V=1)")
            rownames(align.table)<-rep("",times=dim(align.table)[1])

           print.default(align.table,print.gap=2,quote=FALSE,row.names=FALSE)
  } else if(object@Commens==FALSE) cat("\n\tVariables are not commensurated")}

  
    if(object@l_HasExcess==TRUE){
     cat("\n\t* Excess Group")}
  cat("\n\n Results:")
        cat("\n\tDelta Observed                : ",object@ObsDelta)
    
    if(object@DoResamp==FALSE){
        cat("\n\tDelta Expected                : ",object@ExpectDelta)
        cat("\n\tDelta Variance                : ",object@DeltaVar)
        cat("\n\tDelta Skewness                : ",object@DeltaSkew)
      cat("\n\n\tStandardized test statistic   : ", object@StandTestStat)

        cat("\n\tProbability (Pearson Type III)")
        cat("\n\tof a smaller or equal delta   : ",
          object@P_value)} else {
        cat("\n\tProbability (Resample)")
        cat("\n\tof a smaller or equal delta   : ",
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

#################################################################
#### MRPP Print

setMethod("print",
    signature(x = "MRPPObj"),
    function (x, ...)
    {
    cat("\nCall: \n")
     cat(deparse(x@Call,width.cutoff=40L),sep="\n")

  cat("\nGroup Summary:\n")
    if(x@l_HasExcess==TRUE){
     group.names<-c(x@group.names[-c(which(x@da_GpVals==x@d_ExcessVal,arr.ind=TRUE))],
               paste(x@group.names[(which(x@da_GpVals==x@d_ExcessVal,arr.ind=TRUE))],"*",sep=""))
     AvgDist<-c(signif(x@AvgDist,digits=getOption("digits")),"")
     GpSizes<-c(x@GpSizes[-c(which(x@da_GpVals==x@d_ExcessVal,arr.ind=TRUE))],
               paste(x@GpSizes[(which(x@da_GpVals==x@d_ExcessVal,arr.ind=TRUE))],"*",sep=""))
     group.table <- cbind(group.names,GpSizes,AvgDist)
     }  else {
        group.table <- cbind(x@group.names,x@GpSizes,signif(x@AvgDist,digits=getOption("digits")))
       }

    colnames(group.table) <- list("Group Value","Group Size","Group Distance")
    rownames(group.table)=rep("",times=length(x@GpSizes))
    print.default(group.table,print.gap=2,quote=FALSE,row.names=FALSE)
    if(x@l_HasExcess==TRUE){
     cat("\n\t* Excess Group")}
     
  cat("\n\n Results:")
    cat("\nDelta Observed: ",x@ObsDelta)
    cat("\nDelta Expected: ",x@ExpectDelta)
    cat("\nDelta Variance: ",x@DeltaVar)
    cat("\nDelta Skewness: ",x@DeltaSkew)
    cat("\nStandardized test statistic: ", x@StandTestStat)
    cv <- c(.001,.01,.05,.1,1)
    if(x@DoResamp==FALSE){
    cat("\n\tProbability (Pearson Type III) ")
    cat("\n\tof a smaller or equal delta    : ",
      x@P_value)} else {
    cat("\n\tProbability (Resample)")
    cat("\n\tof a smaller or equal delta    : ",
        x@P_value)}

      switch(min(which(cv>=x@P_value,arr.ind=TRUE)),
  cat("***\n"),
  cat("**\n"),
  cat("*\n"),
  cat(".\n"),
  cat("\n")
)
    invisible(x)
    }
)


setMethod("show", "MRPPObj", function(object) {
print(object)
})

# Now setting up just a few generic C style accessor functions as recommended
# in the S4 class documentation

if (!isGeneric("pvalue")) {
  setGeneric("pvalue", function(x){
  standardGeneric("pvalue")
 })
 }
 
setMethod("pvalue", "MRPPObj", function(x) {
 x@P_value
 })
 
if (!isGeneric("ResampVals")) {
  setGeneric("ResampVals", function(x){
  standardGeneric("ResampVals")
 })
 }

setMethod("ResampVals", "MRPPObj", function(x) {
 x@PermVals
 })
 
 

