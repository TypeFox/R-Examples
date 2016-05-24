#############################################################
##### EMRPP Object

setClass("EMRPPObj",
  representation=list(inputData="data.frame",GpSizes="numeric",l_HasExcess="logical",d_ExcessVal="numeric",Hotelling="logical",Commens="logical",
    MaxDist="numeric",Interval="numeric", DistExp="numeric",CForm="numeric",NumGrps="numeric",NumObs="numeric",NumVars="numeric",da_GpVals="numeric",
    da_GroupV="numeric",ObsDelta="numeric",YHot="matrix",P_value="numeric",
    Call="character",group.names="vector",CommAvgDist="vector"),
  prototype=list(inputData=data.frame(),GpSizes=0,l_HasExcess=FALSE,d_ExcessVal=0,Hotelling=FALSE,Commens=TRUE,
    MaxDist=0,Interval=0.0,DistExp=0,CForm=0,NumGrps=0,NumObs=0,NumVars=0,da_GpVals=0,da_GroupV=0,
    ObsDelta=0,YHot=matrix(0.0,0.0,0.0),P_value=0,CommAvgDist=0))


#######################################################
#### EMRPP Summary

setMethod("summary",
    signature(object = "EMRPPObj"),
    function (object, ...)
    {
     cv <- c(.001,.01,.05,.1,1)
     name.vect<-names(object@inputData)

    cat("\tExact Multi-Response Permutation Procedure (EMRPP) \n")
    cat("\nCall: \n")
      cat(paste(strwrap(object@Call,width=70),collapse="\n\t"),"\n")
        if(name.vect[1]!="V1"){
        cat("\n  Grouping Variable : ", name.vect[1])
        cat("\n  Response Variables: ", name.vect[2:length(name.vect)])}

  cat("\n\nSpecification of Analysis:")
                      cat("\n\tNumber of Observations:", object@NumObs, "\n")
                        cat("\tNumber of Groups      :", object@NumGrps, "\n")
                        cat("\tDistance Exponent     :", object@DistExp, "\n")
    if(object@CForm==1) cat("\tWeighting Factor      : n(I)/sum(n(I))=C(I) = ",object@CForm)
    if(object@CForm==2) cat("\tWeighting Factor      : (n(I)-1)/sum(n(I)-1) = C(I) = ",object@CForm)
    if(object@CForm==3) cat("\tWeighting Factor      :  1/sum(1) = C(I) = ",object@CForm)
    if(object@CForm==4) cat("\tWeighting Factor      : (n(I)*(n(I)-1))/sum(n(I)*(n(I)-1)) = C(I) = ",object@CForm)
  if(object@MaxDist!=0) cat("\n\tDistance Cutoff     :",object@MaxDist)
 if(object@Interval!=0) cat("\n\tARC distances used  :",object@Interval, "\n\t     Intervals in unit circle")


  cat("\n\nGroup Summary:\n")
    if(object@l_HasExcess==TRUE){
    
     group.names<-c(object@group.names[-c(which(object@da_GpVals==object@d_ExcessVal,arr.ind=TRUE))],
               paste(object@group.names[(which(object@da_GpVals==object@d_ExcessVal,arr.ind=TRUE))],"*",sep=""))

     GpSizes<-c(object@GpSizes[-c(which(object@da_GpVals==object@d_ExcessVal,arr.ind=TRUE))],
               paste(object@GpSizes[(which(object@da_GpVals==object@d_ExcessVal,arr.ind=TRUE))],"*",sep=""))
     group.table <- cbind(group.names,GpSizes)
     }  else {
        group.table <- cbind(object@group.names,object@GpSizes)
       }

    colnames(group.table) <- list("Group Value","Group Size")
    rownames(group.table)=rep("",times=length(object@GpSizes))
    print.default(group.table,print.gap=2,quote=FALSE,row.names=FALSE)
    if(object@l_HasExcess==TRUE){
     cat("\n\t* Excess Group")}

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
     
  cat("\n\n Results:")
    cat("\n\tDelta Observed             : ",object@ObsDelta)
    cat("\n\tProbability (Exact)")
    cat("\n\tof a smaller or equal delta: ",
      object@P_value)
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
#### EMRPP Print

setMethod("print",
    signature(x = "EMRPPObj"),
    function (x, ...)
    {
    cat("\nCall: \n")
      cat(deparse(x@Call,width.cutoff=40L),sep="\n")

  cat("\n Results:")
    cat("\nDelta Observed: ",x@ObsDelta)
    cv <- c(.001,.01,.05,.1,1)
    cat("\nExact P-value: ",
      x@P_value)
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

setMethod("show", "EMRPPObj", function(object) {
print(object)
})

# a couple accessor functions to follow S4 recommendations
if (!isGeneric("pvalue")) {
  setGeneric("pvalue", function(x){
  standardGeneric("pvalue")
 })
 }

setMethod("pvalue", "EMRPPObj", function(x) {
 x@P_value
 })




