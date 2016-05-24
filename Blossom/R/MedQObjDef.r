#############################################################
##### MEDQ Object

setClass("MEDQObj",
  representation=list(NumGrps="numeric",NumVars="numeric",NumObs="numeric",
                    GrpSizes="vector",inputData="data.frame",NumQuantVals="numeric",quant="vector",
                    MaxGpSize="vector",VariableWInGpMedian="matrix",GpAvgDistToGpMVMedian="vector",
                    GpMedQTolerance="vector",WInGpQuantDist="matrix",
                    ObsDistToGpMedian="matrix",NumIterations="vector",
                    Call="character",group.names="vector"),

   prototype=list(NumGrps=0,NumVars=0,NumObs=0,
                    GrpSizes=0,inputData=data.frame(),NumQuantVals=0,quant=0,
                    MaxGpSize=0,VariableWInGpMedian=matrix(0,0,0),GpAvgDistToGpMVMedian=0,
                    GpMedQTolerance=0,WInGpQuantDist=matrix(0,0,0),
                    ObsDistToGpMedian=matrix(0,0,0),NumIterations=0,
                    group.names=0))
                    


#######################################################
#### MEDQ Summary

setMethod("summary",
    signature(object = "MEDQObj"),
    function (object, ...)
    {
   name.vect<-colnames(object@inputData)

    cat(paste("\t",object@NumVars,"-Dimensional Median and Distance Quantiles \n",sep=""))
    cat("\nCall: \n")
      cat(paste(strwrap(object@Call,width=70),collapse="\n\t"),"\n")

  cat("\nSpecification of Analysis:\n")

    cat("\n\tGrouping Variable         : ", name.vect[ncol(object@inputData)])
    cat("\n\tNumber of Report Variables: ", object@NumVars)
    cat("\n\tReport Variables          : ", paste(name.vect[1:(ncol(object@inputData)-1)],collapse=", "), "\n")
    
    cat("\n\tNumber of Observations    : ", object@NumObs)
    cat("\n\tNumber of Groups          : ", object@NumGrps)
 #Within Group Summaries
 for(i in 1:object@NumGrps){
 if(object@NumGrps>1){cat("\n_________________________________________________\n")

       cat("\nResults for Group Value: ", as.character(object@group.names[i]),"\n")
          cat("\n\tObservations in Group     : ", object@GrpSizes[i])
       }
          cat("\n\tIterations to Solution    : ", object@NumIterations[i])
          cat("\n\tSolution Tolerance        : ", object@GpMedQTolerance[i])
        cat("\n\n",ifelse(object@NumGrps>1,"   Within Group","   "), "Median Coordinates for Variables\n")
        
       group.table<-cbind(name.vect[1:(ncol(object@inputData)-1)],format(object@VariableWInGpMedian[i,],digits=getOption("digits")))
            colnames(group.table) <- list("Variable Name","Multivariate Median Coordinate")
            rownames(group.table)<-rep("       ",times=dim(group.table)[1])
            print.default(group.table,print.gap=2,quote=FALSE,row.names=FALSE)
        cat(paste("\n    ",object@NumVars,"-Dimensional Distance From Median Quantiles",sep=""))
           cat("\n\t",ifelse(object@NumGrps>1,"","Group"), "Average Distance to Multivariate Median: ",object@GpAvgDistToGpMVMedian[i],"\n")
          dist.table<-cbind(object@quant,"",format(object@WInGpQuantDist[i,],digits=getOption("digits")))
            dist.table[object@quant==0,2]="(Minimum)"
            dist.table[object@quant==.5,2]="(Median)"
            dist.table[object@quant==1,2]="(Maximum)"
            colnames(dist.table) <- list("Quantile","","Distance from Median")
            rownames(dist.table)<-rep("       ",times=dim(dist.table)[1])
            print.default(dist.table,print.gap=2,quote=FALSE,row.names=TRUE)


        }
   }

)


setMethod("print", "MEDQObj", function(x) {
summary(x)
})

setMethod("show", "MEDQObj", function(object) {
print(object)
})

if (!isGeneric("Dist2mvm")) {
  setGeneric("Dist2mvm", function(x){
  standardGeneric("Dist2mvm")
 })
 }

setMethod("Dist2mvm", "MEDQObj", function(x) {

Dist2MVM<-as.vector(t(x@ObsDistToGpMedian))
v<-as.vector(rbind(x@GrpSizes,x@MaxGpSize-x@GrpSizes))
v<-rep(rep(c(1,0),times=x@NumGrps),times=v)
Dist2MVM<-Dist2MVM[v!=0]
Dist2MVM<-cbind(x@inputData,Dist2MVM)
if(x@NumGrps==1) Dist2MVM<-Dist2MVM[,-3]
row.names(Dist2MVM)<-NULL
Dist2MVM
 })
