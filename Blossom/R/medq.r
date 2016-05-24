medq<-function(variables,group,data,save.test=FALSE,quant=c(0.0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1.0)){
 Call<-match.call()

       
      if(missing(data)){x<-cbind(variables,group)
      }else{             
         variables<-as.character(Call$variables)
            if(length(variables>1)) variables<-variables[-c(variables=="c")]
         group<-as.character(Call$group)
            if(sum(length(variables)+ length(group))==0) {orig.x<-data[,c(seq(from=2,to=ncol(data),by=1),1)]
             x<-orig.x
            }
            else {
            
            if(is.list(variables)) variables<-unlist(variables)
                  orig.x<-data[,match(c(variables,group),names(data))]
                  x<-orig.x
                  if(length(group)==0)
                    x<-cbind(x,1)
                  }
       }
    
   x<-x[order(x[,ncol(x)]),]
   group.names<-unique(x[,ncol(x)])
   input.data<-as.data.frame(x)
    #remove incomplete cases
    i_NumObs=nrow(x)
       x<-x[complete.cases(x),]
      comp.cases<-nrow(x)
      if(comp.cases/i_NumObs!=1) {warning(paste(i_NumObs-comp.cases,ifelse(i_NumObs-comp.cases==1," case was removed because of a missing value ",
        " cases were removed because of missing values."),sep=""))
           i_NumObs<-nrow(x)   
        }
   
   x[,ncol(x)]<-as.numeric(x[,ncol(x)])
  NumVars=ncol(x)-1
  NumGrps=length(group.names)
  NumObs=nrow(x)
  NumQuantVals=length(quant)
  GrpSizes=as.vector(table(x[,ncol(x)]))
  MaxGpSize<-max(GrpSizes)

NumIterations=rep(0,NumGrps)
VariableWInGpMedian<-matrix(data=0,nrow=NumGrps,ncol=NumVars)
GpAvgDistToGpMVMedian<-rep(0,times=NumGrps)
WInGpQuantDist<-matrix(data=0,nrow=NumGrps,ncol=NumQuantVals)
ObsDistToGpMedian<-matrix(data=0,nrow=NumGrps,ncol=MaxGpSize)
GpMedQTolerance<-rep(0,NumGrps)
WInGrpVarEstimateVal<-matrix(data=0,nrow=NumGrps,ncol=NumQuantVals)
    
iEr=0
x<-as.matrix(x)
storage.mode(x)<-"double"
storage.mode(GrpSizes)<-"integer"
storage.mode(quant)<-"double"
storage.mode(VariableWInGpMedian)<-"double"
storage.mode(GpAvgDistToGpMVMedian)<-"double"
storage.mode(GpMedQTolerance)<-"double"
storage.mode(WInGrpVarEstimateVal)<-"double"
storage.mode(WInGpQuantDist)<-"double"
storage.mode(ObsDistToGpMedian)<-"double"
storage.mode(NumIterations)<-"integer"

   if(NumVars>65536) stop("Number of variables should not exceed 65536 for medq")
   if(NumVars<1) stop("You must have at least one variables for medq")

medqOut<-.Fortran("wrapmedq",
        as.integer(NumGrps),
					as.integer(NumVars),
					as.integer(NumObs),
                    GrpSizes,
                    x,
                    as.integer(NumQuantVals),
                    quant,
                    as.integer(MaxGpSize),
                    VariableWInGpMedian,
                    GpAvgDistToGpMVMedian,
                    GpMedQTolerance,
                    WInGrpVarEstimateVal,
                    WInGpQuantDist,
                    ObsDistToGpMedian,
                    NumIterations,
                    as.integer(iEr))
             
indx<-c(1,GrpSizes[1:(length(GrpSizes)-1)])
for(i in 1:length(GrpSizes)){ 
    if(all(apply(x[sum(indx[1:i]):sum(GrpSizes[1:i]),(1:ncol(x)-1)],2,var)==0)) {
           medqOut[[9]][i,]<-apply(x[sum(indx[1:i]):sum(GrpSizes[1:i]),(1:ncol(x)-1)],2,median)
           medqOut[[10]][i]<-0
           medqOut[[13]][i,]<-0
        }
    }                    

 medqOut<-new("MEDQObj",NumGrps=NumGrps,NumVars=medqOut[[2]],NumObs=medqOut[[3]],
                    GrpSizes=medqOut[[4]],inputData=input.data,NumQuantVals=medqOut[[6]],quant=medqOut[[7]],
                    MaxGpSize=medqOut[[8]],VariableWInGpMedian=medqOut[[9]],GpAvgDistToGpMVMedian=medqOut[[10]],
                    GpMedQTolerance=medqOut[[11]],WInGpQuantDist=medqOut[[13]],
                    ObsDistToGpMedian=medqOut[[14]],NumIterations=medqOut[[15]],
                    Call=deparse(Call,width.cutoff=200L),group.names=group.names)

medqOut
   }



























