formatData <-
function(d,cellLines=NULL,inhibitors=NULL,stimuli=NULL,times=NULL,nodes=NULL,intercept=TRUE,initialIntercept=TRUE,gradients=FALSE) {
  # data is a samples x (4 + nodes) matrix or dataframe
  # column 1 gives the cell line in each sample
  # column 2 gives the inhibitor used in each sample
  # column 3 gives the stimuli used in each sample
  # column 4 gives the time each sample was measured
  # cellLines - a vector that can be used to specify a subset of the cell lines to analyse (default is to use them all)
  # inhibitors - a vector that can be used to specify a subset of the inhibitors to analyse (default is to use them all)
  # stimuli - a vector that can be used to specify a subset of the stimuli to analyse (default is to use them all)
  # times - a vector that can be used to specify a subset of the sample times to analyse as the response (default is to use them all).
  # Note: The entries of d[,4] must be real numbers. Missing values are acceptable and are handled as follows:
  # 1) Missing values in the response are just ignored.
  # 2) For the predictors, if a single timepoint is missing, the predictors are interpolated from the two immediate neighbours.
  # 3) If one of the two immediate neighbours is missing then the response is ignored.
  # 4) UNLESS the predictor in question is for the time zero observation (which is always missing), in which case NA is returned, which is later set to the mean of the predictor during centring. See initialIntercept.
  # nodes - a vector that can be used to specify the indices of a subset of nodes to include in the analysis. Further nodes can be removed from the response in the interventionalInferenceDBN function.
  # intercept - include an intercept parameter in all models?
  # initialIntercept - include an intercept parameter to estimate the level at time zero. Only used if the first sample is included in the response.
  # gradients - If true, changes in concentration are used as the response, rather than the raw values.
  if (length(dim(d))!=2 | dim(d)[2]<5) {stop("d must be a matrix with columns 1 to 4 filled with sample information (cell line, inhibitor, stimuli, time)\n")}
  if (is.null(cellLines)) {cellLines<-unique(d[,1])}
  if (is.null(inhibitors)) {inhibitors<-unique(d[,2])}
  if (is.null(stimuli)) {stimuli<-unique(d[,3])}
  sampleTimes<-sort(unique(as.numeric(d[,4])))
  if (is.null(times)) {times<-sampleTimes}
  if (length(intersect(times,sampleTimes))!=length(times)) {stop("All times must appear in d[,4].\n")}
  sampleTimepoints<-1:length(sampleTimes)
  timepoints<-which(sampleTimes %in% times)
  timeIntervals<-sampleTimes[2:length(sampleTimes)]-sampleTimes[1:(length(sampleTimes)-1)]
  if (is.null(nodes)) {nodes<-1:(dim(d)[2]-4)}
  if (initialIntercept & gradients) {stop("initialIntercept and gradients cannot both be TRUE.\n")}
  if (gradients & 1 %in% timepoints) {timepoints<-setdiff(timepoints,1)}
  dm<-data.matrix(d[,5:(dim(d)[2])])
  y<-matrix(NA,0,length(nodes))
  X0<-matrix(NA,0,intercept+initialIntercept)
  X1<-matrix(NA,0,length(nodes))
  Sigma<-matrix(NA,0,0)
  n<-0
  n.cellLines<-rep(0,length(cellLines))
  n.inhibitors<-rep(0,length(inhibitors))
  n.stimuli<-rep(0,length(stimuli))
  n.timepoints<-rep(0,length(timepoints))
  n.interpolated<-0
  interpolated<-matrix(NA,0,4)
  sampleInfo<-matrix(NA,0,4)
  current.condition<-0
  cond<-NULL
  for (cellLine in cellLines) {
    for (i in inhibitors) {
      for (j in stimuli) {
        current.condition<-current.condition+1
        for (k in timepoints) {
          response<-which(d[,1]==cellLine & d[,2]==i & d[,3]==j & d[,4]==sampleTimes[k])
          if (length(response)>0) {
            if (k==1) {
              predictor<-rep(0,length(nodes))
            } else {
              wh<-which(d[,1]==cellLine & d[,2]==i & d[,3]==j & d[,4]==sampleTimes[k-1])
              if (length(wh)==1) {
                predictor<-dm[wh,nodes]
              } else if (length(wh)>1) {
                predictor<-apply(dm[wh,nodes],2,mean) 
              } else if (!gradients & k>2) {# interpolation
                before<-which(d[,1]==cellLine & d[,2]==i & d[,3]==j & d[,4]==sampleTimes[k-2])
                if (length(before)>1 & length(response)>1) {
                  predictor<-(timeIntervals[k+1]*apply(dm[before,nodes],2,mean)+timeIntervals[k]*apply(dm[response,nodes],2,mean))/(timeIntervals[k+1]+timeIntervals[k])
                } else if (length(before)==1 & length(response)>1) {
                  predictor<-(timeIntervals[k+1]*dm[before,nodes]+timeIntervals[k]*apply(dm[response,nodes],2,mean))/(timeIntervals[k+1]+timeIntervals[k])
                } else if (length(before)>1 & length(response)==1) {
                  predictor<-(timeIntervals[k+1]*apply(dm[before,nodes],2,mean)+timeIntervals[k]*dm[response,nodes])/(timeIntervals[k+1]+timeIntervals[k])
                } else if (length(before)==1 & length(response)==1) {
                  predictor<-(timeIntervals[k+1]*dm[before,nodes]+timeIntervals[k]*dm[response,nodes])/(timeIntervals[k+1]+timeIntervals[k])
                } else {# Give up!
                  predictor<-NULL
                }
              } else {# don't use interpolation and gradients
                predictor<-NULL
              }
            }
            if (!is.null(predictor) & gradients) {
              n<-n+1
              Sigma<-cbind(rbind(Sigma,rep(0,n-1)),rep(0,n))
              if (length(response)==1) {
                y<-rbind(y,(dm[response,nodes]-predictor)/timeIntervals[k-1])
              } else {
                y<-rbind(y,(apply(dm[response,nodes],2,mean)-predictor)/timeIntervals[k-1])
              }
              Sigma[n,n]<-(1/length(response)+1/length(wh))/(timeIntervals[k-1])^2
              if (n>1 && prod(sampleInfo[n-1,]==c(cellLine,i,j,sampleTimes[k-1]))==1) {
                Sigma[n-1,n]<--1/length(wh)/timeIntervals[k-1]/timeIntervals[k-2]
                Sigma[n,n-1]<--1/length(wh)/timeIntervals[k-1]/timeIntervals[k-2]
              }
              X1<-rbind(X1,predictor)
              sampleInfo<-rbind(sampleInfo,c(cellLine,i,j,sampleTimes[k]))
              cond<-c(cond,current.condition)
              n.cellLines[which(cellLines==cellLine)]<-n.cellLines[which(cellLines==cellLine)]+1
              n.inhibitors[which(inhibitors==i)]<-n.inhibitors[which(inhibitors==i)]+1
              n.stimuli[which(stimuli==j)]<-n.stimuli[which(stimuli==j)]+1
              n.timepoints[which(timepoints==k)]<-n.timepoints[which(timepoints==k)]+1
              if (intercept) {X0<-matrix(1,n,1)}
            } else if (!is.null(predictor)) {
              for (r in response) {
                n<-n+1
                y<-rbind(y,dm[r,nodes])
                X1<-rbind(X1,predictor)
                Sigma<-diag(rep(1,n))
                if (k>1 && length(wh)==0) {interpolated<-rbind(interpolated,c(cellLine,i,j,sampleTimes[k]));n.interpolated<-n.interpolated+1}
                sampleInfo<-rbind(sampleInfo,c(cellLine,i,j,sampleTimes[k]))
                cond<-c(cond,current.condition)
                n.cellLines[which(cellLines==cellLine)]<-n.cellLines[which(cellLines==cellLine)]+1
                n.inhibitors[which(inhibitors==i)]<-n.inhibitors[which(inhibitors==i)]+1
                n.stimuli[which(stimuli==j)]<-n.stimuli[which(stimuli==j)]+1
                n.timepoints[which(timepoints==k)]<-n.timepoints[which(timepoints==k)]+1
                if (intercept & initialIntercept & k==1) {
                  X0<-rbind(X0,c(1,1))
                } else if (intercept & initialIntercept & k>1) {
                  X0<-rbind(X0,c(1,0))
                } else if (intercept | (initialIntercept & k==1)) {
                  X0<-rbind(X0,1)
                } else if (initialIntercept & k>1) {
                  X0<-rbind(X0,0)
                }
              }
            }
          }
        }
      }
    }
  }
  row.names(X1)<-NULL
  colnames(interpolated)<-c("Cell line","Inhibitor","Stimuli","Time")
  colnames(sampleInfo)<-c("Cell line","Inhibitor","Stimuli","Time")
  cat("n =",n,"\n")
  cat(current.condition,"conditions:\n")
  for (cellLine in cellLines) {cat("  Cell line",cellLine,":",n.cellLines[which(cellLines==cellLine)],"samples.\n")}
  for (i in inhibitors) {cat("  Inhibitor",i,":",n.inhibitors[which(inhibitors==i)],"samples.\n")}
  for (j in stimuli) {cat("  Stimulus",j,":",n.stimuli[which(stimuli==j)],"samples.\n")}
  for (k in timepoints) {cat("  Time",sampleTimes[k],":",n.timepoints[which(timepoints==k)],"samples.\n")}
  if (n.interpolated>0) {cat(n.interpolated,"predictors produced by interpolation.\n")}
  return(list(y=y,X0=X0,X1=X1,Sigma=Sigma,sampleInfo=sampleInfo,interpolated=interpolated,cond=cond))
}
