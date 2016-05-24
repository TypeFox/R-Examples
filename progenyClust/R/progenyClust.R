progenyClust <-
function(data,FUNclust=kmeans,method='gap',score.invert=F,ncluster=2:10,size=10,iteration=100,repeats=1,nrandom=10,...){
  # parameter checking
  if(!method%in%c('gap','score','both')){
    warning(paste0(method,' is not found. The default method "gap" is used instead.'))
    method='gap'
  }
  if(length(ncluster)<2){
    if(method!='score'){
      warning('The number of clusters for comprison is too small to use the "gap" method. Method "score" is used instead.')
      method='both'
    }
    if(sum(ncluster<2)!=0){
      stop("The number of clusters should be greater than 1.")
    }
  }
  if(repeats<1){
    stop("The number of repeats can't be zero or negative. Please use a positive value.")
  }
  if(iteration<1){
    stop("The number of iterations can't be zero or negative. Please use a positive value.")
  }
  if(method!='gap'){
    if(nrandom<3){
      warning('The number of random datasets might be too small to produce unbiased reference scores.')
      if(nrandom<1){
        stop("The number of random datasets can't be zero or negative. Please use a positive value or use the method 'gap' instead.")
      }
    }
  }
  if(method=='gap' & !all((ncluster[2:length(ncluster)]-ncluster[1:(length(ncluster)-1)])==1) & !all((ncluster[2:length(ncluster)]-ncluster[1:(length(ncluster)-1)])==-1)){
    stop('In order to use the method "gap", ncluster should be a continuous sequence of integers. Otherwise, use method "score" instead.')
  }

  # intialization
  ## matrix to record cluster assignment
  cluster=matrix(0,dim(data)[1],length(ncluster))
  ## matrix to record stability score
  score=matrix(0,repeats,length(ncluster))
  rscore=NULL
  mean.score=NULL
  mean.gap=NULL
  sd.score=NULL
  sd.gap=NULL
  # progeny clustering
  # generate scores for the input data
  for (rep in 1:repeats){
    for(k in 1:length(ncluster)){
      cluster[,k]<-FUNclust(data,ncluster[k],...)$cluster
      probmatrix=matrix(0,(size*ncluster[k]),(size*ncluster[k]))
      for(iter in 1:iteration){
        progeny=matrix(0,size*ncluster[k],dim(data)[2])
        for(c in 1:ncluster[k]){
          for(j in 1:dim(data)[2]){
            progeny[((c-1)*size+1):(c*size),j]=data[cluster[,k]==c,j][sample(sum(cluster[,k]==c),size,replace=T)]
          }
        }
        ## cluster progenies
        pcluster<-FUNclust(progeny,ncluster[k],...)$cluster
        for(i in 1:(size*ncluster[k])){
          for(j in 1:(size*ncluster[k])){
            if(pcluster[i]==pcluster[j]){
              probmatrix[i,j]=probmatrix[i,j]+1
            }
          }
        }
      }
      probmatrix=probmatrix/iteration
      # caculation of true classification & false classification probabilities
      trueprob=0
      falseprob=0
      for(c in 1:ncluster[k]){
        trueprob=trueprob+sum(probmatrix[((c-1)*size+1):(c*size),((c-1)*size+1):(c*size)])
      }
      falseprob=sum(probmatrix)-trueprob
      # caculation of the score
      if(score.invert==T){
        score[rep,k]=(falseprob/(size*(ncluster[k]-1)*size*ncluster[k]))/((trueprob-size*ncluster[k])/((size-1)*size*ncluster[k]))
      }else{
        score[rep,k]=((trueprob-size*ncluster[k])/((size-1)*size*ncluster[k]))/(falseprob/(size*(ncluster[k]-1)*size*ncluster[k]))
      }
    }
  }
  # compute cluster selection criteria
  if(method=='gap'|method=='both'){
    mean.gap=sapply(1:(dim(score)[2]-2),function(x){score[,x+1]-score[,x]+score[,x+1]-score[,x+2]})
    if(!is.null(dim(mean.gap))){
      sd.gap=apply(mean.gap,2,sd)
      mean.gap=apply(mean.gap,2,mean)
    }
  }
  if(method=='score'|method=='both'){
    rscore=matrix(0,nrandom,length(ncluster))
    for(rep in 1:nrandom){
      rdata=data
      for(i in 1:dim(rdata)[2]){
        rdata[,i]=runif(dim(rdata)[1],min(data[,i]),max(data[,i]))
      }
      for(k in 1:length(ncluster)){
        rcluster<-FUNclust(rdata,ncluster[k],...)$cluster
        rprobmatrix=matrix(0,(size*ncluster[k]),(size*ncluster[k]))
        for(iter in 1:iteration){
          rprogeny=matrix(0,size*ncluster[k],dim(rdata)[2])
          for(c in 1:ncluster[k]){
            for(j in 1:dim(rdata)[2]){
              rprogeny[((c-1)*size+1):(c*size),j]=rdata[rcluster==c,j][sample(sum(rcluster==c),size,replace=T)]
            }
          }
          rpcluster<-FUNclust(rprogeny,ncluster[k],...)$cluster
          for(i in 1:(size*ncluster[k])){
            for(j in 1:(size*ncluster[k])){
              if(rpcluster[i]==rpcluster[j]){
                rprobmatrix[i,j]=rprobmatrix[i,j]+1
              }
            }
          }
        }
        rprobmatrix=rprobmatrix/iteration
        rtrueprob=0
        rfalseprob=0
        for(c in 1:ncluster[k]){
          rtrueprob=rtrueprob+sum(rprobmatrix[((c-1)*size+1):(c*size),((c-1)*size+1):(c*size)])
        }
        rfalseprob=sum(rprobmatrix)-rtrueprob
        if(score.invert==T){
          rscore[rep,k]=(rfalseprob/(size*(ncluster[k]-1)*size*ncluster[k]))/((rtrueprob-size*ncluster[k])/((size-1)*size*ncluster[k]))
        }else{
          rscore[rep,k]=((rtrueprob-size*ncluster[k])/((size-1)*size*ncluster[k]))/(rfalseprob/(size*(ncluster[k]-1)*size*ncluster[k]))
        }
      }
    }
    mean.score=apply(score,2,mean)-apply(rscore,2,mean)
    sd.score=apply(score,2,sd)
  }
  # export
  colnames(cluster)=paste0('C',ncluster)
  output=list(cluster=cluster,score=score,random.score=rscore,mean.gap=mean.gap,mean.score=mean.score,sd.gap=sd.gap,sd.score=sd.score,
              call=match.call(),ncluster=ncluster,method=method,score.invert=score.invert)
  class(output)='progenyClust'
  return(output)
}
