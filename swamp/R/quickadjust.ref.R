quickadjust.ref <-
function(g,batches,refbatch){

  if(class(g)!="matrix"){stop("g is not a matrix")}
  if(class(batches)!="factor"){stop("batches is not a factor")}
  if(length(levels(batches))<1.5){stop("batches has to be a factor with at least two levels")}
  if(length(batches)!=ncol(g)){stop("batches has not the same length as ncol(g)")}
  if (any(table(batches)<1.5)){stop("a level of batches has one or zero counts")}
  if(length(refbatch)!=1){stop("refbatch does not have length 1")}
  if(class(refbatch)!="character"){stop("refbatch is not a character")}
  if(any(levels(batches)==refbatch)==FALSE){stop("refbatch is not a level of batches")}
  isna<-which(is.na(batches))
  if (length(isna)>0) {warning(paste("Samples",toString(isna),"will not be adjusted because of NAs in batches"))   }
  
  gafter<-g
  gafter[]<-NA
  adjustervalues<-matrix(ncol=length(levels(batches)),nrow=nrow(g),dimnames=list(rownames(g),levels(batches)))
  refindex<-which(batches==refbatch)
  gref<-g[,refindex]
  rowMedianref<-apply(gref, 1, median, na.rm = T)
  

  for (i in 1:length(levels(batches))){
  index<-which(batches==levels(batches)[i])  
  b1<-g[,index]
  
    b1adj<-b1
    b1adj[]<-NA
    for (j in 1:nrow(b1)){
            adjuster<-rowMedianref[j]/median(b1[j,],na.rm=T)
            adjustervalues[j,i]<-adjuster
            b1adj[j,]<-b1[j,]*adjuster
            }
  gafter[,index]<-b1adj
  }
  return(list(adjusted.data=gafter,scaling.factors=adjustervalues))
  }
