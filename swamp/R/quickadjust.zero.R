quickadjust.zero <-
function(g,batches){

  if(class(g)!="matrix"){stop("g is not a matrix")}
  if(class(batches)!="factor"){stop("batches is not a factor")}
  if(length(levels(batches))<1.5){stop("batches has to be a factor with at least two levels")}
  if(length(batches)!=ncol(g)){stop("batches has not the same length as ncol(g)")}
  if (any(table(batches)<1.5)){stop("a level of batches has one or zero counts")}
  isna<-which(is.na(batches))
  if (length(isna)>0) {warning(paste("Samples",toString(isna),"will not be adjusted because of NAs in batches"))   }
 
  
  gafter<-g
  gafter[]<-NA
  for (i in 1:length(levels(batches))){
  index<-which(batches==levels(batches)[i])
  gb<-g[,index]
  rowMedian<-apply(gb, 1, median, na.rm = T) # median centering not mean centering (more robust)
  k1<-(gb-rowMedian)
  gafter[,index]<-k1
  }
  return(gafter)
  }

