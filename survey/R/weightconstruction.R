##
##  Functions for constructing nonresponse weights
##
##
##  nonresponse(): constructor
##
##  sparseCells(): identify cells with low count, high non-response weight, high total weight.
##
##  neighbours(): find neighbours of specified cells
##
##  joinCells():  combine two cells.
##
##  reweight(): take a nonreponse object and a set of pweights or a survey design, and
##           produce an updates set of weights or design
##
##  weights(): extract nonresponse weights
##
##  After constructing an object with nonresponse(), use sparse() to find
##  cells that need combining, neighbours() to find things to combine them with,
##  and joinCells() to combine them. Rinse and repeat.
##  Use weights() to extract the final weights
##


nonresponse<-function(sample.weights,sample.counts,population){

  if (!all.equal(dimnames(sample.weights),dimnames(sample.counts)))
    stop("Counts and weights have different dimensions")
  if (!all.equal(dimnames(sample.weights),dimnames(population)))
    stop("sample and population dimensions do not match")

  index<-array(1:length(sample.weights),dim=dim(sample.weights),
               dimnames=dimnames(sample.weights))
  
  rval<-list(weights=sample.weights,
             counts=sample.counts,
             population=population,
             call=sys.call(),
             index=index,
             joins=NULL)
  class(rval)<-"nonresponse"
  rval
  
}


print.nonresponse<-function(x,digits=max(3,getOption("digits")-4),...,max.print=5){
  cat("Call: ")
  print(x$call)
  n<-length(x$index)
  nunq<-length(unique(as.vector(x$index)))
  cat(n,"original cells, ",nunq,"distinct cells remaining\n")
  show.all<-nunq<=max.print
  counts<-rowsum(as.vector(x$counts), as.vector(x$index),reorder=FALSE)
  ratios<-rowsum(as.vector(x$population),as.vector(x$index),reorder=FALSE)/rowsum(as.vector(x$weights),as.vector(x$index),reorder=FALSE)
  totals<-rowsum(as.vector(x$population),as.vector(x$index),reorder=FALSE)/counts

  if(length(x$joins)){
    cat("Joins:\n")
    lapply(x$joins, cat,"\n")
  }
  if (show.all){
    cat("Counts: ")
    cat(as.vector(signif(counts,digits)))
    cat("\nNR weights: ")
    cat(as.vector(signif(ratios,digits)))
    cat("\nTotal weights: ")
    cat(as.vector(signif(totals,digits)),"\n")
  } else{ 
    print(summary(data.frame(counts=counts,NRweights=ratios,totalwts=totals)))
  }
  invisible(NULL)
}

"[.nonresponse"<-function(x,i,..){
  expand<-as.integer(as.factor(x$index))
  
  counts<-rowsum(as.vector(x$counts), as.vector(x$index),reorder=FALSE)
  
  ratios<-rowsum(as.vector(x$population),as.vector(x$index),reorder=FALSE)/rowsum(as.vector(x$weights),as.vector(x$index),reorder=FALSE)
  
  totals<-rowsum(as.vector(x$population),as.vector(x$index),reorder=FALSE)/counts

  mm<-matrix(nrow=length(i),ncol=length(dim(x$index)))
  ii<-i-1
  for(j in 1:length(dim(x$index))){
    mm[,j]<-dimnames(x$index)[[j]][(ii %% dim(x$index)[j])+1]
    ii<-ii/(dim(x$index)[j])
  }
  colnames(mm)<-names(dimnames(x$index))
  rownames(mm)<-i
  
  rval<-list(index=as.vector(x$index)[i], names=mm, totals=totals[expand][i],ratios=ratios[expand][i],counts=counts[expand][i])
  rval$call<-sys.call()
  class(rval)<-"nonresponseSubset"
  rval
}

print.nonresponseSubset<-function(x,digits=max(2, getOption("digits")-4),...){
  print(x$call)
  cat("Cells: ")
  cat(x$index)
  cat("\nIndices:\n")
  print(x$names)
  cat("Summary:\n")
  mm<-cbind(NRwt=signif(x$ratios,digits),wt=signif(x$totals,digits),n=x$counts)
  rownames(mm)<-x$index
  prmatrix(mm)
  invisible(NULL)
}

sparseCells<-function(object,count=0,totalweight=Inf,nrweight=1.5){

  expand<-as.integer(as.factor(object$index))

  counts<-rowsum(as.vector(object$counts), as.vector(object$index),reorder=FALSE)

  ratios<-rowsum(as.vector(object$population),as.vector(object$index),reorder=FALSE)/rowsum(as.vector(object$weights),as.vector(object$index),reorder=FALSE)

  totals<-rowsum(as.vector(object$population),as.vector(object$index),reorder=FALSE)/counts

  bad<- (ratios[expand]>nrweight | counts[expand]<count | totals[expand]>totalweight)
  i<-which(bad & !duplicated(as.vector(object$index)))
  if (length(i)==0)
    return(NULL)
  
  d<-dim(object$weights)
  nd<-length(d) 
  dd<-cumprod(d)
  dn<-dimnames(object$weights)
  
  ii<-t(t(outer(i-1,dd,"%%")) %/% c(1,dd[-nd])) +1
  keep<-!duplicated(i)
  ii<-ii[keep,,drop=FALSE]
  i<-i[keep]

  mm<-matrix("",ncol=NCOL(ii),nrow=NROW(ii))
  colnames(mm)<-names(dimnames(object$index))
  rownames(mm)<-i
  
  for(j in seq(length=length(d))){
    mm[,j]<-dn[[j]][ii[,j]]
  }

  rval<-list(index=i, names=mm, totals=totals[expand][i],ratios=ratios[expand][i],
             counts=counts[expand][i],call=sys.call())
  class(rval)<-"nonresponseSubset"
  rval
}


neighbours<-function(index,object){

  d<-dim(object$index)
  if (length(index)==1)
    i<-object$index[index]
  else
    i<-sum(c(index-1,0)*cumprod(c(1,d)))+1
  nd<-length(d)
  
  ## all the joins of that index 
  ii<-which(object$index==object$index[i],arr.ind=TRUE)
  nii<-nrow(ii)
  diffs<-t(rbind(diag(rep(1,nd)),diag(rep(-1,nd))))
  nbours<-matrix(nrow=2*nd*nii,ncol=nd)
  counter<-0
  for(j in 1:nii){
    nbours[counter+1:(2*nd),]<-t(ii[j,]+diffs)
    counter<-counter+2*nd
  }
  keep<-apply(nbours,1, function(x) all(x>0 & x<=d))
  keep<- keep & !duplicated(nbours)
  nbours<-nbours[keep,,drop=FALSE]
  
  nbour.index<-apply(nbours,1, function(x) sum(c(x-1,0)*cumprod(c(1,d)))+1)

  nbour.index<-nbour.index[!(object$index[nbour.index] %in% object$index[index])]
  
  object[nbour.index]
}



weights.nonresponse<-function(object,...){
  w<-rowsum(as.vector(object$population),as.vector(object$index))/rowsum(as.vector(object$weights),as.vector(object$index))
  expand<-as.integer(as.factor(object$index))
  array(w[expand],dim=dim(object$index),dimnames=dimnames(object$index))
}

joinCells<-function(object,a,...){
  if (!is.list(a))
    a<-list(a,...)

  d<-dim(object$index)
  nd<-length(d)
  if (length(a[[1]])>1)
    a<-sapply(a,function(ai) sum(c(ai-1,0)*cumprod(c(1,d)))+1)
  else
    a<-do.call("c",a)
  nd<-length(d)
  
  if(length(a)<2){
    warning("Can't join a single cell")
    return(invisible(object))
  }

  indices<-object$index[a]
  if (length(unique(indices))<2){
    warning("These cells are already joined")
    return(invisible(object))
  }

  object$index[object$index %in% object$index[a]]<-min(object$index[a])
  object$joins<-c(object$joins,list(which(object$index %in% object$index[a])))
  object
  
}
