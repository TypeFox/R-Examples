#distance computation for distance matrices
gromovdist<-function(d1,d2=NULL,type="l1",p=NULL,...){
  typelist<-c("l1","l2","linfinity","lp")
  pos<-pmatch(type,typelist)
  if (is.na(pos)) stop(paste("type",type,"not implemented")) 
  if ((type=="lp")&&!missing(p)){ if ((p<1)|(!is.numeric(p))) stop(paste("value of p",quote(p)," is not admissible"))
  }
  if (hasArg(d2)&&((c1<-class(d1)[1])!=(c2<-class(d2)[1]))) stop(paste("classes for d1, ",c1,", and d2, ",c2,", differ"))
  UseMethod("gromovdist")
}
gromovdist.default<-function(d1,d2=NULL,type="l1",p=NULL,...){
  gromovdist(d1=as.list(d1),d2=as.list(d2),type=type,p=p,...)
}
gromovdist.list<-function(d1,d2=NULL,type="l1",p=NULL,...){
  if (length(table(sapply(d1,FUN=function(o)class(o)[1])))>1L) stop("not all elements of the list have the same type")
  if (!missing(d2)) warning("argument 2 is ignored") 
  ln<-length(d1)
  mmm<-cbind(rep(1:ln,ln),rep(1:ln,each=ln))
  mmm<-mmm[which(mmm[,1]<mmm[,2]),]
  ll<-length(mmm[,1])
  do<-apply(mmm,FUN=function(i)gromovdist(d1=d1[[i[1]]],d2=d1[[i[2]]],type=type,p=p,...),MARGIN=1)
  class(do)<-c("dissimilarity","dist")
  attr(do,"Size")<-ln
  attr(do,"Diag")<-FALSE
  attr(do,"Upper")<-FALSE
  attr(do,"Metric")<-paste0("gromov",type)
  return(do)
}

gromovdist.matrix<-function(d1,d2=NULL,type="l1",p=NULL,...){
  if (missing(d2)) stop("argument 2 missing") 
  if (is.null(rownames(d1))) {rownames(d1)<-1:(dim(d1)[1])
                              colnames(d1)<-1:(dim(d1)[2])
  }                            
  d1lab<-intersect(rownames(d1),colnames(d1))
  if (is.null(d1lab)) stop("no consistent labeling for argument 1")
  if (is.null(rownames(d2))) {rownames(d2)<-1:(dim(d2)[1])
                              colnames(d2)<-1:(dim(d2)[2])
  }                            
  d2lab<-intersect(rownames(d2),colnames(d2))
  if (is.null(d2lab)) stop("no consistent labeling for argument 2")
  comlab<-intersect(d1lab,d2lab)
  if (!is.null(comlab)) {
    dd1<-d1[comlab,comlab]
    dd2<-d2[comlab,comlab]
    typelist<-c("l1","l2","linfinity","lp")
    pos<-pmatch(type,typelist)
    return(do.call(paste0("gromov",typelist[pos[1]],".intern"),list(d1=dd1,d2=dd2,n=length(comlab),p=p,...) ))
  }
  else 
    stop("no common labels")
}

gromovdist.dist<-function(d1,d2=NULL,type="l1",p=NULL,...){
  if (missing(d2)) stop("argument 2 missing")
  gromovdist(d1=as.matrix(d1),d2=as.matrix(d2),type=type,p=p,...)  
}
gromovdist.dissimilarity<-function(d1,d2=NULL,type="l1",p=NULL,...){
  if (missing(d2)) stop("argument 2 missing")
  gromovdist(d1=as.matrix(d1),d2=as.matrix(d2),type=type,p=p,...)  
}
gromovdist.phylo<-function(d1,d2=NULL,type="l1",p=NULL,...){
  if (missing(d2)) stop("argument 2 missing")
  if (is.null(d1$tip.label)) d1$tip.label<-sapply(1:(d1$Nnode+1),paste0) 
  if (is.null(d2$tip.label)) d2$tip.label<-sapply(1:(d2$Nnode+1),paste0) 
  comlab<-intersect(d1$tip.label,d2$tip.label)
  if (!is.null(comlab)) {
    dd1<-cophenetic(d1)[comlab,comlab]
    dd2<-cophenetic(d2)[comlab,comlab]
    typelist<-c("l1","l2","linfinity","lp")
    pos<-pmatch(type,typelist)
    return(do.call(paste0("gromov",typelist[pos[1]],".intern"),list(d1=dd1,d2=dd2,n=length(comlab),p=p,...)))
  }
  else 
    stop("no common labels")
}
gromovdist.multiPhylo<-function(d1,d2=NULL,type="l1",p=NULL,...){
  d1<-c(d1,d2)
  ln<-length(d1)
  mmm<-cbind(rep(1:ln,ln),rep(1:ln,each=ln))
  mmm<-mmm[which(mmm[,1]<mmm[,2]),]
  ll<-length(mmm[,1])
  do<-apply(mmm,FUN=function(i)gromovdist(d1=d1[[i[1]]],d2=d1[[i[2]]],type=type,p=p,...),MARGIN=1)
  class(do)<-c("dissimilarity","dist")
  attr(do,"Size")<-ln
  attr(do,"Diag")<-FALSE
  attr(do,"Upper")<-FALSE
  attr(do,"Metric")<-paste0("gromovphylo",type)
  return(do)
}


gromovdist.igraph<-function(d1,d2=NULL,type="l1",p=NULL,leavesonly=TRUE,...){
  if (missing(d2)) stop("argument 2 missing")
  vv1<-V(d1)
  if (leavesonly) vv1<-vv1[which(degree(d1)==1)]
  dd1<-shortest.paths(d1,v=vv1,to=vv1)
  d1lab<-rownames(dd1)
  vv2<-V(d2)
  if (leavesonly) vv2<-vv2[which(degree(d2)==1)]
  dd2<-shortest.paths(d2,v=vv2,to=vv2)
  d2lab<-rownames(dd2)
  comlab<-intersect(d1lab,d2lab)
  if (!is.null(comlab)) {
    dd1<-dd1[comlab,comlab]
    dd2<-dd2[comlab,comlab]
    typelist<-c("l1","l2","linfinity","lp")
    pos<-pmatch(type,typelist)
    return(do.call(paste0("gromov",typelist[pos[1]],".intern"),list(d1=dd1,d2=dd2,n=length(comlab),p=p,...)))
  }
  else 
    stop("no common labels")
}



