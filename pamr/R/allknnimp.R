"pamr.knnimpute" <-
  function (data, k = 10,rowmax=0.5,colmax=0.8,maxp=1500) 
{
  x <- data$x
  p<-nrow(x)
  col.nas <- drop(rep(1,p)%*%is.na(x))
  if (any(col.nas>colmax*p)) 
    stop(paste("a column has more than",format(round(colmax*100)),"% missing values!"))
  data$x <- knnimp(x,k,maxmiss=rowmax,maxp=maxp)
  data
}
knnimp<-function(x,k=10,maxmiss=0.5,maxp=1500){
  pn<-dim(x)
  dn<-dimnames(x)
  p<-as.integer(pn[1])
  n<-as.integer(pn[2])
  imiss<-is.na(x)
  x[imiss]<-0
  irmiss<-drop(imiss%*%rep(1,n))
  imax<-trunc(maxmiss*n)
  imax<-irmiss>imax
  simax<-sum(imax)
  if(simax>0){
    warning(paste(simax,"rows with more than", format(round(maxmiss*100,1)),"% entries missing;\n",
                  "mean imputation used for these rows"))
    irmiss<-irmiss[!imax]
    imiss.omit<-imiss[imax,,drop=FALSE]
    imiss<-imiss[!imax,]
    xomit<-x[imax,,drop=FALSE]
    x<-x[!imax,]
    discards<-seq(imax)[imax]
    p<-as.integer(p-simax)
  }
  storage.mode(imiss)<-"integer"
  storage.mode(irmiss)<-"integer"
  storage.mode(x)<-"double"
  if(p<=maxp)
    ximp<-knnimp.internal(x,k,imiss,irmiss,p,n,maxp=maxp)
  else
    ximp<-knnimp.split(x,k,imiss,irmiss,p,n,maxp=maxp)
  imiss.new<-is.na(ximp)
  newmiss<-any(imiss.new)
  if( (simax>0) | newmiss ){
    xbar<-mean.miss(x,imiss=imiss)
    if(newmiss)ximp<-meanimp(ximp,imiss.new,xbar)
    if(simax>0){
      xomit<-meanimp(xomit,imiss.omit,xbar)
      xout<-array(0,dim=pn)
      xout[!imax,]<-ximp
      xout[imax,]<-xomit
      ximp<-xout
    }
  }
  dimnames(ximp)<-dn
  ximp
}
           

knnimp.internal<-function(x,k,imiss,irmiss,p,n,maxp=maxp){
  if(p<=maxp){
    junk<-.Fortran("knnimp",
                   x,
                   ximp=x,
                   p,
                   n,
                   imiss=imiss,
                   irmiss,
                   as.integer(k),
                   double(p),
                   double(n),
                   integer(p),
                   integer(n),
                   PACKAGE="pamr")

    ximp<-junk$ximp
### Should we check or iterate?
    ximp[junk$imiss==2]<-NA
    ximp
  }
  else
    knnimp.split(x,k,imiss,irmiss,p,n,maxp=maxp)
}
"knnimp.split" <-
  function(x,k,imiss,irmiss,p,n,maxp){
    junk<-twomeans.miss(x)
    size<-junk$size
    cat("Cluster size",p,"broken into",size,"\n")
    clus<-junk$cluster
    for(i in seq(size)){
      p<-as.integer(size[i])
      index<-clus==i
      x[index,]<-if(p<k)
        meanimp(x[index,])
      else
        knnimp.internal(x[index,],k,imiss[index,],irmiss[index],p,n,maxp)
      cat("Done cluster",size[i],"\n")
    }
    x
  }
mean.miss<-function(x,index=seq(p),imiss=is.na(x)){
  pn<-dim(x)
  p<-pn[1]
  n<-pn[2]
  storage.mode(index)<-"integer"
  x[imiss]<-0
  storage.mode(x)<-"double"
  storage.mode(imiss)<-"integer"
  junk<-  .Fortran("misave",
           x,
           x0=double(n),
           p,
           n,
           imiss0=as.integer(rep(1,n)),
           imiss,
           index,
           as.integer(length(index)),
           PACKAGE="pamr")

  x0<-junk$x0
  x0[junk$imiss0==2]<-NA
x0
}
           
meanimp<-function(x,imiss=is.na(x),xbar=mean.miss(x,imiss=imiss)){
  nr<-nrow(x)
  if(!is.null(nr)&&(nr>1))x[imiss]<-outer(rep(1,nr),xbar)[imiss]
  x
}
                                         
"twomeans.miss" <-
function(x, imiss=is.na(x),imbalance=.2,maxit=5,eps=0.001){
  ### Compute the two-means cluster solution for data with missing
  ### entries
  pn<-dim(x)
  p<-pn[1];n<-pn[2]
  if(missing(imiss))x[imiss]<-0
  storage.mode(imiss)<-"integer"
  starters<-sample(seq(p),2)
  junk<-.Fortran("twomis",
                 x,
                 as.integer(p),
                 as.integer(n),
                 imiss,
                 double(2*n),
                 integer(2*n),
                 as.integer(maxit),
                 as.double(eps),
                 as.integer(starters),
                 cluster=integer(2*p),
                 nsize=integer(2),
                 double(2*p),
                 ratio=double(1),
                 iter=integer(1),
                 integer(p),
                 integer(n),
                PACKAGE="pamr"
               )

  clus=matrix(junk$cluster,ncol=2)
  cluster<-as.list(1:2)
  for(i in 1:2)cluster[[i]]<-clus[seq(junk$nsize[i]),i]
  clus<-rep(1,p)
  clus[cluster[[2]]]<-2
  list(cluster=clus,ratio=junk$ratio,iter=junk$iter,size=junk$nsize)
}
