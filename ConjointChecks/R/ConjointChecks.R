##main function
omni.check<-function(N,n,n.iter,burn=1000,thin=4,CR) {#this checks both single and double cancellation
  n/N->dat
  chain<-list()
  #initialize
  inits<-dat
  inits->foo1->foo2
  for (i in 1:nrow(foo1)) foo1[i,]<-i
  for (i in 1:ncol(foo2)) foo2[,i]<-i
  foo1+foo2->index
  sort(as.numeric(dat))->hold
  counter<-1
  for (i in unique(as.numeric(index))) {
    grep(paste("^",i,"$",sep=""),index)->index2
    stop.here<-length(index2)
    counter:(counter+stop.here-1)->replace
    inits[index2]<-hold[replace]
    counter<-max(replace)+1
  }
  inits->old#chain[[1]]
  like<-function(theta,N,n) {
    #choose(N,n)->x1 #this part cancels!
    1->x1
    n*log(theta)->x2
    (N-n)*log(1-theta)->x3
    sum(x2+x3)
  }
  old.ll<-inits
  for (i in 1:nrow(old.ll)) for (j in 1:ncol(old.ll)) like(inits[i,j],N[i,j],n[i,j])->old.ll[i,j]
  dc.counter<-hands.bl<-hands.tr<-list()
  #iterate
  for (I in 2:n.iter) {
    for (i in 1:nrow(dat)) for (j in 1:ncol(dat)) {
      #############################
      #get left hand
      if (j==1) lh1<-0 else lh1<-old[i,j-1]
      if (i==1) lh2<-0 else lh2<-old[i-1,j]
      #get right hand
      if (j==ncol(dat)) rh1<-1 else rh1<-old[i,j+1]
      if (i==nrow(dat)) rh2<-1 else rh2<-old[i+1,j]
      #now make sure double cancellation needs to hold
      test.1 <- as.logical(old[2,1] < old[1,2])
      test.2 <- as.logical(old[3,2] < old[2,3])
      lh3<-0
      rh3<-1
      if (test.1!=test.2) {
        #dc.counter[[paste(I,i,j)]]<-0
      } else {
        if (test.1 & test.2) {
          if (i==1 & j==3) {
            lh3<-old[3,1]
            #dc.counter[[paste(I,i,j)]]<-"a1"
          }
          if (i==3 & j==1) {
            rh3<-old[1,3]
            #dc.counter[[paste(I,i,j)]]<-"a2"
          }
        }
        if (!test.1 & !test.2) {
          if (i==3 & j==1) {
            lh3<-old[1,3]
            #dc.counter[[paste(I,i,j)]]<-"b1"
          }
          if (i==1 & j==3) {
            rh3<-old[3,1]
            #dc.counter[[paste(I,i,j)]]<-"b2"
          }
        }
      }
      #now work everything out all nice like...
      lh<-max(lh1,lh2,lh3)
      if (rh3>lh) rh<-min(rh1,rh2,rh3) else rh<-min(rh1,rh2)
      if (rh<lh) rh<-1
      #sample new point
      runif(1,lh,rh)->draw
      #
      ## if (i==3 & j==1) hands.bl[[I]]<-c(lh3,rh3,lh,rh)
      ## if (i==1 & j==3) hands.tr[[I]]<-c(lh3,rh3,lh,rh)
      ## if ((i==1 & j==3) | (i==3 & j==1)) {
      ##   print(c(i,j))
      ##   print(c(lh1,lh2,lh3))
      ##   print(c(rh1,rh2,rh3))
      ##   print(c(lh,rh))
      ##   print(old)
      ## }
      #acceptance ratio
      ar<-2
      like(draw,N[i,j],n[i,j])->new.ll
      if (!(old[i,j] %in% 0:1)) exp(new.ll-old.ll[i,j])->ar
      if (ar>runif(1)) {
        draw->old[i,j]
        new.ll->old.ll[i,j]
      }
    }
    if (I>burn & I%%4==0) old->chain[[as.character(I)]]
  }
  hi<-lo<-M<-chain[[1]]
  for (i in 1:3) for (j in 1:3) {
    post<-sapply(chain,function(x) x[i,j])
    quantile(post,CR[1])->lo[i,j]
    quantile(post,CR[2])->hi[i,j]
    mean(post)->M[i,j]    
  }
  list(low=lo,high=hi,mean=M)
}

############################################################
#glorified wrapper
ConjointChecks<-function(N,n,n.3mat=1,par.options=NULL,CR=c(.025,.975),seed=NULL) {
  #N is the number of total tries per cell
  #n is the number correct
  #processing function
  proc.fun<-function(dummy,arg.list) {
    arg.list[[1]]->N
    arg.list[[2]]->n
    arg.list[[3]]->lof
    arg.list[[4]]->CR
    lof[[1]]->omni.check
    #lof[[2]]->chain.2.ci
    #lof[[3]]->compare
    test<-1
    nrow(N)->nrows
    ncol(N)->ncols
    while (test>0) {
      sort(sample(1:nrows,3,replace=FALSE))->rows
      sort(sample(1:ncols,3,replace=FALSE))->cols
      N[rows,cols]->nt
      n[rows,cols]->nc
      nc/nt->dat
      sum (dat==1|dat==0)->test
    }
    omni.check(nt,nc,n.iter=3000,CR=CR)->out
    #chain.2.ci(out)->out
    list(rows,cols,out)->save.dat
    save.dat
  }
  proc.fun_adjacent<-function(dummy,arg.list) {
    dummy[1]->r1
    r1:(r1+2)->rows
    dummy[2]->c1
    c1:(c1+2)->cols
    arg.list[[1]]->N
    arg.list[[2]]->n
    arg.list[[3]]->lof
    arg.list[[4]]->CR
    lof[[1]]->omni.check
    #lof[[2]]->chain.2.ci
    #lof[[3]]->compare
    N[rows,cols]->nt
    n[rows,cols]->nc
    nc/nt->dat
    sum (dat==1|dat==0)->test
    if (test>0) NULL->save.dat else {
      omni.check(nt,nc,n.iter=3000,CR=CR)->out
      #chain.2.ci(out)->out
      list(rows,cols,out)->save.dat
    }
    save.dat
  }
  out<-list()
  n/N->dat
  ifelse(abs(dat-.5)<=.5,TRUE,FALSE)->test
  if (!all(test)) stop("There is a problem with n/N, values not between 0 and 1 (inclusive)")
  #list(omni.check,chain.2.ci,compare)->lof
  list(omni.check)->lof
  #if (!is.null(par.options)) {#sequential last
  if (!require(parallel)) stop("Package 'parallel' not available.")
  if (is.null(par.options$n.workers)) par.options$n.workers<-1
  if (is.null(par.options$type)) {
    par.options$type<-"PSOCK"
  } else {
    if (par.options$type=="MPI") {
      if (!require(snow)) stop("Package 'snow' not available.")
      if (!require(Rmpi)) stop("Package 'Rmpi' not available.")
      #stop("here!")
    }
  }
  list(N,n,lof,CR)->arg.list
  cl<-makeCluster(par.options$n.workers,type=par.options$type)
  if (!is.null(seed)) clusterSetRNGStream(cl,iseed=seed) else clusterSetRNGStream(cl)
  dummy<-list()
  if (n.3mat=="adjacent") {
    nrow(N)->nr
    ncol(N)->nc
    for (i in 1:(nr-2)) for (j in 1:(nc-2)) c(i,j)->dummy[[paste(i,j)]]
    clusterApply(cl,dummy,proc.fun_adjacent,arg.list=arg.list)->out
    sapply(out,is.null)->destroy
    out[!destroy]->out
  } else {
    for (i in 1:n.3mat) dummy[[i]]<-i
    clusterApply(cl,dummy,proc.fun,arg.list=arg.list)->out
  }
  stopCluster(cl)             
  list(N=N,n=n,Checks=out)->out
  #now do some summarizing
  compare<-function(dat,lim) {
    ifelse(dat<lim[[2]] & dat>lim[[1]],0,1)->tab
    tab
  }
  out$N->N
  out$n->n
  out$Checks->out
  (n/N)->dat
  mat.num<-matrix(0,nrow(N),ncol(N))
  mat.den<-matrix(-1,nrow(N),ncol(N))
  for (i in 1:length(out)) {
    out[[i]]->x
    x[[1]]->ro
    x[[2]]->co
    compare(dat[ro,co],x[[3]])->comp
    for (i in 1:3) for (j in 1:3) {
      1+mat.den[ro[i],co[j]]->mat.den[ro[i],co[j]]
      comp[i,j]+mat.num[ ro[i],co[j] ]->mat.num[ro[i],co[j]]
    }
  }
  ifelse(mat.den<0,NA,mat.den)->mat.den
  mat.den+1->mat.den
  mat.num/mat.den->tab
  mean(tab,na.rm=TRUE)->m1
  weight<-N
  sum(tab*weight,na.rm=TRUE)/sum(weight)->m2
  list(n=n,N=N,tab=tab,mean=m1,wt.mean=m2,tab.checked=mat.den)
  new("checks", N=N,n=n,Checks=out,tab=tab,means=list(unweighted=m1,weighted=m2),check.counts=mat.den)
}


