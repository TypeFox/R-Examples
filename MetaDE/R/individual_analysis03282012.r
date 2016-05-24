#-------------------------------------------------------------#
# Functions for individual analysis
# Authors: Jia Li and Xingbin Wang
# Institution: University of pittsburgh
# Date: 03/08/2011
#-------------------------------------------------------------#
library(survival)
#-------------------------------------------------------------#
# ind.analysis
# X: expr object
# method=: individual analysis method 
#          regt
#          modt
#          pairedt
#          r
#          F
#          logrank
# ...
# nperm
# tail:  'low' 'high' 'abs'
# var.equal: T or F for F test only
#---------------------------------------------------------------#

ind.analysis<-function(x,ind.method=c("f","regt","modt","pairedt","pearsonr","spearmanr","F","logrank"),miss.tol=0.3,nperm=NULL,tail,...)
{   
   #------impute the missing values ----------------------#
	x<-MetaDE.impute(x,y=miss.tol)
	check.indmethod(x,ind.method)

  ind.method<-match.arg(ind.method,c("regt","modt","pairedt","pearsonr","spearmanr","F","logrank"),several.ok=TRUE)
  #tail<-check.tail(test)  
   x<-check.exp(x)
  # extract expression matrix and labels
  #dat<-x[[1]]
  #lbl<-x[[2]]   
  #--------test for individual study----#
  K<-length(x) # # of studies
  G<-nrow(x[[1]][[1]]) # # of genes
  P<-stat<-P.perm<-NULL #put p values from each datasets together, each column represents a study
  N<-n<-NULL
  for (k in 1:K)
  {
   y<-x[[k]][[1]]
   l<-x[[k]][[2]]
   ind.res<-switch(ind.method[k],regt={get.p.regt(y,l,nperm=nperm,tail=tail)}
            ,modt={get.p.modt(y,l,nperm=nperm,tail=tail)},F={get.p.F(y,l,nperm=nperm,tail=tail)}
            ,pairedt={get.p.pairedt(y,l,nperm=nperm,tail=tail)},pearsonr={get.p.pearsonr(y,l,nperm=nperm,tail=tail)}
            ,spearmanr={get.p.spearmanr(y,l,nperm=nperm,tail=tail)}
            ,logrank={get.p.logrank(y,l,x[[k]]$censoring.status,nperm=nperm,tail=tail)})

   P<- cbind(P,(ind.res$obs)[,2]) #p value
   stat<-cbind(stat,(ind.res$obs[,1])) # observed stats
   if (!is.null(nperm)) P.perm<-cbind(P.perm,c(ind.res$perm)) #p value from permutation 
   cat("dataset",k,"is done\n")
      nns<-get.sample.label.number(x[[k]][[2]],ind.method[k])
	N<-c(N,nns$N)
	n<-c(n,nns$n)
  }
  if(is.null(names(x))){colnames(stat)<-colnames(P)<-paste("dataset",1:K,sep="")
	}else{colnames(stat)<-colnames(P)<-names(x)}
  if (!is.null(nperm)) colnames(P.perm)<-paste("dataset",1:K,sep="")
  all.res<-list(stat=stat,p=P,bp=P.perm)
  attr(all.res$stat,"nperstudy")<-N
  attr(all.res$stat,"nperlabelperstudy")<-n 
  attr(all.res$stat,"individual.analysis")<-ind.method
  return(all.res)
}

#--function to get the number of sample and lables in each study
get.sample.label.number<-function(lbl,method){
if (method=="logrank"){
  N<-length(lbl) #nperstudy
  n<-NA
 	}else{
  N<-length(lbl) #nperstudy
  n<-c(length(lbl)-sum(lbl),sum(lbl)) #n per label per study
  }
return(list(N=N,n=n))
}
get.sample.label.number2<-function(lbl){
  		N<-length(lbl) #nperstudy
  		n<-table(lbl) #n per label per study
	return(list(N=N,n=n))
	}
#--------------------------------------------------------#
#  perform test on individual study                      #
#--------------------------------------------------------#

#-- moderated t for all genes--------#
cal.modt<-function(y,l) {
    stopifnot(length(l)==ncol(y), is.matrix(y))
    l<-unclass(factor(l))
    n<-table(factor(l))
    ind<-diag(rep(1,length(n)))[l,]
    ym<-y%*%ind%*%diag(1/n) #group mean
    a<-(1/n[1]+1/n[2])/(n[1]+n[2]-2)
    s<-sqrt(a*((y^2%*%ind)%*%diag(1/(n-1))-ym^2%*%diag(n/(n-1)))%*%(n-1)) #variance
    s0<-median(s)
    t<-(ym[,1]-ym[,2])/(s0+s)
    return(t)
}

#-------get p value for mod t using permutation-----------------#

get.p.modt<-function(y,l,nperm,tail){
    rnum<-which(apply(y,1,function(x) !any(is.na(x))))
    stat.obs<-p<-rep(NA,nrow(y))
 
    names(stat.obs)<-names(p)<-rownames(y)

    stat.obs[rnum]<-cal.modt(y[rnum,],l) #calculate stats
    if (!is.null(nperm))
    {
     stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
     rownames(stat.perm)<-rownames(p.perm)<-rownames(y)
     stat.perm[rnum,]<-replicate(nperm,cal.modt(y[rnum,],factor(sample(l))))# get stats from permutation
     p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail) #calculate p value from permuation
     p.perm[rnum,]<-emp(stat.perm[rnum,],tail)  # calculate p value of permutation matrix
     res <- list(obs=cbind(stat.obs,p),perm=p.perm)
    }else
    {
     stop("there're no parametric results for moderated t statistic")   
    }
    return(res)
}


#---- regular two sample t-----#
cal.regt<-function(y,l) {
    stopifnot(length(l)==ncol(y), is.matrix(y))
    l<-unclass(factor(l))
    n<-table(factor(l))
    ind<-diag(rep(1,length(n)))[l,]
    ym<-y%*%ind%*%diag(1/n) #group mean
    a<-(1/n[1]+1/n[2])/(n[1]+n[2]-2)
    s<-sqrt(a*((y^2%*%ind)%*%diag(1/(n-1))-ym^2%*%diag(n/(n-1)))%*%(n-1)) #variance
    t<-(ym[,1]-ym[,2])/s
    return(t)
}

#-------get p value for reg t using permutation-----------------#

get.p.regt<-function(y,l,nperm,tail){
    rnum<-which(apply(y,1,function(x) !any(is.na(x))))
    stat.obs<-p<-q<-rep(NA,nrow(y))
    names(stat.obs)<-names(p)<-rownames(y)
    
    n<-table(factor(l))
    df<-n[1]+n[2]-2

    stat.obs[rnum]<-cal.regt(y[rnum,],l) #calculate stat

   if (!is.null(nperm))
    {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)

    stat.perm[rnum,]<-replicate(nperm,cal.regt(y[rnum,],factor(sample(l))))
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail)
    q[rnum]<-p.adjust(p[rnum],method="BH")

    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
   }else
   {
    if (tail=="low") p[rnum]<-pt(stat.obs[rnum],df)
    if (tail=="high") p[rnum]<-1-pt(stat.obs[rnum],df)
    if (tail=="abs") p[rnum]<-2*(pmin(pt(stat.obs[rnum],df),1-pt(stat.obs[rnum],df)))
    res<-list(obs=cbind(stat.obs,p))
   }
   return(res)
}


#-------- calculate  F stat -------for all genes----------------#
cal.F<-function(y, l, var.equal=TRUE) {
    stopifnot(length(l)==ncol(y), is.matrix(y))
    K <- nlevels(l)
    l<-unclass(factor(l))
    n<-table(factor(l))
    ind<-diag(rep(1,length(n)))[l,] 
    ym<-y%*%ind%*%diag(1/n) #group mean
    y..<-matrix(rowMeans(y),nrow(y),ncol(y))  #grand mean
    y.<-ym[,l]
    dfb <- K - 1 #between group df K-1

    if(var.equal) {
        dfw <- ncol(y)-K # n-K within group
        ## mean sum of squares
        num   <- rowSums((y.-y..)^2) / dfb
        den  <- rowSums(( y - y.)^2) / dfw
        ## F statistic
        f<-num/den
    } else {
        ssb <- (y-y.)^2
        ssm<-ssb%*%ind
        wi <- (1/ssm) %*% diag(n*(n-1))
        u <- rowSums(wi)
        MR <- rowSums((1 - wi/u)^2 %*%diag(1/(n-1)))*1/(K^2-1)
        num <- 1/dfb * rowSums((ym - rowSums(wi*ym)/u)^2 * wi)
        den <- 1+ 2* (K-2)*MR
        f <- num/den 
    }
    return(f)
}

#-------get p value for F test  using permutation-----------------#

get.p.F<-function(y,l,nperm,tail,var.equal=TRUE){
    rnum<-which(apply(y,1,function(x) !any(is.na(x))))

    stat.obs<-p<-q<-rep(NA,nrow(y))
    names(stat.obs)<-names(p)<-rownames(y)

    stat.obs[rnum]<-cal.F(y[rnum,],l,var.equal=var.equal)#statistic

   if (!is.null(nperm))
  {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)

    stat.perm[rnum,]<-replicate(nperm,cal.F(y[rnum,],factor(sample(l)),var.equal=var.equal))#stats from perm
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail) #p value from perm
    #q[rnum]<-p.adjust(p[rnum],method="BH")

    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
   }else
   {
    K <- nlevels(l)
    dfb <- K - 1
    dfw <- ncol(y)-K
    if (tail=="high") p[rnum]<-1-pf(stat.obs[rnum],dfb,dfw)
    if (tail=="abs"|tail=="low") stop("Only upper tail is valid for F test")
    res<-list(obs=cbind(stat.obs,p))
   }
   
   return(res)
}

#------------calculate r statistic (pearson's correlation)for all genes----------------#
#---note: l= is continuous------------#
cal.pearsonr<-function(y,l) {
    stopifnot(length(l)==ncol(y), is.matrix(y))
    n<-length(l)
    num<-n*y%*%l-rowSums(y)*sum(l)
    den<-sqrt(n*rowSums(y^2)-rowSums(y)^2)*sqrt(n*sum(l^2)-sum(l)^2)
    r<-num/den
    r
}

#-------get p value for r using permutation-----------------#

get.p.pearsonr<-function(y,l,nperm,tail){
    rnum<-which(apply(y,1,function(x) !any(is.na(x))))
    stat.obs<-p<-q<-rep(NA,nrow(y))
    names(stat.obs)<-names(p)<-rownames(y)

    stat.obs[rnum]<-cal.pearsonr(y[rnum,],l)#observed stat

   if (!is.null(nperm))
   {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)

    stat.perm[rnum,]<-replicate(nperm,cal.pearsonr(y[rnum,],l)) # stat from perm
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail) #p from perm
    #q[rnum]<-p.adjust(p[rnum],method="BH")
    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
    }else
    {
     n<-length(l)
     t<-stat.obs*sqrt((n-2)/(1-stat.obs^2)) #fisher's transformation
     if (tail=="low") p[rnum]<-pt(t,n-2)
     if (tail=="high") p[rnum]<-1-pt(t,n-2)
     if (tail=="abs") p[rnum]<-2*(pmin(pt(t,n-2),1-pt(t,n-2)))
     res<-list(obs=cbind(stat.obs,p))
    }
   return(res)
}

#------------calculate r statistic (spearman's correlation)for all genes----------------#
#---note: l= is continuous------------#
cal.spearmanr<-function(y,l) {
    stopifnot(length(l)==ncol(y), is.matrix(y))
    n<-length(l)
    y<-t(apply(y,1,rank))
    l<-rank(l)
    num<-n*y%*%l-rowSums(y)*sum(l)
    den<-sqrt(n*rowSums(y^2)-rowSums(y)^2)*sqrt(n*sum(l^2)-sum(l)^2)
    r<-num/den
    r
}

#-------get p value for r using permutation-----------------#
get.p.spearmanr<-function(y,l,nperm,tail){
    rnum<-which(apply(y,1,function(x) !any(is.na(x))))
    stat.obs<-p<-q<-rep(NA,nrow(y))
    names(stat.obs)<-names(p)<-rownames(y)

    stat.obs[rnum]<-cal.spearmanr(y[rnum,],l)#observed stat

   if (!is.null(nperm))
   {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)

    stat.perm[rnum,]<-replicate(nperm,cal.spearmanr(y[rnum,],l)) # stat from perm
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail) #p from perm
    #q[rnum]<-p.adjust(p[rnum],method="BH")
    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
    }else
    {
     n<-length(l)
     t<-stat.obs*sqrt((n-2)/(1-stat.obs^2)) #fisher's transformation
     if (tail=="low") p[rnum]<-pt(t[rnum],n-2)
     if (tail=="high") p[rnum]<-1-pt(t[rnum],n-2)
     if (tail=="abs") p[rnum]<-2*(pmin(pt(t[rnum],n-2),1-pt(t[rnum],n-2)))
     res<-list(obs=cbind(stat.obs,p))
    }
   return(res)
}


#---- paired two sample t-----#
#note note note !!!!!!!!!!!!!!  columns of y should be ordered that can be matched---#
cal.pairedt<-function(y,l) {
    l<-unclass(factor(l))
    n<-table(factor(l))
    nn<-sum(n)
    stopifnot(n[1]==n[2]) # not paired design
    ydiff<-y[,1:n[1]]-y[,(n[1]+1):(n[1]+n[2])]
    den<-sqrt(1/(n[1]-1)*(rowSums(ydiff^2))-1/(n[1]^2-n[1])*(rowSums(ydiff))^2)

    onet<-rowMeans(ydiff)/(den/sqrt(n[1]))
    onet
}

#-------get p value for paired t using permutation-----------------#
get.p.pairedt<-function(y,l,nperm,tail){
    rnum<-which(apply(y,1,function(x) !any(is.na(x))))
    stat.obs<-p<-q<-rep(NA,nrow(y))
 
    names(stat.obs)<-names(p)<-rownames(y)
   
   stat.obs[rnum]<-cal.pairedt(y[rnum,],l)
   if(!is.null(nperm))
   {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)

    stat.perm[rnum,]<-replicate(nperm,cal.pairedt(y[rnum,],factor(perm.lab(l,paired=TRUE))))
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail)
    #q[rnum]<-p.adjust(p[rnum],method="BH")

    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
    }else
    {
     n<-length(l)
     if (tail=="low") p[rnum]<-pt(stat.obs[rnum],n/2-1)
     if (tail=="high") p[rnum]<-1-pt(stat.obs[rnum],n/2-1)
     if (tail=="abs") p[rnum]<-2*(pmin(pt(stat.obs[rnum],n/2-1),1-pt(stat.obs[rnum],n/2-1)))
     res<-list(obs=cbind(stat.obs,p))
    }
   res
}

#-------get logrank z statistic-----------------#
cal.logrank<-function(y,time,event)
{
  get.stat<-function(x,time,event)
  {  
     stat<-summary(coxph(Surv(time,event)~x),method="breslow")$sctest[1]
     stat
  }
   z<-apply(y,1,get.stat,time=time,event=event)
  z
}
cal.p.logrank<-function(y,time,event)
{
  get.p<-function(x,time,event)
  {  
     p<-summary(coxph(Surv(time,event)~x,method="breslow"))$sctest[3]
     p
  }
  p<-apply(y,1,get.p,time=time,event=event)
  p
}

#-------get p value for logrank z using permutation-----------------#
get.p.logrank<-function(y,time,event,nperm,tail)
{
    rnum<-which(apply(y,1,function(x) !any(is.na(x))))
    stat.obs<-p<-pp<-rep(NA,nrow(y))

    stat.obs[rnum]<-cal.logrank(y[rnum,],time=time,event=event)

    names(stat.obs)<-names(p)<-rownames(y)
    if(!is.null(nperm))
    {
    stat.perm<-p.perm<-matrix(NA,nrow=nrow(y),ncol=nperm)
    rownames(stat.perm)<-rownames(p.perm)<-rownames(y)
    stat.perm[rnum,]<-replicate(nperm,cal.logrank(y[rnum,],time=time,event=event))
    p[rnum]<-perm.p(stat.obs[rnum],stat.perm[rnum,],tail)
    #q[rnum]<-p.adjust(p[rnum],method="BH")

    p.perm[rnum,]<-emp(stat.perm[rnum,],tail)
    res <- list(obs=cbind(stat.obs,p),perm=p.perm)
    }else
    {
     p[rnum]<-cal.p.logrank(y[rnum,],time=time,event=event)
     if (tail=="low")  pp[rnum]<-p
     if (tail=="high") pp[rnum]<-1-p
     if (tail=="abs")  pp[rnum]<-2*(pmin(p,1-p))
     res<-list(obs=cbind(stat.obs,p=pp))
    }
   res
}

