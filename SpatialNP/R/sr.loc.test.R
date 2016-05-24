`sr.loc.test` <-
function(X,Y=NULL,g=NULL,score=c("sign","rank"),nullvalue=NULL,cond=FALSE,cond.n=1000,na.action=na.fail,...) 
{
 if (all(is.null(Y),is.null(g))) { #there is only X
  DNAME<-deparse(substitute(X))
  X<-na.action(X)
  X<-as.matrix(X)
  g<-as.factor(rep(1,dim(X)[1]))
 }
 else if(!is.null(Y)) { 
             #there are X and Y
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  if(dim(X)[2]!=dim(Y)[2]) stop("X and Y must have the same number of columns")
  DNAME<-paste(deparse(substitute(X)),"and",deparse(substitute(Y)))
  X<-na.action(X)
  Y<-na.action(Y)
  g<-factor(c(rep(1,dim(X)[1]),rep(2,dim(Y)[1])))
  X<-rbind(X,Y)
 }
 else if(!is.factor(g))            #there is a g but it's not a factor
  stop("g must be a factor or NULL")
 else {                            #there are X and g
  DNAME<-paste(deparse(substitute(X)),"by",deparse(substitute(g)))
  X<-as.matrix(X)
  Xandg<-cbind(g,X)
  Xandg<-na.action(Xandg)
  g<-factor(Xandg[,1])
  X<-as.matrix(Xandg[,-1])
  rm(Xandg)
 }

 n<-dim(X)[1]
 p<-dim(X)[2]
 c<-nlevels(g)
 if(!is.null(nullvalue)) {
  if(length(nullvalue)!=p) 
   stop("'nullvalue' must have length equal to the number of columns of 'X'")
 }
 else nullvalue<-rep(0,p)
 X<-sweep(X,2,nullvalue)
 NVAL<-paste("c(",paste(nullvalue,collapse=","),")",sep="")
 if(c==1) names(NVAL)<-"location" 
  else if(c==2) names(NVAL)<-"difference between group locations"
  else names(NVAL)<-"difference between some group locations"

 score=match.arg(score)

 switch(score,
     "sign"=
     {
      if (c==1)
      {
       METHOD<-"One sample location test using spatial signs"
       scoremat<-spatial.signs(X,center=F)
      }
      else
      {
       METHOD<-"Several samples location test using spatial signs"
       scoremat<-spatial.signs(X)
      }
     },
     "rank"=
     {
      if (c==1) 
      {
       METHOD<-"One sample location test using spatial signed ranks"
       if (p>1) V<-signrank.shape(X)
      }
      else 
      {
       METHOD<-"Several samples location test using spatial ranks"
       if (p>1) V<-rank.shape(X)
      }
      if (p==1) V<-diag(1)
      scoremat<-spatial.rank(X%*%solve(mat.sqrt(V)),shape=FALSE)
      c2<-mean(norm(scoremat)^2)
     })
 
 if (c==1) {
  STATISTIC<-switch(score,
     "sign"=
     { 
      n*p*sum(apply(scoremat,2,mean)^2)
     },
     "rank"=
     {
      
      sums<-pairsum(X)%*%mat.sqrt(solve(V))
      ave<-apply(spatial.signs(rbind(sums,X),center=F,shape=F),2,mean)
      rm(sums)
      n*p*sum(ave^2)/(4*c2)
     })
 } # end c==1
 else { # c != 1
  bar<-numeric(0)
  sizes<-numeric(0)
  for (i in 1:c) {
   bar<-rbind(bar,apply(scoremat[g==levels(g)[i],,drop=F],2,mean))
   sizes<-c(sizes,sum(g==levels(g)[i]))
  }
  STATISTIC<-p*sum(sizes*(norm(bar)^2))/switch(score,"sign"=1,"rank"=c2)
 }

 if (all(cond,score=="sign"))
 {
  Qd<-numeric(0)
  if(c==1) {
   for (i in 1:cond.n) {
    d<-sample(c(-1,1),n,replace=T)
    Qd<-c(Qd,n*p*sum(apply(sweep(scoremat,1,d,"*"),2,mean)^2))
   }
  }
  else {
   for (i in 1:cond.n) {
    gd<-sample(g)
    bar<-numeric(0)
    for (j in 1:c)
     bar<-rbind(bar,apply(scoremat[gd==levels(gd)[j],,drop=F],2,mean))
    Qd<-c(Qd,p*sum(sizes*(norm(bar)^2)))
   }
  }
  PARAMETER<-cond.n
  names(PARAMETER)<-"replications"
  PVAL<-mean(Qd>=STATISTIC)
 }
 else  
 {  
  PVAL<-1-pchisq(STATISTIC,(df<-p*max(1,c-1)))
  PARAMETER<-df
  names(PARAMETER)<-"df"
 }
 ALTERNATIVE<-"two.sided"
 names(STATISTIC)<-"Q.2"
 res<-c(list(statistic=STATISTIC,parameter=PARAMETER,p.value=PVAL,null.value=NVAL,alternative=ALTERNATIVE,method=METHOD,data.name=DNAME))
 class(res)<-"htest"
 return(res)
}

