`sr.indep.test` <- function(X,Y=NULL,g=NULL,score=c("sign","symmsign","rank"),regexp=FALSE,cond=FALSE,cond.n=1000,na.action=na.fail)
{
 if(all(is.null(Y),is.null(g))) stop("Y or g must be given")

 score=match.arg(score)
 
 if(!is.null(g)) {
  if(is.factor(g)) {
   DNAME<-paste(deparse(substitute(X)),"by",deparse(substitute(g)))
   X<-na.action(X)
   Y<-X[,g==levels(g)[2]]
   X<-X[,g==levels(g)[1]]
  }
  else {
   if(is.character(g)) 
    if (regexp) gn<-grep(g,colnames(X))
    else gn<-match(g,colnames(X))
   else gn<-g    # should be is.numeric(g)==TRUE
   gn.char<-paste("c(",paste(g,collapse=","),")",sep="")
   DNAME<-paste(deparse(substitute(X)),"columns",gn.char,"vs. the rest")
   X<-na.action(X)
   Y<-X[,-gn]
   X<-X[,gn]
  }
 }
 else
  {
   DNAME<-paste(deparse(substitute(X)),"and",deparse(substitute(Y)))
   X<-na.action(X)
   if(!is.null(attr(X,"na.action"))) Y<-Y[-(attr(X,"na.action")),]
   Y<-na.action(Y)
   if(!is.null(attr(Y,"na.action"))) X<-X[-(attr(Y,"na.action")),]
  }      

 X<-as.matrix(X)
 Y<-as.matrix(Y) 
 p1<-dim(X)[2]
 p2<-dim(Y)[2]

 n<-dim(X)[1]
 if(dim(Y)[1]!=n) stop("the number of observations in the data sets differ")
 STATISTIC<-switch(score,
	"sign"=
	{
         METHOD<-"Multivariate independence test using spatial signs"
	 SX<-spatial.signs(X,center=T,shape=T)
	 SY<-spatial.signs(Y,center=T,shape=T)
	 ave<-t(SX)%*%SY/n
	 n*p1*p2*mat.norm(ave)^2
	},
	"symmsign"=
	{
         METHOD<-"Multivariate independence test using spatial symmetrized signs"
	 SX<-spatial.signs(pairdiff(X),center=F,shape=T)
	 SY<-spatial.signs(pairdiff(Y),center=F,shape=T)
	 # there are now n*(n-1)/2 rows in these matrices
	 m<-choose(n,2)
	 RX<-spatial.rank(X,shape=T)
	 RY<-spatial.rank(Y,shape=T)
	 cx1<-mean(norm(RX)^2)
	 cx2<-mean(norm(RY)^2)
	 ave<-t(SX)%*%SY/m
	 ((n*p1*p2)/(4*cx1*cx2))*mat.norm(ave)^2
	},
	"rank"=
	{
         METHOD<-"Multivariate independence test using spatial ranks"
	 RX<-spatial.rank(X,shape=T)
	 RY<-spatial.rank(Y,shape=T)
	 ave<-t(RX)%*%RY/n
	 if (p1==1) cx1<-1/3 
	  else cx1<-mean(norm(RX)^2)
	 if (p2==1) cx2<-1/3
	  else cx2<-mean(norm(RY)^2)
	 ((n*p1*p2)/(cx1*cx2))*mat.norm(ave)^2
	})
 names(STATISTIC)<-"Q.2"

 if (cond)
 {
  PARAMETER<-cond.n
  names(PARAMETER)<-"replications"
  statg<-NULL
  for (j in 1:cond.n)
  {
   statg<-c(statg,switch(score,
	"sign"=
	{
	 perm<-sample(n)
	 ave<-t(SX)%*%SY[perm,]/n
   	 n*p1*p2*mat.norm(ave)^2
	},
	"symmsign"=
	{
	 perm<-sample(m)
	 ave<-t(SX)%*%SY[perm,]/m
   	 ((n*p1*p2)/(4*cx1*cx2))*mat.norm(ave)^2
	},
	"rank"=
	{
	 perm<-sample(n)
	 ave<-t(RX)%*%RY[perm,]/n
	 ((n*p1*p2)/(cx1*cx2))*mat.norm(ave)^2
	}))
  }
  PVAL<-mean(statg>=STATISTIC)
 }
 else
 {
  df<-p1*p2
  PARAMETER<-df
  names(PARAMETER)<-"df"
  PVAL<-1-pchisq(STATISTIC,df)
 }
 ALTERNATIVE<-"two.sided"
 NVAL<-0
 names(NVAL)<-"measure of dependence"
 res<-list(statistic=STATISTIC,parameter=PARAMETER,null.value=NVAL,alternative=ALTERNATIVE,method=METHOD,data.name=DNAME,p.value=PVAL)
 class(res)<-"htest"
 res
}

