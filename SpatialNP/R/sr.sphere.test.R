`sr.sphere.test` <-
function(X, score=c("sign","symmsign"), shape=NULL, na.action=na.fail)
{
 X<-na.action(X)
 DNAME=deparse(substitute(X))
 
 score=match.arg(score)
 X<-as.matrix(X)

 p<-dim(X)[2]
 if (p<2) stop("'X' must be at least bivariate")
 n<-dim(X)[1]
 
 if(!is.null(shape))
 {
  if(!is.matrix(shape)) stop("'shape' must be a matrix")
  if(!all(dim(shape)==c(p,p))) stop("dimensions of 'shape' and 'X' do not match")
  X<-X%*%solve(mat.sqrt(shape))
 }

 Cp<-Cpp(p)

 STATISTIC<-switch(score,
  "sign"=
  {
   METHOD="Test of sphericity using spatial signs" 
   S<-spatial.signs(X,F,F)
   S1<-as.vector(t(S)%*%S/n)
   Q1<-(sum((Cp%*%S1)^2))
   gamma<-2/(p*(p+2))
   n*Q1/gamma
  },
  "symmsign"=
  {
   METHOD="Test of sphericity using spatial symmetrized signs" 
   tmp<-Q2internal(X)
   S2<-tmp[1:p^2]
   covmat<-4*(matrix(tmp[-(1:p^2)],ncol=p^2)-S2%*%t(S2))/n
   as.vector(t(Cp%*%S2)%*%gen.inv(covmat)%*%(Cp%*%S2))
  })
 names(STATISTIC)<-"Q.2"
 NVAL<-paste("diag(",paste(p),")",sep="")
 names(NVAL)<-"shape"
 PVAL<-1-pchisq(STATISTIC,(df<-(p+2)*(p-1)/2))
 PARAMETER<-df
 names(PARAMETER)<-"df"
 ALTERNATIVE<-"two.sided"
 res<-list(statistic=STATISTIC,parameter=PARAMETER,p.value=PVAL,null.value=NVAL,alternative=ALTERNATIVE,method=METHOD,data.name=DNAME)
 class(res)<-"htest"
 res
}

