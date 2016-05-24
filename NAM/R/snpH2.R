snpH2=function(gen){
  anyNA = function(x) any(is.na(x))
if(any(is.na(gen))) stop("Missing values not allowed")
if(is.numeric(gen)!=T) stop("Object gen must to be a numeric matrix")

snps = ncol(gen)
obs = nrow(gen)
x = matrix(1:obs,ncol=1)

K = tcrossprod(gen)
K = K/mean(diag(K))
  
# Defining log-REML
loglike<-function(theta){
  lambda<-exp(theta)
  logdt<-sum(log(lambda*delta+1))
  h<-1/(lambda*delta+1)
  yy<-sum(yu*h*yu)
  yx<-matrix(0,q,1)
  xx<-matrix(0,q,q)
  for(i in 1:q){
    yx[i]<-sum(yu*h*xu[,i])
    for(j in 1:q){xx[i,j]=sum(xu[,i]*h*xu[,j])}}
  loglike = -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
  return(-loglike)}
fixed<-function(lambda){
  h<-1/(lambda*delta+1)
  yy<-sum(yu*h*yu)
  yx=timesVec(yu,h,xu,q)
  xx=timesMatrix(xu,h,xu,q,q)
  beta<-qr.solve(xx,yx)
  sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
  sigma2 = as.numeric(sigma2)
  var<-diag((chol2inv(xx))*sigma2)
  stderr<-sqrt(var)
  return(c(beta,stderr,sigma2))}

# Eigendecomposition of K
qq<-eigen(as.matrix(K),symmetric=T)
delta<-qq[[1]]
uu<-qq[[2]]
q<-1
n<-ncol(K)

H2 = c()
theta<-0


pb=txtProgressBar(style=3)
for(i in 1:snps){
  
  y = gen[,i]
  yu<-t(uu)%*%y
  xu<-t(uu)%*%x
  vp<-var(y)
  # Finding lambda through optimization
  parm<-optim(par=theta,fn=loglike,method="L-BFGS-B",lower=0,upper=10)
  lambda<-exp(parm$par)
  # Results
  parmfix<-fixed(lambda)
  Ve<-parmfix[2*q+1]
  Vg<-lambda*Ve
  h2=Vg/(Vg+Ve)
  
  # Saving loop
  H2 = c(H2,h2)
  setTxtProgressBar(pb,i/snps)
};close(pb)

names(H2)=colnames(gen)
class(H2)="H2"
return(H2)
}

# Plot

plot.H2 = function(x,...,chr=NULL){
  anyNA = function(x) any(is.na(x))
  her = as.numeric(x)
  plot(her,xlab="Genome",ylab="Heritability",main="Gene content")
  rect(par("usr")[1],par("usr")[3],
       par("usr")[2],par("usr")[4],
       col=rgb(0.25,0.25,0.25,0.5))
  lines(her,lwd=4,pch=20,type="p",col="blue")
  lines(her,lwd=1,pch=20,type="p",col="lightblue")
  if(is.null(chr)!=T) abline(v=cumsum(chr[-length(chr)])+0.5,lty=3)
}
