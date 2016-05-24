.eigenConst <- function(Sev, p, lambda=12){
	moptv <- 999999
	if( Sev[p] > 0 & Sev[1]/Sev[p] < lambda){
		moptv <- sum(log(Sev) + 1)
		return(Sev)
	} else{
		for(i in 1:p){
			cm <- Sev[i]
			if( cm > 0){
				m <- cm/lambda
				Sevm <- pmin( pmax( Sev, m ), cm)
				moptv_tmp <- sum(log(Sevm) + Sev/Sevm)
				if( moptv_tmp < moptv ){
					moptv <- moptv_tmp
					mopt <- m
				}
			}
			m <- Sev[p-i+1]
			if( m > 0){
				cm <- m*lambda
				Sevm <- pmin( pmax( Sev, m ), cm)
				moptv_tmp <- sum(log(Sevm) + Sev/Sevm)
				if( moptv_tmp < moptv ){
					moptv <- moptv_tmp
					mopt <- m
				}
			}
		}
		Sevm <- pmin( pmax( Sev, mopt ), mopt*lambda)

		return( Sevm )
	}
}


sclust<- function(X,k,V,R,restr.fact=12,tol=1e-4,maxiters=100,maxiters.S=1000, print.it=FALSE) {
	if( missing(X) ) stop("'X' missing")
	if( missing(k) ) stop("'k' missing")
	if( missing(V) ) stop("'V' missing")
	if( missing(R) ) stop("'R' missing")
	
	if(is.data.frame(X) | is.matrix(X))
		X <- data.matrix(X)
	else stop("Data matrix must be of class matrix or data.frame")
	n <- nrow(X)
	p <- ncol(X)
	
	if(is.data.frame(V) | is.matrix(V))
		V <- data.matrix(V)
	if( any(dim(V) != dim(X)) ) stop("'X' and 'V' have non-conforming size")
	epsilon=sum(V==0)/(n*p)

	if( length(unique(R)) != k) stop("Number of cluster labels must be 'k'")
	
  ## init ##

m=matrix(NA,k,p)
Sigma=array(0,c(k,p,p))
D=matrix(NA,n,p)
Dd=Ddtmp=matrix(NA,n,k)
autovalues=matrix(NA,p,k)
U=Sigma
CSmat=array(NA,c(k,n,p))
det=rep(NA,k)
  
  ## init values ##

  pi=prop.table(table(R[R!=0]))

Xt=X
Xt[V==0]=NA 
  for(j in 1:k) {
    m[j,]=apply(Xt[R==j,],2,mean,na.rm=T)
    Sigma[j,,]=var(Xt[R==j,],na.rm=T)
    s=eigen(var(Xt[R==j,],na.rm=T))
    U[j,,]=s$vectors
    autovalues[,j]=s$values
}

autovalues <- matrix(.eigenConst(as.vector( autovalues), p, restr.fact),ncol=k)

    for(j in 1:k) {
      Sigma[j,,]=U[j,,]%*%diag(autovalues[,j])%*%t(U[j,,])
      det[j]=prod(autovalues[,j])
    }

for(j in 1:k) {

Dd[,j]=log(pi[j])+ldmvnorm(Xt,m[j,],Sigma[j,,])

}

Dd=exp(Dd-apply(Dd,1,sumlog))

 lik=0
    for(j in 1:k) {
      lik=lik+sum(Dd[,j]*(log(pi[j])+ldmvnorm(Xt,m[j,],Sigma[j,,])))}

      likold=lik-2*tol
ii=0 
  while(lik-likold>tol & ii < maxiters) {
ii=ii+1

## CES step

iter=0
flag=FALSE
while(iter < maxiters.S & flag==FALSE) {
iter=iter+1
s1=sample(which(V==1),1)
s2=sample(which(V==0),1)

Vc=V
Vc[s1]=0
Vc[s2]=1

Xt=X
Xt[Vc==0]=NA
likcand=0
    for(j in 1:k) {
      likcand=likcand+sum(Dd[,j]*(log(pi[j])+ldmvnorm(Xt,m[j,],Sigma[j,,])))}

if(likcand>lik) {
V=Vc
flag=TRUE
}
}

for(j in 1:k) {

Dd[,j]=log(pi[j])+ldmvnorm(Xt,m[j,],Sigma[j,,])

}

R=apply(Dd,1,which.max)
R[which(apply(V==0,1,all))]=0
Dd[which(apply(V==0,1,all)),]=0

## M step

pi=apply(Dd[which(apply(V!=0,1,any)),],2,sumlog)
pi=exp(pi-sumlog(pi))

Dd=exp(Dd-apply(Dd,1,sumlog))

                Xt <- X
		Xt[V==0] <- NA
	
 for(j in 1:k) {
	XtDd <- sweep(Xt, 1, Dd[,j], "*") 
		VDd <- sweep(V, 1, Dd[,j], "*")
		m[j,] <- colSums(XtDd, na.rm=T)
		m[j,] <- m[j,]/colSums(VDd, na.rm=T)
                Stmp <- matrix(NA, p,p)
		for(h in 1:(p-1)) {
			for(l in (h+1):p) {
	    			Stmp[h,l] <- sum(Dd[,j]*(Xt[,h]-m[j,h])*(Xt[,l]-m[j,l]),na.rm=T)
	    			Stmp[h,l] <- Stmp[h,l]/sum(Dd[,j]*V[,h]*V[,l])
				Stmp[l,h] <- Stmp[h,l]
			}
		}
		for(h in 1:p) Stmp[h,h] <- sum(Dd[,j]*(Xt[,h]-m[j,h])^2,na.rm=T)/sum(Dd[,j]*V[,h])
    s=eigen(Stmp)
    U[j,,]=s$vectors
    autovalues[,j]=s$values
}


autovalues <- matrix(.eigenConst(as.vector( autovalues), p, restr.fact),ncol=k)

    for(j in 1:k) {
      Sigma[j,,]=U[j,,]%*%diag(autovalues[,j])%*%t(U[j,,])
            det[j]=prod(autovalues[,j])
    }

    likold=lik
   lik=0
    for(j in 1:k) {
      lik=lik+sum(Dd[,j]*(log(pi[j])+ldmvnorm(Xt,m[j,],Sigma[j,,])))}
   
if( print.it )cat("iter", ii, "; current lik:", lik, "; change in lik:",lik-likold, "\n")
 
  }

return(list(R=R,pi=pi,mu=m,S=Sigma,V=V,lik=lik,iter=ii))}


