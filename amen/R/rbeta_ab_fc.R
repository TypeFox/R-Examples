#' Gibbs sampling of additive row and column effects and regression coefficient
#' 
#' Simulates from the joint full conditional distribution of (a,b,beta)
#' 
#' 
#' @usage rbeta_ab_fc(Z, Sab, rho, X, mX, mXt, XX, XXt, Xr, Xc, s2 = 1)
#' @param Z n X n (latent) normal relational matrix, with multiplicative
#' effects subtracted out
#' @param Sab row and column covariance
#' @param rho dyadic correlation
#' @param X n x n x p covariate array
#' @param mX design matrix (matricizied version of X)
#' @param mXt dyad-transposed design matrix
#' @param XX regression sums of squares
#' @param XXt crossproduct sums of squares
#' @param Xr row sums for X
#' @param Xc column sums for X
#' @param s2 dyadic variance
#' @return \item{beta}{regression coefficients} \item{a}{additive row effects}
#' \item{b}{additive column effects}
#' @author Peter Hoff
#' @export rbeta_ab_fc
rbeta_ab_fc <-
function(Z,Sab,rho,X,mX,mXt,XX,XXt,Xr,Xc,s2=1) 
{
  ### 
  p<-dim(X)[3] 
  Se<-matrix(c(1,rho,rho,1),2,2)*s2
  iSe2<-mhalf(solve(Se))
  td<-iSe2[1,1] ; to<-iSe2[1,2]
  Sabs<-iSe2%*%Sab%*%iSe2
  tmp<-eigen(Sabs)
  k<-sum(zapsmall(tmp$val)>0 )
  ###

  ###
  mXs<-td*mX+to*mXt                  # matricized transformed X
  XXs<-(to^2+td^2)*XX + 2*to*td*XXt  # sum of squares for transformed X
  Zs<-td*Z+to*t(Z)
  zr<-rowSums(Zs) ; zc<-colSums(Zs) ; zs<-sum(zc) ; n<-length(zr) 
  ###

  ## dyadic and prior contributions  
  if(p>0)
  {
    lb<- crossprod(mXs,c(Zs)) 
    Qb<- XXs + XX/nrow(mXs) 
  }
  ##

  ## row and column reduction
  ab<-matrix(0,nrow(Z),2)
  if(k>0)
  {
    n<-nrow(Z)
    G<-tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
    K<-matrix(c(0,1,1,0),2,2)
    A<-n*t(G)%*%G + diag(k)
    B<-t(G)%*%K%*%G
    iA0<-solve(A)
    C0<- -solve(A+ n*B)%*%B%*%iA0

    iA<-G%*%iA0%*%t(G)
    C<-G%*%C0%*%t(G)

    if(p>0) 
    {    
    Xsr<-td*Xr + to*Xc  # row sums for transformed X
    Xsc<-td*Xc + to*Xr  # col sums for transformed X
    Xss<-colSums(Xsc)  
    lb<- lb - (iA[1,1]*crossprod(Xsr,zr) + iA[2,2]*crossprod(Xsc,zc) +
               iA[1,2]*(crossprod(Xsr,zc) + crossprod(Xsc,zr)) +
               sum(C)*Xss*zs )

    tmp<-crossprod(Xsr,Xsc)
    Qb<- Qb - (iA[1,1]*crossprod(Xsr,Xsr) + iA[2,2]*crossprod(Xsc,Xsc) +
               iA[2,1]*(tmp+t(tmp)) +  sum(C)*Xss%*%t(Xss) ) 
    }
  }
  ##
  if(dim(X)[3]==0){     beta<-numeric(0) }
 
  if(p>0) 
  { 
  V<-solve(Qb)
  m<-V%*%(lb)
  beta<-c(rmvnorm(1,m,V)) 
  }
  ####
 
  #### simulate a, b 
  if(k>0) 
  {
    E<- Zs-Xbeta(td*X+to*aperm(X,c(2,1,3)),beta)
    er<-rowSums(E) ; ec<-colSums(E) ; es<-sum(ec) ; n<-length(er) 
    m<-t(t(crossprod(rbind(er,ec),t(iA0%*%t(G)))) + rowSums(es*C0%*%t(G)) )
    hiA0<-mhalf(iA0)
    e<-matrix(rnorm(n*k),n,k) 
    w<-m+ t( t(e%*%hiA0) - c(((hiA0-mhalf(iA0+n*C0))/n)%*% colSums(e) ) )
    ab<- w%*%t(G)%*%solve(iSe2) 
  }

list(beta=beta,a=ab[,1],b=ab[,2] )  
}
