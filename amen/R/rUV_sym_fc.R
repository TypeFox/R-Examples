#' Gibbs sampling of U and V
#' 
#' A Gibbs sampler for updating the multiplicative effect matrices U and V
#' in the symmetric case. In this case \code{U\%*\%t(V)} is symmetric, so
#' this is parameterized as \code{V=U\%*\%L} where \code{L} is the 
#' diagonal matrix of eigenvalues of \code{U\%*\%t(V)}. 
#' 
#' @usage rUV_sym_fc(E, U, V, s2 = 1, shrink=TRUE)
#' @param E square residual relational matrix
#' @param U current value of U
#' @param V current value of V
#' @param s2 dyadic variance
#' @param shrink adaptively shrink the factors with a hierarchical prior
#' @return \item{U}{a new value of U} \item{V}{a new value of V}
#' @author Peter Hoff
#' @examples
#' 
#' U0<-matrix(rnorm(30,2),30,2) ; V0<-U0%*%diag(c(3,-2)) 
#' E<- U0%*%t(V0) + matrix(rnorm(30^2),30,30) 
#  rUV_sym_fc(E,U0,V0) 
#' rUV_sym_fc 
#' 
#' @export rUV_sym_fc
rUV_sym_fc<-function(E,U,V,s2=1,shrink=TRUE)
{

  R<-ncol(U) ; n<-nrow(U) 
  L<-diag( (V[1,]/U[1,]) ,nrow=R) 
  L[is.na(L)]<-1  # to handle zero start vals

  ## update inverse variance of U
  if(shrink){ivU<-diag( rgamma(R, (2+n)/2 ,(1+apply(U^2,2,sum))/2) ,nrow=R )}
  if(!shrink){ivU<-diag(1/n,nrow=R) }

  ## update each U[i,]
  for(i in rep(sample(1:n),4))
  {
    l<-L%*%( apply(U*E[i,],2,sum) -  U[i,]*E[i,i] )/s2
    iQ<- solve( ( ivU +    L%*%( crossprod(U) - U[i,]%*%t(U[i,]) )%*%L/s2 ) )
    U[i,]<- iQ%*%l + t(chol(iQ))%*%rnorm(R) 
  }

  ## consider MH update - add in later verision

  ## update "eigenvalues"
  for(r in 1:R)
  {
    Er<-E-U[,-r,drop=FALSE]%*%L[-r,-r,drop=FALSE]%*%t(U[,-r,drop=FALSE]) 
    l<- sum( ( Er*(U[,r]%*%t(U[,r])))[upper.tri(Er)] )/s2
    iq<- 1/(1+sum( ( (U[,r]%*%t(U[,r]))^2 )[upper.tri(Er)] )/s2 )
    L[r,r]<-rnorm(1, iq*l,sqrt(iq) ) 
  }

  list(U=U,V=U%*%L) 
}







