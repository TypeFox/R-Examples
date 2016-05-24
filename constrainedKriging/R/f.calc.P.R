f.calc.P <- function(t.bb.covmat, t.pred.designmat, t.cov.beta.coef)
### purpose: calculate the P matrix = t.bb.covmat - var[t(t.Xm) %*% t.beat.gls]
###          of the covariance-matching constrained kriging predictor
###
###          function returns a list with three elements
###              i) P.posdef = logial element if TRUE then matrix P is positive definite
###              ii) P = m x m matrix m = number of polygons
###              iii) P1 = P^0.5 m x m matrix
### arguments:
###            t.bb.covmat = m x m covariance matrix of between the predicton location
###            t.Xm = model matrix of the prediction location
###            t.var.beta.gls = ( t(X) %*% inv.Sigma %*% X)^(-1)
###            t.eps = is the tuning parameter. It will replace those eigen values
###                     when the positive definiteness of  P is not satisfied.
###
### author: Ch. Hofer
### date: 27.3.2007
{

t.P <- t.bb.covmat - tcrossprod( crossprod( t( t.pred.designmat ) , t.cov.beta.coef ), t.pred.designmat )

t.P.eig <- eigen(t.P,symmetric = T)

if( min( t.P.eig$values ) >= 0 )
{
 t.P1 <- crossprod( t( t.P.eig$vectors ), tcrossprod( sqrt( diag( t.P.eig$values, nrow=length( t.P.eig$values ) ) ), t.P.eig$vector ) )
}
else
{
t.P1 <- t.P
t.P1[,] <- NaN
}
return( list( P =  t.P, P1 = t.P1 ) )
} ## end function f.calc.P
