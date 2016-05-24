f.calc.Q <- function(t.skvar, t.pred.covmat.ichol.trans, t.orth.designmat, t.cov.beta.coef)
### purpose: calculate the Q matrix = t.uk.var - var[t(t.Xm) %*% t.beat.gls]
###          of the covariance-matching constrained kriging predictor
###
###          function returns a list with four elements
###              i)   Q.posdef = logial element if TRUE then matrix Q is positive definite
###              ii)  Q = m x m matrix m = number of polygons
###              iii) Q1 = P^0.5 m x m matrix
###              iv)  Q1.inv = inverse matrix of Q1
### arguments:
###            t.SK.var = simple kiriging (co)varianz
###            t.C = n x m matrix with the covariance between the n sampling points and the m
###                  predicton location
###            t.cov.beta.coef = ( t(X) %*% inv.Sigma %*% X)^(-1)
###            t.isupport.covmat  = n x n matrix inverse of the covariance matrix of the n sampling points
###            t.X = design matrix of the samples (# sample points x # beta parameters)-matrix
###            t.eps = is the tuning parameter. It will replace those eigen values
###                     when the positive definiteness of Q is not satisfied.
###
### author: Ch. Hofer
### date: 27.3.2007
{
# # # Q = Var[UK predictor] - VAR[ fitted values]
# # # C' SIGMA^-1 C - C' SIGMA^1 X (X' SIGMA^1 X )^-1 X' SIGMA^-1 C
# # # t.Q = t.SK.var - t(t.C) %*% t.isupport.covmat %*% t.X %*% t.cov.beta.coef %*% t(t.X) %*% t.isupport.covmat %*% t.C
t.aux.I <- crossprod( t( t.pred.covmat.ichol.trans), t.orth.designmat)
t.Q <- t.skvar - crossprod( t( crossprod( t( t.aux.I ), t.cov.beta.coef ) ), t( t.aux.I ) )

### Eigenwerte der Matrix Q
t.Q.eig    <- eigen(t.Q, symmetric=T)

# # # kontrolle ob die Eigenwerte gr?sser als Null sind
if( min( t.Q.eig$values ) > 0 )
{
    t.sqrt.diag.Qeig <- diag( sqrt( t.Q.eig$values ), nrow=length( t.Q.eig$values ) )
    
    t.Q1 <- crossprod( t( t.Q.eig$vectors ), tcrossprod( t( t.sqrt.diag.Qeig ), t.Q.eig$vectors))
    
    t.chol.Q1 <- chol(t.Q1)

t.ichol.Q1 <- forwardsolve(
    t( t.chol.Q1 ),
    diag( nrow( t.Q1 ) )
	)
t.iQ1 <- crossprod( t.ichol.Q1, t.ichol.Q1 )


}
else
{
    t.Q1 <- t.Q
    t.Q1[,] <- 0
    t.iQ1 <- t.Q1
    
}


# Inversere Choleskymatrix L^-1

### Then Q1 inverse.
return( list( Q =  t.Q, Q1 = t.Q1, iQ1 = t.iQ1 ) )
} ## end of function f.calc.Q 
