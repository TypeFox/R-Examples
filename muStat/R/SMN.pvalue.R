`SMN.pvalue` <-
function(pP,qP, pX=0,qX=0, pQ=0,qQ=0, wP=.25,wQ=.25, exact = NULL)
{
  M <- function(P,X,Q, wP,wQ, e = 1, op = "+", f = function(x) x) {
    O3 <- function(X,Y,Z,Op) matrix(outer(outer(X,Y,Op),Z,Op))
    O3(wP^e*f(P), (1-wP-wQ)^e*f(X), wQ^e*f(Q), op)
  }

  exact <-
    (dP<-pP-qP)+(dX<-pX-qX)+(dQ<-pQ-qQ)<100 && is.null(exact) || exact
  if (((nP<-pP+qP)*(nX<-pX+qX)*(nQ<-pQ+qQ)>10^6) || !exact )
    return( 1-pchisq(
      M(dP,dX,dQ, wP,wQ)^2/              # Eq. (1) in Wittkowski (2002)
      M(nP,nX,nQ, wP,wQ, 2),1)[1])       # Eq. (2) in Wittkowski (2002)
  else {
    tb <- cbind(
      M(nP,nX,nQ, wP,wP, 0,"*", function(n) mu.dbinom(0:n, n, .5)),
      M(nP,nX,nQ, wP,wQ, 1,"+", function(n) (0:n)-(n:0) )^2)
    return(1-sum(tb[tb[,2]<c(M(dP,dX,dQ, wP,wQ)^2),1]))
  }
}
