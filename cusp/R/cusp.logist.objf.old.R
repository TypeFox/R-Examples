`cusp.logist.objf.old` <-
function(p, xa, xb, y, alpha = as.matrix(xa) %*% p[1:NCOL(xa)], 
  beta = as.matrix(xb) %*% if(NCOL(xb)>1) {c(1,p[2:NCOL(xb)-1+NCOL(xa)])} else {c(1)}){
    if(length(p)!=NCOL(xa)+NCOL(xb)-1) {
        stop('p should have length ncol(xa)+ncol(xb)-1') 
    }
    X = cbind(1, (alpha/beta^2 > -50.0) / (1 + exp(-alpha/beta^2)))
    yhat <- qr.fitted(qr(X), y)
    rss <- sum((y - yhat)^2)
    attr(rss,'X') = X
    attr(rss,'Y') = y
    attr(rss,'alpha') = alpha
    attr(rss,'beta' ) = beta
    attr(rss,'RSq' ) = 1 - var(yhat)/var(y)
    rss
}

