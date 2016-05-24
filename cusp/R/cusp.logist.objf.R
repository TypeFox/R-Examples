`cusp.logist.objf` <-
function(p, xa, xb, y, alpha = as.matrix(xa) %*% p[1:NCOL(xa)], 
  beta = as.matrix(xb) %*% if(NCOL(xb)>1) {c(1,p[2:NCOL(xb)-1+NCOL(xa)])} else {c(1)}){
    if(sum(diag(var(xa)))<.Machine$double.eps^0.5 && sum(diag(var(xa)))<.Machine$double.eps^0.5 ){
        stop("No variation in predictors")
    }
    X = cbind(1, (alpha/beta^2 > -50.0) / (1 + exp(-alpha/beta^2)))
    Y = if(is.qr(y)) {y} else {qr(cbind(1,y))}
    ssrss <- cusp.subspacerss(predictors = X, dependents = Y); 
#cat(var(X[,2]), p, ssrss$cor, '\n')
    idx <- c(which(zapsmall(ssrss$cor)<1),1)[1] # index of first non-perfect correlation or 1
    rss <- ssrss$rss[idx]
    attr(rss,'X') = X
    attr(rss,'Y') = if(is.qr(y)) { qr.X(y) } else { y }
    attr(rss,'alpha') = alpha
    attr(rss,'beta' ) = beta
    attr(rss,'RSq' ) = ssrss$cor[idx]^2 # c(which(ssrss$cor<1),1)[1]]^2
    rss
}

