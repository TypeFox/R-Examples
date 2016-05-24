"vif" <-
function(obj, digits=5){
Qr <- obj$qr
if (is.null(obj$terms) || is.null(Qr)) 
        stop("invalid 'lm' object:  no terms or qr component")
tt <- terms(obj)
hasintercept <- attr(tt, "intercept") > 0
p <- Qr$rank
if(hasintercept) p1 <- 2:p else p1 <- 1:p
R <- Qr$qr[p1,p1, drop=FALSE]
if(length(p1)>1) R[row(R)>col(R)] <- 0
Rinv <- qr.solve(R)
vv <- apply(Rinv, 1, function(x)sum(x^2))
ss <- apply(R, 2, function(x)sum(x^2))
vif <- ss*vv
signif(vif, digits)
}
