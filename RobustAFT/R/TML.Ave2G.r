TML.Ave2G <-
function(X,y,delta, sigma,sigma.t, mui,mui.t, wgt,cu){ cl <- -cu
n <- length(y); p <- ncol(X); ku <- kc <- 0; nc <- sum(delta==0); zero <- 1e-6
indu  <- (1:n)[delta==1]; indc  <- (1:n)[delta==0]
rs    <- (y-mui)/sigma
if (nc < n) {
 ki  <- wgt[indu]*rs[indu]^2
 ku  <- sum(ki)}
if (nc > 0) {
 muic <- mui[indc]
 muit <- mui.t[indc]
 rci  <- rs[indc]
 ai <- pmax( rci, (sigma.t*cl - muic + muit )/sigma   )
 bi <-            (sigma.t*cu - muic + muit )/sigma
 den <- 1-pnorm(rci)
 ok  <- den > zero
 Iki <- ai*dnorm(ai) - bi*dnorm(bi) + pnorm(bi) - pnorm(ai)
 eki <- rep(0, nc)
 eki[bi > ai & ok] <- Iki[bi > ai & ok]/den[bi > ai & ok]
 kc <- sum(eki)}
(ku + kc)/(n - p)}

