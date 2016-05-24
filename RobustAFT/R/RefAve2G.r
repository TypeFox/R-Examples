RefAve2G <-
function(sigma,Beta,X,y,delta) {
# auxliary functionfor RefSigma()
n <- length(y); p <- ncol(X); xk <- 1.5477
ku <- kc <- 0; nu <- sum(delta)
r  <- as.vector((y-X%*%as.matrix(Beta))/sigma)
if (nu > 0) {
 ru <- r[delta==1]
 ku <- sum(ChiSG(ru)) }
if (nu < n) {
 rc <- r[delta==0]
 kc <- sum(unlist(lapply(rc,IChidnorm,k=xk))) }
(ku+kc)/(n-p)}

