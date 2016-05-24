RefAve2W <-
function(sigma,Beta,X,y,delta) {
# auxliary functionfor RefSigmaW()
n <- length(y); p <- ncol(X); xk <- 1.718; c0 <- -0.1351788
ku <- kc <- 0; nu <- sum(delta)
r  <- as.vector((y-X%*%as.matrix(Beta))/sigma)
if (nu > 0) {
 ru <- r[delta==1]
 ku <- sum(ChiSw(ru-c0)) }
if (nu < n) {
 rc <- r[delta==0]
 kc <- sum(unlist(lapply(rc,IChidlweibul,k=xk))) }
(ku+kc)/(n-p)}

