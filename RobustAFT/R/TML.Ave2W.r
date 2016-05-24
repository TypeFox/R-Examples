TML.Ave2W <-
function(X,y,delta,sigma,sigma.t, mui, mui.t,wgt, cl,cu,const){
n <- length(y); p <- ncol(X); ku <- kc <- 0; zero <- 1.e-6
indu <- (1:n)[delta==1]; indc <- (1:n)[delta==0]; nc <- sum(delta==0)
rs <- (y-mui)/sigma
if (nc < n) {
  ki <- wgt[indu]* (exp(rs[indu])-1)*rs[indu]
  ku <- sum(ki) }
if (nc > 0) {
  muic <- mui[indc]
  muit <- mui.t[indc]
  rci  <- rs[indc]
  ai   <- pmax( rci, (sigma.t*cl-muic+muit)/sigma )
  bi   <-            (sigma.t*cu-muic+muit)/sigma 
  Fo   <- plweibul(rci)
  ok   <- (1-Fo) > zero 
  Iki <- ai*dlweibul(ai)-bi*dlweibul(bi)+plweibul(bi)-plweibul(ai)
  eki  <- rep(0,nc)
  eki[bi > ai & ok] <- Iki[bi > ai & ok]/(1-Fo[bi > ai & ok])
  kc <- sum(eki) }
(ku+kc)/(n-p)}

