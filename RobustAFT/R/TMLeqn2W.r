TMLeqn2W <-
function(d.sigma,d.Beta,rs0,delta,X,cl,cu,const) {
# equation left-hand sides for delta-sigma using a rectangular weight function
p   <- ncol(X); n <- length(rs0); ku <- kc <- 0; zero <- 1e-6
iu  <- delta==1; ic  <- delta==0; nc  <- sum(ic)
if (nc < n ) {
 Xu  <- X[iu,]; r0u <- rs0[iu]
 ru  <- (r0u-Xu%*%d.Beta)/d.sigma
 ki  <- ww(r0u,cl,cu)*(exp(ru)-1)*ru
 ku  <- sum(ki) }
if (nc > 0) {
 Xc  <- X[ic,]; r0c <- rs0[ic]
 rc  <- (r0c-Xc%*%d.Beta)/d.sigma
 Fo  <- plweibul(rc)
 ok  <- (1-Fo) > zero
 ai  <- (pmax(r0c,cl)- Xc%*%d.Beta)/d.sigma
 bi  <- (         cu - Xc%*%d.Beta)/d.sigma
 Iki <- ai*dlweibul(ai)-bi*dlweibul(bi)+plweibul(bi)-plweibul(ai)
 eki <- rep(0,nc)
 eki[bi>ai & ok]  <- Iki[bi>ai & ok]/(1-Fo[bi>ai & ok])
 kc  <- sum(eki)}
(ku+kc)/(n-p)-const}

