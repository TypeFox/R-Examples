"adaptw" <-
function(rso,cl,cu,option){
# adaptive cut-off values
alpha <- NA; n <- length(rso)
if (option=="adaptive"){
 rhos  <- sort(rhow(rso,const=0)); rhocu <- rhow(cu,const=0)
 Po    <- apply(as.matrix(rhos),1,F0w,tol=0.0001,maxit=150) # F0w is a Fortran version of F0w.s
 dsr   <- Discr(rhos,Po,rhocu)
 alpha <- dsr$r.min
 tu2 <- as.numeric(quantile(rhos,probs=alpha)); logtu2 <- log(tu2); upp <- log(tu2+1.2*logtu2) 
 if (tu2 <=1 )             {tl <- tu <- 0}
 if (tu2>1   & tu2 <=1.5)  tl <- uniroot(rhow,lower=-tu2,  upper=0,       const=tu2)$root
 if (tu2>1.5 & tu2 <= 16)  tl <- uniroot(rhow,lower=-tu2,  upper=-tu2+1.5,const=tu2)$root  
 if (tu2 > 16)             tl <- -tu2
 if (tu2>1   & tu2 <= 50)  tu <- uniroot(rhow,lower=logtu2,upper=tu2,     const=tu2)$root
 if (tu2 > 50)             tu <- uniroot(rhow,lower=logtu2,upper=upp,     const=tu2)$root
 tu <- max(tu,cu); tl <- Izero(tu)}
if (option=="fixed"){tu <- cu; tl <- cl}
tlow <- tl; if (tlow < -16) tlow <- -16
Beta <- integrate(tChiww,tlow,tu,args=list(tl=tl,tu=tu))$value
beta <- Beta/(pweibull(exp(tu),shape=1)-pweibull(exp(tl),shape=1))
list(cl=cl,cu=cu,tl=tl,tu=tu,alpha=alpha,beta=beta)}

