"adaptn" <-
function(rso,cl,cu,option){
# adaptive cut-off
alpha  <- NA; n <- length(rso)
if (option=="adaptive"){
 rhos  <- sort(rso^2/2)
 Fo    <- pchisq(2*rhos,df=1) 
 dsr   <- Discr(rhos,Fo,cu^2/2); if (is.na(dsr$j)) return(list(tu=NA))
 alpha <- dsr$r.min
 tu2   <- as.numeric(quantile(rhos,probs=alpha)); tu <- sqrt(2*tu2)
 tu    <- max(tu,cu); tl    <- -tu}
if (option=="fixed"){tu <- cu; tl <- cl}
Beta  <- -tu*dnorm(tu)+tl*dnorm(tl)+pnorm(tu)-pnorm(tl)
beta  <-  Beta/(pnorm(tu)-pnorm(tl))
list(cl=cl,cu=cu,tl=tl,tu=tu,alpha=alpha,beta=beta)}

