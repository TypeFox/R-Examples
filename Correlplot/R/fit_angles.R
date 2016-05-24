fit_angles <-
function(R,ifun="cos",ntrials=10,verbose=FALSE) {
n   <- nrow(R)
aux <- list(R=R,n=n,FUN=ifun)
vmin <- NULL
vtheta <- NULL
for(i in 1:ntrials) {
  x0  <- runif(n-1,min=0.2,max=6) # initial angles
  init <- obj_fun(x0,aux)
  optout <- nlminb(start=x0,objective=obj_fun,scale=1,aux=aux)
  minimum <- optout$objective
  vmin <- c(vmin,minimum)
  theta <- c(0,optout$par)
  vtheta <- rbind(vtheta,theta)
  if(verbose) cat(i,init,minimum,"\n")
}
ind <- which.min(vmin)
theta <- vtheta[ind,]
minim <- vmin[ind]
if(verbose) cat("mininum of objective function found:",minim,"\n")
return(theta)
}

