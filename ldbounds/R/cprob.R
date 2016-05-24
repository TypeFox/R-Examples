"cprob" <-
function(last,nints,ya,yb,i,stdv){
   hlast <- (yb[i-1]-ya[i-1])/nints[i-1]
   grid <- seq(ya[i-1],yb[i-1],length=nints[i-1]+1)
   pupr <- (1-pnorm(yb[i],mean=grid,sd=stdv))*last
   plow <- pnorm(ya[i],mean=grid,sd=stdv)*last
   tqpos <- 0.5*hlast*(2*sum(pupr)-pupr[1]-pupr[length(pupr)]) # This is "trap"
   tqneg <- 0.5*hlast*(2*sum(plow)-plow[1]-plow[length(plow)]) # This is "trap"
   ans <- list(qpos=tqpos,qneg=tqneg)
   return(ans)
 }

