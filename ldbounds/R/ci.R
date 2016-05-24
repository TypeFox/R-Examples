"ci" <-
function(conf,value,t,za,zb){
   zb[length(t)] <- value
   zcrit <- qnorm(1-(1-conf)/2)
   limit <- (value+c(-1,1)*zcrit)/sqrt(t[length(t)])
   target <- c(0,1)*conf+(1-conf)/2
   lim1 <- bisect(t,za,zb,target[1],limit[1],upper=TRUE)
   lim2 <- bisect(t,za,zb,target[2],limit[2],upper=TRUE)
   lim <- list(lower.limit=lim1,upper.limit=lim2)
   return(lim)
 }

