"glan" <-
function(t,za,zb,drft){
   h <- 0.05
   stdv <- sqrt(t-c(0,t[-length(t)])) # These are subroutine "sd"
   sdproc <- sqrt(t)                  # These are subroutine "sd"
   yb <- zb*sdproc-drft*t
   ya <- za*sdproc-drft*t
   nints <- ceiling((yb-ya)/(h*stdv))
   qneg1 <- pnorm(za[1],mean=drft*t[1]/stdv[1])
   qpos1 <- 1-pnorm(zb[1],mean=drft*t[1]/stdv[1])
   cp <- matrix(0,length(t),2)
   cp[1,] <- c(qpos1,qneg1)
   if (length(t) >= 2){
      grid <- seq(ya[1],yb[1],length=nints[1]+1) # These are "first"
      last <- dnorm(grid,mean=0,sd=stdv[1])      # These are "first"
      for (i in 2:length(t)){
         cpr <- cprob(last,nints,ya,yb,i,stdv[i])
         cp[i,] <- c(cpr[[1]],cpr[[2]])
         if (i < length(t)){
            hlast <- (yb[i-1]-ya[i-1])/nints[i-1]                 # These are "other"
            x <- seq(ya[i],yb[i],length=nints[i]+1)               # These are "other"
            last <- fcab(last,nints[i-1],ya[i-1],hlast,x,stdv[i]) # These are "other"
          }
       }
    }
   pr <- sum(cp)
   ans <- list(pr=pr,qpos=cp[,1],qneg=cp[,2])
   return(ans)
 }

