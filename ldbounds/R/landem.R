"landem" <-
function(t,t2,side,iuse,asf,alpha,phi,ztrun){
   h <- 0.05
   zninf <- -8
   tol <- 0.0000001
   stdv <- sqrt(t2-c(0,t2[-length(t2)])) # These are subroutine "sd"
   sdproc <- sqrt(t2)                    # These are subroutine "sd"
   alph <- alphas(iuse,asf,alpha,phi,side,t)
   za <- zb <- ya <- yb <- nints <- rep(0,length(t))
   pd <- alph$pd
   pe <- alph$pe
   if (pd[1]==0){
      zb[1] <- -zninf
      if (zb[1] > ztrun){
         zb[1] <- ztrun
         pd[1] <- side*(1-pnorm(zb[1]))
         pe[1] <- pd[1]
         if (length(t) > 1) pd[2] <- pe[2]-pe[1]
       }
      yb[1] <- zb[1]*stdv[1]
    }
   else if (pd[1] < 1){
      zb[1] <- qnorm(1-pd[1]/side)
      if (zb[1] > ztrun){
         zb[1] <- ztrun
         pd[1] <- side*(1-pnorm(zb[1]))
         pe[1] <- pd[1]
         if (length(t) > 1) pd[2] <- pe[2]-pe[1]
       }
      yb[1] <- zb[1]*stdv[1]
    }
   if (side==1){
      za[1] <- zninf
      ya[1] <- za[1]*stdv[1]
    }
   else if (side != 1){
      za[1] <- -zb[1]
      ya[1] <- -yb[1]
    }
   nints[1] <- ceiling((yb[1]-ya[1])/(h*stdv[1]))
   if (length(t) >= 2){
      grid <- seq(ya[1],yb[1],length=nints[1]+1) # These are "first"
      last <- dnorm(grid,mean=0,sd=stdv[1])      # These are "first"
      for (i in 2:length(t)){
         if ({pd[i] < 0}|{pd[i] > 1}){
            warning("Possible error in spending function.  May be due to truncation.")
            pd[i] <- min(1,pd[i])
            pd[i] <- max(0,pd[i])
          }
         if (pd[i] < tol){
            zb[i] <- -zninf
            if (zb[i] > ztrun){
               zb[i] <- ztrun
               pd[i] <- side*qp(zb[i]*sdproc[i],last,nints[i-1],ya[i-1],yb[i-1],stdv[i])
               pe[i] <- pd[i]+pe[i-1]
               if (i < length(t)) pd[i+1] <- pe[i+1]-pe[i]
             }
            yb[i] <- zb[i]*sdproc[i]
          }
         else if (pd[i]==1) zb[i] <- yb[i] <- 0
         else if ({pd[i] >= tol}&{pd[i] < 1}){
            yb[i] <- bsearch(last,nints,i,pd[i]/side,stdv[i],ya,yb)
            zb[i] <- yb[i]/sdproc[i]
            if (zb[i] > ztrun){
               zb[i] <- ztrun
               pd[i] <- side*qp(zb[i]*sdproc[i],last,nints[i-1],ya[i-1],yb[i-1],stdv[i])
               pe[i] <- pd[i]+pe[i-1]
               if (i < length(t)){
                  pd[i+1] <- pe[i+1]-pe[i]
                }
             }
            yb[i] <- zb[i]*sdproc[i]
          }
         if (side==1){
            ya[i] <- zninf*sdproc[i]
            za[i] <- zninf
          }
         else if (side==2){
            ya[i] <- -yb[i]
            za[i] <- -zb[i]
          }
         nints[i] <- ceiling((yb[i]-ya[i])/(h*stdv[i]))
         if (i < length(t)){
            hlast <- (yb[i-1]-ya[i-1])/nints[i-1]                 # These are "other"
            x <- seq(ya[i],yb[i],length=nints[i]+1)               # These are "other"
            last <- fcab(last,nints[i-1],ya[i-1],hlast,x,stdv[i]) # These are "other"
          }
       }
    }
   ans <- list(lower.bounds=za,upper.bounds=zb,exit.pr=pe,diff.pr=pd,spend=alph$spend)
   return(ans)
 }

