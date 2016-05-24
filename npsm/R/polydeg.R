polydeg <- function(y,x,P,alpha=.05){

     xbar <- mean(x)
     xc <- x - xbar
     xmat <- cbind(xc)
     n <- length(y)

     if(P != 1){
        for(j in 2:P){
            xmat <- cbind(xmat,xc^j)
        }
     }

     ic <- 0
     p <- P
     coll <- rep(0,3)

     while(ic == 0){
        xr <- xmat[,1:(p-1)]
        fitf <- rfit(y~xmat)
        fitr <- rfit(y~xr)
        if(p != 1){
             tst<-drop.test(fitf,fitr)
        } else {
             tst<-drop.test(fitf)
        }
        pv <- tst$p.value
        ftst <- tst$F
        coll <- rbind(coll,c(p,ftst,pv))
        
        if(pv <= alpha){
              ic <- 1
              deg <- p
        }
        if(p == 1){
              ic <- 1
              deg <- 0
        } else {
              xmat <- xr
              p <- p - 1
        }
    }
    coll <- coll [-1,]
    colnames(coll) <- c("Deg","Robust F","p-value")
    list(coll=coll,deg,fitf=fitf)
}
              

     
         
