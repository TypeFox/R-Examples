### function of GCV for doubly penalized BJ method
     gcv.enet <- function(x, y, cens){
     tmp <- expand.grid(y=seq(0.01, 1, by=0.05), x=c(0, 0.001, 0.01, 0.1, 1, 10))
     tmp <- cbind(tmp[,2], tmp[,1])
     nn <- dim(tmp)[1]
     mse.tr <- rep(NA, nn)
     for(i in 1:nn){
     cat("ENET GCV i=",i,"of", nn,"lambda=",tmp[i,1],"s=",tmp[i,2],"\n")
     mse.tr[i] <- bujar(x=x, y=y, cens=cens, learner="enet", lamb=tmp[i,1], s= tmp[i,2], iter.bj=20, max.cycle=1, tol=1e-8, mimpu=NULL, tuning=TRUE, trace=TRUE)$mse.tr
     }
     tmp1 <- which.min(mse.tr)
     lamb <- tmp[tmp1, 1]
     s <- tmp[tmp1, 2]
     cat("\nBJ with enet GCV results: lambda=",lamb,"s=",s,"\n")
     return(list(lambda=lamb, s=s, mode="fraction"))
    } 
