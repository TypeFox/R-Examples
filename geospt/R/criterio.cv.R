assign("criterio.cv",
   function(m.cv){
   m.cv<-as.data.frame(m.cv)
   criterio.mx <- data.frame(matrix(0, ncol=9))
   colnames(criterio.mx) <- c("MPE", "ASEPE", "RMSPE", "MSPE", "RMSSPE", "MAPPE", "CCPE", "R2", "pseudoR2")
   resid.mean <- m.cv[,3]- mean(m.cv[,3])
   criterio.mx[,1] <- mean(m.cv[,4])
   criterio.mx[,2] <- mean((m.cv[,2])^0.5)
   criterio.mx[,3] <- sqrt(sum((m.cv[,4])^2)/nrow(m.cv))
   criterio.mx[,4] <- mean(m.cv[,5])
   criterio.mx[,5] <- sqrt(sum((m.cv[,5])^2)/nrow(m.cv))
   criterio.mx[,6] <- mean(abs(m.cv[,4]/m.cv[,3]))
   criterio.mx[,7] <- cor(m.cv[,1],m.cv[,3])
   criterio.mx[,8] <-  1- sum((m.cv[,4])^2)/sum(resid.mean^2)
   criterio.mx[,9] <- cor(m.cv[,1],m.cv[,3])^2
   criterio.mx
   })
   