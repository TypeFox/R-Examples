assign("idwST",
       function(formula, data, newdata, n.neigh, C, factor.p, progress=TRUE){
         s = cbind(coordinates(data),data["t"]@data)
         s0 = cbind(coordinates(newdata),newdata["t"]@data)
         z = extractFormula(formula, data, newdata)$z
         S <- scale(s)
         dist.newdata <- rdist(S, cbind(standardize(s0[,1],mean(s[,1]),sd(s[,1])),standardize(s0[,2],mean(s[,2]),sd(s[,2])),
                                        C*standardize(s0[,3],mean(s[,3]),sd(s[,3]))))
         Pred <- as.numeric(NA,length= nrow(coordinates(newdata)))
         idw0 <- function(z, dist.newdata, n.neigh, factor.p){
           vec.orden <- order(dist.newdata)
           vc <- vec.orden[1:n.neigh]
           dist.vec.cerca <- dist.newdata[vc]
           Lambda <- dist.vec.cerca^(-factor.p)/sum(dist.vec.cerca^(-factor.p))
           pred <- t(Lambda)%*%z[vc]
           pred
         }
         if(progress)
         pb <- txtProgressBar(min = 0, max = nrow(coordinates(newdata)), char = "=", style = 3)
         for (i in 1:nrow(coordinates(newdata))){
         Pred[i] <- idw0(z, dist.newdata=dist.newdata[,i], n.neigh, factor.p=factor.p)
         if(progress)
         setTxtProgressBar(pb, i)
         }
         if(progress)
         close(pb)
         idw.pred <- data.frame(s0,Pred,NA)
         names(idw.pred) <- c("x","y","t","var1.pred","var1.var")
         idw.pred
       }
)  
