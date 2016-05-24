assign("idwST.tcv",
       function (formula, data, n.neigh, C, factor.p, progress=TRUE)
       {
        s0 = cbind(coordinates(data),data["t"]@data)
        z = extractFormula(formula, data, newdata = data)$z
        pred <- as.numeric(NA,length= length(z))
        if(progress)
        pb <- txtProgressBar(min = 0, max = length(z), char = "=", style = 3)
        for (i in 1:length(z)) {
        pred[i] <- idwST(formula, data[-i, ], newdata = data[i, ], n.neigh, C, factor.p, progress=F)[,4]
        if(progress)
        setTxtProgressBar(pb, i)
        }
        if(progress)
        close(pb)
         idw.pred <- data.frame(pred,NA,z,NA,NA,NA,s0)
         colnames(idw.pred) <- c("var1.pred","var1.var","observed","residual","zscore","fold","x","y","t")
         idw.pred[,6] <- 1:length(z)
         idw.pred[,3]<- z
         idw.pred[,4]<- idw.pred[,3]-idw.pred[,1]
         idw.pred
       }
)
