assign("idwST.cv",
       function (formula, data, n.neigh, C, factor.p, progress=TRUE)
       {
         z = extractFormula(formula, data, newdata = data)$z
         idw.pred <- as.data.frame(matrix(NA, nrow = length(z), ncol = 5))
         colnames(idw.pred) <- c("x", "y", "t", "var1.pred", "var1.var")
         if(progress)
         pb <- txtProgressBar(min = 0, max = length(z), char = "=", style = 3)
         for (i in 1:length(z)) {
           idw.pred[i, 4] <- idwST(formula, data[-i, ], newdata = data[i, ], n.neigh, C, factor.p, progress=F)[, 4]
           if(progress)
           setTxtProgressBar(pb, i)
         }
         if(progress)
         close(pb)
         RMSPE <- sqrt(sum((idw.pred$var1.pred - z)^2)/length(z))
         RMSPE
       }
)