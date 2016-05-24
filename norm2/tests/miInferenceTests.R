library(norm2)

# generate M=25 imputations for cholesterol data
data(cholesterol)
emResult <- emNorm(cholesterol)
set.seed(999)
mcmcResult <- mcmcNorm(emResult, iter=2500, impute.every=100)

M <- 25
est.list <- std.err.list <- as.list( NULL )
est <- std.err <- df.complete <- numeric(3)
names(est) <- names(std.err) <- names(df.complete) <-
   c("Y2 versus Y1", "Y3 versus Y1", "Y3 versus Y2")
for( m in 1:M ){
   # extract the m-th imputed dataset
   y.imp <- data.frame( mcmcResult$imp.list[[m]] )
   # contrast: Y2 versus Y1
   diff <- y.imp$Y2 - y.imp$Y1
   est[1] <- mean( diff )
   std.err[1] <- sqrt( var(diff) / length(diff) )
   df.complete[1] <- length(diff) - 1
   # contrast: Y3 versus Y1
   diff <- y.imp$Y3 - y.imp$Y1
   est[2] <- mean( diff )
   std.err[2] <- sqrt( var(diff) / length(diff) )
   df.complete[2] <- length(diff) - 1
   # contrast: Y3 versus Y2
   diff <- y.imp$Y3 - y.imp$Y2
   est[3] <- mean( diff )
   std.err[3] <- sqrt( var(diff) / length(diff) )
   df.complete[3] <- length(diff) - 1
   # append lists
   est.list[[m]] <- est
   std.err.list[[m]] <- std.err
   }

print( miInference( est.list, std.err.list, df.complete=df.complete ) )
