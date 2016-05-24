####################################
#### cut MCMC chains
####################################
cutChain <- function(res, chain.ind=2:4, burnin, thinning){
  ind <- seq(burnin+1, dim(res[[1]])[2], by = thinning)
  if( length(setdiff(chain.ind, 1:length(res)))>0 ) stop("chain.ind is out of range!")
  lapply(res[chain.ind], function(x, ind){ 
    if(is.matrix(x)){ xx <- x[, ind]
                    } else{xx <- x[ind]}
    xx
    }, ind)
}
####################################
#### mix multiple MCMC chains
####################################
mixChain <- function(res.m.prl){
res.mix <- vector("list", length = length(res.m.prl[[1]]))
names(res.mix) <- names(res.m.prl[[1]])
for(i in 1:length(res.m.prl)){
    for(j in 1:length(res.mix)){
      if( is.matrix(res.m.prl[[i]][[j]]) ){
        res.mix[[j]] <- cbind( res.mix[[j]], res.m.prl[[i]][[j]] )
        } else{
            res.mix[[j]] <- c( res.mix[[j]], res.m.prl[[i]][[j]] )
          }
      }
  
  }
res.mix
}
####################################
#### Plot autocorrelations
####################################
plotACF <- function(S.mcmc, lags = NULL, ...){
for(i in 1:nrow(S.mcmc)){
  temp <- acf(S.mcmc[i,], plot=FALSE, lag.max = lags)
  if(i==1){
    plot(as.vector(temp$lag), as.vector(temp$acf), 
        type="l", xlab="Lag", ylab="Autocorrelation", ...)
  } else lines(temp$lag, temp$acf)
}
}
####################################
#### Find mode 
####################################
findMode <- function(x, ...){
  tmp <- density(x, ...)
  ind <- which(tmp$y==max(tmp$y))
  tmp$x[ind[1]]
}
####################################
#### END
####################################