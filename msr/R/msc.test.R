
#sample from a square crystal in d dimensions with slope gamma and gaussina noise sigma
linear.gaussianNd.sample.full <- function(gamma = 1, sigma = 0.25, k=100, dim = 2){
  d <- matrix(nrow = k, ncol = dim)
  for(i in 1:dim){
    d[,i] <- runif(k)/sqrt(dim)
  }
  fx <- -sqrt(rowSums(d^2))*gamma + rnorm(n=k, sd = sigma)

  cbind(d, fx)
}



#compute p value of introducing crystal i2 in ms2 in crystal i1 in ms1
#ms1 and ms2 should be morse samle complicies on the same data
msc.p <- function(ms1, i1, per, knn=ms1$knn, nruns = 10000){
  l1 <- lm( ms1$y[ms1$crystals==i1] ~ ms1$x[ms1$crystals==i1, ])

  sigma <- sqrt(var(l1$residuals))
  imax <- which.max(l1$fitted)
  imin <- which.min(l1$fitted)
  d <- ms1$x[imax, ] - ms1$x[imin, ]
  d <- sqrt(sum(d^2))
  gamma <- (max(l1$fitted) - min(l1$fitted)) * d
  
  p <- 0
  k <- length(l1$residuals)
  dim <- ncol(ms1$x)
  for( i in 1:nruns){
    d <- linear.gaussianNd.sample.full(gamma = gamma, sigma=sigma, k=k, dim=dim)
    ms <- msc.nn(y = d[, dim+1], x = d[, 1:dim], knn=knn, pLevel = 0)
    if(length(ms$persistence )> 1){
      if( ms$persistence[length(ms$persistence)-1] > per){
        p <- p + 1
      }
    }
  }
  p <- p/nruns
  p 
}



msc.test <- function (ms, nruns = 1000) 
{
   if(ms$nLevels == 1){
     stop()
   }

   pl = length(ms$persistence)
   pvals <- c() 
   for(i in 1:(ms$nLevels-1)){
     ms1 <- ms$mscl[[i]]
     ms2 <- ms$mscl[[i+1]]
     np <- length(ms1$crystalsSize)
     pLevel <- ms$persistence[pl-i]
     pv <- matrix(ncol=2, nrow=np)
     for(pID in 1:np){
       pv[pID, 1] <- pID
       s <- ms1$crystals == pID
       u <- unique(ms2$crystals[s])
       if(length(u) > 1){
         val <- msc.p(ms1, pID, pLevel, knn=ms$knn, nruns = nruns)
         pv[pID,2] <- val
       }
       else{
        pv[pID,2] <- NA
       }
     } 
     pvals[[i]] <- pv
   }
   pvals
}

