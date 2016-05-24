`Nelson.Siegel` <-
function( rate, maturity )
  {
    rate <- try.xts(rate,error=as.matrix)
    if(ncol(rate)==1) rate<-matrix(as.vector(rate),1,nrow(rate))
    pillars.number <- length(maturity)
    lambdaValues <- seq(maturity[1], maturity[ pillars.number ], by=0.5)

    FinalResults <- matrix(0, nrow(rate), 4)
    colnames( FinalResults ) <- c("beta_0","beta_1","beta_2","lambda")
    j <- 1
    while(j <= nrow(rate) )
      {
        InterResults <- matrix(0, length(lambdaValues), 5)
        colnames( InterResults ) <- c("beta0","beta1","beta2","lambda","SSR")
        for( i in 1:length(lambdaValues))
          {
            lambdaTemp <- optimize(.factorBeta2,interval=c(0.001,1),
              maturity=lambdaValues[i],maximum=TRUE)$maximum
            InterEstimation <- .NS.estimator(as.numeric(rate[j,]), maturity, lambdaTemp)
            BetaCoef <- InterEstimation$Par
	    if( BetaCoef[1]>0 & BetaCoef[1]<20)
              {
                SSR <- sum(InterEstimation$Res^2)
                InterResults[i,] <- c(BetaCoef, lambdaTemp, SSR)
              } else
            {
              InterResults[i,] <- c(BetaCoef,lambdaValues[i],1e+5)
            } 
          }
        BestRow <- which.min(InterResults[,5])
        FinalResults[j,] <- InterResults[BestRow,1:4]
        j <- j+1
      }
    reclass( FinalResults, rate )
  }
