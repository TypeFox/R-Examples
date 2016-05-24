`Svensson` <-
function( rate, maturity )
  {

    rate <- try.xts(rate, error=as.matrix)
    if(ncol(rate)==1) rate<-matrix(as.vector(rate),1,nrow(rate))
    pillars.number <- length(maturity)
    Tau1Values <- seq(maturity[1], median(maturity), by=1)
    Tau2Values <- seq(median(maturity), maturity[pillars.number], by=1.5)
    
    FinalResults <- matrix(0, nrow(rate), 6)
    FinalResultsTau2 <- matrix(0, length(Tau1Values), 7)   
    colnames( FinalResults ) <- c("beta_0","beta_1","beta_2","beta_3","tau1","tau2" )
    j <- 1
    while(j <= nrow(rate) )
      {
        InterResultsTau1 <- matrix(0,length(Tau1Values), 7)
        InterResultsTau2 <- matrix(0,length(Tau2Values), 7)
        # colnames( InterResults ) <- c("beta0","beta1","beta2","beta_3","Tau1","Tau2","SSR")
        for( i in 1:length(Tau1Values))
          {
            Tau1Temp <- optimize(.beta2Spot,interval=c(0.001,max(Tau1Values)),maturity=Tau1Values[i],maximum=TRUE)$maximum
            for( a in 1:length(Tau2Values))
              {
                Tau2Temp <- optimize(.beta2Spot,interval=c(0.001,maturity[pillars.number]),maturity=Tau2Values[a],maximum=TRUE)$maximum
                InterEstimation <- .NSS.estimator(as.numeric(rate[j,]), maturity, Tau1Temp, Tau2Temp)
                BetaCoef <- InterEstimation$Par
                SSR <- sum(InterEstimation$Res^2)
                InterResultsTau2[a,] <- c(BetaCoef, Tau1Temp, Tau2Temp, SSR)
              }
            BestRowTau2 <- which.min(InterResultsTau2[,7])
            FinalResultsTau2[i,] <- InterResultsTau2[BestRowTau2,]
          }
        BestRow <- which.min(FinalResultsTau2[,7])
        FinalResults[j,] <- FinalResultsTau2[BestRow,1:6]
        j <- j+1
      }
    reclass( FinalResults, rate )
  }

