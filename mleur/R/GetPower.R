GetPower <-
function(phi, n, NSIM=1000, tests=c("DF", "MLEp", "MLEn", "MCT"),
            noiseDist=c("normal", "t", "stable", "GARCH11"),
            df=5, ALPHA=1.5, BETA=0, alpha=0.2, beta=0.7) 
{
#compute power by simulation
hitsMCT <- hitsMLEn <- hitsMLEp <- hitsDF <- 0
noiseD <- noiseDist[match(noiseDist, noiseDist)[1]]
for (isim in 1:NSIM){
  y <- simar1(phi=phi, n=n, noiseDist=noiseD, df=df, 
        ALPHA=ALPHA, BETA=BETA, alpha=alpha, beta=beta)
  #DF test
  if (!is.na(match("DF", tests))) {
    out <- dftest(y)
    hitsDF <- hitsDF + as.numeric(out$tau<out$criticalValues[2])
  }
  #using rsr methods
  if (!is.na(match("MLEp", tests))) {
    out <- as.vector(mleur(y))
    hitsMLEp <- hitsMLEp + as.numeric(out[1]<out[3])
    }
  if (!is.na(match("MLEn", tests))) {
    out <- as.vector(mleur(y, type="n"))
    hitsMLEn <- hitsMLEn + as.numeric(out[1]<out[3])
    }
  #MCT
  if (!is.na(match("MCT", tests))) {
    out <- mctest(y)
    hitsMCT <- hitsMCT+as.numeric(out<=0.05)
    }
  }
ans <- c(hitsDF, hitsMLEp, hitsMLEn, hitsMCT)/NSIM
ind <- c( !is.na(match("DF", tests)),   !is.na(match("MLEp", tests)),
          !is.na(match("MLEn", tests)), !is.na(match("MCT", tests))
        )
ans <- ans[ind]
names(ans) <- (c("DF", "MLEp", "MLEn", "MCT"))[ind]
list(power=ans, phi=phi, n=n, NSIM=NSIM, noiseDist=noiseD, MOE=1.96*0.5/sqrt(NSIM))
}

