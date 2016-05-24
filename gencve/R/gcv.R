gcv <- function(X, y,  MaxIter = 1000, d = ceiling(length(y)/10), NCores = 1,
             cost = mse, yhat = yhat_lm, libs = character(0),
             seed = "default", ...) {
  stopifnot(is.character(libs))
  n <- length(y)
  stopifnot(is.numeric(y) && n>=10)
  stopifnot(is.numeric(d) && d>=1 && n>d)
  d <- ceiling(d)
#KIter
  KIter <- function(K, ...) {
    avCOST <- 0
    SSQ <- 0
    coryHatyte <- 0
    for (i in 1:K) {
      iTe <- sample(n, size=d)
      iTr <- (1:n)[!(1:n %in% iTe)]
      trdf <- data.frame(X[iTr,], y=y[iTr])
      tedf <- data.frame(X[iTe,], y=y[iTe])
      yHat <- yhat(trdf, tedf, ...)
      COST <- cost(y[iTe], yHat)
      SSQ <- SSQ + (i-1)*(COST-avCOST)*((COST-avCOST)/i)
      avCOST <-  avCOST + (COST-avCOST)/i
      coryHatyte <- coryHatyte + cor(yHat, tedf$y)
    }
    c(avCOST, SSQ, coryHatyte/K)
  }
#
  if (NCores==1) { #NCORES=1 ...
    if (length(libs) > 0) {
    	library(libs, character.only=TRUE)
	  }
    if (!identical(seed, "default")) {
    		set.seed(seed)
	  } #otherwise use random seed
    OUT <- KIter(MaxIter, ...)
    epe <- OUT[1]
    sd_epe <- sqrt(OUT[2]/(MaxIter^2))
    coryHatyte <- OUT[3]
  }  else {#NCores > 1 ...
    Ks <- rep(ceiling(MaxIter/NCores), NCores) #load balancing
    cl <- makeCluster(spec=NCores, type="PSOCK")
    if (!identical(seed, "default")) {
               clusterSetRNGStream(cl, iseed=seed)
	       }
    #Export variables
    clusterExport(cl, list("n", "d", "X", "y"), envir=environment())
    clusterExport(cl, "yhat", envir=environment())
    clusterExport(cl, "cost", envir=environment())
    #Export library or libraries
    if (length(libs) > 0) {
      for (i in 1:length(libs)) {
        libi <- libs[i]
        clusterExport(cl, "libi", envir=environment())
        clusterEvalQ(cl, library(libi, character.only=TRUE))
      }
    }
    #parallel sapply
    OUT <- parSapply(cl, Ks, KIter, ...)
    stopCluster(cl) #stop cluster
    OUT <- matrix(OUT, byrow=TRUE, ncol=3)
    epe <- mean(OUT[,1])
    sd_epe <- sqrt(sum(OUT[,2])/(sum(Ks)^2))
    coryHatyte <- mean(OUT[,3])
  }#end NCores>1
  #
  SNR <- function(y, EPE) {
    MST <- mean((y-mean(y))^2)
    (MST-EPE)/EPE
  }
  snr <- SNR(y, epe)
  r <- (var(y)-epe)/var(y)
  r <- sign(r)*sqrt(abs(r))
  cbind(epe=epe, sd_epe=sd_epe, pcorr=r)
}
