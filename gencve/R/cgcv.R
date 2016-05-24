cgcv <- function(X, y, yh = yh_NN, MaxIter = 1000, d = ceiling(length(y)/10),
                 NCores = 1, libs = character(0), seed = "default", ...) {
  stopifnot(is.character(libs))
  n <- length(y)
  stopifnot(is.factor(y) && n>=10)
  stopifnot(is.numeric(d) && d>=1 && n>d)
  d <- ceiling(d)
#KIter
  KIter <- function(K, ...) {
    avMrate <- 0
    for (i in 1:K) {
      iTe <- sample(n, size=d)
      iTr <- (1:n)[!(1:n %in% iTe)]
      trdf <- data.frame(X[iTr,], y=y[iTr])
      tedf <- data.frame(X[iTe,], y=y[iTe])
      mrate <- yh(trdf, tedf, ...)
      avMrate <-  avMrate + (mrate-avMrate)/i
    }
    avMrate
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
    epe <- OUT
  }  else {#NCores > 1 ...
    Ks <- rep(ceiling(MaxIter/NCores), NCores) #load balancing
    cl <- makeCluster(spec=NCores, type="PSOCK")
    if (!identical(seed, "default")) {
               clusterSetRNGStream(cl, iseed=seed)
	       }
    #Export variables
    clusterExport(cl, list("n", "d", "X", "y"), envir=environment())
    clusterExport(cl, "yh", envir=environment())
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
    epe <- apply(OUT, 1, mean)
  }#end NCores>1
  #
 epe
}
