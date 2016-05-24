errorObsBin <- function(sigma, STAR, parallel=FALSE, corr=c(0,0,0,0,0,0,0)) {

  STARp <- NULL
  if(!is.matrix(STAR)) 
    locSTAR <- matrix(STAR, nrow=1)
  else
    locSTAR <- STAR

  num <- nrow(locSTAR)
  if(parallel) {
      lf <- parallel::mclapply
    } else {
      lf <- lapply
    }

  corrTeff <- corr[1]
  corrFeH <- corr[3]
  corrM <- corr[6]
  corrR <- corr[7]

  STARp <- do.call(rbind, lf(1:num, function(i) {
    if(length(sigma) == 7) {
      mysigma <- c(sigma, sigma)
    } else {
      mysigma <- sigma
    }
    mysigma[4] <- mysigma[4]*locSTAR[i,4];
    mysigma[5] <- mysigma[5]*locSTAR[i,5];
    mysigma[6] <- mysigma[6]*locSTAR[i,6];
    mysigma[7] <- mysigma[7]*locSTAR[i,7];
    mysigma[11] <- mysigma[11]*locSTAR[i,13];
    mysigma[12] <- mysigma[12]*locSTAR[i,14];
    mysigma[13] <- mysigma[13]*locSTAR[i,15];
    mysigma[14] <- mysigma[14]*locSTAR[i,16];
    
    COV <- diag(mysigma^2)
    COV[is.na(COV)] <- 0
    
    # corr Teff
    COV[1,8] <- COV[8,1] <- corrTeff * mysigma[1]*mysigma[8];
    # corr [Fe/H]
    COV[3,10] <- COV[10,3] <- corrFeH * mysigma[3]*mysigma[10];
    # corr M
    COV[6,13] <- COV[13,6] <- corrM * mysigma[6]*mysigma[13];
    # corr R
    COV[7,14] <- COV[14,7] <- corrR * mysigma[7]*mysigma[14];
    
    tmpSTARp <- mvrnorm(1, locSTAR[i,c(1:7, 10:16)], COV);
    return(c(tmpSTARp, locSTAR[i,c(6:9, 15:18, 3, 12)])) } ))

  colnames(STARp)[c(15:24)] <- c("origM1", "origR1", "logage1", "pcAge1", "origM2", "origR2", "logage2", "pcAge2", "origFeH1", "origFeH2")
  STARp <- STARp[,c(1:7, 17:18, 8:14, 21:22, 15:16, 19:20, 23:24)]
  return(STARp)
}



estimateBin <- function(data, STAR, sigma, thr, sel, parallel=FALSE) {

  # currently unused
  pw <- 4
  res <- 0

  # max age difference (in Gyr)
  tsp <- 0.01
  
  if(!is.matrix(data))
    stop("the recovery grid must be of \"matrix\" class")
  if(dim(data)[2] != 9)
    stop("uncorret number of columns in the data grid")

  if(!is.matrix(STAR)) 
    locSTAR <- matrix(STAR, nrow=1)
  else
    locSTAR <- STAR

  if(length(sigma) == 7) {
    mysigma <- c(sigma, sigma)
  } else {
    mysigma <- sigma
  }

  dim <- nrow(locSTAR)
  if(parallel) {
    lf <- parallel::mclapply
  } else {
    lf <- lapply
  }
  res <- lf(1:dim, function(i) {
    val <- .Call("lik4bin", data, STAR[i,], mysigma, thr, sel, pw, res, tsp)
    if(!is.null(val)) {
      val <- c(i, val, c(STAR[i,c(19:20,8:9,3,21:22,17:18,12,23:24)]));
      return(val)}})
  res <- do.call(rbind, res)
  if(is.null(res)) return(res)
  res <- as.data.frame(res, row.names=1:dim)
  colnames(res) <- c("id", "M1", "R1", "estlogage1", "age1", "lik1", "num1", "M2", "R2", "estlogage2", "age2", "lik2", "num2", "M1b", "R1b", "estlogage1b", "age1b", "M2b", "R2b", "estlogage2b", "age2b", "numbin", "r", "trueM1", "trueR1", "logAge1", "pcAge1", "FeH1", "trueM2", "trueR2", "logAge2", "pcAge2", "FeH2", "trueFeH1", "trueFeH2")

  res <- res[,c(1,2,3,5,8,9,11,14,15,17,18,19,21,23)]
  
  estage <- 0.5*(res$age1b + res$age2b)
  estage12 <- 0.5*(res$age1 + res$age2)
  res$ageBin <- estage
  res$ageMean <- estage12


  return(res)
}


block <- function(grid) {
  nr <- nrow(grid)
  nm <- length(unique(grid[,6]))
  np <- 110 # track points
  nb <- nr/(nm*np)
  block <- gl(nb, np*nm)
  return(block)
}

sampleBinStar <- function(n, grid, block, restrict=TRUE, parallel=FALSE) {

  tspread <- 1e-7
  
  nrow <- nrow(grid)
  nm <- length(unique(grid[,6]))
  np <- 110  # track points
  nb <- nrow/(nm*np)
  dimblock <- np*nm

  if(restrict) {
    res <- (1e-9*10^grid[,8]) < 14
  } else {
    res <- rep(TRUE, nrow)
  }
  
  sel <- sample((1:nrow)[res], n, replace=TRUE)
  star1 <- grid[sel,]
  wblock <- as.numeric(block[sel])

  time <- 10^grid[,8]
  
  if(parallel) {
    lf <- parallel::mclapply
  } else {
    lf <- lapply
  }

  # secondary stars at same initial [Fe/H]
  # max age difference  1/tspread yr
  star2 <- do.call(rbind, lf(1:n, function(i) {
    wb <- wblock[i]
    selgrid <- ((wb-1)*np*nm+1):(wb*np*nm)
    usegrid <- grid[selgrid,]
    timegrid <- time[selgrid]
    seltime <- which(abs(tspread*(10^star1[i,8]- timegrid)) < 1)
    if(length(seltime) > 1) {
      sels <- sample(seltime, 1)
    } else {
      sels <- seltime
    }
    return(usegrid[sels,])
  }))

  
  # luminosity computation
  coeff <- 3.3338e7
  L1 <- star1[,7]^2*star1[,1]^4/coeff^2
  L2 <- star2[,7]^2*star2[,1]^4/coeff^2

  # output data frame
  res <- cbind(star1, star2)

  # primary stars have higher luminosity...
  inv <- which(L2 > L1)
  res[inv,] <- res[inv, c(10:18,1:9)]

  return(res)
}
