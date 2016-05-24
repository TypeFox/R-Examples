errorObs <- function(sigma, STAR, parallel=FALSE, corr=0) {

  STARp <- NULL
  if(!is.matrix(STAR)) 
    locSTAR <- matrix(STAR, nrow=1)
  else
    locSTAR <- STAR

  num <- nrow(locSTAR)
  nc <- ncol(locSTAR)
  
  if(nc < 9) {
    s <- matrix(NA, nrow=num, ncol=9)
    s[,1:nc] <- locSTAR
    locSTAR <- s
  }
    
  if(parallel) {
      lf <- parallel::mclapply
    } else {
      lf <- lapply
    }
  STARp <- do.call(rbind, lf(1:num, function(i) {
    mysigma <- sigma;
    # sigma[4,5] are expressed as % ...
    mysigma[4] <- mysigma[4]*locSTAR[i,4];
    mysigma[5] <- mysigma[5]*locSTAR[i,5];
    # seismic parameters correlation 
    COV <- diag(mysigma^2)
    COV[4,5] <- mysigma[4]*mysigma[5]*corr;
    COV[5,4] <- mysigma[4]*mysigma[5]*corr;
    # sample the observation from multivariate normal distribution
    # assuming COV covariance matrix
    tmpSTARp <- MASS::mvrnorm(1, locSTAR[i,1:7], COV);
    return(c(tmpSTARp, locSTAR[i,c(6:9,3)])) } ))

  colnames(STARp)[c(8:9,11:12)] <- c("origM", "origR", "pcAge", "origFeH")
  return(STARp)

}

estimate <- function(data, STAR, sigma, thr, sel, parallel=FALSE) {
# main function: estimate M, R, and age 

  if(!is.matrix(data))
    stop("the recovery grid must be of \"matrix\" class")
  if(dim(data)[2] != 9)
    stop("uncorret number of columns in the data grid")
  
  if(!is.matrix(STAR))
    locSTAR <- matrix(STAR, nrow=1)
  else
    locSTAR <- STAR

  if(is.unsorted(data[,1]))
    stop("the recovery grid must be sorted on effective temperature")
  
  dim <- dim(locSTAR)[1]

  lf <- ifelse(parallel, parallel::mclapply, lapply)
    
  res <- lf(1:dim, function(i) {
    val <- .Call("lik4", data, locSTAR[i,], sigma, thr, sel)

    if(!is.null(val)) {
      val <- c(i, val, locSTAR[i,c(8:12)])
      return(val) }
  })
    
  res <- do.call(rbind, res)
  res <- as.data.frame(res, row.names=1:dim)
    
  if(nrow(res) > 0) {
    colnames(res) <- c("id", "M", "R", "estlogage", "age", "num", "trueM", "trueR", "logAge", "pcAge", "FeH")
    res <- res[,c(1,2,3,5)]
  }
  
#  res$errorM <- (res$M - res$trueM)/res$trueM
#  res$errorR <- (res$R - res$trueR)/res$trueR
 
  return(res)
}



sampleStar <- function(n, grid, restrict=TRUE) {

  if(!is.matrix(grid))
    stop("the recovery grid must be of \"matrix\" class")
  
  nrow <- nrow(grid)
  # select only models with age < 14 Gyr
  if(restrict) {
    res <- (1e-9*10^grid[,8]) < 14
  } else {
    res <- rep(TRUE, nrow)
  }
  sel <- sample((1:nrow)[res], n, replace=TRUE)

  return(grid[sel,])
}
