describe.r <-
function(x, na.rm=TRUE, tr=.2, type=3) {
  if(sum(sum(x > 1, na.rm=T) + sum(x < -1, na.rm=T)) > 0) {stop("All values in 'x' must be correlation coefficients between -1.00 and 1.00.")}
  if(sum(sum(x==1, na.rm=T) + sum(x==-1, na.rm=T) > 0)) {warning("'x' contains values of 1.00 or -1.00 resulting in an infinite value.")}
  get.valid <- function(x) {
    sum(!is.na(x))
  }
  get.miss <- function(x) {
    sum(is.na(x))
  }
  z <- fisherz(x)
  if(is.null(dim(x)[2])) {
    vars <- 1
    valid <- get.valid(z)
    miss <- get.miss(z)
    stats <- matrix(NA, nrow=1, ncol=11)
    stats[1,1] <- mean(z, na.rm=na.rm)
    stats[1,2] <- sd(z, na.rm=na.rm)
    stats[1,3] <- median(z, na.rm=na.rm)
    stats[1,4] <- mean(z, na.rm=na.rm, trim=tr)
    stats[1,5] <- mad(z, na.rm=na.rm)
    stats[1,6] <- min(z, na.rm=na.rm)
    stats[1,7] <- max(z, na.rm=na.rm)
    stats[1,8] <- stats[1,7] - stats[1,6]
    stats[1,9] <- skew(z, na.rm=na.rm)
    stats[1,10] <- kurtosi(z, na.rm=na.rm, type=type)
    stats[1,11] <- stats[1,2] / sqrt(valid)
  }
  else {
    vars <- 1:ncol(z)
    valid <- apply(z, 2, get.valid)
    miss <- apply(z, 2, get.miss)
    stats <- matrix(NA, nrow=ncol(z), ncol=11)
    stats[,1] <- apply(z, 2, mean, na.rm=na.rm)
    stats[,2] <- apply(z, 2, sd, na.rm=na.rm)
    stats[,3] <- apply(z, 2, median, na.rm=na.rm)
    stats[,4] <- apply(z, 2, mean, na.rm=na.rm, trim=tr)
    stats[,5] <- apply(z, 2, mad, na.rm=na.rm)
    stats[,6] <- apply(z, 2, min, na.rm=na.rm)
    stats[,7] <- apply(z, 2, max, na.rm=na.rm)
    stats[,8] <- stats[,7] - stats[,6]
    stats[,9] <- apply(z, 2, skew, na.rm=na.rm)
    stats[,10] <- apply(z, 2, kurtosi, na.rm=na.rm, type=type)
    stats[,11] <- stats[,2] / sqrt(valid)
  }
  out <- data.frame(vars, valid, miss, fisherz2r(stats))
  colnames(out) <- c("var", "n", "miss", "mean", "sd", "median", "trimmed", "mad", "min", "max", "range", "skew", "kurtosis", "se")
  return(out)
}
