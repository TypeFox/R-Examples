"jsmurf" <-
function(y, x = 1:length(y), x0 = 2 * x[1] - x[2], q, alpha = 0.05, r = 4e3, lengths = 2^(floor(log2(length(y))):floor(log2(max(length(param$kern) + 1, 1 / param$param$cutoff)))), param, rm.out = FALSE, jumpint = confband, confband = FALSE)
{
  # compute standard deviation
  sdi <- sdrobnorm(y, lag = length(param$kern) + 1)

  # determine and remove outliers
  if(rm.out) {
    cutmin <- quantile(y, 1e-4) + qnorm(1 / length(y)) * sdi
    cutmax <- quantile(y, 1 - 1e-4) - qnorm(1 / length(y)) * sdi
    outlier <- y < cutmin | y > cutmax
    if(any(outlier)) {
      ycut <- y[-neighbours(which(outlier), 1:length(y), length(param$kern) - 1)]
      xcut <- x[-neighbours(which(outlier), 1:length(y), length(param$kern) - 1)]
    } else {
      ycut <- y
      xcut <- x
    }
  } else {
    ycut <- y
    xcut <- x
  }
  
  # compute quantile
  if(missing(q)) {
    if(is.null(r)) stop("q or r need to be specified!")
    q <- kMRC.quant(1 - alpha, length(y), r, param$kern, lengths)
  } else {
    alpha <- NA
  }
  
  # compute bounds
  bs <- bounds.MRC(ycut, q = q, lengths = lengths, family = 'gaussKern', param = param)
  
  # compute bounded solution
  allblocks <- stepbound.default(ycut, bounds = bs, x = xcut, x0 = x0, family = 'gaussKern', param = param, refit = FALSE, jumpint = jumpint, confband = confband)

  if(!is.na(alpha)) attr(allblocks, "alpha") <- alpha
  attr(allblocks, "q") <- q
  attr(allblocks, "sd") <- sdi
  allblocks
}
