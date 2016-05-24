vario.func <- function (x, y, glob.mean, glob.sd, glob.N, is.multivar=FALSE,
                        type=c("semivar", "cov", "pearson", "spearman", "kendall", "moran", "geary")) {
  types=c("semivar", "cov", "pearson", "spearman", "kendall", "moran", "geary")
  type=match.arg(tolower(type), types)
  
  if (!is.multivar) {    
    if (type=="semivar")
      results=(1/(2*length(x)))*sum((x-y)^2)
    else if (type=="cov")
      results=(1/length(x))*(sum((x-glob.mean[1])*(y-glob.mean[2])))
    else if (type=="moran") {
      results=(1/length(x))*(sum((x-glob.mean[1])*(y-glob.mean[2])))/
        (glob.sd[1]*glob.sd[2])
    }
    else if (type=="geary") {
      results=((1/(2*length(x)))*(sum((x-y)^2)))/((1/(glob.N-1))*(glob.sd[1]*glob.sd[2]*glob.N))
    }
    else{
      stop("Error: variogram type must be one of semivar, cov, geary, or moran for univariate datasets")
    }
  }
  else {
    if (type=="cov") {
      results=cov(x, y, use = "pairwise.complete.obs")
    }
    else if (type %in% c("pearson", "spearman", "kendall")) {
      results=mean(cor(t(cbind(x, y)), method=type, use = "pairwise.complete.obs"))
    }
    else {
      stop("Error: variogram type must be cov, pearson, spearman, or kendall for multivariate datasets")
    }
  }
  return (as.numeric(results))
}
