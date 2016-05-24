TFMLazyScore <- function(mat, pvalue, bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                         type=c("PFM", "PWM"), granularity=1e-5){
  type <- match.arg(type)
  if(length(pvalue) != 1L){
    stop("pvalue must be length of 1")
  }
  if(pvalue > 1 || pvalue < 0){
    stop("pvalue must be between 0 and 1")
  }
  if(granularity <= 0){
    stop("granularity must be larger than 0")
  }
  bg <- normargPriorParams(bg)
  if(type == "PFM"){
    mat <- normargMat(mat)
  }
  score <- .Call("lazyScore", mat, pvalue, bg, type, granularity, 
                 PACKAGE="TFMPvalue")
  return(score)
}

