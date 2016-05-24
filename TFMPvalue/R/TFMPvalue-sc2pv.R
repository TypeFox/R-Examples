

TFMsc2pv <- function(mat, score, bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                    type=c("PFM", "PWM")){
  ## TODD: valiate these mat, bg
  type <- match.arg(type)
  if(length(score) != 1L){
    stop("score must be length of 1")
  }
  bg <- normargPriorParams(bg)
  if(type == "PFM"){
    mat <- normargMat(mat)
  }
  pvalue <- .Call("sc2pv", mat, score, bg, type, PACKAGE="TFMPvalue")
  return(pvalue)
}

