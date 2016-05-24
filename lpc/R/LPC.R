LPC <- function(x,y, type="regression", soft.thresh=NULL, u=NULL,censoring.status=NULL){
  call <- match.call()
  dat <- list(x=x,y=y,censoring.status=censoring.status)
  CheckLPCFormat(dat,type,soft.thresh,u)
  if(type=="regression") scores <- quantitative.func(dat$x, dat$y, .05)$tt
  else if(type=="survival") scores <- cox.func(dat$x, dat$y, dat$censoring.status, .05)$tt
  else if(type=="two class") scores <- ttest.func(dat$x, dat$y, .05)$tt
  else if(type=="multiclass") scores <- multiclass.func(dat$x, dat$y, .05)$tt
  if(is.null(u))  u <- svd(t(scale(t(dat$x), center=TRUE, scale=FALSE)))$u
    coefs <- lm((scores-mean(scores)) ~ u+0)$coef
    if(is.null(soft.thresh)) soft.thresh <- GetSoftThresh(dat, u, type, .999*max(abs(coefs))) 
    lpc.coefs <- NULL
    for(i in 1:length(coefs)) lpc.coefs <- c(lpc.coefs, sign(coefs[i])*max(0, abs(coefs[i])-soft.thresh))
    lpcscores <- mean(scores) + u%*%lpc.coefs
  output <- list(lpcscores=lpcscores, tscores=scores, soft.thresh=soft.thresh, coefs=lpc.coefs,call=call)
  class(output) <- "lpcobj"
  return(output)
}  
