#'@rdname lrt
#'@method LRT DRM
#'@export


LRT.DRM <-
function(object, splitcrit="score",...){
  
  if(is.character(splitcrit) && splitcrit == "score"){
    sc  <- rowSums(object$data)
    scm <- ifelse(sc > median(sc), 1,0)
  }
  else{
    if(!is.vector(splitcrit)){stop("Error: split criterium has to be a vector!", call. = FALSE)}
    scm <- splitcrit
  }  
  sp_dat <- split(as.data.frame(object$data), as.factor(scm), drop=FALSE)
  
  sp_res <- lapply(sp_dat, function(dat) DRM(dat))
  sp_val    <- sapply(sp_res, function(ex) ex$logLikelihood)
  sp_npar  <- sapply(sp_res, function(ex) length(ex$estpar))
  
  emp_Chi2 <- -2*(object$logLikelihood - sum(sp_val))
  df       <- sum(sp_npar) - length(object$estpar)
  pval     <- 1-pchisq(emp_Chi2, df)
  
  itempar_split <- sapply(sp_res, function(re) list(re$itempar*(-1)))
  itemse_split <- sapply(sp_res, function(re) list(re$itempar_se))
  
  res_lrt <- list(emp_Chi2 = emp_Chi2, df=df, pval=pval, itempar=itempar_split, item_se=itemse_split)
  class(res_lrt) <- "aLR"
  res_lrt
}
