stablecouple <- function(mt, selfprefmethod = "customized", similarity = TRUE, 
                initselfpref = if(similarity) 1e-07 else median(mt[lower.tri(mt)]), 
                prob = 0.5){
  METHODS <- c("two-means", "mean", "quantile", "customized")
  selmethod <- pmatch(selfprefmethod, METHODS)
  if (is.na(selmethod)) stop("Invalid procedure to set the values of the main diagonal")
  if (selmethod == -1)  stop("Ambiguous procedure to set the values of the main diagonal")
  stopifnot(is.matrix(mt))
  stopifnot(isSymmetric(mt))
  stopifnot(is.logical(similarity))
  if(similarity) {fun <- "which.max"; subs <- ">="} else {fun <- "which.min"; subs <- "<="}
  nsp <- nrow(mt)
  if(!is.null(initselfpref)) {
    stopifnot(is.numeric(initselfpref) & is.vector(initselfpref))
    diag(mt) <- rep(initselfpref, length.out = nsp)
  }
  if(any(is.na(mt))) stop(paste("Please, check your matrix for misleading NA values.", 
                          "Replace them by numeric values"))
  if(selmethod == 3) {
    prob <- prob[1]
    if(prob > 1 | prob < 0) 
      stop("You must provide a value between 0 and 1 for the argument 'prob'")
  }     
  if(selmethod != 4) {
    i <- 0
    ref <- apply(mt, 1, function(x) {
                                     i <- i + 1
                                    return(x[do.call(subs, list(x, diag(mt)[i]))])
                                    }) 
  }
  kfast <- function(vec){
    refcut <- unique(vec)
    if(length(refcut)< 3) return(mean(refcut))
    aux <- unlist(tapply(vec, kmeans(vec, 2)[[1]], range))
    if(aux[2] < aux[3]) return(mean(aux[2:3]))
    else return(mean(aux[c(1, 4)])) 
  }
  if(selmethod == 1) diag(mt) <- unlist(lapply(ref, kfast))
  if(selmethod == 2) diag(mt) <- unlist(lapply(ref, mean)) 
  if(selmethod == 3) diag(mt) <- unlist(lapply(ref, quantile, probs = prob)) 
  selfpref <- diag(mt)
  valref <- stpairs <- array(NA, nsp)
  while(any(is.na(stpairs))){
    hi <- do.call(fun, list(mt))
    nf <- hi%%nsp
    nc <-  ceiling(hi/nsp)
    if(nf == 0) nf <- nsp
    stpairs[c(nf, nc)] <- c(nc, nf)
    valref[c(nf, nc)] <- mt[nf, nc]
    mt[c(nf,nc),] <- NA
    mt[,c(nf,nc)] <- NA
  }
  return(list(stpairs = stpairs, valref = valref, selfpref = selfpref))
}



