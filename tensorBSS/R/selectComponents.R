selectComponents <-
function(x, first=2, last=2){
  
  r <- length(dim(x)) - 1
  n <- dim(x)[r+1]
  dimprod <- prod(dim(x)[-(r+1)])
  
  xDev <- apply(x, 1:r, sd)
  xScaled <- sweep(x, 1:r, xDev^(-1), '*')
  xScaled_vec <- tensorVectorize(xScaled)

  kappa <- apply(xScaled_vec^4, 1, mean)
  indiceshi <- order(kappa, decreasing = TRUE)[1:first]
  indiceslo <- order(kappa, decreasing = FALSE)[1:last]

  indices <- rep(0, length(kappa))
  if(first >= 0.5){
    for(i in 1:first){
      indices[indiceshi[i]] <- 1
    }
  }
  if(last >= 0.5){
    for(i in 1:last){
      indices[indiceslo[i]] <- 1
    }
  }
  indexnames <- which(array(indices, dim=dim(x)[1:r]) != 0, arr.ind = T)
  # comps <- order(kappa, decreasing = TRUE)[c(1:first, (length(kappa)-last+1):length(kappa))]
  
  results <- NULL
  kappas <- NULL
  for(i in 1:nrow(indexnames)){
    results <- cbind(results, eval(parse(text=paste("x[", paste(as.character(indexnames[i, ]), sep=" ", collapse=","), ", ]", sep=""))))
    # Not optimal...
    kappas[i] <- mean(eval(parse(text=paste("xScaled[", paste(as.character(indexnames[i, ]), sep=" ", collapse=","), ", ]", sep="")))^4)
  }

  results <- results[, order(kappas, decreasing = TRUE)]
  colnames(results) <- apply(indexnames[order(kappas, decreasing = TRUE), , drop = FALSE], 1, function(a) paste(as.character(a), sep=" ", collapse=","))

  return(results)
}

