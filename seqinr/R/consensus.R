consensus <- function(matali, method = c( "majority", "threshold", "IUPAC", "profile"), threshold = 0.60,
  warn.non.IUPAC = FALSE, type = c("DNA", "RNA")){

  if(inherits(matali, "alignment")) matali <- as.matrix(matali)
  if(!is.matrix(matali)) stop("matrix or alignment object expected")
  if(storage.mode(matali) != "character") stop("matrix of characters expected")
  
  method <- match.arg(method)
  
  if(method == "IUPAC"){
    type <- match.arg(type)
    res <- apply(matali, 2, bma, warn.non.IUPAC = warn.non.IUPAC, type = type)
    names(res) <- NULL
    return(res)
  }
  
  if(method == "majority"){
    majority <- function(x) names(which.max(table(x)))
    res <- apply(matali, 2, majority)
    names(res) <- NULL
    return(res)
  }
  
  if(method == "profile"){
    obsvalue <- levels(factor(matali))
    nrow <- length(obsvalue)
    row.names(matali)<-NULL
    res <- apply(matali, 2, function(x) table(factor(x, levels = obsvalue)))
    return(res)
  }
  
  if(method == "threshold"){
    profile <- consensus(matali, method = "profile")
    profile.rf <- apply(profile, 2, function(x) x/sum(x))
    res <- rownames(profile.rf)[apply(profile.rf, 2, which.max)]
    res <- ifelse(apply(profile.rf, 2, max) >= threshold, res, NA)
    names(res) <- NULL
    return(res)
  }
}

con <- consensus
