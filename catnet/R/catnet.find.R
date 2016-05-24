#########################################################################
# Categorical Network Class Methods
# Find-members

setMethod("cnFind", "catNetworkEvaluate",
          function(object, complexity=0, alpha=0, factor=1) {
            if(!is(object, "catNetworkEvaluate"))
               stop("catNetworkEvaluate object should be specified.")
            if(length(object@nets) < 1)
              return(NULL)
            if(complexity>0) {
              complx <- sapply(object@nets, function(pnet) pnet@complexity)
              idx <- which(complx == complexity)
              if(length(idx) == 0) {
                for(i in 1:length(complx))
                  if(complx[i] > complexity)
                    break
                idx <- i
              }
              return(object@nets[[idx]])
            }
            if(alpha == "BIC") 
              alpha <- -1
            if(alpha == "AIC") 
              alpha <- -2
            if(factor <= 0) {
              factor <- 1
              warning("factor set to 1")
            }
            if(alpha==-2) {##AIC
              ff <- factor
              score <- sapply(object@nets, function(pnet) {                
                cmplx <- sapply(1:pnet@numnodes, function(node) as.integer(nodeComplexity(pnet, as.integer(node))))
                alpha.n <- ff/pnet@nodeSampleSizes
                return(pnet@likelihood - sum(alpha.n*cmplx))
              })
            }
            if(alpha==-1) {##BIC
              ff <- 0.5*factor
              score <- sapply(object@nets, function(pnet) {                
                cmplx <- sapply(1:pnet@numnodes, function(node) as.integer(nodeComplexity(pnet, as.integer(node))))
                alpha.n <- ff*log(pnet@nodeSampleSizes)/pnet@nodeSampleSizes
                return(pnet@likelihood - sum(alpha.n*cmplx))
              })
            }
            if(alpha!=-1&&alpha!=-2) {
              score <- sapply(object@nets, function(pnet) {
                ff <- factor/pnet@numnodes
                cmplx <- sapply(1:pnet@numnodes, function(node) as.integer(nodeComplexity(pnet, as.integer(node))))              
                alpha.n <- ff*exp(-alpha*log(pnet@nodeSampleSizes))
                return(pnet@likelihood - sum(alpha.n*cmplx))
              })
            }
            idx <- which(score == max(score))
            if(length(idx)!=1)
              return(NULL)
            return(object@nets[[idx]])
            })

setMethod("cnFind", "list",
          function(object, complexity=0, alpha=0, factor=1) {
            idx <- NULL
            idx <- sapply(object, function(pnet) {
              if(is(pnet, "catNetwork"))
                return(pnet@complexity==complexity)
              else
                return(FALSE)
            })
            if(is.null(idx))
              return(NULL)
            if(sum(idx) == 0) {
              idx <- sapply(object, function(pnet) {
                if(is(pnet, "catNetwork"))
                  return(abs(cnComplexity(pnet)-complexity))
              else
                return(Inf)
              })
              cc <- idx
              idx[cc==min(cc)] <- TRUE
              idx[cc!=min(cc)] <- FALSE              
            }
            idx <- which(idx==max(idx))
            id <- max(idx)
            return(object[[id]])
            })

setMethod("cnFindAIC", "catNetworkEvaluate", function(object) {
  objectlist <- object@nets
  if(length(objectlist) < 1)
    return(NULL)
  liststr <- ""
  maxobj <- objectlist[[1]]
  numsamples <- object@numsamples
  maxaic <- -Inf
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    curaic <- numsamples*object@likelihood - object@complexity
    if(maxaic < curaic) {
      maxaic <- curaic
      maxobj <- object
    }
  }
  return(maxobj)
})

setMethod("cnFindAIC", "list", function(object, numsamples) {
  if(length(object) < 1)
    return(NULL)
  numsamples <- as.integer(numsamples)
  if(numsamples < 1)
    stop("numsamples should be greater than 0")
  objectlist <- object
  liststr <- ""
  maxobj <- objectlist[[1]]
  maxaic <- -Inf
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    curaic <- numsamples*object@likelihood - object@complexity
    if(maxaic < curaic) {
      maxaic <- curaic
      maxobj <- object
    }
  }
  return(maxobj)
})

setMethod("cnFindBIC", "catNetworkEvaluate", function(object) {
  objectlist <- object@nets
  if(length(objectlist) < 1)
    return(NULL)
  liststr <- ""
  maxobj <- objectlist[[1]]
  maxbic <- -Inf
  numsamples <- object@numsamples
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    curbic <- numsamples*object@likelihood - 0.5*object@complexity*log(numsamples)
    if(maxbic < curbic) {
      maxbic <- curbic
      maxobj <- object
    }
  }
  return(maxobj)
})

setMethod("cnFindBIC", "list", function(object, numsamples) {
  if(length(object) < 1)
    return(NULL)
  numsamples <- as.integer(numsamples)
  if(numsamples < 1)
    stop("numsamples should be greater than 0")
  objectlist <- object
  liststr <- ""
  maxobj <- objectlist[[1]]
  maxbic <- -Inf
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    curbic <- numsamples*object@likelihood - 0.5*object@complexity*log(numsamples)
    if(maxbic < curbic) {
      maxbic <- curbic
      maxobj <- object
    }
  }
  return(maxobj)
})


setMethod("cnFindKL", "catNetworkEvaluate", function(object) {
  objectlist <- object@nets
  if(length(objectlist) < 1)
    return(NULL)
  liststr <- ""
  maxobj <- objectlist[[1]]
  maxbic <- -Inf
  numsamples <- object@numsamples
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    curbic <- numsamples*object@likelihood - cnKLComplexity(object)
    if(maxbic < curbic) {
      maxbic <- curbic
      maxobj <- object
    }
  }
  return(maxobj)
})

setMethod("cnFindKL", "list", function(object, numsamples) {
  if(length(object) < 1)
    return(NULL)
  numsamples <- as.integer(numsamples)
  if(numsamples < 1)
    stop("numsamples should be greater than 0")
  objectlist <- object
  liststr <- ""
  maxobj <- objectlist[[1]]
  maxbic <- -Inf
  for(object in objectlist) {
    if(!is(object, "catNetwork"))
      next
    curbic <- numsamples*object@likelihood - cnKLComplexity(object)
    if(maxbic < curbic) {
      maxbic <- curbic
      maxobj <- object
    }
  }
  return(maxobj)
})

