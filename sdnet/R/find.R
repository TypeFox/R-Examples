#########################################################################
# Categorical Network Class Methods
# Find-members

setMethod("cnCatnetFromDagEvaluate", "dagEvaluate",
  function(object, index) {

  cnobj <- .Call("ccnCatnetFromDagEvaluate", 
                    object, index,
                    PACKAGE="sdnet")
  if(is.null(cnobj))
    return(NULL)
  return(cnobj)

  ncats <- as.integer(object@maxcats)
  npars <- as.integer(object@maxpars)
  ## create catNetwork without prob table
  cnobj <- new("catNetwork", as.integer(object@numnodes), npars, ncats, pDefault=NULL)
  cnobj@nodes <- object@nodes
  cnobj@pars <- vector("list", as.integer(object@numnodes))
  objpars <- as.integer(object@numPars[[index]])
  cbuff <- NULL
  k <- 1
  while(k <= length(objpars)) {
    dn <- as.numeric(objpars[k])
    if(dn < 0) 
      dn <- as.numeric(0x100000000) + as.numeric(dn)
    n1 <- as.integer(dn/as.integer(0x1000000))
    dn <- dn-n1*as.numeric(0x1000000)
    n2 <- as.integer(dn/as.integer(0x10000))
    dn <- dn-n2*as.numeric(0x10000)
    n3 <- as.integer(dn/as.integer(0x100))
    n4 <- dn-n3*as.numeric(0x100)
    cbuff <- c(cbuff,n4)
    cbuff <- c(cbuff,n3)
    cbuff <- c(cbuff,n2)
    cbuff <- c(cbuff,n1)
    k <- k + 1
  }
  numpars <- NULL
  k <- 1
  i <- 1
  while(k <= length(cbuff)) {
    n <- cbuff[k]
    if(n >= 128) {
      n <- 256 - n
      k <- k + 1
      numpars <- c(numpars, rep(cbuff[k], n))
    }
    else
      numpars <- c(numpars, n)
    k <- k + 1
  }
  numpars <- numpars[1:object@numnodes]
  for(i in 1:object@numnodes) { 
    nodeParSize <- numpars[i]
    if(nodeParSize > 0) {
      cnobj@pars[[i]] <- object@parSlots[[i]][((nodeParSize-1)*npars+1):((nodeParSize-1)*npars+nodeParSize)]
      names(cnobj@pars[[i]]) <- object@nodes[cnobj@pars[[i]]]
    }
    cnobj@cats[[i]] <- 1:ncats 
    cnobj@nodecomplx[i] <- object@parComplx[[i]][nodeParSize+1]
    cnobj@nodelik[i] <- object@parLogliks[[i]][nodeParSize+1]
  }
  cnobj@complx <- object@complx[index]
  cnobj@loglik <- object@loglik[index]
  cnobj@nodeSampleSizes <- NA
  return(cnobj)
})

setMethod("cnFind", "dagEvaluate", 
          function(object, complx=0, alpha=0, factor=1) {
            if(complx>0) {
              idx <- which(object@complx == complx)
              if(length(idx) == 0) {
                for(i in 1:length(object@complx))
                  if(object@complx[i] > complx)
                    break
                idx <- i
              }
              return(cnCatnetFromDagEvaluate(object, idx))
            }
            if(alpha == "BIC") 
              alpha <- -1
            if(alpha == "AIC") 
              alpha <- -2
            if(factor <= 0) {
              factor <- 1
              warning("factor set to 1")
            }

            nets <- lapply(1:object@numDags, function(idx) {
              return(cnCatnetFromDagEvaluate(object, idx))
            })
            
            if(alpha==-2) {##AIC
              score <- sapply(nets, function(pnet) {                
                cmplx <- sapply(1:pnet@numnodes, function(node) as.integer(nodeComplexity(pnet, as.integer(node))))
                alpha.n <- 1/pnet@nodeSampleSizes
                return(pnet@loglik - sum(alpha.n*cmplx))
              })
            }
            if(alpha==-1) {##BIC
              score <- sapply(nets, function(pnet) {                
                cmplx <- sapply(1:pnet@numnodes, function(node) as.integer(nodeComplexity(pnet, as.integer(node))))
                alpha.n <- 0.5*log(pnet@nodeSampleSizes)/pnet@nodeSampleSizes
                return(pnet@loglik - sum(alpha.n*cmplx))
              })
            }
            if(alpha!=-1&&alpha!=-2) {
              score <- sapply(nets, function(pnet) {
                cmplx <- sapply(1:pnet@numnodes, function(node) as.integer(nodeComplexity(pnet, as.integer(node))))              
                alpha.n <- factor*exp(-alpha*log(pnet@nodeSampleSizes))
                return(pnet@loglik - sum(alpha.n*cmplx))
              })
            }
            idx <- which(score == max(score))
            if(length(idx)!=1)
              return(NULL)
            return(nets[[idx]])
          })

setMethod("cnFind", "catNetworkEvaluate", 
          function(object, complx=0, alpha=0, factor=1) {
            if(length(object@nets) < 1)
              return(NULL)
            if(complx>0) {
              netcomplx <- sapply(object@nets, function(pnet) pnet@complx)
              idx <- which(netcomplx == complx)
              if(length(idx) == 0) {
                for(i in 1:length(netcomplx))
                  if(netcomplx[i] > complx)
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
              score <- sapply(object@nets, function(pnet) {                
                cmplx <- sapply(1:pnet@numnodes, function(node) as.integer(nodeComplexity(pnet, as.integer(node))))
                alpha.n <- 1/pnet@nodeSampleSizes
                return((pnet@loglik - sum(alpha.n*cmplx))/pnet@numnodes)
              })
            }
            if(alpha==-1) {##BIC
              score <- sapply(object@nets, function(pnet) {                
                cmplx <- sapply(1:pnet@numnodes, function(node) as.integer(nodeComplexity(pnet, as.integer(node))))
                alpha.n <- 0.5*log(pnet@nodeSampleSizes)/pnet@nodeSampleSizes
                return((pnet@loglik - sum(alpha.n*cmplx))/pnet@numnodes)
              })
            }
            if(alpha!=-1&&alpha!=-2) {
              score <- sapply(object@nets, function(pnet) {
                cmplx <- sapply(1:pnet@numnodes, function(node) as.integer(nodeComplexity(pnet, as.integer(node))))              
                alpha.n <- factor*exp(-alpha*log(pnet@nodeSampleSizes))
                return((pnet@loglik - sum(alpha.n*cmplx))/pnet@numnodes)
              })
            }
            idx <- which(score == max(score))
            if(length(idx)!=1)
              return(NULL)
            return(object@nets[[idx]])
            })

setMethod("cnFind", "list",
          function(object, complx=0, alpha=0, factor=1) {
            idx <- NULL
            idx <- sapply(object, function(pnet) {
              if(is(pnet, "catNetwork"))
                return(pnet@complx==complx)
              else
                stop("List of catNetwork's is expected")
            })
            if(is.null(idx))
              return(NULL)
            if(sum(idx) == 0) {
              idx <- sapply(object, function(pnet) {
                if(is(pnet, "catNetwork"))
                  return(abs(cnComplexity(pnet)-complx))
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
    curaic <- (numsamples*object@loglik - object@complx)/object@numnodes
    if(maxaic < curaic) {
      maxaic <- curaic
      maxobj <- object
    }
  }
  return(maxobj)
})

setMethod("cnFindAIC", "dagEvaluate", function(object) {
  if(object@numDags < 1)
    return(NULL)
  maxk <- -1
  maxaic <- -Inf
  numsamples <- object@numsamples
  for(k in 1:object@numDags) {
    curaic <- numsamples*object@loglik[k] - object@complx[k]
    if(maxaic < curaic) {
      maxaic <- curaic
      maxk <- k
    }
  }
  cnobj <- cnCatnetFromDagEvaluate(object, maxk)
  return(cnobj)
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
    curaic <- (numsamples*object@loglik - object@complx)/object@numnodes
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
    curbic <- (numsamples*object@loglik - 0.5*object@complx*log(numsamples))/object@numnodes
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
    curbic <- (numsamples*object@loglik - 0.5*object@complx*log(numsamples))/object@numnodes
    if(maxbic < curbic) {
      maxbic <- curbic
      maxobj <- object
    }
  }
  return(maxobj)
})

setMethod("cnFindBIC", "dagEvaluate", function(object) {
  if(object@numDags < 1)
    return(NULL)
  maxk <- -1
  maxbic <- -Inf
  numsamples <- object@numsamples
  logNumSamples <- log(object@numsamples)
  for(k in 1:object@numDags) {
    curbic <- numsamples*object@loglik[k] - 0.5*object@complx[k]*logNumSamples
    if(maxbic < curbic) {
      maxbic <- curbic
      maxk <- k
    }
  }
  cnobj <- cnCatnetFromDagEvaluate(object, maxk)
  return(cnobj)
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
    curbic <- numsamples*object@loglik - cnKLComplexity(object)
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
    curbic <- numsamples*object@loglik - cnKLComplexity(object)
    if(maxbic < curbic) {
      maxbic <- curbic
      maxobj <- object
    }
  }
  return(maxobj)
})

