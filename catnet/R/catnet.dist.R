#########################################################################
# Categorical Network Class Methods
# Distance between networks

setMethod("initialize", "catNetworkDistance", 
          function(.Object, hamm=0, hammexp=0, tp = 0, fp=0, fn=0, order = 0, markov = 0, KLdist=0) {
            .Object@hamm <- hamm
            .Object@hammexp <- hammexp
            .Object@tp <- tp
            .Object@fp <- fp
            .Object@fn <- fn
            .Object@sp <- 0
            .Object@sn <- 0
            .Object@fscore <- 0
            .Object@skel.tp <- 0
            .Object@skel.fp <- 0
            .Object@skel.fn <- 0
            .Object@order.fp <- order
            .Object@order.fn <- order
            .Object@markov.fp <- markov
            .Object@markov.fn <- markov
            .Object@KLdist <- KLdist
            return(.Object)
            })

setMethod("show", "catNetworkDistance",
          function(object) {
            if(is(object, "catNetworkDistance"))
              if(object@fscore<1) 
                str <- sprintf(
" Edges:
                  TP = %d, 
                  FP = %d, 
	          FN = %d,
             F-score = %f, 
\n Hamming:
             (FP+FN) = %d, 
                 exp = %d,
\n Skeleton:
                  TP = %d, 
                  FP = %d, 
	          FN = %d,
\n Order:         
                  FP = %d,
                  FN = %d,
\n Markov blanket:
                  FP = %d, 
                  FN = %d \n",
                  object@tp, 
                  object@fp,
                  object@fn,
                  object@fscore, 
                  object@hamm,
                  object@hammexp,
                  object@skel.tp, 
                  object@skel.fp,
                  object@skel.fn,
                  object@order.fp,
                  object@order.fn,
                  object@markov.fp,
                  object@markov.fn)
            else
                str <- sprintf(
" Edges:
                  TP = %d, 
                  FP = %d, 
	          FN = %d,
             F-score = %f\n",
                  object@tp, 
                  object@fp,
                  object@fn,
                  object@fscore)
            cat(str, "\n")
            return(str)
            })

setMethod("initialize", "catNetworkEvaluate", 
          function(.Object, nnodes, numsamples, nnets) {
            .Object@numnodes <- nnodes
            .Object@numsamples <- numsamples
            .Object@nets <- vector("list", nnets)
            .Object@complexity <- rep(NA, nnets)
            .Object@loglik <- rep(NA, nnets)
            .Object@KLdist <- rep(NA, nnets)
            .Object@hamm <- rep(NA, nnets)
            .Object@hammexp <- rep(NA, nnets)
            .Object@tp <- rep(NA, nnets)
            .Object@fp <- rep(NA, nnets)
            .Object@fn <- rep(NA, nnets)
            .Object@sp <- rep(NA, nnets)
            .Object@sn <- rep(NA, nnets)
            .Object@fscore <- rep(NA, nnets)
            .Object@skel.tp <- rep(NA, nnets)
            .Object@skel.fp <- rep(NA, nnets)
            .Object@skel.fn <- rep(NA, nnets)
            .Object@order.fp <- rep(NA, nnets)
            .Object@order.fn <- rep(NA, nnets)
            .Object@markov.fp <- rep(NA, nnets)
            .Object@markov.fn <- rep(NA, nnets)
            .Object@time <- 0
            return(.Object)
            })

setMethod("show", "catNetworkEvaluate",
          function(object) {
            if(is(object, "catNetworkEvaluate"))
              str <- sprintf(
" Number of nodes    = %d, 
 Sample size        = %d,
 Number of networks = %d
 Processing time    = %.3f\n",                             
                  object@numnodes, 
                  object@numsamples,
                  length(object@nets),
                  object@time)
            cat(str, "\n")
            return(str)
            })

setMethod("cnPlot", "catNetworkEvaluate",
          function(object, file) {

            if(length(object@loglik) > 0 && length(object@complexity) > 0 &&
               (length(object@tp) == 0 || is.na(object@tp[1]))) {
              par(mfrow=c(1,1))
              plot(object@complexity, object@loglik, xlab="complexity", ylab="log(likelihood)", lty=1, 
                   main=paste(object@numsamples, " samples, ", object@numnodes, " nodes.", sep=""))
            }            
            else if(length(object@loglik) > 0 && length(object@complexity) > 0 &&
                    length(object@tp) > 0 && !is.na(object@tp[1]) && 
                    length(object@hamm) > 0 && !is.na(object@hamm[1]) &&
                    length(object@markov.fp) > 0 && !is.na(object@markov.fp[1]) && 
                    length(object@hammexp) > 0 && !is.na(object@hammexp[1])) {
              par(mfrow=c(3,2))
              plot(object@complexity, object@loglik, xlab="complexity", ylab="log(likelihood)", lty=1, 
                   main=paste(object@numsamples, " samples, ", object@numnodes, " nodes.", sep=""))
              ##plot(object@complexity, object@KLdist,
              ##     xlab="complexity", ylab="KL-dist", lty=1, 
              ##     main="Probability Distance")
              plot(object@complexity, object@tp,
                   xlab="complexity", ylab="TP", lty=1,
                   main="True Positives Directed Edges")
              plot(object@complexity, object@hamm,
                   xlab="complexity", ylab="Hamming", lty=1, 
                   main="Parent Matrix Distance")              
              plot(object@complexity, object@markov.fp + object@markov.fn,
                   xlab="complexity", ylab="Markov", lty=1, 
                   main="Markov Neighbor Distance")
              plot(object@complexity, object@fp,
                   xlab="complexity", ylab="FP", lty=1, 
                   main="False Positive Directed Edges")
              plot(object@complexity, object@fn,
                   xlab="complexity", ylab="FN", lty=1,
                   main="False Nagative Directed Edges")
            }
            else if(length(object@loglik) > 0 && length(object@complexity) > 0 &&
                    length(object@tp) > 0 && !is.na(object@tp[1]) &&
                    length(object@fp) > 0 && !is.na(object@fp[1]) && 
                    length(object@hamm) > 0 && !is.na(object@hamm[1]) ) {
              par(mfrow=c(2,2))
              xx <- object@complexity
              plot(xx, object@loglik, xlab="complexity", ylab="log(likelihood)", lty=1, 
                   main=paste(object@numsamples, " samples, ", object@numnodes, " nodes.", sep=""))
              plot(xx, object@hamm,
                   xlab="complexity", ylab="hamm", lty=1, 
                   main="Hamming Distance")
              plot(xx, object@tp,
                   xlab="complexity", ylab="TP", lty=1,
                   main="True Positive Directed Edges")
              plot(xx, object@fp,
                   xlab="complexity", ylab="FP", lty=1,
                   main="False Positive Directed Edges")
              ##plot(xx, object@skel.tp,
              ##     xlab="complexity", ylab="skeleton TP", lty=1,
              ##     main="True Positive Edges")
            }
            else if(length(object@nets) > 0) {
              loglik <- sapply(object@nets, function(net) net@likelihood)
              complx <- sapply(object@nets, function(net) net@complexity)
              par(mfrow=c(1,1))
              plot(complx, loglik, xlab="complexity", ylab="log(likelihood)", lty=1, 
                   main=paste(object@numsamples, " samples, ", object@numnodes, " nodes.", sep=""))              
            }
            
            })

setMethod("cnProcTime", "catNetworkEvaluate",
          function(object) {
          return(object@time)
        })

setMethod("cnCompare", c("catNetwork","catNetwork"),
          function(object1, object2, extended = TRUE) {
            return(.compare(object1, object2, extended))
          })

setMethod("cnCompare", c("catNetwork","matrix"),
          function(object1, object2, extended = TRUE) {
            return(.compare(object1, object2, extended))
          })

setMethod("cnCompare", c("catNetwork","list"),
          function(object1, object2, extended = TRUE) {
            return(findNetworkDistances(object1, as.integer(1), object2, extended))
          })

setMethod("cnCompare", c("catNetwork","catNetworkEvaluate"),
          function(object1, object2, extended = TRUE) {
            return(findNetworkDistances(object1, object2@numsamples, object2@nets, extended))
          })

.compare <- function(object1 ,object2, extended = TRUE) {

  if(!is(object1, "catNetwork"))
    stop("catNetwork should be specified")
  numnodes <- object1@numnodes

  if(is(object2, "catNetwork")) {
    if(object1@numnodes != object2@numnodes)
      stop("Only networks with the same number of nodes can be compared")
    if(prod(tolower(object1@nodes) == tolower(object2@nodes)) == 0) {
      norder <- order(cnNodes(object1))
      object1 <- cnReorderNodes(object1, norder)
      norder <- order(cnNodes(object2))
      object2 <- cnReorderNodes(object2, norder)
    }
    if(prod(tolower(object1@nodes) == tolower(object2@nodes)) == 0)
      stop("Only networks with the same nodes can be compared")
    numnodes2 <- object2@numnodes
    mpartrue <- cnMatParents(object1)
    mpar <- cnMatParents(object2)
  }
  else if(is.matrix(object2) && dim(object2)[1] == dim(object2)[2]) {
    numnodes2 <- dim(object2)[1]
    if(numnodes != numnodes2)
      stop("The number of nodes should be equal")
    if(!is.null(rownames(object2)) &&
       prod(tolower(object1@nodes) == tolower(rownames(object2))) == 0)
      stop("Only objects with the same nodes can be compared")
    mpartrue <- cnMatParents(object1)
    mpar <- matrix(as.integer(object2 > 0), nrow=numnodes)
  }
  else
    stop("No valid second object is specified")
  
  out <- new("catNetworkDistance")

  ## KL-distance
  ##out@KLdist <- cnMarginalKLdist(object1, object2)

  out@tp <- sum(mpartrue==1 & mpar==1)
  out@fn <- sum(mpartrue == 1 & mpar == 0)
  out@fp <- sum(mpartrue == 0 & mpar == 1)

  tn <- sum(mpartrue==0 & mpar==0)
  out@sp <- 1
  if(tn+out@fp>0)
    out@sp <- tn/(tn+out@fp)
  out@sn <- 1
  if(out@tp+out@fn>0)
    out@sn <- out@tp/(out@tp+out@fn)
  out@fscore <- 2*out@sp*out@sn/(out@sp+out@sn)

  out@skel.tp <- sum((t(mpartrue)+mpartrue)==1 & (t(mpar)+mpar)==1)/2
  out@skel.fn <- sum((t(mpartrue)+mpartrue) == 1 & (t(mpar)+mpar) == 0)/2
  out@skel.fp <- sum((t(mpartrue)+mpartrue) == 0 & (t(mpar)+mpar) == 1)/2
  
  out@hamm <- sum(mpartrue!=mpar)

  out@hammexp <- NA
  out@order.fp <- NA
  out@order.fn <- NA
  out@markov.fp <- NA
  out@markov.fn <- NA

  out@KLdist <- NA

  if(!extended)
	return(out)
  
  if(numnodes > 512) {
    cat("The network is large and calculating the exponential distance may be very slow. Continue? ('y' or 'n')\n")
    if(scan("", what="character", nmax=1, quiet=TRUE) != "y" ) 
      return(out)
  }
  
  m1 <- mpartrue
  sum1 <- m1
  m2 <- mpar
  sum2 <- m2
  ##fd <- 0
  k <- 1
  while((sum(m1) > 0 || sum(m2) > 0) && k <= numnodes) {
    k <- k+1
    ##fd <- fd + sum(m1 != m2)
    m1 <- m1%*%mpartrue
    sum1 <- sum1 + m1
    m2 <- m2%*%mpar
    sum2 <- sum2 + m2
  }
  out@hammexp <- sum(sum1!=sum2)

  out@order.fp <- sum((t(sum1) > 0) & (sum2 > 0))
  out@order.fn <- sum((sum1 > 0) & (t(sum2) > 0))
  ## it's better
  #out@order.fp <- sum((sum1 > 0) & (sum2 == 0))
  #out@order.fn <- sum((sum1 == 0) & (sum2 > 0))
  
  out@markov.fp <- 0
  out@markov.fn <- 0
  for(nnode in 1:numnodes) {
    for(n1 in 1:(numnodes-1)) {
      for(n2 in (n1+1):numnodes) {
        if(mpar[nnode, n1] == 1 && mpar[nnode, n2] == 1 &&
           (mpartrue[nnode, n1] != 1 || mpartrue[nnode, n2] != 1)) {
          out@markov.fp <- out@markov.fp + 1
        }
        if(mpartrue[nnode, n1] == 1 && mpartrue[nnode, n2] == 1 &&
           (mpar[nnode, n1] != 1 || mpar[nnode, n2] != 1)) {
          out@markov.fn <- out@markov.fn + 1
          ##cat(nnode, ": ", n1, ", ", n2, "\n")
        }
      }
    }
  }
  
  return(out)
}


findNetworkDistances <- function(object, numsamples, nets, extended = TRUE) {

  nnets <- length(nets)
  if(nnets==0)
    stop("No networks are given")
  for(i in 1:nnets)
    if(is.null(nets[[i]]))
      break
  nnets <- i
  out <- new("catNetworkEvaluate", object@numnodes, numsamples, nnets)
  out@nets <- nets

  numnodes <- object@numnodes
  
  ## KL-distance
  ##out@KLdist <- cnMarginalKLdistList(object, nets)

  norder <- order(cnNodes(object))
  if(prod(norder == seq(1, numnodes)) == 0)
    object <- cnReorderNodes(object, norder)
  
  mpartrue <- cnMatParents(object)
  m1 <- vector("list", numnodes)
  m1[[1]] <- mpartrue
  sum1 <- m1[[1]]
  k <- 1
  while(k < numnodes) {
    k <- k + 1
    m1[[k]] <- m1[[k-1]]%*%mpartrue
    sum1 <- sum1 + m1[[k]]
  }
  
  for(i in 1:nnets){

    ##cat(nets[[i]]@nodes, "\n")
    if(numnodes != nets[[i]]@numnodes)
      stop("numnodes != nets[[i]]@numnodes")

    if(prod(tolower(object@nodes) == tolower(nets[[i]]@nodes)) == 0) {      
      norder <- order(cnNodes(nets[[i]]))
      if(prod(norder == seq(1, numnodes)) == 0)
        nets[[i]] <- cnReorderNodes(nets[[i]], norder)
    }
    
    if(prod(tolower(object@nodes) == tolower(nets[[i]]@nodes)) == 0)
      stop("Only networks with the same nodes can be compared.")
    
    mpar <- cnMatParents(nets[[i]])
    out@complexity[i] <- nets[[i]]@complexity
    out@loglik[i] <- nets[[i]]@likelihood

    out@tp[i] <- sum(mpartrue==1 & mpar==1)
    out@fn[i] <- sum(mpartrue == 1 & mpar == 0)
    out@fp[i] <- sum(mpartrue == 0 & mpar == 1)

    tn <- sum(mpartrue==0 & mpar==0)
    out@sp[i] <- 1
    if(tn+out@fp[i]>0)
      out@sp[i] <- tn/(tn+out@fp[i])
    out@sn[i] <- 1
    if(out@tp[i]+out@fn[i]>0)
      out@sn[i] <- out@tp[i]/(out@tp[i]+out@fn[i])
    out@fscore[i] <- 2*out@sp[i]*out@sn[i]/(out@sp[i]+out@sn[i])
    
    out@skel.tp[i] <- sum((t(mpartrue)+mpartrue)==1 & (t(mpar)+mpar)==1)/2
    out@skel.fn[i] <- sum((t(mpartrue)+mpartrue) == 1 & (t(mpar)+mpar) == 0)/2
    out@skel.fp[i] <- sum((t(mpartrue)+mpartrue) == 0 & (t(mpar)+mpar) == 1)/2
    
    out@hamm[i] <- sum(mpartrue!=mpar)
  
    out@hammexp[i] <- NA
    out@order.fp[i] <- NA
    out@order.fn[i] <- NA
    out@markov.fp[i] <- NA
    out@markov.fn[i] <- NA

    out@KLdist[i] <- NA

    if(!extended)
	next

    m2 <- mpar
    sum2 <- m2
    ##fd <- sum(m1[[1]] != m2)
    k <- 1
    while(sum(m2) > 0 && k <= numnodes) {
      k <- k+1
      ##fd <- fd + sum(m1[[k]] != m2)
      m2 <- m2%*%mpar
      sum2 <- sum2 + m2
    }
    out@hammexp[i] <- sum(sum1!=sum2)

    out@order.fp[i] <- sum((t(sum1) > 0) & (sum2 > 0))
    out@order.fn[i] <- sum((sum1 > 0) & (t(sum2) > 0))
    ## it's better
    ##out@order.fp[i] <- sum((sum1 > 0) & (sum2 == 0))
    ##out@order.fn[i] <- sum((sum1 == 0) & (sum2 > 0))
    
    out@markov.fp[i] <- 0
    out@markov.fn[i] <- 0
    for(nnode in 1:numnodes) {
      for(n1 in 1:(numnodes-1)) {
        for(n2 in (n1+1):numnodes) {
          if(mpar[nnode, n1] == 1 && mpar[nnode, n2] == 1 &&
             (mpartrue[nnode, n1] != 1 || mpartrue[nnode, n2] != 1)) {
            out@markov.fp[i] <- out@markov.fp[i] + 1
          }
          if(mpartrue[nnode, n1] == 1 && mpartrue[nnode, n2] == 1 &&
             (mpar[nnode, n1] != 1 || mpar[nnode, n2] != 1))
            out@markov.fn[i] <- out@markov.fn[i] + 1
        }
      }
    }

    ##cat(out@markov.fn[i], "\n")
    
  }

  return(out)
}


setMethod("cnParHist", "catNetworkEvaluate",
          function(object) {

            numnodes <- object@numnodes

            bfirst <- TRUE
            mhist <- matrix(rep(0, numnodes*numnodes), nrow=numnodes)
            for(bnet in object@nets) {
              if(!is(bnet, "catNetwork"))
                next
              if(bnet@numnodes != numnodes)
                next
              if(bfirst) {
                rownames(mhist) <- bnet@nodes
                colnames(mhist) <- bnet@nodes
                bfirst <- FALSE
              }
              mhist <- mhist + cnMatParents(bnet)
            }
            
            return(mhist)
        })

setMethod("cnParHist", "list",
          function(object) {

            numnodes <- 0
            
            for(bnet in object) {
              if(!is(bnet, "catNetwork"))
                next
              if(numnodes < 1) {
                numnodes <- bnet@numnodes
                mhist <- matrix(rep(0, numnodes*numnodes), nrow=numnodes)
                rownames(mhist) <- bnet@nodes
                colnames(mhist) <- bnet@nodes
              }
              if(bnet@numnodes != numnodes)
                next
              mhist <- mhist + cnMatParents(bnet)
            }
            
            return(mhist)
        })
