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
            .Object@pr <- 0
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

setMethod("cnCompare", c("catNetwork","catNetwork"),
          function(object1, object2, extended = FALSE) {
            return(.compare(object1, object2, extended))
          })

setMethod("cnCompare", c("catNetwork","matrix"),
          function(object1, object2, extended = FALSE) {
            return(.compare(object1, object2, extended))
          })

setMethod("cnCompare", c("catNetwork","list"),
          function(object1, object2, extended = FALSE) {
            return(findNetworkDistances(object1, as.integer(1), object2, extended))
          })

setMethod("cnCompare", c("catNetwork","catNetworkEvaluate"),
          function(object1, object2, extended = FALSE) {
            return(findNetworkDistances(object1, object2@numsamples, object2@nets, extended))
          })

.compare <- function(object1 ,object2, extended = FALSE) {

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
  out@pr <- 1
  if(out@tp+out@fp>0)
    out@pr <- out@tp/(out@tp+out@fp)
  out@sp <- 1
  if(tn+out@fp>0)
    out@sp <- tn/(tn+out@fp)
  out@sn <- 1
  if(out@tp+out@fn>0)
    out@sn <- out@tp/(out@tp+out@fn)
  out@fscore <- 0
    if(out@pr+out@sn > 0)
  out@fscore <- 2*out@pr*out@sn/(out@pr+out@sn)

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


findNetworkDistances <- function(object, numsamples, nets, extended = FALSE) {

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
    out@loglik[i] <- nets[[i]]@loglik

    out@tp[i] <- sum(mpartrue==1 & mpar==1)
    out@fn[i] <- sum(mpartrue == 1 & mpar == 0)
    out@fp[i] <- sum(mpartrue == 0 & mpar == 1)

    tn <- sum(mpartrue==0 & mpar==0)
    out@pr[i] <- 1
    if(out@tp[i]+out@fp[i]>0)
      out@pr[i] <- out@tp[i]/(out@tp[i]+out@fp[i])
    out@sp[i] <- 1
    if(tn+out@fp[i]>0)
      out@sp[i] <- tn/(tn+out@fp[i])
    out@sn[i] <- 1
    if(out@tp[i]+out@fn[i]>0)
      out@sn[i] <- out@tp[i]/(out@tp[i]+out@fn[i])
    out@fscore[i] <- 0
    if(out@pr[i]+out@sn[i] > 0)
      out@fscore[i] <- 2*out@pr[i]*out@sn[i]/(out@pr[i]+out@sn[i])
    
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


