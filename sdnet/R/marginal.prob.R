#########################################################################
# Categorical Network Class Methods
# Marginal Probability Calculations

# Kullback-Leibler distance between two conditional probs at node [idroot]
# with pars set [ppars]
.nodeKLdist <- function(idroot, ppars, pcatlist, idx, problist1, problist2) {
  if(length(problist1) != length(problist2))
    stop("Incompatible probability lists while calling nodeKLdist.")
  if(is.null(ppars) || length(idx) < 1) {
    if(length(pcatlist[[idroot]]) < 2) {
       return(0)
    }
    probs <- problist1
    ps <- sum(probs)
    if(ps > 0)
      probs <- probs/ps
    probs[problist2==0] <- 0
    problist2[problist2==0] <- 1
    problist1[problist1==0] <- 1
    return(sum(probs*log(problist1/problist2)))
  }
  idnode <- ppars[idx[1]]
  poutlist <- sapply(seq(1,length(pcatlist[[idnode]])),
                     function(cat, idroot, ppars, pcatlist, idx, problist1, problist2)
                     .nodeKLdist(idroot, ppars, pcatlist, idx,
                                problist1[[cat]], problist2[[cat]]),
                     idroot, ppars, pcatlist, idx[-1], problist1, problist2
                     )
  return(sum(poutlist))
}

.nodeKLdistProbs <- function(idx, problist1, problist2) {
  if(length(problist1) != length(problist2))
    stop("Incompatible probability lists while calling nodeKLdist.")
  if(length(idx) < 1) {
    probs <- problist1
    ps <- sum(probs)
    if(ps > 0)
      probs <- probs/ps
    probs[problist2==0] <- 0
    problist2[problist2==0] <- 1
    problist1[problist1==0] <- 1
    return(sum(probs*log(problist1/problist2)))
  }
  poutlist <- sapply(1:length(problist1),
                     function(i, idx, problist1, problist2)
                     .nodeKLdistProbs(idx, problist1[[i]], problist2[[i]]),
                     idx[-1], problist1, problist2
                     )
  return(sum(poutlist))
}

.KLdistCatnets <- function(object1, object2) {
  ldist <- sapply(1:object1@numnodes, function(i)
         .nodeKLdist(i, object1@pars[[i]], object1@categories[[i]],
                    seq(1:length(object1@pars[[i]])),
                    object1@probs[[i]], object2@probs[[i]])
         )
  return(sum(ldist))
}

.KLdistProbs <- function(pars, problist1, problist2) {
  if(length(pars) < 1)
    return(0)
  ##cat(length(pars),"\n")
  ldist <- sapply(1:length(pars),
                  function(i, p1, p2) {
                    ##cat(i, length(pars[[i]]), length(problist2[[i]]),"\n")
                    if(length(pars[[i]]) > 0)
                      .nodeKLdistProbs(1:length(pars[[i]]), p1[[i]], p2[[i]])
                    else
                      return(0)
                  },
                  problist1, problist2
                  )
  return(sum(ldist))
}

setMethod("cnCondKLdist", "catNetwork",
          function(object1, object2,...) {
            if(!is(object1, "catNetwork") || !is(object2, "catNetwork"))
              stop("catNetwork object is required.")
            if(object1@numnodes != object2@numnodes)
              stop("Number of nodes should be equal.")
            blist <- sapply(1:object1@numnodes, function(i, o1, o2)
                            return(length(o1@pars[[i]]) == length(o2@pars[[i]]) &&
                                   sum(o1@pars[[i]] == o2@pars[[i]]) == length(o1@pars[[i]])),
                            object1, object2
                            )
            if(sum(blist) < object1@numnodes)
              stop("Incompatible parent sets.")
            return(.KLdistCatnets(object1, object2))
          })

#########################################################################

nodeAncestors <- function(idroot, ppars) {
  if(length(ppars[[idroot]]) < 1) {
    return(NULL)
  }
  ancset <- c(ppars[[idroot]])
  out <- lapply(ppars[[idroot]],
           function(i, ppars) {
               nodeAncestors(i, ppars)
             },
          ppars)
  if(!is.null(out))
    for(s in out)
      if(!is.null(s))
        ancset <- c(ancset, s)
  i <- 2
  while(i <= length(ancset)) {
    if(i < 2)
      next
    cc <- ancset[i]
    if(sum(ancset[1:(i-1)]==cc) > 0){
      ancset <- ancset[-i]
    }
    else
      i <- i + 1
  }
    
  return(ancset)
}

.probTreeAddLeaf <- function(ptree, leafproblist, leafpars, idtree, idleafpars, pcatlist) {
  if(length(idtree) < 1) {
    ##cat("idtree: ", idtree,"\n")
    ##cat("idleafpars: ", idleafpars,"\n")
    if(length(idleafpars)>0){
      cat(idleafpars)
      stop("Length(idleafpars) should be zero.")
    }
    ## tree is a scalar, leafproblist is a vector
    if(is.null(ptree))
       return(leafproblist)
    return(ptree*leafproblist)
  }
  treenode <- idtree[1]
  if(length(idleafpars) > 0)
    leafnode <- leafpars[idleafpars[1]]
  else
    leafnode <- 0
  ##cat("Nodes: ", idtree,"; ", leafnode, idleafpars, "\n")
  if(length(idtree) > 0 && treenode == leafnode)
    poutlist <- lapply(seq(1, length(pcatlist[[leafnode]])),
                       function(cat, ptree, leapproblist, leafpars, idtree, idleafpars, pcatlist)
                       .probTreeAddLeaf(ptree[[cat]], leafproblist[[cat]], leafpars, idtree, idleafpars, pcatlist), 
                       ptree, leafproblist, leafpars, idtree[-1], idleafpars[-1], pcatlist
                       )
  else
    poutlist <- lapply(seq(1, length(pcatlist[[treenode]])),
                       function(cat, ptree, leapproblist, leafpars, idtree, idleafpars, pcatlist)
                       .probTreeAddLeaf(ptree[[cat]], leafproblist, leafpars, idtree, idleafpars, pcatlist), 
                       ptree, leafproblist, leafpars, idtree[-1], idleafpars, pcatlist
                       )
  return(poutlist)
}


## 
.probTreeToMatrix <- function(ptree, idx, pcatlist, prob, offset) {
  if(length(idx) < 1) {
    prob[offset] <- ptree
    return(prob)
  }
  idnode <- idx[1]
  nodeoff <- 1
  for(i in 1:length(pcatlist)) ## number of nodes == length(pcatlist)
    if(i > idnode)
      nodeoff <- nodeoff*length(pcatlist[[i]])
  ##cat(idnode, nodeoff, length(pcatlist), "\n")
  for(cat in 1:length(pcatlist[[idnode]])) {
    off <- offset + nodeoff*(cat-1)
    prob <- .probTreeToMatrix(ptree[[cat]], idx[-1], pcatlist, prob, off)
  }
  return(prob)
}

.nodeMarginalProb <- function(idnode, pars, probs, categories) {

  if(idnode > length(pars) || is.null(pars[[idnode]])) {
    return(probs[[idnode]])
  }
    
  nodesOrder <- cnOrder(pars)

  nodepar <- nodeAncestors(idnode, pars)

  if(length(nodepar)>15){
    warning("Parent set is too big; returns default (0.5, 0.5).")
    return(rep(1/length(categories[[idnode]]), length(categories[[idnode]])))
  }

  ##cat("marginalNodeProb" , nodepar, "Order: ", nodesOrder, "\n")

  parorder <- rep(0, length(nodepar))
  for(i in 1:length(nodepar))
    parorder[i] <- which(nodesOrder==nodepar[i])
  ##cat(parorder, "\n")
  parorder <- order(parorder)
  ##cat(parorder, "\n")

  ptree <- NULL  
  idtree <- NULL

  ## find the joint probability of [nodepar,idnode]
  
  parordered <- rep(0, length(nodepar)+1)
  parcatlist <- vector("list", length(nodepar)+1)
  for(i in 1:length(nodepar)) {
    nnode <- nodepar[parorder[i]]
    parordered[i] <- nnode
    parcatlist[[i]] <- categories[[nnode]]
    
    idleafpars <- NULL
    if(length(pars[[nnode]]) > 0)
      idleafpars <- 1:length(pars[[nnode]])

    ##cat("add ", i, ":", nnode, " par: " ,pars[[nnode]], "\n")
    
    ptree <- .probTreeAddLeaf(ptree,
                              probs[[nnode]],
                              pars[[nnode]],
                              idtree,
                              idleafpars, 
                              categories)
    idtree <- c(idtree, nnode)
  }
  i <- i + 1
  nnode <- idnode
  parordered[i] <- nnode
  parcatlist[[i]] <- categories[[nnode]]
  idleafpars <- NULL
  if(length(pars[[nnode]]) > 0)
    idleafpars <- 1:length(pars[[nnode]])
  ptree <- .probTreeAddLeaf(ptree,
                           probs[[nnode]],
                           pars[[nnode]],
                           idtree,
                           idleafpars, 
                           categories)

  ##cat(parordered,"\n")
  ## list the joint prob in a matrix
  n <- 1
  for(i in 1:length(nodepar))
    n <- n*length(categories[[parordered[i]]])
  n <- n*length(categories[[idnode]])
  jointprob <- rep(0, n)
  jointprob <- .probTreeToMatrix(ptree, 1:(length(nodepar)+1), parcatlist, jointprob, 1)

  ncat <- length(categories[[idnode]])
  margprob <- rep(0, ncat)
  for(i in 1:ncat) {
    ids <- seq(i, length(jointprob), ncat)
    margprob[i] <- sum(jointprob[ids])
  }
  return(margprob)
}

setMethod("cnMarginalKLdist", signature("catNetwork", "catNetwork"), 
          function(object1, object2,...) {
            if(object1@numnodes != object2@numnodes)
              stop("Number of nodes should be equal.")
            dist <- 0
            for(i in 1:object1@numnodes) {
              p1 <- .nodeMarginalProb(i, object1@pars, object1@probs, object1@categories)
              p2 <- .nodeMarginalProb(i, object2@pars, object2@probs, object2@categories)
              probs <- p1
              probs[p2==0] <- 0
              p2[p2==0] <- 1
              p1[p1==0] <- 1
              dist <- dist + sum(probs*log(p1/p2))
                
            }
            return(dist)
          })

setMethod("cnNodeMarginalProb", signature("catNetwork"), 
          function(object, node) {
            if(is.character(node))
              node <- which(object@nodes == node)
            if(node < 1 || node > object@numnodes)
              stop("Invalid node ", node)
            vmarg <- .Call("ccnMarginalProb", 
                           object,
                           as.integer(node), 
                           PACKAGE="sdnet")
            return(vmarg)
            ##return(.nodeMarginalProb(nnode, object@pars, object@probs, object@categories))
          })

setMethod("cnMarginalKLdistList", signature("catNetwork", "list"),
          function(object1, object2list,...) {
            dist <- rep(0, length(object2list))

            cat("Marginal KL-distance...")
            
            p1list <- vector("list", object1@numnodes)
            for(i in 1:object1@numnodes)
              p1list[[i]] <- .nodeMarginalProb(i, object1@pars, object1@probs, object1@categories)
             
            for(j in 1:length(object2list)) {
              object2 <- object2list[[j]]
              if(is.null(object2))
                next
              if(object1@numnodes != object2@numnodes)
                stop("Number of nodes should be equal.")
              dist[j] <- 0
              for(i in 1:object1@numnodes) {
                p2 <- .nodeMarginalProb(i, object2@pars, object2@probs, object2@categories)
                p1 <- p1list[[i]]
                probs <- p1
                probs[p2==0] <- 0
                p2[p2==0] <- 1
                p1[p1==0] <- 1
                dist[j] <- dist[j] + sum(probs*log(p1/p2))
              }
            }

            cat("Done\n")
            return(dist)
          })
