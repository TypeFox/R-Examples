
#########################################################################
# Categorical Network Class Methods
# Random Sampling

##setMethod("cnSamples", c("catNetwork"),   
##	function(object, numsamples=1, output="frame", as.index=FALSE) {
##	  if(!is.numeric(numsamples) && !is.integer(numsamples))
##            stop("The number of samples should be integer")
##	  numsamples <- as.integer(numsamples)
##	  data <- sapply(1:numsamples, function(x) genRandomCats(object))
##	  data <- matrix(data, ncol=numsamples)
##          if(!as.index) {
##            for(i in 1:object@numnodes) {
##              cats <- object@cats[[i]]
##              for(j in 1:numsamples) {
##                data[i, j] <- cats[as.integer(data[i,j])]
##              }
##            }
##          }
##	  rownames(data)<-object@nodes
##	  if(!missing(output) && output=="matrix")
##            return(data)
##          return(as.data.frame(t(data)))
##	})

setMethod("cnSamples", c("catNetwork"),
	function(object, numsamples=1, pert=NULL, output="frame", as.index=FALSE, naRate=0) {
	  if(!is.numeric(numsamples) && !is.integer(numsamples))
            stop("The number of samples should be integer")
	  numsamples <- as.integer(numsamples)
          if(!is.null(pert)) {
            if(output=="frame" && !is.vector(pert))
              pert <- t(pert)
            if(is.vector(pert)) {
              pert <- sapply(1:numsamples, function(j) pert)
            }
            if(ncol(pert) != numsamples || nrow(pert) != object@numnodes)
              stop("The perturbation size should be ", object@numnodes, " by ", numsamples)
            for(i in 1:object@numnodes)
              if(is.character(pert[i,]))
                pert[i,] <- sapply(pert[i,], function(cc) {
                  id <- which(object@cats[[i]]==cc)
                  if(length(id) > 0)
                    return(id[1])
                  else
                    return(0)
                  })
          }

          if(naRate < 0)
            naRate <- 0
          if(naRate>1)
            naRate <- 1
          
          fast <- TRUE
          if(fast) {
            data <- .Call("ccnSamples", 
                          object, as.integer(numsamples), pert, as.numeric(naRate), 
                          PACKAGE="sdnet")
          }
          else {
            if(!is.null(pert))
              data <- sapply(1:numsamples, function(j) {
                genRandomCatsPert(object, pert[, j])
              })
            else
              data <- sapply(1:numsamples, function(j) genRandomCats(object))
          }

	  data <- matrix(data, ncol=numsamples)
          rownames(data)<-object@nodes
          if(!as.index) {
            for(i in 1:object@numnodes) {
              cats <- object@cats[[i]]
              for(j in 1:numsamples) {
                data[i, j] <- cats[as.integer(data[i,j])]
              }
            }
          }
          else {
            for(i in 1:object@numnodes)
                data[i, ] <- as.integer(data[i,])           
          }
          if(output!="matrix") 
            data <- as.data.frame(t(data))
          return(data)
	})

# recursive probability search
findProbSlot <- function(idroot, ppars, pcatlist, idx, problist, psample) {
  if(is.null(ppars) || length(idx) < 1) {
    return(problist)
  }
  idnode <- ppars[idx[1]]
  if(psample[idnode] < 1 || psample[idnode] > length(pcatlist[[idnode]]))
    stop("Wrong category\n")
  return(findProbSlot(idnode, ppars, pcatlist, idx[-1], problist[[psample[idnode]]], psample))
}

genRandomCats <- function(object) {  
  data <- rep(-1, length(object@nodes))  
  rpath <- orderNodesDescend(object@pars)
  for(i in (1:length(rpath))) {
    ##cat(i,", ")
    nnode <- rpath[i]
    if(length(object@cats[[nnode]]) < 1){
      data[nnode] <- 1
      next
    }
    if(is.null(object@pars[[nnode]])) {
      r <- runif(1,0,1)
      icat <- 1
      rcat <- object@probs[[nnode]][icat]
      while(rcat < r) {
        icat <- icat + 1
        rcat <- rcat + object@probs[[nnode]][icat]
      }
      data[nnode] <- icat
      next
    }
    plist <- findProbSlot(nnode, object@pars[[nnode]], object@cats,
                          seq(1,length(object@pars[[nnode]])), object@probs[[nnode]], data)
    r <- runif(1,0,1)
    icat <- 1
    rcat <- plist[icat]
    while(rcat < r) {
      icat <- icat + 1
      rcat <- rcat + plist[icat]
    }
    data[nnode] <- icat
    ##cat("\n")
  }
  return(data) 
}

genRandomCatsPert <- function(object, pert) {
  data <- rep(-1, length(object@nodes))  
  rpath <- orderNodesDescend(object@pars)
  for(i in (1:length(rpath))) {
    ##cat(i,", ")
    nnode <- rpath[i]
    if(pert[nnode] >= 1 && pert[nnode] <= length(object@cats[[nnode]])) {
      data[nnode] <- pert[nnode]
      next
    }
    if(length(object@cats[[nnode]]) < 1){
      data[nnode] <- 1
      next
    }
    if(is.null(object@pars[[nnode]])) {      
      r <- runif(1,0,1)
      icat <- 1
      rcat <- object@probs[[nnode]][icat]
      while(rcat < r) {
        icat <- icat + 1
        rcat <- rcat + object@probs[[nnode]][icat]
      }
      data[nnode] <- icat
      next
    }
    plist <- findProbSlot(nnode, object@pars[[nnode]], object@cats,
                   seq(1,length(object@pars[[nnode]])), object@probs[[nnode]], data)
    r <- runif(1,0,1)
    icat <- 1
    rcat <- plist[icat]
    while(rcat < r) {
      icat <- icat + 1
      rcat <- rcat + plist[icat]
    }
    data[nnode] <- icat
    ##cat("\n")
  }
  return(data) 
}

# creates a tree-list of zeroes
# for a node [idroot] with vector [ppars] of pars and [pcatlist] of cats
# idx is a recursion controlling index set
initSampleProb <- function(idroot, ppars, pcatlist, idx) {
  if(is.null(ppars) || length(idx) < 1) {
    return(rep(0, length(pcatlist[[idroot]])))
  }
  idnode <- ppars[idx[1]]
  lapply(seq(1, length(pcatlist[[idnode]])),
         function(cat, idroot, ppars, pcatlist, idx)
         initSampleProb(idroot, ppars, pcatlist, idx),
         idroot, ppars, pcatlist, idx[-1])
}

# for a node [idroot] with pars [ppars] and
# probability tree-list [problist], 
# update the conditional probability 
# based on a sample [psample] by incrementing the coresponding frequency value
updateSampleProb <- function(idroot, ppars, pcatlist, idx, problist, psample) {
  if(is.null(ppars) || length(idx) < 1) {
    if(is.na(psample[idroot]))
      return(problist)
    if(psample[idroot] > length(pcatlist[[idroot]]))
      stop("Wrong category\n")
    problist[psample[idroot]] <- problist[psample[idroot]] + 1
    return(problist)
  }
  idnode <- ppars[idx[1]]
  poutlist <- lapply(seq(1,length(pcatlist[[idnode]])),
           function(cat, idroot, ppars, pcatlist, idx, problist, psample) {
             if(is.na(psample[idnode]) || psample[idnode] != cat)
               return(problist[[cat]])
             else
               return(updateSampleProb(idroot, ppars, pcatlist, idx[-1], problist[[cat]], psample))
             },
           idroot, ppars, pcatlist, idx, problist, psample)
  return(poutlist)
}

# here [psample] contains the samples cats for the nodes listed in [ppars]
# while [pcounts] gives corresponding frequences
# update is for all cats simultaneously
setNodeSampleProb <- function(idroot, ppars, pcatlist, idx, problist, data) {
  if(is.null(ppars) || length(idx) < 1) {
    if(!is.null(dim(data)) && dim(data)[2] > 0)
      problist <- sapply(1:length(pcatlist[[idroot]]),
                         function(ncat) {                           
                           return(sum(data[idroot,]==ncat))
                         })
    else
      problist <- sapply(1:length(pcatlist[[idroot]]), function(cat) return(0))
    return(problist)
  }
  idnode <- ppars[idx[1]]
  poutlist <- lapply(1:length(pcatlist[[idnode]]),
      function(catidx, idroot, ppars, pcatlist, idx, problist, data, pcounts) {
      if(!is.null(dim(data)) && dim(data)[1] >= idnode){
        psubsample <- data[,data[idnode,] == catidx]
        ## beware: psubsample may bacome a vector
        if(is.null(dim(psubsample)) && length(psubsample) > 0)
          psubsample <- as.matrix(psubsample, rows=length(psubsample))
      }
      else
        psubsample <- NULL
      return(setNodeSampleProb(idroot, ppars, pcatlist, idx[-1], problist[[catidx]], psubsample))
    },
                     idroot, ppars, pcatlist, idx, problist, data)
  return(poutlist)
}

normalizeProb <- function(net) {
  if(!is(net,"catNetwork"))
    stop("net should be catNetwork.")
  if(length(net@pars) < net@numnodes)
      net@pars <- c(net@pars, vector("list", net@numnodes - length(net@pars)))
  if(length(net@pars) != length(net@probs))
     stop("length(net@pars) != length(net@probs)")
  for(i in 1:length(net@probs)) {
    idx <- NULL
    if(length(net@pars[[i]]) > 0)
      idx <- 1:length(net@pars[[i]])
    net@probs[[i]] <- normalizeProbSlot(i, net@pars[[i]], net@cats, idx, net@probs[[i]])
  }
  return(net)
}

# based on a sample [psample] with pars [ppars] and
# probability tree-list [problist], 
# normalize the conditional probability tree-list at [idroot] so that it sums to 1
normalizeProbSlot <- function(idroot, ppars, pcatlist, idx, problist) {
  if(is.null(ppars) || length(idx) < 1) {
    if(length(problist) != length(pcatlist[[idroot]]) || length(pcatlist[[idroot]]) < 1) {
       return(problist)
    }
    ps <- sum(problist)
    if(!is.na(ps) && is.numeric(ps)) {
      if(ps > 0)
        problist <- problist/ps
      else {
        ## set neutral probability
        nn <- length(problist)
        avg <- 1 / nn
        problist <- rep(avg, nn)
        if(nn > 1)
          problist[nn] <- 1 - sum(problist[1:(nn-1)])
      }
    }
    return(problist)
  }
  idnode <- ppars[idx[1]]
  poutlist <- lapply(seq(1,length(pcatlist[[idnode]])),
                     function(cat) {
                       normalizeProbSlot(idroot, ppars, pcatlist, idx[-1], problist[[cat]])
                     })
  return(poutlist)
}

# calculates the loglik at node [idroot] with pars [ppars] and
# conditional probability tree-list [problist]
# which is assumed to contain non-normalized frequences
# multinomial model is assumed 
nodeCondLoglik <- function(idroot, ppars, pcatlist, idx, problist) {
  if(is.null(ppars) || length(idx) < 1) {
    if(length(pcatlist[[idroot]]) < 1) {
       return(0)
    }
    probs <- problist
    ps <- sum(probs)

    ####################################################
    ## Next condition determines whether parent sets with
    ## non-sample-populated slots should be discarded
    ##if(ps <= 0)
    ##  return(-Inf)
    if(ps > 0) {
      ps <- 1 / ps
      probs <- probs * ps
    }
    ####################################################
    probs[probs==0] <- 1
    return(sum(problist*log(probs)))
  }
  idnode <- ppars[idx[1]]
  poutlist <- sapply(seq(1,length(pcatlist[[idnode]])),
                     function(cat, idroot, ppars, pcatlist, idx, problist)
                     nodeCondLoglik(idroot, ppars, pcatlist, idx[-1], problist[[cat]]),
                     idnode, ppars, pcatlist, idx, problist
                     )
  return(sum(poutlist))
}

.setSampleProb <- function(object, data) {
  if(!is(object, "catNetwork"))
    stop("Object should be catNetwork.")
  
  if(dim(data)[1] != object@numnodes)
    stop("Incompatible sample dimensions.\n")
  numsamples <- dim(data)[2]
  if(numsamples < 1)
    stop("No samples\n")

  for(nnode in (1:object@numnodes)) {
    object@probs[[nnode]] <- initSampleProb(nnode, 
                                                object@pars[[nnode]],
                                                object@cats,
                                                seq(1,length(object@pars[[nnode]])))
  }
  
  for(j in (1:numsamples)) {
    ps <- data[,j]
    for(nnode in 1:object@numnodes) {
      ## increment frequency      
      object@probs[[nnode]] <- updateSampleProb(nnode, object@pars[[nnode]], object@cats,
                                                        seq(1,length(object@pars[[nnode]])), object@probs[[nnode]], ps)
    }
    
    for(nnode in (1:object@numnodes)) {
      object@probs[[nnode]] <- normalizeProbSlot(nnode, object@pars[[nnode]], object@cats,
                                                         seq(1,length(object@pars[[nnode]])), object@probs[[nnode]])
    }
  }
  
  return(object)
}


setMethod("cnSetProb", "catNetwork",
          function(object, data, pert=NULL, nodeCats = NULL, softmode = FALSE) {
            
            if((!is.matrix(data) && !is.data.frame(data)) || nrow(data)*ncol(data)<1)
              stop("'data' should be a matrix or data frame of cats")

            if(is.data.frame(data)) {
              data <- as.matrix(t(data))
              if(!is.null(pert)) { 
                if(!is.data.frame(pert))
                  stop("Perturbations should be a data frame")
                pert <- as.matrix(t(pert))
              }
            }
            
            if(softmode) {
              if(!is.numeric(data) || is.integer(data))
                stop("Numeric data is expected in soft mode")
              if(is.null(nodeCats) || !is.list(nodeCats))
                stop("data is numeric; nodeCats should be specified")
              object@cats <- nodeCats
              object@maxcats <- max(sapply(nodeCats, function(cl) length(cl)))
              rownames <- names(nodeCats)
              numnodes <- length(nodeCats)
              numsamples <- ncol(data)
            }
            else  {
              ## call this if you want to keep objects' cats
              ##r <- .categorizeSample(data, pert, object)
              ## or that if you want to change them
              r <- .categorizeSample(data, pert, object=NULL, nodeCats=nodeCats, ask=TRUE)
              data <- r$data
              pert <- r$pert
            
              if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
                stop("The number of nodes in the object and data should be equal.")
              
              rownames <- rownames(data)
              if(length(rownames) != object@numnodes)
                stop("The data rows should be named after the nodes of the object.")

              ## objects' cats has to be reset for they might be different from the original if set by nodeCats
              object@cats <- r$cats
              object@maxcats <- r$maxcats
              
              numnodes <- dim(data)[1]
              numsamples <- dim(data)[2]
            }

            if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
              norder <- order(rownames)
              data <- data[norder,]
              pert <- pert[norder,]
              rownames <- rownames(data)
              norder <- order(object@nodes)
              object <- cnReorderNodes(object, norder)
            }

            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The row names should correspond to the object nodes.")

            ## we don't need this actually
            ## specify object's cats for the search
            catIndices <- NULL
            if(!is.null(nodeCats)) {
              ## the cats are explicitly given and may not come from the data
              catIndices <- lapply(1:object@numnodes, function(i) 1:length(object@cats[[i]]))
            }
            

	    if(TRUE) {
		## needs to create a probability list
		object@probs <- lapply(seq(1, numnodes),
                  function(parid) {
                      setDefaultProb(parid, object@pars[[parid]], object@cats,
                                    seq(1, length(object@pars[[parid]])))
                  })
		newobject <- .Call("ccnSetProb", 
                      	object, data, pert, 
                      	PACKAGE="sdnet")
                ## awkward but necessary
		newobject@nodes <- object@nodes
	    }
	    else {
              newobject <- .setSampleProb(object, data)
            }

            newobject@cats <- object@cats
            newobject@maxcats <- object@maxcats
            newobject@complx <- cnComplexity(object)
            names(newobject@probs) <- object@nodes
            
	    return(newobject)
          }
)

