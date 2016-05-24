#########################################################################
# Categorical Network Class Methods  
 
## generate a catNetwork object 
setMethod("initialize", "catNetwork",  
          function(.Object, ...) { 
            if(!validObject(.Object)) 
              stop("Error creating the object.") 
            .Object@objectName <- "catNetwork" 
            return(.Object) 
          }) 
 
setMethod("initialize", "catNetwork",  
          function(.Object, numnodes, maxpars = 0, numcats = 2, p.default=FALSE, p.delta1=0.01, p.delta2=0.01, ... ) { 
 
            .Object@objectName <- "catNetwork" 
 
            if(nargs() >= 2 && is.numeric(numnodes) && is.numeric(maxpars) && is.numeric(numcats)) 
              return(genRandomCatnet(.Object, numnodes, maxpars, numcats, p.default, p.delta1, p.delta2)) 
 
            return(edges2catnet(.Object, c(1), edges = vector("list", 1))) 
          }) 
 
cnRandomCatnet <- function(numnodes, maxpars, numcats, p.delta1=0.01, p.delta2=0.01) { 
  return(new("catNetwork", as.integer(numnodes), as.integer(maxpars), as.integer(numcats),
             as.logical(FALSE), p.delta1, p.delta2 )) 
}
 
#cnCatnetFromGraph <- function(graph, numcats = 2, cats = NULL, probs = NULL) { 
#  maxcats <- numcats 
#  for(cat in cats) 
#    if(maxcats < length(cat)) 
#      maxcats <- length(cat) 
#  object <- new("catNetwork", numnodes=1, maxpars=0, maxcats=2) 
#  object <- graph2catnet(object, graph, as.integer(maxcats), cats, probs) 
#  return(object) 
#} 
 
cnCatnetFromEdges <- function(nodes, edges, numcats = 2) { 
  numnodes <- length(nodes) 
  if(numnodes < 1) 
    stop("'At least one node should be specified") 
  if(length(edges) != numnodes) 
    stop("'edges' should be a list with the length of 'nodes'") 
  if(numcats < 2) { 
    numcats <- 2 
    warning("'numcats' is set to 2") 
  } 
   
  maxpars <- 0 
  for(nodeedges in edges) 
    if(maxpars < length(nodeedges)) 
      maxpars <- length(nodeedges) 
   
  object <- new("catNetwork", numnodes=numnodes, as.numeric(maxpars), as.integer(numcats)) 
  object <- edges2catnet(object, nodes, edges, as.integer(numcats)) 
  return(object) 
} 
 
cnCatnetFromSif <- function(file, numcats = 2) { 
  nodes <- NULL
  fsif <- read.table(file)
  if(ncol(fsif) < 3)
    fsif <- read.csv(file,header=FALSE)
  if(ncol(fsif)<2)
    stop("file with at least 2 columns should be given")
  if(ncol(fsif) == 3) {
    col1 <- as.character(fsif[,1])
    col2 <- as.character(fsif[,3])
  }
  else {
    col1 <- as.character(fsif[,1])
    col2 <- as.character(fsif[,2])
  }
  for(s in as.character(col1)) {
    if(length(which(nodes==s)) > 0) 
      next 
    nodes <- c(nodes, s)
  }
  for(s in as.character(col2)) {
    if(length(which(nodes==s)) > 0) 
      next 
    nodes <- c(nodes, s)
  }
  numnodes <- length(nodes) 
  edges <- vector("list", numnodes) 
  for(n in 1:dim(fsif)[1]) { 
    k <- which(nodes==col1[n]) 
    m <- as.character(col2[n]) 
    if(length(which(edges[[k]]==m)) > 0) 
      next 
    edges[[k]] <- c(edges[[k]], m) 
  }

  object <- new("catNetwork", numnodes=length(nodes)) 
  object <- edges2catnet(object, nodes, edges, as.integer(numcats)) 
  return(object) 
} 
 
cnNew <- function(nodes, cats, pars, probs = NULL, p.delta1=0.01, p.delta2=0.01, dagonly = FALSE) {

  if(length(nodes) < 1) 
    stop("At least 1 node is needed") 
 
  if(length(nodes) != length(pars)) 
    stop("Incompatible parent list") 
 
  if(length(nodes) != length(cats)) 
    stop("Incompatible category list") 
   
  maxcats <- 0 
  for(cat in cats)
    if(!is.null(cat) && maxcats < length(cat))
      maxcats <- length(cat) 
 
  maxpars <- 0
  i <- 1
  for(par in pars) {
    if(!is.null(par) && length(par) > 0 && par != "") {
      par <- sapply(par, function(pp)
             if(is.character(pp))
               return (which(nodes == pp)[1])
             else
               return(as.integer(pp))
             )
      par <- par[par>=1 & par<= length(nodes)]
      pars[[i]] <- as.integer(par)
      if(maxpars < length(par))
        maxpars <- length(par)
    }
    i <- i + 1
  }
 
  object <- new("catNetwork", length(nodes)) 
  object@nodes <- nodes 
  object@pars <- pars 
  object@cats <- cats
  object@maxcats <- as.integer(maxcats)
  object@maxpars <- as.integer(maxpars)
  object@probs <- vector("list", length(nodes))

  object@nodecomplx <- sapply(1:object@numnodes, function(x) nodeComplexity(object, x)) 
  object@complx <- as.integer(sum(object@nodecomplx)) 

  if(!dagonly) { 
    if(!is.null(probs)) { 
      if(length(nodes) != length(probs)) 
        stop("Incompatible probability list") 
      object@probs <- probs 
    }
    else {
      object@probs <-
       lapply(seq(1, object@numnodes), function(parid) {
          if(length(object@pars[[parid]]) > 0)
            setRandomProb(parid, object@pars[[parid]], object@cats,
                        seq(1, length(object@pars[[parid]])), p.delta1, p.delta2)
          else
            setRandomProb(parid, object@pars[[parid]], object@cats, NULL, p.delta1, p.delta2)
        })
    }
  }
   
  object@loglik <- 0 
  object@nodelik <- NA
  object@nodeSampleSizes <- NA

  if(!dagonly)
    if(!validCatNetwork(object, FALSE)) 
      stop("Invalid network") 
  
  return(object)   
} 
 
setMethod("cnAddNode", c("catNetwork", "character", "character", "ANY", "ANY"), 
          function(object, node, cats, pars = NULL, prob = NULL) { 
 
            object@numnodes <- as.integer(object@numnodes + 1) 
             
            if(length(which(object@nodes == node))>0) { 
              node <- paste("N", object@numnodes ,sep="") 
            } 
            object@nodes <- c(object@nodes, node) 
             
            if(length(cats) < 2) { 
              warning("The node should have at least two cats") 
              cats <- paste("C", 1:2 ,sep="") 
            } 
            object@cats[[object@numnodes]] <- cats 
            if(length(cats) > object@maxcats) 
              object@maxcats <-  length(cats) 
 
            if(is.null(pars)) { 
              warning("Set random pars") 
              idx <- sample(1:(object@numnodes-1)) 
              pars <- idx[1:object@maxpars] 
            } 
            object@pars[[object@numnodes]] <- pars 
            if(length(pars) > object@maxpars) 
               object@maxpars <-  length(pars) 
 
            if(is.null(prob)) { 
              warning("Set default conditional probability") 
              prob <- setDefaultProb(object@numnodes, pars, object@cats, 1:length(pars)) 
            } 
            object@probs[[object@numnodes]] <- prob 
 
            return(object) 
            }) 
 
setMethod("show", "catNetwork", 
          function(object) { 
            if(is(object, "catNetwork")) 
              cat("A catNetwork object with ", object@numnodes, " nodes, ", 
                  object@maxpars, " pars, ", object@maxcats, 
                  " cats,\n Likelihood = ", object@loglik, 
                  ", Complexity = ", object@complx, ".\n") 
            }) 
 
setMethod("cnNumNodes", "catNetwork", function(object) length(object@nodes)) 
 
genRandomCatnet <- function(object, numnodes, maxpars, maxcats, defaultprob, p.delta1, p.delta2) { 
 
  nodes <- sapply(seq(1,numnodes), function(i) paste("N", i ,sep=""))             
  pars <- genRandomParents(numnodes, maxpars) 
   
  cats <- lapply(seq(1:numnodes), function(i, maxcats) 
                 paste("C", seq(1:maxcats), sep=""), 
                 maxcats) 

  if(p.delta1<0) p.delta1 <- 0
  if(p.delta2<0) p.delta2 <- 0
  while(p.delta1+p.delta2>=0.5) {
    p.delta1 <- p.delta1/2
    p.delta2 <- p.delta2/2
  }
  
  probs <- lapply(seq(1, numnodes), 
                  function(parid) { 
                    if(defaultprob) 
                      setDefaultProb(parid, pars[[parid]], cats, 
                                    seq(1, length(pars[[parid]]))) 
                    else 
                      setRandomProb(parid, pars[[parid]], cats, 
                                    seq(1, length(pars[[parid]])), p.delta1, p.delta2) 
                  }) 
 
  object@numnodes <- as.integer(numnodes) 
  object@nodes <- nodes 
  object@meta <- "" 
  object@maxpars <- as.integer(maxpars) 
  object@pars <- pars 
  object@cats <- cats 
  object@maxcats <- as.integer(maxcats) 
  object@probs <- probs 
 
  object@nodecomplx <- sapply(1:as.integer(numnodes), function(x) nodeComplexity(object, x)) 
  object@complx <- as.integer(sum(object@nodecomplx)) 
   
  object@loglik <- 0 
  object@nodelik <- NA
  object@nodeSampleSizes <- NA 
   
  return(object)   
} 
 
# returns vector   
setMethod("cnNodes", c("catNetwork", "missing"),  
          function(object) { 
            which <- seq(1:object@numnodes) 
            cnNodes(object, which) 
          }) 

setMethod("cnNodes", c("catNetwork", "vector"),  
          function(object, which) { 
            if(is.null(which)) 
              which <- seq(1:object@numnodes) 
            as.character(object@nodes[which]) 
          }) 
 
setMethod("cnEdges", c("catNetwork", "missing"),  
          function(object) { 
            which <- seq(1:object@numnodes) 
            cnEdges(object, which) 
          }) 

setMethod("cnEdges", c("catNetwork", "character"),  
          function(object, which) { 
            id <- sapply(which, function(node) which(object@nodes == node))
            if(length(id) < 1)
              return(NULL)
            cnEdges(object, id) 
          })

setMethod("cnEdges", c("catNetwork", "vector"),  
	function(object, which) { 
          if(object@numnodes < 1) 
            return(NULL) 
          if(is.null(which)) 
            which <- seq(1, object@numnodes) 
	  onodes <- cnNodes(object)   
	  plist <- vector("list", length(which)) 
          for(j in (1:length(which))) 
            names(plist)[[j]] <- onodes[which[j]] 
	  for(j in (1:object@numnodes)) { 
	    parlist <- object@pars[[j]]   
	    if(length(parlist) < 1) next 
	    for(i in (1:length(parlist))) {   
		idx <- which(which == parlist[i]) 
                if(length(idx) < 1) next 
                for(jj in 1:length(idx)) 
                  plist[[idx[jj]]] <- c(plist[[idx[jj]]], onodes[j]) # letter nodes    
	    }  
	  } 
          if(length(plist)==0) 
            return(plist) 
          i <- 1 
          while(i<=length(plist)) 
            if(is.null(plist[[i]])) 
              plist <- plist[-i] 
            else 
              i <- i + 1 
          return(plist) 
	})   
   
 
setMethod("cnMatEdges", "catNetwork",  
          function(object) { 
            mout <- NULL 
            if(object@numnodes < 1) 
              return(mout) 
            for(n in 1:object@numnodes) { 
              if(length(object@pars[[n]]) < 1) 
                next 
              for(j in object@pars[[n]]) { 
                mout <- c(mout, object@nodes[j], object@nodes[n]) 
              } 
            } 
            if(is.null(mout)) 
              return(NULL) 
            return(matrix(mout, ncol=2, byrow=TRUE)) 
          }) 
 
# returns list 
setMethod("cnParents", c("catNetwork", "missing"),  
          function(object) { 
            which <- seq(1:object@numnodes) 
            cnParents(object, which) 
          }) 

setMethod("cnParents", c("catNetwork", "character"),  
          function(object, which) { 
            id <- sapply(which, function(node) which(object@nodes == node))
            if(length(id) < 1)
              return(NULL)
            cnParents(object, id) 
          })

setMethod("cnParents", c("catNetwork", "vector"),  
          function(object, which) { 
            if(is.null(which)) 
              which <- seq(1, object@numnodes) 
            plist <- lapply(which, function(n) {
              n <- as.integer(n)
              if(n<1 || n>object@numnodes)
                return(NULL)
              if(is.null(object@pars[[n]])) 
                return(NULL) 
              as.character(object@nodes[object@pars[[n]]]) 
            })
            if(length(plist)==0) 
              return(plist)            
            if(!is.list(plist))
              return(as.vector(plist))
            pnames <- object@nodes[which]
            plist <- setNames(plist, pnames) 
            i <- 1 
            while(i<=length(plist)) 
              if(is.null(plist[[i]])) 
                plist <- plist[-i] 
              else 
                i <- i + 1 
            return(plist) 
          }) 
 
matParents <- function(object, nodeorder) { 
            n <- object@numnodes 
            mat <- matrix(rep(0, n*n), n, n) 
            for(j in 1:length(object@pars)) { 
              plist <- object@pars[[j]] 
              if(length(plist) == 0) 
                next 
              for(i in plist) { 
                mat[nodeorder[j], nodeorder[i]] <- 1 
              } 
            } 
            rownames(mat) <- object@nodes 
            colnames(mat) <- object@nodes 
            return(mat) 
            } 
 
## generates the binary matrix of pairwise node directed connections 
setMethod("cnMatParents", c("catNetwork", "missing"),  
          function(object) { 
            n <- object@numnodes 
            nodeorder <- seq(1,n) 
            return(matParents(object, nodeorder)) 
          })

setMethod("cnMatParents", c("catNetwork", "vector"),  
          function(object, nodeorder) { 
            n <- object@numnodes 
            if(missing(nodeorder) || is.null(nodeorder)) 
              nodeorder <- seq(1,n) 
            return(matParents(object, nodeorder)) 
          }) 
  
setMethod("cnOrder", "catNetwork",    
	function(object) { 
          cnOrderNodes(object@pars) 
	})  
 
setMethod("cnOrder", "list",    
	function(object) { 
          cnOrderNodes(object) 
	}) 
 
validCatNetwork <- function(obj, quietly=FALSE) {
  res = TRUE
  onodes <- obj@nodes   
  nnodes <- length(onodes)   
  opars <- obj@pars   
  omaxpars <- obj@maxpars   
  ocats <- obj@cats   
  omaxcats <- obj@maxcats   
  oprob <- obj@probs   

  if(nnodes < 1 || omaxcats < 2) {   
    if(nnodes < 1 && !quietly)   
      cat("At least one node is needed.\n")   
    if(omaxcats < 2 && !quietly)   
      cat("At least two cats are needed.\n")   
    res <- FALSE   
    return(res)   
  }   
  for(i in (1:length(onodes))) {   
    if(length(opars[[i]]) > omaxpars){   
      omaxpars <- length(opars[[i]])   
    }   
    if(length(ocats[[i]]) > omaxcats){   
      omaxcats <- length(ocats[[i]])   
    }   
    while(length(opars[[i]]) == 0 && i < length(opars))   
      i <- i + 1   
    if(length(opars[[i]]) == 0 && i >= length(opars))   
      break 
    nn <- 1 
    for(j in (1:length(opars[[i]]))) {   
      ipar <- opars[[i]][j]   
      if(ipar < 1 || ipar > nnodes) {   
         if(!quietly)   
           cat("Wrong parent set: node=", i,", parent=", j, ".\n")   
        res <- FALSE   
        return(res)   
      }   
      #Check DAG condition; Note that the parent set is an increasing sequence   
      if(ipar == i)  {   
        if(!quietly)   
          cat("Wrong network, not a DAG: self-reference for node ", i, ".\n")   
        res <- FALSE   
        return(res)   
      }   
      if(ipar < i && length(opars[[ipar]])) {   
        for(ii in (1:length(opars[[ipar]])) ) {   
          jj <- opars[[ipar]][ii]   
          if(jj == i) {   
            if(!quietly)   
              cat("Wrong network, not a DAG: nodes ", i, " and ", ipar, " are mutually connected.\n")   
            res <- FALSE   
            return(res)   
          }   
          if(jj > i)   
            break;   
        }   
      }   
      nn <- nn*length(ocats[[j]]) 
    } 
  } 
 
  for(j in 1:nnodes) { 
    if(!is.null(oprob[[j]]) && !checkProbSet(j, opars[[j]], ocats, seq(1,length(opars[[j]])), oprob[[j]])) { 
      if(!quietly)   
        cat("Wrong probability for node ", j, "\n")   
      res <- FALSE   
      return(res)   
    }   
  } 
 
  # update 
  obj@maxpars <- omaxpars   
  obj@maxcats <- omaxcats 

  if(!isDAG(onodes, opars)) {   
    if(!quietly)   
      cat("Wrong network: not a DAG.\n")   
    res <- FALSE   
    return(res)   
  } 
   
  return(res)    
}   
 
nodeComplexity <- function(object, nnode) { 
  ll <- sapply(object@pars[[nnode]], function(i) length(object@cats[[i]])) 
  if(length(ll)>0) 
    return(prod(ll)*(length(object@cats[[nnode]])-1)) 
  else 
    return(length(object@cats[[nnode]])-1) 
} 
 
parentSetComplexity <- function(numnodecats, pars, cats) { 
  ll <- sapply(pars, function(i) length(cats[[i]])) 
  if(length(ll)>0) 
    return(prod(ll)*(numnodecats-1)) 
  else 
    return(numnodecats-1) 
} 

condNonUnifComplexity <- function(idroot, ppars, pcatlist, probs, idx) {
  if(is.null(ppars) || length(idx) < 1) {
    nodep <- probs
    nodep <- nodep[nodep!=0]
    return(as.numeric(length(pcatlist[[idroot]])-1)) 
  }
  idnode <- ppars[idx[1]]
  if(length(idx)==1) {
    vv <- probs[[1]]
    for(cat in 2:length(pcatlist[[idnode]])) {
      vv <- vv + probs[[cat]]
    }
    mv <- vv / length(pcatlist[[idnode]])
    vv <- (probs[[1]]-mv)^2
    for(cat in 2:length(pcatlist[[idnode]])) {
      vv <- vv + (probs[[cat]]-mv)^2
    }
    if(sum(vv)<=1e-16)
      return(0)
    return(as.numeric(length(pcatlist[[idnode]])*(length(pcatlist[[idroot]])-1)))
  }
  ccmplx <- sapply(seq(1, length(pcatlist[[idnode]])),
         function(cat) condNonUnifComplexity(idroot, ppars, pcatlist, probs[[cat]], idx[-1]))
  return(sum(ccmplx))
}

setMethod("cnComplexity", signature("catNetwork"), function(object, node=NULL, include.unif=TRUE) {
  if(is.null(node)) { 
    pc <- sapply(1:object@numnodes, function(x) nodeComplexity(object, x)) 
    return(as.integer(sum(pc))) 
  }
  if(is.character(node))
    node <- sapply(node, function(str) which(object@nodes == str))
  cmplx <- 0
  if(!is.numeric(node)) { 
    pc <- sapply(1:object@numnodes, function(x) {
      if(include.unif)
        return(nodeComplexity(object, x))
      else {
        idx <- 1:length(object@pars[[x]])
        return(condNonUnifComplexity(x,object@pars[[x]], object@cats, object@probs[[x]], idx))
      }
    })
    cmplx <- as.integer(sum(pc))
  } 
  else {
    if(include.unif)
      cmplx <- sapply(node, function(i) as.integer(nodeComplexity(object, as.integer(i))))
    else {
      cmplx <- sapply(node, function(i) {
        idx <- 1:length(object@pars[[i]])
        condNonUnifComplexity(i,object@pars[[i]], object@cats, object@probs[[i]], idx)
      })
    }
  }
  return(cmplx)
}) 

setMethod("cnSubNetwork", "catNetwork",  
function(object, nodeIndices, indirectEdges = FALSE) { 

  if(is.character(nodeIndices)) {
    nodeIndices <- sapply(nodeIndices, function(node) which(object@nodes == node))
    if(length(nodeIndices) < 1)
      return(object)
  }
  
  nodeIndices <- nodeIndices[nodeIndices<=object@numnodes] 
  nodes <- object@nodes[nodeIndices] 
  numnodes <- length(nodes) 
  if(numnodes < 1) 
    stop("Can't create a network without nodes.") 
 
  pars <- vector("list", numnodes)   
  cats <- vector("list", numnodes) 
   
  maxpars <- 0 
  maxcats <- 0 
  for(j in 1:numnodes) { 
    i <- nodeIndices[j] 
    cats[[j]] <- object@cats[[i]] 
    if(indirectEdges==TRUE) { 
      nodepars <- object@pars[[i]] 
      while(!is.null(nodepars)) { 
        pool <- NULL 
        for(par in nodepars) { 
          id <-  which(nodeIndices==par) 
          if(length(id) > 0 && 
             ## prevent repeating 
             length(which(pars[[j]] == id)) <= 0) 
            pars[[j]] <- c(pars[[j]], id) 
          else 
            pool <- c(pool, object@pars[[par]]) 
        } 
        nodepars <- pool 
      } 
    } 
    else { 
      for(par in object@pars[[i]]) { 
        id <- which(nodeIndices==par) 
        if(length(id) > 0) 
          pars[[j]] <- c(pars[[j]], id) 
      } 
    } 
    if(maxpars < length(pars[[j]])) 
      maxpars <- length(pars[[j]]) 
    if(maxcats < length(cats[[j]])) 
      maxcats <- length(cats[[j]])
  } 

  newnet <- cnNew(nodes, cats, pars, probs = NULL, dagonly = TRUE)
  ##newnet <- new("catNetwork", numnodes, as.integer(maxpars), as.integer(maxcats)) 
  ##newnet@nodes <- nodes 
  ##newnet@pars <- pars 
  ##newnet@cats <- cats
  problist <- vector("list", numnodes)
  newnet@probs <- problist 

  pc <- sapply(1:newnet@numnodes, function(x) nodeComplexity(newnet, x)) 
  newnet@complx <- as.integer(sum(pc))
      
  return(newnet) 
}) 
 
 
setProbSlot <- function(idroot, ppars, pcatlist, idx, problist, catvec, probslot) { 
  if(is.null(ppars) || length(idx) < 1) { 
    if(length(problist) != length(probslot)) 
      stop("length(problist) != length(probslot)") 
    problist <- probslot 
    return(problist) 
  } 
  idnode <- ppars[idx[1]] 
  cat <- catvec[idx[1]] 
  if(cat < 1 || cat > length(pcatlist[[idnode]])) 
    stop("Wrong category\n") 
  return(setProbSlot(idnode, ppars, pcatlist, idx[-1], problist[[cat]], catvec, probslot)) 
} 
 
 
reorderNodeProb <- function(idroot, ppars, cats, catvec, idx, prob, 
                            nodeIndices, newparents, newcategories) { 
  if(is.null(ppars) || length(idx) < 1) { 
    newroot <- nodeIndices[idroot] 
    cat("newroot = ", newroot, "\n") 
    newidx <- NULL 
    if(length(newparents) > 0) 
      newidx <- seq(1, length(newparents)) 
    newprob <- initSampleProb(newroot, newparents, newcategories, newidx) 
    newprob <- setProbSlot(newroot, newparents, newcategories, 
                           newidx, newprob, catvec, prob) 
    return(newprob) 
  } 
  idnode <- ppars[idx[1]] 
  lapply(seq(1, length(cats[[idnode]])), 
         function(cat, idx) 
         reorderNodeProb(idroot, ppars, cats, c(catvec, cat), idx, prob[[cat]], 
                         nodeIndices, newparents, newcategories), 
         idx[-1]) 
} 
 
setMethod("cnReorderNodes", c("catNetwork", "vector"),  
function(object, nodeIndices) {

  if(is.character(nodeIndices)) {
    nodeIndices <- sapply(nodeIndices, function(node) which(object@nodes == node))
    if(length(nodeIndices) < 1)
      return(object)
  }

  nodeIndices <- nodeIndices[nodeIndices<=object@numnodes] 
  if(length(nodeIndices) != object@numnodes) { 
    warning("length(nodeIndices) != object@numnodes") 
    return(object) 
  } 
  nodeIndicesInvert <- nodeIndices 
  for(i in 1:object@numnodes) { 
    nodeIndicesInvert[i] = which(nodeIndices == i) 
  } 
   
  newnodes <- object@nodes[nodeIndices] 
  pars <- vector("list", object@numnodes) 
  cats <- vector("list", object@numnodes) 
  probs <- vector("list", object@numnodes) 
 
  for(i in 1:object@numnodes) { 
    if(length(object@pars[[nodeIndices[i]]]) > 0) {
      pars[[i]] <- nodeIndicesInvert[object@pars[[nodeIndices[i]]]]
      names(pars[[i]]) <- newnodes[pars[[i]]]
    }
    cats[[i]] <- object@cats[[nodeIndices[i]]] 
  }

  for(i in 1:object@numnodes) { 
    ##probs[[nodeIndices[i]]] <- reorderNodeProb(i, object@pars[[i]], object@cats, 
    ##                                      catvec=NULL, 
    ##                                      seq(1, length(object@pars[[i]])), 
    ##                                      object@probs[[i]], 
    ##                                      nodeIndices, pars, cats) 
    probs[[i]] <- object@probs[[nodeIndices[i]]] 
  } 
 
  object@nodes <- newnodes 
  object@cats <- cats 
  object@pars <- pars 
  object@probs <- probs 
   
  return(object) 
}) 

cnSetSeed <- function(seed) {
  ## set the seed both in R and in the C-library
  .Call("ccnSetSeed", seed, PACKAGE="sdnet")
  set.seed(seed)
}
