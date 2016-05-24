#########################################################################
# Directed Acyclic Graph Routines

# generate a DAG pars list
genRandomParents<- function(numnodes, maxparents, fixedparsize = 0) {
  idx <- sample(seq(1, numnodes))
  pars <- vector("list", numnodes)
  for(i in (2:numnodes)) {
    if(fixedparsize)
      npars <- 2*maxparents
    else
      npars <- floor(2*maxparents*runif(1,0,1) + 0.5)
    if(npars <= maxparents)
      next;
    npars <- npars - maxparents
    if(npars > i - 1)
      npars <- i - 1
    pars[[idx[i]]] <- idx[sort(sample(1:(i-1), npars))]
  }
  ## very important
  ## make it compatible with nodesOrder function
  norder <- orderNodesDescend(pars)
  for(i in norder) {
    if(length(pars[[i]]) < 2)
      next;
    for(j in 2:length(pars[[i]]))
      for(k in 1:(j-1))
        if(which(norder == pars[[i]][k]) > which(norder == pars[[i]][j]) ) {
          x <- pars[[i]][k]
          pars[[i]][k] <- pars[[i]][j]
          pars[[i]][j] <- x
        }
  }
  return(pars)
}


cnOrderParentList <- function(lpars) findRootPath(lpars)
  
findRootPath <- function(lpars) {
  if(!is.list(lpars))
    stop("Specify a list of parent sets.")
  rpath <- NULL  
  nn <- length(lpars)  
  rflags <- rep(0, nn)  
  brepeat <- TRUE  
  while(brepeat) {  
    iroot <- 0  
    for(i in (1:nn)) {  
      # find a root node  
      if(rflags[i] == 0) {  
        if(i > length(lpars)) {  
          broot <- TRUE  
        }  
        else {  
          pp <- lpars[[i]]  
          broot <- TRUE  
          if(length(pp) > 0) {  
            for(j in (1:length(pp))) {  
              #cat(i, j,pp[j],rflags[pp[j]], "\n")  
              if(pp[j] > 0 && rflags[pp[j]] == 0){  
                broot <- FALSE  
                break  
              }  
            }  
          }  
        }  
        #cat(broot,"\n")  
        if(broot){  
          iroot <- i  
          break  
        }  
      }  
    }      
    if(iroot > 0) {  
      rflags[iroot] <- 1  
	rpath <- c(rpath, iroot)  
    } else  
      break  
    if(sum(rflags) == nn)  
      brepeat <- FALSE  
  }  
  return(rpath)  
} 


cnOrderNodes <- function(lpars) orderNodesDescend(lpars)

# extract a topological order of a DAG
# [lpars] is a list of parent sets
orderNodesDescend <- function(lpars) {
  if(!is.list(lpars))
    stop("Specify a list of parent sets.")
  ## remove 2-loops from CPDAGs
  lpars <- lapply(1:length(lpars), function(i) {
    if(length(lpars[[i]]) < 1)
      return(NULL)
    pp <- NULL
    for(np in lpars[[i]]) {
      if(length(lpars[[np]])<1 || sum(which(lpars[[np]]==i)) < 1)
        pp <- c(pp, np)
    }
    return(pp)
  })
  
  rpath <- NULL  
  nn <- length(lpars)  
  rflags <- rep(0, nn)
  brepeat <- TRUE  
  while(brepeat) {
    brepeat <- FALSE
    for(i in (1:nn)) {
      # find a root node  
      if(rflags[i] == 0) {  
        pp <- lpars[[i]]
        ##cat(pp, "\n")
        broot <- TRUE  
        if(length(pp) > 0) {  
          for(j in (1:length(pp))) {  
            ##cat(i, j, pp[j], rflags[pp[j]], "\n")  
            if(pp[j] > 0 && rflags[pp[j]] == 0){  
              broot <- FALSE  
              break  
            }  
          }  
        }  
        if(broot){  
          rpath <- c(rpath, i)
          rflags[i] <- 1
          brepeat <- TRUE
        }  
      }  
    }
  }
  return(rpath)  
}  

  
isDAG <- function(lnodes, lpars) {  
  nn <- length(lnodes)  
  rflags <- rep(0, nn)  
  brepeat <- TRUE  
  while(brepeat) {  
    iroot <- 0  
    for(i in (1:nn)) {  
      # find a root node  
      if(rflags[i] == 0) {  
        if(i > length(lpars)) {  
          broot <- TRUE  
        }  
        else {  
          pp <- lpars[[i]]  
          broot <- TRUE  
          if(length(pp) > 0) {  
            for(j in (1:length(pp))) {  
              
              if(is(pp[j], "character"))
                ppj <- which(lnodes==pp[j])
              else
                ppj <- pp[j]
              if(ppj > 0 && rflags[ppj] == 0){  
                broot <- FALSE  
                break  
              }  
            }  
          }  
        }  
        if(broot){  
          iroot <- i  
          break  
        }
      }
    }
    if(iroot > 0)  
      rflags[iroot] <- 1  
    else  
      break  
    if(sum(rflags) == nn)  
      brepeat <- FALSE  
  }  
  return(!brepeat)  
}

############################################################################
## converting DAG to CPDAG; equivalence class reporesentation
## ORDER-EDGES and LABEL-EDGES algorithms as given in Chickering ...

orderEdges <- function(object) {
  if(!is(object, "catNetwork"))
    stop("catNetwork object is expected.")
  nodeorder <- orderNodesDescend(object@pars)
  #cat(nodeorder)
  order <- sapply(1:object@numnodes, function(k, pars) {
    rep(0, length(pars[[k]]))
    }, object@pars)
  i <- 1
  while(1) {
    ## find the lowest ordered node y that has an unordered edge incident into it
    y <- 0
    for(k in 1:object@numnodes) {
      if(length(order[[k]]) == 0)
        next
      for(kk in 1:length(order[[k]]))
        if(order[[k]][kk] == 0)
          ## edge: pars[[k]][kk] -> k
          if(y <= 0 || nodeorder[y] > nodeorder[k])
            y <- k
    }
    #cat("y = ", y,"\n")
    if(y <= 0)
      break
    ## find the highest odered node x for which with edge x->y not ordered 
    x <- 0
    kkx <- 0
    for(kk in 1:length(order[[y]]))
      if(order[[y]][kk] == 0)
        if(x <= 0 || nodeorder[x] < nodeorder[object@pars[[y]][kk]]) {
          x <- object@pars[[y]][kk]
          kkx <- kk
        }
    if(x <= 0)
      break
    #cat("order ", y, kkx, i, "\n")
    order[[y]][kkx] <- i
    i <- i + 1
  }
  return(order)
}

labelEdges <- function(object) {
  if(!is(object, "catNetwork"))
    stop("catNetwork object is expected.")
  order <- orderEdges(object)
  label <- sapply(1:object@numnodes, function(k, pars) {
    rep(-1, length(pars[[k]]))
    }, object@pars)
  ## label = 1 for compelled and 0 for reversible edge
  while(1) {
    ## find the lowest ordered edge x->y that is unknown (-1)
    miny <- 0
    minx <- 0
    for(k in 1:object@numnodes) {
      if(length(order[[k]]) == 0)
        next
      for(kk in 1:length(order[[k]]))
        if(label[[k]][kk] == -1)
          if(miny <= 0 || order[[miny]][minx] > order[[k]][kk]) {
            miny <- k
            minx <- kk
          }
    }
    if(miny <= 0)
      break
    y <- miny
    x <- object@pars[[miny]][minx]
    ##cat(x, "->", y, minx, miny, "\n")
    bnextcycle <- FALSE
    ## check all compelled edges w->x
    if(length(label[[x]]) > 0) {
      for(w in 1:length(label[[x]]))
        if(label[[x]][w] == 1) {
          ## if w is a parent of y
          id = which(object@pars[[y]] == w)
          if(length(id) == 1) {
            ## label w->y as compelled
            ##cat(y, object@pars[[y]][id[1]], " compelled\n")
            label[[y]][id[1]] <- 1
          }
          else {
            ## label all pars of y as compelled
            for(j in 1:length(label[[y]])){
              ##cat(y, object@pars[[y]][j], " compelled\n")
              label[[y]][j] <- 1
              ## IMPORTANT: break the cycle
              bnextcycle <- TRUE
            }
          }
        }
    }
    if(bnextcycle)
      next
    breverse <- 1
    if(length(label[[y]]) > 0) {
      for(j in 1:length(label[[y]])) {
        if(j == minx)
          next
        z <- object@pars[[y]][j]
        if(sum(object@pars[[x]] == z) == 0) {
          ## label all pars of y as compelled
          breverse <- 0
          ## label x->y as compelled
          ##cat(y, object@pars[[y]][minx], " compelled\n")
          label[[y]][minx] <- 1
          ## and all unknowns incident to y also
          for(j in 1:length(label[[y]]))
            if(label[[y]][j] == -1) {
              ##cat(y, object@pars[[y]][j], " compelled\n")
              label[[y]][j] <- 1
            }
        }
      }
    }
    if(breverse) {
      ## label x->y and all unknown edges incident to y as reversible
      ##cat(y, object@pars[[y]][minx], " reversible\n")
      label[[y]][minx] <- 0
      for(j in 1:length(label[[y]])) {
        if(label[[y]][j] == -1) {
          ##cat(y, object@pars[[y]][j], " reversible\n")
          label[[y]][j] <- 0
        }
      }
    }
  }
  
  return(label)
}

setMethod("dag2cpdag", "catNetwork",   
	function(object) {
          labels <- labelEdges(object)
          cpdag <- object
          for(k in 1:object@numnodes) {
            if(length(labels[[k]]) == 0)
              next
            for(kk in 1:length(labels[[k]])) 
              if(labels[[k]][kk] == 0) {
                x <- object@pars[[k]][kk]
                y <- k
                ## add y->x edge
                cat("add ", object@nodes[y], "->", object@nodes[x], "\n")
                cpdag@pars[[x]] <- c(cpdag@pars[[x]], y)
              }
          }
          return(cpdag)
	})  

setMethod("cnMarParents", "catNetwork",   
	function(object, flags=NULL) {
          n <- object@numnodes
	  if(!is.null(flags) && length(flags) != n)
            flags <- NULL
          an <- object@pars
          head <- an
          brep <- TRUE
          while(brep) {
            brep <- FALSE
            newhead <- vector("list",n)
            for(nh in 1:n) {
              vh <- head[[nh]]
              if(length(vh)>0) {
                pars <- NULL
                for(nvh in vh) {
                  for(pp in object@pars[[nvh]]) {
                    if(sum(an[[nh]]==pp)>0)
                      next
                    an[[nh]] <- c(an[[nh]], pp)
                    pars <- c(pars, pp)
                  }
                }
                if(!is.null(pars)) {
                  newhead[[nh]] <- pars
                  brep <- TRUE
                }
              }
            }
            head <- newhead
          }
          if(!is.null(flags)) {
            head <- lapply(1:n, function(i) if(flags[i]) return(an[[i]]) else return(NULL))
            an <- head
          }
          
          b <- vector("list",n)
          for(k in 1:n) {
            for(i in 1:n) {
              if(k==i || sum(an[[k]]==i)>0 || sum(object@pars[[i]]==k)>0) 
                next
              for(pp in an[[i]]) {
                if(sum(b[[k]]==pp)>0)
                  next
                b[[k]] <- c(b[[k]], pp)
              }
            }
          }
          
          r <- vector("list",n)
          for(k in 1:n) {
            for(i in an[[k]]) {
              if(sum(b[[k]]==i)>0)
                next
              ins <- FALSE
              for(pp in an[[i]])
                if(sum(b[[k]]==pp)>0) {
                  ins <- TRUE
                  break
                }
              if(ins)
                next
              r[[k]] <- c(r[[k]],i)
            }
          }
          
          names(r) <- object@nodes
          r <- lapply(r, function(rr) object@nodes[rr])
          return(r)
	})  
  
