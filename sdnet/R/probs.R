#########################################################################
# Categorical Network Class Methods
# Probability Calculations

listProbSet <- function(idroot, ppars, pcatlist, idx, problist, strin) { 
  if(is.null(ppars) || length(idx) < 1) { 
    if(length(pcatlist[[idroot]]) != length(problist)) {
      ##cat(idroot, ",   ", length(pcatlist[[idroot]]), ", ", length(problist), "\n")
      warning("Wrong probability slot") 
      return("") 
    } 
    strout <- sapply(seq(1, length(pcatlist[[idroot]])), function(j, strin, pcatlist, problist) { 
      paste("[", strin, "]", pcatlist[j], "  ", problist[j], "\n", sep="") 
    }, strin, pcatlist[[idroot]], problist) 
    return(paste(strout, sep="", collapse="")) 
  } 
  idnode <- ppars[idx[1]] 
  idx <- idx[-1] 
  strout <- sapply(seq(1,length(pcatlist[[idnode]])), 
    function(cat, strin) { 
      listProbSet(idroot, ppars, pcatlist, idx, problist[[cat]], paste(strin, pcatlist[[idnode]][cat], sep=" ")) 
    }, strin) 
  return(paste(strout, sep="", collapse="")) 
}

listProbTable <- function(idroot, ppars, pcatlist, idx, problist, pars) { 
  if(is.null(ppars) || length(idx) < 1) { 
    if(length(pcatlist[[idroot]]) != length(problist)) {
      warning("Wrong probability slot") 
      return("") 
    }
    return(c(pars, problist))
  } 
  idnode <- ppars[idx[1]] 
  idx <- idx[-1] 
  rout <- NULL
  for(cat in 1:length(pcatlist[[idnode]])) { 
    rout <- rbind(rout, listProbTable(idroot, ppars, pcatlist, idx, problist[[cat]], c(pars, as.character(pcatlist[[idnode]][cat]))))
  }
  rout <- as.table(rout)
  return(rout)
} 
 
setMethod("cnPlotProb", c("catNetwork"),  
          function(object, which=NULL) { 
            if(is.null(which)) 
              which <- seq(1, object@numnodes)
            if(is.character(which))
              which <- sapply(which, function(node) which(object@nodes == node))
            str <- sapply(which, function(n) { 
              paste(paste("Node[", object@nodes[n], "], Parents: ", sep=""), 
                    paste(object@nodes[object@pars[[n]]], sep=",", collapse=" "), "\n",  
                    listProbSet(n, object@pars[[n]], object@cats, 
                                seq(1,length(object@pars[[n]])), 
                                object@probs[[n]], ""), 
                    collapse="", sep="") 
            }) 
            cat(str)
          })

setMethod("cnProb", c("catNetwork"),  
          function(object, which=NULL) { 
            if(is.null(which)) 
              which <- seq(1, object@numnodes)
            if(is.character(which))
              which <- sapply(which, function(node) which(object@nodes == node))
            ltab <- lapply(which,  function(n) { 
              ntab <- listProbTable(n, object@pars[[n]], object@cats, 
                                    seq(1,length(object@pars[[n]])), 
                                    object@probs[[n]], NULL)
              if(length(object@pars[[n]]) > 0) {
                colnames(ntab) <- c(object@nodes[object@pars[[n]]], object@cats[[n]])
                for(j in (length(object@pars[[n]])+1):ncol(ntab))
                  ntab[,j] <- as.numeric(ntab[,j])
              }
              else {                      
                ntab <- as.numeric(ntab)
                names(ntab) <- object@cats[[n]]
              }
              return(ntab)
            })
            names(ltab) <- object@nodes[which]
            return(ltab)
          }) 

checkProbSet <- function(idroot, ppars, pcatlist, idx, problist) {
  if(length(ppars) == 0 || length(ppars) < 1 || length(idx) < 1) {
    if(is.null(problist) || length(pcatlist[[idroot]]) != length(problist)) {
      return(FALSE)
    }
    if(!is.nan(sum(problist)) || sum(problist) == -Inf)
      return(TRUE)
    if(sum(problist<0) > 0) {
      warning("Probability slot with negative values", "\n")
      return(FALSE)
    }
    if(as.integer(sum(problist)) > 1) {
      warning("Probability slot with sum ", sum(problist), "\n")
      return(FALSE)
    }
    return(TRUE)
  }
  idnode <- ppars[idx[1]]
  res <- sapply(seq(1,length(pcatlist[[idnode]])),
                function(cat)
    checkProbSet(idroot, ppars, pcatlist, idx[-1], problist[[cat]]))
  if(sum(res) < length(pcatlist[[idnode]]))
    return(FALSE)
  return(TRUE)
}

## Attach a new leaf-node with a parent set [leafpars] and conditional probability [leafproblist]
## to a tree-list [ptree] that includes all nodes less in order than the new node
## THE PARENTS IN [leafpars] HAS THE SAME ORDER AS THAT IN [ptree],
## the latter is assured if the catNetwork is created by genRandomParents

# recursive probability assignment
setRandomProb <- function(idroot, ppars, pcatlist, idx, delta1=0.01, delta2=0.01) {
  if(is.null(ppars) || length(ppars) < 1 || length(idx) < 1) {
    if(length(pcatlist[[idroot]]) < 1) {
      return(NULL)
    }
    ii <- 0
    while(ii<10*(1+floor(1/(0.5-delta1-delta2)))^length(pcatlist[[idroot]])) {
      plist <- sapply(pcatlist[[idroot]], function(x) runif(1))
      plist <- plist/sum(plist)
      if(sum(plist<delta1 | plist>1-delta1)+sum(plist>0.5-delta2 & plist<0.5+delta2) == 0)
        break
      ii <- ii + 1
    }
    plist <- sapply(seq(1,length(plist)), function(n, plist){
      plist[n] <- floor(1000*plist[n])/1000
      }, plist)
    plist[1] <- 1 - sum(plist[-1])
    #poutlist <- c(poutlist, plist)
    return(as.vector(plist))
  }
  else {
    id <- ppars[idx[1]]
    #cat(ppars[idx], id, "\n")
    poutlist <- lapply(pcatlist[[id]],
           function(cat) setRandomProb(idroot, ppars, pcatlist, idx[-1], delta1, delta2))
  }
}

# recursive probability assignment
setDefaultProb <- function(idroot, ppars, pcatlist, idx) {
  if(is.null(ppars) || length(ppars) < 1 || length(idx) < 1) {
    if(length(pcatlist[[idroot]]) < 1) {
      return(NULL)
    }
    plist <- rep(1, length(pcatlist[[idroot]]))
    plist <- plist/sum(plist)
    plist <- sapply(seq(1,length(plist)), function(n, plist){
      plist[n] <- floor(100*plist[n])/100
      }, plist)
    plist[1] <- 1 - sum(plist[-1])
    return(as.vector(plist))
  }
  else {
    id <- ppars[idx[1]]
    ##cat(id, ": ",idx, "\n")
    poutlist <- lapply(pcatlist[[id]],
           function(cat, idroot, ppars, pcatlist, idx)
           setDefaultProb(idroot, ppars, pcatlist, idx),
           idroot, ppars, pcatlist, idx[-1])
  }
}
