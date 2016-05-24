#########################################################################
# Categorical Network Class Methods
# Prediction

nodePredict <- function(idroot, ppars, idx, problist, psample) {
  if(is.null(ppars) || length(idx) < 1) {
    idcat <- which(problist == max(problist))
    ##cat(problist, ", ", idcat, "\n")
    if(length(idcat) < 1)
	  stop("Probability list is broken")
    if(length(idcat) > 1) {
      ##warning("remi")
      idcat <- idcat[1]
    }
    return(idcat)
  }
  idnode <- ppars[idx[1]]
  if(idnode > length(psample))
    stop("Wrong sample")
  cat <- psample[idnode]
  ##cat(idnode, ": ", cat, "\n")
  if(is.na(cat))
    stop("Insufficient sample")
  return(nodePredict(idroot, ppars, idx[-1], problist[[cat]], psample))
}

samplePredict <- function(object, data) {
  if(dim(data)[1] != object@numnodes)
    stop("Incompatible sample dimensions.\n")
  numSamples <- dim(data)[2]
  if(numSamples < 1)
    stop("No samples\n")
  
  nodeOrder <- cnOrder(object)
  
  for(j in (1:numSamples)) {
    ps <- as.integer(data[,j])
    for(i in 1:object@numnodes) {
      nnode <- nodeOrder[i]
      if(!is.na(ps[nnode]))
        next
      ps[nnode] <- nodePredict(nnode, object@pars[[nnode]],
                               seq(1,length(object@pars[[nnode]])), object@probs[[nnode]], ps)
    }
    
    data[,j] <- ps
  }
  return(data)
}


setMethod("cnPredict", c("catNetwork"), 
          function(object, data) {

            if(!is.matrix(data) && !is.data.frame(data))
              stop("data should be a matrix or data frame of node cats")
            asframe <- FALSE
	    if(is.data.frame(data)) {
	      data <- as.matrix(t(data))
              asframe <- TRUE
            }

            if(length(dim(data)) == 2 && dim(data)[1] != object@numnodes)
	      stop("the number of nodes in 'object' and 'data' should be equal")

	    rownames <- rownames(data)
            if(length(rownames) != object@numnodes)
              stop("The data rows should be named after the nodes of the object.")

            if(prod(tolower(rownames) == tolower(object@nodes)) == 0) {
              norder <- order(rownames)
              data <- data[norder,]
              rownames <- rownames(data)
              norder <- order(object@nodes)
              object <- cnReorderNodes(object, norder)
            }

            if(prod(tolower(rownames) == tolower(object@nodes)) == 0)
              stop("The row names should correspond to the object nodes.")

            for(i in 1:dim(data)[1]) {

              noderow <- data[i, !is.na(data[i,])]
              if(length(noderow) > 0 && is.integer(noderow) &&
                 min(noderow) >= 1 &&
                 max(noderow) <= length(object@cats[[i]])) {
                for(j in 1:length(data[i,])) {
                  if(is.na(data[i,j]))
                    next
                  data[i,j] <- as.integer(data[i,j])
                }
              }
              else {                
                for(j in 1:length(data[i,])) {
                  if(is.na(data[i,j]))
                    next                  
                  c <- which(object@cats[[i]] == data[i,j])
                  if(length(c) < 1 || c < 1 || c > length(object@cats[[i]])) {
                    warning("Category ", data[i,j], " for node ", i, " sample ", j, " is illigal - reset to NA")
                    c <- NA
                  }
                  data[i,j] <- c
                }
              }
            }

            newdata <- samplePredict(object, data)

            for(j in (1:dim(data)[2])) {
              ps <- as.integer(newdata[,j])
              for(i in 1:dim(data)[1])
                ps[i] <- object@cats[[i]][as.integer(ps[i])]
              newdata[,j] <- ps
            }

            if(asframe)
              return(as.data.frame(t(newdata)))
            return(newdata)
          })

