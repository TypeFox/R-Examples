#########################################################################
# Clustering Routines

setMethod("cnCluster", "catNetwork", 
function(object) {

  norder <- cnOrder(object)
  cl <- vector("list", object@numnodes)
  nodecl <- rep(0, object@numnodes)
  icl <- 1
  cl[[1]] <- norder[1]
  nodecl[norder[1]] <- 1
  for(j in norder[2:object@numnodes]) {
    for(k in object@pars[[j]]) {
      cl[[nodecl[k]]] <- c(cl[[nodecl[k]]], j)
      nodecl[j] <- nodecl[k]
    }
    if(nodecl[j] < 1) {
      icl <- icl + 1
      cl[[icl]] <- j
      nodecl[j] <- icl
    }
  }
  for(icl in length(cl):1)
    if(length(cl[[icl]]) < 2)
       cl[icl] <- NULL

  if(length(cl) > 1) {
    b <- TRUE
    while(b) {
      b <- FALSE
      for(i in length(cl):2) {
        ##cat("i = ", i,"\n")
        if(length(cl[[i]]) < 1)
          next
        for(j in (i-1):1) {
          if(length(cl[[j]]) < 1)
            next
          for(k in cl[[j]]) {
            if(length(which(cl[[i]] == k)) > 0) {
              b <- TRUE
              cl[[j]] <- c(cl[[j]], cl[[i]])
              cl[[i]] <- NULL
              ##if(length(cl[[j]])>0)
              ##  cat(j,": ", cl[[j]], "\n")
              break
            }
          }
          if(b)
            break
        }
      }
    }
  }

  for(k in 1:length(cl)) {
    cc <- cl[[k]]
    for(i in length(cc):2)
      if(length(which(cc[1:(i-1)] == cc[i])))
        cc <- cc[-i]
    cl[[k]] <- cc
  }
  
  for(i in 1:length(cl)) {
    names(cl[[i]]) <- object@nodes[cl[[i]]]
  }
  
  return(cl)
})

cnClusterMI <- function(data, pert=NULL, threshold=0) {
  if(!is.numeric(threshold))
    stop("The threshold parameter should be numeric")
  mentropy <- cnEntropy(data, pert)
  numnodes <- nrow(mentropy)
  for(i in 1:numnodes)
    mentropy[i,] <- mentropy[i,i]-mentropy[i,]
  threshold <- threshold[threshold>=0]
  ##cat(threshold,"\n")
  cls.list <- vector("list", length(threshold))
  it <- 0
  for(t in threshold) {
    mm <- mentropy
    cls <- vector("list", numnodes)
    icls <- 0
    ncls <- rep(0, numnodes)
    while(1) {
      ff <- max(mm)
      if(ff < t)
        break
      k <- which(mm == ff)[1]
      j <- floor((k-1)/numnodes)
      i <- k - numnodes*j
      j <- j + 1
      mm[i,j] <- -1
      ##cat("(i,j) = ", i,", ", j, " (", ff, ") \n")
      if(ncls[i] > 0 && ncls[j] > 0 && ncls[i] != ncls[j]) {
        ## unite
        ##cat(ncls[i], ", ", ncls[j], "\n")
        ##cat("predi : ")
        ##for(s in 1:icls)
        ##  cat(cls[[s]], "\n")
        kk <- ncls[[j]]
        cls[[ncls[i]]] <- c(cls[[ncls[i]]], cls[[ncls[j]]])
        ncls[cls[[ncls[j]]]] <- ncls[i]
        cls[[kk]] <- NULL        
        icls <- icls - 1
        ## all indices after kk decreases by 1
        for(ss in kk : icls)
          ncls[cls[[ss]]] <- ss
        ##cat("sled : ")
        ##for(s in 1:icls) {
        ##  cat(cls[[s]], "\n")
        ##  cat(ncls[cls[[s]]], "\n")
        ##}
        next
      }
      if(ncls[i] > 0 && ncls[j] == 0) {
        ## add j to i
        ncls[j] <- ncls[i]
        cls[[ncls[i]]] <- c(cls[[ncls[i]]], j)
        next
      }
      if(ncls[j] > 0 && ncls[i] == 0) {
        ## add i to j
        ncls[i] <- ncls[j]
        cls[[ncls[j]]] <- c(cls[[ncls[j]]], i)
        next
      }
      if(ncls[i] == 0 && ncls[j] == 0) {
        ## add i and j in a new cluster
        icls <- icls + 1
        cls[[icls]] <- c(i,j)
        ncls[j] <- icls
        ncls[i] <- icls
        next
      }
    }
    if(icls < 1)
      next
    it <- it + 1
    cls.list[[it]] <- lapply(1:icls, function(i) cls[[i]])
  }
  if(it < 1)
    return(NULL)
  return(lapply(1:it, function(i) cls.list[[i]]))
}

setMethod("cnClusterSep", "catNetwork", 
function(object, data, pert=NULL) {

  numnodes <- object@numnodes
  mm <- cnEntropy(data, pert)
  ## form the joint entropy matrix
  for(i in 1:numnodes)
    mm[i,] <- mm[i,i] - mm[i,]

  sep.list <- vector("list", 5)
  names(sep.list) <- c("clusters", "tmats", "in.entropy", "out.entropy", "connected")

  cls1 <- cnCluster(object)
  sep.list[[1]] <- cls1

  in.entropy <- rep(0,length(cls1))
  out.entropy <- rep(0,length(cls1))
  tmats <- vector("list", length(cls1))
  connected <- rep(FALSE, length(cls1))
  
  for(nc in 1:length(cls1)) {
    in.entropy[nc] <- Inf
    out.entropy[nc] <- -Inf
    for(i in cls1[[nc]]) {
      max.in <- -Inf
      for(j in cls1[[nc]]) {
        if(j == i)
          next
        if(max.in < mm[i,j]) {
          max.in <- mm[i,j]
        }
      }
      if(in.entropy[nc] > max.in)
        in.entropy[nc] <- max.in
      max.out <- -Inf
      for(j in 1:object@numnodes) {
        if(length(which(cls1[[nc]] == j))>0)
          next
        if(max.out < mm[i,j]) {
          max.out <- mm[i,j]
        }
      }
      if(out.entropy[nc] < max.out)
        out.entropy[nc] <- max.out
    }
    clen <- length(cls1[[nc]])
    mat <- matrix(rep(1, clen^2), nrow=clen)
    for(i in 1:(clen-1)) {
      for(j in (i+1):clen) {
        if(j == i)
          next
        if(mm[cls1[[nc]][i],cls1[[nc]][j]] < in.entropy[nc]-1e-8) {
          mat[i,j] <- 0
          mat[j,i] <- 0
        }
      }
    }
    rownames(mat) <- cls1[[nc]]
    colnames(mat) <- cls1[[nc]]
    tmats[[nc]] <- mat

    m <- mat
    brep <- TRUE
    while(brep) {
      brep <- FALSE
      for(i1 in 1:(clen-1))
        for(j1 in 1:(clen-1))
          for(i2 in (i1+1):clen)
            for(j2 in (j1+1):clen) {
              f <- m[i1,j1] + m[i1,j2] + m[i2,j1] + m[i2,j2]
              if(f == 3) {
                m[i1,j1] <- 1
                m[i1,j2] <- 1
                m[i2,j1] <- 1
                m[i2,j2] <- 1
                brep <- TRUE
              }
            }
    }
    connected[nc] <- prod(m)>0
  }
  sep.list[[2]] <- tmats
  sep.list[[3]] <- in.entropy
  sep.list[[4]] <- out.entropy
  sep.list[[5]] <- connected
  return(sep.list)
})

