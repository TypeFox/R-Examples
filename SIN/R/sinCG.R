sinCG <- function(blocks, S, n, type="AMP", holm=TRUE){
  if(!is.blocks(blocks, dim(S)[1])){
    return("blocks is not a valid block structure over the variables!")
  }
  ## blocks is now a list of vectors 
  sinAMP <- function(blocks, S, n){
    complete.order <- unlist(blocks)
    old.order <- dimnames(S)[[1]]
    S <- S[complete.order,complete.order]
    p <- dim(S)[1]
    q <- length(blocks)
    pvals <- diag(rep(1,p))
    labels <- dimnames(S)[[1]]
    dimnames(pvals) <- list(labels, labels)
    ## Do block 1 explicitly (since no parents)
    b <- cumsum(unlist(lapply(blocks, length)))
    corr.part <- -cov2cor(solve(S[1:b[1],1:b[1]]))
    pvals[1:b[1],1:b[1]] <- simpvalueMx(corr.part,n-b[1]-1,p)
    ## Loop thru blocks 2 to q
    if(q==1){
      return(zapsmall(pvals[old.order,old.order]))
    }
    else{
      for(i in 2:q){
        ## within block i
        corr.part <- -cov2cor(solve(S[1:b[i],1:b[i]]))
        pvals[(b[i-1]+1):b[i],(b[i-1]+1):b[i]] <-
          simpvalueMx(corr.part[(b[i-1]+1):b[i],(b[i-1]+1):b[i]],n-b[i]-1,p)
        ## block i to strict past
        for(w in (b[i-1]+1):b[i]){
          corr.part <- -cov2cor(solve(S[c(1:b[i-1], w),c(1:b[i-1], w)]))
          pvals[1:b[i-1],w] <-
            pvals[w,1:b[i-1]] <-
              simpvalueVec(corr.part[1:b[i-1],b[i-1]+1],
                           n-(b[i-1]+1)-1,p)
        }
      }
      return(zapsmall(pvals[old.order,old.order]))
    }
  }
  ## blocks is now a list of vectors 
  sinLWF <- function(blocks, S, n){
    complete.order <- unlist(blocks)
    old.order <- dimnames(S)[[1]]
    S <- S[complete.order,complete.order]
    p <- dim(S)[1]
    q <- length(blocks)
    pvals <- diag(rep(1,p))
    labels <- dimnames(S)[[1]]
    dimnames(pvals) <- list(labels, labels)
    ## Do block 1 explicitly (since no parents)
    b <- cumsum(unlist(lapply(blocks, length)))
    corr.part <- -cov2cor(solve(S[1:b[1],1:b[1]]))
    pvals[1:b[1],1:b[1]] <- simpvalueMx(corr.part,n-b[1]-1,p)
    ## Loop thru blocks 2 to q
    if(q==1){
      return(zapsmall(pvals[old.order,old.order]))
    }
    else{
      for(i in 2:q){
        ## block i vs past i
        corr.part <- -cov2cor(solve(S[1:b[i],1:b[i]]))
        pvals[(b[i-1]+1):b[i],1:b[i]] <-
          simpvalueMx(corr.part,n-b[i]-1,p)[(b[i-1]+1):b[i],1:b[i]]
        pvals[1:b[i],(b[i-1]+1):b[i]] <- t(pvals[(b[i-1]+1):b[i],1:b[i]])
      }
      return(zapsmall(pvals[old.order,old.order]))
    }
  }
  ## Call one of the above functions now
  if(type=="AMP"){
    pvals <- sinAMP(blocks,S,n)
  }
  else{
    if(type=="LWF"){
    pvals <- sinLWF(blocks,S,n)
    }
    else{
      return("type must be AMP or LWF (as string)!")
    }
  }
  if(holm==TRUE){
    return(holm(pvals))
  } else{
    return(pvals)
  }
}

 
plotCGpvalues <- function(blocks, pvals, legend=TRUE, legendpos=NULL){
  if(!is.blocks(blocks, dim(pvals)[1])){
    return("blocks is not a valid block structure over the variables!")
  }
  ## two functions to order p-values and create labels
  vecCGpvalues <- function(blocks, pvals){
    p <- dim(pvals)[1]
    q <- length(blocks)
    b <- cumsum(unlist(lapply(blocks, length)))
    pvec <- c()
    oo <- unlist(blocks)
    pvals <- pvals[oo,oo]
    ## Block 1 first
    if(b[1]>=2){
      for(j in 1:(b[1]-1)){
        pvec <- c(pvec, pvals[j,(j+1):b[1]])
      }
    }
    ## Blocks 2 thru q
    for(i in 2:q){
      ## directed edges 
      for(j in (b[i-1]+1):b[i]){
        for(k in 1:b[i-1]){
          pvec <- c(pvec, pvals[k,j])
        }
      }
      ## undirected edges in Block i
      if(b[i]-b[i-1]>=2){
        for(j in (b[i-1]+1):(b[i]-1)){
          pvec <- c(pvec, pvals[j,(j+1):b[i]])
        }
      }
    }
    return(pvec)
  }
  
  createCGlabels <- function(blocks, pvals){
    p <- dim(pvals)[1]
    q <- length(blocks)
    b <- cumsum(unlist(lapply(blocks, length)))
    labels <- c()
    ## Block 1 first
    if(b[1]>=2){
      for(j in 1:(b[1]-1)){
        for(k in (j+1):b[1]){
          labels <-
            c(labels,
              paste(letters[1],j,"-",letters[1],k, sep=""))
        }
      }
    }
    ## Blocks 2 thru q
    for(i in 2:q){
      ## directed edges 
      for(j in (b[i-1]+1):b[i]){
        for(l in 1:b[1]){
          labels <-
            c(labels,
              paste(letters[1],l,"->",letters[i],j-b[i-1], sep=""))
        }
        if(i >=3){
          for(k in 2:(i-1)){
            for(l in (b[k-1]+1):b[k]){
              labels <-
                c(labels,
                  paste(letters[k],l-b[k-1],"->",letters[i],j-b[i-1], sep=""))
            }
          }
        }
      }
      ## undirected edges in Block i
      if(b[i]-b[i-1]>=2){
        for(j in (b[i-1]+1):(b[i]-1)){
          for(k in (j+1):b[i]){
            labels <-
              c(labels,
                paste(letters[i],j-b[i-1],"-",letters[i],k-b[i-1], sep=""))
          }
        }
      }
    }
    return(labels)
  }
  ## actual plotting 
  par(mar=c(7,5,2,2)+0.1)
  CGlab <- createCGlabels(blocks, pvals)
  CGpvals <- vecCGpvalues(blocks, pvals)
  temp <- length(CGlab)
  plot(as.factor(1:temp), CGpvals, type="n",
       ylab="P-value", xlab="", axes=FALSE, ylim=c(0,1), cex.lab=1.2, las=2)
  title(xlab = "Edge", line = 5, cex.lab=1.2)
  axis(1, at=1:temp, labels=CGlab[1:temp], las=2)
  axis(2, at=seq(0,1,by=0.1), las=1)
  temp2 <- sapply(seq(0,1,by=0.1), abline, 0, lty="dotted", col="grey")
  plot(as.factor(1:temp), CGpvals, axes=FALSE, add=TRUE) 
  box()
  plotlabels <- c()
  for(i in 1:length(blocks)){
    for(j in 1:length(blocks[[i]])){
      plotlabels <- c(plotlabels, paste(letters[i],j,sep="" ))
    }
  }
  oo <- unlist(blocks)
  for(i in 1:length(plotlabels)){
    plotlabels[i] <- paste(plotlabels[i],dimnames(pvals[oo,oo])[[1]][i], sep="  ")
  }
  if(legend==TRUE){
    if(is.null(legendpos)){
      legend(temp+0.6,1, x.intersp=-0.3, 
             plotlabels, bg="white", xjust=1, yjust=1)
    }
    else{
      x <- legendpos[1]
      y <- legendpos[2]
      legend(x,y,  x.intersp=-0.3,
             plotlabels, bg="white", xjust=1, yjust=1)
    }
  }
}
