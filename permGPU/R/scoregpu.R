scoregpu = function(y, event, markers, test, B=0, stand=TRUE, pval=TRUE, index=FALSE, scale=FALSE){
  error=FALSE
  if(!is.element(test,c("cox","npcox"))){
    print("This test has not been implemented yet")
    error = TRUE
  }

  gid = rownames(markers)
  n = ncol(markers)
  K = nrow(markers)

  if (scale)
    markers = scale(t(markers), scale=FALSE)
  else
    markers = t(markers)

  if(!error){
    if (index){
      ind = 1
      output = rep(0, (K+n)*(B+1))
    }
    else{
      output=rep(0,K*(B+1))
      ind = 0
    }
    out <- .C("scoregpu", as.single(markers), as.single(y), as.single(event), as.integer(n), as.integer(K), as.character(test), as.integer(B), as.integer(ind), as.integer(as.numeric(stand)), as.integer(as.numeric(pval)), as.single(output), NAOK=TRUE, PACKAGE="permGPU")
    out = out[[11]]
    stat = out[1:(K*(B+1))]
    stat = matrix(stat, K, B+1) 
    colnames(stat) = paste("P", seq(0,B), sep="")
    rownames(stat) = gid
    
    if (index){
      indexmat = out[(K*(B+1)+1):((K+n)*(B+1))]
      indexmat = matrix(indexmat, n, B+1)
      rownames(indexmat) = paste("N", seq(1,n), sep="")
      colnames(indexmat) = paste("P", seq(0,B), sep="")
    }

    if(index)
      out = list(stat=stat, index=indexmat)
    else
      out = stat
  }
  else
    {
      out = NULL
    }

  return(out)
  }
