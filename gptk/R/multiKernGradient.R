multiKernGradient <-
function (kern, x, x2, covGrad) {
  if ( is.list(x) ) {
    if ( (nargs()>3) & !is.list(x2) )
      stop("Time course information is not matched in List format.")

    arg1 <- list()
    arg2 <- list()
    dim1 <- array()
    dim2 <- array()
    for ( i in seq(length=kern$numBlocks) ) {
      dim1[i] <- dim(as.array(x[[i]]))[1]
      arg1[[i]] <- x[[i]]
      if ( nargs()>3 ) {
        dim2[i] <- dim(as.array(x2[[i]]))[1]
        arg2[[i]] <- x2[[i]]
      } else {
        dim2[i] <- dim1[i]
        covGrad <- x2
        arg2[[i]] <- arg1[[i]]
      }
    }

    g <- array(0, dim(kern$paramGroups)[1])
    startVal <- 1
    endVal <- 0

    for ( i in seq(length=kern$numBlocks) ) {
      endVal <- endVal + kern$comp[[i]]$nParams
      startOne <- sum(dim1[seq(length.out=(i-1))])+1
      endOne <- sum(dim1[seq(length.out=i)])
      startThree <- sum(dim2[seq(length.out=(i-1))])+1
      endThree <- sum(dim2[seq(length.out=i)])

      if ( nargs()>3 ) {
        g[startVal:endVal] <- .multiKernGradientBlock(kern, arg1[[i]], arg2[[i]], covGrad[startOne:endOne, startThree:endThree], i, i)
      } else {
        g[startVal:endVal] <- .multiKernGradientBlock(kern, arg1[[i]], covGrad[startOne:endOne, startThree:endThree], i, i)
      }

      startVal2 <- 1
      endVal2 <- 0
      
      for ( j in seq(length.out=(i-1)) ) {
        endVal2 <- endVal2 + kern$comp[[j]]$nParams
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- sum(dim2[seq(length.out=(j-1))])+1
          endTwo <- sum(dim2[seq(length.out=j)])

          gList <- .multiKernGradientBlock(kern, arg1[[i]], arg2[[j]], covGrad[startOne:endOne, startTwo:endTwo], i, j)

          g[startVal:endVal] <- g[startVal:endVal] + 2*gList$g1
          g[startVal2:endVal2] <- g[startVal2:endVal2] + 2*gList$g2
        }
        startVal2 <- endVal2 + 1
      }
      startVal <- endVal + 1
    }

  } else {   # non-list x
    dim1 <- dim(as.array(x))[1]
    arg1 <- x
    if ( nargs() > 3 ) {
      dim2 <- dim(as.array(x2))[1]
      arg2 <- x2
    } else {
      dim2 <- dim1
      covGrad <- x2
      arg2 <- arg1
    }

    g <- array(0, dim(kern$paramGroups)[1])
    startVal <- 1
    endVal <- 0
    for ( i in seq(length=kern$numBlocks) ) {
      endVal <- endVal + kern$comp[[i]]$nParams
      startOne <- (i-1)*dim1 + 1
      endOne <- i*dim1
      if ( nargs() > 3 ) {
        g[startVal:endVal] <- .multiKernGradientBlock(kern, arg1, arg2, covGrad[startOne:endOne, ((i-1)*dim2+1):(i*dim2)], i, i)
      } else {
        g[startVal:endVal] <- .multiKernGradientBlock(kern, arg1, covGrad[startOne:endOne, ((i-1)*dim2+1):(i*dim2)], i, i)
      }
      
      startVal2 <- 1
      endVal2 <- 0
      
      for ( j in seq(length=(i-1)) ) {
        endVal2 <- endVal2 + kern$comp[[j]]$nParams
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- (j-1)*dim2 + 1
          endTwo <- j*dim2

          gList <- .multiKernGradientBlock(kern, arg1, arg2, covGrad[startOne:endOne, startTwo:endTwo], i, j)

          g1 <- gList$g1
          g2 <- gList$g2

          if ( nargs() > 3 ) {
            startThree <- (j-1)*dim1 + 1
            endThree <- j*dim1
            gList <- .multiKernGradientBlock(kern, arg2, arg1, t(covGrad[startThree:endThree, startTwo:endTwo]), j, i)

            g3 <- gList$g1
            g4 <- gList$g2
            g[startVal:endVal] <- g[startVal:endVal] + g1 + g4
            g[startVal2:endVal2] <- g[startVal2:endVal2] + g2 + g3
          } else {
            g[startVal:endVal] <- g[startVal:endVal] + 2*g1
            g[startVal2:endVal2] <- g[startVal2:endVal2] + 2*g2           
          }
        }
        startVal2 <- endVal2 + 1
      }
      startVal <- endVal + 1
    }
  }

  g <- (g %*% kern$paramGroups)[1,]
  return (g)
}
