multiKernCompute <-
function (kern, x, x2=x) {
  if ( is.list(x) ) {
    if ( length(x) != kern$numBlocks )
      stop ("Time information is not matched among blocks!")

    dim1 <- array(0, dim=length(x))
    dim2 <- array(0, dim=length(x))
    
    for ( i in seq(length=kern$numBlocks) ) {
      dim1[i] <- dim(as.array(x[[i]]))[1]
      if ( nargs()>2 ) {
        if ( length(x) != length(x2) )
          stop ("Time information is not matched within the block!")
        dim2[i] <- dim(as.array(x2[[i]]))[1]
      } else {
        dim2[i] <- dim1[i]
      }
    }

    K <- matrix(0, sum(dim1), sum(dim2))

    for ( i in seq(length=kern$numBlocks) ) {
      startOne <- sum(dim1[seq(length.out=(i-1))])+1
      endOne <- sum(dim1[seq(length.out=i)])
      startThree <- sum(dim2[seq(length.out=(i-1))])+1
      endThree <- sum(dim2[seq(length.out=i)])

      if ( nargs()<3 ) {
        K[startOne:endOne, startThree:endThree] <- .multiKernComputeBlock(kern, i, i, x[[i]])
      } else {
        K[startOne:endOne, startThree:endThree] <- .multiKernComputeBlock(kern, i, i, x[[i]], x2[[i]])
      }

      for ( j in seq(length.out=(i-1)) )
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- sum(dim2[seq(length.out=(j-1))])+1
          endTwo <- sum(dim2[seq(length.out=j)])

          if ( nargs()<3 ) {
            K[startOne:endOne, startTwo:endTwo] <- .multiKernComputeBlock(kern, i, j, x[[i]], x[[j]])
            K[startTwo:endTwo, startOne:endOne] <- t(K[startOne:endOne, startTwo:endTwo])
          } else {
            K[startOne:endOne, startTwo:endTwo] <- .multiKernComputeBlock(kern, i, j, x[[i]], x2[[j]])
            startFour <- sum(dim1[seq(length.out=(j-1))])+1
            endFour <- sum(dim1[seq(length.out=j)])
            K[startFour:endFour, startThree:endThree] <- t(.multiKernComputeBlock(kern, j, i, x2[[i]], x[[j]]))
          }
        }
    }
  } else {
                                        # non-cell part
    dim1 = dim(as.array(x))[1]
    
    if ( nargs() > 2 ) {
      dim2 = dim(as.array(x2))[1]
    } else {
      dim2 = dim1;
    }
    
    K <- matrix(0, kern$numBlocks*dim1, kern$numBlocks*dim2)
    
    for ( i in seq(length=kern$numBlocks) ) {
      startOne <- (i-1)*dim1 + 1
      endOne <- i*dim1
      startThree <- (i-1)*dim2 + 1
      endThree <- i*dim2
      if ( nargs() < 3 ) {
        K[startOne:endOne, startThree:endThree] <- .multiKernComputeBlock(kern, i, i, x)
      } else {
        K[startOne:endOne, startThree:endThree] <- .multiKernComputeBlock(kern, i, i, x, x2)
      }

      for ( j in seq(length=(i-1)) ) {
        if ( !is.na(kern$block[[i]]$cross[j]) ) {
          startTwo <- (j-1)*dim2 + 1
          endTwo <- j*dim2
          if ( nargs() < 3 ) {
            K[startOne:endOne, startTwo:endTwo] <- .multiKernComputeBlock(kern, i, j, x)
          } else {
            K[startOne:endOne, startTwo:endTwo] <- .multiKernComputeBlock(kern, i, j, x, x2)
          }
          if ( nargs()< 3 ) {
            K[startTwo:endTwo, startOne:endOne] <- t(K[startOne:endOne, startTwo:endTwo])
          } else {
            startFour <- (j-1)*dim1 + 1
            endFour <- j*dim1
            K[startFour:endFour, startThree:endThree] <- t(.multiKernComputeBlock(kern, j, i, x2, x))
          }
        }
      }
    }
  }

  return (K)
}
