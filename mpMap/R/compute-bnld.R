compute.bnld <-
function(mpcross1, mpcross2, rmat)
{
  output <- list()
  nmrk1 <- ncol(mpcross1$founders)
  nmrk2 <- ncol(mpcross2$founders)

  # pairs of loci
  pairs <- expand.grid(1:nmrk1, (nmrk1+1):(nmrk1+nmrk2))

  n.pairs <- nrow(pairs)

  # order by rows
  pairs <- pairs[order(pairs[,1], pairs[,2]), ]

  # combine two sets of data
  finals <- cbind(mpcross1$id, mpcross1$finals, mpcross2$finals)
  founders <- cbind(mpcross1$founders, mpcross2$founders)

  pedigree <- mpcross1$pedigree 

  rpairs <- CR_calcLD(finals, founders, pedigree, pairs, rmat)

  output$ld <- list()
  output$ld$W <- matrix(data=rpairs[,1], nrow=nmrk1, ncol=nmrk2, byrow=TRUE)
  output$ld$LewontinD <- matrix(data=rpairs[,2], nrow=nmrk1, ncol=nmrk2, byrow=TRUE)
  output$ld$delta2 <- matrix(data=rpairs[,3], nrow=nmrk1, ncol=nmrk2, byrow=TRUE)
  output$ld$r2 <- matrix(data=rpairs[,4], nrow=nmrk1, ncol=nmrk2, byrow=TRUE)

  return(output)
}

