combine_ld <-
function(mpcross1, mpcross2, rmat)
{
  output <- list()

  ldb <- compute.bnld(mpcross1, mpcross2, rmat)

  # now combine the four pieces
  nmrk1 <- ncol(mpcross1$founders)
  nmrk2 <- ncol(mpcross2$founders)
  ntot <- nmrk1+nmrk2

  output$W <- matrix(nrow=ntot, ncol=ntot)
  output$LewontinD <- matrix(nrow=ntot, ncol=ntot)
  output$delta2 <- matrix(nrow=ntot, ncol=ntot)
  output$r2 <- matrix(nrow=ntot, ncol=ntot)

  output$W[1:nmrk1, 1:nmrk1] <- mpcross1$ld$W
  output$LewontinD[1:nmrk1, 1:nmrk1] <- mpcross1$ld$LewontinD
  output$delta2[1:nmrk1, 1:nmrk1] <- mpcross1$ld$delta2
  output$r2[1:nmrk1, 1:nmrk1] <- mpcross1$ld$r2

  output$W[1:nmrk1, (nmrk1+1):ntot] <- (ldb$W)
  output$LewontinD[1:nmrk1, (nmrk1+1):ntot] <- (ldb$LewontinD)
  output$delta2[1:nmrk1, (nmrk1+1):ntot] <- (ldb$delta2)
  output$r2[1:nmrk1, (nmrk1+1):ntot] <- (ldb$r2)

  output$W[(nmrk1+1):ntot, 1:nmrk1] <- t(ldb$W)
  output$LewontinD[(nmrk1+1):ntot, 1:nmrk1] <- t(ldb$LewontinD)
  output$delta2[(nmrk1+1):ntot, 1:nmrk1] <- t(ldb$delta2)
  output$r2[(nmrk1+1):ntot, 1:nmrk1] <- t(ldb$r2)

  output$W[(nmrk1+1):ntot, (nmrk1+1):ntot] <- mpcross2$ld$W
  output$LewontinD[(nmrk1+1):ntot, (nmrk1+1):ntot] <- mpcross2$ld$LewontinD
  output$delta2[(nmrk1+1):ntot, (nmrk1+1):ntot] <- mpcross2$ld$delta2
  output$r2[(nmrk1+1):ntot, (nmrk1+1):ntot] <- mpcross2$ld$r2

  colnames(output$W) <- rownames(output$W) <- colnames(output$LewontinD) <- rownames(output$LewontinD) <- colnames(output$r2) <- rownames(output$r2) <- colnames(output$delta2) <- rownames(output$delta2) <- c(colnames(mpcross1$finals), colnames(mpcross2$finals))

  return(output)
}

