combine_rf <-
function(mpcross1, mpcross2, r)
{
  output <- list()

  rfb <- compute_bnrf(mpcross1, mpcross2, r)

  # now combine the three pieces
  nmrk1 <- ncol(mpcross1$founders)
  nmrk2 <- ncol(mpcross2$founders)
  ntot <- nmrk1+nmrk2

  output$rgrid <- r
  output$theta <- matrix(nrow=ntot, ncol=ntot)
  output$lod <- matrix(nrow=ntot, ncol=ntot)
  output$lkhd <- matrix(nrow=ntot, ncol=ntot)

  output$theta[1:nmrk1, 1:nmrk1] <- mpcross1$rf$theta
  output$lkhd[1:nmrk1, 1:nmrk1] <- mpcross1$rf$lkhd
  output$lod[1:nmrk1, 1:nmrk1] <- mpcross1$rf$lod

  output$theta[1:nmrk1, (nmrk1+1):ntot] <- (rfb$theta)
  output$lkhd[1:nmrk1, (nmrk1+1):ntot] <- (rfb$lkhd)
  output$lod[1:nmrk1, (nmrk1+1):ntot] <- (rfb$lod)

  output$theta[(nmrk1+1):ntot, 1:nmrk1] <- t(rfb$theta)
  output$lkhd[(nmrk1+1):ntot, 1:nmrk1] <- t(rfb$lkhd)
  output$lod[(nmrk1+1):ntot, 1:nmrk1] <- t(rfb$lod)

  output$theta[(nmrk1+1):ntot, (nmrk1+1):ntot] <- mpcross2$rf$theta
  output$lkhd[(nmrk1+1):ntot, (nmrk1+1):ntot] <- mpcross2$rf$lkhd
  output$lod[(nmrk1+1):ntot, (nmrk1+1):ntot] <- mpcross2$rf$lod

  # modify the parts corresponding to the second set of data
#  mpcross2$rf$lkhdgrid[,1:2] <- mpcross2$rf$lkhdgrid[,1:2]+nmrk1
#  mpcross2$rf$lodgrid[,1:2] <- mpcross2$rf$lodgrid[,1:2]+nmrk1
#  rfb$lkhdgrid[,2] <- rfb$lkhdgrid[,2]+nmrk1
#  rfb$lodgrid[,2] <- rfb$lodgrid[,2]+nmrk1
   # not working correctly at this time; need to fix up.
#  output$lkhdgrid <- rbind(mpcross1$rf$lkhdgrid, rfb$lkhdgrid, mpcross2$rf$lkhdgrid) 
#  output$lodgrid <- rbind(mpcross1$rf$lodgrid, rfb$lodgrid, mpcross2$rf$lodgrid)
#  output$markers <- c(mpcross1$rf$markers, mpcross2$rf$markers+nmrk1)
#  output$mrk.names <- c(mpcross1$rf$mrk.names, mpcross2$rf$mrk.name)

  colnames(output$theta) <- rownames(output$theta) <- colnames(output$lod) <- rownames(output$lod) <- colnames(output$lkhd) <- rownames(output$lkhd) <- c(colnames(mpcross1$finals), colnames(mpcross2$finals))

  return(output)
}

