compute_bnrf <-
function(mpcross1, mpcross2, r, grid=FALSE)
{
  output <- list()
  nmrk1 <- ncol(mpcross1$founders)
  nmrk2 <- ncol(mpcross2$founders)
  n.finals <- length(mpcross1$id)

  # pairs of loci
  pairs <- expand.grid(1:nmrk1, (nmrk1+1):(nmrk1+nmrk2))

  n.pairs <- nrow(pairs)

  # order by rows
  pairs <- pairs[order(pairs[,1], pairs[,2]), ]

  # combine two sets of data
  finals <- cbind(mpcross1$id, mpcross1$finals, mpcross2$finals)
  founders <- cbind(mpcross1$founders, mpcross2$founders)

  pedigree <- mpcross1$pedigree 

  rpairs <- CR_estrf(finals, founders, pedigree, pairs, r)

  if (grid) output$lkhdgrid <- cbind(pairs, rpairs)

  maxlkhd <- apply(rpairs, 1, max)
  minlkhd <- apply(rpairs, 1, min)

  same <- ((maxlkhd-minlkhd)<.000001)
  tm <- r[apply(rpairs, 1, which.max)]
  tm[which(same==TRUE)] <- NA

  output$theta <- matrix(data=tm, nrow=nmrk1, ncol=nmrk2, byrow=TRUE)
 
  theta.5 <- which(r==0.5)
  if (grid) output$lodgrid <- cbind(pairs, rpairs-rpairs[,theta.5])
  lm <- apply(rpairs, 1, max)
  output$lkhd <- matrix(data=lm, nrow=nmrk1, ncol=nmrk2, byrow=TRUE)
  lm <- apply(rpairs-rpairs[,theta.5], 1, max)
  lm[which(same==TRUE)] <- NA
  output$lod <- matrix(data=lm, nrow=nmrk1, ncol=nmrk2, byrow=TRUE)
  output$rgrid <- r

  return(output)
}

