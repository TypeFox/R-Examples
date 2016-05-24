overlap <-
function(niche.par, nreps, nprob, alpha = 0.95,
                    species.names, norm.redraw = TRUE) {
  niso <- ncol(niche.par[[1]]$mu) # number of isotopes
  nspec <- length(niche.par) # number of species
  nlevels <- length(alpha) # number of levels
  if(missing(species.names)) species.names <- names(niche.par)
  # temporary variables
  mu <- matrix(NA, niso, nspec)
  x <- array(NA, dim = c(niso, nprob, nspec))
  C <- array(NA, dim = c(niso, niso, nspec))
  Zsq <- array(NA, dim = c(nprob, nspec-1, nlevels))
  # constants
  qlevels <-  qchisq(alpha, df = niso)
  qlevels <- array(rep(qlevels, each = nprob*(nspec-1)),
                   dim = c(nprob, nspec-1, nlevels))
  if(!norm.redraw) z <- matrix(rnorm(niso*nprob), niso, nprob)
  # output
  over <- array(NA, dim = c(nspec, nspec, nreps, nlevels))
  dimnames(over) <- list(species.names, species.names, NULL,
                         paste0(alpha*100, "%"))
  names(dimnames(over)) <- c("Species A", "Species B", "", "alpha")
  # subsample the parameters posteriors of each niche nreps times
  ind <- sapply(niche.par, function(mv) {
    sample(nrow(mv$mu), nreps, replace = TRUE)
  })
  # calculate each overlap
  for(jj in 1:nreps) {
    # get means, t(chol(variance)), and draws for each species
    for(ii in 1:nspec) {
      mu[,ii] <- niche.par[[ii]]$mu[ind[jj,ii],]
      C[,,ii] <- t(chol(niche.par[[ii]]$Sigma[,,ind[jj,ii]]))
      if(norm.redraw) z <- matrix(rnorm(niso*nprob), niso, nprob)
      x[,,ii] <- C[,,ii] %*% z + mu[,ii]
    }
    # check whether the draw of each species is within the ellipse of every other
    for(ii in 1:nspec) {
      Z <- backsolve(C[,,ii], matrix(x[,,-ii], nrow = niso)-mu[,ii], upper.tri = FALSE)
      Zsq[] <- colSums(Z^2)
      over[-ii,ii,jj,] <- colMeans(Zsq < qlevels)
    }
  }
  if(dim(over)[4] == 1) over <- over[,,,1]
  over
}
