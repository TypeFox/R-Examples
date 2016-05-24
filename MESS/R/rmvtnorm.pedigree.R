rmvtnorm.pedigree <- function(n=1, pedigree, h2=0, c2=0, d2=0) {

# simulate.pedigree

  if (! "pedigree" %in% class(pedigree))
    stop("requires a single pedigree as input")
  
  if ((sum(h2, c2, d2)>=1) | min(h2, c2, d2) <0)
    stop("input parameters wrong. Must me non-negative and sum to less than 1")

  # Twins are automatically handled by the kinship function from kinship2 if coded correctly in the
  # pedigree
  
  if (d2!=0)
    stop("cannot handle d2 != 0 yet")

  # Could be programmed easily for non-inbred pedigrees using the fact that Delta7 = phi_km phi_ln + phi_kn phi_lm
  # (see Lange p. 80).

  pedsize <- length(pedigree$id)
  
  # VCmat per block
  vcmat <- (1-h2-c2)*diag(pedsize) + h2*2*kinship(pedigree) + c2
  t(rmvnorm(n, mean=rep(0, pedsize), sigma=vcmat))
}

rmvt.pedigree <- function(n=1, pedigree, h2=0, c2=0, d2=0, df=1) {

# simulate.pedigree

  if (! "pedigree" %in% class(pedigree))
    stop("requires a single pedigree as input")
  
  if ((sum(h2, c2, d2)>=1) | min(h2, c2, d2) <0)
    stop("input parameters wrong. Must me non-negative and sum to less than 1")

  # Twins are automatically handled by the kinship function from kinship2 if coded correctly in the
  # pedigree
  
  if (d2!=0)
    stop("cannot handle d2 != 0 yet")

  # Could be programmed easily for non-inbred pedigrees using the fact that Delta7 = phi_km phi_ln + phi_kn phi_lm
  # (see Lange p. 80).

  pedsize <- length(pedigree$id)
  
  # VCmat per block
  vcmat <- (1-h2-c2-d2)*diag(pedsize) + h2*2*kinship(pedigree) + c2
  t(rmvt(n, df=df, sigma=vcmat))
}


#rmvnorm.matlist <- function(n, mean=0, 
#                            sigma, varlist, sigmaI=1) {

#  if (is.matrix(varlist))
#    varlist <- list(varlist)
#  if (!is.list(varlist))
#    stop("varlist must be a list")

#  if (length(sigma) != length(varlist))
#    stop("lengths of sigma and varlist must match")

  # Check that all matrices have thesame dimensions?
#  size <- nrow(varlist[[1]])

#  vcmat <- sigmaI*diag(size)

#  sapply(1:length(sigma), function(i) {vcmat <- vcmat + sigma[i]*varlist[[i]]} )
  
#  for (i in 1:length(sigma)) {
#    lapply(varlist)
#  }
#  as.vector(rmvnorm(1, mean=rep(0, pedsize), sigma=vcmat))

  
#}
