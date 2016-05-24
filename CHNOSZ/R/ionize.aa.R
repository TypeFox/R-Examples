# CHNOSZ/ionize.aa.R
# rewritten ionization function 20120526 jmd

ionize.aa <- function(aa, property="Z", T=25, P="Psat", pH=7, ret.val=NULL, suppress.Cys=FALSE) {
  # calculate the additive ionization property of proteins with amino acid 
  # composition in aa as a function of vectors of T, P and pH;
  # property if NULL is the net charge, if not NULL is one of the subcrt() properties
  # or "A" to calculate A/2.303RT for the ionization reaction
  # ret.val can be 'pK', 'alpha' or 'aavals' to get these values;
  # T, P and pH should be same length, or larger a multiple of the smaller
  lmax <- max(c(length(T), length(P), length(pH)))
  T <- rep(T, length.out=lmax)
  P <- rep(P, length.out=lmax)
  pH <- rep(pH, length.out=lmax)
  # turn pH into a matrix with as many columns as ionizable groups
  pH <- matrix(rep(pH, 9), ncol=9)
  # turn charges into a matrix with as many rows as T,P,pH conditions
  charges <- c(-1, -1, -1, 1, 1, 1, -1, 1, -1)
  charges <- matrix(rep(charges, lmax), nrow=lmax, byrow=TRUE)
  # the rownumbers of the ionizable groups in thermo$obigt
  neutral <- c("[Cys]", "[Asp]", "[Glu]", "[His]", "[Lys]", "[Arg]", "[Tyr]", "[AABB]", "[AABB]")
  charged <- c("[Cys-]","[Asp-]","[Glu-]","[His+]","[Lys+]","[Arg+]","[Tyr-]","[AABB+]","[AABB-]")
  ineutral <- info(neutral, "aq")
  icharged <- info(charged, "aq")
  # we'll only call subcrt() with the unique pressure/temperature combinations
  pTP <- paste(T, P)
  dupPT <- duplicated(pTP)
  # what property are we after
  sprop <- c("G", property)
  if(property %in%  c("A", "Z")) sprop <- "G"
  sout <- subcrt(c(ineutral, icharged), T=T[!dupPT], P=P[!dupPT], property=sprop)$out
  # the G-values
  Gs <- sapply(sout, function(x) x$G)
  # keep it as a matrix even if we have only one unique T, P-combo
  if(length(pTP[!dupPT])==1) Gs <- t(Gs)
  # now the Gibbs energy difference for each group
  DG <- Gs[, 10:18, drop=FALSE] - Gs[, 1:9, drop=FALSE]
  # build a matrix with one row for each of the (possibly duplicated) T, P values
  uPT <- unique(pTP)
  DG <- t(sapply(pTP, function(x) DG[match(x, uPT), , drop=FALSE]))
  # the pK values (-logK) 
  DG <- DG * charges
  pK <- apply(DG, 2, function(x) convert(x, "logK", T=convert(T, "K")))
  # keep it as a matrix even if we have only one T, P-combo
  if(lmax==1) pK <- t(pK)
  if(identical(ret.val, "pK")) {
    colnames(pK) <- charged
    return(pK)
  }
  # now to calculate alpha! - degrees of formation of the charged groups
  alpha <- 1 / (1 + 10 ^ (charges * (pH - pK)))
  # suppress cysteine ionization if requested
  if(suppress.Cys) alpha[, 1] <- 0
  if(identical(ret.val, "alpha")) return(alpha)
  # now to calculate the properties of the ionizable groups - can be charges, 
  # the chemical affinities of the ionization reactions,
  # or another property from subcrt()
  if(identical(property, "Z")) aavals <- charges
  else if(identical(property, "A")) aavals <- - charges * (pH - pK)
  else {
    # it's not charge, so compile it from the subcrt output
    # the property-values
    icol <- match(property, colnames(sout[[1]]))
    aavals <- sapply(sout, function(x) x[,icol])
    # keep it as a matrix even if we have only one unique T, P-combo
    if(length(pTP[!dupPT])==1) aavals <- t(aavals)
    # build a matrix with one row for each of the (possibly duplicated) T, P values
    aavals <- t(sapply(pTP, function(x) aavals[match(x, uPT), , drop=FALSE], USE.NAMES=FALSE))
    # the property difference for each group
    aavals <- aavals[, 10:18, drop=FALSE] - aavals[, 1:9, drop=FALSE]
  }
  if(identical(ret.val, "aavals")) {
    colnames(aavals) <- charged
    return(aavals)
  }
  # the contribution from each group to the ionization property of the protein
  aavals <- aavals * alpha
  # now work with 'aa'; the next line is so that a missing argument shows
  # the name of this function in the error message
  aa <- aa
  # the columns where we find the counts of ionizable groups
  iionize <- match(c("Cys", "Asp", "Glu", "His", "Lys", "Arg", "Tyr", "chains", "chains"), colnames(aa))
  aa <- as.matrix(aa[, iionize])
  # add it all up
  out <- apply(aa, 1, function(x) {
    aavals %*% x
  })
  # keep it as a matrix even if we have only one T, P-combo
  if(lmax==1) out <- t(out)
  rownames(out) <- rownames(aavals)
  # that's all folks!
  return(out)
}
