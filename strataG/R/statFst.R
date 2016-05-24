#' @rdname popStructStat
#' @export
#' 
statFst <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {   
  if(is.null(strata.mat)) g <- g[, , strataNames(g)]
  if(ploidy(g) == 1 | nStrata(g) == 1) {
    g@sequences <- NULL
    result <- statPhist(g, nrep = nrep, strata.mat = strata.mat, keep.null = keep.null)
    result$stat.name <- "Fst"
    return(result)
  }
  
  if(nStrata(g) == 1) {
    return(list(
      stat.name = "Fst", 
      result = c(estimate = NA, p.val = NA),
      null.dist = NULL
    ))
  }
  
  strata.mat <- .checkStrataMat(strata.mat, g, nrep)
  
  result <- statFst_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    strata.mat, ploidy(g)
  )
  
  return(.formatResult(result, "Fst", keep.null))
   
#   # returns numerator and denominator sums for
#   #   theta-w calculation for alleles at each locus (page 1363)
#   locus.sums <- matrix(0, nrow = 2, ncol = nLoc(g))
#   for(i in 1:ncol(locus.sums)) {
#     loc.mat <- matrix(g@loci[, i], ncol = ploidy(g))
#     
#     # get rid of NAs
#     to.use <- apply(cbind(strata, loc.mat), 1, function(x) !any(is.na(x)))
#     if(sum(to.use) == 0) next
#     strata.vec <- strata[to.use]
#     loc.mat <- loc.mat[to.use, ]
#     
#     # identify unique alleles and return zeroes if locus is fixed
#     alleles.vec <- unique(c(loc.mat))
#     if(length(alleles.vec) < 2) next
#     
#     # calculate variables constant for locus
#     nvec <- table(strata.vec)
#     nvec <- nvec[nvec > 0]
#     r <- length(nvec)
#     if(r < 2) next
#     nbar <- mean(nvec)
#     rnbar <- r * nbar
#     nc <- (rnbar - (sum(nvec ^ 2) / rnbar)) / (r - 1)
#     
#     # get allele frequencies in each population
#     allele.freq <- prop.table(table(c(loc.mat), rep(strata.vec, ploidy(g))), 2)
#     
#     # get proportion heterozygosity for each population (rows) and allele (cols)
#     is.het <- apply(loc.mat, 1, function(x) {
#       alleles.vec %in% x & length(unique(x)) > 1
#     })
#     rownames(is.het) <- alleles.vec
#     p.het <- do.call(rbind, tapply(1:ncol(is.het), strata.vec, function(i) {
#       rowMeans(is.het[, i])
#     }, simplify = FALSE))
#     
#     # create matrix of Va, Vb, Vc for all alleles
#     varcomp.mat <- matrix(nrow = 3, ncol = length(alleles.vec))
#     rownames(varcomp.mat) <- c("Va", "Vb", "Vc")
#     colnames(varcomp.mat) <- alleles.vec
#     for(a in alleles.vec) {      
#       pbar <- sum(nvec * allele.freq[a, ]) / rnbar
#       hbar <- sum(nvec * p.het[, a]) / rnbar
#       # s = numerator inside summation for s2 calculation
#       s <- nvec * (allele.freq[a,] - pbar) ^ 2
#       s2 <- sum(s) / (r - 1) / nbar
#       
#       # used in equations for Va and Vb
#       inner.term <- (pbar * (1 - pbar)) - ((r - 1) * s2 / r)
#       nbar.m1 <- nbar - 1
#       
#       # Va (between strata) - Eqn. 2
#       Va1 <- 0.25 * hbar
#       Va2 <- s2 - (inner.term - Va1) / nbar.m1
#       Va <- nbar * Va2 / nc
#       
#       # Vb (between individuals within strata) - Eqn. 3
#       Vb1 <- ((2 * nbar - 1) * hbar) / 4 / nbar
#       Vb2 <- inner.term - Vb1
#       Vb <- nbar * Vb2 / nbar.m1
#       
#       # Vc (between gametes within individuals) - Eqn. 4
#       Vc <- 0.5 * hbar
#       
#       varcomp.mat[, a] <- c(Va, Vb, Vc)
#     } # rows are variance components (1=Va, 2=Vb, 3=Vc), columns are alleles
#     
#     # calculate three parameters from Eqn. 1
#     # F <- 1 - sum(varcomp.mat[3,]) / sum(varcomp.mat)
#     # theta <- sum(varcomp.mat[1,] / sum(varcomp.mat)
#     # f <- 1 - sum(varcomp.mat[3,]) / sum(varcomp.mat[c(2,3),])
#     
#     locus.sums[, i] <- c(sum(varcomp.mat["Va", ]), sum(varcomp.mat))
#   } # row 1 is sum of Va, row 2 is sum of Va+Vb+Vc
#   
#   # return W & C Fst estimate (theta-w) from Eqn. 10
#   est <- sum(locus.sums[1, ], na.rm = TRUE) / sum(locus.sums[2, ], na.rm = TRUE)
#   
#   if(is.nan(est)) est <- NA
#   c(Fst = est)
}


#' @rdname popStructStat
#' @export
#' 
statFstPrime <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {  
  if(is.null(strata.mat)) g <- g[, , strataNames(g)]
  if(ploidy(g) == 1 | nStrata(g) == 1) {
    return(list(
      stat.name = "F'st", 
      result = c(estimate = NA, p.val = NA),
      null.dist = NULL
    ))
  }
  
  strata.mat <- .checkStrataMat(strata.mat, g, nrep)
  
  result <- statFstPrime_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    strata.mat, ploidy(g)
  )
  
  return(.formatResult(result, "F'st", keep.null))
  
  
#   loci.max <- apply(loci.fst, 2, function(x) {
#     new.locus <- paste(strata, x, sep = ".")
#     new.locus[is.na(x) | is.na(strata)] <- NA
#     as.numeric(factor(new.locus)) - 1
#   })
#   
#   result <- statFst_C(loci.fst, strata, ploidy) / statFst_C(loci.max, strata, ploidy)
#   
#   
  
#   
#   fst.max.g <- g
#   for(i in 1:ncol(fst.max.g@loci)) {
#     new.locus <- paste(strata.vec, fst.max.g@loci[, i], sep = ".")
#     new.locus[is.na(g@loci[, i]) | is.na(strata.vec)] <- NA
#     fst.max.g@loci[, i] <- factor(new.locus)
#   }
#   est <- statFst(g, strata, ...) / statFst(fst.max.g, strata, ...)
#     
#   names(est) <- "F'st"
#   est
}