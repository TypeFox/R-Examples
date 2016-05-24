#' @title Simulation function for the RepeatABEL package.
#'
#' @description  The function takes a GenABEL object as input and generates simulated phenotypic values for related individuals having repeated obserevations.
#' @param formula.FixedEffects A formula including the name of the simulated variable as response, and cofactors as fixed effects.
#' @param genabel.data A GenABEL object of class gwaa.data.
#' @param n.obs A vector including the number of observations per individual. The length of n.obs must be equal to the number if individuals in genabel.data.
#' @param SNP.eff The size of a simulated SNP.effect.
#' @param SNP.nr The SNP genotype that the SNP effect is simulated on. SNP.nr=i is the i:th SNP.
#' @param beta The simulated fixed effects. Must be equal to the number of cofactors simulated (including the intercept term).
#' @param VC A vector of length 3 including the simulated variances of the polygenic effect, permanent environmental effect and residuals, respectively.
#' @param GRM  An optional input where the Genetic Relationship Matrix can be given. Otherwise it is computed using the GenABEL package.
#' @param sim.gamma A logical parameter specifying whether the residuals shuld be simulated from a gamma distribution or not. If specified as TRUE then residuals are drawn from a gamma distribution with variance equal to the residual variance specified in \code{VC[3]} 

#' @return Returns a data frame including the simulated phenotypic values, cofactors and IDs.
#' @author Lars Ronnegard

#' @examples
#' data(gen.data)
#'  #Simulate 4 observations per individual
#'  set.seed(1234)
#'  Phen.Sim <- simulate_PhenData(y ~ age, genabel.data=gen.data, 
#'                  n.obs=rep(4, nids(gen.data)), SNP.eff=1, SNP.nr=1000, VC=c(1,1,1))
#'  GWAS1 <- rGLS(y ~ age, genabel.data = gen.data, phenotype.data = Phen.Sim)
#'  plot(GWAS1, main="Simulated Data Results")
#'  
simulate_PhenData <-
function(formula.FixedEffects = y~1, genabel.data, n.obs, SNP.eff = NULL, SNP.nr = NULL, beta = NULL, VC = c(1,1,1), GRM = NULL, sim.gamma = FALSE){
    #library(GenABEL)
  if (length(SNP.eff) != length(SNP.nr)) stop("The number of elements in SNP.eff and SNP.nr must be the same")
  if (is.null(SNP.eff)) {
      SNP.eff <- 0
      SNP.nr <- 1
  }
  SNP <- as.double(genabel.data@gtdata[,SNP.nr])
  rownames(SNP) <- NULL
  SNP <- SmoothSNPmatrix(SNP)
  n <- nids(genabel.data)
  if (length(n.obs) != n) stop("The number of individuals and the length of n.obs must be the same")
  N <- sum(n.obs)
  id1 <- rep(genabel.data@phdata$id, n.obs)
  id2 <- genabel.data@phdata$id
  #Construct incidence matrix for repeated observations
  indx <- numeric(N)
  for (i in 1:n) {
      indx <- indx + i * (id1 %in% id2[i])
  }
  genabel.data@phdata$temp.dummy <- 0
  n.names <- length(colnames(genabel.data@phdata))
  colnames(genabel.data@phdata) <- c(colnames(genabel.data@phdata)[-n.names], all.vars(formula.FixedEffects)[1])
  ###########
  if (is.null(GRM)) {
    autosomalMarkers <- which(chromosome(genabel.data) != "X")
    GRM <- compute.GRM(genabel.data[ , snpnames(genabel.data)[autosomalMarkers]])
  }
  eig <- eigen(GRM)
  if (max(diag(GRM)) > 1.6) print("There seems to be highly inbred individuals in your data")
  if (min(eig$values < -0.5)) print("The genetic relationship matrix is far from positive definite")
  non_zero.eigenvalues <- eig$values>(1e-6) #Put numerically small eigenvalues to zero
  eig$values[!non_zero.eigenvalues] <- 0
  Z.GRM <- ( eig$vectors %*% diag(sqrt(eig$values)) )[indx, ]
  rownames(Z.GRM) <- NULL
  ############
  X0 <- model.matrix(formula.FixedEffects, data = genabel.data@phdata)
  row.names(X0) <- NULL
  if (is.null(beta)) beta <- rep(0, ncol(X0))
  if (ncol(X0) != length(beta)) stop("The length of beta must be the same as the number columns in the cofactor design matrix")
  X <- X0[indx, 1:ncol(X0)]
  if (class(X)=="numeric") X <- matrix(X, length(X), 1)
  sigma2a <- VC[1]
  sigma2p <- VC[2]
  sigma2e <- VC[3]
  a <- rnorm(n, 0, sqrt(sigma2a))
  p <- rnorm(n, 0, sqrt(sigma2p))
  e <- rnorm(N, 0, sqrt(sigma2e))
  if (sim.gamma) e <- rgamma(N, sigma2e) - sigma2e #Follows a gamma distribution of variance sigma2e and mean=0
  y0 <- X %*% beta + Z.GRM %*% a + p[indx] + e
  if (length(SNP.nr)==1) y <- y0 + SNP[indx,] * SNP.eff
  if (length(SNP.nr)>1) y <- y0 + SNP[indx,] %*% SNP.eff
  id <- rep(genabel.data@phdata$id, n.obs)
  for (i in 1:length(all.vars(formula.FixedEffects))) {
    var.name <- all.vars(formula.FixedEffects)[i]
    if (i == 1) xx <- y
    if (i > 1) xx <- cbind(xx, rep(genabel.data@phdata[ , names(genabel.data@phdata) %in% var.name], n.obs))
  }
  colnames(xx) <- all.vars(formula.FixedEffects)
  Phen.Data <- data.frame(xx, id)
  return(Phen.Data)
}
