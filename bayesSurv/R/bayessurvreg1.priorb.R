####################################################
#### AUTHOR:     Arnost Komarek                 ####
####             (2004)                         ####
####                                            ####
#### FILE:       bayessurvreg1.priorb.R         ####
####                                            ####
#### FUNCTIONS:  bayessurvreg1.priorb           ####
####################################################

### ======================================
### bayessurvreg1.priorb
### ======================================
## Subfunction for 'bayessurvreg1' function
##  -> just to make it more readable
##
## Manipulation with proposal and prior specification
##  for random effects
##
bayessurvreg1.priorb <- function(prior.b, nrandom, ncluster, toler.chol)
{
  thispackage = "bayesSurv"
  #thispackage = NULL

  if (!nrandom){
    priori <- c(0, 0, 0, 0, 0, 0, 0, 0)
    priord <- c(0, 0, 0, 0, 0)
    
    pdi.b <- list(double = priord, integer = priori)
    attr(pdi.b, "prior.b") <- list()    
    return(pdi.b)
  }

  if (length(prior.b) == 0) inprior <- "arnost"
  else                      inprior <- names(prior.b)
  
  ## Type of prior for D
  ## =====================
  tmp <- match("prior.D", inprior, nomatch=NA)
  if (is.na(tmp)) prior.b$prior.D <- "inv.wishart"
  priorD <-  pmatch(prior.b$prior.D, table = c("inv.wishart", "sduniform"), nomatch = 0)
  if (priorD == 0){
    prior.b$prior.D <- "inv.wishart"
    warning("Non-matching prior.b$prior.D changed to inv.wishart")
    priorD <- 1    
  }
  if (priorD != 1 & nrandom > 1){
    prior.b$prior.D <- "inv.wishart"
    warning("Dimension of the random effect > 1: prior.b$prior.D changed to inv.wishart")
  }    
  priorD <- priorD - 1     ## now: 0 = inv.wishart, 1 = sduniform
  
  ## Type of proposal etc.
  ## =====================
  tmp <- match("type.upd", inprior, nomatch=NA)
  if (is.na(tmp)) prior.b$type.upd <- "gibbs"
  if (length(prior.b$type.upd) != 1) stop("Incorrect prior.b$type.upd parameter.")
  typeUpd <- pmatch(prior.b$type.upd, table = c("random.walk.metropolis", "adaptive.metropolis", "gibbs"), nomatch = 0)
  if (typeUpd == 2 | typeUpd == 0){
    warning("Update type for random effects changed to 'gibbs'.")
    typeUpd <- 3
  }    
  typeUpd <- typeUpd - 1            ## now: 0 = random walk, 1 = adaptive, 2 = gibbs

  if (typeUpd == 2){
    prior.b$blocks <- NULL
    prior.b$weight.unif <- NULL
    prior.b$half.range.unif <- NULL

    nBlocks <- 1
    nInBlock <- nrandom
    lcovparLV <- 0.5*nInBlock*(nInBlock+1)
    indBlockLV <- 1:nrandom
    covparLV <- numeric(lcovparLV)
    halfRangeUnif <- numeric(nrandom)
    weightUnif <- numeric(nBlocks)    
  }
  else{
    tmp <- match("blocks", inprior, nomatch=NA)
    if (is.na(tmp)) stop("Blocks for random effects must be specified.")
    inblocks <- names(prior.b$blocks)
    tmp <- match("ind.block", inblocks, nomatch=NA)
    if (is.na(tmp)) stop("Blocks for random effects must be specified.")
    tmp <- match("cov.prop", inblocks, nomatch=NA)
    if (is.na(tmp)) stop("Covariance matrices for proposal for random effects must be specified.")

    nBlocks <- length(prior.b$blocks$ind.block)
    if (nBlocks != length(prior.b$blocks$cov.prop)) stop("Not all random effects blocks have defined a covariance matrix for the proposal")

    nInBlock <- sapply(prior.b$blocks$ind.block, length)
    if (sum(nInBlock) != nrandom) stop("Some random effects are not assigned to any block or to more blocks.")

    indBlockLV <- unlist(prior.b$blocks$ind.block)
    if (length(indBlockLV) != nrandom) stop("Some b parameters are assigned to more blocks")
    if (sum(indBlockLV %in% 1:nrandom) != nrandom) stop("Incorrect prior.b$ind.block parameter.")

    covparLV <- unlist(prior.b$blocks$cov.prop)
    lcovparLV <- sum(0.5*nInBlock*(nInBlock+1))
    if (sum(is.na(covparLV))) stop("Incorrect prior.b$cov.prop parameter.")
    if (length(covparLV) != lcovparLV) stop("Incorrect prior.b$cov.prop parameter.")

    tmp <- match("weight.unif", inprior, nomatch=NA)
    if (is.na(tmp)) prior.b$weight.unif <- rep(0.5, nBlocks)
    if (sum(is.na(prior.b$weight.unif))) stop("Incorrect prior.b$weight.unif parameter.")
    if (length(prior.b$weight.unif) != nBlocks) stop("Incorrect prior.b$weight.unif parameter.")  
    prior.b$weight.unif[prior.b$weight.unif < 0] <- 0
    prior.b$weight.unif[prior.b$weight.unif > 1] <- 1  
    weightUnif <- prior.b$weight.unif

    tmp <- match("half.range.unif", inprior, nomatch=NA)
    if (is.na(tmp)) stop("half.range.unif for random effects must be specified.")
    if (sum(is.na(prior.b$half.range.unif))) stop("Incorrect prior.b$half.range.unif parameter.")
    if (length(prior.b$half.range.unif) != nrandom) stop("Incorrect prior.b$half.range.unif parameter.")  
    halfRangeUnif <- prior.b$half.range.unif
  }

  ## Parameters of prior for covariance matrix D
  ## ===========================================
  tmp <- match("df.D", inprior, nomatch=NA)
  if (is.na(tmp)) prior.b$df.D <- nrandom + 2
  tmp <- match("scale.D", inprior, nomatch=NA)
  if (is.na(tmp)){
    if (priorD == 1) stop("Upper limit of the uniform prior for std. dev. of the random effect must be given.")
    prior.b$scale.D <- 0.002*diag(nrandom)
    prior.b$scale.D <- prior.b$scale.D[lower.tri(prior.b$scale.D, diag = TRUE)]
  }
  if (priorD == 0){
    if (prior.b$df.D <= nrandom - 1) stop("Too low prior degrees of freedom for D matrix.")
    indD <- 0:(nrandom - 1)
    diagI <- (indD * (2*nrandom - indD + 1)) / 2
    lD <- (nrandom * (nrandom + 1))/2
    if (length(prior.b$scale.D) != lD) stop("Incorrect prior scale matrix for D matrix.")
    cholD <- .C("cholesky", as.double(prior.b$scale.D), rank = integer(1), as.integer(nrandom), as.integer(diagI), as.double(toler.chol),
                PACKAGE = thispackage)
    if (cholD$rank < nrandom) stop("Prior scale matrix for D matrix is not positive definite.")
  }
  if (priorD == 1){
    prior.b$scale.D <- prior.b$scale.D[1]
    if (prior.b$scale.D[1] <= 0) stop("Upper limit of the uniform prior for std. dev. of the random effect must be positive.")
    prior.b$df.D <- 1/(prior.b$scale.D*prior.b$scale.D)
  }    

  priordD <- c(prior.b$df.D, prior.b$scale.D)
  if (priorD == 0) names(priordD) <- c("df.D", paste("scale.D", 1:lD))
  if (priorD == 1) names(priordD) <- c("1/B.sq", "B")  

  ## Put everything to long vectors
  ##  for indBlockLV, do R --> C++ transformation  
  ## =============================================
  priori <- c(nrandom, ncluster, priorD, typeUpd, nBlocks, nInBlock, lcovparLV, indBlockLV - 1)
  priord <- c(priordD, covparLV, halfRangeUnif, weightUnif)

  pdi.b <- list(double = priord, integer = priori)
  attr(pdi.b, "prior.b") <- prior.b

  return(pdi.b)
}  
