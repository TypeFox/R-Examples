####################################################
#### AUTHOR:     Arnost Komarek                 ####
####             (2004)                         ####
####                                            ####
#### FILE:       bayessurvreg1.priorBeta.R      ####
####                                            ####
#### FUNCTIONS:  bayessurvreg1.priorBeta        ####
####################################################

### ======================================
### bayessurvreg1.checkStore
### ======================================
## Subfunction for 'bayessurvreg1' function
##  -> just to make it more readable
##
## Manipulation with the prior and proposal specification
##  for the regression parameters and means of random effects
##
##
## factors ..... a vector of length nX with 0 on places of columns which are not factors
##               and with 1, 2, ... on places of factors
##               * all columns of X corresponding to one factor have same index
##               * the first factor has teh highest number
## n.factors ... number of factor variables
##
bayessurvreg1.priorBeta <- function(prior.beta, nX, indb, factors, n.factors, n.in.factors)
{
  if(length(prior.beta) == 0) inprior <- "arnost"
  else                        inprior <- names(prior.beta)
  
  if (!nX){
    priori <- c(0, 0, 0, 0, 0, 0, 0)
    priord <- c(0, 0, 0, 0, 0, 0, 0, 0)
    
    pdi.beta <- list(double = priord, integer = priori)
    attr(pdi.beta, "prior.beta") <- list()    
    return(pdi.beta)
  }    

  tmp <- match("mean.prior", inprior, nomatch=NA)
  if(is.na(tmp)) stop("Prior means for betas must be specified.")
  tmp <- match("var.prior", inprior, nomatch=NA)
  if(is.na(tmp)) stop("Prior vars for betas must be specified.")
  if (length(prior.beta$mean.prior) != nX) stop("Incorrect length of a vector of prior means for betas.")
  if (length(prior.beta$var.prior) != nX) stop("Incorrect length of a vector of prior vars for betas.")  
  if (sum(is.na(prior.beta$mean.prior))) stop("Prior means for betas must not be missing.")
  if (sum(is.na(prior.beta$var.prior))) stop("Prior vars for betas must not be missing.")  
  if (sum(prior.beta$var.prior <= 0)) stop("Prior vars for betas must be all positive.")    

  tmp <- match("blocks", inprior, nomatch=NA)
  if (is.na(tmp)) tmp2 <- NA
  else            tmp2 <- match("ind.block", names(prior.beta$blocks), nomatch=NA)
  if (is.na(tmp) | is.na(tmp2)){

    ## One block with fixed effects, one block with means of random effects
    prior.beta$blocks <- list()
    fixed <- (1:nX)[indb == -1]
    random <- (1:nX)[indb > 0]

    prior.beta$blocks$ind.block <- list()
    prior.beta$blocks$cov.prop <- list()
    nBlocks <- 0
    if (length(fixed) > 0){
      nBlocks <- nBlocks + 1
      prior.beta$blocks$ind.block[[nBlocks]] <- fixed
      prior.beta$blocks$cov.prop[[nBlocks]] <- numeric((length(fixed) * (1 + length(fixed)))/2)
    }
    if (length(random) > 0){
      nBlocks <- nBlocks + 1
      prior.beta$blocks$ind.block[[nBlocks]] <- random
      prior.beta$blocks$cov.prop[[nBlocks]] <- numeric((length(random) * (1 + length(random)))/2)
    }      
    prior.beta$type.upd <- rep("gibbs", nBlocks)
    inprior <- names(prior.beta)    
    
    ## ===== OLDER VERSION =====
    ## Each parameter in one block, except factors
    ## cov.prop = diag(prior variance) for each block
##    prior.beta$blocks <- list()
##    nBlocks <- length(n.in.factors)
##    nInBlock <- n.in.factors
##    cumn <- c(0, cumsum(nInBlock))

##    prior.beta$blocks$ind.block <- list()
##    prior.beta$blocks$cov.prop <- list()
##    for (i in 1:nBlocks){
##      prior.beta$blocks$ind.block[[i]] <- (cumn[i]+1):cumn[i+1]
##      covmat <- diag(prior.beta$var.prior[(cumn[i]+1):cumn[i+1]], nrow = nInBlock[i])      
##      prior.beta$blocks$cov.prop[[i]] <- covmat[lower.tri(covmat, diag = TRUE)]
##    }    

  }

  if (is.na(tmp)) tmp2 <- NA
  else            tmp2 <- match("cov.prop", names(prior.beta$blocks), nomatch=NA)  
  if (is.na(tmp2)){
    prior.beta$blocks$cov.prop <- list()
  }     
  
  nBlocks <- length(prior.beta$blocks$ind.block)
  tmp <- match("type.upd", inprior, nomatch=NA)
  if(is.na(tmp)) prior.beta$type.upd <- rep("gibbs", nBlocks)
  if (length(prior.beta$type.upd) != nBlocks) stop("Incorrect prior.beta$type.upd parameter.")

  typeUpd <- pmatch(prior.beta$type.upd, table = c("random.walk.metropolis", "adaptive.metropolis", "gibbs"), nomatch = 0, duplicates.ok = TRUE)
  typeUpd[typeUpd == 0] <- 3        ## no matching ==> update using gibbs
  typeUpd <- typeUpd - 1            ## now: 0 = random walk, 1 = adaptive, 2 = gibbs
  
  for (i in 1:nBlocks){
    if (typeUpd[i] == 2){
      lb <- length(prior.beta$blocks$ind.block[[i]])
      prior.beta$blocks$cov.prop[[i]] <- numeric((lb*(1+lb))/2)
    }      
  }    
  
  if (nBlocks != length(prior.beta$blocks$cov.prop)) stop("Not all beta blocks have defined a covariance matrix for the proposal")

  nInBlock <- sapply(prior.beta$blocks$ind.block, length)
  if (sum(nInBlock) != nX) stop("Some beta parameters are not assigned to any block or to more blocks.")

  indBlockLV <- unlist(prior.beta$blocks$ind.block)
  if (length(indBlockLV) != nX) stop("Some beta parameters are assigned to more blocks")
  if (sum(indBlockLV %in% 1:nX) != nX) stop("Incorrect prior.beta$ind.block parameter.")

  covparLV <- unlist(prior.beta$blocks$cov.prop)
  lcovparLV <- sum(0.5*nInBlock*(nInBlock+1))
  if (sum(is.na(covparLV))) stop("Incorrect prior.beta$cov.prop parameter.")
  if (length(covparLV) != lcovparLV) stop("Incorrect prior.beta$cov.prop parameter.")

  tmp <- match("mean.sampled", inprior, nomatch=NA)
  if(is.na(tmp)) prior.beta$mean.sampled <- rep(0, nX)
  if (length(prior.beta$mean.sampled) != nX) stop("Incorrect prior.beta$mean.sampled parameter.")
  if (sum(is.na(prior.beta$mean.sampled))) stop("Incorrect prior.beta$mean.sampled parameter.")
  meanSampled <- prior.beta$mean.sampled

  tmp <- match("eps.AM", inprior, nomatch=NA)
  if(is.na(tmp)) prior.beta$eps.AM <- rep(0.05, nBlocks)
  eps <- prior.beta$eps.AM
  if (length(eps) != nBlocks) stop("Incorrect prior.beta$eps.AM parameter.")
  if (sum(is.na(eps))) stop("Incorrect prior.beta$eps.AM parameter.")
  if (sum(eps < 0)) stop("Incorrect prior.beta$eps.AM parameter.")  

  tmp <- match("sd.AM", inprior, nomatch=NA)
  if(is.na(tmp)) prior.beta$sd.AM <- (2.4*2.4)/(1:max(nInBlock))
  sdNum <- prior.beta$sd.AM
  if (length(sdNum) < max(nInBlock)) stop("Too short prior.beta$sd.AM parameter.")
  sdNum <- sdNum[1:max(nInBlock)]
  if (sum(is.na(sdNum))) stop("Incorrect prior.beta$sd.AM parameter.")
  if (sum(sdNum < 0)) stop("Incorrect prior.beta$sd.AM parameter.")    

  tmp <- match("weight.unif", inprior, nomatch=NA)
  if(is.na(tmp)) prior.beta$weight.unif <- rep(0.5, nBlocks)
  if (sum(is.na(prior.beta$weight.unif))) stop("Incorrect prior.beta$weight.unif parameter.")
  if (length(prior.beta$weight.unif) != nBlocks) stop("Incorrect prior.beta$weight.unif parameter.")  
  prior.beta$weight.unif[prior.beta$weight.unif < 0] <- 0
  prior.beta$weight.unif[prior.beta$weight.unif > 1] <- 1  
  weightUnif <- prior.beta$weight.unif

  tmp <- match("half.range.unif", inprior, nomatch=NA)
  if(is.na(tmp)) prior.beta$half.range.unif <- 0.5*sqrt(12*prior.beta$var.prior)
  if (sum(is.na(prior.beta$half.range.unif))) stop("Incorrect prior.beta$half.range.unif parameter.")
  if (length(prior.beta$half.range.unif) != nX) stop("Incorrect prior.beta$half.range.unif parameter.")  
  halfRangeUnif <- prior.beta$half.range.unif
  
  ## Put everything to long vectors
  ##  for indBlockLV, do R --> C++ transformation
  priord <- c(prior.beta$mean.prior, prior.beta$var.prior, meanSampled, halfRangeUnif, covparLV,  weightUnif, eps, sdNum)
  priori <- c(nBlocks, nX, nInBlock, max(nInBlock), lcovparLV, indBlockLV - 1, typeUpd)

  pdi.beta <- list(double = priord, integer = priori)
  attr(pdi.beta, "prior.beta") <- prior.beta

  return(pdi.beta)  
}  
