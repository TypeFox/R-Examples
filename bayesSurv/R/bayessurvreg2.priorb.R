################################################
#### AUTHOR:     Arnost Komarek             ####
####             (2005)                     ####
####                                        ####
#### FILE:       bayessurvreg2.priorb.R     ####
####                                        ####
#### FUNCTIONS:  bayessurvreg2.priorb       ####
################################################

### ======================================
### bayessurvreg2.priorb
### ======================================
## Manipulation with the prior specification for the random effects
## and their covariance matrix
## - version for bayessurvreg2 (AFT model with G-spline error and normal random effects)
##
## 30/01/2005
## =======================================================================================
bayessurvreg2.priorb <- function(prior.b, init, design)
{
  thispackage = "bayesSurv"
  #thispackage = NULL
  
  if (design$nrandom){

    if (length(prior.b) == 0) inprior <- "arnost"
    else                      inprior <- names(prior.b)
    if(length(init) == 0) ininit <- "arnost"
    else                  ininit <- names(init)

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
    if (priorD != 1 & design$nrandom > 1){
      prior.b$prior.D <- "inv.wishart"
      warning("Dimension of the random effect > 1: prior.b$prior.D changed to inv.wishart")
    }    
    priorD <- priorD - 1     ## now: 0 = inv.wishart, 1 = sduniform

    ## Parameters of prior for covariance matrix D
    ## ===========================================
    tmp <- match("df.D", inprior, nomatch=NA)
    if (is.na(tmp)) prior.b$df.D <- design$nrandom + 2
    tmp <- match("scale.D", inprior, nomatch=NA)
    if (is.na(tmp)){
      if (priorD == 1) stop("Upper limit of the uniform prior for std. dev. of the random effect must be given.")
      prior.b$scale.D <- 0.002*diag(design$nrandom)
      prior.b$scale.D <- prior.b$scale.D[lower.tri(prior.b$scale.D, diag = TRUE)]
    }
    if (priorD == 0){
      if (prior.b$df.D <= design$nrandom - 1) stop("Too low prior degrees of freedom for D matrix.")
      indD <- 0:(design$nrandom - 1)
      diagI <- (indD * (2*design$nrandom - indD + 1)) / 2
      lD <- (design$nrandom * (design$nrandom + 1))/2
      if (length(prior.b$scale.D) != lD) stop("Incorrect prior scale matrix for D matrix.")
      cholD <- .C("cholesky", as.double(prior.b$scale.D), rank = integer(1),
                              as.integer(design$nrandom), as.integer(diagI), as.double(1e-10),
                PACKAGE = thispackage)
      if (cholD$rank < design$nrandom) stop("Prior scale matrix for D matrix is not positive definite.")
    }
    if (priorD == 1){
      prior.b$scale.D <- prior.b$scale.D[1]
      if (prior.b$scale.D[1] <= 0) stop("Upper limit of the uniform prior for std. dev. of the random effect must be positive.")
      prior.b$df.D <- 1/(prior.b$scale.D*prior.b$scale.D)
    }

    ## Initial values of the covariance matrix D
    ## =========================================
    tmp <- match("D", ininit, nomatch=NA)
    if(is.na(tmp)){                     
      init$D <- diag(design$nrandom)[lower.tri(diag(design$nrandom), diag = TRUE)]   ## identity matrix as initial D
    }
    else{
      if (length(init$D) < 0.5*design$nrandom*(1+design$nrandom)) stop("Incorrect init$D parameter supplied.")
      init$D <- init$D[1:(0.5*design$nrandom*(1+design$nrandom))]
    }
    if (sum(is.na(init$D))) stop("Incorrect init$D parameter supplied.")

    ## Create vectors for CovMatrix constructor
    ## ========================================
    DparmI <- c(design$nrandom, priorD)
    DparmD <- c(init$D, prior.b$df.D, prior.b$scale.D)
    names(DparmI) <- c("nrow", "type.prior")
    if (priorD == 0){
      Dtemp <- diag(design$nrandom)
      Drows <- row(Dtemp)[lower.tri(row(Dtemp), diag = TRUE)]
      Dcols <- col(Dtemp)[lower.tri(col(Dtemp), diag = TRUE)]      
      names(DparmD) <- c(paste("D.", Drows, ".", Dcols, sep = ""), "df.D", paste("scale.D.", Drows, ".", Dcols, sep = ""))
    }      
    if (priorD == 1){
      names(DparmD) <- c("D", "1/B.sq", "B")
    }    
    
    ## Initial values of the random effects b
    ## ========================================    
    tmp <- match("b", ininit, nomatch=NA)
    if(is.na(tmp)){
      tmp2 <- match("beta", ininit, nomatch=NA)
      if (is.na(tmp2)) stop("At least init$beta should be known at this moment")
      bb <- init$beta[design$indb > 0]
      if (design$randomInt) bb <- c(0, bb)
      init$b <- rep(bb, design$ncluster)
    } 
    else{
      if (length(init$b) == 0){
        tmp2 <- match("beta", ininit, nomatch=NA)
        if (is.na(tmp2)) stop("At least init$beta should be known at this moment")
        bb <- init$beta[design$indb > 0]
        if (design$randomInt) bb <- c(0, bb)
        init$b <- rep(bb, design$ncluster)
      } 
      else{
        if (length(init$b) < design$nrandom*design$ncluster) stop("Incorrect init$b parameter supplied.")
        init$b <- init$b[1:(design$nrandom*design$ncluster)]
      }
    }  
    if (sum(is.na(init$b))) stop("Incorrect init$b parameter supplied.")

    ## Create vectors for RandomEff constructor
    ## ========================================
    bparmI <- c(0, design$nrandom, design$ncluster, design$nwithin)    
    bparmD <- as.vector(init$b)
    names(bparmI) <- c("type.prior", "nRandom", "nCluster", paste("Cl", 1:design$ncluster, sep=""))
    names(bparmD) <- paste("b.", rep(1:design$nrandom, design$ncluster), ".",
                           rep(1:design$ncluster, rep(design$nrandom, design$ncluster)), sep="")
  }
  else{
    bparmI <- c(0, 0, design$n, rep(1, design$n))
    bparmD <- rep(0, design$n)
    names(bparmI) <- c("type.prior", "nRandom", "nCluster", paste("Cl", 1:design$n, sep=""))
    names(bparmD) <- paste("b", 1:design$n, sep="")

    DparmI <- c(0, 0)
    DparmD <- c(0, 0, 0)
    names(DparmI) <- c("nrow", "type.prior")
    names(DparmD) <- c("D", "df.D", "scale.D")

    init$b <- numeric(0)
    init$D <- numeric(0)
    prior.b <- list(pior.D = character(0), df.D = numeric(0), scale.D = numeric(0))
  }

  toreturn <- list(bparmI=bparmI, bparmD=bparmD, DparmI=DparmI, DparmD=DparmD)
  attr(toreturn, "init.b") <- init$b
  attr(toreturn, "init.D") <- init$D
  attr(toreturn, "prior.b") <- prior.b
  
  return(toreturn)  
}

