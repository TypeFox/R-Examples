####################################################
#### AUTHOR:     Arnost Komarek                 ####
####             (07/12/2006)                   ####
####                                            ####
#### FILE:       bayessurvreg3.priorinitNb.R    ####
####                                            ####
#### FUNCTIONS:  bayessurvreg3.priorinitNb      ####
####################################################

### ======================================
### bayessurvreg3.priorinitNb
### ======================================
bayessurvreg3.priorinitNb <- function(priorinit.Nb, init, init2, design, design2, doubly)
{
  if (!doubly) stop("Doubly censored data are expected if priorinit.Nb is specified")
  if (design$nrandom != 1 || design2$nrandom != 1) stop("Only random intercepts are allowed if priorinit.Nb is specified")
  if (design$ncluster != design2$ncluster) stop("Onset and time-to-event parts of the model have to have the same number of clusters if priorinit.Nb is specified")
  if (sum(abs(design$nwithin - design$nwithin))) stop("Onset and time-to-event parts of the model have to have the same number of observations within each clusters if priorinit.Nb is specified")
  
  ##### Initial values of random effects + some related parameters #####
  ##### ========================================================== #####
  if(length(init) == 0) ininit <- "arnost"
  else                  ininit <- names(init)
  if(length(init2) == 0) ininit2 <- "arnost"
  else                   ininit2 <- names(init2)

  ## Initial values of the ONSET random intercept
  tmp <- match("b", ininit, nomatch=NA)
  if(is.na(tmp)){
    init$b <- rep(0, design$ncluster)
  }
  else{
    if (length(init$b) == 0){
      init$b <- rep(0, design$ncluster)
    } 
    else{
      if (length(init$b) < design$ncluster) stop("Incorrect init$b parameter supplied.")
      init$b <- init$b[1:design$ncluster]
    }
  }  
  if (sum(is.na(init$b))) stop("Incorrect init$b parameter supplied.")

  ## Initial values of the TIME-TO-EVENT random intercept
  tmp <- match("b", ininit2, nomatch=NA)
  if(is.na(tmp)){
    init2$b <- rep(0, design2$ncluster)
  }
  else{
    if (length(init2$b) == 0){
      init2$b <- rep(0, design2$ncluster)
    } 
    else{
      if (length(init2$b) < design2$ncluster) stop("Incorrect init2$b parameter supplied.")
      init2$b <- init2$b[1:design2$ncluster]
    }
  }  
  if (sum(is.na(init2$b))) stop("Incorrect init2$b parameter supplied.")

  ## Create vectors for RandomEff32::RE initializer
  ## ==============================================
  bparmI <- c(0, 1, design$ncluster, design$nwithin)    
  bparmD <- as.vector(init$b)
  names(bparmI) <- c("type.prior", "nRandom", "nCluster", paste("Cl", 1:design$ncluster, sep=""))
  names(bparmD) <- paste("d.", 1:design$ncluster, sep="")        

  bparmI2 <- c(0, 1, design2$ncluster, design2$nwithin)    
  bparmD2 <- as.vector(init2$b)
  names(bparmI2) <- c("type.prior", "nRandom", "nCluster", paste("Cl", 1:design$ncluster, sep=""))
  names(bparmD2) <- paste("b.", 1:design2$ncluster, sep="")
  

  ##### Covariance matrix of random effects and its prior #####
  ##### ================================================= #####
  if(length(priorinit.Nb) == 0) inprinit <- "arnost"
  else                          inprinit <- names(priorinit.Nb)

  ### Initial value of the covariance matrix of random effects
  iinit <- match("init.D", inprinit, nomatch=NA)
  if(is.na(iinit)) stop("priorinit.Nb$init.D must be given")
  Dmat <- priorinit.Nb$init.D  
  if (length(Dmat) == 3){
    Dtri <- Dmat
  }
  else{
    if (is.null(dim(Dmat))) stop("priorinit.Nb$init.D must be a matrix")
    if (nrow(Dmat) != 2 || ncol(Dmat) != 2) stop("priorinit.Nb$init.D must be a matrix 2 x 2")
    Dtri <- Dmat[lower.tri(Dmat, diag=TRUE)]
  }  

  ### Degrees of freedom of the Wishart prior
  idf <- match("df.Di", inprinit, nomatch=NA)
  if(is.na(idf)) stop("priorinit.Nb$df.Di must be given")
  df <- priorinit.Nb$df.Di[1]
  if (df <= 1) stop("priorinit.Nb$df.Di must be > 1")
  
  ### Scale matrix of the Wishart prior
  iS <- match("scale.Di", inprinit, nomatch=NA)
  if(is.na(iS)) stop("priorinit.Nb$scale.Di must be given")
  Smat <- priorinit.Nb$scale.Di 
  if (length(Smat) == 3){
    Stri <- Smat
  }
  else{
    if (is.null(dim(Smat))) stop("priorinit.Nb$scale.Di must be a matrix")
    if (nrow(Smat) != 2 || ncol(Smat) != 2) stop("priorinit.Nb$scale.Di must be a matrix 2 x 2")
    Stri <- Smat[lower.tri(Smat, diag=TRUE)]
  }  

  wdim <- 2
  Imat <- diag(wdim)
  rowsI <- row(Imat)[lower.tri(row(Imat), diag=TRUE)]
  colsI <- col(Imat)[lower.tri(col(Imat), diag=TRUE) ] 
  naam <- paste("(", rowsI, ".", colsI, ")", sep="")  
  parD <- c(Dtri, df, Stri)
  names(parD) <- c(paste("D", naam, sep=""), "df", paste("S", naam, sep=""))


  reffdi  <- list(bparmI=bparmI,  bparmD=bparmD,  GsplI=0, GsplD=parD, specification=2, r=0)
  reffdi2 <- list(bparmI=bparmI2, bparmD=bparmD2, GsplI=0, GsplD=0,    specification=2, r=0)
  toreturn <- list(reffdi=reffdi, reffdi2=reffdi2)
  attr(toreturn, "priorinit.Nb") <- priorinit.Nb
  attr(toreturn, "init") <- init
  attr(toreturn, "init2") <- init2  
  
  return(toreturn)
}

