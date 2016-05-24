##################################################
#### AUTHOR:     Arnost Komarek               ####
####             (07/12/2006)                 ####
####                                          ####
#### FILE:       bayessurvreg3.checkrho.R     ####
####                                          ####
#### FUNCTIONS:  bayessurvreg3.checkrho       ####
##################################################

### ======================================
### bayessurvreg3.checkrho
### ======================================
bayessurvreg3.checkrho <- function(rho, doubly)
{
  if (!doubly){
    rho <- list(type.update = "fixed.zero", init=0, sigmaL=0.1)
  }  
  
  if (missing(rho)){
    rho <- list(type.update = "fixed.zero", init=0, sigmaL=0.1)
  }  
  
  if(length(rho) == 0) inrho <- "arnost"
  else                 inrho <- names(rho)

  itype <- match("type.update", inrho, nomatch=NA)
  if(is.na(itype)) rho$type.update <- "fixed.zero"
  Type <- pmatch(rho$type.update, table=c("fixed.zero", "normal.around.mode", "langevin"), nomatch=NA)
  if (is.na(Type)) stop("rho$type.update was not given")

  if (Type == "fixed.zero"){
    rho$init <- 0
    rho$sigmaL <- 0.1
  }
  else{
    iinit <- match("init", inrho, nomatch=NA)
    if(is.na(iinit)) rho$init <- 0
    if (rho$init <= -1 || rho$init >= 1) stop("initial rho must lie between -1 and 1")

    if (Type == "normal.around.mode"){
      rho$sigmaL <- 0.1
    }
    else{
      sinit <- match("sigmaL", inrho, nomatch=NA)
      if(is.na(sinit)) stop("scale parameter for the Langevin update of rho must be given")
      if (rho$sigmaL <= 0) stop("scale parameter for the Langevin update of rho must be positive")
    }      
  }  

  rho$typeI <- Type - 2  ## => 0 = normal.around.mode,  1 = langevin
  return(rho)
}

