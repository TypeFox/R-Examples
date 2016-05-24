##
##  PURPOSE:   Argument manipulation for NMixRelabel functions
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    26/02/2010 (by taking sub-code originally included in NMixRelabel.NMixMCMC function)
##
##  FUNCTIONS:  NMixRelabelAlgorithm
##
## ================================================================================================

## *************************************************************
## NMixRelabelAlgorithm
## *************************************************************
##
NMixRelabelAlgorithm <- function(type=c("mean", "weight", "stephens"), par, dim)
{
  type <- match.arg(type)
  if (type == "mean"){
    Ctype <- 1
    if (missing(par)) par <- 1
    if (!is.numeric(par)) stop("par must be a number")
    par <- par[1]
    iparam <- par
    if (iparam <= 0 | iparam > dim) stop(paste("par must be between 1 and ", dim, sep=""))

    relabel <- list(type=type, par=iparam)      ## resulting re-labeling
    iparam <- iparam - 1                        ## R -> C indexing
  }else{
    if (type == "weight"){
      Ctype <- 2
      par <- 1
      iparam <- 0
      relabel <- list(type=type, par=1)         ## resulting re-labeling
    }
    else{
      if (type == "stephens"){
        Ctype <- 3
        if (missing(par)) par <- list(type.init = "identity", par = 1, maxiter = 50)
        if (!is.list(par)) stop("par must be a list")
        inpar <- names(par)
        ind_type.init  <- match("type.init", inpar, nomatch=NA)
        ind_par        <- match("par",       inpar, nomatch=NA)
        ind_maxiter    <- match("maxiter",   inpar, nomatch=NA)

        if (is.na(ind_type.init)) par$type.init <- "identity"
        Ctype.init <- pmatch(par$type.init, table=c("identity", "mean", "weight"), nomatch=NA) - 1
        if (is.na(Ctype.init)) stop("unknown par$type.init option supplied")

        if (Ctype.init == 0){               ### par$init == identity
          par$par <- 1
        }else{
          if (Ctype.init == 1){             ### par$init == mean
            if (is.na(ind_par)) par$par <- 1
            if (!is.numeric(par$par)) stop("par$par must be a number")
            par$par <- par$par[1]
            if (par$par <= 0 | par$par > dim) stop(paste("par$par must be between 1 and ", dim, sep=""))            
          }else{
            if (Ctype.init == 2){           ### par$init == weight
              par$par <- 1
            }  
          }  
        }

        if (is.na(ind_maxiter)) par$maxiter <- 50
        if (!is.numeric(par$maxiter)) stop("par$maxiter must be a number")        
        par$maxiter <- par$maxiter[1]
        if (par$maxiter <= 0) stop("par$maxiter must be positive")

        iparam <- c(Ctype.init, par$par - 1, par$maxiter, 1)              ### iparam[4] determines whether search or transportation problem should
                                                                          ### be used in step2 (0 = transportation problem, 1 = search)
        names(iparam) <- c("type.init", "par", "maxiter", "type.step2")
        relabel <- list(type=type, par=par)                               ## resulting re-labeling
      }  
    }
  }

  RET <- list(relabel=relabel, Ctype=Ctype, iparam=iparam)
  return(RET)
}

