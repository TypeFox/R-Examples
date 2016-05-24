#################################################
#### AUTHOR:     Arnost Komarek              ####
####             (2005)                      ####
####                                         ####
#### FILE:       bayesBisurvreg.priorInit.R  ####
####                                         ####
#### FUNCTIONS:  bayesBisurvreg.priorInit    ####
#################################################

### ======================================
### bayesBisurvreg.priorInit
### ======================================
## Subsection of bayesBisurvreg
## -> just to make it more readable
##
## 14/01/2005
##
## dim ................. dimension of the response (1 or 2)
## design, design2 ..... having the structure as returned by function bayessurvreg.design!!!
##
## REMARK: both design$Y and design2$Y contain not-transformed times (i.e. t and not y)
##
bayesBisurvreg.priorInit <- function(dim, prior, init, design, mcmc.par, prior2, init2, design2, mcmc.par2, doubly)
{
  thispackage = "bayesSurv"
  #thispackage = NULL
  
  if (dim < 1 | dim > 2) stop("dim parameter must be either 1 or 2")
  
  ### Onset and times-to-event
  ### ========================
  if(length(init) == 0) ininit <- "arnost"
  else                  ininit <- names(init)
  if(length(init2) == 0) ininit2 <- "arnost"
  else                   ininit2 <- names(init2)
  
  tmp <- match("y", ininit, nomatch=NA)
  if(is.na(tmp)) init.y <- numeric(0)
  else           init.y <- init$y
  tmp <- match("y", ininit2, nomatch=NA)
  if(is.na(tmp)) init2.y <- numeric(0)
  else           init2.y <- init2$y
  
  yy <- give.init.y2(init.y, init2.y, dim, design, design2, doubly)    

  
  ### Prior/inits for the onset time/event time related quantities
  ### ============================================================
    ## G-spline:
    ## ---------
  gspl <- give.init.Gspline(prior, init, mcmc.par, dim)
  init <- attr(gspl, "init")
  prior <- attr(gspl, "prior")
  mcmc.par <- attr(gspl, "mcmc.par")
  if(length(init) == 0) ininit <- "arnost"
  else                  ininit <- names(init)
  
    ## Compute residuals:
    ## ------------------
  if (!design$nX){
    resid <- yy$init.y
  }
  else{
    tmp <- match("beta", ininit, nomatch=NA)
    if (is.na(tmp)) stop("init$beta should already be known...")
    resid <- matrix(as.vector(t(yy$init.y)) - design$X %*% init$beta, ncol=dim, byrow=TRUE)    
  }

    ## Initial allocations:
    ## --------------------
  tmp <- match("r", ininit, nomatch=NA)
  if(is.na(tmp)) init$r <- numeric(0)
  init$r <- give.init.r(init$r, resid, dim, prior$K, init$gamma, init$sigma, prior$c4delta, init$intercept, init$scale)
  
    ## Recalculate r into the vector of length n with entries from {1,...,total_length}
  if (dim == 1)   rsm <- init$r + prior$K[1] + 1
  else{
    if (dim == 2) rsm <- (init$r[,2]+prior$K[2])*(2*prior$K[1]+1) + init$r[,1]+prior$K[1]+1
    else          stop("Unimplemented dimension appeared in bayesBisurvreg.priorInit()")
  }  

  init$y <- yy$init.y
  
  ### Prior/inits for the time-to-event in the case of doubly censoring
  ### ==================================================================
  if (doubly){

      ## G-spline:
      ## ---------
    gspl2 <- give.init.Gspline(prior2, init2, mcmc.par2, dim)
    init2 <- attr(gspl2, "init")
    prior2 <- attr(gspl2, "prior")
    mcmc.par2 <- attr(gspl2, "mcmc.par")
    if(length(init2) == 0) ininit2 <- "arnost"
    else                   ininit2 <- names(init2)
        
      ## Compute residuals:
      ## ------------------
    if (!design2$nX){
      resid2 <- yy$init2.y
    }
    else{
      tmp <- match("beta", ininit2, nomatch=NA)
      if (is.na(tmp)) stop("init2$beta should already be known...")
      resid2 <- matrix(as.vector(t(yy$init2.y)) - design2$X %*% init2$beta, ncol=dim, byrow=TRUE)
    }

      ## Initial allocations:
      ## --------------------
    tmp <- match("r", ininit2, nomatch=NA)
    if(is.na(tmp)) init2$r <- numeric(0)
    init2$r <- give.init.r(init2$r, resid2, dim, prior2$K, init2$gamma, init2$sigma, prior2$c4delta, init2$intercept, init2$scale)
  
      ## Recalculate r into the vector of length n with entries from {1,...,total_length}
    if (dim == 1)   rsm2 <- init2$r + prior2$K[1] + 1
    else{
      if (dim == 2) rsm2 <- (init2$r[,2]+prior2$K[2])*(2*prior2$K[1]+1) + init2$r[,1]+prior2$K[1]+1
      else          stop("Unimplemented dimension appeared in bayesBisurvreg.priorInit()")
    }

    init2$y <- yy$init2.y
  }
  else{    ## not doubly:
    gspl2 <- list(Gparmi = 0, Gparmd = 0, specification = 1)
    rsm2 <- 0
  }

  toreturn <- list(Gparmi = gspl$Gparmi, Gparmd = gspl$Gparmd, 
                   y = t(yy$init.y), r = rsm,
                   Gparmi2 = gspl2$Gparmi, Gparmd2 = gspl2$Gparmd,
                   y2 = t(yy$init2.y), r2 = rsm2,
                   iter = init$iter, specification = c(gspl$specification, gspl2$specification),
                   y.left=yy$y1.left, y.right=yy$y1.right, status=yy$status1,
                   t2.left=yy$t2.left, t2.right=yy$t2.right, status2=yy$status2)

  attr(toreturn, "init")     <- init
  attr(toreturn, "prior")    <- prior
  attr(toreturn, "mcmc.par") <- mcmc.par
  attr(toreturn, "init2")     <- if (doubly) init2     else list()
  attr(toreturn, "prior2")    <- if (doubly) prior2    else list()
  attr(toreturn, "mcmc.par2") <- if (doubly) mcmc.par2 else list()

  return(toreturn)  
}  
