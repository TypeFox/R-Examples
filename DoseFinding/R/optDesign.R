## optimal designs for model-fitting

## calculate gradient of model and gradient of TD
calcGrads <- function(fmodels, doses, weights,
                      Delta, off, scal, direction,
                      designCrit){
  modgrad <- TDgrad <- nPar <- vector("list", modCount(fmodels, fullMod=TRUE))
  z <- 1
  for(nam in names(fmodels)){
    pars <- fmodels[[nam]]
    if(is.matrix(pars)){
      for(i in 1:nrow(pars)){
        modgrad[[z]] <- t(gradCalc(nam, pars[i,], doses, off=off, scal=scal)*sqrt(weights))
        if(designCrit != "Dopt")
          TDgrad[[z]] <- calcTDgrad(nam, pars[i,], Delta, direction, off, scal)
        nPar[[z]] <- nPars(nam)
        z <- z+1
      }
    } else {
      modgrad[[z]] <- t(gradCalc(nam, pars, doses, off=off, scal=scal)*sqrt(weights))
      if(designCrit != "Dopt")
        TDgrad[[z]] <- calcTDgrad(nam, pars, Delta, direction, off, scal)
      nPar[[z]] <- nPars(nam)
      z <- z+1        
    }
  }
  modgrads <- do.call("c", modgrad)
  TDgrad <- do.call("c", TDgrad)
  nPar <- do.call("c", nPar)

  list(modgrads=modgrads, TDgrad=TDgrad, nPar=nPar)
}


## returns the number of parameters (needed for C call)
nPars <- function(mods){
  builtIn <- c("linlog", "linear", "quadratic", 
               "emax", "exponential", "logistic", 
               "betaMod", "sigEmax")
  ind <- match(mods, builtIn)
  if(any(is.na(ind))){
    stop(mods[which(is.na(ind))], " model not allowed in optDesign")
  }
  c(2,2,3,3,3,4,4,4)[ind]
}

## function which calls different optimizers
callOptim <- function(func, method, nD, control, lowbnd, uppbnd){
  ## actual optimizer
  if(method == "nlminb"){ # nlminb and optim run on transformed values
    res <- nlminb(getStart(nD), objective=func, control = control,
                  lower=rep(0, nD), upper=rep(pi, nD))
  } else if(method == "Nelder-Mead"){
    res <- optim(getStart(nD), fn=func, control = control)
  } else if(method == "solnp"){ # no need for transformed values for solnp
    avail <- requireNamespace("Rsolnp", quietly = TRUE)
    if(!avail)
      stop("Need suggested package Rsolnp for this calculation to use solnp optimizer")
    ## get starting value (need feasible starting value for solnp)
    ## try whether equal allocation is feasible
    eq <- rep(1/nD, nD)
    if(all(eq > lowbnd+0.001) & all(eq < uppbnd-0.001)){
      strt <- eq
    } else {
      slb <- sum(lowbnd)
      sub <- sum(uppbnd)
      gam <- (1-slb)/(sub-slb)
      strt <- lowbnd+gam*(uppbnd-lowbnd)
    }
    eqfun <- function(x, ...){
      sum(x)
    }
    con <- list(trace = 0)
    con[(namc <- names(control))] <- control
    res <- Rsolnp::solnp(strt, fun=func, eqfun=eqfun, eqB=1,
                         control = con, LB = lowbnd, UB = uppbnd)
  } 
  res
}

## transforms from unconstrained values R^k into constrained
## values in S^k = {w|sum_i w_i=1 and w_i >= 0}
transTrig <- function(y, k){
  a <- numeric(k)  
  if(k == 2){
    a[1] <- sin(y[1])^2
  } else {
    a[1:(k-1)] <- sin(y)^2
    a[2:(k-1)] <- a[2:(k-1)]*cumprod(cos(y[1:(k-2)])^2)
  }
  a[k] <- prod(cos(y[1:(k-1)])^2)
  a
}

## identity function
idtrans <- function(y, k){
  y
}

## calculate uniform design but on R^k scale
## (inverse of transTrig at uniform design)
getStart <- function(k){
  y <- numeric(k-1)
  eq <- 1/k
  y[1] <- asin(sqrt(eq))
  for(j in 2:(k-1)){
    y[j] <- asin(sqrt(eq/prod(cos(y[(1:j)-1])^2)))
  }
  y
}

## function called in the optimization (design criterion is
## implemented in C and called "critfunc")
optFunc <- function(x, xvec, pvec, nD, probs, M, n, nold, bvec, designCrit,
                    trans, standInt){
  xtrans <- do.call("trans", list(x, nD))
  res <- .C("critfunc", xvec, pvec, nD, probs, M, xtrans, n,
            nold, double(16), as.double(1e-15), bvec, designCrit, standInt,
            double(1), PACKAGE = "DoseFinding")
  res[[14]]
}

## user visible function calling all others
optDesign <- function(models, probs, doses,
                      designCrit = c("Dopt", "TD", "Dopt&TD", "userCrit"),
                      Delta, standDopt = TRUE, weights,
                      nold = rep(0, length(doses)),  n,
                      control=list(), 
                      optimizer = c("solnp", "Nelder-Mead", "nlminb", "exact"),
                      lowbnd = rep(0, length(doses)), uppbnd = rep(1, length(doses)),
                      userCrit, ...){
  if(!missing(models)){
    if(!inherits(models, "Mods"))
      stop("\"models\" needs to be of class Mods")
    direction <- attr(models, "direction")
    off <- attr(models, "off")
    scal <- attr(models, "scal")
    if(missing(doses))
      doses <- attr(models, "doses")
  } else {
    if(missing(userCrit))
      stop("either \"models\" or \"userCrit\" need to be specified")
    if(missing(doses))
      stop("For userCrit one always needs to specify doses")
  }
  ## check arguments
  designCrit <- match.arg(designCrit)
  optimizer <- match.arg(optimizer)
  if(missing(n)){
    if(optimizer == "exact")
      stop("need to specify sample size via n argument")
    if(any(nold > 0))
      stop("need to specify sample size for next cohort via n argument")
    n <- 1 ## value is arbitrary in this case
  } else {
    if(length(n) > 1)
      stop("n needs to be of length 1")
  }
  if(missing(Delta)){
    if(designCrit %in% c("TD", "Dopt&TD"))
      stop("need to specify target difference \"Delta\"")
  } else {
    if(Delta <= 0)
      stop("\"Delta\" needs to be > 0, if curve decreases use \"direction = decreasing\"")
  }
  if(missing(weights)){
    weights <- rep(1, length(doses))
  } else {
    if(length(weights) != length(doses))
      stop("weights and doses need to be of equal length")
  }
  if(length(lowbnd) != length(doses))
    stop("lowbnd needs to be of same length as doses")
  if(length(uppbnd) != length(doses))
    stop("uppbnd needs to be of same length as doses")
  if(any(lowbnd > 0) | any(uppbnd < 1)){
    if(optimizer != "solnp" & optimizer != "exact")
      stop("only optimizers solnp or exact can handle additional constraints on allocations")
  }
  if(sum(lowbnd) > 1)
    stop("Infeasible lower bound specified (\"sum(lowbnd) > 1\"!)")
  if(sum(uppbnd) < 1)
    stop("Infeasible upper bound specified (\"sum(lowbnd) < 1\"!)")
  if(!is.logical(standDopt))
    stop("standDopt needs to contain a logical value")
  standInt <- as.integer(standDopt) # use standardized or non-stand. D-optimality
  nD <- length(doses)
  if(designCrit == "TD" | designCrit == "Dopt&TD"){ # check whether TD exists in (0,max(dose))
    if(length(unique(direction)) > 1)
      stop("need to provide either \"increasing\" or \"decreasing\" as direction to optDesign, when TD optimal designs should be calculated")
    direction <- unique(direction)
    tdMods <- TD(models, Delta, "continuous", direction)
    tdMods[tdMods > max(doses)] <- NA
    if(any(is.na(tdMods)))
      stop("TD does not exist for ",
           paste(names(tdMods)[is.na(tdMods)], collapse=", " ), " model(s)")
  }
  if(designCrit == "Dopt" | designCrit == "Dopt&TD"){ # check whether Fisher matrix can be singular
    np <- nPars(names(models))
    if(max(np) > length(doses))
      stop("need at least as many dose levels as there are parameters to calculate Dopt design.")
  } 

  ## use transformation for Nelder-Mead and nlminb
  if(is.element(optimizer, c("Nelder-Mead", "nlminb"))){ 
    transform <- transTrig
  } else {
    transform <- idtrans
  }
  
  if(designCrit != "userCrit"){ # prepare criterion function
    ## check arguments
    if(abs(sum(probs)-1) > sqrt(.Machine$double.eps)){
      stop("probs need to sum to 1")
    }
    ## prepare criterion function
    lst <- calcGrads(models, doses, weights,
                     Delta, off, scal, direction, designCrit)
    ## check for invalid values (NA, NaN and +-Inf)
    checkInvalid <- function(x)
      any(is.na(x)|(is.nan(x)|!is.finite(x)))
    grInv <- checkInvalid(lst$modgrads)
    MvInv <- ifelse(designCrit != "Dopt", checkInvalid(lst$TDgrad), FALSE)
    if(grInv | MvInv)
      stop("NA, NaN or +-Inf in gradient or bvec")
    ## prepare arguments before passing to C
    M <- as.integer(length(probs))
    if(M != length(lst$nPar))
      stop("probs of wrong length")
    if(length(lst$modgrads) != length(doses)*sum(lst$nPar))
      stop("Gradient of wrong length.")
    if(length(nold) != nD)
      stop("Either nold or doses of wrong length.")
    nD <- as.integer(nD)
    p <- as.integer(lst$nPar)
    intdesignCrit <- match(designCrit, c("TD", "Dopt", "Dopt&TD"))
    objFunc <- function(par){
      optFunc(par, xvec=as.double(lst$modgrads),
              pvec=as.integer(p), nD=nD, probs=as.double(probs),
              M=M, n=as.double(n), nold = as.double(nold),
              bvec=as.double(lst$TDgrad), trans = transform,
              standInt = standInt,designCrit = as.integer(intdesignCrit))
    }
  } else { # user criterion
    if(missing(userCrit))
      stop("need design criterion in userCrit when specified")
    if(!is.function(userCrit))
      stop("userCrit needs to be a function")
    objFunc <- function(par){
      par2 <- do.call("transform", list(par, nD))
      userCrit((par2*n+nold)/(sum(nold)+n), doses, ...)
    }
  }

  ## perform actual optimization
  if(optimizer != "exact"){ # use callOptim function
    res <- callOptim(objFunc, optimizer, nD, control, lowbnd, uppbnd)
    if(optimizer == "Nelder-Mead" | optimizer == "nlminb"){ # transform results back
      des <- transTrig(res$par, length(doses))
      if(optimizer == "Nelder-Mead"){
        crit <- res$value
      } else {
        crit <- res$objective
      }
    }
    if(optimizer == "solnp"){ # no need to transform back
      des <- res$pars
      crit <- res$values[length(res$values)]
    }
    if(res$convergence){
      message("Message: algorithm indicates no convergence, the 'optimizerResults'
               attribute of the returned object contains more details.")
    }
  } else { # exact criterion (enumeration of all designs)
    ## enumerate possible exact designs
    con <- list(maxvls1 = 1e6, maxvls2 = 1e5, groupSize = 1)
    con[(namc <- names(control))] <- control    
    mat <- getDesMat(n, nD, lowbnd, uppbnd,
                     con$groupSize, con$maxvls1, con$maxvls2)
    designmat <- sweep(mat*n, 2, nold, "+")
    res <- sweep(designmat, 2, n+sum(nold), "/")
    ## evaluate criterion function
    if(designCrit != "userCrit"){
      critv <- calcCrit(res, models, probs, doses,
                        designCrit, Delta, standDopt,
                        weights, nold, n)
    } else {
      critv <- apply(res, 1, objFunc)
    }
    des <- mat[which.min(critv),]
    crit <- min(critv)
  }
  out <- list(crit = crit, design = des, doses = doses, n = n,
              nold = nold, designCrit = designCrit)
  attr(out, "optimizerResults") <- res
  class(out) <- "DRdesign"
  out
}

calcCrit <- function(design, models, probs, doses, 
                     designCrit = c("Dopt", "TD", "Dopt&TD"),
                     Delta, standDopt = TRUE, weights,
                     nold = rep(0, length(doses)), n){
  if(!inherits(models, "Mods"))
    stop("\"models\" needs to be of class Mods")
  off <- attr(models, "off")
  scal <- attr(models, "scal")
  if(missing(doses))
    doses <- attr(models, "doses")  
  ## extract design
  if(inherits(design, "DRdesign"))
    design <- design$design
  if(!is.numeric(design))
    stop("design needs to be numeric")
  if(!is.matrix(design))
    design <- matrix(design, ncol = length(design))
  if(ncol(design) != length(doses))
    stop("design and doses should be of the same length")      
  if(any(abs(rowSums(design)-1) > 0.001))
    stop("design needs to sum to 1")
  if(missing(n)){
    n <- 1 # value arbitrary
  } else {
    if(length(n) > 1)
      stop("n needs to be of length 1")
  }
  if(missing(weights)){
    weights <- rep(1, length(doses))
  } else {
    if(length(weights) != length(doses))
      stop("weights and doses need to be of equal length")
  }
  designCrit <- match.arg(designCrit)
  if(missing(Delta) & substr(designCrit, 1, 3) == "TD")
    stop("need to specify clinical relevance parameter")
  direction <- attr(models, "direction")

  if(designCrit == "TD" | designCrit == "Dopt&TD"){ # check whether TD exists in (0,max(dose))
    if(length(unique(direction)) > 1)
      stop("need to provide either \"increasing\" or \"decreasing\" as direction to optDesign, when TD optimal designs should be calculated")
    direction <- unique(direction)
    tdMods <- TD(models, Delta, "continuous", direction)
    tdMods[tdMods > max(doses)] <- NA
    if(any(is.na(tdMods)))
      stop("TD does not exist for ",
           paste(names(tdMods)[is.na(tdMods)], collapse=", " ), " model(s)")
  }
  if(designCrit == "Dopt" | designCrit == "Dopt&TD"){ # check whether Fisher matrix can be singular
    np <- nPars(names(models))
    if(max(np) > length(doses))
      stop("need more dose levels to calculate Dopt design.")
  }
  if(!is.logical(standDopt))
    stop("standDopt needs to contain a logical value")
  standInt <- as.integer(standDopt)
  lst <- calcGrads(models, doses, weights, Delta, off, scal,
                   direction, designCrit)
  ## check for invalid values (NA, NaN and +-Inf)
  checkInvalid <- function(x)
    any(is.na(x)|(is.nan(x)|!is.finite(x)))
  grInv <- checkInvalid(lst$modgrads)
  MvInv <- ifelse(designCrit != "Dopt", checkInvalid(lst$TDgrad), FALSE)
  if(grInv | MvInv)
    stop("NA, NaN or +-Inf in gradient or bvec")
  ## prepare for input into C
  M <- as.integer(length(probs))
  nD <- as.integer(length(doses))
  if(M != length(lst$nPar))
    stop("Probs of wrong length")
  if(length(lst$modgrads) != length(doses)*sum(lst$nPar))
    stop("Gradient of wrong length.")
  
  if(length(nold) != nD)
    stop("Either nold or doses of wrong length.")
  p <- as.integer(lst$nPar)
  intdesignCrit <- match(designCrit, c("TD", "Dopt", "Dopt&TD"))
  res <- numeric(nrow(design))
  ## check for sufficient number of design points
  iter <- 1:nrow(design)
  design0 <- sweep(design, 2, nold, "+")
  count <- apply(design0, 1, function(x) sum(x > 0.0001))
  ind <- count < max(p[probs > 0])
  if(any(ind)){
    iter <- iter[!ind]
    res[ind] <- NA
    if(all(is.na(res)))
      warning("need at least as many dose levels in the design as parameters in the model")
  }
  for(i in iter){
    res[i] <- optFunc(design[i,], xvec=as.double(lst$modgrads),
                      pvec=as.integer(p), nD=nD, probs=as.double(probs),
                      M=M, n=as.double(n), nold = as.double(nold),
                      bvec=as.double(lst$TDgrad), trans = idtrans,
                      standInt = standInt, designCrit = as.integer(intdesignCrit))
  }
  res
}

## print designs
print.DRdesign <- function(x, digits = 5, eps = 0.001, ...){
  nam <- switch(x$designCrit,
                "TD" = "TD",
                "Dopt" = "D",
                "Dopt&TD" = "TD and D mixture",
                "userCrit" = "userCrit")
  cat("Calculated", nam, "- optimal design:\n")
  ind <- x$design > eps
  vec <- x$design[ind]
  names(vec) <- x$doses[ind]
  print(round(vec, digits = digits))
}

## auxiliary function for efficient rounding
which.is.max <- function (x){
    y <- seq_along(x)[x == max(x)]
    if (length(y) > 1L) 
        sample(y, 1L)
    else y
}

## efficient rounding (see Pukelsheim (1993), Ch. 12)
rndDesign <- function(design, n, eps = 0.0001){

  if(missing(n))
    stop("total sample size \"n\" needs to be specified")
  n <- round(n) # ensure n is an integer (at least numerically)
  if(inherits(design, "DRdesign")){
    design <- design$design
  }
  if(!inherits(design, "numeric"))
    stop("design needs to be a numeric vector.")
  zeroind <- design < eps
  if(any(zeroind)){
    design <- design[!zeroind]/sum(design[!zeroind])
  }
  l <- sum(!zeroind)
  nn <- ceiling((n-0.5*l)*design)
  while(sum(nn)!=n){
    if(sum(nn)<n){
      indmin <- which.is.max(-nn/design)
      nn[indmin] <- nn[indmin]+1
    } else {
      indmax <- which.is.max((nn-1)/design)
      nn[indmax] <- nn[indmax]-1
    }
  }
  if(any(zeroind)){
    out <- numeric(length(design))
    out[zeroind] <- 0
    out[!zeroind] <- nn
    return(out)
  } else {
    nn
  }
}


getCompositions <- function(N, M){
  nC <- choose(N+M-1, M-1)
  lst <- .C("getcomp", comp=integer(nC*M), integer(M-1),
            as.integer(N), as.integer(M-1), as.integer(nC),
            PACKAGE = "DoseFinding")
  matrix(lst$comp, byrow = TRUE, nrow = nC)
}


## calculate all possible compositions of n patients to nDoses groups
## (assuming a certain block-size) upper and lower bounds on the
## allocations can also be specified
getDesMat <- function(n, nDoses, lowbnd = rep(0, nDoses), 
                      uppbnd = rep(1, nDoses), groupSize,
                      maxvls1, maxvls2){
  if(n %% groupSize)
    stop("n needs to be divisible by groupSize")
  nG <- n/groupSize
  combn <- choose(nG+nDoses-1,nDoses-1)
  if(combn > maxvls1)
    stop(combn, " (unrestricted) combinations, increase maxvls1 in control 
         argument if this calculation should be performed")

  desmat <- getCompositions(nG, nDoses)/nG
 
  if(any(lowbnd > 0) | any(uppbnd < 1)){
    comp <- matrix(lowbnd, byrow = TRUE, ncol = nDoses, nrow=nrow(desmat))
    LindMat <- desmat >= comp
    comp <- matrix(uppbnd, byrow=TRUE, ncol = nDoses, nrow=nrow(desmat))
    UindMat <- desmat <= comp
    ind <- rowSums(LindMat*UindMat) == nDoses
    desmat <- desmat[ind,]
    if(nrow(desmat) == 0)
      stop("no design is compatible with bounds specified in lowbnd and uppbnd")
  }
  if(nrow(desmat) > maxvls2)
    stop(nrow(desmat), " combinations, increase maxvls2 in control argument if
         this calculation should be performed")
  desmat
}

## plot method for design objects
plot.DRdesign <- function(x, models, lwdDes = 10, colDes = rgb(0,0,0,0.3), ...){
  if(missing(models))
    stop("need object of class Mods to produce plot")
  plot(models, ...)
  layoutmat <- trellis.currentLayout()
  nc <- ncol(layoutmat)
  nr <- nrow(layoutmat)
  total <- sum(layoutmat > 0)
  z <- 1
  for(i in 1:nc){
    for(j in 1:nr){
      if(z > total)
        break
      trellis.focus("panel", i, j)
      args <- trellis.panelArgs()
      miny <- min(args$y)
      maxy <- max(args$y)
      dy <- maxy-miny
      for(k in 1:length(x$doses)){
        yy <- c(0,(x$design*dy)[k])+miny
        xx <- rep(x$doses[k],2)
        panel.xyplot(xx, yy, type="l", col = colDes, lwd = lwdDes)
      }
      z <- z+1
      trellis.unfocus()
    }
  }
}
