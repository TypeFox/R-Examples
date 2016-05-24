## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  

# source("rf.R")
# see getNset.R for the variable .methods



RFboxcox <- function(data, boxcox, vdim=1, inverse=FALSE, ignore.na=FALSE) {
  if (missing(boxcox)) boxcox <- .Call("get_boxcox")
  if (any(is.na(boxcox)) && !ignore.na)
    stop("non-finte values in Box-Cox transformation")
  if (!all(is.finite(boxcox))) return(data)
  if (is.list(data)) {
    for (i in 1:length(data))
      data[[i]] <- RFboxcox(data[[i]], boxcox=boxcox, vdim=vdim,
                            inverse=inverse, ignore.na=ignore.na)
    return(data)
  }
  Data <- data + 0
  .Call("BoxCox_trafo", as.double(boxcox), as.double(Data), as.integer(vdim),
        as.logical(inverse));
  return(Data)
}


RFlinearpart <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                         data, distances, dim, set=0, ...) {
  Reg <- MODEL_USER  
  relax <- isFormulaModel(model)
  RFoptOld <- internal.rfoptions(..., RELAX=relax)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  model <- list("linearpart", PrepareModel2(model))

  rfInit(model=model, x=x, y=y, z=z, T=T, grid=grid,
         distances=distances, dim=dim, reg = Reg, dosimulate=FALSE)

  .Call("get_linearpart", Reg, as.integer(set))
}




predictGauss <- function(Reg,
                         Reg.predict= if (missing(model)) Reg else
                         RFopt$register$predict_register,
                         model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                         data, distances, dim, kriging_variance, ...) {
  relax <- isFormulaModel(model)
  RFoptOld <- internal.rfoptions(..., RELAX=relax)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  if (missing(model) || is.null(model)) {
    if (Reg != Reg.predict) stop("'Reg.predict' should not be given")
    model <- list("predict", register=Reg)
  } else {
    if (Reg == Reg.predict) stop("'Reg.predict' must be different from 'Reg'")
     model <- list("predict", PrepareModel2(model), register=Reg)
  }

 # str(model); xxxx
  
  rfInit(model=model, x=x, y=y, z=z, T=T, grid=grid,
         distances=distances, dim=dim, reg = Reg.predict, dosimulate=FALSE)

  .Call("EvaluateModel", double(0), as.integer(Reg.predict),
        PACKAGE="RandomFields")
}

setvector <- function(model, preceding, len, factor) {
  if (model[[1]] == ZF_SYMBOLS_PLUS) {
    return(c(ZF_SYMBOLS_PLUS,
             lapply(model[-1], setvector, preceding=preceding, len=len,
                    factor=factor)))
  }
  if (found <- model[[1]] == ZF_SYMBOLS_C) {
    if (length(model) != len + 1) stop("bug")
    model <- c(model[1], rep(0, preceding), model[-1])
    if (!missing(factor)) model <- list(ZF_SYMBOLS_MULT, model, factor)
  } else if (model[[1]] == ZF_SYMBOLS_MULT) {
    for (i in 2:length(model)) {
      if (found <- model[[i]][[1]]==ZF_SYMBOLS_C) {
        if (length(model[[i]]) != len + 1) stop("bug")
        model[[i]] <- c(ZF_SYMBOLS_C, rep(0, preceding), model[[i]][-1])
        if (!missing(factor)) model[[length(model) + 1]] <- factor
        break
      }
    }
  }
  if (!found) {
    bind <- c(list(ZF_SYMBOLS_C), rep(0, preceding), rep(1, len))
    if (model[[1]] == ZF_SYMBOLS_MULT)
      model <- c(ZF_SYMBOLS_MULT, list(bind), model[-1])
    else model <- list(ZF_SYMBOLS_MULT, bind, model)
    if (!missing(factor)) model[[length(model) + 1]] <- factor
  }
  return(model)
}

splittingC <- function(model, preceding, factor) {
  const <- sapply(model[-1],
                  function(m) (is.numeric(m) && !is.na(m)) ||
                  (m[[1]] == ZF_SYMBOLS_CONST && !is.na(m[[2]]))
                  )
  if (all(const)) {
    model <- c(model[1], if (preceding > 0) rep(0, preceding), model[-1])
    return(list(ZF_SYMBOLS_MULT, model, if (!missing(factor)) list(factor)))
  }
  for (i in 2:length(model)) {
    vdim <- preceding + (if (i==2) 0 else GetDimension(model[[i-1]]))
  #  Print(i, length(model), vdim)
    m  <- ReplaceC(model[[i]])
    L <- GetDimension(m)
  #  Print(m, L, preceding, vdim)
    model[[i]] <- setvector(m, preceding = vdim, len = L, factor=factor)
   # Print(model[[i]])
  }
  model[[1]] <- ZF_SYMBOLS_PLUS
  names(model) <- NULL
  L <- GetDimension(model[[length(model)]])
 # Print(L)
  model <- SetDimension(model, L)
}


GetDimension <- function(model){
  if (model[[1]] == ZF_SYMBOLS_PLUS) {
    return(max(sapply(model[-1], GetDimension)))
  }
  if (model[[1]] == ZF_SYMBOLS_C) return(length(model) - 1)
  else if (model[[1]] == ZF_SYMBOLS_MULT) {
    L <- sapply(model, function(m) if (m[[1]]==ZF_SYMBOLS_C) length(m)-1 else 1)
    return(max(L))
  }
  return(1)
}


SetDimension <- function(model, L){
  if (model[[1]] == ZF_SYMBOLS_PLUS) {
    return(c(ZF_SYMBOLS_PLUS, lapply(model[-1], SetDimension, L=L)))
  }
  if (model[[1]] == ZF_SYMBOLS_C) {
    if (length(model) <= L) for (i in length(model) : L) model[[i+1]] <- 0
    names(model) <- c("", letters[1:L])
  } else if (model[[1]] == ZF_SYMBOLS_MULT) {
    for (i in 2:length(model)) {
#      Print(model[[i]], L, model[[i]][[1]]==ZF_SYMBOLS_C)
      if (model[[i]][[1]]==ZF_SYMBOLS_C) {
        if (length(model[[i]]) <= L)
          for (j in length(model[[i]]) : L) model[[i]][[j+1]] <- 0
        names(model[[i]]) <- c("", letters[1:L])
      }
    }
  }
  return(model)
}



SplitC <- function(model, factor) {
#  Print("SplitC", model)
  model <- splittingC(model, 0, factor)
  L <- GetDimension(model)
  return(SetDimension(model, L))
}


AnyIsNA <- function(model) {
#  Print(model)
  if (is.list(model)) {
    for (i in 1:length(model)) if (AnyIsNA(model[[i]])) return(TRUE)
    return(FALSE)
  } else return((is.numeric(model) || is.logical(model)) && any(is.na(model)))
}


ReplaceC <- function(model) {
 # Print("Replace", model)
  if (model[[1]] == ZF_SYMBOLS_PLUS) {
    for (i in 2:length(model)) model[[i]] <- ReplaceC(model[[i]])
  } else if (model[[1]] == ZF_SYMBOLS_MULT) {
    if (length(model) <= 2) {
      stop("here, products must have at least 2 factors")
    }
    cs <- sapply(model, function(m) m[[1]] == ZF_SYMBOLS_C)
    s <- sum(cs)
    if (s > 0) {
      if (s > 1)
        stop("multiplication with '", ZF_SYMBOLS_C, "' may happen only once")
      cs <- which(cs)
      C <- model[[cs]]      
      if (!AnyIsNA(model[-cs])) {
#        Print(model[-cs], model[-cs])
        return(SplitC(C, factor=model[-cs]))
      }
    }      
  } else if (model[[1]] == ZF_SYMBOLS_C) {
    return(SplitC(model))
  }
  return(model)
}

initRFlikelihood <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                             data, distances, dim, likelihood,
                             estimate_variance = NA,
                              Reg, ignore.trend=FALSE, ...) {

  if (!missing(likelihood)) ## composite likelihood etc
    stop("argument 'likelihood' is a future feature, not programmed yet")


  model <- PrepareModel2(model, ...)
  model <- ReplaceC(model); ## multivariates c() aufdroeseln
  
  model <- list("loglikelihood", model, data = data,
                estimate_variance=estimate_variance,
                betas_separate = FALSE, ignore_trend=ignore.trend)

  return (rfInit(model=model, x=x, y=y, z=z, T=T, grid=grid,
                 distances=distances, dim=dim, reg = Reg, dosimulate=FALSE))

}


RFlikelihood <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                         data, distances, dim, likelihood,
                         estimate_variance = NA,
                         ...) {
  relax <- isFormulaModel(model)
  RFoptOld <-
    if (missing(likelihood)) internal.rfoptions(..., RELAX=relax)
    else internal.rfoptions(likelihood=likelihood, ..., RELAX=relax)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  Reg <- RFopt$register$likelihood_register
    
  initRFlikelihood(model=model, x=x, y = y, z = z, T=T, grid=grid,
                   data=data,
                   distances=distances, dim=dim, likelihood=likelihood,
                   estimate_variance = estimate_variance,
                   Reg = Reg, ...)  

  likeli <- .Call("EvaluateModel", double(0),  Reg, PACKAGE="RandomFields")
  info <- .Call("get_likeliinfo", Reg)
  globalvariance <- info$estimate_variance
  where <- 1 + globalvariance
  param <- likeli[-1:-where]
  if (length(param) > 0) names(param) <- info$betanames
  
  return(list(loglikelihood = likeli[1], likelihood = exp(likeli[1]),
              global.variance = if (globalvariance) likeli[where] else NULL,
              parameters = param
              )
         )
}


rfdistr <- function(model, x, q, p, n, dim=1, ...) {
  ## note: * if x is a matrix, the points are expected to be given row-wise
  ##       * if not anisotropic Covariance expects distances as values of x
  ##
  ## here, in contrast to Covariance, nonstatCovMatrix needs only x

  RFoptOld <- internal.rfoptions(..., RELAX=isFormulaModel(model))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  if (!missing(n) && n>10 && RFopt$internal$examples_reduced) {
    message("number of simulations reduced")
    n <- 10
  }
    
  model<- list("Distr", PrepareModel2(model), dim=as.integer(dim));
  if (!missing(x)) {
    model$x <- if (is.matrix(x)) t(x) else x
  }
  if (!missing(q)) {
    model$q <- if (is.matrix(q)) t(q) else q
  }
  if (!missing(p)) {
    model$p <- if (is.matrix(p)) t(p) else p
  }
  if (!missing(n)) {
    model$n <- n
  }

#  Print(model)
  
  rfInit(model=model, x=matrix(0, ncol=dim, nrow=1),
         y=NULL, z=NULL, T=NULL, grid=FALSE, reg = MODEL_USER,
         dosimulate=FALSE, old.seed=RFoptOld[[1]]$general$seed)

  res <-  .Call("EvaluateModel", double(0), as.integer(MODEL_USER),
                PACKAGE="RandomFields")

  if (RFoptOld[[2]]$general$returncall) attr(res, "call") <-
    as.character(deparse(match.call(call=sys.call(sys.parent()))))
  attr(res, "coord_system") <- c(orig=RFoptOld[[2]]$coords$coord_system,
                                 model=RFoptOld[[2]]$coords$new_coord_system)
  return(res)
}

RFdistr <- function(model, x, q, p, n, dim=1, ...) {
   rfdistr(model=model, x=x, q=q, p=p, n=n, dim=dim, ...)
}
RFddistr <- function(model, x, dim=1,...) {
  if (hasArg("q") || hasArg("p") || hasArg("n")) stop("unknown argument(s)");
  rfdistr(model=model, x=x, dim=dim,...)
}
RFpdistr <- function(model, q, dim=1,...) {
  if (hasArg("x") || hasArg("p") || hasArg("n")) stop("unknown argument(s)");
  rfdistr(model=model, q=q, dim=dim,...)
}
RFqdistr <- function(model, p, dim=1,...){
  if (hasArg("x") || hasArg("q") || hasArg("n")) stop("unknown argument(s)");
  rfdistr(model=model, p=p, dim=dim,...)
}
RFrdistr <- function(model, n, dim=1,...) {
  if (hasArg("x") || hasArg("q") || hasArg("p")) stop("unknown argument(s)");
  rfdistr(model=model, n=n, dim=dim,...)
}



rfeval <- function(model, x, y = NULL, z = NULL, T=NULL, grid=NULL,
                  distances, dim, ..., 
                  ##                  dim = ifelse(is.matrix(x), ncol(x), 1),
                  fctcall=c("Cov", "CovMatrix", "Variogram",
                    "Pseudovariogram", "Fctn")) {

  ## note: * if x is a matrix, the points are expected to be given row-wise
  ##       * if not anisotropic Covariance expects distances as values of x
  ##
  ## here, in contrast to Covariance, nonstatCovMatrix needs only x

   # print("rfeval in rf.R")  #Berreth
    
  RFoptOld <- internal.rfoptions(xyz_notation=2*(length(y)!=0 && !is.matrix(y)),
                                 ..., RELAX=isFormulaModel(model))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  
  fctcall <- match.arg(fctcall)

  if (fctcall != "CovMatrix" && !missing(distances) && !is.null(distances)) {
    if(missing(dim) || length(dim) != 1) {
      warning("'dim' not given or not of length 1, hence set to 1");
      dim <- 1
    }
    if (length(y) != 0 || length(z) != 0 || length(T) != 0 ||
        (!missing(grid) && length(grid) != 0))
      stop("if distances are given 'y', 'z', 'T', 'grid' may not be given")
    x <- (if (is.matrix(distances)) distances else
          cbind(distances, matrix(0, nrow=length(distances), ncol = dim - 1)))
    distances <- NULL
    dim <- NULL
    grid <- FALSE
  }
  
  p <- list(fctcall, PrepareModel2(model));

  # Print(p, x, distances)
 
  rfInit(model=p, x=x, y=y, z=z, T=T, grid=grid,
         distances=distances, dim=dim, reg = MODEL_USER, dosimulate=FALSE,
         old.seed=RFoptOld[[1]]$general$seed)

  res <- .Call("EvaluateModel", double(0), as.integer(MODEL_USER),
               PACKAGE="RandomFields")
  
  if (RFoptOld[[2]]$general$returncall) attr(res, "call") <-
    as.character(deparse(match.call(call=sys.call(sys.parent()))))
  attr(res, "coord_system") <- .Call("GetCoordSystem", as.integer(MODEL_USER),
              RFoptOld[[2]]$coords$coord_system,
              RFoptOld[[2]]$coords$new_coord_system)
   return(res)
}


RFcov <- function(model, x, y = NULL, z = NULL, T=NULL, grid,
                        distances, dim, ...) {
  rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid,
                distances=distances, dim=dim, ..., fctcall="Cov")
}


RFcovmatrix <- function(model, x, y = NULL, z = NULL, T=NULL, grid,
                        distances, dim, ...) {
  rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid,
                distances=distances, dim=dim, ..., fctcall="CovMatrix")
}

RFvariogram <- function (model, x, y=NULL,  z = NULL, T=NULL, grid,
                        distances, dim,...){
  rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid,
                distances=distances, dim=dim, ..., fctcall="Variogram")
}

RFpseudovariogram <- function(model, x, y=NULL,  z = NULL, T=NULL, grid,
                        distances, dim,...){
   rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid,
                distances=distances, dim=dim, ..., fctcall="Pseudovariogram")
}

RFfctn <- function(model, x, y=NULL,  z = NULL, T=NULL, grid,
                   distances, dim,...) {
   # print("RFfctn in rf.R") #Berreth
  rfeval(model=model, x=x, y=y, z=z, T=T, grid=grid,
                distances=distances, dim=dim, ..., fctcall="Fctn")
}

RFcalc <- function(model) {
  if (is.numeric(model)) return(model)
  RFfctn(model, 0,
         coord_system="cartesian", new_coord_system="keep", spConform = FALSE)
}


######################################################################
######################################################################


rfDoSimulate <- function(n = 1, reg, spConform) {
  stopifnot(length(n) == 1, n>0, is.finite(n))
  RFopt <- RFoptions()
  if (missing(spConform)) spConform <- RFopt$general$spConform

  if (RFopt$gauss$paired && (n %% 2 != 0))
    stop("if paired, then n must be an even number")

  info <- RFgetModelInfo(RFopt$registers$register, level=3)

  len <- info$loc$len
  vdim <- info$vdim
  total <- info$loc$totpts
  if (is.null(total) || total <= 0)
    stop("register ", RFopt$registers$register, " does not look initialized")

  error <- integer(1)

  result <- .Call("EvaluateModel", as.double(n), as.integer(reg), #userdefined,
                  PACKAGE="RandomFields")
  
  if (!spConform) return(result)
  
  prep <- prepare4RFspDataFrame(model=NULL, info=info, RFopt=RFopt)
  attributes(result)$varnames <- prep$names$varnames
  
  res2 <- conventional2RFspDataFrame(result,
                                     coords=prep$coords,
                                     gridTopology=prep$gridTopology,
                                     n=n,
                                     vdim=prep$vdim,
                                     T=info$loc$T,
                                     vdim_close_together
                                     =RFopt$general$vdim_close_together)
  return(res2)
}



rfInit <- function(model, x, y = NULL, z = NULL, T=NULL, grid=FALSE,
                   distances, dim, reg, dosimulate=TRUE, old.seed=NULL) {
 
                                        ##print("rfInit in rf.R") #Berreth
  stopifnot(xor(missing(x), #|| length(x)==0,
                missing(distances) || length(distances)==0))

  RFopt <- RFoptions() 
  if (!is.na(RFopt$general$seed)) {
    allequal <- all.equal(old.seed, RFopt$general$seed)
    allequal <- is.logical(allequal) && allequal
    if (dosimulate && RFopt$general$printlevel >= PL_IMPORTANT &&
        (is.null(old.seed) || (!is.na(old.seed) && allequal)
         )
        ) {
       message("NOTE: simulation is performed with fixed random seed ",
               RFopt$general$seed,
               ".\nSet 'RFoptions(seed=NA)' to make the seed arbitrary.")
     }
    set.seed(RFopt$general$seed)
  }
  ##  if (missing(x) || length(x) == 0) stop("'x' not given")

  new <- C_CheckXT(x, y, z, T, grid=grid, distances=distances, dim=dim,
                 y.ok=!dosimulate)
 
  return(.Call("Init", as.integer(reg), model, new, NAOK=TRUE, # ok
               PACKAGE="RandomFields"))
}




    
RFsimulate <- function (model, x, y = NULL, z = NULL, T = NULL, grid=NULL,
                        distances, dim, data, given = NULL, err.model,
                        n = 1,  ...) {  
  mc <- as.character(deparse(match.call()))

### preparations #########################################################  
  if (!missing(distances) && length(distances)  > 0)
    RFoptOld <- internal.rfoptions(xyz_notation=length(y)!=0,
                                   expected_number_simu=n, ..., 
                                   general.spConform = FALSE,
                                   RELAX=isFormulaModel(model))
  else 
    RFoptOld <- internal.rfoptions(xyz_notation=length(y)!=0,
                                   expected_number_simu=n, ...,
                                   RELAX=isFormulaModel(model))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  if (n>2 && RFopt$internal$examples_reduced) {
    message("number of simulations reduced")
    n <- 2
  }
  
  cond.simu <- !missing(data) && !is.null(data) 
  reg <- RFopt$registers$register

  ### simulate from stored register ########################################
  mcall <- as.list(match.call(expand.dots=FALSE))
  if (length(mcall)==1 ||
      length(mcall)==2 && !is.null(mcall$n) ||
      length(mcall)==3 && !is.null(mcall$n) && "..." %in% names(mcall)) {
    if (cond.simu) {
      stop("repeated performance of conditional simulation not programmed yet")
    } else {
      # userdefined <- GetParameterModelUser(PrepareModel2(model, ...))
      res <- rfDoSimulate(n=n, reg=reg, spConform=RFopt$general$spConform
                          #userdefined=userdefined
                          )
      if (RFopt$general$returncall) attr(res, "call") <- mc
      attr(res, "coord_system") <- .Call("GetCoordSystem", reg,
                                     RFopt$coords$coord_system,
                                     RFopt$coords$new_coord_system)
      return(res)
    }
  }
  
  ### preparations #########################################################
  stopifnot(!missing(model) && !is.null(model))

  model.orig <- model
  model <- PrepareModel2(model, ...)
  err.model <-  if (missing(err.model)) NULL else PrepareModel2(err.model, ...)


  ### conditional simulation ###############################################
  if (cond.simu) {
    #Print(RFoptions()$general)
    if (isSpObj(data)) data <- sp2RF(data)
    stopifnot(missing(distances) || is.null(distances))
    res <- switch(GetProcessType(model),
                  RPgauss = 
                  rfCondGauss(model=model.orig, x=x, y=y, z=z, T=T,
                              grid=grid, n=n, data=data, given=given,
                              err.model=err.model,
                              ## next line to make sure that this part
                              ## matches with predictGauss
                              predict_register = MODEL_PREDICT,
                              ...),
                  stop(GetProcessType(model),
                       ": conditional simulation of the process not programmed yet")
                  )
  } else { ## unconditional simulation ####
    if(!is.null(err.model))
      warning("error model is unused in unconditional simulation")

    rfInit(model=list("Simulate",
               setseed=eval(parse(text="quote(set.seed(seed=seed))")),
               env=.GlobalEnv, model), x=x, y=y, z=z, T=T,
           grid=grid, distances=distances, dim=dim, reg=reg,
           old.seed=RFoptOld[[1]]$general$seed)
    if (n < 1) return(NULL)
    res <- rfDoSimulate(n=n, reg=reg, spConform=FALSE)
   } # end of uncond simu


  ## output: RFsp   #################################
  if ((!cond.simu || (!missing(x) && length(x) != 0)) ## not imputing
      && RFopt$general$spConform) {
    info <- RFgetModelInfo(if (cond.simu) MODEL_PREDICT else reg, level=3)
    if (length(res) > 1e7) {
      message("Too big data set (more than 1e7 entries) to allow for 'spConform=TRUE'. So the data are returned as if 'spConform=FALSE'")
      return(res)
    }
    
    prep <- prepare4RFspDataFrame(model.orig, info, x, y, z, T,
                                  grid, data, RFopt)
    attributes(res)$varnames <- prep$names$varnames

#    Print(res, prep, n, info)
    
    res <- conventional2RFspDataFrame(data=res, coords=prep$coords,
                                      gridTopology=prep$gridTopology,
                                      n=n,
                                      vdim=prep$vdim,
                                      T=info$loc$T,
                                      vdim_close_together=
                                      RFopt$general$vdim_close_together)
    if (is.raster(x)) {
      res <- raster::raster(res)
      raster::projection(res) <- raster::projection(x)
    }
  }
  
  if (RFopt$general$returncall) attr(res, "call") <- mc
  attr(res, "coord_system") <- .Call("GetCoordSystem",
                                     as.integer(reg),
                                     RFopt$coords$coord_system,
                                     RFopt$coords$new_coord_system)
  return(res)
}


