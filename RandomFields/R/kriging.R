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


## source("modelling.R")




FinImputIntern <- function(data, simu, coords, coordnames, data.col, vdim,
                           spConform, fillall=FALSE) {
  n <- length(data) / (vdim * coords$restotal)
  #Print(data, all, tail(all$simu), spConform);
  if (is(data, "RFsp")) {
    if (spConform) {
      data@data[ , ] <- as.vector(simu)
      return(data)
    } else {
      values <- as.matrix(data@data)
      values[is.na(values) | fillall] <- simu
      return(cbind(coordinates(data), values))
    }
  } else { ## not RFsp
    #Print("for testing")
    if (coords$grid) {
      ## to do
      stop("not programmed yet")
    } else {
      ##  coords <- all$x
      colnames(coords$x) <- coordnames
      
      values <- data[, data.col]
      values[is.na(values) | fillall] <- simu
      
      if (!spConform)  return(cbind(coords$x, values))
      
      tmp.all <- conventional2RFspDataFrame(data=values, coords=coords$x,
                                            gridTopology=NULL,
                                            n=n, vdim=vdim,
                                            vdim_close_together=FALSE)
      if (is(tmp.all, "RFspatialPointsDataFrame"))
        try(tmp.all <- as(tmp.all, "RFspatialGridDataFrame"), silent=TRUE)
      if (is(tmp.all, "RFpointsDataFrame"))
        try(tmp.all <- as(tmp.all, "RFgridDataFrame"), silent=TRUE)
    }
    return(tmp.all)
  }
}


FinishImputing <- function(data, simu, Z, spConform, fillall) {
  ## to do: grid

  if (is.list(data)) {
    for (i in 1:length(data))
      data[[i]] <- FinImputIntern(data=data[[i]], simu=simu[[i]],
                                  coords=Z$coord[[i]], coordnames=Z$coordnames,
                                  data.col=Z$data.col, vdim=Z$vdim,
                                  spConform = spConform, fillall=fillall)
    return(data)
  }

  return(FinImputIntern(data=data[[1]], simu=simu, coords=Z$coord[[1]],
                        coordnames=Z$coordnames, data.col=Z$data.col, 
                        vdim=Z$vdim, spConform=spConform, fillall=fillall))

}



ExpandGrid <- function(x) {
  #### ACHTUNG! ZWINGENDE REIHENFOLGE
  if (x$grid) { # 0
    x$x <-
      as.matrix(do.call(expand.grid,
                        lapply(apply(cbind(x$x, x$T), 2,
                                     function(x) list(seq(x[1],by=x[2],length.out=x[3]))), function(x) x[[1]])))
  } else if (x$Zeit) {
    dim.x <- if (is.vector(x$x)) c(length(x$x), 1) else dim(x$x)
    x$x <- cbind(matrix(rep(t(x$x), times=x$T[3]),
                            ncol=dim.x[2], byrow=FALSE),
                     rep(seq(x$T[1], by=x$T[2],
                             length.out=x$T[3]), each=dim.x[1]))
  }
  if (length(x$y) > 0) stop("no expansion within a kernel definition")
#  x$y <- double(0) #1 
  x$T <- double(0) #2
  x$grid <- FALSE  #3
#  x$spatialdim <- ncol(x$x) #4
  x$Zeit <- FALSE           #5
#  x$dist.given <- FALSE      #6
  x$restotal <- nrow(x$x)   #7
  x$l <- x$restotal         #8
  return(x)
}



rfPrepareData <- function(model, x, y=NULL, z=NULL, T=NULL,
                          distances=NULL, dim, grid,
                          data, given=NULL, 
                          RFopt, reg, err.model = NULL,
                          ...) {  

  if (!missing(distances) && length(distances)>0)
    stop("option distances not programmed yet.")

#  Print(model=model, data=data, given=given, T, ...)
  missing.x <- missing(x) || length(x) == 0
  imputing <- missing.x && length(distances) == 0
   krige <- model <- PrepareModel2(model)
  if (!is.null(err.model)) {
    linpart <- RFlinearpart(model=err.model, new$x, set=1)
    if (length(linpart$X) > 0 || any(linpart$Y != 0))
      stop("a trend is not allowed for the error model.")
    krige <- list(ZF_SYMBOLS_PLUS, PrepareModel2(err.model, ...), krige)
  }

  if (length(given) == 0) {
    ## so either within the data or the same the x-values
    Z <- StandardizeData(model=model, data=data, RFopt=RFopt, ...)
    if (Z$matrix.indep.of.x.assumed) {
      
      if (missing.x) stop("coordinates cannot be detected")
      Z <- StandardizeData(model=model, x=x, y=y, z=z, T=T, RFopt=RFopt,
                              distances=distances, dim=dim, grid=grid,
                              data=data, ...)
    }
  } else   {
    Z <- StandardizeData(model=model, data=data, x=given, RFopt=RFopt, ...)
  }
  
  #str(Z)
 
  if (length(Z$data) != 1) stop("exactly one data set must be given.")
  dimdata <- base::dim(Z$data[[1]])
  Z$data[[1]] <- as.double(Z$data[[1]])
  repet <- Z$repetitions
  new.dimdata <- c(prod(dimdata) / repet, repet)
  base::dim(Z$data[[1]]) <- new.dimdata
  data.na <- is.na(Z$data[[1]])
  data.na.var <- rowSums(data.na)
  base::dim(Z$data[[1]]) <- dimdata
  base::dim(data.na.var) <- c(length(data.na.var) / Z$vdim , Z$vdim)
  data.na.loc <- rowSums(data.na.var > 0) > 0
  any.data.na <- any(data.na.loc)
  split <- any(data.na.var > 0 & data.na.var != repet)

  
  if (any.data.na && Z$coord[[1]]$dist.given)
    stop("missing values not programmed yet for given distancs")

  if (imputing) {
    if (Z$vdim > 1) stop("imputing does not work in the multivariate case")
    if (repet == 1)  {
      if (RFopt$krige$fillall || !any.data.na) {
        data.na <- rep(TRUE, length(data.na)) ## nur
    ## um Daten im Endergebnis einzutragen
        new <- Z$coord[[1]]
      } else {
        new <- ExpandGrid(Z$coord[[1]])
        new$x <- new$x[data.na.loc, , drop=FALSE]
      }
    } else new <- NULL
  } else {
    new <- CheckXT(x, y, z, T, grid=grid, distances=distances, dim=dim)
    if (Z$tsdim != new$spatialdim + new$Zeit)
      stop("coodinate dimension of locations with and without data, ",
           "respectively, do not match.")
  }
   
  if (any.data.na) {
    if (RFopt$general$na_rm_lines && (!imputing || repet==1)) {
      Z$data[[1]] <-  Z$data[[1]][!data.na.loc, , drop=FALSE]
      Z$coord[[1]] <- ExpandGrid(Z$coord[[1]])
      Z$coord[[1]]$x <- Z$coord[[1]]$x[!data.na.loc,  , drop=FALSE]
    } else if (split) {
      data <- list()
      for (i in 1:repet) {
        dim(Z$data[[1]]) <- c(length(Z$data)[[1]] / (Z$vdim*repet), Z$vdim, repet)
        data[[i]] <-  Z$data[[1]][ , , i, drop=FALSE]
      }
      Z$data <- data
    }
  }

   

   return(list(Z=Z, X=new, krige = krige, model=model,
              imputing=imputing, data.na = if (imputing) data.na))
}




RFinterpolate <- function(model, x, y=NULL, z=NULL, T=NULL, grid=NULL,
                          distances, dim, data, given=NULL,
                          err.model=NULL, ignore.trend=FALSE, ...) {
  if (!missing(distances) && length(distances) > 0) stop("'distances' not programmed yet.")
    
  opt <- list(...)
  i <- pmatch(names(opt), c("MARGIN"))
  opt <- opt[is.na(i)]
  
  RFoptOld <- do.call("internal.rfoptions", c(opt, RELAX=isFormulaModel(model)))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  boxcox <- .Call("get_boxcox")


  ## eingabe wird anstonsten auch als vdim_close erwartet --
  ## dies ist nocht nicht programmiert! Ausgabe ist schon programmiert
  ## CondSimu ist auch noch nicht programmiert
  if (RFopt$general$vdim_close_together)
    stop("'vdim_close_together' must be FALSE")

  
  reg <- MODEL_KRIGE
  return.variance <- RFopt$krige$return_variance

  all <- rfPrepareData(model=model, x=x, y=y, z=z, T=T,
                       distances=distances, dim=dim, grid=grid,
                       data=data, given=given, RFopt=RFopt,
                       reg=reg, err.model = err.model,
                       ...)
  #Print(all);
  
  imputing <- all$imputing
  tsdim <- as.integer(all$Z$tsdim)
  repet <- as.integer(all$Z$repetitions)
  vdim <- all$Z$vdim
  if (!imputing) {
    coordnames.incl.T <-
      c(if (!is.null(all$Z$coordnames)) all$Z$coordnames else
        paste(ZF_GENERAL_COORD_NAME[1], 1:all$Z$spatialdim, sep=""),
        if (all$Z$Zeit) ZF_GENERAL_COORD_NAME[2] else NULL)
    if (all$X$grid) {
      coords <- list(x=NULL, T=NULL)
      xgr <- cbind(all$X$x, all$X$T)
      colnames(xgr) <- coordnames.incl.T
      gridTopology <- sp::GridTopology(xgr[1, ], xgr[2, ], xgr[3, ])
      ## bis 3.0.70 hier eine alternative
    } else {
      coords <- list(x=all$X$x, T=all$X$T)
      ## wenn bei gegeben unklar was zu tun ist. Ansonsten
      if (length(coords$T) == 0)  colnames(coords$x) <- coordnames.incl.T
      gridTopology <- NULL
    }
  }
  nx <- all$X$restotal
  dimension <-
    if (all$X$grid) c(if (length(all$X$x) > 0) all$X$x[3, ],
                      if (length(all$X$T) > 0) all$X$T[3]) else nx # to do:grid
  newdim <- c(dimension, if (vdim>1) vdim, if (repet>1) repet)

  if (imputing && return.variance) {
    return.variance <- FALSE
    warning("with imputing, currently the variance cannot be returned")
  }


  if (length(all$Z$data) > 1) {
    Res <- array(dim=c(nx, vdim, repet))
    for (i in 1:length(all$Z$data)) {
       Res[ , , i] <-
        RFinterpolate(model=model, x=x, y=y, z=z, T=T, grid=grid,
                      distances=distances, dim=dim,
                      data = all$Z$data[[i]], given = all$Z$coord,
                      err.model=err.model, ...,
                      spConform = FALSE, return.variance=FALSE)      
    }
    dim(Res) <- c(nx * vdim, repet)
  } else { ## length(all$Z$data) == 1   
    exact <- RFopt$general$exact
    maxn <- RFopt$krige$locmaxn
    ngiven <- as.integer(all$Z$coord[[1]]$restotal) ## number of given points
    split <- RFopt$krige$locsplitn[1]
    split <- ngiven > maxn || (!is.na(exact) && !exact && ngiven > split)

    data <- RFboxcox(all$Z$data[[1]])
    .Call("set_boxcox", c(Inf, 0))
 
    if (imputing) {
      Res <- data
    } else {
      Res <- matrix(nrow=nx, ncol=repet * vdim)
    }
    
    if (split) {
      ## to do:
      all$X <- ExpandGrid(all$X) ## to  do
      all$Z$coord[[1]] <- ExpandGrid(all$Z$coord[[1]]) ## to  do
      
      ## neighbourhood kriging !
      if (!is.na(exact) && exact)
        stop("number of conditioning locations too large for an exact result.")
      if (ngiven > maxn && is.na(exact) &&
          RFopt$general$printlevel>=PL_IMPORTANT)
        message("performing neighbourhood kriging")

      stop("neighbourhood kriging currently not programmed")

      ## calculate the boxes for the locations where we will interpolate
      idx <- GetNeighbourhoods(Z=Z,
                               X=all$X, ## given locations; to do: grid
                               splitfactor=RFopt$krige$locsplitfactor,
                               maxn=RFopt$krige$locmaxn,
                               split_vec = RFopt$krige$locsplitn,
                               )
      totalparts <- length(idx[[2]])
      
      if (totalparts > 1) RFoptions(general.pch="")
      pr <- totalparts > 1 && RFopt$general$pch != "" &&RFopt$general$pch != " "

      for (p in 1:totalparts) {
        stopifnot((Nx <- as.integer(length(idx[[3]][[p]]))) > 0)
        if (pr && p %% 5==0) cat(RFopt$general$pch)
        givenidx <- unlist(idx[[1]][idx[[2]][[p]]])
       if (ignore.trend) 
         initRFlikelihood(all$krige, Reg=reg, grid=FALSE,
                         x=all$Z$coord[[1]]$x[givenidx,  , drop=FALSE],
                         data=data[givenidx, , drop=FALSE],
                         ignore.trend = ignore.trend)
        else
          RFlikelihood(all$krige, Reg=reg, grid=FALSE,
                       x=all$Z$coord[[1]]$x[givenidx,  , drop=FALSE],
                       data=data[givenidx, , drop=FALSE],
                       likelihood_register = reg)
        res <- predictGauss(Reg=reg, model=all$model,
                            x = all$X[idx[[3]][[p]], ], grid = FALSE,
                            kriging_variance=FALSE) 
        
        if (imputing) {
          ## TO DO : idx[[3]] passt nicht, da sowohl fuer Daten
          ##         als auch coordinaten verwendet wird. Bei repet > 1
          ##         ist da ein Problem -- ueberpruefen ob repet=1
          
          where <- all$data.na[idx[[3]][[p]]]  ## to do:grid
          isNA <- is.na(Res[where, ])
          Res[where, ][isNA] <- res[isNA]        
        } else {
          Res[idx[[3]][[p]], ] <- res
        }
      } ## for p in totalparts
      if (pr) cat("\n")
    } else { ## not split
      if (ignore.trend) 
        initRFlikelihood(all$krige, Reg=reg, x=all$Z$coord, data=data,
                         ignore.trend = ignore.trend)
      else RFlikelihood(all$krige, x=all$Z$coord, data=data,
                        likelihood_register = reg)
      Res <- predictGauss(Reg=reg, model=all$model, x=all$X,
                          kriging_variance=FALSE)

    #  Print(Res)
      
     }
  }## !is.list(Z)
  Z <- all$Z ## achtung! oben kann sich noch all$Z aendern!
  X <- all$X

  #  Print(newdim, Res, vdim, repet, dimension, X$grid, nx, X)
#  print(Res)
 #  Print(newdim)

  
  
  if (length(newdim)>1)  base::dim(Res) <- newdim else Res <- as.vector(Res)

  if (return.variance && length(newdim <- c(dimension, if (vdim>1) vdim)) > 1)
    base::dim(sigma2) <- newdim  
  if (!is.null(Z$varnames)) attributes(Res)$varnames <- Z$varnames

  spConform <- RFopt$general$spConform

  Res <- RFboxcox(data=Res, boxcox = boxcox, inverse=TRUE)

  if (!spConform && !imputing) {
    if (vdim > 1 && RFopt$general$vdim_close_together) {
      Resperm <- c(length(dimension)+1, 1:length(dimension),
                   if(repet>1) length(dimension)+2)      
      Res <- aperm(Res, perm=Resperm)
      
      if (return.variance)
        sigma2 <- aperm(sigma2, perm=Resperm[1:(length(dimension)+1)])
    }
    if (return.variance) Res <- list(estim = Res, var = sigma2)
    #class(Res) <- "RandomFieldsReturn"    
    return(Res)
  }


  if (imputing) {
    ## Print(data, 1, spConform, Res) 
    Res <- FinishImputing(data=data, simu=Res, Z=Z, spConform=spConform,
                          fillall = RFopt$krige$fillall) ## to do : grid
    if (return.variance){
      var <- FinishImputing(data=data, simu=sigma2, Z=Z, spConform=spConform,
                            fillall = RFopt$krige$fillall)# to do : grid
      if (spConform) {
        names(var@data) <- paste("var.", names(var@data), sep="")
        Res@.RFparams$has.variance <- TRUE
        Res <-  cbind(Res, var)
      } else Res <- list(Res, var)
    }
 #   print(Res)
    return(Res)
  } else {
  
    Res <- conventional2RFspDataFrame(Res, coords=coords$x,
                                      gridTopology=gridTopology,
                                      n=repet, vdim=vdim, T = coords$T,
                                      vdim_close_together =
                                      RFopt$general$vdim_close_together)
  
    if (return.variance){
      var <- conventional2RFspDataFrame(sigma2, coords=coords$x,
                                        gridTopology=gridTopology,
                                        n=1, vdim=vdim, T = coords$T,
                                        vdim_close_together =
                                        RFopt$general$vdim_close_together)
      names(var@data) <- paste("var.", names(var@data), sep="")
      Res@.RFparams$has.variance <- TRUE
      Res <-  cbind(Res, var)
    }
  }
  
#  Res@.RFparams$krige.method <-
#    c("Simple Kriging", "Ordinary Kriging", "Kriging the Mean",
#      "Universal Kriging", "Intrinsic Kriging")[krige.meth.nr]
  

  ## Res@.RFparams$var <- sigma2 ## sehr unelegant.
  ## * plot(Res) sollte zwei Bilder zeigen
  ## * var(Res) sollte sigma2 zurueckliefern
  ## * summary(Res) auch summary der varianz, falls vorhanden
  ## * summary(Res) auch die Kriging methode

  if (is.raster(x)) {
    Res <- raster::raster(Res)
    raster::projection(Res) <- raster::projection(x)
  }
   
  return(Res)
}




rfCondGauss <- function(model, x, y=NULL, z=NULL, T=NULL, grid, n=1,
                        data,   # first coordinates, then data
                        given=NULL, ## alternative for coordinates of data
                        err.model=NULL, ...) { # ... wegen der Variablen  
  dots <- list(...)
  if ("spConform" %in% names(dots)) dots$spConform <- NULL

  RFoptOld <- internal.rfoptions(..., RELAX=isFormulaModel(model)) 
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  boxcox <- .Call("get_boxcox")
  cond.reg <- RFopt$registers$register
 
  all <- rfPrepareData(model=model, x=x, y=y, z=z, T=T, grid=grid,
                       data=data, given=given, RFopt=RFopt,
                       reg=MODEL_KRIGE, err.model = err.model, ...)
  Z <- all$Z
  X <- all$X
  simu.grid <- X$grid
  tsdim <- Z$tsdim
  vdim <- Z$vdim

  data <- RFboxcox(Z$data[[1]])
  .Call("set_boxcox", c(Inf, 0))

  if (all$Z$repetitions != 1)
     stop("conditional simulation allows only for a single data set")
  
  txt <- "kriging in space time dimensions>3 where not all the point ly on a grid is not possible yet"
  ## if 4 dimensional then the last coordinates should ly on a grid

  ## now check whether and if so, which of the given points belong to the
  ## points where conditional simulation takes place
  coord <- ExpandGrid(Z$coord[[1]])
  simu <- NULL
  if (simu.grid) {
    xgr <- cbind(X$x, X$T)
    ind <- 1 + (t(coord$x) - xgr[1, ]) / xgr[2, ] 
    index <- round(ind)
    outside.grid <-
      apply((abs(ind-index)>RFopt$general$gridtolerance) | (index<1) |
            (index > 1 + xgr[3, ]), 2, any)

    if (any(outside.grid)) {
      ## at least some data points are not on the grid:
      ## simulate as if there is no grid
      simu.grid <- FALSE
      ll <- NULL ##  otherwise check will give a warning
      l <- ncol(xgr)

      if (l>3) stop(txt)
      xx <- if (l==1) ## dim x locations
             matrix(seq(from=xgr[1], by=xgr[2], len=xgr[3]),
                        nrow=1)
            else eval(parse(text=paste("t(expand.grid(",
                            paste("seq(from=xgr[1,", 1:l, 
                                  "], by=xgr[2,", 1:l,
                                  "], len=xgr[3,", 1:l, "])", collapse=","),
                         "))")))  
      ll <- eval(parse(text=paste("c(",
                   paste("length(seq(from=xgr[1,", 1:l, 
	                 "], by=xgr[2,", 1:l, 
		         "], len=xgr[3,", 1:l, "]))",
                         collapse=","),
                   ")")))

      new.index <- rep(0,ncol(index))
      ## data points that are on the grid, must be registered,
      ## so that they can be used as conditioning points of the grid
      if (!all(outside.grid)) {
        new.index[!outside.grid] <- 1 +
          colSums((index[, !outside.grid, drop=FALSE]-1) *
                  cumprod(c(1, ll[-length(ll)])))
      }
      index <- new.index
      new.index <- NULL
    } else {  
      ## data points are all lying on the grid
     
      simu <- do.call(RFsimulate, args=c(list(model=all$krige,
                                      x=X$x, # y=y, z=z,
                                      T=X$T,
                                      grid=X$grid,
                                      n=n, 
                                      register=cond.reg,
                                      seed = NA),
                                      dots, list(spConform=FALSE)))
      ## for all the other cases of simulation see, below
      index <- t(index)
      if (is.vector(simu)) dim(simu) <- c(length(simu), 1)
      else {
        d <- dim(simu)
        last <- d + 1 - (1 : ((vdim > 1) + (n > 1)))
       if (!is.matrix(simu)) dim(simu) <- c(prod(d[-last]), prod(d[last]))
      }
   }
  } else { ## not simu.grid
    xx <- t(X$x)  ## dim x locations
   
    ## the next step can be pretty time consuming!!!
    ## to.do: programme it in C
    ##
    ## identification of the points that are given twice, as points to
    ## be simulated and as data points (up to a tolerance distance !)
    ## this is important in case of nugget effect, since otherwise
    ## the user will be surprised not to get the value of the data at
    ## that point
    one2ncol.xx <- 1:ncol(xx)
    index <- apply(coord$x, 1, function(u){
      i <- one2ncol.xx[colSums(abs(xx - u)) < RFopt$general$gridtolerance]
      if (length(i)==0) return(0)
      if (length(i)==1) return(i)
      return(NA)
    })
  }

  
  if (!simu.grid) {
    ## otherwise the simulation has already been performed (see above)
    tol <- RFopt$general$gridtolerance * nrow(xx)
    if (any(is.na(index)))
      stop("identification of the given data points is not unique - `tol' too large?")
    if (any(notfound <- (index==0))) {
      index[notfound] <- (ncol(xx) + 1) : (ncol(xx) + sum(notfound))
    }
    
    xx <- rbind(t(xx), coord$x[notfound, , drop=FALSE])
    simu <- do.call(RFsimulate,
                    args=c(list(model=all$krige, x=xx, grid=FALSE, n=n,
                        register = cond.reg, seed = NA), dots,
                        spConform=FALSE, examples_reduced = FALSE))
    if (is.vector(simu)) dim(simu) <- c(length(simu), 1)
    #Print(simu, X$restotal, index)
    #print(simu)
#    Print(c(list(model=all$krige, x=xx, grid=FALSE, n=n,
 #                       register = cond.reg, seed = NA), dots,
  #          spConform=FALSE, examples_reduced = FALSE));
   # print(tail(simu));
    
    rm("xx")
  }

  if (is.null(simu)) stop("random field simulation failed")
  
  simu.given <- simu[index, ]
  simu <- as.vector(simu[1:X$restotal, ]) # as.vector is necessary !! Otherwis
##                                       is not recognized as a vector


#  Print(simu, simu.given, simu.grid, index, as.vector(data),  simu.given, X)

 
  ## to do: als Naeherung bei UK, OK:
  ## kriging(data, method="A") + simu - kriging(simu, method="O") !
  stopifnot(length(X$y)==0, length(X$z)==0)
  simu <- simu + RFinterpolate(x=X, model=model,
                               err.model = err.model,
                               register=MODEL_KRIGE,
                               given = coord,
                               data = as.vector(data) - simu.given,
                               spConform=FALSE, ignore.trend = TRUE)
  simu <- RFboxcox(data=simu, boxcox = boxcox, inverse=TRUE)
  
  if (all$imputing) {
    return(FinishImputing(data=Z$data[[1]], simu=simu, Z=Z,
                          spConform=RFopt$general$spConform,
                          fillall = RFopt$krige$fillall))
  }
  
  return(simu)
  
}
