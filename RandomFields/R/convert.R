
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



PrepareModel2 <- function(model, ..., x=NULL) {
  if (missing(model) || is.null(model)) stop("'model' must be given.")
  method <- "ml"
  if (class(model) == "RF_fit") model <- model[[method]]$model
  else if (class(model) == "RFfit") model <- model[method]

  m <- parseModel(model, ..., x=x)
  
  if (notplus <- !(m[[1]] %in% ZF_PLUS)) m <- list(ZF_SYMBOLS_PLUS, m)
     
  for (i in 2:length(m)) {
    if ((m[[i]][[1]] %in% ZF_MIXED) && length(m[[i]]$X)==1 &&
        is.numeric(m[[i]]$X) && m[[i]]$X==1 && !is.null(m[[i]]$b)) {        
      m[[i]] <- list(ZF_TREND[2], mean=m[[i]]$b)
      if (RFoptions()$general$printlevel > PL_IMPORTANT)
        message(paste("The '1' in the mixed model definition has been replaced by '", ZF_TREND[1], "(mean=", m[[i]]$mean, ")'.", sep=""))
      }
  }

  if (notplus) m <- m[[2]]
  class(m) <- "RM_model"
  return(m)

#  if (class(model) != "formula") {
#    if (is.list(model)) return(model)
#    else stop("model of unknown form -- maybe you have used an obsolete definition. See ?RMmodel for the model definition")
#  }
#  return(listmodel)
}


PrepareModel <-  function(model, param, trend=NULL, 
                          nugget.remove=TRUE, method=NULL) {
  ## any of the users model definition (standard, nested, list) for the
  ## covariance function is transformed into a standard format, used
  ## especially in the c programs
  ##
  ## overwrites in some situation the simulation method for nugget.
  ## allows trend to be NA (or any other non finite value  -- is not checked!)
  ## trend has not been implemented yet!

  if (is(model, ZF_MODEL))
    stop("models of class ZF_MODEL cannot be combined with obsolete RandomFields functions")

  if (!is.null(method)) stop("to give method in PrepareModel is obsolete")
  
  if (!is.null(trend))      
     if (!is.numeric(trend) || length(trend)!=1)
        stop("in the obsolete setting, only constant mean can used")
 
  if (is.list(model) && is.character(model[[1]]) &&
      (is.null(names(model)) || names(model)[[1]]=="")) {
     if (!missing(param) && !is.null(param))
        stop("param cannot be given in the extended definition")
     if (is.null(trend)) return(model)
     trend <- list(ZF_TREND[2], mean=trend)
     if (model[[1]] %in% ZF_PLUS) return(c(model, list(trend)))
     else return(list(ZF_SYMBOLS_PLUS, model, trend))
  }
    
  printlevel <- RFoptions()$general$printlevel
  STOP <- function(txt) {
    if (printlevel>=PL_ERRORS) {
      cat("model: ")
      if (!missing.model) Print(model) else cat(" missing.\n") #
      cat("param: ")
      if (!missing.param) Print(param) else cat(" missing.\n") #
      cat("trend: ")
      Print(trend) # 
     }
    stop("(in PrepareModel) ", txt, call.=FALSE)
  }
   
  transform <- function(model) {
    if (!is.list(model)) {
      STOP("some elements of the model definition are not lists")
    }
    m <- list(DOLLAR[1], var=model$v)
    lm <- length(model) - 3 # var, scale/aniso, name
    if (!is.null(model$a)) m$aniso <- model$a else m$scale <- model$scale
##    model <- c(model, if (!is.null(model$a))
##               list(aniso=model$a) else list(scale=model$s)) ## ???

    
    if (!is.na(p <- pmatch("meth", names(model), duplicates.ok=TRUE))) {
      if (printlevel>=PL_ERRORS)  Print(p, model) #
      stop("method cannot be given with the model anymore. It must be given as a parameter to the function. See 'RFoptions' and 'RFsimulate'")
     }
  
    if (!is.null(model$me))
      stop("'mean' seems to be given within the inner model definitions"); 
    if (!is.character(model$m)) {
       stop("'model' was not given extacly once each odd number of list entries or additional unused list elements are given.")
    }
    m1 <- list(model$m)
    if (!is.null(model$k)) {
      lm <- lm - 1
      if (length(model$k) != 0)
        for (i in 1:length(model$k)) {
          eval(parse(text=paste("m1$k", i, " <- model$k[", i, "]", sep="")))
      }
    }
    if (lm != 0) {
      if (printlevel>=PL_ERRORS) Print(lm, model) #
      stop("some parameters do not fit")
    }
    m <- c(m, list(m1))

    return(m)
    
  } # end transform

  op.list <- c(ZF_SYMBOLS_PLUS, ZF_SYMBOLS_MULT)  ## if others use complex list definition !
  missing.model <- missing(model)
  missing.param <- missing(param) || is.null(param)

  if (missing.param && is.null(model$param)) { ## full model
    if (RFoptions()$internal$warn_oldstyle)
      warning("the sequential list format is depreciated.")
    if (missing.model || (length(model)==0)) model <- list()
    else if (!is.list(model))
      STOP("if param is missing, model must be a list of lists (or a list in the extended notation)")
    if (is.null(trend) + is.null(model$mean) + is.null(model$trend)<2)
      STOP("trend/mean is given twice")
    if (!is.null(model$mean)) trend <- model$mean else
    if (!is.null(model$trend)) trend <- model$trend else trend <- NULL

    model$trend <- model$mean <- NULL
    ## the definition might be given at a deeper level as element
    ## $model of the list:
    if (is.list(model$model)) {
      if (!is.list(model$model[[1]]))
        STOP("if param is missing, the model$model must be a list of lists")
      model <- model$model
    }
    if (length(model)==0) { ## deterministic      
      return(if (is.null(trend)) NULL else list(ZF_TREND[2], mean=trend))
    }
    if (length(model) %% 2 !=1) STOP("list for model definition should be odd")
    if (length(model)==1)
      return(if (is.null(trend) ||
                 is.numeric(trend) && length(trend)==1 && !is.na(trend)&&trend==0)
             transform(model[[1]]) 
             else list(ZF_SYMBOLS_PLUS, transform(model[[1]]),
                       list(ZF_TREND[2], mean=trend)));

    op <- pmatch(c(model[seq(2, length(model), 2)], recursive=TRUE),
                 op.list, duplicates.ok=TRUE) - 1
    if (!all(is.finite(op))) STOP("operators are not all allowed; see the extended list definition for extensions")
    model <- model[seq(1, length(model), 2)]

    plus <- which(op==0)
    if (length(plus) == 0) {
      m <- list("*", lapply(model, transform))
    } else {
      plus <- c(0, plus, length(op)+1)
      m <- list(ZF_SYMBOLS_PLUS)
      for (i in 1:(length(plus) - 1)) {
        m[[i+1]] <-
          if (plus[i] + 1 == plus[i+1]) transform(model[[plus[i] + 1]])
          else list(ZF_SYMBOLS_MULT,
                    lapply(model[(plus[i] + 1) : plus[i+1]], transform))
      }
    }
   model <- m
  } else { ## standard definition or nested model
    if (missing.param) { ## a simple list of the model and the
      ##                    parameters is also possible
      if (is.null(param <- model$p)) STOP("is.null(model$param)")
      stopifnot(is.null(trend) || is.null(model$trend))
      if (is.null(trend)) trend <- model$trend
      if (!is.null(model$mean)) {
        if (!is.null(trend)) STOP("mean and trend given twice")
        trend <- model$mean
      }
      model <- model$model
    }
    stopifnot(is.character(model), length(model)==1)
    if (is.matrix(param)) { ## nested
      if (nrow(param) == 1)
        return(PrepareModel(model=model, param=c(param[1], 0, param[-1]),
                            trend=trend))
      name <- model
      model <- list(ZF_SYMBOLS_PLUS)#, method=method)
      for (i in 1:nrow(param)) {
        model <- c(model,
                   if (is.na(param[i, 2]) || param[i, 2] != 0)
                   list(list(DOLLAR[1], var=param[i, 1], scale=param[i, 2],
                             if (ncol(param) >2) list(name, k=param[i,-1:-2])
                             else list(name)))
                   else list(list(DOLLAR[1], var=param[i,1],
                                  list(ZF_NUGGET[2]))))
      }
    } else if (is.vector(param)) {  ## standard, simple way
      ## falls trend gegeben, dann ist param um 1 Komponente gekuerzt
      if (is.null(trend)) {
        trend <- param[1]
        param <- param[-1]
      } else message("It is assumed that no mean is given so that the first component of param is the variance")
      if (model == ZF_NUGGET[2]) {
        model <- transform(list(model=model, var=sum(param[1:2]), scale=1)) 
      } else {
        if  (length(param) > 3)
          model <- transform(list(model=model, var=param[1], scale=param[3],
                                  k=param[-1:-3]))
        else 
          model <- transform(list(model=model, var=param[1], scale=param[3]))
        if (is.na(param[2]) || param[2] != 0 || !nugget.remove) {# nugget
          model <- list(ZF_SYMBOLS_PLUS,
                        model,
                        transform(list(model=ZF_NUGGET[2], var=param[2], scale=1)))
        }
        ## if (!is.null(method)) model <- c(model, method=method) ## doppelt
      }
    } else stop("unknown format")  # end nested/standard definition
  }

  return(if (is.null(trend) ||
             is.numeric(trend) && length(trend)==1 &&  !is.na(trend) &&trend==0)
         return(model)
         else if (model[[1]] %in% ZF_PLUS)
                 c(model, list(list(ZF_TREND[2], mean=trend)))
         else list(ZF_SYMBOLS_PLUS, model, list(ZF_TREND[2], mean=trend)))
}


seq2grid <- function(x, name, grid, warn_ambiguous, gridtolerance) {
  xx <- matrix(nrow=3, ncol=length(x))
  step0 <- rep(FALSE, length(x))
  gridnotgiven <- missing(grid) || length(grid) == 0
  
  for (i in 1:length(x)) {
    if (length(x[[i]]) == 1) {
      xx[,i] <- c(x[[i]], 0, 1)
      next
    }
    step <- diff(x[[i]])
    if (step[1] == 0.0) {
      
      ok <- step0[i] <- all(step == 0.0)      
    } else {
      ok <- max(abs(step / step[1] - 1.0)) <= gridtolerance
    }

    if (!ok) {
      if (gridnotgiven) return(FALSE)
      if (!TRUE)
        Print(i, x[[i]][1:min(100, length(x[[i]]))], #
              step[1:min(100,length(step))],
              range(diff(step[1:min(100,length(step))])))
      stop("Different grid distances detected, but the grid must ",
           "have equal distances in each direction -- if sure that ",
           "it is a grid, increase the value of 'gridtolerance' which equals ",
           gridtolerance,".\n")
    }

    xx[,i] <- c(x[[i]][1], step[1], if (step0[i]) 1 else length(x[[i]]))
  }

  if (FALSE && gridnotgiven && warn_ambiguous && length(x) > 1) {
    RFoptions(internal.warn_ambiguous = FALSE)
    message("Ambiguous interpretation of coordinates. Better give 'grid=TRUE' explicitly. (This message appears only once per session.)")
  }

  if (any(step0)) {
    if (all(step0)) {
      if (gridnotgiven) return(FALSE)
      else stop("Within a grid, the coordinates must be distinguishable")
    } else {
      if (gridnotgiven && warn_ambiguous) {
        RFoptions(internal.warn_ambiguous = FALSE)
        warning("Interpretation as degenerated grid. Better give 'grid' explicitely. (This warning appears only once per session.)")
      }
    }
  }

  return(xx)
}

CheckXT <- function(x, y=NULL, z=NULL, T=NULL, grid, distances=NULL,
                    dim=NULL, # == spatialdim!
                    length.data,
                    y.ok = FALSE, 
                    printlevel = RFoptions()$general$printlevel){

  ## do not pass anything on "..." ! --- only used for internal calls
  ## when lists are re-passed

  ## converts the given coordinates into standard formats
  ## (one for arbitrarily given locations and one for grid points)
  #print("CheckXT in convert.R")#Berreth
  
  if (!missing(x)) {
    if (is(x, "CheckXT")) return(x)  
    if (is.list(x))  { 
      if (!is.list(x[[1]])) return(do.call("CheckXT", x))
      L <- list()
      for (i in 1:length(x)) {
        L[[i]] <-
          if (is(x[[i]], "CheckXT")) x[[i]] else do.call("CheckXT", x[[i]])
      }
      if (length(x) > 1) {    
        if (!all(diff(sapply(L, function(x) x$Zeit)) == 0) ||
          !all(diff(sapply(L, function(x) x$spatialdim)) == 0))
          stop("all sets must have the same dimension")
        if (!all(diff(sapply(L, function(x) x$dist.given)) == 0))
        stop("either all the sets must be based on distances or none")
      }
      class(L) <- "CheckXT"
      return(L)
    }
  }

  
  RFopt <- RFoptions()
  curunits <- RFopt$coords$coordunits
  newunits <-  RFopt$coords$new_coordunits
  coord_system <-  RFopt$coords$coord_system
  new_coord_system <-  RFopt$coords$new_coord_system
  ex.red <- RFopt$internal$examples_reduced
  
  if (!missing(distances) && !is.null(distances)) { ## length==0 OK!
    stopifnot(is.matrix(distances) || (!missing(dim) && !is.null(dim)),
              (missing(grid) || length(grid) == 0),
              missing(x) || is.null(x),
              length(y)==0,
              length(z)==0,
              length(T)==0)
    
    if (coord_system != new_coord_system && new_coord_system != "keep")      
      stop("coordinate systems differ")
    
    if (is.list(distances)) {
      L <- list()
      for (i in 1:length(distances))
        L[[i]] <- do.call("CheckXT", list(distances=distances[[i]], dim=dim))
       class(L) <- "CheckXT"
      return(L)
    }
        
    if (class(distances) == "dist") {
      x <- as.vector(distances)
      len <- length(distances)
    } else if (is.matrix(distances) || is.vector(distances)) {
      if (is.matrix(distances)) {        
        len <- nrow(distances)
        if (is.null(dim)) dim = ncol(distances)
        else if (dim != ncol(distances))
          stop("matrix of distances does not fit the given dimension")
      } else {
        len <- length(distances)
        if (is.null(dim))
          stop("dim is not given although 'distances' are used")
      }
      x <- distances
    } else {
      stop("'distances' not of required format.")
    }

    if (ex.red && len > ex.red^2 / 2) {
      LEN <- as.integer(ex.red)
      len <- as.integer(LEN * (LEN - 1) / 2)
      x <- if (is.matrix(x)) x[1:len ,] else x[1:len]
    } else {
      LEN <- as.integer(1e-9 + 0.5 * (1 + sqrt(1 + 8 * len)))
      if (LEN * (LEN-1) / 2 != len) LEN <- NaN
    }

    ## keep exactly the sequence up to 'distances'
    if (storage.mode(x) != "double") storage.mode(x) <- "double"
    L <- list(x = as.matrix(x), #0
              y = double(0),   #1
                T= double(0),  #2
              grid = FALSE, #3
              spatialdim=as.integer(dim),#4
              Zeit=FALSE, #5
              dist.given = TRUE, #6
              restotal = LEN, ## number of points
              l = LEN, ## ?? physical length??
              coordunits = curunits,
              new_coordunits = newunits
              )
    class(L) <- "CheckXT"
    return(L)
  }

  stopifnot(!missing(x))
  if (is(x, "RFsp") || isSpObj(x)) {
    return(CheckXT(x=coordinates(x), y=y, z=z, T=T, grid=grid,
                   distances=distances, dim=dim, length.data=length.data,
                   y.ok=y.ok, printlevel=printlevel))
  }    
  if (is.raster(x)) x <- as(x, 'GridTopology')
 
  if ((missing(grid) || length(grid) == 0) && !missing(length.data)) {
    new <-  try(CheckXT(x=x, y=y, z=z, T=T, grid=TRUE, distances=distances,
                        dim=if (!missing(dim)) dim,
                        length.data = length.data, y.ok =y.ok,
                        printlevel = printlevel
                        ), silent=TRUE)
    if (grid <- (class(new) != "try-error")) {
      ratio <- length.data / new$restotal

      if (grid <- ratio == as.integer(ratio)) {
        if (printlevel>=PL_IMPORTANT && new$spatialdim > 1)
          message("Grid detected. If it is not a grid, set grid=FALSE.\n")
      }
    }
    return(if (grid) new else {
      CheckXT(x, y, z, T, grid=FALSE, distances,
              if (!missing(distances) && length(distances) > 0) dim=1,
              length.data = length.data,
              printlevel = printlevel) }
           )
  } # if (missing(grid) && !missing(length.data))


  gridtriple <- FALSE

  if (is.GridTopology <- is(x, "GridTopology")){
    x <- rbind(x@cellcentre.offset,
               x@cellsize,
               x@cells.dim)
    if ((missing(grid) || length(grid) == 0)) grid <- TRUE else stopifnot(grid)
    gridtriple <- TRUE
  }
  ##else {
  ##  is.GridTopology <- FALSE
  ##}

  
  if (is.data.frame(x)) {
    if (ncol(x)==1) x <- as.vector(x) else x <- as.matrix(x)
  }
  
  stopifnot(length(x) != 0)
#  stopifnot(all(unlist(lapply(as.list(x), FUN=function(li) is.numeric(li))))) ## wann benoetigt???

  stopifnot(is.numeric(x))# um RFsimulte(model, data) statt data=data abzufangen
  
  
#  stopifnot(all(is.finite(x)), all(is.finite(y)), all(is.finite(z))) ; s.u. unlist
 
     
  if (is.matrix(x)) {
    if (!is.numeric(x)) stop("x is not numeric.")
    if (length(z)!=0) stop("If x is a matrix, then z may not be given")
    if (length(y)!=0) {
      if (!y.ok) stop("If x is a matrix, then y may not be given")
      if (length(T)!=0)
        stop("If x is a matrix and y is given, then T may not be given")
      if (!is.matrix(y) || ncol(y) != ncol(x) ||
          nrow(x)==3 && nrow(y)!=3 && ((missing(grid) || length(grid) == 0) ||
                                grid))
        stop("y does not match x (it must be a matrix)")
    }

    if (coord_system == COORD_SYS_NAMES[coord_auto + 1] && ncol(x) >= 2
        && ncol(x) <= 3 && !is.null(n <- dimnames(x)[[2]])) {
      if (any(idx <- earth_coordinate_names(n))) {
        if (length(idx) == 2 && !all(idx == 1:2))
          stop("earth coordinates not in order longitude/latitude")
        cur <- curunits[1]
        newunits <- RFopt$coords$new_coordunits
        curunits <- RFopt$coords$coordunits
        curunits[1:2] <- ZF_EARTHCOORD_NAMES[1:2]
        if (newunits[1] == "") newunits[1] <-  UNITS_NAMES[units_km + 1]
        newunits[2:3] <- newunits[1]                
        if (RFopt$internal$warn_coordinates)
          message("\n\nNOTE: current units are ",
                  if (cur=="") "not given and" else paste("'", cur, "', but"),
                  " earth coordinates detected:\n",
                  "earth coordinates will be transformed into units of '",
                  newunits[1],
                  "'.\nIn particular, the values of all scale parameters of ",
                  "any model defined\nin R^3 (currently all models!) are ",
                  "understood in units of '", newunits[1],
                  "'.\nChange options 'coord_system' and/or 'units' if ",
                  "necessary.\n(This message appears only once per session.)\n")
        coord_system <- COORD_SYS_NAMES[earth + 1]
        RFoptions(coords.coord_system = coord_system,
                  coords.coordunits = curunits,
                  coords.new_coordunits = newunits,
                  internal.warn_coordinates=FALSE)

      } else {
         RFoptions(coords.coord_system =  COORD_SYS_NAMES[cartesian + 1])
      }
    }    
   
    spatialdim <- ncol(x)
    len <- nrow(x)
    if (spatialdim==1 && len != 3 && (missing(grid) || length(grid) == 0)) {
      if (length(x) <= 2) grid <- TRUE
      else {
        dx <- diff(x)
        grid <- max(abs(diff(dx))) < dx[1] * RFopt$general$gridtolerance
      }
    } # else {

    if ((missing(grid) || length(grid) == 0) &&
        any(apply(x, 2, function(z) (length(z) <= 2) || max(abs(diff(diff(z))))
                  > RFopt$general$gridtolerance))) {
      grid <- FALSE
    }

    if ((missing(grid) || length(grid) == 0) || !is.logical(grid)) {
      grid <- TRUE
      if (spatialdim > 1 && RFopt$internal$warn_ambiguous) {
        RFoptions(internal.warn_ambiguous = FALSE)
        warning("Ambiguous interpretation of the coordinates. Better give the logical parameter 'grid=TRUE' explicitely. (This warning appears only once per session.)")
      }
    }

    if (grid && !is.GridTopology) {
      if (gridtriple <- len==3) {
        if (printlevel >= PL_SUBIMPORTANT && RFopt$internal$warn_oldstyle) {
          message("x was interpreted as a gridtriple; the new gridtriple notation is:\n  1st row of x is interpreted as starting values of sequences,\n  2nd row as step,\n 3rd row as number of points (i.e. length),\n  in each of the ", ncol(x), " directions.")
        } 
      } else len <- rep(len, times=spatialdim)   # Alex 8.10.2011
    }

    if (grid && !gridtriple) {
      ## list with columns as list elements -- easier way to
      ## do it??
      x <- lapply(apply(x, 2, list), function(r) r[[1]])
      if (length(y) != 0) y <- lapply(apply(y, 2, list), function(r) r[[1]])
    }
  } else { ## x, y, z given separately
    if (length(y)==0 && length(z)!=0) stop("y is not given, but z")
    xyzT <- list(x=if (!missing(x)) x, y=y, z=z, T=T)
    for (i in 1:4) {
      if (!is.null(xyzT[[i]]) && !is.numeric(xyzT[[i]])) {
        if (printlevel>PL_IMPORTANT) 
          message(names(xyzT)[i],
                  " not being numeric it is converted to numeric")
        assign(names(xyzT)[i], as.numeric(xyzT[[i]]))
      }
    }
    remove(xyzT)
    spatialdim <- 1 + (length(y)!=0) + (length(z)!=0)
    if (spatialdim==1 && ((missing(grid) || length(grid) == 0) || !grid)) {
      ## ueberschreibt Einstellung des Nutzers im Falle d=1
      if (length(x) <= 2) newgrid <- TRUE
      else {
        dx <- diff(x)
        newgrid <- max(abs(diff(dx))) < dx[1] * RFopt$general$gridtolerance
      }
      if ((missing(grid) || length(grid) == 0)) grid <- newgrid
      else if (xor(newgrid, grid) && RFopt$internal$warn_on_grid) {
        RFoptions(internal.warn_on_grid = FALSE)
        message("coordinates", if (grid) " do not",
                " seem to be on a grid, but grid = ", grid)
      }
    }
    len <- c(length(x), length(y), length(z))[1:spatialdim]
    
    if (!(missing(grid) || length(grid) == 0) && !grid) { ## sicher nicht grid, ansonsten ausprobieren
      if (any(diff(len) != 0)) stop("some of x, y, z differ in length")
      x <- cbind(x, y, z)
      ## make a matrix out of the list
      len <- len[1]
    } else {
      if ((missing(grid) || length(grid) == 0) && any(len != len[1]))
        grid <- TRUE
      x <- list(x, y, z)[1:spatialdim]
    }
    y <- z <- NULL ## wichtig dass y = NULL ist, da unten die Abfrage
  }  ## end of x, y, z given separately 
  
  if (!all(is.finite(unlist(x)))) {
    stop("coordinates are not all finite")
  }


  if ((missing(grid) || length(grid) == 0) || grid) {
    if (gridtriple) {
      if (len != 3)
        stop("In case of simulating a grid with option gridtriple, exactly 3 numbers are needed for each direction")
      lr <- x[3,] # apply(x, 2, function(r) length(seq(r[1], r[2], r[3])))
      ##x[2,] <- x[1,] + (lr - 0.999) * x[3,] ## since own algorithm recalculates
      ##                               the sequence, this makes sure that
      ##                               I will certainly get the result of seq
      ##                               altough numerical errors may occurs
      restotal <- prod(x[3, ])
      if (length(y)!=0 && !all(y[3,] == x[3,]))
        stop("the grids of x and y do not match ")        
    } else {     
      xx <- seq2grid(x, "x",  grid,
                     RFopt$internal$warn_ambiguous, RFopt$general$gridtolerance)
     if (length(y)!=0) {
        yy <- seq2grid(y, "y", grid,
                       RFopt$internal$warn_ambiguous,
                       RFopt$general$gridtolerance)
        if (xor(is.logical(xx), is.logical(yy)) ||
            (!is.logical(xx) && !all(yy[3,] == xx[3,])))
          stop("the grids for x and y do not match")      
      }
      if (missing(grid) || length(grid) == 0) grid <- !is.logical(xx)       
      if (grid) {
        x <- xx
        if (length(y) != 0) y <- yy
        restotal <- prod(len)
        len <- 3
      } else {
        x <- sapply(x, function(z) z)
        if (length(y) != 0) y <- sapply(y, function(z) z)
      }
    }
    if (grid && any(x[3, ] <= 0))
      stop(paste("step must be postive. Got as steps",
                 paste(x[3,], collapse=",")))
    ##if (len == 1) stop("Use grid=FALSE if only a single point is simulated")
  }
 
  if (!grid) {
    restotal <- nrow(x)
    if (length(y)==0) {
      if (restotal < 200 && any(as.double(dist(x)) == 0)) {
        d <- as.matrix(dist(x))
        diag(d) <- 1
        idx <-  which(as.matrix(d) ==0)
        if (printlevel>PL_ERRORS)
          Print(x, dim(d), idx , cbind( 1 + ((idx-1)%% nrow(d)), #
                                       1 + as.integer((idx - 1)  / nrow(d))) )
        warning("locations are not distinguishable")
      }
      ## fuer hoehere Werte con total ist ueberpruefung nicht mehr praktikabel
    }
  }

  if (coord_system == "earth") {
    # if (ncol(x) > 4) stop("earth coordinates have maximal 3 components")
    opt <- RFoptions()$coords ## muss nochmals neu sein
    global.units <- opt$new_coordunits[1]
    if (global.units[1] == "") global.units <- "km"
   
    Raumdim <- ncol(x) #if (grid) ncol(x) else
    new_is_cartesian <- new_coord_system %in% CARTESIAN_SYSTEMS
    if (new_is_cartesian) {
      if (sum(idx <- is.na(opt$zenit))) {
         zenit <- (if (grid) x[1, 1:2] + x[2, 1:2] * (x[3, 1:2] - 1) else 
                   if (opt$zenit[!idx] == 1) colMeans(x[, 1:2]) else
                   if (opt$zenit[!idx] == Inf) colMeans(apply(x[, 1:2], 2, range)) else
                   stop("unknown value of zenit"))
         RFoptions(zenit = zenit)
      }

      code <- switch(new_coord_system,
                     "cartesian" = CARTESIAN_COORD,
                     "gnomonic" = GNOMONIC_PROJ,
                     "orthographic" = ORTHOGRAPHIC_PROJ,
                     stop("unknown projection method")
                     )
       x <- RFfctn(RMtrafo(new=code), x, grid=grid, 
                   coords.new_coordunits=global.units,
                   coords.new_coord_system = "keep")
      
       if (length(y) != 0)         
         y <- RFfctn(RMtrafo(new=code), y, grid=grid, 
                   coords.new_coordunits=global.units,
                   coords.new_coord_system = "keep")
     
      if (new_coord_system == "cartesian") {
        Raumdim <- max(3, Raumdim)
        spatialdim <- Raumdim
      }
      dim(x) <- c(length(x) /Raumdim, Raumdim)
      #x <- t(x)

      ## never try to set the following lines outside the 'if (new_coord_system'
      ## as in case of ..="keep" none of the following lines should be set
      RFoptions(coords.coord_system = 
                if (new_is_cartesian) "cartesian" else new_coord_system)
      grid <- FALSE
    } else if (!(new_coord_system %in% c("keep", "sphere", "earth"))) {
      warning("unknown new coordinate system")
    }
  }

  if (Zeit <- length(T)!=0) {
    Ttriple <- length(T) == 3;
    if (length(T) <= 2) Tgrid <- TRUE
      else {
        dT <- diff(T)
        Tgrid <- max(abs(diff(dT))) < dT[1] * RFopt$general$gridtolerance
      }
    if (is.na(RFopt$general$Ttriple)) {
      if (Ttriple && Tgrid)
        stop("ambiguous definition of 'T'. Set RFoptions(Ttriple=TRUE) or ",
             "RFoptions(Ttriple=FALSE)")
      if (!Ttriple && !Tgrid) stop("'T' does not have a valid format")
    } else if (RFopt$general$Ttriple) {
      if (!Ttriple)
        stop("'T' is not given in triple format 'c(start, step, length)'")
      Tgrid <- FALSE
    } else {
      if (!Tgrid) stop("'T' does not define a grid")
      Ttriple <- FALSE
    }
    if (Tgrid)
      T <- as.vector(seq2grid(list(T), "T", Tgrid,
                              RFopt$internal$warn_ambiguous,
                              RFopt$general$gridtolerance))
    restotal <- restotal * T[3]
  }

  if (!missing(dim) && !is.null(dim) && spatialdim != dim) {
    stop("'dim' should be given only when 'distances' are given. Here, 'dim' contradicts the given coordinates.")
  }

  if (ex.red) {
    if (grid) {
      x[3, ] <- pmin(x[3, ], ex.red)
      if (length(y) > 0) y[3, ] <- pmin(y[3, ], ex.red)
      restotal <- as.integer(prod(x[3, ]))
    } else {
      len <- restotal <- as.integer(min(nrow(x), ex.red^spatialdim))
      x <- x[1:len, , drop=FALSE]
      if (length(y) > 0) y <- y[1:len, , drop=FALSE]
    }
    
    if (Zeit) {
      T[3] <- min(T[3], 3)
      restotal <- as.integer(restotal * T[3])
    }
  }
 
   
  ## keep exactly the sequence up to 'grid'
  if (length(x) > 0) {
    if (storage.mode(x) != "double") storage.mode(x) <- "double"
  } else x <- double(0)
  if (length(y) > 0) {
    if (storage.mode(y) != "double") storage.mode(y) <- "double"
  } else y <- double(0)

  L <- list(x=x, #0
            y=y, #1
            T=as.double(T), #2
            grid=as.logical(grid), #3
            spatialdim=as.integer(spatialdim), #4
            Zeit=Zeit, #5
            dist.given=FALSE, #6
            restotal=as.integer(restotal), ## 7, nr of locations
            l=as.integer(len),             ## 8, physical "length/rows" of input
            coordunits = curunits,  #9
            new_coordunits = newunits) #10
  class(L) <- "CheckXT"
  return(L)  
}


trafo.to.C_CheckXT <- function(new) {
  if (is.list(new[[1]])) {
    for(i in 1:length(new)) {
      if (length(new[[i]]$x)>0 && !new[[i]]$grid) new[[i]]$x = t(new[[i]]$x)
      if (length(new[[i]]$y)>0 && !new[[i]]$grid) new[[i]]$y = t(new[[i]]$y)
    }
  } else {
    if (length(new$x)>0 && !new$grid) new$x = t(new$x)
    if (length(new$y)>0 && !new$grid) new$y = t(new$y)
  }
  new
}


C_CheckXT <- function(x, y=NULL, z=NULL, T=NULL, grid, distances=NULL,
                      dim=NULL, # == spatialdim!
                      length.data,
                      y.ok = FALSE, 
                      printlevel = RFoptions()$general$printlevel){
 
  neu <- CheckXT(x=x, y=y, z=z, T=T, grid=grid, distances=distances,
            dim=dim, length.data=length.data, y.ok=y.ok, 
            printlevel = printlevel)
  return(trafo.to.C_CheckXT(neu))
}
   
    

RFearth2cartesian <- function(coord, units=NULL, system = "cartesian",
                              grid=FALSE) {
  if (is.character(system)) system <- pmatch(system, ISONAMES) - 1
  stopifnot(system %in%
            c(CARTESIAN_COORD, GNOMONIC_PROJ, ORTHOGRAPHIC_PROJ))
  if (is.null(units)) {
    global.units <- RFoptions()$coords$new_coordunits[1]
    units <- if (global.units[1] == "") "km" else global.units 
  }
  if (!is.matrix(coord)) coord <- t(coord)
  res <- RFfctn(RMtrafo(new=system), coord, grid=grid,
                coords.new_coord_system = "keep",
                coords.new_coordunits=units,
                coords.coord_system="earth")
  dimnames(res) <- list(NULL, c("X", "Y", "Z", "T")[1:ncol(res)])
  return(res)
}

RFearth2dist <- function(coord, units=NULL, system="cartesian",
                         grid=FALSE, ...) {
  if (is.character(system)) system <- pmatch(system, ISONAMES) - 1
  stopifnot(system %in%
            c(CARTESIAN_COORD, GNOMONIC_PROJ, ORTHOGRAPHIC_PROJ))
  if (is.null(units)) {
    global.units <- RFoptions()$coords$new_coordunits[1]
    units <- if (global.units[1] == "") "km" else global.units 
  }
  if (!is.matrix(coord)) coord <- t(coord)
  z <- RFfctn(RMtrafo(new=system), coord, grid=grid,
                coords.new_coord_system = "keep",
                coords.new_coordunits=units,
                coords.coord_system="earth")
  return(dist(z, ...))
}




## used by RFratiotest, fitgauss, Crossvalidation, likelihood-ratio,  RFempir
StandardizeData <- function(model,
                            x, y=NULL, z=NULL, T=NULL,
                            grid, data, distances=NULL, 
                            RFopt, mindist_pts=2,
                            dim=NULL, allowFirstCols=TRUE, vdim = NULL, ...) {
  #if (missing(x)) Print(data, T) else Print(data, T, x)
  RFoptions(internal.examples_reduced=FALSE)

  #Print(data); if (!missing(x)) print(x); Print(missing(x),  y, z, T, missing(dim), missing(grid), missing(distances))
  
  if (missing(dim)) dim <- NULL
  if (missing(grid)) grid <- NULL
  
  dist.given <- !missing(distances) && length(distances)>0
 
  matrix.indep.of.x.assumed <- FALSE  
  rangex <- neu <- gridlist <- RFsp.coord <-  gridTopology <- data.RFparams <-
    mindist <- data.col <- NULL
  if (missing(data)) stop("missing data")
  missing.x <- missing(x)  
  
  if (isSpObj(data)) data <- sp2RF(data)
  if (isRFsp <- is(data, "RFsp") || (is.list(data) && is(data[[1]], "RFsp"))){    
    ## ||(is.list(data) && is(data[[1]], "RFsp")))
    if ( (!missing.x && length(x)!=0) || length(y)!=0   || length(z) != 0 ||
       length(T) !=  0 || dist.given || length(dim)!=0 || length(grid) != 0)
      stop("data object already contains information about the locations. So, none of 'x' 'y', 'z', 'T', 'distance', 'dim', 'grid' should be given.")
    if (!is.list(data)) data <- list(data)
    sets <- length(data)
    x <- RFsp.coord <- gridTopology <- data.RFparams <- vector("list", sets)
    
    if (!is.null(data[[1]]@.RFparams)) {
      if (length(vdim) > 0) stopifnot( vdim == data[[1]]@.RFparams$vdim)
      else vdim <- data[[1]]@.RFparams$vdim
    }
   
    dimdata <- NULL
    dimensions <- (if (isGridded(data[[1]])) data[[1]]@grid@cells.dim
                   else nrow(data[[1]]@data))

    dimensions <- c(dimensions, data[[1]]@.RFparams$vdim)
    for (i in 1:length(data)) {
      xi <- list()
      xi$grid <- isGridded(data[[i]])
      compareGridBooleans(grid, xi$grid)
      
      data[[i]] <- selectDataAccordingFormula(data[[i]], model=model)
      data.RFparams[[i]] <- data[[i]]@.RFparams
      gridTopology[[i]] <- if (xi$grid) data[[i]]@grid else NULL
      RFsp.coord[[i]] <- if (!xi$grid) data[[i]]@coords else NULL
      
      dimensions <- if (xi$grid) data[[i]]@grid@cells.dim else nrow(data[[i]]@data)
      dimensions <- c(dimensions, data[[i]]@.RFparams$vdim)
      if (RFopt$general$vdim_close_together) dimensions <- rev(dimensions)
      dimdata <- rbind(dimdata, c(dimensions, data[[i]]@.RFparams$n))      
      tmp <- RFspDataFrame2conventional(data[[i]])
      
      xi$x <- tmp$x
      if (!is.null(tmp$T)) xi$T <- tmp$T
      data[[i]] <- as.matrix(tmp$data)   
      x[[i]] <- xi
   }

    idx <- if (RFopt$general$vdim_close_together) 1 else length(dimensions)
    if (all(dimdata[, idx] == 1))
      dimdata <- dimdata[, -idx, drop=FALSE]
    if (all(dimdata[, ncol(dimdata)] == 1)) # repet
      dimdata <- dimdata[, -ncol(dimdata), drop=FALSE]
    
  } else { # !isRFsp
    ## dimdata wird spaeter bestimmt
   
    if (dist.given) {
      stopifnot(missing(x) || length(x)==0, length(y)==0, length(z)==0)
      if (!is.list(distances)) {
        distances <- list(distances)
        if (is.list(data))
          stop("if list of data is given then also for distances ")
        data <- list(as.matrix(data))
      } else if (!is.list(data)) {
        stop("if list of distances is given then also for data ")
        if (length(data) != length(distances))
          stop("length of distances does not match length of data")
      }      
      for (i in 1:length(distances)) {
        if (any(is.na(data)))
          stop("missing data are not allowed if distances are used.")
      }

      stopifnot(missing(T) || length(T)==0)
      if (is.matrix(distances[[1]])) {

        dimensions <- sapply(distances, nrow)
        spatialdim <- tsdim <- xdimOZ <- dimensions[1]
        if (length(dim) > 0 && dim != spatialdim)
          stop("unclear specification of the distances: either the distances is given as a vector or distance vectors should given, where the number of rows matches the spatial dimension")
        lcc <- sapply(distances, function(x) 0.5 * (1 + sqrt(1 + 8 * ncol(x))) )
        if (!all(diff(dimensions) == 0))
          stop("sets of distances show different dimensions")
        range_distSq <- function(M) range(apply(M, 2, function(z) sum(z^2)))
        rangex <- sqrt(range(sapply(distances, range_distSq)))
      } else {        
        xdimOZ <- 1L
        spatialdim <- tsdim <- as.integer(dim)      
        lcc <- sapply(distances, function(x) if (is.matrix(x)) -1
                                        else 0.5 * (1 + sqrt(1 + 8* length(x))))
        rangex <- range(sapply(distances, range))
      }
#      Print(mindist, rangex, RFopt$nugget$tol)      
      mindist <- min(rangex)
      if (is.na(mindist)) mindist <- 1 ## nur 1 pkt gegeben, arbitraerer Wert
      if (mindist <= RFopt$nugget$tol) {
        if (!RFopt$general$allowdistanceZero)
          stop("distance with value 0 identified -- use allowdistanceZero=T?")
        mindist <- 1e-15 * (RFopt$nugget$tol == 0) + 2 * RFopt$nugget$tol

        for (i in 1:length(distances))
          if (is.vector(distances[[i]]))
            distances[[i]][distances[[i]] == 0] <- mindist
          else distances[[i]][1, apply(distances[[i]], 2,
                                       function(z) sum(z^2))] <- mindist
      }

      len <- as.integer(lcc)
      if (any(len != lcc)) stop("number of distances not of form k(k-1)/2")
      neu <- CheckXT(distances=distances, dim = spatialdim) 
      coordunits <- RFopt$coords$coordunits
      Zeit <- FALSE      
    } else { ## distances not given
       if (is.data.frame(data) || !is.list(data)) {
 #         Print(missing(x), x, data, is.data.frame(data), !is.list(data))
        if (!missing(x) && is.list(x) && !is.data.frame(x) &&
             (length(x$grid)==0 || length(x$restot)==0))
           stop("either both coordinates and data must be lists or none")
         data <- list(data)
      }

      sets <- length(data)
      for (i in 1:sets) {
        if (is.data.frame(data[[i]]) || is.vector(data[[i]]))
          data[[i]] <- as.matrix(data[[i]])
      }
      sel <- try(selectAccordingFormula(data[[1]], model), silent=TRUE)
      if (is(sel, "try-error")) sel <- NULL
      
       if (missing(x)) { ## dec 2012: matrix.indep.of.x.assumed        
        if (!is.null(dnames <- colnames(data[[1]]))) {
          if ((!any(is.na(xi <- RFopt$coord$coordnames))) ||
              (length(xi <- earth_coordinate_names(dnames)) == 2) ||
              (length(xi <- cartesian_coordinate_names(dnames)) > 0) ||
              (length(xi <- general_coordinate_names(dnames)) > 0)  ) {
            x <- list()
            for (i in 1:sets) {
              xx <- data[[i]][ , xi, drop=FALSE]
              storage.mode(xx) <- "numeric"
              x[[i]] <- list(x=xx, grid = FALSE)
              if (length(sel) == 0) sel <- -xi 
            }
          }         
        }

        if (missing(x)) { ## if still missing
          data.col <- try(data.columns(data[[1]], xdim=dim,
                                       force=allowFirstCols,
                                       halt=!allowFirstCols))
         
          x <- list()
          if (is(data.col, "try-error")) {
            if (length(sel) > 0){
              for (i in 1:sets) {
                x[[i]] <- data[[i]][ , !sel, drop=FALSE]
                storage.mode(x[[i]]) <- "numeric"
              }
              if (length(dim) == 0) {
                warning("better give 'dim' explicitely.")
              }
              if (length(dim) > 0 && ncol(x[[i]]) != dim)
                stop("'dim' does not match the recognized coordindates")
            } else {
              sel <- TRUE
              data.col <- NULL
              matrix.indep.of.x.assumed <- TRUE
              for (i in 1:sets) {
                x[[i]] <- 1:nrow(as.matrix(data[[i]]))
                storage.mode(x[[i]]) <- "numeric"
              }
            }
### x[1] <- 0 ## so no grid !  ## why forbidding ?? 15.5.2015
          } else {
            for (i in 1:sets) {
              xx <- data[[i]][, data.col$x, drop=FALSE]
               storage.mode(xx) <- "numeric"
              x[[i]] <- list(x=xx, grid=FALSE)
              if (length(sel) == 0) sel <- data.col$data
            }
          }
        }

        for (i in 1:sets) {
          data[[i]] <- data[[i]][ , sel, drop=FALSE]
          storage.mode(data[[i]]) <- "numeric"
        }
        
      } ## xgiven; KEIN ELSE, auch wenn nachfolgend z.T. gedoppelt wird

    
      if (is.data.frame(x)) x <- as.matrix(x)   
      if (is.list(x)) {
        if (length(y)!=0 || length(z)!=0 || length(T)!=0)
          stop("if x is alist 'y', 'z', 'T' may not be given")       

       if (!is.list(x[[1]])) {
         if (length(data) == 1) x <- list(x)
         else stop("number of sets of 'x' and 'data' differ")
       }
       
        
       } else {
        x <- list(x=x)
        if (length(y)!=0) {
          stopifnot(!is.list(y))
          x$y <- y
        }
        if (length(z)!=0) {
          stopifnot(!is.list(z))
          x$z <- z
        }
        if (length(T)!=0) {
          stopifnot(!is.list(T))
          x$T <- T
        }
        if (!is.null(grid))
          x$grid <- grid
        if (!is.list(data)) data <- list(as.matrix(data))
        x <- list(x)          
      }
      ##}
    } # ! distance
    sets <- length(data)
    dimdata <- matrix(nrow=sets, ncol=length(base::dim(data[[1]])))
    for (i in 1:sets) dimdata[i, ] <- base::dim(data[[i]])
    
  } # !isRFsp

 
  if (!dist.given) { ##   x coordinates, not distances

     
    neu <- CheckXT(x=x) #, y=y, z=z, T=T, grid=grid, distances=distances,
# dim=dim, length) # , length.data=length(data[[i]]), printlevel = 0)

    if (!is.list(neu[[1]])) neu <- list(neu)
    
    coordunits<- neu[[1]]$coordunits
    spatialdim <- as.integer(neu[[1]]$spatialdim)
    Zeit <- neu[[1]]$Zeit
    tsdim <- as.integer(spatialdim + Zeit)

    len <- sapply(neu, function(x)
                  (if (x$grid) prod(x$x[3, ]) else nrow(x$x)) *
                  (if (Zeit) x$T[3] else 1))

    getrange <- function(x)
      if (x$grid) rbind(x$x[1, ], x$x[1, ] + x$x[2, ] * (x$x[3, ] - 1))
      else apply(x$x, 2, range)
    rangex <- sapply(neu, getrange)

    ## falls mehrere datasets:
    if (ncol(x[[1]]$x) > 1 || is.null(x[[1]]$dist.given) || !x[[1]]$dist.given){
      rangex <- t(rangex) 
      base::dim(rangex) <- c(length(rangex) / spatialdim, spatialdim)
    }
    rangex <- apply(rangex, 2, range)

    getmindistSq <- function(x) {
      if (x$grid) sum(x$x[2,]^2)
      else if (nrow(x$x) < 2) NA
      else if (nrow(x$x) <= mindist_pts) min(dist(x$x))
      else min(dist(x$x[sample(nrow(x$x), mindist_pts), ]))
    }

    if (Zeit && any(sapply(neu, function(x) x$T[2]) <= RFopt$nugget$tol))
      stop("step of time component smaller than nugget tolerance 'tol'")
    
    if (any(sapply(neu, function(x) x$grid && any(x$x[2, ]<=RFopt$nugget$tol))))
      stop("step of some spatial component smaller than nugget tolerance 'tol'")

    zaehler <- 0

    repeat {
      mindist <- sqrt(min(sapply(neu, getmindistSq)))      
      if (is.na(mindist)) mindist <- 1 ## nur 1 pkt gegeben, arbitraerer Wert
      if (mindist <= RFopt$nugget$tol) {
        if (!RFopt$general$allowdistanceZero)
          stop("Distance with value 0 identified -- use allowdistanceZero=T?")
        if ((zaehler <- zaehler + 1) > 10)
          stop("unable to scatter point pattern")
        for (i in 1:length(neu)) if (!neu[[i]]$grid)
          neu[[i]]$x <- neu[[i]]$x + rnorm(length(neu[[i]]$x), 0,
                                           10 * RFopt$nugget$tol)
      } else break;
    }
  
    xdimOZ <- ncol(neu[[1]]$x)
  }
 
  if (length(dim) > 0) stopifnot(dim == tsdim)
  varnames <- try(colnames(data[[1]]))

## geht x[[1]]$x immer gut ??
#  Print(missing(x), neu)
  names <- GetDataNames(model=model,
                        coords=if (missing(x)) NULL else x[[1]]$x,
                        locinfo=neu[[1]]) #ohne data!

  if (is.null(names$varnames))
    names$varnames <-
      if (class(varnames) == "try-error") NULL else varnames

  restotal <- sapply(neu, function(x) x$restotal)
  ldata <- sapply(data, length)


  
  if (length(vdim) == 0) {
    if (all(sapply(data, function(x) is.vector(x) || ncol(x) == 1)))
      vdim <- 1
    else if (!missing(model)) {
      vdim <- rfInit(list("Cov", PrepareModel2(model=model, ...,
                                               x=trafo.to.C_CheckXT(neu))),
                     x=x, y=y, z=z, T=T, grid=grid, distances=distances,
                     dim=dim, reg=MODEL_AUX, dosimulate=FALSE)[1]
    } else vdim <- NA
  }
    
  repetitions <- as.integer(ldata / (restotal * vdim))
  #  Print(data, ldata, repetitions, restotal, vdim, neu, dist.given)
  if (!is.na(vdim) && any(ldata != repetitions * restotal * vdim))
    stop("mismatch of data dimensions")
  
  RFoptions(internal.examples_reduced=RFopt$internal$examples_red)
  return(list(
      ## coord = expandiertes neu # #
      model = if (missing(model)) NULL else
              PrepareModel2(model, ..., x=trafo.to.C_CheckXT(neu)),
      orig.model = if (missing(model)) NULL else model,
      data=data, dimdata=dimdata, isRFsp = isRFsp,
      RFsp.coord = RFsp.coord,
      coord = neu,
      dist.given=dist.given,
      gridTopology = gridTopology,
      data.RFparams = data.RFparams,
      spatialdim=spatialdim,
      tsdim=tsdim,
      rangex = as.matrix(rangex),
      coordunits=coordunits,
      Zeit = Zeit,
      matrix.indep.of.x.assumed = matrix.indep.of.x.assumed,
      len = len,
      mindist = mindist,
      xdimOZ = xdimOZ,
      vdim = vdim,
      coordnames=names$coordnames,
      varnames=if (length(names$varnames)==0) ""
                   else names$varnames,   
      data.col = data.col,
      repetitions = repetitions
      ))
}

