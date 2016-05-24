 
#   stop("big data sets currently not allowed") 
#    printlevel <- mle.methods <- lsq.methods <- recall <- sdvar <- general <- TRUE




## to do: grid
GetNeighbourhoods <- function(model, Z, X,
                              splitfactor, maxn, split_vec, shared=FALSE,
                              consider.correlations = TRUE) {

  model.nr <- MODEL_AUX
  rfInit(model, Z$coord, dosimulate=FALSE, reg=model.nr)
  
  lc <- Z$len
  maxn <- as.integer(maxn)        ## max number of points, including neighbours
  minimum <- as.integer(split_vec[1])## min. number of points in a neighbourhood
  splitn <- as.integer(split_vec[2]) ## number of location when split up
  maximum <- as.integer(split_vec[3])## maximum number of points when still
  ##                                    neighbours of neighbours are included.
  ##                         Note that, mostly, an additional box is included.
  locfactor <- as.integer(splitfactor * 2 + 1) ## total number of
  ##    neighbouring boxes in each direction, including the current box itself
  xdimOZ <- Z$xdimOZ
  tsdim <- Z$tsdim
  newZ <- list()
 
  if (Z$dist.given) { ## to do
   stop("space splitting for 'distances' not programmed yet")
   if (is.vector(dist)) {
      ### nur reinkopiert
      j <- 1:lc
      composite <- list()
      li <- 0
      maxi <-  (splitfactor[2] - 1) / Z$vdim 
        # mini <-  splitfactor[3] / Z$vdim
      while (length(j) >= maxi) {
        distpos <- (j[1]-1) * (lc - j[1] / 2) + 1
          distidx <- distpos : (distpos + (lc-j[1]-1))
        
        locs <- (j[1]+1):lc
        locidx <- which(!is.na(pmatch(locs, j)))
        locidx <- c(locs[locidx], j[1])            
        li <- li + 1
        composites[[li]] <- locidx
        j <- j[is.na(pmatch(j, locidx))]
      }
      if (length(j) > 0) {
        li <- li + 1
        composite[[li]] <- j
      }
      
      ## kann jetzt noch optimiert werden hinsichtlich schaetzqualitaet durch
      ## zusammenfassen
    } else {
    }
   
   ## return(result)
   return(newZ)
  }

  newZ <- list()
  for (set in 1:lc) {
    n <- as.integer(dim(Z$coord[[set]])[2])
    coord <- Z$coord[[set]]
    
    if (consider.correlations) {
      natsc <- .C("MultiDimRange", as.integer(model.nr),
                  as.integer(set),
                  natsc = double(xdimOZ),
                  PACKAGE="RandomFields")$natsc
    } else natsc <- rep(1, xdimOZ)  

    data <- Z$data[[i]]    
    if (coord$grid) {
      stop("space splitting for 'grid' not programmed yet")
      X <- cbind(coord$x, coord$T)
      dim(data) <- c(X[3, ], dim(data)[2])
      pts <- natsc / X[2, ]
      total <- prod(pts)
      factor <-(total / splitn) ^ (1/tsdim)
      pts <- pts / factor
      blocks <- round(X[3,] / pts)
      blocks[blocks < 1] <- 1
      old <- FALSE
      while(prod(ceiling(X[3, ] / blocks)) > maximum) {
        stopifnot(any(!old))
        idx <- blocks == max(blocks[!old])        
        old <- old | idx
        blocks[idx] <- blocks[idx] + 1
      }
      minpts <- trunc(X[3, ] / blocks)
      remaining <- X[3, ] - minpts * blocks
      blocknpts <- (blocks - remaining) * minpts
      listseq <- combi <- list()
      for (i in 1:tsdim) {
        combi[[i]] <- c(-1, 1)
        listseq[[i]] <- 1:blocknpts[i]
      }
      combi <- do.call(expand.grid, combi)
      for (j in 1:nrow(combi)) {
        L <- list(data)
        idx <- combi[j, ] 
        for (i in 1:tsdim) L[[i+1]] <- idx[i] * listseq[[i]]
        
        if (all(remaining[idx < 0] > 0)) {
          newZ[[length(newZ) + 1]] <- Z[[i]]
          newZ[[length(newZ)]]$coords[3, ] <- minpts + (idx < 0)
          newZ[[length(newZ)]]$data <- matrix(ncol=ncol(data), do.call("[", L))
        }
      }
    } else {    
      u <- numeric(xdimOZ)
      for (i in 1:xdimOZ) {
        u[i] <- length(unique(coord$x[i,]))
      }
      
      Range <- apply(coord, 1, range)
      rd <- apply(Range, 2, diff) 
      len <- pmax(1e-10 * max(rd), rd * (1 + 1e-10))
      units <- pmax(1, len * natsc)
      nDsplitn <- n / splitn
      
      ## * "gerechte" Aufteilung in alle Richtungen waere nDsplitn
      ## * die Richtung in die viele units sind, soll eher aufgespalten werden
      ## * ebenso : wo viele Werte sind eher aufspalten
      idx <- (nDsplitn / prod(units * u))^{1/xdimOZ} * units * u > 0.5 
      reddim <- sum(idx)
      units <- units[idx]
      zaehler <- 1
      parts <- rep(1, xdimOZ)
      OK <- integer(1)
      
      repeat {
        parts[idx] <- (nDsplitn / prod(units))^{1/reddim} *
          locfactor * zaehler * units * Z$vdim
        parts <- as.integer(ceiling(parts))
        
        ## zuordnung der coordinaten_Werte zu den jeweiligen "parts"
        ## given ist liegend
        coord.idx <- floor((coord$x - Range[1,]) / (len / parts))
        
        cumparts <- cumprod(parts)
        totparts <- as.integer(cumparts[length(cumparts)])
        Ccumparts <- as.integer(c(1, cumparts))
        cumparts <- Ccumparts[-length(Ccumparts)]
        
        ## elms.in.boxes <- integer(totparts)  
        ## neighbours <- integer(totparts)
        cumidx <- as.integer(colSums(coord.idx * cumparts))
        
        elms.in.boxes <- .Call("countelements", cumidx, n, totparts,
                               PACKAGE="RandomFields")
        
        neighbours <- .Call("countneighbours", xdimOZ, parts, locfactor,
                            Ccumparts, elms.in.boxes, PACKAGE="RandomFields")
      
        ## if there too many points within all the neighbours, then split
        ## into smaller boxes
        zaehler <- zaehler * 2
      
        ## image(neighbours, zlim=c(0:(prod(parts)-1)))
        if (!is.null(neighbours)) break;
      }
    
      l <- list()
      l[[1]] <- .Call("getelements", cumidx, xdimOZ, n, Ccumparts,
                      elms.in.boxes, PACKAGE="RandomFields")
      l1len <- sapply(l[[1]], length)
    }
    
    if (length(X$x) > 0) {
      if (X$grid) {
        stop("not programmed yet")
      } else { 
        ## now calculate the boxes for the locations where we will interpolate
        i <- pmax(0, pmin(parts-1,
                          floor((t(coord$x) - Range[1,]) / (len / parts))))
        dim(i) <- rev(dim(coord$x))
        i <- as.integer(colSums(i * cumparts))
      
        res.in.boxes <- .C("countelements", i, nrow(coord$x), totparts,
                         PACKAGE="RandomFields")
        
        notzeros <- res.in.boxes > 0
        l[[3]] <-
          .Call("getelements", i, xdimOZ, as.integer(nrow(coord$x)),
                Ccumparts, res.in.boxes, PACKAGE="RandomFields")[notzeros]
        ## TO DO : idx[[3]] passt nicht, da sowohl fuer Daten
        ##          als auch coordinaten verwendet wird. Bei repet > 1
        ##         ist da ein Problem -- ueberpruefen ob repet=1
        
        
        ll <- .Call("getneighbours", xdimOZ, parts, locfactor, Ccumparts,
                    neighbours, PACKAGE="RandomFields")[notzeros]
        less <-
          sapply(ll, function(x) sum(elms.in.boxes[x]) < minimum) | !shared
        ##                  if !shared then all(less)==TRUE
      
        if (any(less)) {
          not.considered.yet <- sapply(l[[1]], length) > 0   
          newll <- ll
          for (i in which(less)) {
            current <- ll[[i]]
            elements <- sum(elms.in.boxes[current] *
                            (shared | not.considered.yet[current]))# number of pts in a neighbourhood
            while (elements < minimum) {
              new <- unique(unlist(ll[current])) # neighbours of neighbours, but not
              new <- new[which(is.na(pmatch(new, current)))]# neighbours themselves
              nn <- elms.in.boxes[new] * (shared | not.considered.yet[new]) # how many pts are in each of these boxes?
              ordr <- order(nn)
              new <- new[ordr]
              nn <- nn[ordr]
              cs <- elements + cumsum(nn)
              smaller <- sum(cs <= maximum) ## now, check which neighbours of
              ## the neigbours can be included in the list of neighbours of i
              ## to increase the number of points in the kriging neighbourhood
              if (smaller == 0) break; ## none
              if (smaller == length(cs) || cs[smaller] >= minimum ||
                  cs[smaller+1] > maxn) {
                if ( (elements <- cs[length(cs)]) <= maxn ) {            
                  current <- c(current, new)            
                } else {
                  current <- c(current, new[1:smaller])
                  elements <- cs[smaller]
                }
                if (smaller != length(cs)) break
              } else {
                ## smaller < length(cs) && cs[smaller] < minimum &&
                ## cs[smaller+1]<=maxn
                ## i.e., include the next one, but there is no chance to include
                ## more within the rules.
                elements <- cs[smaller+1]
                current <- c(current, new[1:(smaller+1)])
                break;
              }
            }
            current <- current[l1len[current] > 0]
            if (!shared) current <- current[not.considered.yet[current]]
            newll[[i]] <- current
            not.considered.yet[current] <- FALSE                            
          }
          newll <- newll[sapply(newll, length) > 0]
          l[[2]] <- newll
        } else l[[2]] <- ll
      } ## locations to be estimated not on grid
    } ## locations to be estimated
  } ## for
  return(if (shared) l else lapply(l[[2]], function(x) unlist(l[[1]][x])))
}

GetComposites <- function(Z, cliquesize) {
  stopifnot(cliquesize == 2)
  return(Z)
}


BigDataSplit <- function(Z, RFopt) {
  fit <- RFopt$fit
  method <- fit$likelihood   
  if (is.na(method <- pmatch(method, RC_LIKELIHOOD_NAMES)))
    stop("unknown value for 'likelihood'.")
  method <- RC_LIKELIHOOD_NAMES[method + 1] ## + 1  ok??
  
  if (method == "full" ||
      (method %in% c("auto", "tesselation") && all(Z$len <= fit$max_neighbours)
       )) return(Z)
  if (method == "auto") {
    method <- "tesselation"
    if (RFopt$general$printlevel>=PL_IMPORTANT)
      message("Too many locations to use standard estimation methods.\n",
              "Hence an approximative methods is used. However, it is *not* ",
              "ensured\n",
              "that they will do well in detecting long memory dependencies.")
  }

  stop("not programmed yet")

  model <- list("RFfctn", Z$model)


  stop("Neigbour kriging/selection not programmed yet")

  if (method == "tesselation") {
    if (any(diff(fit$cliquesize) <= 0)) stop("in case of 'tesselation', 'cliquesize' must contain three increasing numbers")
    return(GetNeighbourhoods(model, Z=Z, 
                             splitfactor=fit$splitfactor_neighbours,
                             maxn=fit$max_neighbours,
                             split_vec=fit$cliquesize,
                             shared=FALSE)
           )
  } else if (method == "composite") {
    if (any(diff(fit$cliquesize) != 0)) stop("in case of 'composite', 'cliquesize' must be a single value")
    return(GetComposites(Z=Z, cliquesize = fit$cliquesize))
  } else stop("unknown 'likelihood' value")  
}
