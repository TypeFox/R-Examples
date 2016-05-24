"build.hessian" <-
  function(xyz, pfc.fun, fc.weights=NULL,
           sequ=NULL, sse=NULL, ss.bonds=NULL, ... )  {
    
    if(missing(xyz))
      stop("build.hessian: 'xyz' coordinates must be provided")
    
    if (!is.function(pfc.fun))
      stop("build.hessian: 'pfc.fun' must be a function")

    ## Coordinates
    xyz = as.xyz(xyz)
    if(nrow(xyz)>1)
      xyz=xyz[1,,drop=FALSE]
    xyz=as.vector(xyz)
    
    xyz    <- matrix(xyz, ncol=3, byrow=TRUE)
    natoms <- nrow(xyz)

    ## Check provided weight matrix
    if(!is.null(fc.weights)) {
      if(!is.matrix(fc.weights))
        stop("'fc.weights' must be a numeric matrix")
      
      if((nrow(fc.weights) != natoms) ||
         (ncol(fc.weights) != natoms) )
        stop("'fc.weights' must be numeric matrix with dimensions NxN")
    }

    build.submatrix <- function(xyz, natoms, 
                                fc.weights=NULL, 
                                ssdat=NULL, ...) {
     
      ## Full Hessian
      Hsm <- matrix(0, ncol=3*natoms, nrow=3*natoms)
      
      ## Indices relating atoms and columns in the sub-hessian
      col.inds <- seq(1, ncol(Hsm), by=3)

      ## Weight indices
      inds <- rep(1:natoms, each=3)

      ## Convenient indices for accessing the hessian
      inds.x <- seq(1, natoms*3, by=3)
      inds.y <- inds.x+1
      inds.z <- inds.x+2
            
      for ( i in 1:natoms ) {
        ## Calculate difference vectors and force constants
        diff.vect <- t(t(xyz) - xyz[i,])
        
        ##dists <- apply(diff.vect, 1, function(x) sqrt(sum(x**2)))
        dists <- sqrt(rowSums(diff.vect**2))  ## quicker !

        ## pfc.fun takes a vector of distances
        if("ssdat" %in% names(formals( pfc.fun )) &&
           "atom.id"     %in% names(formals( pfc.fun )) )
          force.constants <- pfc.fun(dists, atom.id=i, ssdat=ssdat, ...)
        else 
          force.constants <- pfc.fun(dists, ...)
          
        ## Scale the force constants
        if(!is.null(fc.weights)) {
          force.constants <- force.constants * fc.weights[i,]
        }
        
        force.constants <- (-1) * force.constants / (dists**2)

        ## since we divide on zero, ensure no Inf values
        force.constants[i] <- 0
        diff.vect[i,] <- 0

        ## Hessian elements
        dxx <- diff.vect[,1] * diff.vect[,1] * force.constants
        dyy <- diff.vect[,2] * diff.vect[,2] * force.constants
        dzz <- diff.vect[,3] * diff.vect[,3] * force.constants
        
        dxy <- diff.vect[,1] * diff.vect[,2] * force.constants
        dxz <- diff.vect[,1] * diff.vect[,3] * force.constants
        dyz <- diff.vect[,2] * diff.vect[,3] * force.constants

        ## Place the elements
        m <- col.inds[i]
        
        ## Off-diagonals 
        Hsm[inds.x, m   ] <- dxx
        Hsm[inds.y, m+1 ] <- dyy
        Hsm[inds.z, m+2 ] <- dzz

        Hsm[inds.y, m   ] <- dxy
        Hsm[inds.z, m   ] <- dxz
        
        Hsm[inds.x, m+1 ] <- dxy
        Hsm[inds.z, m+1 ] <- dyz
        
        Hsm[inds.x, m+2 ] <- dxz
        Hsm[inds.y, m+2 ] <- dyz

        ## Diagonal super elements
        Hsm[inds.x[i], m] <- sum(Hsm[inds.x, m]) * (-1)
        Hsm[inds.y[i], m] <- sum(Hsm[inds.y, m]) * (-1)
        Hsm[inds.z[i], m] <- sum(Hsm[inds.z, m]) * (-1)

        Hsm[inds.x[i], m+1] <- sum(Hsm[inds.x, m+1]) * (-1)
        Hsm[inds.y[i], m+1] <- sum(Hsm[inds.y, m+1]) * (-1)
        Hsm[inds.z[i], m+1] <- sum(Hsm[inds.z, m+1]) * (-1)

        Hsm[inds.x[i], m+2] <- sum(Hsm[inds.x, m+2]) * (-1)
        Hsm[inds.y[i], m+2] <- sum(Hsm[inds.y, m+2]) * (-1)
        Hsm[inds.z[i], m+2] <- sum(Hsm[inds.z, m+2]) * (-1)

      }
      return(Hsm)
    }
    
 
     
    ## Sequence-structure properties in a list
    ssdat <- list()

    ## Sequence
    if(!is.null(sequ)) {
      ssdat$seq <- sequ
    } else {
      ssdat$seq <- NULL
    }

    ## SSE
    if(!is.null(sse)) {
      if(is.null(sse$call$resno) | is.null(sse$call$full)) {
        use.sse <- FALSE
      }
      else {
        if((sse$call$resno=="F" | sse$call$resno=="FALSE") &
           (sse$call$full=="T"  | sse$call$full=="TRUE"  ))
          use.sse <- TRUE
        else
          use.sse <- FALSE
      }
      
      ssdat$sse <- sse
      if(!use.sse)
        warning("sse argument ignored: use 'dssp' with 'resno=FALSE' and 'full=TRUE'")
    }
    else {
      use.sse <- FALSE
      ssdat$sse <- NULL
    }
    
    ## SS-bonds
    if(!is.null(ss.bonds)) {
      if(class(ss.bonds)!="matrix")
        stop("ss.bonds must be two a column matrix")
      if(ncol(ss.bonds)!=2)
        stop("ss.bonds must be two a column matrix")
      ssdat$ss.bonds <- rbind(ss.bonds, ss.bonds[,2:1])
    }
    else {
      ssdat$ss.bonds <- NULL
    }

    ## Identify bridge pairs
    if(use.sse) {
      ## Helix 1-4 interactions
      ssdat$helix14      <- sse.bridges(sse, type="helix", hbond=TRUE, energy.cut=-1.0)
      ssdat$helix14      <- rbind(ssdat$helix14, ssdat$helix14[,2:1])
      
      ## Beta bridges
      ssdat$beta.bridges <- sse.bridges(sse, type="sheet", hbond=TRUE, energy.cut=-1.0)
      ssdat$beta.bridges <- rbind(ssdat$beta.bridges, ssdat$beta.bridges[,2:1])
    }      
    else {
      ssdat$helix14      <- NULL
      ssdat$beta.bridges <- NULL
    }

    H <- build.submatrix(xyz=xyz, natoms=natoms,
                         fc.weights=fc.weights,
                         ssdat=ssdat, ... )
    
    return(H)
  }
