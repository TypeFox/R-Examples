correct.degfree <- function(object)
{ daeTolerance <- get("daeTolerance", envir=daeEnv)
#check that degrees of freedom are correct
  check.df <- FALSE
  if (!is.projector(object))
    stop("Must supply a valid object of class projector")
  if (is.na(object@degfree))
    stop("Degrees of freedom are missing. Can use degfree<- to set.")
  e <- eigen(object@.Data, symmetric=T, only.values=T)
  nonzero.e <- e$values[e$values > daeTolerance[["eigen.tol"]]]
	dflen <- length(nonzero.e)
  if ((object@degfree - dflen) > daeTolerance[["eigen.tol"]])
      stop("Degrees of freedom of projector are not correct")   
  check.df <- TRUE
  check.df
 }

validProjector <- function(object)
{ daeTolerance <- get("daeTolerance", envir=daeEnv)
#checks that a matrix is a square projection operator
  Q <- object@.Data
  isproj <- TRUE
  if (nrow(Q) != ncol(Q))
     isproj <- "Matrix is not square"
     else #if (!is.allzero(t(Q)-Q))
          if (!isTRUE(all.equal(t(Q), Q, tolerance=daeTolerance[["element.tol"]])))
              isproj <- "Matrix is not symmetric"
          else 
          { #if (!is.allzero(Q%*%Q - Q))
            if (!isTRUE(all.equal(Q%*%Q, Q, tolerance=daeTolerance[["element.tol"]])))
              isproj <- "Matrix is not idempotent"
          }
   isproj
}

projector <- function(Q)
{ daeTolerance <- get("daeTolerance", envir=daeEnv)
  p <- new("projector", .Data=Q)
  validity <- validProjector(p)
  if (class(validity) == "character")
     stop(validity)
  e <- eigen(Q, symmetric=T, only.values=T)
  nonzero.e <- e$values[abs(1 - e$values) < 0.9]
  dflen <- length(nonzero.e)
  p@degfree=dflen
  p
}

setClass("projector", representation("matrix", degfree = "integer"), prototype(degfree=as.integer(NA)))
setAs(from="projector", to="matrix", 
      def=function(from){m <- from@.Data; m})
setValidity("projector", validProjector, where=".GlobalEnv")

is.projector <- function(object)
{ inherits(object, "projector") & validObject(object)
}

degfree <- function(object)
{ if (!inherits(object, "projector"))
    stop("Must supply an object of class projector")
  object@degfree
}

"degfree<-" <- function(object, value)
#A function to replace supplied or computed the degrees of freedom of the projector object
{ if (!is.projector(object))
    stop("Must assign to a valid object of class projector")
  if (length(value) == 1)
    object@degfree <- as.integer(value)
  else
  { daeTolerance <- get("daeTolerance", envir=daeEnv)
    e <- eigen(object@.Data, symmetric=T, only.values=T)
	  nonzero.e <- e$values[e$values > daeTolerance[["eigen.tol"]]]
	  object@degfree <- length(nonzero.e)
  }
  object
}

print.projector <- function(x, ...)
{ if (!inherits(x, "projector"))
    stop("Must supply an object of class projector")
  print(as(x, "matrix"))
  cat("degfree: ",x@degfree,"\n")
  invisible(x)
}
setMethod("show", "projector", function(object) print.projector(object))

proj2.decomp <- function(...)
{ .Deprecated(new = "proj2.eigen", package = "dae")
  invisible()
}

proj2.eigen <- function(Q1, Q2)
{ #A procedure to compute the eigenvalues and eigenvectors for the decomposition 
  #of Q1 pertaining to Q2  i.e. the common eigenvalues of Q1Q2Q1 and Q2Q1Q2 and the 
  #eigenvectors of Q1 in their joint dcomposition.
  #They are stored in a list with elements named efficiencies and eigenvectors
  if (!is.projector(Q1) | !is.projector(Q2))
    stop("Must supply valid objects of class projector")
  if (nrow(Q1) != nrow(Q2))
    stop("Matrices not conformable.")
  daeTolerance <- get("daeTolerance", envir=daeEnv)
  Q121 <- Q1 %*% Q2 %*% Q1
	eff <- eigen(Q121, symmetric=T)
	nonzero.eff <- eff$values[eff$values > daeTolerance[["eigen.tol"]]]
	r <- length(nonzero.eff)
	if (r==0)
	{ nonzero.eff <- 0
	  nonzero.eigen <- NULL
	}
	else
	{ nonzero.eigen <- eff$vectors[,1:r]
	  nonzero.eigen <- (abs(nonzero.eigen) > daeTolerance[["eigen.tol"]])* nonzero.eigen
	}
	list(efficiencies = nonzero.eff, eigenvectors = nonzero.eigen)
}

proj2.efficiency <- function(Q1, Q2)
{ #A procedure to compute the canonical efficiency factors (eigenvalues) in the joint 
  #decomposition of Q1 and Q2 
  proj.Q1Q2 <- proj2.eigen(Q1, Q2)
  nonzero.eff <- proj.Q1Q2$efficiencies
	nonzero.eff
}

decomp.relate <- function(decomp1, decomp2)
{ #A procedure to examine the the relationship between the eigenvectors in 
  #decomp1 and decomp2
  #decomp is a list produced by proj2.eigen or proj2.combine
  daeTolerance <- get("daeTolerance", envir=daeEnv)
  relmat <- crossprod(decomp1$eigenvectors, decomp2$eigenvectors)
  relmat <- (abs(relmat) > daeTolerance[["element.tol"]])* relmat
  dimnames(relmat) <- list(as.character(round(decomp1$efficiencies,4)), as.character(round(decomp2$efficiencies,4)))
  relmat
}

proj2.ops <- function(...)
{ .Deprecated(new = "proj2.combine", package = "dae")
  invisible()
}


"proj2.combine" <- function(Q1, Q2)
{ #A procedure to compute the Residual operator for P remove Q when P and Q are nonorthogonal.
  #  Corresponding projection operator for Q in P is also obtained.
  n <- nrow(Q1)
  if (n != nrow(Q2))
    stop("Matrices not conformable.")
  isproj <- is.projector(Q1) & is.projector(Q2)
  Qconf <- projector(matrix(0, nrow = n, ncol = n))
  Qres <- Q1
  Eff.Q1.Q2 <- 0
  eigenvec <- NULL
  if (degfree(Q1) > 0)
  { #compute efficiencies
    decomp <- proj2.eigen(Q1, Q2)
    Eff.Q1.Q2 <- decomp$efficiencies 
    if (length(Eff.Q1.Q2) == 1 & Eff.Q1.Q2[1]==0) #check matrices are orthogonal
    { warning("Matrices are orthogonal.")
    }
    else
    { daeTolerance <- get("daeTolerance", envir=daeEnv)
      EffUnique.Q1.Q2 <- remove.repeats(Eff.Q1.Q2)
      K <- length(EffUnique.Q1.Q2)
      #check for all confounded (i.e. eff = 1)
      if (K == 1 & EffUnique.Q1.Q2[1] == 1 & length(Eff.Q1.Q2) == degfree(Q2))
      { Qconf <- projector(Q2)
        Qres <- projector(Q1 - Q2)
      } else if(length(Eff.Q1.Q2) == degfree(Q1)) # all of Q1 is confounded by Q2
      { Qconf <- projector(Q1)
        Qres <- projector(matrix(0, nrow = nrow(Q1), ncol = ncol(Q1)))
      } else      #compute projection operators for partially confounded case
      { Qres <- Q1
        Q121 <- Q1 %*% Q2 %*% Q1
        for(eff in EffUnique.Q1.Q2[1:K])
        { Qres <- Qres %*% (Q1 - (Q121/eff))
          Qres <- (Qres + t(Qres))/2  #force symmetry because know that it must be
        }  
        eff <- eigen(Qres, symmetric=T)
        eff <- eff$values[eff$values > daeTolerance[["eigen.tol"]]]
        Qres <- projector(Qres)
        Qconf <- projector(Q1 - Qres)
      }
    }
    eigenvec <- decomp$eigenvectors
  }
  list(efficiencies = Eff.Q1.Q2, eigenvectors=eigenvec, Qconf = Qconf, Qres = Qres)
}
