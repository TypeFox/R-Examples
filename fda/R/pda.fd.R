pda.fd  <-  function(xfdlist, bwtlist=NULL, awtlist=NULL, ufdlist=NULL,
                     nfine=501, returnMatrix=FALSE)
{
#  PDA_FD computes the basis function expansions of the
#  estimates of the coefficient functions a_k(t) and b_j(t)
#  in the possibly nonhomogeneous linear differential operator
#
#    Lx(t) =
#       b_0(t)x(t) + b_1(t)Dx(t) + ... + b_{M-1}D^{M-1}x(t) + D^M x(t)
#       - a_1(t)u_1(t) - ... - a_k(t)u_K(t)
#
#  of order M = DIFEORDER that minimizes in a least squares sense the residual
#  functions f(t) = Lx(t).
#
#  The J equations may be of different orders and have different numbers
#     of forcing functions.  In the rather complicated description of the
#     arguments below, we use M_j to stand for the order of the jth equation.
#
#  If (DIFEORDER = 0, PDALIST fits the varying coefficient or pointwise
#  linear model using the functions x(t) as dependent variables and
#  the forcing functions u(t) as indep}ent variables.  In this case,
#  there must be at least one forcing function.
#
#  The functions x(t) are in functional data object XFDOBJ.
#  The forcing functions u_k(t) are in functional data object UFDOBJ.
#  The coefficient functions for u_k(t) and x(t) are expanded in terms of the
#  basis functions specified in AWTLIST and BWTLIST, respectively.
#
#
#  Arguments:
#  XFDLIST     list vector of length J of functional data objects for the
#                 J functions whose derivatives define the DIFE.
#  UFDLIST     list of length J whose components are themselves lists.
#                  The jth component of this list contains a list vector of
#                  length K_j of functional data objects for the
#                 independent variables or u-variables forcing the equation
#                  for the jth variable.
#              In the special univariate case where there is only a single
#                  DIFE, UFDLIST can have components which are the
#                  functional data objects for the forcing functions, rather
#                  than being a list with a single component list containing
#                  these functions.
#  BWTLIST     list vector of length J, the jth component of which is a list
#                   vector of length J, each component which is itself a list
#                   vector of length M_k equal to the order of the kth
#                   differential equation in the containing list.
#                   Each component of these lists within lists within a list
#                   is a functional parameter object defining a weighting
#                   coefficient function.
#              that is, BWTLIST is three levels or layers of lists, the top
#                   level a single list of length J, the second level a
#                   set of J lists, each corresponding to an equation, and the
#                   third level containing lists of length M_k containing
#                   coefficient functions defining the contribution of the
#                   m_kth derivative of variable k to the equation for
#                   variable j.
#              that is, if the orders of all the equations were the same
#                   and R supported list arrays, their dimensions would be
#                   J, J and M.
#              However, in the special case of J = 1, where only M
#                   coefficients are required, BWTLIST may be a simple
#                   list vector of length M.
#  AWTLIST     list of length J whose components are themselves lists.
#                 The jth component of this list contains a list vector of
#                  length K_j of functional parameter objects for the
#                  coefficient functions a_{jk}(t) multipling the
#                  corresponding forcing function in UFDLIST.
#              In the special univariate case where there is only a single
#                  DIFE, AWTLIST can have components which are the
#                  functional parameter objects for the forcing functions,
#                  rather than being a list with a single component list
#                  containing these functional parameter objects.
#  NFINE       number of sampling points for numerical integration, set by
#                  default to 501, but adjusted as required to define a mesh
#                  that is sufficient fine as to give a satisfactory
#                  approximation to an integral.

#  The value in each component of lists (within lists) XFDLIST and UFDLIST is a
#      scalar FD object.
#  The value in each component of lists (within lists) AWTLIST AND BWTLIST is a
#      scalar FDPAR object.

#  Returns:
#  BWTLIST     list structure identical to that of the argument BWTLIST
#  RESFDLIST   FD object for residual functions.
#  AWTLIST     list structure identical to that of the argument AWTLIST
#  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
#               from a call to function BsplineS.  See this function for
#               enabling this option.

#  last modified 16 April 2014 by Jim Ramsay

#  check dimensions of the lists

# check XFDLIST

  if (inherits(xfdlist, "fd")) xfdlist = list(xfdlist)

  if (!inherits(xfdlist, "list")) stop(
		"XFDLIST is neither a list or a FD object")

  nvar <- length(xfdlist)

#  ----------------------------------------------------------------
#     For efficiency, there are two versions of this code:
#     one for a single variable, and another for multiple variables.
#  ----------------------------------------------------------------

  if (nvar == 1) {

#  ----------------------------------------------------------------
#                   Single variable case
#  ----------------------------------------------------------------

    difeorder <- length(bwtlist)
    difeordp1 <- difeorder + 1

    xfdobj <- xfdlist[[1]]
    xbasis <- xfdobj$basis
    xcoef  <- xfdobj$coefs
    xrange <- xbasis$rangeval

#  check the dimensions of UFDLIST and AWTLIST and get number of forcing
#    functions NFORCE

    if (is.null(ufdlist) | is.null(awtlist)) {
      nforce  <- 0
    } else {
      if (inherits(ufdlist[[1]], "list")) {
    #  UFDLIST is a list with a single component that is a list.
    #  convert to a list of length NFORCE.
        nforce <- length(ufdlist[[1]])
        temp <- vector("list", nforce)
        for (iu in 1:nforce) temp[[iu]] <- ufdlist[[1]][[iu]]
        ufdlist <- temp
      } else {
        nforce <- length(ufdlist)
      }
      if (inherits(awtlist[[1]], "list")) {
    #  AWTLIST is a list with a single component that is a list.
    #  convert to a list of length NFORCE.
        if (length(awtlist[[1]]) != nforce)
              stop("The length of AWTLIST is incorrect.")
        temp <- vector("list", nforce)
        for (iu in 1:nforce) temp[[iu]] <- awtlist[[1]][[iu]]
        awtlist <- temp
      } else {
        if (length(awtlist) != nforce)
            stop("The length of AWTLIST is incorrect.")
      }
    }

    #  check to see if there is anything to estimate

    if (difeorder == 0 && nforce == 0)
      stop("There are no coefficient functions to estimate.")

    ncurve  <- dim(xcoef)[2]
    nbasmax <- xbasis$nbasis

    #  check UFDLIST and AWTLIST

    if (nforce > 0) {
      errorwrd <- FALSE
      for (iu in 1:nforce) {
        if (!inherits(ufdlist[[iu]], "fd")) {
          print(paste("UFDLIST[[",iu,
                      "]] is not a functional data object.",sep=""))
          errorwrd <- TRUE
        } else {
          ufdi   <- ufdlist[[iu]]
          urange <- ufdi$basis$rangeval
          #  check that urange is equal to xrange
	    if (any(urange != xrange)) {
            print(paste(
             "XRANGE and URANGE are not identical for UFDLIST[[",
                        iu,"]].",sep=""))
            errorwrd <- TRUE
          }
        }
        afdPari <- awtlist[[iu]]
    	afdi    <- afdPari$fd
    	if (!inherits(afdi, "fd")) {
            print(paste(
	       "AFDI is not a functional data object for AWTLIST[[",
                        iu,"]].",sep=""))
            errorwrd <- TRUE
        } else {
          basisi <- afdi$basis
    	    if (any(basisi$rangeval != urange)) {
              print(paste("Ranges are incompatible for AWTLIST[[",
                          iu,"]].",sep=""))
              errorwrd <- TRUE
          }
    	    nbasmax <- max(c(nbasmax,basisi$nbasis))
        }
      }
      if (errorwrd) stop("")
    }

    #  check BWTLIST

    #  convert to a single-layer list if necessary

    if (inherits(bwtlist[[1]],"list")) {
      temp <- vector("list",difeorder)
      for (j in 1:nvar) {
        if (inherits(bwtlist[[1]][[j]], "list")) {
          bwtlist[[1]][[j]] <- bwtlist[[1]][[j]][[1]]
        }
        temp[[j]] <- bwtlist[[1]][[j]]
      }
      bwtlist <- temp
    }

    #  check the components
    errorwrd <- FALSE
    for (j in 1:difeorder) {
      if (!is.null(bwtlist[[j]])) {
        bfdParj <- bwtlist[[j]]
        if (!inherits(bfdParj,"fdPar")) {
          print(paste(
            "BWTLIST[[",j,"]] is not a functional parameter object.",sep=""))
          errorwrd <- TRUE
        }  else {
          bfdj <- bfdParj$fd
          if (!inherits(bfdj, "fd")) {
            print(paste(
              "BFDJ in BWTLIST[[",j,"]] is not a functional data object.",
              sep=""))
            errorwrd <- TRUE
          } else {
            basisj <- bfdj$basis
            if (any(basisj$rangeval != xrange)) print(paste(
                    "Ranges are incompatible for BWTLIST[[",j,"]].",sep=""))
          }
        }
        nbasmax <- max(c(nbasmax,basisj$nbasis))
      }
    }
    if (errorwrd) stop("")

    #  Set up sampling values to be used in numerical integration
    #    and set up matrix of basis values.  The number of sampling
    #  NFINE is here set to a usually workable value if too small.

    if (nfine < 5*nbasmax) nfine <- 5*nbasmax

    deltax <- (xrange[2]-xrange[1])/(nfine-1)
    tx     <- seq(xrange[1],xrange[2],deltax)

    #  set up  YARRAY to hold values of x functions and their derivatives

    yarray <- array(0,c(nfine,ncurve,difeordp1))
    for (j in 1:difeordp1) yarray[,,j] <- eval.fd(tx, xfdobj, j-1, returnMatrix)

    #  set up  UARRAY to hold values of u functions

    if (nforce > 0) {
        uarray <- array(0,c(nfine,ncurve,nforce))
        for (iu in 1:nforce)
            uarray[,,iu] <- eval.fd(tx, ufdlist[[iu]], returnMatrix=returnMatrix)
    }

    #  set up array YPROD to hold mean of products of values in YARRAY

    yprod <- array(0,c(nfine,difeordp1,difeordp1))
    for (j1 in 1:difeordp1) for (j2 in 1:j1) {
        if (ncurve == 1) yprodval <- yarray[,1,j1]*yarray[,1,j2]
        else             yprodval <- apply(yarray[,,j1]*yarray[,,j2],1,mean)
        yprod[,j1,j2] <- yprodval
        yprod[,j2,j1] <- yprodval
    }

    #  set up array YUPROD to hold mean of u-variables u times
    #    x functions and their derivatives

    if (nforce > 0) {
        yuprod <- array(0,c(nfine, nforce, difeordp1))
        for (iu in 1:nforce) {
            for (j1 in 1:difeordp1) {
                if (ncurve == 1) {
                    yuprodval <- yarray[,1,j1]*uarray[,1,iu]
                } else {
                    yuprodval <- apply(yarray[,,j1]*uarray[,,iu],1,mean)
                }
                yuprod[,iu,j1] <- yuprodval
            }
        }
    }

    #  set up array UPROD to hold mean of products of u-variables u

    if (nforce > 0) {
        uprod <- array(0,c(nfine, nforce, nforce))
        for (iu in 1:nforce) for (ju in 1:iu) {
            if (ncurve == 1) uprodval <- uarray[,1,iu]*uarray[,1,ju]
            else             uprodval <- apply(uarray[,,iu]*uarray[,,ju],1,mean)
            uprod[,iu,ju] <- uprodval
            uprod[,ju,iu] <- uprodval
        }
    }

    #  set up an index array and some arrays of 1's

    onesn <- rep(1,nfine)

    #  set up array to hold coefficients for basis expansions

    if (nforce > 0) {
      aarray <- matrix(0,nfine,nforce)
    } else {           
      aarray <- NULL
    }

    barray <- matrix(0,nfine,difeorder)

    #  --------------  beginning of loop through variables  -------------------

    #  get number of coefficients to be estimated for this equation

    # loop through u-variables

    neqns  <- 0

    if (nforce > 0) {
	for (iu in 1:nforce) {
            if (!is.null(awtlist[[iu]])) {
                afdPari <- awtlist[[iu]]
                if (afdPari$estimate)
                    neqns <- neqns + afdPari$fd$basis$nbasis
            }
	}
    }

    # loop through x functions and their derivatives

    for (j1 in 1:difeorder) {
        if (!is.null(bwtlist[[j1]])) {
            bfdParj <- bwtlist[[j1]]
            if (bfdParj$estimate)
                neqns <- neqns + bfdParj$fd$basis$nbasis
        }
    }

    if (neqns < 1) stop(
                        "Number of equations to solve is not positive.")

    #  set up coefficient array and right side array for linear equation

    cmat   <- matrix(0,neqns, neqns)
    dmat   <- matrix(0,neqns, 1)

    #  evaluate default weight functions for this variable

    if (nforce > 0) {
        for (iu in 1:nforce) {
            if (!is.null(awtlist[[iu]])) {
                afdPari     <- awtlist[[iu]]
                aarray[,iu] <- eval.fd(tx, afdPari$fd,returnMatrix=returnMatrix)
            }
        }
    }

    for (j1 in 1:difeorder) {
        if (!is.null(bwtlist[[j1]])) {
            bfdParj     <- bwtlist[[j1]]
            barray[,j1] <- eval.fd(tx, bfdParj$fd, returnMatrix=returnMatrix)
        }
    }

    #  loop through equations,
    #    corresponding to rows for CMAT and DMAT

    #  loop through equations for u-variables

    mi12 <- 0
    if (nforce > 0) {
        for (iu1 in 1:nforce) {
            if (!is.null(awtlist[[iu1]])) {
                afdPari1   <- awtlist[[iu1]]
                if (afdPari1$estimate) {
                    abasisi1    <- afdPari1$fd$basis
                    abasismati1 <- getbasismatrix(tx, abasisi1, 
                                                  returnMatrix=returnMatrix)
                    mi11 <- mi12 + 1
                    mi12 <- mi12 + abasisi1$nbasis
                    indexi1 <- mi11:mi12
        #  DMAT entry for u-variable
                    weighti1 <- -yuprod[,iu1,difeordp1]
                    dmat[indexi1] <-
                        trapzmat(abasismati1,onesn,deltax,weighti1)
        # add terms corresponding to x-derivate weights
        # that are not estimated
                    for (j1 in 1:difeorder) {
                        bfdParij <- bwtlist[[j1]]
                        if (!bfdParij$estimate) {
                            weightij <- -yuprod[,iu1,j1]
                            dmat[indexi1] <- dmat[indexi1] +
                                trapzmat(abasismati1, barray[,j1],
                                         deltax, weightij)
                        }
                    }
        #  loop through weight functions to be estimated,
        #    corresponding to columns for CMAT
        #  begin with u-variables
                    mi22 <- 0
                    for (iu2 in 1:nforce) {
                        if (!is.null(awtlist[[iu2]])) {
                            afdPari2   <- awtlist[[iu2]]
                            if (afdPari2$estimate) {
                                abasisi2    <- afdPari2$fd$basis
                                abasismati2 <- getbasismatrix(tx, abasisi2,
                                                      returnMatrix=returnMatrix)
                                weighti2    <- uprod[,iu1,iu2]
                                Cprod  <- trapzmat(abasismati1, abasismati2,
                                                   deltax, weighti2)
                                mi21 <- mi22 + 1
                                mi22 <- mi22 + abasisi2$nbasis
                                indexi2 <- mi21:mi22
                #  coefficient matrix CMAT entry
                                cmat[indexi1,indexi2] <- Cprod
                            }
                        }
                    }
        #  remaining columns:
        #    loop through u-variable -- x-derivative pairs
                    mij22 <- mi22
                    for (j2 in 1:difeorder) {
                        if (!is.null(bwtlist[[j2]])) {
                            bfdParj2     <- bwtlist[[j2]]
                            if (bfdParj2$estimate) {
                                bbasisij2    <- bfdParj2$fd$basis
                                bbasismatij2 <- getbasismatrix(tx, bbasisij2,
                                                      returnMatrix=returnMatrix)
                                weightij12   <- -yuprod[,iu1,j2]
                                Cprod <- trapzmat(abasismati1,bbasismatij2,
                                                  deltax,weightij12)
                                mij21 <- mij22 + 1
                                mij22 <- mij22 + bbasisij2$nbasis
                                indexij2  <- mij21:mij22
                                cmat[indexi1,indexij2] <- Cprod
                            }
                        }
                    }
        #  add roughness penalty matrix to diagonal entries
                    lambdai1 <- afdPari1$lambda
                    if (lambdai1 > 0) {
                        Lfdobj <- afdPari1$Lfd
                        penmat <- lambdai1*eval.penalty(abasisi1, Lfdobj)
                        cmat[indexi1,indexi1] <- cmat[indexi1,indexi1] + penmat
                    }
                }
            }
        }
    }

    #  loop through equations for x-derivatives

    mij12 <- mi12
    for (j1 in 1:difeorder) {
        if (!is.null(bwtlist[[j1]])) {
            bfdParj1 <- bwtlist[[j1]]
            if (bfdParj1$estimate) {
                bbasisij1    <- bfdParj1$fd$basis
                bbasismatij1 <- getbasismatrix(tx,bbasisij1, 
                                               returnMatrix=returnMatrix)
                mij11 <- mij12 + 1
                mij12 <- mij12 + bbasisij1$nbasis
                indexij1 <- mij11:mij12
                #  DMAT entry for u-variable -- x-derivative pair
                weightij1 <- yprod[,j1,difeordp1]
                dmat[indexij1] <-
                    trapzmat(bbasismatij1,onesn,deltax,weightij1)
                #  add terms corresponding to forcing functions
                #  with unestimated coefficients
                if (nforce > 0) {
                  for (iu in 1:nforce) {
                    if (!is.null(awtlist[[iu]])) {
                      afdPari <- awtlist[[iu]]
                      if (!afdPari$estimate) {
                        weightijk <- -yuprod[,iu,j1]
                        dmat[indexij1] <- dmat[indexij1] +
                            trapzmat(bbasisij1, aarray[,iu],deltax, weightijk)
                      }
                    }
                  }
                }
            }
        }
        #  first columns of CMAT: u-variable entries
        mi22 <- 0
        if (nforce > 0) {
          for (iu2 in 1:nforce) {
            if (!is.null(awtlist[[iu2]])) {
              afdPari2 <- awtlist[[iu2]]
              if (afdPari2$estimate) {
                abasisi2    <- afdPari2$fd$basis
                abasismati2 <- getbasismatrix(tx, abasisi2, 
                                              returnMatrix=returnMatrix)
                weighti2    <- -yuprod[,iu2,j1]
                Cprod <- trapzmat(bbasismatij1,abasismati2,deltax,weighti2)
                mi21 <- mi22 + 1
                mi22 <- mi22 + abasisi2$nbasis
                indexi2 <- mi21:mi22
                cmat[indexij1,indexi2] <- Cprod
              }
            }
          }
        }
        #  remaining columns: x-derivative pairs
        mij22 <- mi22
        for (j2 in 1:difeorder) {
          if (!is.null(bwtlist[[j2]])) {
            bfdParj2  <- bwtlist[[j2]]
            bbasisij2 <- bfdParj2$fd$basis
            if (bfdParj2$estimate) {
              bbasismatij2 <- getbasismatrix(tx, bbasisij2, 
                                             returnMatrix=returnMatrix)
              weightij22   <- yprod[,j1,j2]
              Cprod <- trapzmat(bbasismatij1,bbasismatij2,deltax,weightij22)
              mij21 <- mij22 + 1
              mij22 <- mij22 + bbasisij2$nbasis
              indexij2 <- mij21:mij22
              cmat[indexij1,indexij2] <- Cprod
            }
          }
        }
        # add roughness penalty matrix to diagonal entries
        lambdaj1 <- bfdParj1$lambda
        if (lambdaj1 > 0) {
          Lfdobj <- bfdParj1$Lfd
          penmat <- lambdaj1*eval.penalty(bbasisij1, Lfdobj)
          cmat[indexij1,indexij1] <- cmat[indexij1,indexij1] + penmat
        }
    }

    #  --------------  end of loop through variables  -------------------

    # solve for coefficients of basis expansions

    dvec <- -symsolve(cmat,dmat)

    #  set up u-function weight functions

    mi2 <- 0
    if (nforce > 0) {
        for (iu in 1:nforce) {
            if (!is.null(awtlist[[iu]])) {
                afdPari <- awtlist[[iu]]
                if (afdPari$estimate) {
                    mi1 <- mi2 + 1
                    mi2 <- mi2 + afdPari$fd$basis$nbasis
                    indexi <- mi1:mi2
                    afdPari$fd$coefs <- as.matrix(dvec[indexi])
                    awtlist[[iu]] <- afdPari
                }
            }
        }
    }

    #  set up X-function derivative weight functions

    mij2 <- mi2
    for (j in 1:difeorder) {
      if (!is.null(bwtlist[[j]])) {
        bfdParj <- bwtlist[[j]]
        if (bfdParj$estimate) {
          mij1 <- mij2 + 1
          mij2 <- mij2 + bfdParj$fd$basis$nbasis
          indexij <- mij1:mij2
          bfdParj$fd$coefs <- as.matrix(dvec[indexij])
          bwtlist[[j]] <- bfdParj
        }
      }
    }

    #  set up residual list RESFDLIST

    #  initialize with highest order derivative for this variable
    resmat  <- eval.fd(tx, xfdobj, difeorder, returnMatrix)
    #  add contributions from weighted u-functions
    if (nforce > 0) {
      onesncurve <- rep(1,ncurve)
      for (iu in 1:nforce) {
	if (!is.null(awtlist[[iu]])) {
          afdPari  <- awtlist[[iu]]
          aveci    <- as.vector(eval.fd(tx, afdPari$fd, 
                                        returnMatrix=returnMatrix))
          umati    <- eval.fd(tx, ufdlist[[iu]], returnMatrix=returnMatrix)
          aumati   <- outer(aveci,onesncurve)*umati
          resmat   <- resmat - aumati
        }
      }
    }
    #  add contributions from weighted x-function derivatives
    for (j in 1:difeorder) {
        if (!is.null(bwtlist[[j]])) {
            bfdParj <- bwtlist[[j]]
            bmatij  <- as.vector(eval.fd(tx, bfdParj$fd, 
                                         returnMatrix=returnMatrix))
            xmatij  <- eval.fd(tx, xfdobj, j-1, returnMatrix)
            resmat  <- resmat + bmatij*xmatij
        }
    }
    #  set up the functional data object
    resbasis <- xbasis
    resfd    <- smooth.basis(tx, resmat, resbasis)$fd
    resfdnames      <- xfdobj$fdnames
    resfdnames[[2]] <- "Residual function"
    resfdnames[[3]] <- "Residual function value"
    resfd$fdnames   <- resfdnames
    resfdlist       <- list(resfd)

    #  ----------------------------------------------------------------
    #                   End of single variable case
    #  ----------------------------------------------------------------

  } else {

    #  ----------------------------------------------------------------
    #                   Multiple variable case
    #  ----------------------------------------------------------------

    #  check the dimensions of UFDLIST and AWTLIST

    if (is.null(ufdlist) || is.null(awtlist)) {
      awtlist <- NULL
    } else {
      if (length(ufdlist) != nvar)
        stop(paste("The length of UFDLIST",
                   " does not match that of XFDLIST."))
      errorwrd = FALSE
      for (j in 1:nvar) {
        if (!is.null(ufdlist[[j]])) {
          nforce <- length(ufdlist[[j]])
          if (length(awtlist[[j]]) != nforce) {
            print(paste("The length of AWTLIST[[",j,
                        "]] is incorrect.",sep=""))
            errorwrd = TRUE
          }
        }
      }
      if (errorwrd) stop("")
    }

    #  check the dimensions of BWTLIST

    if (length(bwtlist) != nvar) stop("Length of BWTLIST is incorrect.")
    errorwrd = FALSE
    for (ivar in 1:nvar) {
      if (length(bwtlist[[ivar]]) != nvar) {
          print(paste("The length of BWTLIST[[",ivar,
                      "]] is incorrect.",sep=""))
          errorwrd = TRUE
      }
    }
    if (errorwrd) stop("")

    #  check XFDLIST and extract NCURVE and XRANGE

    xfd1       <- xfdlist[[1]]
    xcoef1     <- xfd1$coefs
    xbasis1    <- xfd1$basis
    xrange1    <- xbasis1$rangeval
    ncurve     <- dim(xcoef1)[2]
    resfdnames <- xfd1$fdnames

    errorwrd = FALSE
    for (ivar in 1:nvar) {
      xfdi    <- xfdlist[[ivar]]
      xcoefi  <- xfdi$coefs
      xbasisi <- xfdi$basis
      xrangei <- xbasisi$rangeval
      ncurvei <- dim(xcoefi)[2]
      if (!inherits(xfdi, "fd")) {
        print(paste("XFDLIST[[",ivar,
                    "]] is not a functional data object.",sep=""))
        errorwrd = TRUE
      } else {
        if (any(xrangei != xrange1)) {
          print("Ranges are incompatible for XFDLIST.")
          errorwrd = TRUE
        }
        if (ncurvei != ncurve) {
            print("Number of curves is incompatible for XFDLIST.")
            errorwrd = TRUE
        }
      }
    }
    if (errorwrd) stop("")

    nbasmax <- xbasis1$nbasis
    #  This will be the maximum number of basis functions

    #  check compatibility of UFDLIST and AWTLIST

      if (!(is.null(ufdlist) || is.null(awtlist))) {
        urange <- ufdlist[[1]]$basis$rangeval
        errorwrd <- FALSE
        for (ivar in 1:nvar) {
          if (!is.null(ufdlist[[ivar]])) {
            for (iu in 1:length(ufdlist[[ivar]])) {
              ufdiviu <- ufdlist[[ivar]][[iu]]
              if (!inherits(ufdiviu, "fd")) {
                print(paste("UFDLIST[[",ivar,",",iu,
                            "]] is not a functional data object.",
                            sep=""))
                errorwrd <- TRUE
              }
              if (any(ufdiviu$basis$rangeval != urange)) {
                  print("Ranges are incompatible for UFDLIST.")
                  errorwrd <- TRUE
              }
              awtfdPari <- awtlist[[ivar]][[iu]]
              if (!inherits(awtfdPari, "fdPar")) {
                print(paste("AWTFDPAR[[",ivar,"]][[",iu,
                      "]] is not a functional parameter object.",sep=""))
                errorwrd <- TRUE
              }
              afdi   <- awtfdPari$fd
              basisi <- afdi$basis
              if (any(basisi$rangeval != urange)) {
                  print("Ranges are incompatible for AWTLIST.")
                  errorwrd <- TRUE
              }
              nbasmax <- max(c(nbasmax,basisi$nbasis))
            }
            if (errorwrd) stop("")
          }
        }
      }

      #  check BWTLIST

      errorwrd <- FALSE
      for (ivar1 in 1:nvar) {
        for (ivar2 in 1:nvar) {
          difeorder <- length(bwtlist[[ivar1]][[ivar2]])
          for (j in 1:difeorder) {
            if (!is.null(bwtlist[[ivar1]][[ivar2]][[j]])) {
              bfdPari1i2j <- bwtlist[[ivar1]][[ivar2]][[j]]
              if (!inherits(bfdPari1i2j, "fdPar")) {
                  print(paste("BWTLIST[[",ivar1, ",",ivar2, ",",j,
                      "]] is not a functional parameter object.",sep=""))
                  errorwrd = TRUE
              }
              basisi1i2j <- bfdPari1i2j$fd$basis
              if (any(basisi1i2j$rangeval != xrange1)) {
                  print(paste("Ranges are incompatible for BWTLIST[[",
                              ivar1,"]][[",ivar2,"]][[",
                              j,"]]",sep=""))
                  errorwrd <- TRUE
              }
              nbasmax <- max(c(nbasmax,basisi1i2j$nbasis))
            }
          }
        }
    }
    if (errorwrd) stop("")

    #  set up sampling values to be used in numerical integration
    #    and set up matrix of basis values.  The number of sampling
    #  NFINE is here set to a usually workable value if too small.

    if (nfine < 5*nbasmax) nfine <- 5*nbasmax

    deltax <- (xrange1[2]-xrange1[1])/(nfine-1)
    tx     <- seq(xrange1[1],xrange1[2],deltax)

    #  set up  YARRAY to hold values of x functions and their derivatives

    yarray <- vector("list", 0)
    for (ivar in 1:nvar) {
        difeorder <- length(bwtlist[[ivar]][[ivar]])
        difeordp1 <- difeorder + 1
        yarray[[ivar]] <- array(0,c(nfine,ncurve,difeordp1))
        for (j in 1:difeordp1){
            yj <- eval.fd(tx, xfdlist[[ivar]], j-1, returnMatrix=returnMatrix)
            yarray[[ivar]][,,j] <- as.matrix(yj)
        }
    }

    #  set up  UARRAY to hold values of u functions

    if (!is.null(ufdlist)) {
      uarray <- vector("list", nvar)
      for (ivar in 1:nvar) {
        if (is.null(ufdlist[[ivar]])) {
          uarray[[ivar]] <- NULL
        } else {
          nforce <- length(ufdlist[[ivar]])
          uarray[[ivar]] <- vector("list", nforce)
          for (iu in 1:nforce)
                    uarray[[ivar]][[iu]] <- matrix(0,nfine,ncurve)
        }
      }
      for (ivar in 1:nvar) {
        if (!is.null(ufdlist[[ivar]])) {
          nforce <- length(ufdlist[[ivar]])
          for (iu in 1:nforce)
            uarray[[ivar]][[iu]] <- eval.fd(tx, ufdlist[[ivar]][[iu]],
                                                returnMatrix=returnMatrix)
        }
      }
    }

    #  set up array YPROD to hold mean of products of values in YARRAY

    yprod <- vector("list", nvar)
    for (i1 in 1:nvar) yprod[[i1]] <- vector("list", nvar)
    for (i1 in 1:nvar) {
      difeord1p1 <- length(bwtlist[[i1]][[i1]]) + 1
      for (i2 in 1:nvar) {
            difeord2p1 <- length(bwtlist[[i2]][[i2]]) + 1
            yprod[[i1]][[i2]] <- array(0,c(nfine,difeord2p1,difeord2p1))
      }
    }

    for (i1 in 1:nvar) {
      difeord1p1 <- length(bwtlist[[i1]][[i1]]) + 1
      for (j1 in 1:difeordp1) {
        for (i2 in 1:nvar) {
          difeord2p1 <- length(bwtlist[[i2]][[i2]]) + 1
          for (j2 in 1:difeord2p1) {
            if (ncurve == 1) {
              yprodval <-       yarray[[i1]][,1,j1]*yarray[[i2]][,1,j2]
            } else {
              yprodval <- apply(yarray[[i1]][,,j1]*yarray[[i2]][,,j2],1,mean)
            }
            yprod[[i1]][[i2]][,j1,j2] <- yprodval
          }
        }
      }
    }

    #  set up array YUPROD to hold mean of u-variables u times
    #    x functions and their derivatives

    if (!is.null(ufdlist)) {
      yuprod <- vector("list", nvar)
      for (i1 in 1:nvar) {
        if (!is.null(ufdlist[[i1]])) {
          nforce <- length(ufdlist[[i1]])
          if (nforce > 0) {
            yuprod[[i1]] <- vector("list", nforce)
            for (iu in 1:nforce) {
              difeordp1 <- length(bwtlist[[i1]][[i1]]) + 1
              yuprod[[i1]][[iu]] <- matrix(0,nfine,difeordp1)
            }
          }
        }
      }
      onesncurve <- rep(1,ncurve)
      for (i1 in 1:nvar) {
        if (!is.null(ufdlist[[i1]])) {
          nforce <- length(ufdlist[[i1]])
          if (nforce > 0) {
            difeordp1 <- length(bwtlist[[i1]][[i1]]) + 1
            for (iu in 1:nforce) {
              for (j1 in 1:difeordp1) {
                if (ncurve == 1) {
                  yuprodval <- yarray[[i1]][,1,j1]*uarray[[i1]][[iu]]
                } else {
                  yuprodval <- apply(yarray[[i1]][,,j1]*
                                outer(uarray[[i1]][[iu]],onesncurve),1,mean)
                }
                yuprod[[i1]][[iu]][,j1] <- yuprodval
              }
            }
          }
        }
      }
    }

    #  set up array UPROD to hold mean of products of u-variables u

    if (!is.null(ufdlist)) {
      uprod <- vector("list", nvar)
      for (ivar in 1:nvar) {
        nforce <- length(ufdlist[[ivar]])
        if (nforce > 0) {
          uprod[[ivar]] <- array(0,c(nfine, nforce, nforce))
          for (iu in 1:nforce) for (ju in 1:iu) {
            uprodval <- uarray[[ivar]][[iu]]*uarray[[ivar]][[ju]]
            uprod[[ivar]][,iu,ju] <- uprodval
            uprod[[ivar]][,ju,iu] <- uprodval
          }
        }
      }
    }

    #  set up an index array and some arrays of 1"s

    onesn <- rep(1,nfine)

    #  set up array to hold coefficients for basis expansions

    #  --------------  beginning of loop through variables  -------------------

    for (ivar in 1:nvar) {

      #  get number of coefficients to be estimated for this equation

      neqns  <- 0

      # loop through u-variables  if required

      if (is.null(ufdlist) || is.null(ufdlist[[ivar]])) {
        nforce <- 0
      } else {
        nforce <- length(ufdlist[[ivar]])
        if (nforce > 0) {
          for (iu in 1:nforce) {
            afdPari <- awtlist[[ivar]][[iu]]
            if (afdPari$estimate)
                neqns <- neqns + afdPari$fd$basis$nbasis
          }
        }
      }

      # loop through x functions and their derivatives

      for (i2 in 1:nvar) {
        difeorder <- length(bwtlist[[ivar]][[i2]])
        for (j2 in 1:difeorder) {
          if (!is.null(bwtlist[[ivar]][[i2]][[j2]])) {
            bfdParij <- bwtlist[[ivar]][[i2]][[j2]]
            if (bfdParij$estimate) neqns <- neqns + bfdParij$fd$basis$nbasis
          }
        }
      }
      if (neqns < 1)  stop("Number of equations to solve is not positive.")

      #  set up coefficient array and right side array for linear equation

      cmat   <- matrix(0,neqns, neqns)
      dmat   <- matrix(0,neqns, 1)

      #  evaluate default weight functions for this variable

      if (nforce > 0) {
        aarray <- matrix(0,nfine,nforce)
        for (iu in 1:nforce) {
          if (!is.null(awtlist[[ivar]][[iu]])) {
            afdPari <- awtlist[[ivar]][[iu]]
            aarray[,iu] <- eval.fd(tx, afdPari$fd, returnMatrix=returnMatrix)
          }
        }
      }

      barray <- vector("list", nvar)
      for (i in 1:nvar) {
        difeorder <- length(bwtlist[[ivar]][[i]])
        barray[[i]] <- array(0,c(nfine,nvar,difeorder))
        for (j in 1:difeorder) {
          if (!is.null(bwtlist[[ivar]][[i]][[j]])) {
            bfdParij <- bwtlist[[ivar]][[i]][[j]]
            bfdij    <- bfdParij$fd
            barray[[i]][,,j] <- as.matrix(eval.fd(tx, bfdij),
                                                    returnMatrix=returnMatrix)
          }
        }
      }

      #  loop through equations,
      #    corresponding to rows for CMAT and DMAT

      #  loop through equations for u-variables

      mi12 <- 0
      if (nforce > 0) {
        for (iu1 in 1:nforce) {
          if (!is.null(awtlist[[ivar]][[iu1]])) {
            afdPari1   <- awtlist[[ivar]][[iu1]]
            if (afdPari1$estimate) {
              abasisi1    <- afdPari1$fd$basis
              abasismati1 <- getbasismatrix(tx, abasisi1, returnMatrix=returnMatrix)
              mi11 <- mi12 + 1
              mi12 <- mi12 + abasisi1$nbasis
              indexi1 <- mi11:mi12
              #  DMAT entry for u-variable
              weighti1 <- -yuprod[[ivar]][[iu1]][,difeordp1]
              dmat[indexi1] <- trapzmat(abasismati1,onesn,deltax,weighti1)
              #  add terms corresponding to x-derivative weights
              #  that are not estimated
              for (i in 1:nvar) {
                difeorder <- length(bwtlist[[ivar]][[i]])
                for (j in 1:difeorder) {
                  bfdParij <- bwtlist[[ivar]][[i]][[j]]
                  if (!bfdParij$estimate) {
                    weightij <- -yuprod[[ivar]][[iu1]][,j]
                    dmat[indexi1] <- dmat[indexi1] +
                          trapzmat(abasismati1, barray[[ivar]][,,j],
                                   deltax, weightij)
                  }
                }
              }
              #  loop through weight functions to be estimated,
              #    corresponding to columns for CMAT
              #  begin with u-variables
              mi22 <- 0
              for (iu2 in 1:nforce) {
                if (!is.null(awtlist[[ivar]][[iu2]])) {
                  afdPari2   <- awtlist[[ivar]][[iu2]]
                  if (afdPari2$estimate) {
                    abasisi2    <- afdPari2$fd$basis
                    abasismati2 <- getbasismatrix(tx, abasisi2,
                                                returnMatrix=returnMatrix)
                    weighti2    <- uprod[[ivar]][,iu1,iu2]
                    Cprod       <- trapzmat(abasismati1, abasismati2,
                                              deltax, weighti2)
                    mi21 <- mi22 + 1
                    mi22 <- mi22 + abasisi2$nbasis
                    indexi2 <- mi21:mi22
                    #  coefficient matrix CMAT entry
                    cmat[indexi1,indexi2] <- Cprod
                  }
                }
              }
              #  remaining columns:
              #    loop through u-variable -- x-derivative pairs
              mij22 <- mi22
              for (i2 in 1:nvar) {
                if (!is.null(bwtlist[[ivar]][[i2]])) {
                  difeorder <- length(bwtlist[[ivar]][[i2]])
                  for (j2 in 1:difeorder) {
                    bfdParij2   <- bwtlist[[ivar]][[i2]][[j2]]
                    if (bfdParij2$estimate) {
                      bbasisij2    <- bfdParij2$fd$basis
                      bbasismatij2 <- getbasismatrix(tx, bbasisij2,
                                                     returnMatrix=returnMatrix)
                      weightij12   <- -yuprod[[i2]][[iu1]][,j2]
                      Cprod        <- trapzmat(abasismati1,bbasismatij2,
                                           deltax,weightij12)
                      mij21 <- mij22 + 1
                      mij22 <- mij22 + bbasisij2$nbasis
                      indexij2  <- mij21:mij22
                      cmat[indexi1,indexij2] <- Cprod
                    }
                  }
                }
              }
              #  add roughness penalty matrix to diagonal entries
              lambdai1 <- afdPari1$lambda
              if (lambdai1 > 0) {
                Lfdobj <- afdPari1$Lfd
                penmat <- lambdai1*eval.penalty(abasisi1,Lfdobj)
                cmat[indexi1,indexi1] <- cmat[indexi1,indexi1] + penmat
              }
            }
          }
        }
      }

      #  end of loop through equations for u-variables

      #  loop through equations for x-derivatives

      mij12 <- mi12
      for (i1 in 1:nvar) {
        difeorder1 <- length(bwtlist[[ivar]][[i1]])
        for (j1 in 1:difeorder1) {
          if (!is.null(bwtlist[[ivar]][[i1]][[j1]])) {
            bfdParij1 <- bwtlist[[ivar]][[i1]][[j1]]
            if (bfdParij1$estimate) {
              bbasisij1    <- bfdParij1$fd$basis
              bbasismatij1 <- getbasismatrix(tx, bbasisij1, returnMatrix=returnMatrix)
              mij11 <- mij12 + 1
              mij12 <- mij12 + bbasisij1$nbasis
              indexij1 <- mij11:mij12
              #  DMAT entry for u-variable -- x-derivative pair
              weightij1 <- yprod[[i1]][[ivar]][,j1,difeordp1]
              dmat[indexij1] <- trapzmat(bbasismatij1,onesn,
                                       deltax,weightij1)
              #  add terms corresponding to forcing functions
              #  with unestimated coefficients
              if (nforce > 0) {
                for (iu in 1:nforce) {
                  if (!is.null(awtlist[[ivar]][[iu]])) {
                    afdPari <- awtlist[[ivar]][[iu]]
                    if (!afdPari$estimate) {
		                  weightijk <- yprod[,ivar,iu,j1]
		                  dmat[indexij1] <- dmat[indexij1] +
		                    trapzmat(bbasismatij1,aarray[,iu],deltax,weightijk)
		                }
		              }
	              }
              }
              #  first columns of CMAT: u-variable entries
              mi22 <- 0
              if (nforce > 0) {
                for (iu2 in 1:nforce) {
                  if (!is.null(awtlist[[ivar]][[iu2]])) {
                    afdPari2  <- awtlist[[ivar]][[iu2]]
                    if (afdPari2$estimate) {
                      abasisi2    <- afdPari2$fd$basis
                      abasismati2 <- getbasismatrix(tx, abasisi2, returnMatrix=returnMatrix)
                      weighti2    <- -yuprod[[i1]][[iu2]][,j1]
                      Cprod <- trapzmat(bbasismatij1,abasismati2,deltax,weighti2)
                      mi21 <- mi22 + 1
                      mi22 <- mi22 + abasisi2$nbasis
                      indexi2 <- mi21:mi22
                      cmat[indexij1,indexi2] <- Cprod
                    }
                  }
                }
              }
              #  remaining columns: x-derivative pairs
              mij22 <- mi22
              for (i2 in 1:nvar) {
                difeorder2 <- length(bwtlist[[ivar]][[i2]])
                for (j2 in 1:difeorder2) {
                  if (!is.null(bwtlist[[ivar]][[i2]][[j2]])) {
                    bfdParij2 <- bwtlist[[ivar]][[i2]][[j2]]
                    bbasisij2    <- bfdParij2$fd$basis
                    bbasismatij2 <- getbasismatrix(tx, bbasisij2, returnMatrix=returnMatrix)
                    weightij22   <- yprod[[i1]][[i2]][,j1,j2]
                    Cprod <- trapzmat(bbasismatij1,bbasismatij2,deltax,weightij22)
                    if (bfdParij2$estimate) {
                      mij21 <- mij22 + 1
                      mij22 <- mij22 + bbasisij2$nbasis
                      indexij2 <- mij21:mij22
                      cmat[indexij1,indexij2] <- Cprod
                    }
                  }
                }
              }
              #  add roughness penalty terms to diagonal entries
              lambdaij1 <- bfdParij1$lambda
              if (lambdaij1 > 0) {
	              Lfdobj <- bfdParij1$Lfd
	              penmat <- lambdaij1*eval.penalty(bbasisij1,Lfdobj)
	              cmat[indexij1,indexij1] <- cmat[indexij1,indexij1] +penmat
              }
            }
          }
        }
      }

      #  end of loop through x derivatives

      dvec <- -symsolve(cmat,dmat)

      #  set up u-function weight functions

      mi2 <- 0
      if (nforce > 0) {
        for (iu in 1:nforce) {
          if (!is.null(awtlist[[ivar]][[iu]])) {
            afdPari <- awtlist[[ivar]][[iu]]
            if (afdPari$estimate) {
              mi1 <- mi2 + 1
              mi2 <- mi2 + afdPari$fd$basis$nbasis
              indexi <- mi1:mi2
              afdPari$fd$coefs <- as.matrix(dvec[indexi])
              awtlist[[ivar]][[iu]] <- afdPari
            }
          }
        }
      }

      #  set up X-function derivative weight functions

      mij2 <- mi2
      for (i1 in 1:nvar) {
        difeorder <- length(bwtlist[[ivar]][[i1]])
        for (j1 in 1:difeorder) {
          if (!is.null(bwtlist[[ivar]][[i1]][[j1]])) {
            bfdParij <- bwtlist[[ivar]][[i1]][[j1]]
            if (bfdParij$estimate) {
              mij1 <- mij2 + 1
              mij2 <- mij2 + bfdParij$fd$basis$nbasis
              indexij <- mij1:mij2
              bfdParij$fd$coefs <- as.matrix(dvec[indexij])
              bwtlist[[ivar]][[i1]][[j1]] <- bfdParij
            }
          }
        }
      }
    }

    #  --------------  end of loop through variables  -------------------

    #  set up residual list RESFDLIST

    resfdlist <- vector("list", nvar)

    for (ivar in 1:nvar) {
      difeorder <- length(bwtlist[[ivar]][[ivar]])
      xfdi      <- xfdlist[[ivar]]
      resbasis  <- xfdi$basis
      #  initialize with highest order derivative for this variable
      resmat    <- eval.fd(tx, xfdi, difeorder, returnMatrix)
      #  add contributions from weighted u-functions
      onesncurve <- rep(1,ncurve)
      if (!is.null(ufdlist)) {
        nforce <- length(ufdlist[[ivar]])
        if (nforce > 0) {
    	    for (iu in 1:nforce) {
    	      if (!is.null(awtlist[[ivar]][[iu]])) {
	            afdPari  <- awtlist[[ivar]][[iu]]
              amati    <- as.vector(eval.fd(tx, afdPari$fd, 
                                                  returnMatrix=returnMatrix))
              umati    <- eval.fd(tx, ufdlist[[ivar]][[iu]], 
                                                  returnMatrix=returnMatrix)
	            if (ncurve == 1) { 
                aumati <- amati*umati
	            } else {            
                aumati <- outer(amati,onesncurve)*umati
              }
              resmat   <- resmat - aumati
            }
    	    }
        }
      }
      #  add contributions from weighted x-function derivatives
      for (i1 in 1:nvar) {
        difeorder <- length(bwtlist[[ivar]][[i1]])
        for (j1 in 1:difeorder) {
          if (!is.null(bwtlist[[ivar]][[i1]][[j1]])) {
            bfdParij <- bwtlist[[ivar]][[i1]][[j1]]
            bfdij    <- bfdParij$fd
            bvecij   <- as.vector(eval.fd(tx, bfdij, returnMatrix=returnMatrix))
            if (ncurve == 1) {
              bmatij <- bvecij
            }  else  {
              bmatij <- outer(bvecij,onesncurve)
            }
            xmatij <- eval.fd(tx, xfdlist[[i1]], j1-1, returnMatrix)
            resmat <- resmat + bmatij*xmatij
          }
        }
      }
      #  set up the functional data object
      resfdi            <- smooth.basis(tx, resmat, resbasis)$fd
      resfdnames        <- xfdi$fdnames
      resfdnames[[2]]   <- "Residual function"
      resfdnames[[3]]   <- "Residual function value"
      resfdlist[[ivar]] <- resfdi
    }

    #  ----------------------------------------------------------------
    #                   End of multiple variable case
    #  ----------------------------------------------------------------

  }

  pdaList <- list(bwtlist=bwtlist, resfdlist=resfdlist, awtlist=awtlist)
  class(pdaList) <- 'pda.fd'
  pdaList
}

# ---------------------------------------------------------------------------

trapzmat <- function(X,Y,delta=1,wt=rep(1,n)) {
#TRAPZMAT integrates the products of two matrices of values
#   using the trapezoidal rule, assuming equal spacing
#  X is the first  matrix of values
#  Y is the second matrix of values
#  DELTA is the spacing between argument values (one by default)
#  WT is a vector of weights (ones by default)
#
#  XtWY is a matrix of integral estimates, number of rows equal to
#  number of col of X, number of cols equal to number of cols of Y

	X <- as.matrix(X)
	Y <- as.matrix(Y)
	
	n <- dim(X)[1]

	if (dim(Y)[1] != n) {
    	stop("X and Y do not have same number of rows.")
	}

	if (length(wt) != n) {
    	stop("X and WT do not have same number of rows.")
	}

	if (delta <= 0) {
    	stop("DELTA is not a positive value.")
	}

	wt[c(1,n)] <- wt[c(1,n)]/2
	wt <- wt*delta

	X <- X*outer(wt,rep(1,dim(X)[2]))
	XtWY <- crossprod(X,Y)
	return(XtWY)
}

