# CHNOSZ/wjd.R  
# Gibbs energy minimization and supporting functions
# 20111117 jmd

wjd <- function(
  # example problem definition
  # the formula matrix: composition of the species with elements H, N, O
  A = matrix(c(
    1,2,2,0,0,1,0,0,0,1,
    0,0,0,1,2,1,1,0,0,0,
    0,0,1,0,0,0,1,1,2,1),ncol=3,
    dimnames=list(NULL,c("H","N","O"))),
  # the free energies, G0/RT, at 3500 K
  G0.RT = c(
    -10.021,-21.096,-37.986,-9.846,-28.653,
    -18.918,-28.032,-14.640,-30.594,-26.111),
  # initial solution, a positive set of values (numbers of moles)
  Y = c(0.1,0.35,0.5,0.1,0.35,0.1,0.1,0.1,0.1,0.1),
  # pressure in atmospheres
  P = 51,
  # number of values of lambda tested at each iteration
  nlambda = 101,
  # maximum number of iterations
  imax = 10,
  # the free energy change, as a fraction of the total system energy
  # in the current step, below which the iterations will stop
  Gfrac = 1e-7
) {
  # Gibbs energy minimization
  # using steepest descent approach of
  # White, Johnson and Dantzig, 1958
  # J. Chem. Phys. 28(5), 751-755
  # doi:10.1063/1.1744264

  # trying to make it faster than version coded 2005-02-10
  # notation notes: i for species, j for elements
  # use G0 here instead of F0 to stand for standard Gibbs free energy

  ## most common incorrect call is to have non-positive starting mole numbers
  if(any(Y <= 0)) stop("all mole numbers of initial solution (Y) must be positive")

  ## begin function definitions
  # Eq. 18 - coefficients in the system of unknowns
  # A: formula matrix
  # Y: set of mole numbers
  coeff <- function(A,Y) {
    # m - number of elements (columns)
    m <- ncol(A)
    # n - number of species (rows)
    n <- nrow(A)
    # our matrix of r's is m x m
    R <- matrix(0,nrow=m,ncol=m)
    for(j in 1:m) {
      for(k in j:m) {
        # Eq. 17 - r coefficients 
        r <- sum(A[,j] * A[,k] * Y)
        # fill in the matrix
        R[j,k] <- r
        R[k,j] <- r
      }
    }
    # Eq. 4 - total number of (atomic weights?) of element j
    B <- Y %*% A
    coeff <- cbind(R,t(B))
    B <- c(B,0)
    coeff <- rbind(coeff,B)
    return(coeff)
  }
  # Eq. 15 - free energies of the species
  f.Y <- function(Y,C) return(Y * (C + log(Y/sum(Y))))
  # Eq. 18 - rhs of matrix equation
  rhs <- function(A,f.Y) {
    # m - number of elements (columns)
    m <- ncol(A)
    # assemble the sum for each element
    rhs <- sapply(1:m, function(j) {
      sum(A[,j] * f.Y)
    })
    # final value is sum over species
    rhs <- c(rhs,sum(f.Y))
    return(rhs)
  }
  # Eq. 14 - the resulting mole number vector
  X <- function(A,Y,f.Y,mults) {
    # m - number of elements (columns)
    m <- ncol(A)
    # n - number of species (rows)
    n <- nrow(A)
    # third term: the summation over elements, for each species
    X3 <- sapply(1:n, function(i) {
      sum(mults[1:m] * A[i,]) * Y[i]
    })
    # second term: ratio of mole vectors
    # (cf. Eq. 19)
    X2 <- Y * (tail(mults,1) + 1)
    # first term: negative of f.Y
    X1 <- -f.Y
    # put them together
    return(X1 + X2 + X3)
  }
  ## end function definitions

  # now set up up the calculation
  # keep the initial solution around
  Y.0 <- Y
  # following Eq. 2
  # G0.RT: standard Gibbs free energies
  # P: pressure in atmospheres
  C <- G0.RT + log(P)
  # initialize iteration counter
  i <- 0
  # initialize system free energy output
  F.Y <- numeric()
  # initialize fractional distance change output
  lambda <- numeric()
  # initialize free energy change output
  Ffrac <- numeric()
  # initialize mass balance results
  elements <- t(A) %*% Y

  # we iterate 
  repeat {
    # determine f.Y by Eq. 15
    f.Y.1 <- f.Y(Y,C)
    # compute system free energy
    F.Y <- c(F.Y, sum(f.Y.1))
    # don't surpass the maximum number of iterations
    i <- i + 1
    if(i > imax) break
    # stop if the last iteration changed the free energy by less than required
    if(i > 1) {
      d.F.Y <- diff(tail(F.Y,2))
      Ffrac <- c(Ffrac, abs(d.F.Y / tail(F.Y,1)))
      if(tail(Ffrac,1) < Gfrac) break
    }
    # set up the system of equations
    coeff.1 <- coeff(A,Y)
    rhs.1 <- rhs(A,f.Y.1)
    # solve the system to get the Lagrange multipliers
    mults <- solve(coeff.1,rhs.1)
    # get the new (possibly negative) mole values
    X.1 <- X(A,Y,f.Y.1,mults)
    # what are the directional changes
    D <- X.1 - Y
    # lambda is the fractional amount we go along that direction;
    # WJD58 give two constraints but no specific exploration procedure.
    # first constraint is that all mole numbers are positive. 
    if(any(X.1 < 0)) {
      # find the lowest value of lambda where any mole
      # number becomes zero (the species is zapped, not allowed!)
      lam <- -Y/D
      lam[lam < 0] <- 1
      lamzap <- min(lam)
      lastlam <- nlambda - 1
    } else {
      # all mole numbers are positive; we can take it all the way
      lamzap <- 1
      lastlam <- nlambda
    }
    # let's explore lambda between 0 and lamzap
    # (including lamzap if it's 1)
    lams <- seq(0,lamzap,length.out=nlambda)[1:lastlam]
    # second constraint is that the derivative of free energy 
    # doesn't go positive ... Eq. 20
    d.f.Y <- function(lambda,Y,D,C) {
      d.f.Y <- sum(D * (C + log( (Y + lambda * D) / (sum(Y) + lambda * sum(D)) )))
      return(d.f.Y)
    }
    # what are the free energy derivatives
    d.f.Y.1 <- sapply(lams,d.f.Y,Y=Y,D=D,C=C)
    # if any are positive, exclude those lambdas
    lams[d.f.Y.1 > 0] <- 0
    # take the highest lambda
    lambda <- c(lambda,max(lams))
    # we now have lambda, so we can calculate a new value for Y
    Y <- Y + tail(lambda,1) * D
    # it might be wise to check the mass balance
    elements <- cbind(elements,t(A) %*% Y)
    # next iteration
  }

  # the result is in 'X' to be consistent with notation of WJD58
  X <- Y
  # the initial solution is in 'Y'
  Y <- Y.0
  # return problem definition and results
  return(list(A=A,G0.RT=G0.RT,Y=Y,P=P,X=X,F.Y=F.Y,lambda=lambda,Ffrac=Ffrac,elements=elements))
}

element.potentials <- function(w, plot.it=FALSE, iplot=1:ncol(w$A)) {
  # calculate the chemical potentials of the elements 
  # from the output of wjd(), using all or some of the combinations
  # of species that are compositionally independent
  # 20111126 jmd
  # put the species in abundance order
  oX <- order(w$X)
  # the mole fractions, formulas, and energies of the species in this order
  X <- w$X[oX]
  A <- w$A[oX,]
  G0.RT <- w$G0.RT[oX] + log(w$P)
  # get the combinations of species that are compositionally independent
  ic <- invertible.combs(A)
  # a function to calculate chemical potentials of the elements for the ith combination of species
  mu <- function(i) {
    myA <- A[ic[i,],]
    # chemical potentials (/RT) of the species: G0/RT + ln(mole fraction)
    myB <- (G0.RT + log(X/sum(X)))[ic[i,]]
    # chemical potentials of the elements
    myX <- solve(myA,myB)
  }
  # run the calculation over all combinations
  ep <- t(sapply(1:nrow(ic),mu))
  # keep names of the elements
  colnames(ep) <- colnames(w$A)
  # to make a plot
  if(plot.it) {
    par(mfrow=c(length(iplot),1))
    for(i in iplot) {
      ylab <- as.expression(substitute(mu[x]/RT,list(x=colnames(ep)[i])))
      plot(ep[,i],xlab="species combination",ylab=ylab)
      title(main=paste("max difference (range) =",format(diff(range(ep[,i])),digits=2)))
    }
  }
  return(ep)
}

is.near.equil <- function(w, tol=0.01, quiet=FALSE) {
  # given the output of wjd(), make a simple test for equilibrium
  # that is, that the chemical potentials of the elements are nearly
  # the same when calculated using different sets of species in the system
  ep <- element.potentials(w)
  # stop if we don't have at least two combinations
  if(nrow(ep) < 2) stop("can not test for equilibrium because species abundances are determined")
  # equilibrium unless proven guilty
  ine <- TRUE
  for(i in 1:ncol(ep)) if(diff(range(ep[,i])) > tol) ine <- FALSE
  if(!ine & !quiet) {
    # talk about the differences in chemical potentials
    epdiff <- abs(apply(apply(ep, 2, range), 2, diff))
    imax <- which.max(epdiff)
    msgout("is.near.equil: solution has variation of ", epdiff[imax], " in mu/RT of ", names(epdiff)[imax], "\n")
  }
  return(ine)
}

guess <- function(
  A = matrix(c(
    1,2,2,0,0,1,0,0,0,1,
    0,0,0,1,2,1,1,0,0,0,
    0,0,1,0,0,0,1,1,2,1),ncol=3,
    dimnames=list(NULL,c("H","N","O"))),
  B = c(2,1,1), method="stoich", minX=0.001, iguess=1, ic=NULL
){
  # given the elemental stoichiometries of a set of species (A)
  # and the number of moles of elements (B)
  # find moles of species that satisfy mass balance and are all positive
  # generally this will be one of the solutions of an underdetermined system

  # first of all, we can't do anything if all there are no elements
  if(all(B==0)) stop("there are zero moles of all elements")

  # if method="central" get central solution using limSolve package  20120919
  if(identical(method, "central")) {
    if(!"limSolve" %in% row.names(installed.packages())) {
      msgout("guess: skipping 'central' method as limSolve package is not available\n")
    } else {
      # the inequality constraints for moles of species
      G <- diag(nrow(A))
      # minX is the minimum mole number we will accept
      H <- rep(minX, nrow(A))
      # get a solution
      X <- limSolve::xranges(E=t(A), F=B, G=G, H=H, central=TRUE, full=TRUE)[, "central"]
      return(X)
    }
  }

  if(identical(method, "stoich")) {
    # if method="stoich" use a stoichiometric approach: 20111231 jmd
    # - select one of the (many) species combinations (ic) that
    #   make a square, invertible stoichiometric matrix (the "variable" species)
    # - assign equal mole numbers to all the "other" species (Xother),
    #   such that any element has at most max.frac fraction of the desired amount (B)
    #   (max.frac is scanned from 0.01 to 0.99)
    # - calculate the mole numbers of the stoichiometry-setting species
    #   that give the desired elemental composition; accept the provisional
    #   solution if all numbers are positive

    # arguments:
    # A - the stoichiometric matrix
    # B - the moles of elements
    # iguess - which provisional guess to return (NULL for all)
    # ic - which specific combination of species to test (NULL for all)

    # get the various combinations of species that are
    # stoichiometrically independent
    combs <- invertible.combs(A)
    # we will potentially process all of them unless a specific one is identified
    if(is.null(ic)) ic <- 1:nrow(combs)
    # a counter to keep track of the provisional guesses
    iprov <- 0
    # where to store the output if we want all guesses
    out <- list()
    for(i in ic) {
      # which species are the variable ones
      ivar <- combs[i,]
      # moles of elements for one mole of all of the other species
      Bother.1 <- colSums(A[-ivar, , drop=FALSE])
      # which element is proportionally most highly represented w.r.t. the desired composition
      imax <- which.max(Bother.1/B)
      # max.frac - the highest fraction of contribution to moles of elements by the "other" species
      for(max.frac in c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)) {
        # the number of moles of all "other" species that give max.frac of any element
        Xother <- max.frac/(Bother.1/B)[imax]
        # moles of elements for this number of moles of all the other species
        Bother <- Bother.1 * Xother
        # moles of elements that are left for the variable species
        Bvar <- B - Bother
        # now solve for the number of moles of the variable species
        Xvar <- solve(t(A[ivar,]),Bvar)
        # stop the search if we found a positive solution
        if(all(Xvar > 0)) break
      }
      # put them together
      X <- numeric(nrow(A))
      X[-ivar] <- Xother
      X[ivar] <- Xvar
      # add names
      names(X) <- rownames(A)
      # if all the moles are positive, this is a provisional
      # guess, otherwise make the result NA
      if(any(Xvar <= 0)) X <- NA
      else iprov <- iprov + 1
      # return the result if we're at the correct guess number
      if(is.null(iguess)) out <- c(out,list(X))
      else if(iprov==iguess) return(X)
    }
    # if we're here, we should return all guesses, or 
    # make an error (the requested guess number doesn't exist)
    if(is.null(iguess) & iprov > 0) return(out)
    else {
      if(is.null(iguess)) iguess <- "[ALL]"
      stop(paste("you asked for guess number ",iguess,
        " but there are only ",iprov,
        " that satisfy all stoichiometric constraints",sep=""))
    }
  }

  # if we're here we didn't find a guessing method
  stop("no method found")
}

run.wjd <- function(ispecies, B=NULL, method="stoich", Y=run.guess(ispecies, B, method),
  P=1, T=25, nlambda=101, imax=10, Gfrac=1e-7, tol=0.01) {
  ### set up a Gibbs energy minimization
  ### using compositions and standard Gibbs energies of species
  ### from database in CHNOSZ  20120101 jmd
  ## get the stoichiometric matrix for the species
  A <- i2A(ispecies)
  ## assemble the standard molal Gibbs energies of the species
  s <- subcrt(ispecies, P=P, T=T, property="G", exceed.Ttr=TRUE)
  G0 <- sapply(1:length(s$out), function(i) s$out[[i]]$G)
  G0.RT <- G0/get("thermo")$opt$R/convert(T, "K")
  ## if Y is provided use that as initial guess
  if(!missing(Y)) {
    # giving both Y and B is not allowed
    if(!is.null(B)) stop("Y and B can not both be provided")
    # the length of Y must be equal to number of species
    if(length(Y) != nrow(A)) stop("Y must have same length as number of species")
    # a single guess
    w <- wjd(A, G0.RT, Y, P=P, nlambda=nlambda, imax=imax, Gfrac=Gfrac)
  } else {
    # if we're using method "central" there is only one guess
    if(method=="central") {
      w <- wjd(A, G0.RT, Y, P=P, nlambda=nlambda, imax=imax, Gfrac=Gfrac)
    } else {
      # for method "stoich" loop over all the guesses created by run.guess
      Y <- Y[!is.na(Y)]
      for(i in 1:length(Y)) {
        w <- wjd(A, G0.RT, Y[[i]], P=P, nlambda=nlambda, imax=imax, Gfrac=Gfrac)
        if(is.near.equil(w, tol=tol)) {
          msgout("run.wjd: got within tolerance on initial solution ", i, " of ", length(Y), "\n")
          break
        }
        if(i==length(Y)) msgout("run.wjd: tested ", length(Y), " initial solutions\n")
      }
    }
    # only return a near equilibrium solution
    if(!is.near.equil(w, tol=tol)) {
      stop(paste("couldn't find a solution within mu/RT tolerance of", tol))
    }
  }
  return(w)
}

run.guess <- function(ispecies, B=NULL, method="stoich", iguess=NULL) {
  ## run guess() using species from database  20120612
  # get the stoichiometric matrix for the species
  A <- i2A(ispecies)
  # we need B
  if(is.null(B)) stop(paste("please provide B (numbers of moles of ", paste(colnames(A), collapse=", ")  , ")", sep=""))
  # if B is a formula, turn it into a vector
  if(is.character(B)) {
    # zero formula with elements in same order as A
    zero <- paste(colnames(A), "0", collapse="", sep="")
    # sum of zero and B
    mB <- makeup(c(zero, B), sum=TRUE)
    # then turn it into a vector
    B <- as.numeric(unlist(mB))
  } else B <- as.vector(B)
  # setup initial guess
  Y <- guess(A, B, method=method, iguess=iguess)
  # take away NA guesses
  Y <- Y[!is.na(Y)]
  return(Y)
}

equil.potentials <- function(w, tol=0.01, T=25) {
  ## return the average of the element.potentials, only if w is.near.equil  20120613
  if(!is.near.equil(w, tol=tol)) return(NULL)
  else return(colMeans(element.potentials(w)) * get("thermo")$opt$R * convert(T, "K"))
}
