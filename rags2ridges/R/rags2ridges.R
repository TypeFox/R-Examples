################################################################################
################################################################################
################################################################################
##
## Name: rags2ridges
## Authors: Carel F.W. Peeters, Anders E. Bilgrau, and Wessel N. van Wieringen
##
## Maintainer: Carel F.W. Peeters
##             Statistics for Omics Research Unit
##             Dept. of Epidemiology & Biostatistics
##             VU University medical center
##             Amsterdam, the Netherlands
## Email:	     cf.peeters@vumc.nl
##
## Version: 2.0
## Last Update:	12/10/2015
## Description:	Ridge estimation for high-dimensional precision matrices
##              Includes supporting functions for (integrative) graphical modeling
##
## Code files: rags2ridges.R      >> master file/core module
##             rags2ridgesFused.R >> fused module
##             rags2ridgesMisc.R  >> miscellaneous module
##             rags2ridges.cpp    >> C++ work horses
##
## Publications:
##   [1] van Wieringen, W.N. & Peeters, C.F.W. (2015).
##       "Ridge Estimation of Inverse Covariance Matrices from High-Dimensional
##       Data", arXiv:1403.0904v3 [stat.ME].
##   [2] van Wieringen, W.N. & Peeters, C.F.W. (2015).
##       "Application of a New Ridge Estimator of the Inverse Covariance Matrix
##       to the Reconstruction of Gene-Gene Interaction Networks", in: di Serio,
##       C., Lio, P., Nonis, A., and Tagliaferri, R. (Eds.) Computational
##       Intelligence Methods for Bioinformatics and Biostatistics. Lecture
##       Notes in Computer Science, vol. 8623. Springer, pp. 170-179.
##   [3] Bilgrau, A.E., Peeters, C.F.W., Eriksen, P.S., Boegsted, M., &
##       van Wieringen, W.N. (2015).
##       "Targeted Fused Ridge Estimation of Inverse Covariance Matrices from
##       Multiple High-Dimensional Data Classes", arXiv:1509.07982v1 [stat.ME].
## 	 [4] Peeters, C.F.W., van de Wiel, M.A., & van Wieringen, W.N. (2015).
##       "The Spectral Condition Number Plot for Regularization Parameter
##       Determination".
##
################################################################################
################################################################################
################################################################################


################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module A: rags2ridges Core
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

##------------------------------------------------------------------------------
##
## Hidden support functions
##
##------------------------------------------------------------------------------

.trace <- function(M){
  ##############################################################################
  # - Internal function to compute the trace of a matrix
  # - Faster support function (as opposed to 'matrix.trace') when input M is
  #   already forced to 'matrix'
  # - M > matrix input
  ##############################################################################

  return(sum(diag(M)))
}



.is.int <- function(x, tolerance = .Machine$double.eps){
  ##############################################################################
  # - Logical function that checks if a number is an integer within machine
  #   precision
  # - x         > input number
  # - tolerance > tolerance threshold for determining integer quality
  ##############################################################################

  #return(all.equal(x, as.integer(x), ...)) # Faster(?), safer
  return(abs(x - round(x)) < tolerance)
}



.LL <- function(S, P){
  ##############################################################################
  # - Function that computes the value of the (negative) log-likelihood
  # - S > sample covariance matrix
  # - P > precision matrix (possibly regularized inverse of covariance or
  #       correlation matrix)
  ##############################################################################

  LL <- -log(det(P)) + sum(S*P) #.trace(S %*% P)
  return(LL)
}



.Frobenius <- function(X) {
  ##############################################################################
  # - Function computing squared Frobenius norm - the sum of the squarred
  #   entries.
  # - X > A numeric
  ##############################################################################

  return(sum(X^2))
}



.FrobeniusLoss <- function(O, P){
  ##############################################################################
  # - Function computing Frobenius loss
  # - O > Estimated (possibly regularized) precision matrix
  # - P > True (population) precision matrix
  ##############################################################################

  return(.Frobenius(O - P))
}



.RelativeFrobeniusLoss <- function(O, P) {
  ##############################################################################
  # - Function computing the frobenius loss relative to the (frobenious) norm
  #   of P.
  # - O > Estimated (possibly regularized) precision matrix
  # - P > True (population) precision matrix
  ##############################################################################

  return(.FrobeniusLoss(O, P)/.Frobenius(P))
}



.QuadraticLoss <- function(O, C){
  ##############################################################################
  # - Function computing Quadratic loss
  # - O > Estimated (possibly regularized) precision matrix
  # - C > True (population) covariance matrix
  ##############################################################################

  return((sum(((O %*% C - diag(ncol(O))))^2)))
}



.eigShrink <- function(dVec, lambda, const = 0){
  ##############################################################################
  # - Function that shrinks the eigenvalues in an eigenvector
  # - Shrinkage is that of rotation equivariant alternative ridge estimator
  # - Main use is in avoiding expensive matrix square root when choosing a
  #   target that leads to a rotation equivariant version of the alternative
  #   ridge estimator
  # - dVec   > numeric vector containing the eigenvalues of a matrix S
  # - lambda > penalty parameter
  # - const  > a constant, default = 0
  ##############################################################################

  Evector <- (dVec - lambda * const)
  return(sqrt(lambda + Evector^2/4) + Evector/2)
}



.ridgeSi <- function(S, lambda, type = "Alt", target = default.target(S)){
  ##############################################################################
  # - Hidden function that calculates Ridge estimators of a covariance matrix
  # - Function is mirror image main routine 'ridgeS'
  # - Main use is to circumvent (unnecessary) inversion (especially in
  #   'conditionNumberPlot' function)
  # - S       > sample covariance matrix
  # - lambda  > penalty parameter (choose in accordance with type of Ridge
  #             estimator)
  # - type    > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - Alt     > van Wieringen-Peeters alternative ridge estimator of a
  #             covariance matrix
  # - ArchI   > Archetypal I ridge estimator of a covariance matrix
  # - ArchII  > Archetypal II ridge estimator of a covariance matrix
  # - target  > target (precision terms) for Type I estimators,
  #             default = default.target(S)
  #
  # - NOTES:
  # - When type = "Alt" and target is p.d., one obtains the
  #   van Wieringen-Peeters type I estimator
  # - When type = "Alt" and target is null-matrix, one obtains the
  #   van Wieringen-Peeters type II estimator
  # - When target is not the null-matrix it is expected to be p.d. for the vWP
  #   type I estimator
  # - The target is always expected to be p.d. in case of the archetypal I
  #   estimator
  # - When type = "Alt" and target is null matrix or of form c * diag(p), a
  #   rotation equivariant estimator ensues. In these cases the expensive
  #   matrix square root can be circumvented
  ##############################################################################

  # Alternative estimator
  if (type == "Alt"){
    if (all(target == 0)){
      Spectral  <- eigen(S, symmetric = TRUE)
      Eigshrink <- .eigShrink(Spectral$values, lambda)
      C_Alt     <- Spectral$vectors %*% diag(Eigshrink) %*% t(Spectral$vectors)
      colnames(C_Alt) = rownames(C_Alt) <- colnames(S)
      return(C_Alt)
    } else if (all(target[!diag(nrow(target))] == 0) &
                 (length(unique(diag(target))) == 1)){
      varPhi    <- unique(diag(target))
      Spectral  <- eigen(S, symmetric = TRUE)
      Eigshrink <- .eigShrink(Spectral$values, lambda, const = varPhi)
      C_Alt     <- Spectral$vectors %*% diag(Eigshrink) %*% t(Spectral$vectors)
      colnames(C_Alt) = rownames(C_Alt) <- colnames(S)
      return(C_Alt)
    } else {
      D     <- (S - lambda * target)
      C_Alt <- D/2 + sqrtm((D %*% D)/4 + lambda * diag(nrow(S)))
      return(C_Alt)
    }
  }

  # Archetypal I
  if (type == "ArchI"){
    C_ArchI <- (1-lambda) * S + lambda * solve(target)
    return(C_ArchI)
  }

  # Archetypal II
  if (type == "ArchII"){
    C_ArchII <- S + lambda * diag(nrow(S))
    return(C_ArchII)
  }
}



.pathContribution <- function(sparseP, path, detSparseP){
  ##############################################################################
  # - Function calculating the contribution of a path to the covariance between
  #   begin and end node
  # - sparseP    > sparse precision/partial correlation matrix
  # - path       > path between two nodes (start and end node)
  # - detSparseP > determinant of 'sparseP'
  ##############################################################################

  if (length(path) < nrow(sparseP)){
    return((-1)^(1+length(path)) *
             prod(sparseP[cbind(path[-1], path[-length(path)])]) *
             det(sparseP[-path, -path]) / detSparseP)
  }
  if (length(path) == nrow(sparseP)){
    return((-1)^(1+length(path)) *
             prod(sparseP[cbind(path[-1], path[-length(path)])]) / detSparseP)
  }
}



.path2string <- function(path){
  ##############################################################################
  # - Function converting a numeric or character vector into a single string
  # - path > path between two nodes (start and end node)
  ##############################################################################

  pName <- sprintf("%s--%s", path[1], path[2])
  if (length(path) > 2){
    for (w in 3:(length(path))){ pName <- sprintf("%s--%s", pName, path[w]) }
  }
  return(pName)
}



.pathAndStats <- function(Gt, node1t, node2t, nei1t, nei2t, P0t, detP0t,
                          pathNames){
  ##############################################################################
  # - Function determining shortest paths (and their contribution) between
  #   node 1 and node 2
  # - It does so via the neighborhoods 1 and 2
  # - Gt        > graphical object
  # - node1t    > start node of the path
  # - node2t    > end node of the path
  # - nei1t     > neighborhood around the start node
  # - nei2t     > neighborhood around the end node
  # - p0t       > sparse precision/partial correlation matrix
  # - detP0t    > determinant of 'p0t'
  # - pathNames > named path represented as string
  ##############################################################################

  pathsTemp <- list()
  pathStatsTemp <- numeric()
  for (v1 in 1:length(nei1t)){
    for (v2 in 1:length(nei2t)){
      slhNo1Ne1 <- get.all.shortest.paths(Gt, node1t, nei1t[v1])$res
      slhNe1Ne2 <- get.all.shortest.paths(Gt, nei1t[v1], nei2t[v2])$res
      slhNe2No2 <- get.all.shortest.paths(Gt, nei2t[v2], node2t)$res
      for (uNo1Ne1 in 1:length(slhNo1Ne1)){
        for (uNe1Ne2 in 1:length(slhNe1Ne2)){
          for (uNe2No2 in 1:length(slhNe2No2)){
            fullPath <- c(slhNo1Ne1[[uNo1Ne1]], slhNe1Ne2[[uNe1Ne2]][-1],
                          slhNe2No2[[uNe2No2]][-1])
            if (length(unique(fullPath)) == (length(fullPath))){
              pName <- .path2string(fullPath)
              if (!(pName %in% c(pathNames, rownames(pathStatsTemp)))){
                pathsTemp[[length(pathsTemp)+1]] <- fullPath
                pathStatsTemp <-
                  rbind(pathStatsTemp,
                        c(length(fullPath)-1,
                          .pathContribution(P0t, fullPath, detP0t)))
                rownames(pathStatsTemp)[nrow(pathStatsTemp)] <- pName
              }
            }
          }
        }
      }
    }
  }
  return(list(paths=pathsTemp, pathStats=pathStatsTemp))
}



.cvl <- function(lambda, Y,  target = default.target(covML(Y)), type = "Alt"){
  ##############################################################################
  # - Function that calculates a cross-validated negative log-likelihood score
  #   for single penalty value
  # - lambda > value penalty parameter
  # - Y      > (raw) Data matrix, variables in columns
  # - target > target (precision terms) for Type I estimators,
  #            default = default.target(covML(Y))
  # - type   > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  ##############################################################################

  slh <- numeric()
  for (i in 1:nrow(Y)){
    S   <- covML(Y[-i, ])
    slh <- c(slh, .LL(t(Y[i, , drop = FALSE]) %*% Y[i, , drop = FALSE],
                      ridgeP(S, lambda, target = target, type = type)))
  }
  return(mean(slh))
}



.lambdaNullDist <- function(i, Y, id, lambdaMin, lambdaMax,
                            lambdaInit, target, type){
  ##############################################################################
  # - Function that determines the optimal value of the penalty parameter for a
  #   single permutation
  # - Optimal penalty determined using the 'optPenalty.LOOCVauto' function
  # - i          > number of permutations; passed by nPerm in
  #                'GGMblockNullPenalty'
  # - Y          > (raw) Data matrix, variables in columns
  # - id         > indicator variable for the two blocks of the precision matrix
  # - lambdaMin  > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax  > maximum value penalty parameter (dependent on 'type')
  # - lambdaInit > initial value for lambda for starting optimization
  # - target     > target (precision terms) for Type I estimators
  # - type       > must be one of {"Alt", "ArchI", "ArchII"}
  ##############################################################################

  reshuffle    <- sample(1:nrow(Y), nrow(Y))
  Y[, id == 1] <- Y[reshuffle, id == 1]
  return(optPenalty.LOOCVauto(Y, lambdaMin, lambdaMax, lambdaInit,
                              target = target, type = type)$optLambda)
}



.blockTestStat <- function(i, Y, id, lambda, target, type){
  ##############################################################################
  # - Function that calculates logratio statistic for block independence
  # - i      > number of permutations; passed by nPerm in 'GGMblockTest'
  # - Y      > (raw) Data matrix, variables in columns
  # - id     > indicator variable for the two blocks of the precision matrix
  # - lambda > value penalty parameter
  # - target > target (precision terms) for Type I estimators
  # - type   > must be one of {"Alt", "ArchI", "ArchII"}
  ##############################################################################

  reshuffle    <- sample(1:nrow(Y), nrow(Y))
  Y[, id == 1] <- Y[reshuffle, id == 1]
  S <- solve(ridgeP(covML(Y), lambda = lambda, target = target, type = type))
  return(log(det(S[id == 0, id == 0])) +
           log(det(S[id == 1, id == 1])) - log(det(S)))
}




##------------------------------------------------------------------------------
##
## Support functions
##
##------------------------------------------------------------------------------

symm <- function(M){
  ##############################################################################
  # - Large objects that are symmetric sometimes fail to be recognized as such
  #   by R due to rounding under machine precision. This function symmetrizes
  #   for computational purposes matrices that are symmetric in numeric ideality
  # - M > symmetric (in numeric ideality) square matrix
  ##############################################################################

  # Dependencies
  # require("base")

  if (!is.matrix(M)){
    stop("M should be a matrix")
  }
  else if (nrow(M) != ncol(M)){
    stop("M should be a square matrix")
  }
  else {
    # Symmetrize
    Msym <- (M + t(M))/2

    # Return
    return(Msym)
  }
}



adjacentMat <- function(M, diag = FALSE){
  ##############################################################################
  # - Function that transforms a real matrix into an adjacency matrix
  # - Intended use: Turn sparsified precision matrix into an adjacency matrix
  #   for undirected graph
  # - M    > (sparsified precision) matrix
  # - diag > logical indicating if the diagonal elements should be retained
  ##############################################################################

  # Dependencies
  # require("base")

  if (!is.matrix(M)){
    stop("M should be a matrix")
  }
  else if (nrow(M) != ncol(M)){
    stop("M should be square matrix")
  }
  else {
    # Create adjacency matrix
    AM <- M
    AM[AM != 0] <- 1
    diag(AM) <- 0

    if (diag){
      diag(AM) <- 1
    }

    # Return
    return(AM)
  }
}



covML <- function(Y){
  ##############################################################################
  # - function that gives the maximum likelihood estimate of the covariance
  # - Y   > (raw) data matrix, assumed to have variables in columns
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")

  if (!is.matrix(Y)){
    stop("Input (Y) should be a matrix")
  }
  else {
    Ys  <- scale(Y, center = TRUE, scale = FALSE)
    Sml <- crossprod(Ys)/nrow(Ys)  # (t(Ys) %*% Ys)/nrow(Ys)
    return(Sml)
  }
}



covMLknown <- function(Y, covMat = NULL, corMat = NULL,
                       corType = "none", varType = "none", nInit = 100){
  ##############################################################################
  # - Maximum likelihood estimation of the covariance matrix, with various
  #   types of assumptions on the structure of this matrix.
  # - Y	      > (raw) data matrix, assumed to have variables in columns
  # - covMat  > A positive-definite covariance 'matrix'. When specified, the
  #             to-be-estimated covariance matrix is assumed to be proportional
  #             to the specified covariance matrix. Hence, only a constant needs
  #             to be estimated
  # - corMat  > A positive-definite correlation 'matrix'. When specified, the
  #             to-be-estimated covariance matrix is assumed to have this
  #             correlation structure. Hence, only the variances need to
  #		          be estimated.
  # - corType > character that is either "none" (no structure on the correlation
  #             among variate assumed) or "equi" (variates are equi-correlated).
  # - varType	> character that is either "none" (no structure on the marginal
  #             variances of the variates assumed) or "common" (variates have
  #             equal marginal variances).
  # - nInit   > Maximum number of iterations for likelihood maximization
  #             when corType='equi'.
  #
  # NOTES:
  # - Currently no dependencies. In the future chunks of code may be C-ed
  #   through Rcpp.
  # - Future version should a.o. also allow a first order autoregressive
  #   correlation assumption.
  ##############################################################################

  # center data
  Ys <- scale(Y, center = TRUE, scale = FALSE)

  # equicorrelated covariance: parameter estimation
  if (corType=="equi" & varType !="none"){

    # initial variance estimate
    sds <- sqrt(apply(Ys^2, 2, mean))
    rho <- 0
    llNew <- 10^(-10)
    for (k in 1:nInit){
      # update current log-likelihood
      llPrev <- llNew

      # estimate rho
      Sml <- diag(1/sds) %*% ((t(Ys) %*% Ys)/nrow(Ys)) %*% diag(1/sds)
      a1 <- sum(diag(Sml))
      a2 <- (ncol(Ys)-1) * a1 - sum(Sml)
      p <- nrow(Sml)
      minLoglikEqui <- function(rho, p, a1, a2){
        (p-1) * log(1-rho) +  log((p-1) * rho + 1) + ((1-rho) * ((p-1) * rho + 1))^(-1) * (a1 + a2 * rho) }
      rho <- optim(par=0.1, fn=minLoglikEqui, method="Brent", lower=-1/(p-1)
                   + .Machine$double.eps, upper=1-.Machine$double.eps, p=ncol(Ys),
                   a1=a1, a2=a2)$par
      llNew <- minLoglikEqui(rho, p, a1, a2)
      if (abs(llNew - llPrev) < 0.0001){ break }

      # estimate variance(s)
      Sml <- matrix(rho, ncol(Ys), ncol(Ys))
      diag(Sml) <- 1
      if (varType!="common"){
        V <- eigen(Sml)
        D <- abs(V$values)
        V <- V$vectors
        Vinner <- eigen(diag(1/sqrt(D)) %*% t(V) %*% ((t(Ys) %*% Ys)/nrow(Ys))
                        %*% V %*% diag(1/sqrt(D)))
        Vinner$values <- sqrt(abs(Re(Vinner$values)))
        sds <- Re(diag(V %*% diag(sqrt(D)) %*% Vinner$vectors
                       %*% diag(Vinner$values) %*% t(Vinner$vectors)
                       %*% diag(sqrt(D)) %*% t(V)))
      }
      if (varType=="common"){
        sds <- rep(sum(diag(solve(Sml) %*% (t(Ys) %*% Ys)/nrow(Ys))) /
                     ncol(Ys), ncol(Ys))
      }
    }

    # equicorrelated covariance: matrix construction
    Sml <- matrix(rho, ncol(Ys), ncol(Ys))
    diag(Sml) <- 1
    Sml <- diag(sds) %*% Sml %*% diag(sds)
    return(Sml)
  }
  if (!is.null(covMat)){
    # covariance known up to a constant
    c <- sum(diag(solve(covMat) %*% (t(Ys) %*% Ys)/nrow(Ys))) / ncol(Ys)
    Sml <- c * covMat
    return(Sml)
  }
  if (!is.null(corMat)){
    # correlation matrix known, variances unknown
    V <- eigen(corMat)
    D <- abs(V$values)
    V <- V$vectors
    Vinner <- eigen(diag(1/sqrt(D)) %*% t(V) %*% ((t(Ys) %*% Ys)/nrow(Ys))
                    %*% V %*% diag(1/sqrt(D)))
    Vinner$values <- sqrt(abs(Re(Vinner$values)))
    sds <- Re(diag(V %*% diag(sqrt(D)) %*% Vinner$vectors
                   %*% diag(Vinner$values) %*% t(Vinner$vectors)
                   %*% diag(sqrt(D)) %*% t(V)))
    Sml <- diag(sds) %*% corMat %*% diag(sds)
    return(Sml)
  }
  if (!is.null(covMat) & !is.null(corMat) & corType=="none" & varType=="none"){
    return(covML(Ys))
  }
}



evaluateS <- function(S, verbose = TRUE){
  ##############################################################################
  # - Function evualuating various properties of an input matrix
  # - Intended use is to evaluate the various properties of what is assumed to
  #   be a covariance matrix
  # - Another use is to evaluate the various properties of a (regularized)
  #   precision matrix
  # - S       > sample covariance/correlation matrix or (regularized) precision
  #             matrix
  # - verbose > logical indicating if output should be printed on screen
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")

  if (!is.matrix(S)){
    stop("S should be a matrix")
  }
  else if (nrow(S) != ncol(S)){
    stop("S should be a square matrix")
  }
  else if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  else {
    Sproperties <- list()

    # Is S symmetric?
    Sproperties$symm <- isSymmetric(S)

    # Are eigenvalues S real and positive?
    evs                   <- eigen(S)$values
    Sproperties$realEigen <- all(Im(evs) == 0)
    Sproperties$posEigen  <- all(evs > 0)

    # Is S diagonally dominant?
    Sproperties$diagDom <- all(abs(cov2cor(S)[upper.tri(S)]) < 1)

    # Trace and determinant S
    Sproperties$trace <- sum(diag(S))
    Sproperties$det   <- det(S)

    # Spectral condition number S
    Sproperties$condNumber <- abs(max(evs) / min(evs))

    if (verbose){
      cat("Properties of input matrix:\n")
      cat("----------------------------------------\n")
      cat("       symmetric : ", Sproperties$symm, "\n", sep="")
      cat("eigenvalues real : ", Sproperties$realEigen, "\n", sep="")
      cat(" eigenvalues > 0 : ", Sproperties$posEigen, "\n", sep="")
      cat("  diag. dominant : ", Sproperties$diagDom, "\n\n", sep="")
      cat("           trace : ", round(Sproperties$trace, 5), "\n", sep="")
      cat("     determinant : ", round(Sproperties$det, 5), "\n", sep="")
      cat(" l2 cond. number : ", round(Sproperties$condNumber, 5), "\n", sep="")
      cat("----------------------------------------\n")
    }

    # Return
    return(Sproperties)
  }
}



pcor <- function(P, pc = TRUE){
  ##############################################################################
  # - Function computing partial correlation/standardized precision matrix from
  #   a precision matrix
  # - P  > precision matrix (possibly regularized inverse of covariance or
  #        correlation matrix)
  # - pc > logical indicating if the partial correlation matrix should be
  #        computed
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")

  if (!is.matrix(P)){
    stop("P should be a matrix")
  }
  else if (!isSymmetric(P)){
    stop("P should be a symmetric matrix")
  }
  else {
    # Compute partial correlation matrix
    if (pc){
      P       <- -P
      diag(P) <- -diag(P)
      Pcor    <- cov2cor(P)
      return(Pcor)
    }

    # Compute standardized precision matrix
    else {
      SP <- cov2cor(P)
      return(SP)
    }
  }
}



default.target <- function(S, type = "DAIE", fraction = 1e-04, const){
  ##############################################################################
  # - Function that generates a (data-driven) default target for usage in
  #   ridge-type shrinkage estimation
  # - The target that is generated is to be understood in precision terms
  # - See function 'ridgeS'
  # S        > sample covariance/correlation matrix
  # type     > character determining the type of default target;
  #            default = "DAIE" (see notes below)
  # fraction > fraction of largest eigenvalue below which an eigenvalue is
  #            considered zero. Only when type = "DAIE".
  # const    > numeric constant that represents the partial variance.
  #            Only when type = "DCPV"
  #
  # Notes:
  # - The van Wieringen-Peeters type I estimator and the archetypal I estimator
  #   utilize a p.d. target
  # - DAIE: diagonal average inverse eigenvalue
  #   Diagonal matrix with average of inverse nonzero eigenvalues of S as
  #   entries
  # - DIAES: diagonal inverse average eigenvalue S
  #   Diagonal matrix with inverse of average of eigenvalues of S as entries
  # - DUPV: diagonal unit partial variance
  #   Diagonal matrix with unit partial variance as entries (identity matrix)
  # - DAPV: diagonal average partial variance
  #   Diagonal matrix with average of inverse variances of S as entries
  # - DCPV: diagonal constant partial variance
  #   Diagonal matrix with constant partial variance as entries. Allows one to
  #   use other constant than [DAIE, DUPV, DAPV, and in a sense Null]
  # - DEPV: diagonal empirical partial variance
  #   Diagonal matrix with the inverse of variances of S as entries
  # - Null: Null matrix
  #   Matrix with only zero entries
  # - All but DEPV and Null lead to rotation equivariant alternative and
  #   archetype I ridge estimators
  # - Null also leads to a rotation equivariant alternative Type II estimator
  ##############################################################################

  # Dependencies
  # require("base")

  if (!is.matrix(S)){
    stop("Input (S) should be a matrix")
  }
  else if (!isSymmetric(S)){
    stop("Input (S) should be a symmetric matrix")
  }
  else if (class(type) != "character"){
    stop("Input (type) is of wrong class")
  }
  else if (!(type %in% c("DAIE", "DIAES", "DUPV", "DAPV",
                         "DCPV", "DEPV", "Null"))){
    stop("Input (type) should be one of {'DAIE', 'DIAES', 'DUPV', 'DAPV', ",
         "'DCPV', 'DEPV', 'Null'}")
  }
  else {
    # Compute and return a default target matrix
    # Diagonal matrix with average of inverse nonzero eigenvalues of S as
    # entries
    if (type == "DAIE"){
      if (class(fraction) != "numeric"){
        stop("Input (fraction) is of wrong class")
      } else if (length(fraction) != 1){
        stop("Length input (fraction) must be one")
      } else if (fraction < 0 | fraction > 1){
        stop("Input (fraction) is expected to be in the interval [0,1]")
      } else {
        Eigs   <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
        const  <- mean(1/(Eigs[Eigs >= Eigs[1]*fraction]))
        target <- const * diag(ncol(S))
      }
    }

    # Diagonal matrix with inverse of average of eigenvalues of S as entries
    if (type == "DIAES"){
      Eigs   <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
      const  <- 1/mean(Eigs)
      target <- const * diag(ncol(S))
    }

    # Diagonal matrix with unit partial variance as entries
    if (type == "DUPV"){
      target <- diag(ncol(S))
    }

    # Diagonal matrix with average empirical partial variances as entries
    if (type == "DAPV"){
      apv    <- mean(1/diag(S))
      target <- apv * diag(ncol(S))
    }

    # Diagonal matrix with constant partial variance as entries
    if (type == "DCPV"){
      if (class(const) != "numeric"){
        stop("Input (const) is of wrong class")
      } else if (length(const) != 1){
        stop("Length input (const) must be one")
      } else if (const <= 0 | const > .Machine$double.xmax){
        stop("Input (const) is expected to be in the interval (0, Inf)")
      } else {
        target <- const * diag(ncol(S))
      }
    }

    # Diagonal matrix with empirical partial variances as entries
    if (type == "DEPV"){
      target <- diag(1/diag(S))
    }

    # Null matrix
    if (type == "Null"){
      target <- matrix(0, ncol(S), nrow(S))
    }

    # Return
    colnames(target) = rownames(target) <- rownames(S)
    return(target)
  }
}




##------------------------------------------------------------------------------
##
## Function for Ridge Estimation of the Precision matrix
##
##------------------------------------------------------------------------------

ridgeS <- function(S, lambda, type = "Alt", target = default.target(S)){
  ##############################################################################
  # - Function that calculates Ridge estimators of a precision matrix
  # - S       > sample covariance matrix
  # - lambda  > penalty parameter (choose in accordance with type of Ridge
  #             estimator)
  # - type    > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - Alt     > van Wieringen-Peeters alternative ridge estimator of a precision
  #             matrix
  # - ArchI   > Archetypal I ridge estimator of a precision matrix
  # - ArchII  > Archetypal II ridge estimator of a precision matrix
  # - target  > target (precision terms) for Type I estimators,
  #             default = default.target(S)
  #
  # - NOTES:
  # - When type = "Alt" and target is p.d., one obtains the
  #   van Wieringen-Peeters type I estimator
  # - When type = "Alt" and target is null-matrix, one obtains the
  #   van Wieringen-Peeters type II est.
  # - When target is not the null-matrix it is expected to be p.d. for the
  #   vWP type I estimator
  # - The target is always expected to be p.d. in case of the archetypal I
  #   estimator
  # - When type = "Alt" and target is null matrix or of form c * diag(p), a
  #   rotation equivariant estimator ensues. In these cases the expensive
  #   matrix square root can be circumvented
  ##############################################################################

  # Dependencies
  # require("base")
  # require("expm")

  # Deprecation warning
  warning("This function is deprecated. Please use ridgeP instead.")

  if (!isSymmetric(S)) {
    stop("S should be a symmetric matrix")
  }
  else if (lambda <= 0) {
    stop("lambda should be positive")
  }
  else if (!(type %in% c("Alt", "ArchI", "ArchII"))){
    stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
  }
  else{
    # Calculate Ridge estimator
    # Alternative estimator
    if (type == "Alt"){
      if (!isSymmetric(target)){
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("S and target should be of the same dimension")
      } else if (!all(target == 0) &
                   any(eigen(target, symmetric = TRUE,
                             only.values = TRUE)$values <= 0)){
        stop("When target is not a null-matrix it should be p.d. for this ",
             "type of ridge estimator")
      } else if (all(target == 0)){
        Spectral  <- eigen(S, symmetric = TRUE)
        Eigshrink <- .eigShrink(Spectral$values, lambda)
        P_Alt     <- solve(Spectral$vectors %*% diag(Eigshrink) %*%
                             t(Spectral$vectors))
        colnames(P_Alt) = rownames(P_Alt) <- colnames(S)
        return(P_Alt)
      } else if (all(target[!diag(nrow(target))] == 0) &
                   (length(unique(diag(target))) == 1)){
        varPhi    <- unique(diag(target))
        Spectral  <- eigen(S, symmetric = TRUE)
        Eigshrink <- .eigShrink(Spectral$values, lambda, const = varPhi)
        P_Alt     <- solve(Spectral$vectors %*% diag(Eigshrink) %*%
                             t(Spectral$vectors))
        colnames(P_Alt) = rownames(P_Alt) <- colnames(S)
        return(P_Alt)
      } else {
        E     <- (S - lambda * target)
        P_Alt <- solve(E/2 + sqrtm((E %*% E)/4 + lambda * diag(nrow(S))))
        return(P_Alt)
      }
    }

    # Archetypal I
    if (type == "ArchI"){
      if (lambda > 1){
        stop("lambda should be in (0,1] for this type of Ridge estimator")
      } else if (!isSymmetric(target)){
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("S and target should be of the same dimension")
      } else if (any(eigen(target, symmetric = TRUE,
                           only.values = TRUE)$values <= 0)){
        stop("Target should always be p.d. for this type of ridge estimator")
      } else {
        P_ArchI <- solve((1-lambda) * S + lambda * solve(target))
        return(P_ArchI)
      }
    }

    # Archetypal II
    if (type == "ArchII"){
      P_ArchII <- solve(S + lambda * diag(nrow(S)))
      return(P_ArchII)
    }
  }
}



ridgeP <- function(S, lambda, type = "Alt", target = default.target(S)){
  ##############################################################################
  # - Function that calculates ridge estimators of a precision matrix
  # - Version that uses the RcppArmadillo implementation
  # - S       > sample covariance matrix
  # - lambda  > penalty parameter (choose in accordance with type of Ridge
  #             estimator)
  # - type    > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - Alt     > van Wieringen-Peeters alternative ridge estimator of a precision
  #             matrix
  # - ArchI   > Archetypal I ridge estimator of a precision matrix
  # - ArchII  > Archetypal II ridge estimator of a precision matrix
  # - target  > target (precision terms) for Type I estimators,
  #             default = default.target(S)
  #
  # - NOTES:
  # - When type = "Alt" and target is p.d., one obtains the
  #   van Wieringen-Peeters type I estimator
  # - When type = "Alt" and target is null-matrix, one obtains the
  #   van Wieringen-Peeters type II estimator
  # - When target is not the null-matrix it is expected to be p.d. for the
  #   vWP type I estimator
  # - The target is always expected to be p.d. in case of the archetypal I
  #   estimator
  # - When type = "Alt" and target is null matrix or of form c * diag(p), a
  #   rotation equivariant estimator ensues. In these cases the expensive
  #   matrix square root can be circumvented
  ##############################################################################

  if (!isSymmetric(S)) {
    stop("S should be a symmetric matrix")
  }
  else if (lambda <= 0) {
    stop("lambda should be positive")
  }
  else if (!(type %in% c("Alt", "ArchI", "ArchII"))){
    stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
  }
  else{
    # Calculate Ridge estimator
    # Alternative estimator
    if (type == "Alt"){
      if (!isSymmetric(target)) {
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]) {
        stop("S and target should be of the same dimension")
      } else {
        P_Alt <- .armaRidgeP(S, target, lambda)
      }
      dimnames(P_Alt) <- dimnames(S)
      return(P_Alt)
    }

    # Archetypal I
    if (type == "ArchI"){
      if (lambda > 1){
        stop("lambda should be in (0,1] for this type of Ridge estimator")
      } else if (!isSymmetric(target)){
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("S and target should be of the same dimension")
      } else if (any(eigen(target, symmetric = TRUE,
                           only.values = TRUE)$values <= 0)){
        stop("Target should always be p.d. for this type of ridge estimator")
      } else {
        P_ArchI <- solve((1-lambda) * S + lambda * solve(target))
        return(P_ArchI)
      }
    }

    # Archetypal II
    if (type == "ArchII"){
      P_ArchII <- solve(S + lambda * diag(nrow(S)))
      return(P_ArchII)
    }
  }
}




##------------------------------------------------------------------------------
##
## Functions for Penalty Parameter selection
##
##------------------------------------------------------------------------------

optPenalty.LOOCV <- function(Y, lambdaMin, lambdaMax, step, type = "Alt",
                             target = default.target(covML(Y)),
                             output = "light", graph = TRUE, verbose = TRUE) {
  ##############################################################################
  # - Function that selects the optimal penalty parameter by leave-one-out
  #   cross-validation
  # - Y           > (raw) Data matrix, variables in columns
  # - lambdaMin   > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax   > maximum value penalty parameter (dependent on 'type')
  # - step        > determines the coarseness in searching the grid
  #                 [lambdaMin, lambdaMax]
  # - type        > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - target      > target (precision terms) for Type I estimators,
  #                 default = default.target(covML(Y))
  # - output      > must be one of {"all", "light"}, default = "light"
  # - graph       > Optional argument for visualization optimal penalty
  #                 selection, default = TRUE
  # - verbose     > logical indicating if intermediate output should be printed
  #                 on screen
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")
  # require("graphics")

  if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  if (verbose){
    cat("Perform input checks...", "\n")
  }
  if (!is.matrix(Y)){
    stop("Input (Y) should be a matrix")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("Input (lambdaMin) must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("Input (lambdaMin) must be positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("Input (lambdaMax) must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("Input (lambdaMax) must be larger than lambdaMin")
  }
  else if (class(step) != "numeric"){
    stop("Input (step) is of wrong class")
  }
  else if (!.is.int(step)){
    stop("Input (step) should be integer")
  }
  else if (step <= 0){
    stop("Input (step) should be a positive integer")
  }
  else if (!(output %in% c("all", "light"))){
    stop("Input (output) should be one of {'all', 'light'}")
  }
  else if (class(graph) != "logical"){
    stop("Input (graph) is of wrong class")
  }
  else {
    # Set preliminaries
    LLs     <- numeric()
    lambdas <- seq(lambdaMin, lambdaMax, len = step)

    # Calculate CV scores
    if (verbose) {
      cat("Calculating cross-validated negative log-likelihoods...\n")
    }
    for (k in 1:length(lambdas)){
      slh <- numeric()
      for (i in 1:nrow(Y)){
        S   <- covML(Y[-i,])
        slh <- c(slh, .LL(t(Y[i,,drop = F]) %*% Y[i,,drop = F],
                          ridgeP(S, lambdas[k], type = type, target = target)))
      }

      LLs <- c(LLs, mean(slh))
      if (verbose){cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")}
    }

    # Visualization
    optLambda <- min(lambdas[which(LLs == min(LLs))])
    if (graph){
      if (type == "Alt"){Main = "Alternative ridge estimator"}
      if (type == "ArchI"){Main = "Archetypal I ridge estimator"}
      if (type == "ArchII"){Main = "Archetypal II ridge estimator"}
      plot(log(lambdas), type = "l", LLs, axes = FALSE,
           xlab = "ln(penalty value)", ylab = "LOOCV neg. log-likelihood",
           main = Main)
      axis(2, ylim = c(min(LLs),max(LLs)), col = "black", lwd = 1)
      axis(1, col = "black", lwd = 1)
      abline(h = min(LLs), v = log(optLambda), col = "red")
      legend("topright",
             legend = c(paste("min. LOOCV neg. LL: ", round(min(LLs),3),sep=""),
                        paste("Opt. penalty: ", optLambda, sep = "")), cex = .8)
    }

    # Return
    S <- covML(Y)
    if (output == "all"){
      return(list(optLambda = optLambda,
                  optPrec = ridgeP(S, optLambda, type = type, target = target),
                  lambdas = lambdas, LLs = LLs))
    }
    if (output == "light"){
      return(list(optLambda = optLambda,
                  optPrec = ridgeP(S, optLambda, type = type, target = target)))
    }
  }
}



optPenalty.aLOOCV <- function(Y, lambdaMin, lambdaMax, step, type = "Alt",
                              target = default.target(covML(Y)),
                              output = "light", graph = TRUE, verbose = TRUE) {
  ##############################################################################
  # - Function that selects the optimal penalty parameter by approximate
  #   leave-one-out cross-validation
  # - Y           > (raw) Data matrix, variables in columns
  # - lambdaMin   > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax   > maximum value penalty parameter (dependent on 'type')
  # - step        > determines the coarseness in searching the grid
  #                 [lambdaMin, lambdaMax]
  # - type        > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - target      > target (precision terms) for Type I estimators,
  #                 default = default.target(covML(Y))
  # - output      > must be one of {"all", "light"}, default = "light"
  # - graph       > Optional argument for visualization optimal penalty
  #                 selection, default = TRUE
  # - verbose     > logical indicating if intermediate output should be printed
  #                 on screen
  ##############################################################################

  # Dependencies
  # require("base")
  # require("graphics")

  if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  if (verbose){
    cat("Perform input checks...", "\n")
  }
  if (!is.matrix(Y)){
    stop("Input (Y) should be a matrix")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("Input (lambdaMin) must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("Input (lambdaMin) must be positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("Input (lambdaMax) must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("Input (lambdaMax) must be larger than lambdaMin")
  }
  else if (class(step) != "numeric"){
    stop("Input (step) is of wrong class")
  }
  else if (!.is.int(step)){
    stop("Input (step) should be integer")
  }
  else if (step <= 0){
    stop("Input (step) should be a positive integer")
  }
  else if (!(output %in% c("all", "light"))){
    stop("Input (output) should be one of {'all', 'light'}")
  }
  else if (class(graph) != "logical"){
    stop("Input (graph) is of wrong class")
  }
  else {
    # Set preliminaries
    S       <- covML(Y)
    n       <- nrow(Y)
    lambdas <- seq(lambdaMin, lambdaMax, len = step)
    aLOOCVs <- numeric()

    # Calculate approximate LOOCV scores
    if (verbose){cat("Calculating approximate LOOCV scores...", "\n")}
    if (type == "Alt" & all(target == 0)){
      if (!isSymmetric(target)){
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("Covariance matrix based on data input (Y) and target should be ",
             "of the same dimension")
      } else {
        Spectral <- eigen(S, symmetric = TRUE)
        Vinv     <- solve(Spectral$vectors)
        for (k in 1:length(lambdas)){
          Eigshrink <- .eigShrink(Spectral$values, lambdas[k])
          P         <- t(Vinv) %*% diag(1/Eigshrink) %*% Vinv
          nLL       <- .5 * .LL(S, P)
          isum      <- numeric()

          for (i in 1:n){
            S1   <- t(Y[i,,drop = FALSE]) %*% Y[i,,drop = FALSE]
            isum <- c(isum, sum((solve(P) - S1) * (P %*% (S - S1) %*% P)))
          }

          aLOOCVs <- c(aLOOCVs, nLL + 1/(2 * n^2 - 2 * n) * sum(isum))
          if (verbose){cat(paste("lambda = ", lambdas[k], " done\n", sep = ""))}
        }
      }
    } else if (type == "Alt" & all(target[!diag(nrow(target))] == 0) &
                 (length(unique(diag(target))) == 1)) {
      if (!isSymmetric(target)){
        stop("Shrinkage target should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("Covariance matrix based on data input (Y) and target should be ",
             "of the same dimension")
      } else if (any(diag(target) <= 0)){
        stop("Input (target) should be p.d.")
      } else {
        varPhi   <- unique(diag(target))
        Spectral <- eigen(S, symmetric = TRUE)
        Vinv     <- solve(Spectral$vectors)
        for (k in 1:length(lambdas)){
          Eigshrink <- .eigShrink(Spectral$values, lambdas[k], const = varPhi)
          P         <- t(Vinv) %*% diag(1/Eigshrink) %*% Vinv
          nLL       <- .5 * .LL(S, P)
          isum      <- numeric()

          for (i in 1:n){
            S1   <- t(Y[i,,drop = FALSE]) %*% Y[i,,drop = FALSE]
            isum <- c(isum, sum((solve(P) - S1) * (P %*% (S - S1) %*% P)))
          }

          aLOOCVs <- c(aLOOCVs, nLL + 1/(2 * n^2 - 2 * n) * sum(isum))
          if (verbose){cat(paste("lambda = ", lambdas[k], " done\n", sep = ""))}
        }
      }
    } else {
      for (k in 1:length(lambdas)){
        P    <- ridgeP(S, lambdas[k], type = type, target = target)
        nLL  <- .5 * .LL(S, P)
        isum <- numeric()

        for (i in 1:n){
          S1   <- t(Y[i,,drop = FALSE]) %*% Y[i,,drop = FALSE]
          isum <- c(isum, sum((solve(P) - S1) * (P %*% (S - S1) %*% P)))
        }

        aLOOCVs <- c(aLOOCVs, nLL + 1/(2 * n^2 - 2 * n) * sum(isum))
        if (verbose){cat(paste("lambda = ", lambdas[k], " done\n", sep = ""))}
      }
    }

    # Visualization
    optLambda <- min(lambdas[which(aLOOCVs == min(aLOOCVs))])
    if (graph){
      if (type == "Alt"){Main = "Alternative ridge estimator"}
      if (type == "ArchI"){Main = "Archetypal I ridge estimator"}
      if (type == "ArchII"){Main = "Archetypal II ridge estimator"}
      plot(log(lambdas), type = "l", aLOOCVs, axes = FALSE,
           xlab = "ln(penalty value)",
           ylab = "Approximate LOOCV neg. log-likelihood", main = Main)
      axis(2, ylim = c(min(aLOOCVs),max(aLOOCVs)), col = "black", lwd = 1)
      axis(1, col = "black", lwd = 1)
      abline(h = min(aLOOCVs), v = log(optLambda), col = "red")
      legend("topright",
             legend = c(paste("min. approx. LOOCV neg. LL: ",
                              round(min(aLOOCVs),3), sep = ""),
                        paste("Opt. penalty: ", optLambda, sep = "")), cex = .8)
    }

    # Return
    if (output == "all"){
      return(list(optLambda = optLambda,
                  optPrec = ridgeP(S, optLambda, type = type, target = target),
                  lambdas = lambdas, aLOOCVs = aLOOCVs))
    }
    if (output == "light"){
      return(list(optLambda = optLambda,
                  optPrec = ridgeP(S, optLambda, type = type, target = target)))
    }
  }
}



optPenalty.LOOCVauto <- function (Y, lambdaMin, lambdaMax,
                                  lambdaInit = (lambdaMin + lambdaMax)/2,
                                  target = default.target(covML(Y)),
                                  type = "Alt") {
  ##############################################################################
  # - Function that determines the optimal value of the penalty parameter by
  #   application of the Brent algorithm to the (leave-one-out) cross-validated
  #   log-likelihood
  # - Y          > (raw) Data matrix, variables in columns
  # - lambdaMin  > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax  > maximum value penalty parameter (dependent on 'type')
  # - lambdaInit > initial value for lambda for starting optimization
  # - target     > target (precision terms) for Type I estimators,
  #                default = default.target(covML(Y))
  # - type       > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")

  # input checks
  if (!is.matrix(Y)){
    stop("Input (Y) if of wrong class")
  }
  else if (sum(is.na(Y)) != 0){
    stop("Input matrix (Y) should not contain missings")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("Input (lambdaMin) must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("Input (lambdaMin) must be strictly positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("Input (lambdaMax) must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("Input (lambdaMax) must be larger than lambdaMin")
  }
  else if (class(lambdaInit) != "numeric"){
    stop("Input (lambdaInit) is of wrong class")
  }
  else if (length(lambdaInit) != 1){
    stop("Input (lambdaInit) must be a scalar")
  }
  else if (lambdaInit <= lambdaMin){
    stop("Input (lambdaInit) must be larger than input (lambdaMin)")
  }
  else if (lambdaMax <= lambdaInit){
    stop("Input (lambdaInit) must be smaller than input (lambdaMax)")
  }
  else {
    # determine optimal value of ridge penalty parameter
    optLambda <- optim(lambdaInit, .cvl, method = "Brent", lower = lambdaMin,
                       upper = lambdaMax, Y = Y, target = target,
                       type = type)$par

    # Return
    return(list(optLambda = optLambda,
                optPrec = ridgeP(covML(Y), optLambda,
                                 type = type, target = target)))
  }
}



conditionNumberPlot <- function(S, lambdaMin, lambdaMax, step, type = "Alt",
                                target = default.target(S), norm = "2",
                                digitLoss = FALSE, rlDist = FALSE,
                                vertical = FALSE, value, main = TRUE,
                                nOutput = FALSE, verbose = TRUE){
  ##############################################################################
  # - Function that visualizes the spectral condition number against the
  #   regularization parameter
  # - Can be used to heuristically determine the (minimal) value of the penalty
  #   parameter
  # - The ridge estimators operate by shrinking the eigenvalues
  # - This is especially the case when targets are used that lead to rotation
  #   equivariant estimators
  # - Maximum shrinkage (under rotation equivariance) implies that all
  #   eigenvalues will be equal
  # - Ratio of maximum and minimum eigenvalue of P can then function as
  #   a heuristic
  # - It's point of stabilization can give an acceptable value for the penalty
  # - The ratio boils down to the (spectral) condition number of a matrix
  # - S         > sample covariance/correlation matrix
  # - lambdaMin > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax > maximum value penalty parameter (dependent on 'type')
  # - step      > determines the coarseness in searching the grid
  #               [lambdaMin, lambdaMax].The steps on the grid are equidistant
  #               on the log scale
  # - type      > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - target    > target (precision terms) for Type I estimators,
  #               default = default.target(S)
  # - norm      > indicates the norm under which the condition number is to be
  #               estimated. Default is the L2-norm. The L1-norm can be
  #               (cheaply) approximated
  # - digitLoss > logical indicating if the approximate loss in digits of
  #               accuracy should also be plotted. Default = FALSE
  # - rlDist    > logical indicating if relative distance to set of singular
  #               matrices should also be plotted. Default = FALSE
  # - vertical  > optional argument for visualization vertical line in graph
  #               output, default = FALSE. Can be used to indicate the value of,
  #               e.g., the optimal penalty as indicated by some routine. Can be
  #               used to assess if this optimal penalty will lead to a
  #               well-conditioned estimate
  # - value     > indicates constant on which to base vertical line when
  #               vertical = TRUE
  # - main      > logical indicating if plot should contain type of estimator as
  #               main title
  # - nOutput   > logical indicating if numeric output should be given (lambdas
  #               and condition numbers)
  # - verbose   > logical indicating if intermediate output should be printed on
  #               screen
  ##############################################################################

  # Dependencies
  # require("base")
  # require("graphics")
  # require("Hmisc")
  # require("sfsmisc")

  if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  if (verbose){
    cat("Perform input checks...", "\n")
  }
  if (!is.matrix(S)){
    stop("S should be a matrix")
  }
  else if (!isSymmetric(S)){
    stop("S should be a symmetric matrix")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("lambdaMin must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("lambdaMin must be positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("lambdaMax must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("lambdaMax must be larger than lambdaMin")
  }
  else if (class(step) != "numeric"){
    stop("Input (step) is of wrong class")
  }
  else if (!.is.int(step)){
    stop("step should be integer")
  }
  else if (step <= 0){
    stop("step should be a positive integer")
  }
  else if (!(type %in% c("Alt", "ArchI", "ArchII"))){
    stop("type should be one of {'Alt', 'ArchI', 'ArchII'}")
  }
  else if (!isSymmetric(target)){
    stop("Shrinkage target should be symmetric")
  }
  else if (dim(target)[1] != dim(S)[1]){
    stop("S and target should be of the same dimension")
  }
  else if (type == "Alt" & !all(target == 0) &
           any(eigen(target, symmetric = TRUE, only.values = T)$values <= 0)){
    stop("When target is not a null-matrix it should be p.d.")
  }
  else if (type == "ArchI" & lambdaMax > 1){
    stop("lambda should be in (0,1] for this type of Ridge estimator")
  }
  else if (type == "ArchI" &
           any(eigen(target, symmetric = TRUE, only.values = T)$values <= 0)){
    stop("Target should be p.d.")
  }
  else if (!(norm %in% c("2", "1"))){
    stop("norm should be one of {'2', '1'}")
  }
  else if (class(digitLoss) != "logical"){
    stop("Input (digitLoss) is of wrong class")
  }
  else if (class(rlDist) != "logical"){
    stop("Input (rlDist) is of wrong class")
  }
  else if (digitLoss & rlDist){
    stop("Only one of 'digitLoss' and 'rlDist' may be TRUE")
  }
  else if (class(vertical) != "logical"){
    stop("Input (vertical) is of wrong class")
  }
  else if (class(main) != "logical"){
    stop("Input (main) is of wrong class")
  }
  else if (class(nOutput) != "logical"){
    stop("Input (nOutput) is of wrong class")
  }
  else {
    # Set preliminaries
    lambdas <- lseq(lambdaMin, lambdaMax, length = step)
    condNR  <- numeric()

    if (norm == "2"){
      # Calculate spectral condition number ridge estimate on lambda grid
      if (verbose){cat("Calculating spectral condition numbers...", "\n")}
      if (type == "Alt" & all(target == 0)){
        Spectral <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
        for (k in 1:length(lambdas)){
          Eigshrink <- .armaEigShrink(Spectral, lambdas[k])
          condNR[k] <- as.numeric(max(Eigshrink)/min(Eigshrink))
        }
      } else if (type == "Alt" & all(target[!diag(nrow(target))] == 0)
                 & (length(unique(diag(target))) == 1)){
        varPhi   <- unique(diag(target))
        Spectral <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
        for (k in 1:length(lambdas)){
          Eigshrink <- .armaEigShrink(Spectral, lambdas[k], cons = varPhi)
          condNR[k] <- as.numeric(max(Eigshrink)/min(Eigshrink))
        }
      } else {
        if (type == "Alt"){
          for (k in 1:length(lambdas)){
            P         <- .armaRidgePAnyTarget(S, target = target, lambdas[k])
            Eigs      <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
            condNR[k] <- as.numeric(max(Eigs)/min(Eigs))
          }
        }
        if (type != "Alt"){
          for (k in 1:length(lambdas)){
            P         <- .ridgeSi(S, lambdas[k], type = type, target = target)
            Eigs      <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
            condNR[k] <- as.numeric(max(Eigs)/min(Eigs))
          }
        }
      }
    }

    if (norm == "1"){
      # Calculate approximation to condition number under 1-norm
      if (verbose){cat("Approximating condition number under 1-norm...", "\n")}
      if (type == "Alt"){
        for (k in 1:length(lambdas)){
          P         <- .armaRidgeP(S, target = target, lambdas[k])
          condNR[k] <- as.numeric(1/rcond(P, norm = "O"))
        }
      }
      if (type != "Alt"){
        for (k in 1:length(lambdas)){
          P         <- .ridgeSi(S, lambdas[k], type = type, target = target)
          condNR[k] <- as.numeric(1/rcond(P, norm = "O"))
        }
      }
    }

    # Visualization
    if (verbose){cat("Plotting...", "\n")}
    if (main){
      if (type == "Alt"){Main = "Alternative ridge estimator"}
      if (type == "ArchI"){Main = "Archetypal I ridge estimator"}
      if (type == "ArchII"){Main = "Archetypal II ridge estimator"}
    }
    if (!main){Main = " "}
    if (norm == "2"){Ylab = "spectral condition number"}
    if (norm == "1"){Ylab = "condition number under 1-norm"}
    if (digitLoss | rlDist){par(mar = c(5,4,4,5)+.1)}
    plot(log(lambdas), type = "l", condNR, axes = FALSE, col = "blue4",
         xlab = "ln(penalty value)", ylab = Ylab, main = Main)
    axis(2, ylim = c(0,max(condNR)), col = "black", lwd = 1)
    axis(1, col = "black", lwd = 1)
    minor.tick(nx = 10, ny = 0, tick.ratio = .4)
    par(xpd = FALSE)
    if (digitLoss){
      dLoss <- floor(log10(condNR))
      par(new = TRUE)
      plot(log(lambdas), dLoss, axes = FALSE, type = "l", col = "green3",
           xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      axis(4, col = "black", lwd = 1)
      mtext("Loss in digits of accuracy", side = 4, line = 3)
      legend("top", col=c("blue4","green3"), lty = 1, legend =
               c("Condition number", "floor(log10(Condition number))"), cex = .8)
    }
    if (rlDist){
      RlDist <- 1/condNR
      par(new = TRUE)
      plot(log(lambdas), RlDist, axes = FALSE, type = "l", col = "green3",
           xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      axis(4, col = "black", lwd = 1)
      mtext("relative distance to singular matrix", side = 4, line = 3)
      legend("top", col=c("blue4","green3"), lty = 1, legend =
               c("Condition number", "Relative distance"), cex = .8)
    }
    if (vertical){
      if (missing(value)){
        stop("Need to specify input (value)")
      } else if (class(value) != "numeric"){
        stop("Input (value) is of wrong class")
      } else if (length(value) != 1){
        stop("Input (value) must be a scalar")
      } else if (value <= 0){
        stop("Input (value) must be positive")
      } else {
        abline(v = log(value), col = "red")
      }
    }

    # Possible output
    if (nOutput){
      return(list(lambdas = lambdas, conditionNumbers = condNR))
    }
  }
}




##------------------------------------------------------------------------------
##
## Functions for Block Independence Testing and Mutual Information
##
##------------------------------------------------------------------------------

if(getRversion() >= "2.15.1") utils::globalVariables("rags2ridges")

GGMblockNullPenalty <- function(Y, id, nPerm = 25, lambdaMin, lambdaMax,
                                lambdaInit = (lambdaMin+lambdaMax)/2,
                                target = default.target(covML(Y)),
                                type = "Alt", ncpus = 1){
  ##############################################################################
  # - Function generating the distribution of the penalty parameter
  # - It does so under the null hypothesis of block independence
  # - The optimal value of the penalty parameter is determined under multiple
  #   permutations
  # - Optimal penalty determined using the 'optPenalty.LOOCVauto' function
  # - Y          > (raw) Data matrix, variables in columns
  # - id         > indicator variable for the two blocks of the precision matrix
  # - nPerm      > desired number of permutations
  # - lambdaMin  > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax  > maximum value penalty parameter (dependent on 'type')
  # - lambdaInit > initial value for lambda for starting optimization
  # - target     > target (precision terms) for Type I estimators,
  #                default = default.target(covML(Y))
  # - type       > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - ncpus      > desired number of cpus
  #
  # Notes:
  # - Dependency on 'snowfall' when ncpus > 1
  ##############################################################################

  # Dependencies
  # require("base")
  # require("snowfall")

  if (!is.matrix(Y)){
    stop("Input (Y) should be a matrix")
  }
  else if (sum(is.na(Y)) != 0){
    stop("Input (Y) should not contain missings")
  }
  else if (class(id) != "numeric" & class(id) != "integer"){
    stop("Input (id) is of wrong class")
  }
  else if (!(all(unique(id) %in% c(0, 1)))){
    stop("Input (id) has unlawful entries")
  }
  else if (length(id) != ncol(Y)){
    stop("Column dimension input (Y) does not match with length input (id)")
  }
  else if (class(nPerm) != "numeric" & class(nPerm) != "integer"){
    stop("Input (nPerm) is of wrong class")
  }
  else if (!.is.int(nPerm)){
    stop("Input (nPerm) is expected to be a (numeric) integer")
  }
  else if (nPerm <= 0){
    stop("Input (nPerm) must be strictly positive")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("Input (lambdaMin) must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("Input (lambdaMin) must be strictly positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("Input (lambdaMax) must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("Input (lambdaMax) must be larger than input (lambdaMin)")
  }
  else if (class(lambdaInit) != "numeric"){
    stop("Input (lambdaInit) is of wrong class")
  }
  else if (length(lambdaInit) != 1){
    stop("Input (lambdaInit) must be a scalar")
  }
  else if (lambdaInit <= lambdaMin){
    stop("Input (lambdaInit) must be larger than input (lambdaMin)")
  }
  else if (lambdaMax <= lambdaInit){
    stop("Input (lambdaInit) must be smaller than input (lambdaMax)")
  }
  else if (class(ncpus) != "numeric" & class(ncpus) != "integer"){
    stop("Input (ncpus) is of wrong class")
  }
  else if (!.is.int(ncpus)){
    stop("Input (ncpus) is expected to be a (numeric) integer")
  }
  else if (ncpus <= 0){
    stop("Input (ncpus) must be strictly positive")
  }
  else {
    # Determine null distribution
    if (ncpus == 1){
      lambdaNull <- sapply(1:nPerm, .lambdaNullDist, Y = Y, id = id,
                           lambdaMin = lambdaMin, lambdaMax = lambdaMax,
                           lambdaInit = lambdaInit, target = target,
                           type = type)
    }
    if (ncpus > 1){
      sfInit(parallel = TRUE, cpus = ncpus)
      sfLibrary(rags2ridges, verbose = FALSE)
      lambdaNull <- sfSapply(1:nPerm, .lambdaNullDist, Y = Y, id = id,
                             lambdaMin = lambdaMin, lambdaMax = lambdaMax,
                             lambdaInit = lambdaInit, target = target,
                             type = type)
      sfStop()
    }

    # Return
    return(lambdaNull)
  }
}



GGMblockTest <- function (Y, id, nPerm = 1000, lambda,
                          target = default.target(covML(Y)), type = "Alt",
                          lowCiThres = 0.1, ncpus = 1, verbose = TRUE) {
  ##############################################################################
  # - Function performing a permutation test for block structure in the
  #   precision matrix
  # - The setting is a high-dimensional one
  # - Y          > (raw) Data matrix, variables in columns
  # - id         > indicator variable for the two blocks of the precision matrix
  # - nPerm      > desired number of permutations, default = 1000
  # - lambda     > the penalty parameter employed in the permutation test
  # - target     > target (precision terms) for Type I estimators,
  #                default = default.target(covML(Y))
  # - type       > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - lowCiThres > A value between 0 and 1, determining speed of efficient
  #                p-value calculation
  # - ncpus      > desired number of cpus
  # - verbose    > logical indicating if progress/output should be printed on
  #                screen
  #
  # Notes:
  # - Dependency on 'snowfall' when ncpus > 1
  # - When verbose = TRUE, also graphical output is given: a histogram of the
  #   null-distribution
  # - The value for 'lambda' ideally stems from 'GGMblockNullPenalty'
  # - If the probability of a p-value being below 'lowCiThres' is smaller than
  #   0.001 (meaning: the test is unlikely to become significant), the
  #   permutation analysis is terminated and a p-value of unity (1) is reported
  ##############################################################################

  # Dependencies
  # require("base")
  # require("snowfall")
  # require("graphics")

  if (!is.matrix(Y)){
    stop("Input (Y) should be a matrix")
  }
  else if (sum(is.na(Y)) != 0){
    stop("Input (Y) should not contain missings")
  }
  else if (class(id) != "numeric" & class(id) != "integer"){
    stop("Input (id) is of wrong class")
  }
  else if (!(all(unique(id) %in% c(0, 1)))){
    stop("Input (id) has unlawful entries")
  }
  else if (length(id) != ncol(Y)){
    stop("Column dimension input (Y) does not match with length input (id)")
  }
  else if (class(nPerm) != "numeric" & class(nPerm) != "integer"){
    stop("Input (nPerm) is of wrong class")
  }
  else if (!.is.int(nPerm)){
    stop("Input (nPerm) is expected to be a (numeric) integer")
  }
  else if (nPerm <= 0){
    stop("Input (nPerm) must be strictly positive")
  }
  else if (class(lambda) != "numeric"){
    stop("Input (lambda) is of wrong class")
  }
  else if (length(lambda) != 1){
    stop("Input (lambda) must be a scalar")
  }
  else if (lambda <= 0){
    stop("Input (lambda) must be strictly positive")
  }
  else if (class(lowCiThres) != "numeric"){
    stop("Input (lowCiThres) is of wrong class")
  }
  else if (length(lowCiThres) != 1){
    stop("Input (lowCiThres) must be a scalar")
  }
  else if (lowCiThres <= 0 | lowCiThres >= 1){
    stop("Input (lowCiThres) must be in the interval (0,1)")
  }
  else if (class(ncpus) != "numeric" & class(ncpus) != "integer"){
    stop("Input (ncpus) is of wrong class")
  }
  else if (!.is.int(ncpus)){
    stop("Input (ncpus) is expected to be a (numeric) integer")
  }
  else if (ncpus <= 0){
    stop("Input (ncpus) must be strictly positive")
  }
  else if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  else {
    # Observed test statistics
    S     <- solve(ridgeP(covML(Y), lambda = lambda,
                          target = target, type = type))
    llObs <- log(det(S[id == 0, id == 0, drop = FALSE])) +
      log(det(S[id == 1, id == 1, drop = FALSE])) - log(det(S))

    # Define steps at which the possibility of significance should be evaluated
    steps <- sort(unique(c(0, 25, 50, 100, 150, 200,
                           seq(from = 250, to = 2750, by = 250),
                           seq(from = 3000, to = 10000, by = 500),
                           seq(from = 11000, to = 50000, by = 1000), nPerm)))
    steps <- steps[steps <= nPerm]

    # Generate null distribution
    nullDist <- numeric()
    if (ceiling(ncpus) > 1){
      sfInit(parallel = TRUE, cpus = ncpus)
      sfLibrary(rags2ridges, verbose = FALSE)
    }

    for (j in 1:(length(steps) - 1)){
      if (verbose){
        cat(paste(steps[j], " of ", steps[length(steps)],
                  " permutations done, and counting...\n", sep = ""))
      }
      if (ncpus == 1){
        nullDistPart <- sapply(c((steps[j] + 1):steps[j + 1]), .blockTestStat,
                               Y = Y, id = id, lambda = lambda, target = target,
                               type = type)
      }
      if (ncpus > 1){
        nullDistPart <- sfSapply(c((steps[j] + 1):steps[j + 1]), .blockTestStat,
                                 Y = Y, id = id, lambda = lambda,
                                 target = target, type = type)
      }
      nullDist <- c(nullDist, nullDistPart); rm(nullDistPart); gc()
      pVal     <- sum(nullDist >= as.numeric(llObs))/steps[j + 1]
      pBound   <- pVal - sqrt(pVal * (1 - pVal)/steps[j + 1]) * 3.09
      significanceUnlikely <- (pBound > lowCiThres)
      if (significanceUnlikely){pVal <- 1; break}
      if (verbose){cat(paste(steps[j + 1], "of", steps[length(steps)],
                             " permutations done", sep = " "), "\n")}
    }

    if (ncpus > 1){sfStop()}

    # Generating on screen (graphical) output
    if (verbose){
      # Visual summary of test results
      xlims     <- c(min(c(llObs, nullDist)), max(c(llObs, nullDist)))
      histFreqs <- hist(nullDist, n = sqrt(sum(nPerm))+1, plot = FALSE)$counts
      hist(nullDist, xlim = xlims, n = sqrt(sum(nPerm))+1, col = "blue",
           border = "lightblue", xlab = "null statistics",
           ylab = "frequency", main = "Histogram of null distribution")
      lines(rep(llObs, max(histFreqs)), 0.9 * (c(1:max(histFreqs))-1),
            col = "red", lwd = 2)
      text(quantile(c(nullDist, llObs), probs = 0.05), 0.95 * max(histFreqs),
           paste("p-value:", pVal))
      text(llObs, 0.95 * max(histFreqs), "test stat.")

      # Summary of test results
      remark <-
        ifelse(significanceUnlikely,
               "resampling terminated prematurely due to unlikely significance",
               "none")
      cat("\n")
      cat("Likelihood ratio test for block independence\n")
      cat("----------------------------------------\n")
      cat("-> number of permutations : ", nPerm, "\n", sep="")
      cat("-> test statistic         : ", round(llObs, digits = 3), "\n",sep="")
      cat("-> p-value                : ", round(pVal, digits = 3), "\n", sep="")
      cat("-> remark                 : ", remark, "\n", sep="")
      cat("----------------------------------------\n")
      cat("\n")
    }

    # Return
    return(list(statistic = llObs, pvalue = pVal, nulldist = nullDist,
                nperm = nPerm, remark = remark))
  }
}



GGMmutualInfo <- function(S, split1){
  ##############################################################################
  # - Function that calculates the mutual information between two exhaustive and
  #   mutually exclusive splits of normal p-variate random variable.
  # - S      > Sample Covariance matrix.
  # - split1 > A numeric indicating the variates (by column number) forming the
  #            first split. The second split is automatically formed from its
  #            complement.
  #
  # NOTES:
  # - No dependencies at current
  ##############################################################################

  if (!is.matrix(S)){
    stop("Input (S) should be a matrix")
  }
  else if (!isSymmetric(S)){
    stop("Input (S) should be a symmetric matrix")
  }
  else if (length(split1) < 1 & length(split1) > nrow(S)-1){
    stop("Input (split1) is of wrong length.")
  }
  else {
    # mutual information
    MI <- log(det(S[-split1,-split1])) -
          log(det(S[-split1,-split1] - S[-split1,split1,drop=FALSE] %*%
                  solve(S[split1,split1,drop=FALSE])
                  %*% S[split1,-split1,drop=FALSE]))
    return(MI)
  }
}




##------------------------------------------------------------------------------
##
## Test for Vanishing Partial Correlations
##
##------------------------------------------------------------------------------

sparsify <- function(P, threshold = c("absValue", "localFDR", "top"),
                     absValueCut = .25, FDRcut = .9,
                     top = 10, output = "heavy", verbose = TRUE){
  ##############################################################################
  # - Function that sparsifies/determines support of a partial correlation
  #   matrix
  # - Support can be determined by absolute value thresholding or by local FDRs
  #   thresholding
  # - One can also choose to threshold based on the top X of absolute partial
  #   correlations
  # - Local FDR operates on the nonredundant non-diagonal elements of a partial
  #   correlation matrix
  # - Function is to some extent a wrapper around certain 'fdrtool' functions
  # - P           > (possibly shrunken) precision matrix
  # - threshold   > signifies type of thresholding
  # - absValueCut > cut-off for partial correlation elements selection based on
  #                 absolute value thresholding.
  #                 Only when threshold = 'absValue'. Default = .25
  # - FDRcut      > cut-off for partial correlation element selection based on
  #                 local FDR thresholding
  #                 Only when threshold = 'localFDR'. Default = .9
  # - top         > partial correlation element selection based on retainment
  #                 'top' number of absolute partial correlation.
  #                 Only when threshold = 'top'. Default = 10
  # - output      > must be one of {"heavy", "light"}, default = "heavy"
  # - verbose     > logical indicating if intermediate output should be printed
  #                 on screen. Only when threshold = 'localFDR'. Default = TRUE.
  #
  # NOTES:
  # - Input (P) may be the partial correlation matrix or the standardized
  #   precision matrix. These are identical up to the signs of off-diagonal
  #   elements. Either can be used as it has no effect on the thresholding
  #   operator and the ensuing sparsified result.
  # - Input (P) may also be the unstandardized precision matrix. The function
  #   converts it to the partial correlation matrix
  # - The function evaluates if the input (P) is a partial
  #   correlation/standardized precision matrix or an unstandardized precision
  #   matrix. If the input amounts to the latter both the sparsified partial
  #   correlation matrix and the corresponding sparsified precision matrix are
  #   given as output (when output = "heavy"). Otherwise, the ouput consists of
  #   the sparsified partial correlation/standardized precision matrix.
  # - When output = "light", only the (matrix) positions of the zero and
  #   non-zero elements are returned.
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")
  # require("fdrtool")

  if (!is.matrix(P)){
    stop("Input (P) should be a matrix")
  }
  else if (!isSymmetric(P)){
    stop("Input (P) should be a symmetric matrix")
  }
  else if (!evaluateS(P, verbose = FALSE)$posEigen){
    stop("Input (P) is expected to be positive definite")
  }
  else if (missing(threshold)){
    stop("Need to specify type of sparsification ('absValue' or 'localFDR' ",
         "or 'top')")
  }
  else if (!(threshold %in% c("absValue", "localFDR", "top"))){
    stop("Input (threshold) should be one of {'absValue', 'localFDR', 'top'}")
  }
  else if (!(output %in% c("light", "heavy"))){
    stop("Input (output) should be one of {'light', 'heavy'}")
  }
  else {
    # Obtain partial correlation matrix
    if (all(length(unique(diag(P))) == 1 & unique(diag(P)) == 1)){
      stan = TRUE
      PC  <- P
    } else {
      stan = FALSE
      PC  <- symm(pcor(P))
    }

    # Number of nonredundant elements
    NR <- (ncol(P)*(ncol(P) - 1))/2

    # Obtain sparsified matrix
    if (threshold == "top"){
      if (class(top) != "numeric"){
        stop("Input (top) is of wrong class")
      } else if (length(top) != 1){
        stop("Input (top) must be a scalar")
      } else if (!.is.int(top)){
        stop("Input (top) should be a numeric integer")
      } else if (top <= 0){
        stop("Input (top) must be strictly positive")
      } else if (top >= NR){
        stop("Input (top) must be smaller than the number of nonredundant ",
             "off-diagonal elements of the input matrix P")
      } else {
        absValueCut <- sort(abs(PC[upper.tri(PC)]),
                            decreasing = TRUE)[ceiling(top)]
        threshold   <- "absValue"
      }
    }

    if (threshold == "absValue"){
      if (class(absValueCut) != "numeric"){
        stop("Input (absValueCut) is of wrong class")
      } else if (length(absValueCut) != 1){
        stop("Input (absValueCut) must be a scalar")
      } else if (absValueCut <= 0 | absValueCut >= 1){
        stop("Input (absValueCut) must be in the interval (0,1)")
      } else {
        PC0 <- PC
        PC0[!(abs(PC0) >= absValueCut)] <- 0
        if (!stan){
          P0 <- P
          P0[PC0 == 0] <- 0
        }
      }
    }

    if (threshold == "localFDR"){
      if (class(FDRcut) != "numeric"){
        stop("Input (FDRcut) is of wrong class")
      } else if (length(FDRcut) != 1){
        stop("Input (FDRcut) must be a scalar")
      } else if (FDRcut <= 0 | FDRcut >= 1){
        stop("Input (FDRcut) must be in the interval (0,1)")
      } else if (class(verbose) != "logical"){
        stop("Input (verbose) is of wrong class")
      } else {
        lFDRs <- 1 - fdrtool(PC[upper.tri(PC)], "correlation",
                             plot = verbose, verbose = verbose)$lfdr
        PC0   <- diag(nrow(PC))
        PC0[lower.tri(PC0)] <- 1
        zeros <- which(PC0 == 0, arr.ind = TRUE)
        zeros <- zeros[which(lFDRs <= FDRcut),]
        PC0   <- PC
        PC0[zeros] <- 0
        PC0[cbind(zeros[,2], zeros[,1])] <- 0
        if (!stan){
          P0 <- P
          P0[PC0 == 0] <- 0
        }
      }
    }

    # Return
    NNZ <- length(which(PC0[upper.tri(PC0)] != 0))
    cat("- Retained elements: ", NNZ, "\n")
    cat("- Corresponding to", round(NNZ/NR,4) * 100,"% of possible edges \n")
    cat(" \n")

    if (output == "heavy"){
      if (stan){
        colnames(PC0) = rownames(PC0) <- colnames(P)
        return(PC0)
      }
      if (!stan){
        colnames(PC0) = rownames(PC0) <- colnames(P)
        colnames(P0)  = rownames(P0)  <- colnames(P)
        return(list(sparseParCor = PC0, sparsePrecision = P0))
      }
    }
    if (output == "light"){
      return(list(zeros = which(PC0 == 0, arr.ind = TRUE),
                  nonzeros = which(PC0 != 0, arr.ind = TRUE)))
    }
  }
}




##------------------------------------------------------------------------------
##
## Functions for Loss/Entropy/Fit Evaluation
##
##------------------------------------------------------------------------------

loss <- function(E, T, precision = TRUE, type = c("frobenius", "quadratic")){
  ##############################################################################
  # - Function evualuating various loss functions on the precision
  # - E         > Estimated (possibly regularized) precision matrix
  # - T         > True (population) covariance or precision matrix
  # - precision > Logical indicating if T is a precision matrix (when TRUE)
  # - type      > character indicating which loss function is to be used
  ##############################################################################

  if (!is.matrix(E)){
    stop("Input (E) is of wrong class")
  }
  else if (!isSymmetric(E)){
    stop("E should be a symmetric matrix")
  }
  else if (!is.matrix(T)){
    stop("Input (T) is of wrong class")
  }
  else if (!isSymmetric(T)){
    stop("T should be a symmetric matrix")
  }
  else if (dim(E)[1] != dim(T)[1]){
    stop("E and T should be of the same dimension")
  }
  else if (class(precision) != "logical"){
    stop("Input (precision) is of wrong class")
  }
  else if (missing(type)){
    stop("Need to specify loss type ('frobenius' or 'quadratic')")
  }
  else if (!(type %in% c("frobenius", "quadratic"))){
    stop("type should be one of {'frobenius', 'quadratic'}")
  }
  else {
    # Frobenius loss
    if (type == "frobenius"){
      if (precision)  {loss <- .FrobeniusLoss(E, T)}
      if (!precision) {loss <- .FrobeniusLoss(E, solve(T))}
    }

    # Quadratic loss
    if (type == "quadratic"){
      if (precision)  {loss <- .QuadraticLoss(E, solve(T))}
      if (!precision) {loss <- .QuadraticLoss(E, T)}
    }

    # Return
    return(loss)
  }
}



KLdiv <- function(Mtest, Mref, Stest, Sref, symmetric = FALSE){
  ##############################################################################
  # - Function that calculates the Kullback-Leibler divergence between two
  #   normal distributions
  # - Mtest     > mean vector approximating m.v. normal distribution
  # - Mref      > mean vector 'true'/reference m.v. normal distribution
  # - Stest     > covariance matrix approximating m.v. normal distribution
  # - Sref      > covariance matrix 'true'/reference m.v. normal distribution
  # - symmetric > logical indicating if original symmetric version of KL div.
  #               should be calculated
  ##############################################################################

  # Dependencies
  # require("base")

  if (class(Mtest) != "numeric"){
    stop("Input (Mtest) is of wrong class")
  }
  else if (class(Mref) != "numeric"){
    stop("Input (Mref) is of wrong class")
  }
  else if (length(Mtest) != length(Mref)){
    stop("Mtest and Mref should be of same length")
  }
  else if (!is.matrix(Stest)){
    stop("Input (Stest) is of wrong class")
  }
  else if (!is.matrix(Sref)){
    stop("Input (Sref) is of wrong class")
  }
  else if (!isSymmetric(Stest)){
    stop("Stest should be symmetric")
  }
  else if (!isSymmetric(Sref)){
    stop("Sref should be symmetric")
  }
  else if (dim(Stest)[1] != length(Mtest)){
    stop("Column and row dimension of Stest should correspond to length Mtest")
  }
  else if (dim(Sref)[1] != length(Mref)){
    stop("Column and row dimension of Sref should correspond to length Mref")
  }
  else if (class(symmetric) != "logical"){
    stop("Input (symmetric) is of wrong class")
  }
  else {
    # Evaluate KL divergence
    KLd <- (sum(diag(solve(Stest) %*% Sref)) +
              t(Mtest - Mref) %*% solve(Stest) %*% (Mtest - Mref) -
              nrow(Sref) - log(det(Sref)) + log(det(Stest)))/2

    # Evaluate (original) symmetric version KL divergence
    if (symmetric){
      KLd <- KLd + (sum(diag(solve(Sref) %*% Stest)) +
                      t(Mref - Mtest) %*% solve(Sref) %*% (Mref - Mtest) -
                      nrow(Sref) - log(det(Stest)) + log(det(Sref)))/2
    }

    # Return
    return(as.numeric(KLd))
  }
}



evaluateSfit <- function(Phat, S, diag = FALSE, fileType = "pdf", nameExt = "",
                         dir = getwd()){
  ##############################################################################
  # - Function aiding the visual inspection of the fit of the estimated
  #   (possibly regularized) precision matrix vis-a-vis the sample
  #   covariance matrix
  # - Phat     > (regularized) estimate of the precision matrix
  # - S        > sample covariance matrix
  # - diag     > logical determining treatment diagonal elements for plots
  # - fileType > signifies filetype of output
  # - nameExt  > character giving extension of default output names.
  #              Circumvents overwriting of output when working in single
  #              directory.
  # - dir      > specifies the directory in which the visual output is stored
  ##############################################################################

  # Dependencies
  # require("base")
  # require("graphics")

  if (!is.matrix(Phat)){
    stop("Input (Phat) should be a matrix")
  }
  else if (!isSymmetric(Phat)){
    stop("Input (Phat) should be a symmetric matrix")
  }
  else if (all(diag(Phat) == 1)){
    stop("Input (Phat) should be a nonstandardized precision matrix")
  }
  else if (!is.matrix(S)){
    stop("Input (S) should be a matrix")
  }
  else if (!isSymmetric(S)){
    stop("Input (S) should be a symmetric matrix")
  }
  else if (all(diag(S) == 1)){
    stop("Input (S) should be a nonstandardized covariance matrix")
  }
  else if (class(diag) != "logical"){
    stop("Input (diag) is of wrong class")
  }
  else if (missing(fileType)){
    stop("Need to specify type of output file ('pdf' or 'eps')")
  }
  else if (!(fileType %in% c("pdf", "eps"))){
    stop("fileType should be one of {'pdf', 'eps'}")
  }
  else if (class(nameExt) != "character"){
    stop("Input (nameExt) is of wrong class")
  }
  else if (class(dir) != "character"){
    stop("Specify directory for output as 'character'")
  }
  else {
    # Obtain estimated covariance matrix
    Shat <- solve(Phat)


    print("Visualizing covariance fit")
    # plot 1: QQ-plot of covariances
    if (fileType == "pdf"){
      pdf(paste(dir, "QQplot_covariances_", nameExt, ".pdf"))
    }
    if (fileType == "eps"){
      setEPS(); postscript(paste(dir, "QQplot_covariances_", nameExt, ".eps"))
    }
    if (diag){
      cObs <- as.numeric(S[upper.tri(S, diag = TRUE)])
      cFit <- as.numeric(Shat[upper.tri(Shat, diag = TRUE)])
    }
    if (!diag){
      cObs <- as.numeric(S[upper.tri(S)])
      cFit <- as.numeric(Shat[upper.tri(Shat)])
    }
    op <- par(pty = "s")
    qqplot(x = cObs, y = cFit, pch = 20, xlab = "sample covariances",
           ylab = "fits", main = "QQ-plot, covariances")
    lines(seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
          seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
          col = "grey", lty = 2)
    par(op); dev.off()

    # plot 2: Comparison of covariances by heatmap
    if (fileType == "pdf"){
      pdf(paste(dir, "heatmap_covariances_", nameExt, ".pdf"))
    }
    if (fileType == "eps"){
      setEPS(); postscript(paste(dir, "heatmap_covariances_", nameExt, ".eps"))
    }
    op  <- par(pty = "s")
    slh <- S
    slh[lower.tri(slh)] <- Shat[lower.tri(Shat)]
    gplot <- edgeHeat(slh, diag = diag, legend = FALSE, main = "Covariances")
    print(gplot); par(op); dev.off()


    print("Visualizing correlation fit")
    # plot 3: QQ-plot of correlations
    if (fileType == "pdf"){
      pdf(paste(dir, "QQplot_correlations_", nameExt, ".pdf"))
    }
    if (fileType == "eps"){
      setEPS(); postscript(paste(dir, "QQplot_correlations_", nameExt, ".eps"))
    }
    if (diag){
      cObs <- as.numeric(cov2cor(S)[upper.tri(S, diag = TRUE)]);
      cFit <- as.numeric(cov2cor(Shat)[upper.tri(Shat, diag = TRUE)])
    }
    if (!diag){
      cObs <- as.numeric(cov2cor(S)[upper.tri(S)]);
      cFit <- as.numeric(cov2cor(Shat)[upper.tri(Shat)])
    }
    op <- par(pty = "s")
    qqplot(x = cObs, y = cFit, pch = 20, xlab = "sample correlations",
           ylab = "fits", main = "QQ-plot, correlations")
    lines(seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
          seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
          col = "grey", lty = 2)
    par(op); dev.off()

    # plot 4: Comparison of correlations by heatmap
    if (fileType == "pdf"){
      pdf(paste(dir, "heatmap_correlations_", nameExt, ".pdf"))
    }
    if (fileType == "eps"){
      setEPS(); postscript(paste(dir, "heatmap_correlations_", nameExt, ".eps"))
    }
    op  <- par(pty = "s")
    slh <- cov2cor(S)
    slh[lower.tri(slh)] <- cov2cor(Shat)[lower.tri(Shat)]
    gplot <- edgeHeat(slh, diag = diag, legend = FALSE, main = "Correlations")
    print(gplot); par(op); dev.off()


    print("Visualizing partial correlation fit")
    # If the sample covariance matrix non-singular,
    # also evaluate partial correlation fit
    if (evaluateS(S, verbose = FALSE)$posEigen){

      # plot 5: QQ-plot of partial correlations
      if (fileType == "pdf"){
        pdf(paste(dir, "QQplot_partCorrelations_", nameExt, ".pdf"))
      }
      if (fileType == "eps"){
        setEPS();
        postscript(paste(dir, "QQplot_partCorrelations_", nameExt, ".eps"))
      }
      if (diag){
        cObs <- as.numeric(pcor(solve(S))[upper.tri(S)], diag = TRUE);
        cFit <- as.numeric(pcor(Phat)[upper.tri(Phat, diag = TRUE)])}
      if (!diag){
        cObs <- as.numeric(pcor(solve(S))[upper.tri(S)]);
        cFit <- as.numeric(pcor(Phat)[upper.tri(Phat)])
      }
      op <- par(pty = "s")
      qqplot(x = cObs, y = cFit, pch = 20, xlab = "sample partial correlations",
             ylab = "fits", main = "QQ-plot, partial correlations")
      lines(seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
            seq(min(cFit, cObs), max(cFit, cObs), length.out = 100),
            col = "grey", lty = 2)
      par(op); dev.off()

      # plot 6: Comparison of partial correlations by heatmap
      if (fileType == "pdf"){
        pdf(paste(dir, "heatmap_partCorrelations_", nameExt, ".pdf"))
      }
      if (fileType == "eps"){
        setEPS();
        postscript(paste(dir, "heatmap_partCorrelations_", nameExt, ".eps"))
      }
      op  <- par(pty = "s")
      slh <- pcor(solve(S))
      slh[lower.tri(slh)] <- pcor(Phat)[lower.tri(Phat)]
      gplot <- edgeHeat(slh, diag = diag, legend = FALSE,
                        main = "Partial correlations")
      print(gplot); par(op); dev.off()

    } else {
      print(paste("sample covariance matrix is singular:",
                  "partial correlation fit not visualized"))
    }
  }
}




##------------------------------------------------------------------------------
##
## Functions for Visualization
##
##------------------------------------------------------------------------------

ridgePathS <- function (S, lambdaMin, lambdaMax, step, type = "Alt",
                        target = default.target(S), plotType = "pcor",
                        diag = FALSE, vertical = FALSE, value, verbose = TRUE){
  ##############################################################################
  # - Function that visualizes the regularization path
  # - Regularization path may be visualized for (partial) correlations,
  #   covariances and precision elements
  # - S         > sample covariance/correlation matrix
  # - lambdaMin > minimum value penalty parameter (dependent on 'type')
  # - lambdaMax > maximum value penalty parameter (dependent on 'type')
  # - step      > determines the coarseness in searching the grid
  #               [lambdaMin, lambdaMax]
  # - type      > must be one of {"Alt", "ArchI", "ArchII"}, default = "Alt"
  # - target    > target (precision terms) for Type I estimators,
  #               default = default.target(S)
  # - plotType  > specificies the elements for which the regularization path is
  #               to be visualized.
  #               Must be one of {"pcor", "cor", "cov", "prec"},
  #               default = "pcor"
  # - diag      > logical indicating if the diagonal elements should be retained
  #               for plotting, default = FALSE.
  # - vertical  > optional argument for visualization vertical line in graph
  #               output, default = FALSE
  #               Can be used to indicate the value of, e.g., the optimal
  #               penalty as indicated by some
  #               routine. Can be used to assess the whereabouts of this optimal
  #               penalty along the regularization path.
  # - value     > indicates constant on which to base vertical line when
  #               vertical = TRUE
  # - verbose   > logical indicating if intermediate output should be printed
  #               on screen
  ##############################################################################
  # Dependencies
  # require("base")

  if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  if (verbose){
    cat("Perform input checks...", "\n")
  }
  if (!is.matrix(S)){
    stop("input (S) should be a matrix")
  }
  if (!isSymmetric(S)){
    stop("Input (S) should be a covariance matrix")
  }
  else if (class(lambdaMin) != "numeric"){
    stop("Input (lambdaMin) is of wrong class")
  }
  else if (length(lambdaMin) != 1){
    stop("Input (lambdaMin) must be a scalar")
  }
  else if (lambdaMin <= 0){
    stop("Input (lambdaMin) must be positive")
  }
  else if (class(lambdaMax) != "numeric"){
    stop("Input (lambdaMax) is of wrong class")
  }
  else if (length(lambdaMax) != 1){
    stop("Input (lambdaMax) must be a scalar")
  }
  else if (lambdaMax <= lambdaMin){
    stop("lambdaMax must be larger than lambdaMin")
  }
  else if (class(step) != "numeric") {
    stop("Input (step) is of wrong class")
  }
  else if (!.is.int(step)){
    stop("Input (step) should be integer")
  }
  else if (step <= 0){
    stop("Input (step) should be a positive integer")
  }
  else if (class(plotType) != "character") {
    stop("Input (plotType) is of wrong class")
  }
  else if (!(plotType %in% c("pcor", "cor", "cov", "prec"))){
    stop("Input (plotType) should be one of {'pcor', 'cor', 'cov', 'prec'}")
  }
  else if (length(nchar(plotType)) != 1){
    stop("Input (plotType) should be exactly one of {'pcor', 'cor', 'cov', ",
         "'prec'}")
  }
  if (class(diag) != "logical") {
    stop("Input (diag) is of wrong class")
  }
  else if (class(vertical) != "logical"){
    stop("Input (vertical) is of wrong class")
  }
  else {
    # Set preliminaries
    lambdas  <- seq(lambdaMin, lambdaMax, len = step)
    YforPlot <- numeric()

    # Calculate paths
    if (verbose){cat("Calculating...", "\n")}
    if (type == "Alt" & all(target == 0)){
      if (!isSymmetric(target)){
        stop("Input (target) should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("Inputs ('S' and 'target') should be of the same dimension")
      } else {
        Spectral <- eigen(S, symmetric = TRUE)
        Vinv     <- solve(Spectral$vectors)
        for (k in 1:length(lambdas)){
          Eigshrink <- .eigShrink(Spectral$values, lambdas[k])
          P         <- t(Vinv) %*% diag(1/Eigshrink) %*% Vinv
          if (plotType=="pcor"){
            YforPlot <- cbind(YforPlot, pcor(symm(P))[upper.tri(P)])
          }
          if (plotType=="prec"){
            YforPlot <- cbind(YforPlot, P[upper.tri(P, diag = diag)])
          }
          if (plotType=="cov") {
            YforPlot <- cbind(YforPlot, solve(P)[upper.tri(P, diag = diag)])
          }
          if (plotType=="cor") {
            YforPlot <- cbind(YforPlot, cov2cor(solve(P))[upper.tri(P)])
          }
          if (verbose){
            cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")
          }
        }
      }
    } else if (type == "Alt" & all(target[!diag(nrow(target))] == 0) &
                 (length(unique(diag(target))) == 1)){
      if (!isSymmetric(target)){
        stop("Input (target) should be symmetric")
      } else if (dim(target)[1] != dim(S)[1]){
        stop("Inputs ('S' and 'target') should be of the same dimension")
      } else if (any(diag(target) <= 0)){  ## THIS TEST IS NOT CORRECT I GUESS?
        stop("Input (target) should be p.d.")
      } else {
        varPhi   <- unique(diag(target))
        Spectral <- eigen(S, symmetric = TRUE)
        Vinv     <- solve(Spectral$vectors)
        for (k in 1:length(lambdas)){
          Eigshrink <- .eigShrink(Spectral$values, lambdas[k], const = varPhi)
          P         <- t(Vinv) %*% diag(1/Eigshrink) %*% Vinv
          if (plotType=="pcor"){
            YforPlot <- cbind(YforPlot, pcor(symm(P))[upper.tri(P)])
          }
          if (plotType=="prec"){
            YforPlot <- cbind(YforPlot, P[upper.tri(P, diag = diag)])
          }
          if (plotType=="cov") {
            YforPlot <- cbind(YforPlot, solve(P)[upper.tri(P, diag = diag)])
          }
          if (plotType=="cor") {
            YforPlot <- cbind(YforPlot, cov2cor(solve(P))[upper.tri(P)])
          }
          if (verbose){
            cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")
          }
        }
      }
    } else {
      for (k in 1:length(lambdas)){
        P <- ridgeP(S, lambdas[k], type = type, target = target)
        if (plotType=="pcor"){
          YforPlot <- cbind(YforPlot, pcor(symm(P))[upper.tri(P)])
        }
        if (plotType=="prec"){
          YforPlot <- cbind(YforPlot, P[upper.tri(P, diag = diag)])
        }
        if (plotType=="cov") {
          YforPlot <- cbind(YforPlot, solve(P)[upper.tri(P, diag = diag)])
        }
        if (plotType=="cor") {
          YforPlot <- cbind(YforPlot, cov2cor(solve(P))[upper.tri(P)])
        }
        if (verbose){
          cat(paste("lambda = ", lambdas[k], " done", sep = ""), "\n")
        }
      }
    }

    # Visualize
    if (plotType=="cor") {ylabel <- "penalized correlation"}
    if (plotType=="cov") {ylabel <- "penalized covariances"}
    if (plotType=="pcor"){ylabel <- "penalized partial correlation"}
    if (plotType=="prec"){ylabel <- "penalized precision elements"}
    if (type == "Alt"){Main = "Alternative ridge estimator"}
    if (type == "ArchI"){Main = "Archetypal I ridge estimator"}
    if (type == "ArchII"){Main = "Archetypal II ridge estimator"}

    plot(YforPlot[1,] ~ log(lambdas), axes = FALSE, xlab = "ln(penalty value)",
         ylab = ylabel, main = Main, col = "white",
         ylim = c(min(YforPlot), max(YforPlot)))
    for (k in 1:nrow(YforPlot)){
      lines(YforPlot[k, ] ~ log(lambdas), col = k, lty = k)
    }
    axis(2, col = "black", lwd = 1)
    axis(1, col = "black", lwd = 1)
    if (vertical){
      if (missing(value)){
        stop("Need to specify input (value)")
      } else if (class(value) != "numeric"){
        stop("Input (value) is of wrong class")
      } else if (length(value) != 1){
        stop("Input (value) must be a scalar")
      } else if (value <= 0){
        stop("Input (value) must be positive")
      } else {
        abline(v = log(value), col = "red")
      }
    }
  }
}



if (getRversion() >= "2.15.1") utils::globalVariables(c("X1", "X2", "value"))

edgeHeat <- function(M, lowColor = "blue", highColor = "red", textsize = 10,
                     diag = TRUE, legend = TRUE, main = ""){
  ##############################################################################
  # - function that visualizes precision matrix as a heatmap
  # - can be used to assess (visually) the performance of set of graphical
  #   modeling techniques
  # - M         > Precision matrix
  # - lowColor  > determines color scale in the negative range, default = "blue"
  # - highColor > determines color scale in the positive range, default = "red"
  # - textsize  > set textsize row and column labels, default = 10
  # - diag      > logical determining treatment diagonal elements M. If FALSE,
  #               then the diagonal elements are given the midscale color of
  #               white; only when M is a square matrix
  # - legend    > optional inclusion of color legend, default = TRUE
  # - main      > character specifying the main title, default = ""
  ##############################################################################

  # Dependencies
  #require("ggplot2")
  #require("reshape")

  if (!is.matrix(M)){
    stop("Supply 'M' as matrix")
  }
  else if (class(lowColor) != "character"){
    stop("Input (lowColor) is of wrong class")
  }
  else if (length(lowColor) != 1){
    stop("Length lowColor must be one")
  }
  else if (class(highColor) != "character"){
    stop("Input (highColor) is of wrong class")
  }
  else if (length(highColor) != 1){
    stop("Length highColor must be one")
  }
  else if (class(textsize) != "numeric"){
    stop("Input (textsize) is of wrong class")
  }
  else if (length(textsize) != 1){
    stop("Length textsize must be one")
  }
  else if (textsize <= 0){
    stop("textsize must be positive")
  }
  else if (class(diag) != "logical"){
    stop("Input (diag) is of wrong class")
  }
  else if (class(legend) != "logical"){
    stop("Input (legend) is of wrong class")
  }
  else if (class(main) != "character"){
    stop("Input (main) is of wrong class")
  }
  else {
    # Put matrix in data format
    if (nrow(M) == ncol(M) & !diag) {diag(M) <- 0}
    Mmelt    <- melt(M)
    Mmelt$X1 <-
      factor(as.character(Mmelt$X1), levels = unique(Mmelt$X1), ordered = TRUE)
    Mmelt$X2 <-
      factor(as.character(Mmelt$X2), levels = unique(Mmelt$X2), ordered = TRUE)

    # Visualize
    if (legend){
      ggplot(Mmelt, aes(X2, X1, fill = value)) + geom_tile() +
        scale_fill_gradient2("", low = lowColor,  mid = "white",
                             high = highColor, midpoint = 0) +
        theme(axis.ticks = element_blank()) +
        theme(axis.text.y = element_text(size = textsize)) +
        theme(axis.text.x =
                element_text(angle = -90, vjust = .5, size = textsize)) +
        xlab(" ") + ylab(" ") +
        ylim(rev(levels(Mmelt$X1))) +
        ggtitle(main)
    } else {
      ggplot(Mmelt, aes(X2, X1, fill = value)) + geom_tile() +
        scale_fill_gradient2("", low = lowColor,  mid = "white",
                             high = highColor, midpoint = 0) +
        theme(axis.ticks = element_blank()) +
        theme(axis.text.y = element_text(size = textsize)) +
        theme(axis.text.x =
                element_text(angle = -90, vjust = .5, size = textsize)) +
        xlab(" ") + ylab(" ") +
        ylim(rev(levels(Mmelt$X1))) +
        ggtitle(main) +
        theme(legend.position = "none")
    }
  }
}



Ugraph <- function(M, type = c("plain", "fancy", "weighted"),
                   lay = layout.circle, Vsize = 15, Vcex = 1,
                   Vcolor = "orangered", VBcolor = "darkred",
                   VLcolor = "black", prune = FALSE, legend = FALSE,
                   label = "", Lcex = 1.3, PTcex = 4, cut = .5,
                   scale = 10, pEcolor = "black", nEcolor = "grey",
                   main = ""){
  ##############################################################################
  # - Function that visualizes the sparsified precision matrix as an undirected
  #   graph
  # - Function is partly a wrapper around certain 'igraph' functions
  # - M       > (Possibly sparsified) precision matrix
  # - type    > graph type: 'plain' gives plain undirected graph. 'fancy' gives
  #             undirected graph in which dashed lines indicate negative partial
  #             correlations while solid lines indicate positive partial
  #             correlations, and in which black lines indicate strong edges.
  #             'weighted' gives an undirected graph in which edge thickness
  #             indicates the strenght of the partial correlations. Grey lines
  #             then indicate negative partial correlations while black lines
  #             represent positive partial correlations.
  # - lay     > determines layout of the graph. All layouts in 'layout{igraph}'
  #             are accepted. Default = layout.circle, giving circular layout.
  # - Vsize   > gives vertex size, default = 15
  # - Vcex    > gives size vertex labels, default = 1
  # - Vcolor  > gives vertex color, default = "orangered"
  # - VBcolor > gives color of the vertex border, default = "darkred"
  # - VLcolor > gives color of the vertex labels, default = "black"
  # - prune   > logical indicating if vertices of degree 0 should be removed
  # - legend  > optional inclusion of color legend, default = FALSE
  # - label   > character label for the endogenous variables, default = "";
  #             only when legend = TRUE
  # - Lcex    > scaling legend box, default = 1.3; only when legend = TRUE
  # - PTcex   > scaling node in legend box, default = 4; only when legend = TRUE
  # - cut     > cut-off for indication of strong edge, default = .5; only when
  #             type = "fancy"
  # - scale   > scale factor for visualizing strenght of edges, default = 10;
  #             only when type = "weighted"
  # - pEcolor > gives edge color for edges tied to positive precision elements,
  #             default = "black"; only when type = "weighted"
  # - nEcolor > gives edge color for edges tied to negative precision elements,
  #             default = "grey"; only when type = "weighted"
  # - main    > character specifying heading figure, default = ""
  ##############################################################################

  # Dependencies
  # require("igraph")
  # require("reshape")

  if (!is.matrix(M)){
    stop("M should be a matrix")
  }
  else if (nrow(M) != ncol(M)){
    stop("M should be square matrix")
  }
  else if (missing(type)){
    stop("Need to specify graph type ('plain' or 'fancy' or 'weighted')")
  }
  else if (!(type %in% c("plain", "fancy", "weighted"))){
    stop("type should be one of {'plain', 'fancy', 'weighted'}")
  }
  else if (class(Vsize) != "numeric"){
    stop("Input (Vsize) is of wrong class")
  }
  else if (length(Vsize) != 1){
    stop("Length Vsize must be one")
  }
  else if (Vsize <= 0){
    stop("Vsize must be positive")
  }
  else if (class(Vcex) != "numeric"){
    stop("Input (Vcex) is of wrong class")
  }
  else if (length(Vcex) != 1){
    stop("Length Vcex must be one")
  }
  else if (Vcex <= 0){
    stop("Vcex must be positive")
  }
  else if (class(Vcolor) != "character"){
    stop("Input (Vcolor) is of wrong class")
  }
  else if (length(Vcolor) != 1){
    stop("Length Vcolor must be one")
  }
  else if (class(VBcolor) != "character"){
    stop("Input (VBcolor) is of wrong class")
  }
  else if (length(VBcolor) != 1){
    stop("Length VBcolor must be one")
  }
  else if (class(VLcolor) != "character"){
    stop("Input (VLcolor) is of wrong class")
  }
  else if (length(VLcolor) != 1){
    stop("Length VLcolor must be one")
  }
  else if (class(prune) != "logical"){
    stop("Input (prune) is of wrong class")
  }
  else if (class(legend) != "logical"){
    stop("Input (legend) is of wrong class")
  }
  else if (class(main) != "character"){
    stop("Input (main) is of wrong class")
  }
  else {
    # Preliminaries
    AM <- adjacentMat(M)
    GA <- graph.adjacency(AM, mode = "undirected")
    if (prune){GA <- delete.vertices(GA, which(degree(GA) < 1))}

    # Plain graph
    if (type == "plain"){
      plot(GA, layout = lay, vertex.size = Vsize, vertex.label.family = "sans",
           vertex.label.cex = Vcex, vertex.color = Vcolor,
           vertex.frame.color = VBcolor,
           vertex.label.color = VLcolor, main = main)
    }

    # Fancy graph
    if (type == "fancy"){
      if (class(cut) != "numeric"){
        stop("Input (cut) is of wrong class")
      } else if (length(cut) != 1){
        stop("Length cut must be one")
      } else if (cut <= 0){
        stop("cut must be positive")
      } else {
        Names <- colnames(M)
        colnames(M) = rownames(M) <- seq(1, ncol(M), by = 1)
        Mmelt <- melt(M)
        Mmelt <- Mmelt[Mmelt$X1 > Mmelt$X2,]
        Mmelt <- Mmelt[Mmelt$value != 0,]
        E(GA)$weight <- Mmelt$value
        E(GA)$color  <- "grey"
        E(GA)[E(GA)$weight < 0]$style <- "dashed"
        E(GA)[E(GA)$weight > 0]$style <- "solid"
        E(GA)[abs(E(GA)$weight) > cut]$color <- "black"
        plot(GA, layout = lay, vertex.size = Vsize,
             vertex.label.family = "sans", vertex.label.cex = Vcex,
             vertex.color = Vcolor, vertex.frame.color = VBcolor,
             vertex.label.color = VLcolor,
             edge.color = E(GA)$color, edge.lty = E(GA)$style, main = main)
      }
    }

    # Weighted graph
    if (type == "weighted"){
      if (class(scale) != "numeric"){
        stop("Input (scale) is of wrong class")
      } else if (length(scale) != 1){
        stop("Length scale must be one")
      } else if (scale <= 0){
        stop("scale must be positive")
      } else if (class(pEcolor) != "character"){
        stop("Input (pEcolor) is of wrong class")
      } else if (length(pEcolor) != 1){
        stop("Length pEcolor must be one")
      } else if (class(nEcolor) != "character"){
        stop("Input (nEcolor) is of wrong class")
      } else if (length(nEcolor) != 1){
        stop("Length nEcolor must be one")
      } else {
        Names <- colnames(M)
        colnames(M) = rownames(M) <- seq(1, ncol(M), by = 1)
        Mmelt <- melt(M)
        Mmelt <- Mmelt[Mmelt$X1 > Mmelt$X2,]
        Mmelt <- Mmelt[Mmelt$value != 0,]
        E(GA)$weight <- Mmelt$value
        E(GA)[E(GA)$weight < 0]$color <- nEcolor
        E(GA)[E(GA)$weight > 0]$color <- pEcolor
        plot(GA, layout = lay, vertex.size = Vsize,
             vertex.label.family = "sans", vertex.label.cex = Vcex,
             vertex.color = Vcolor, vertex.frame.color = VBcolor,
             vertex.label.color = VLcolor,
             edge.color = E(GA)$color, edge.width = scale*abs(E(GA)$weight),
             main = main)
      }
    }

    # Legend
    if (legend){
      if (class(label) != "character"){
        stop("Input (label) is of wrong class")
      } else if (length(label) != 1){
        stop("Length label must be one")
      } else if (class(Lcex) != "numeric"){
        stop("Input (Lcex) is of wrong class")
      } else if (length(Lcex) != 1){
        stop("Length Lcex must be one")
      } else if (Lcex <= 0){
        stop("Lcex must be positive")
      } else if (class(PTcex) != "numeric"){
        stop("Input (PTcex) is of wrong class")
      } else if (length(PTcex) != 1){
        stop("Length PTcex must be one")
      } else if (PTcex <= 0){
        stop("PTcex must be positive")
      } else{
        legend("bottomright", label, pch = 20, col = Vcolor,
               cex = Lcex, pt.cex = PTcex)
      }
    }
  }
}




##------------------------------------------------------------------------------
##
## Functions for Topology Statistics
##
##------------------------------------------------------------------------------

GGMnetworkStats <- function(sparseP, as.table = FALSE){
  ##############################################################################
  # - Function that calculates various network statistics from a sparse matrix
  # - Input matrix is assumed to be a sparse precision of partial correlation
  #   matrix
  # - The sparse precision matrix is taken to represent a conditional
  #   independence graph
  # - sparseP  > sparse precision/partial correlation matrix
  # - as.table > logical indicating if output should be returned as table;
  #              default = FALSE
  #
  # - NOTES (network statistics produced):
  # - Node degree
  # - Betweenness centrality
  # - Closeness centrality
  # - Eigenvalue centrality
  # - Number of negative edges for each node
  # - Number of positive edges for each node
  # - Assessment if network/graph is chordal (triangulated)
  # - Mutual information of each variate with all other variates
  # - Variance of each variate (based on inverse sparsified precision matrix)
  # - Partial variance of each variate (= 1 when input matrix is partial
  #   correlation matrix)
  # - Future versions of this function may include additional statistics
  #
  # - REFERENCE:
  # - Newman, M.E.J. (2010), "Networks: an introduction",
  #   Oxford University Press
  ##############################################################################

  # Dependencies
  # require("base")
  # require("igraph")

  if (!is.matrix(sparseP)){
    stop("Input (sparseP) should be a matrix")
  }
  else if (!isSymmetric(sparseP)){
    stop("Input (sparseP) should be a symmetric matrix")
  }
  else if (!evaluateS(sparseP, verbose = FALSE)$posEigen){
    stop("Input (sparseP) is expected to be positive definite")
  }
  else if (class(as.table) != "logical"){
    stop("Input (as.table) is of wrong class")
  }
  else{
    # Some warnings
    if (all(sparseP != 0)){
      warning("Given input (sparseP) implies a saturated conditional ",
              "independence graph")
    }
    if (all(sparseP[!diag(nrow(sparseP))] == 0)){
      warning("Given input (sparseP) implies an empty conditional ",
              "independence graph")
    }

    # Obtain corresponding sample covariance matrix
    pvars <- 1/diag(sparseP)
    S     <- solve(sparseP)

    # Calculate nodes' mutual information
    MI <-
      unlist(lapply(1:nrow(S),
                    function(j, S){
                      log(det(S[-j,-j])) -
                        log(det(S[-j,-j] - S[-j,j,drop=FALSE] %*%
                                  S[j,-j,drop=FALSE]/S[j,j]))
                    }, S = S))
    names(MI) <- colnames(sparseP)

    # Signs of edges
    diag(sparseP) <- 0
    nPos <- apply(sign(sparseP), 2, function(Z){ sum(Z == 1) })
    nNeg <- apply(sign(sparseP), 2, function(Z){ sum(Z == -1) })

    # Adjacency to graphical object
    AM  <- adjacentMat(sparseP)
    CIG <- graph.adjacency(AM, mode = "undirected")

    # Return
    if (as.table){
      networkStats <- cbind(degree(CIG), betweenness(CIG), closeness(CIG),
                            evcent(CIG)$vector, nNeg, nPos, MI, diag(S), pvars)
      colnames(networkStats) <-
        c("degree", "betweenness", "closeness", "eigenCentrality", "nNeg",
          "nPos", "mutualInfo", "variance", "partialVar")
      return(networkStats)
    }
    if (!as.table){
      return(list(degree = degree(CIG), betweenness = betweenness(CIG),
                  closeness = closeness(CIG),
                  eigenCentrality = evcent(CIG)$vector,
                  nNeg = nNeg, nPos = nPos, chordal = is.chordal(CIG)$chordal,
                  mutualInfo = MI, variance = diag(S), partialVar = pvars))
    }
  }
}



GGMpathStats <- function(P0, node1, node2, neiExpansions = 2, verbose = TRUE,
                         graph = TRUE, nrPaths = 2, lay = layout.circle,
                         nodecol = "skyblue", Vsize = 15, Vcex = .6,
                         VBcolor = "darkblue", VLcolor = "black",
                         all.edges = TRUE, prune = TRUE, legend = TRUE,
                         scale = 1, Lcex = .8, PTcex = 2, main = ""){
  ##############################################################################
  # - Function that expresses the covariance between a pair of variables as a
  #   sum of path weights
  # - The sum of path weights is based on the shortest paths connecting the pair
  #   in an undirected graph
  # - P0            > sparse precision/partial correlation matrix
  # - node1         > start node of the path
  # - node2         > end node of the path
  # - neiExpansions > a numeric determining how many times the neighborhood
  #                   around the start and end node should be expanded in the
  #                   search for shortest paths between the node pair.
  #                   Default = 2
  # - verbose       > logical indicating if output should also be printed on
  #                   screen. Default = TRUE
  # - graph         > Optional argument for visualization strongest paths,
  #                   default = TRUE
  # - nrPaths	      > indicates the number of paths with the highest
  #                   contribution to the marginal covariance
  #                   between the indicated node pair (node1 and node2) to be
  #                   visualized/highlighted;
  #                   only when graph = TRUE
  # - lay           > determines layout of the graph. All layouts in
  #                   'layout{igraph}' are accepted. Default = layout.circle,
  #                    giving circular layout; only when graph = TRUE
  # - nodecol       > gives color of node1 and node2; only when graph = TRUE
  # - Vsize   	    > gives vertex size, default = 15; only when graph = TRUE
  # - Vcex    	    > gives size vertex labels, default = .6; only when
  #                   graph = TRUE
  # - VBcolor     	> gives color of the vertex border, default = "darkblue";
  #                   only when graph = TRUE
  # - VLcolor      	> gives color of the vertex labels, default = "black";
  #                   only when graph = TRUE
  # - all.edges     > logical indicating if edges other than those implied by
  #                   the 'nrPaths' paths between
  #                   node1 and node2 should also be visualized. Default = TRUE;
  #                   only when graph = TRUE
  # - prune         > logical indicating if vertices of degree 0 should be
  #                   removed. Default = TRUE; only when graph = TRUE
  # - legend        > optional inclusion of color legend, default = TRUE; only
  #                   when graph = TRUE
  # - scale         > scale factor for visualizing strenght of edges,
  #                   default = 1. It is a relative scaling
  #                   factor, in the sense that the edges implied by the
  #                   'nrPaths' paths between node1 and node2 have edge
  #                   thickness that is twice this scaling factor (so it is
  #                   a scaling factor vis-a-vis the unimplied edges); only
  #                   when all.edges = TRUE
  # - Lcex          > scaling legend box, default = .8; only when legend = TRUE
  # - PTcex         > scaling node in legend box, default = 2; only when
  #                   legend = TRUE
  # - main          > character specifying heading figure, default = ""
  #
  # - NOTES:
  # - As in Jones & West (2005), paths whose weights have an opposite sign to
  #   the marginal covariance (between endnodes of the path) are referred to
  #   as 'moderating paths' while paths whose weights have the same sign as the
  #   marginal covariance are referred to as 'mediating' paths
  ##############################################################################

  # Dependencies
  # require("base")
  # require("igraph")
  # require("reshape")

  if (!is.matrix(P0)){
    stop("Input (P0) should be a matrix")
  }
  else if (!isSymmetric(P0)){
    stop("Input (P0) should be a symmetric matrix")
  }
  else if (!evaluateS(P0, verbose = FALSE)$posEigen){
    stop("Input (P0) is expected to be positive definite")
  }
  else if (class(node1) != "numeric"){
    stop("Input (node1) is of wrong class")
  }
  else if (length(node1) != 1){
    stop("Length input (node1) must be 1")
  }
  else if (!.is.int(node1)){
    stop("Input (node1) should be a numeric integer")
  }
  else if (node1 < 1){
    stop("Input (node1) cannot be zero or negative")
  }
  else if (node1 > ncol(P0)){
    stop("Input (node1) cannot exceed the number of variables in P0")
  }
  else if (class(node2) != "numeric"){
    stop("Input (node2) is of wrong class")
  }
  else if (length(node2) != 1){
    stop("Length input (node2) must be 1")
  }
  else if (!.is.int(node2)){
    stop("Input (node2) should be a numeric integer")
  }
  else if (node2 < 1){
    stop("Input (node2) cannot be zero or negative")
  }
  else if (node2 > ncol(P0)){
    stop("Input (node2) cannot exceed the number of variables in P0")
  }
  else if (node1 == node2){
    stop("Inputs (node1 and node2) cannot be equal")
  }
  else if (class(neiExpansions) != "numeric"){
    stop("Input (neiExpansions) is of wrong class")
  }
  else if (length(neiExpansions) != 1){
    stop("Length input (neiExpansions) must be 1")
  }
  else if (!.is.int(neiExpansions)){
    stop("Input (neiExpansions) should be a numeric integer")
  }
  else if (neiExpansions < 1){
    stop("Input (neiExpansions) cannot be zero or negative")
  }
  else if (class(graph) != "logical"){
    stop("Input (graph) is of wrong class")
  }
  else if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  else {
    # Some warnings
    if (all(P0 != 0)){
      warning("Given input (P0) implies a saturated conditional independence ",
              "graph")
    }
    if (all(P0[!diag(nrow(P0))] == 0)){
      warning("Given input (P0) implies an empty conditional independence ",
              "graph")
    }

    # Precision associated graph
    colnames(P0) <- colnames(P0, do.NULL = FALSE)
    Names  <- colnames(P0)
    adjMat <- adjacentMat(P0)
    colnames(adjMat) <- 1:nrow(P0)
    rownames(adjMat) <- 1:nrow(P0)
    G <- graph.adjacency(adjMat, mode = "undirected")

    # Is there a connection between the specified nodes?
    pathExits <- is.finite(shortest.paths(G, node1, node2))

    # If nodes are connected, evaluate path contributions
    if (!pathExits){
      stop("provided node pair is not connected.")
    } else {
      # Determinant of precision
      detP0 <- det(P0)

      # Objects to be returned
      paths <- list()
      pathStats <- numeric()

      # Shortest paths between the nodes
      slh <- get.all.shortest.paths(G, node1, node2)$res
      for (u in 1:length(slh)){
        fullPath <- slh[[u]]
        pName <- .path2string(fullPath)
        paths[[length(paths)+1]] <- fullPath
        pathStats <- rbind(pathStats, c(length(slh[[u]])-1,
                                        .pathContribution(P0, fullPath, detP0)))
        rownames(pathStats)[nrow(pathStats)] <- pName
      }
      nei1 <- node1
      nei2 <- node2

      for (u in 1:neiExpansions){
        # Consider longer paths between the nodes
        nei1temp <- nei1
        nei2temp <- nei2
        for (v1 in 1:length(nei1)){
          nei1temp <- c(nei1temp, neighbors(G, nei1[v1]))
        }
        for (v2 in 1:length(nei2)){
          nei2temp <- c(nei2temp, neighbors(G, nei2[v2]))
        }
        nei1 <- setdiff(unique(nei1temp), node1)
        nei2 <- setdiff(unique(nei2temp), node2)
        slh  <- .pathAndStats(G, node1, node2, nei1, nei2, P0, detP0,
                              rownames(pathStats))
        pathStats <- rbind(pathStats, slh$pathStats)
        paths <- c(paths, slh$paths)
      }

      # Wrap up
      pNames <- rownames(pathStats)
      paths  <- paths[order(abs(pathStats[,2]), decreasing=TRUE)]
      pathStats <- matrix(pathStats[order(abs(pathStats[,2]),
                                          decreasing=TRUE),], ncol = 2)
      names(paths) = rownames(pathStats) <- pNames
      colnames(pathStats) <- c("length", "contribution")
      if (verbose | graph){covNo1No2 <- solve(P0)[node1, node2]}

      # Summary
      if (verbose){
        covNo1No2expl <- sum(pathStats[,2])

        # Reformat results
        statsTable <- data.frame(cbind(rownames(pathStats), pathStats[,1],
                                       round(pathStats[,2], 5)))
        rownames(statsTable) <- NULL
        colnames(statsTable) <- c("path", "length", "contribution")

        # Print results on screen
        cat("Covariance between node pair :", round(covNo1No2, 5), "\n")
        cat("----------------------------------------\n")
        print(statsTable, quote=FALSE)
        cat("----------------------------------------\n")
        cat("Sum path contributions       :", round(covNo1No2expl, 5), "\n")

      }

      # Visualize
      if (graph){
        if (class(nrPaths) != "numeric"){
          stop("Input (nrPaths) is of wrong class")
        } else if (length(nrPaths) != 1){
          stop("Length input (nrPaths) must be one")
        } else if (!.is.int(nrPaths)){
          stop("Input (nrPaths) should be a numeric integer")
        } else if (nrPaths <= 0){
          stop("Input (nrPaths) must be a strictly positive integer")
        } else if (nrPaths > length(paths)){
          stop("Input (nrPaths) cannot exceed the total number of paths ",
               "discerned")
        } else if (class(nodecol) != "character"){
          stop("Input (nodecol) is of wrong class")
        } else if (length(nodecol) != 1){
          stop("Length input (nodecol) must be one")
        } else if (class(Vsize) != "numeric"){
          stop("Input (Vsize) is of wrong class")
        } else if (length(Vsize) != 1){
          stop("Length input (Vsize) must be one")
        } else if (Vsize <= 0){
          stop("Input (Vsize) must be strictly positive")
        } else if (class(Vcex) != "numeric"){
          stop("Input (Vcex) is of wrong class")
        } else if (length(Vcex) != 1){
          stop("Length input (Vcex) must be one")
        } else if (Vcex <= 0){
          stop("Input (Vcex) must be strictly positive")
        } else if (class(VBcolor) != "character"){
          stop("Input (VBcolor) is of wrong class")
        } else if (length(VBcolor) != 1){
          stop("Length input (VBcolor) must be one")
        } else if (class(VLcolor) != "character"){
          stop("Input (VLcolor) is of wrong class")
        } else if (length(VLcolor) != 1){
          stop("Length input (VLcolor) must be one")
        } else if (class(all.edges) != "logical"){
          stop("Input (all.edges) is of wrong class")
        } else if (class(prune) != "logical"){
          stop("Input (prune) is of wrong class")
        } else if (class(legend) != "logical"){
          stop("Input (legend) is of wrong class")
        } else if (class(main) != "character"){
          stop("Input (main) is of wrong class")
        } else {
          # Preliminaries
          AM <- adjacentMat(P0)
          GA <- graph.adjacency(AM, mode = "undirected")
          colnames(P0) = rownames(P0) <- seq(1, ncol(P0), by = 1)
          Mmelt <- melt(P0)
          Mmelt <- Mmelt[Mmelt$X1 > Mmelt$X2,]
          Mmelt <- Mmelt[Mmelt$value != 0,]

          # Determine if path is mediating or moderating and color accordingly
          for (k in 1:nrPaths){
            Path <- unlist(paths[k])
            if (sign(covNo1No2) == sign(pathStats[k,2])){COL = "green"}
            if (sign(covNo1No2) != sign(pathStats[k,2])){COL = "red"}
            for (i in 1:(length(Path) - 1)){
              if (Path[i] > Path[i + 1]){
                tempX1 <- Path[i]
                tempX2 <- Path[i + 1]
              } else {
                tempX1 <- Path[i + 1]
                tempX2 <- Path[i]
              }
              row <- which(Mmelt$X1 == tempX1 & Mmelt$X2 == tempX2)
              E(GA)[row]$color <- COL
            }
          }

          # Coloring nodes
          V(GA)$color <- "white"
          V(GA)$color[node1] <- nodecol
          V(GA)$color[node2] <- nodecol

          # Produce graph
          if (all.edges){
            if (class(scale) != "numeric"){
              stop("Input (scale) is of wrong class")
            } else if (length(scale) != 1){
              stop("Length input (scale) must be one")
            } else if (scale <= 0){
              stop("Input (scale) must be strictly positive")
            } else {
              E(GA)[is.na(E(GA)$color)]$color <- "grey"
              E(GA)$weight <- 1
              E(GA)[E(GA)$color == "green"]$weight <- 2
              E(GA)[E(GA)$color == "red"]$weight   <- 2
              if (prune){GA <- delete.vertices(GA, which(degree(GA) < 1))}
              plot(GA, layout = lay, vertex.size = Vsize,
                   vertex.label.family = "sans", vertex.label.cex = Vcex,
                   edge.width = scale*abs(E(GA)$weight),
                   vertex.color = V(GA)$color, vertex.frame.color = VBcolor,
                   vertex.label.color = VLcolor, main = main)
            }
          } else {
            if (prune){GA <- delete.vertices(GA, which(degree(GA) < 1))}
            plot(GA, layout = lay, vertex.size = Vsize,
                 vertex.label.family = "sans", vertex.label.cex = Vcex,
                 vertex.color = V(GA)$color, vertex.frame.color = VBcolor,
                 vertex.label.color = VLcolor, main = main)
          }

          # Legend
          if (legend){
            if (class(Lcex) != "numeric"){
              stop("Input (Lcex) is of wrong class")
            } else if (length(Lcex) != 1){
              stop("Length input (Lcex) must be one")
            } else if (Lcex <= 0){
              stop("Input (Lcex) must be strictly positive")
            } else if (class(PTcex) != "numeric"){
              stop("Input (PTcex) is of wrong class")
            } else if (length(PTcex) != 1){
              stop("Length input (PTcex) must be one")
            } else if (PTcex <= 0){
              stop("Input (PTcex) must be strictly positive")
            } else {
              legend("bottomright", c("mediating path", "moderating path"),
                     lty=c(1,1), col = c("green", "red"), cex = Lcex,
                     pt.cex = PTcex)
            }
          }
        }
      }

      # Return
      Numeric    <- rownames(adjMat)
      VarName    <- Names
      identifier <- data.frame(Numeric, VarName)
      return(list(pathStats = pathStats, paths = paths,
                  Identifier = identifier))
    }
  }
}




##------------------------------------------------------------------------------
##
## Wrapper function
##
##------------------------------------------------------------------------------

fullMontyS <- function(Y, lambdaMin, lambdaMax,
                       target = default.target(covML(Y)), dir = getwd(),
                       fileTypeFig = "pdf", FDRcut = .9, nOutput = TRUE,
                       verbose = TRUE){
  ##############################################################################
  # - Function that forms a wrapper around the rags2ridges functionalities
  # - Invokes functionalities to get from data to graph and topology summaries
  # - Y           > (raw) Data matrix, variables in columns
  # - lambdaMin   > minimum value penalty parameter
  # - lambdaMax   > maximum value penalty parameter
  # - target      > target (precision terms) for Type I estimators,
  #                 default = default.target(covML(Y))
  # - dir         > specifies the directory in which the (visual) output is
  #                 stored
  # - fileTypeFig > signifies filetype of visual output; Should be one of
  #                 {"pdf", "eps"}
  # - FDRcut      > cut-off for partial correlation element selection based on
  #                 local FDR thresholding. Default = .9.
  # - nOutput     > logical indicating if numeric output should be given
  # - verbose     > logical indicating if intermediate output should be printed
  #                 on screen
  #
  # - NOTES:
  # - Always uses the alternative ridge estimator by van Wieringen and Peeters
  #   (2015)
  # - Always uses LOOCV by Brent for optimal penalty parameter determination
  # - Always uses support determination by local FDR thresholding
  # - Visualizes the network by 'Ugraph' with type = "fancy". Network on basis
  #   partial correlations
  # - Network statistics calculated on sparsified partial correlation network
  # - There are no elaborate input checks as these are all covered by the
  #   indvidual functions invoked
  ##############################################################################

  # Dependencies
  # require("base")
  # require("stats")
  # require("graphics")
  # require("Hmisc")
  # require("fdrtool")
  # require("igraph")
  # require("reshape")
  # require("utils")

  if (class(dir) != "character"){
    stop("Specify directory for output (dir) as 'character'")
  }
  else if (!(fileTypeFig %in% c("pdf", "eps"))){
    stop("Input (fileTypeFig) should be one of {'pdf', 'eps'}")
  }
  else if (class(nOutput) != "logical"){
    stop("Input (nOutput) is of wrong class")
  }
  else if (class(verbose) != "logical"){
    stop("Input (verbose) is of wrong class")
  }
  else {
    # In case one does not know:
    cat("Output files are stored in the following directory:", dir, "\n")
    cat("\n")

    # Determine optimal penalty and precision
    if (verbose){cat("Progress:", "\n")}
    if (verbose){cat("Determining optimal penalty value...", "\n")}
    optimal <- optPenalty.LOOCVauto(Y, lambdaMin = lambdaMin,
                                    lambdaMax = lambdaMax,
                                    lambdaInit = (lambdaMin+lambdaMax)/2,
                                    target = target)

    # Condition number plot
    if (verbose){cat("Generating condition number plot...", "\n")}
    if (fileTypeFig == "pdf"){pdf(paste(dir, "Condition_Number_Plot.pdf"))}
    if (fileTypeFig == "eps"){
      setEPS()
      postscript(paste(dir, "Condition_Number_Plot.eps"))
    }
    conditionNumberPlot(covML(Y), lambdaMin = lambdaMin, lambdaMax = lambdaMax,
                        step = 100000, target = target, vertical = TRUE,
                        value = optimal$optLambda, main = FALSE, verbose =FALSE)
    dev.off()

    # Sparsify the precision matrix
    if (verbose){cat("Determining support...", "\n")}
    PC0 <- sparsify(optimal$optPrec, threshold = "localFDR", FDRcut = FDRcut,
                    output = "heavy", verbose = FALSE)$sparseParCor

    # Visualize the network
    if (verbose){cat("Visualizing network...", "\n")}
    if (fileTypeFig == "pdf"){pdf(paste(dir, "Network.pdf"))}
    if (fileTypeFig == "eps"){setEPS(); postscript(paste(dir, "Network.eps"))}
    Ugraph(PC0, type = "fancy", Vsize = ncol(PC0)/10, Vcex = ncol(PC0)/130)
    dev.off()

    # Calculate network statistics
    if (verbose){cat("Calculating network statistics...", "\n")}
    Stats <- GGMnetworkStats(PC0, as.table = TRUE)
    capture.output(print(Stats), file = paste(dir, "Network_Statistics.txt"))

    # Done
    if (verbose){cat("DONE!", "\n")}

    # Return
    if (nOutput){
      return(list(optLambda = optimal$optLambda, optPrec = optimal$optPrec,
                  sparseParCor = PC0, networkStats = Stats))
    }
  }
}






################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module B: rags2ridges Fused
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

# See R/rags2ridgesFused.R






################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module C: rags2ridges Miscellaneous
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

# See R/rags2ridgesMisc.R






################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module D: rags2ridges DCMG
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

# In development






################################################################################
################################################################################
##------------------------------------------------------------------------------
##
## Module E: rags2ridges Robust
##
##------------------------------------------------------------------------------
################################################################################
################################################################################

# In development


