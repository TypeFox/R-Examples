#############################################
#### AUTHOR:    Arnost Komarek           ####
####            (2003)                   ####
####                                     ####
#### FILE:      smoothSurvReg.control.R  ####
####                                     ####
#### FUNCTIONS: smoothSurvReg.control    ####
#############################################

### ====================================================================================
### smoothSurvReg.coontrol: More options for smoothSurvReg function
### ====================================================================================
## est.c
## est.scale
## maxiter ............ maximal number of Newton-Raphson iterations
## firstiter .......... number of the first iteration (useful when not starting
##                      iterations from the beginning)
## rel.tolerance ...... tolerance for the convergence (norm of the appropriate score vector)
## toler.chol ......... tolerance for the Cholesky decomposition to detect
##                      non positive definite matrices
## toler.eigen ........ value used in eigen value decomposition to replace
##                      eigen values too close to zero or negative eigen values
## maxhalf ............ maximal number of step-halving steps before reporting
##                      non-convergence
## debug .............. do I want to debug?
## info ............... do I want to print some information during the iterations
## lambda.use.......... tuning parameter for the penalty used in a given fit
## sdspline ........... standard deviation of one basis spline
## difforder .......... order of the difference used in the penalty
## dist.range ......... approximate range of knots
## by.knots ........... distance between the two knots
## knots .............. knots
## nsplines
## last.three ......... indeces of the knots which are to be functions of the remaining ones
##                      (applicable if c's are estimated)
smoothSurvReg.control <- function(
                              est.c = TRUE,
                              est.scale = TRUE,
                              maxiter = 200,
                              firstiter = 0,
                              rel.tolerance = 5e-5,
                              toler.chol = 1e-15,
                              toler.eigen = 1e-3,
                              maxhalf = 10,
                              debug = 0,
                              info = TRUE,
                              lambda.use = 1.0,
                              sdspline = NULL,
                              difforder = 3,
                              dist.range = c(-6, 6),
                              by.knots = 0.3,
                              knots = NULL,
                              nsplines = NULL,
                              last.three = NULL
                           )
{

  if (length(dist.range) != 2) stop("Invalid 'dist.range' ")
  if (dist.range[2] < dist.range[1]) stop("Invalid 'dist.range' ")
  if ((dist.range[1] >= 0) || (dist.range[2] <= 0)) stop("Invalid 'dist.range' ")
  if (lambda.use < 0) stop("Penalty 'lambda.use' parameter has to be positive ")
  if (by.knots <= 0) stop("Distance between the two knots must be positive ")

  ## Knots
  if (is.null(knots)){
     middle.knot <- 0
     between.knots <- by.knots
     knots1 <- seq(0, dist.range[2], by = by.knots)
     knots2 <- seq(0, dist.range[1], by = -by.knots)
     knots2 <- knots2[-1]                ## remove the first zero
     knots2 <- knots2[order(knots2)]
     knots <- c(knots2, knots1)
     nknot <- length(knots)
     nsplines <- nknot
  }
  else {
     nknot <- length(knots)
     if ((est.c) && (length(unique(knots)) != nknot)) stop("All knots have to be distinct ")
     knots <- knots[order(knots)]
     nsplines <- nknot
  }

  if (nsplines <= 3) est.c <- FALSE

  ## SD spline
  if (nsplines > 1) between.knots <- max(knots[2:nsplines] - knots[1:(nsplines-1)])
  else              between.knots <- 3/2
  if (is.null(sdspline))  sdspline <- (2/3) * between.knots
  else if (sdspline <= 0) sdspline <- (2/3) * between.knots
  if ((sdspline >= 1) && (nsplines > 1) && est.c){
      warning("sdspline higher than 1 changed into 0.9 ")
      sdspline <- 0.9
  }

  ## last.three
  if (est.c){     # there are at least 4 knots
     if (is.null(last.three)){
        which.zero <- which.min(abs(knots))
        if (which.zero > 1 && which.zero < nsplines)
           last.three <- c(which.zero, which.zero - 1, which.zero + 1)
        else
           if (which.zero == 1) last.three <- 1:3
           else                 last.three <- nsplines:(nsplines-2)
     }
     if (length(last.three) != 3) stop("Incorrect 'last.three' parameter ")
     if (sum(last.three %in% (1:nsplines)) != 3) stop("Incorrect 'last.three' parameter ")
     if (length(unique(last.three)) != 3) stop("Incorrect 'last.three' parameter ")
     if (abs(knots[last.three[2]]) <  1e-4)
            stop("Zero reference knot[last.three[2]] ")
     if (abs(knots[last.three[2]]-knots[last.three[3]]) <  1e-4)
            stop("Too close reference knots[last.three[2]] and knots[last.three[3]] ")
     if (abs(1 - sdspline + knots[last.three[2]]*knots[last.three[3]]) < 1e-4)
            stop("Badly conditioned reference knots[last.three[2]] and knots[last.three[3]] ")
  }
  else
     last.three <- 1:3

  ## Order of the difference in the penalty
  if (est.c){
    if ((difforder < 0) || (difforder > nsplines - 1))
         stop("'difforder' has to be non-negative and smaller than 'nsplines'  \nDefault value of 'difforder' is 2!")
  }
  else
    difforder <- 0

  ## Check some values
  if (maxiter < 0) stop("'maxiter' has to be non-negative ")
  if (rel.tolerance <= 0) stop("'rel.tolerance' has to be positive ")
  if (toler.chol <= 0) stop("'toler.chol' has to be positive ")
  if (toler.eigen <= 0) stop("'toler.eigen' has to be positive ")
  if (maxhalf < 0) stop("'maxhalf' has to be a non-negative integer ")

  return(list(est.c = est.c,
              est.scale = est.scale,
              maxiter = maxiter,
              firstiter = firstiter,
              rel.tolerance = rel.tolerance,
              toler.chol = toler.chol,
              toler.eigen = toler.eigen,
              maxhalf = maxhalf,
              debug = debug,
              info = info,
              lambda.use = lambda.use,
              sdspline = sdspline,
              difforder = difforder,
              knots = knots,
              nsplines = nsplines,
              last.three = last.three
         ))
}

