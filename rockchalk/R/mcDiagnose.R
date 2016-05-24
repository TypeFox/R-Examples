##' Estimate leave-one-variable-out regressions
##'
##' This is a convenience for analysis of multicollinearity in regression.
##'
##' @param model a fitted regression model
##' @return a list of fitted regressions, one for each omitted variable.
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @export
lmAuxiliary <- function(model){
  dat <- as.data.frame(model.matrix(model))
  ## ivnames <- attr(delete.response(terms(model)), "term.labels")
  ## previous does not work with transforms like poly
  hasIntercept <- attr(terms(model), "intercept")
  if (hasIntercept) dat <- dat[ , -1] # removes intercept. assumes intercept in column 1
  ivnames <- colnames(dat)
  cat("The following auxiliary models are being estimated and returned in a list:\n")

  results <- list()
  colnames(dat) <- gsub("`","",  colnames(dat)) #rm ticks on varnames
  for (i in ivnames) {
    fmla <- paste( "`",i,"`", " ~ ." , sep = "") #add ticks in fmla
    fmla <- gsub("``","`",  fmla) # remove double ticks
    lmcall <- call("lm", formula(fmla, data = dat), quote(dat))
    results[[ i ]] <- maux <- eval(lmcall)
    print(formula(maux))
  }
 results
}
NULL

##' retrieves estimates of the coefficient of determination from a list of regressions
##'
##' Asks each regression model in a list for a summary and then reports the R-squares.
##' @param auxRegs a list of fitted regression objects
##' @return a numeric vector of the same length as auxRegs.
##' @author Paul E. Johnson <pauljohn@@ku.edu>
getAuxRsq <- function(auxRegs){
  auxRsq <- numeric(length=length(auxRegs))
  j <- 0
  for ( i in auxRegs ){
    j <- j + 1
    sumry <- summary(i)
    auxRsq[j] <- sumry$r.squared
  }

  names(auxRsq) <- names(auxRegs)

  invisible(auxRsq)
}
NULL


##' Converts the R-square to the variance inflation factor
##'
##' calculates vif = 1/(1-R-square)
##'
##' @param rsq a vector of real values, presumably fitted R-squares
##' @return a vector of vif estimates
##' @author Paul E. Johnson <pauljohn@@ku.edu>
getVIF <- function(rsq){
  vif <- 1/(1-rsq)

  invisible(vif)
}

##' Calculates the delta R-squares, also known as squared
##' semi-partial correlation coefficients.
##'
##' The change in the R-square when a variable is removed from a
##' regression is called delta R-square. It is sometimes suggested as
##' a way to determine whether a variable has a substantial effect on
##' an outcome. This is also known as the squared semi-partial correlation
##' coefficient.
##'
##' @export
##' @importFrom stats drop1
##' @param model a fitted regression model
##' @return a vector of estimates of the delta R-squares
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @examples
##' dat1 <- genCorrelatedData(N=250, means=c(100,100),
##' sds=c(30,20), rho=0.0,  stde = 7, beta=c(1.1, 2.4, 4.1, 0))
##' m1 <- lm(y ~ x1 + x2, data=dat1)
##' getDeltaRsquare(m1)
##' ## more problematic in presence of collinearity
##' dat2 <- genCorrelatedData(N=250, means=c(100,100),
##' sds=c(30,20), rho=0.6,  stde = 7, beta=c(1.1, 2.4, 4.1, 0))
##' m2 <- lm(y ~ x1 + x2, data=dat2)
##' getDeltaRsquare(m2)
getDeltaRsquare <- function(model){
  modeldrop1 <- drop1(model)
  RSS <- modeldrop1[1, "RSS"] ##Residual Sum of Squares
  deltaSS <- modeldrop1[ , "Sum of Sq"]
  SST = sum((model$model[ , 1] - mean(model$model[, 1]))^2)
  deltaRsquare <- deltaSS/SST
  names(deltaRsquare) <- row.names(modeldrop1)
  omit <- is.na(deltaRsquare)
  deltaRsquare <- deltaRsquare[-omit]
  cat("The deltaR-square values: the change in the R-square
      observed when a single term is removed.", fill = TRUE)
  cat("Same as the square of the 'semi-partial correlation coefficient'", fill = TRUE)
  print(as.data.frame(deltaRsquare))
  invisible(deltaRsquare)
}
NULL

##' Multi-collinearity diagnostics
##'
##' Conducts a series of checks for multicollinearity.
##' @param model a fitted regression model
##' @return a list of the "auxiliary regressions" that were
##' fitted during the analysis
##' @export
##' @importFrom stats model.matrix
##' @author Paul E. Johnson <pauljohn@@ku.edu>
##' @example inst/examples/mcDiagnose-ex.R
mcDiagnose <- function(model){
  if (any(is.na(coef(model)))) stop("There are redundant variables in the model. Fix the specification before diagnosing multicollinearity")
  auxRegs <- lmAuxiliary(model)
  auxRsq <- getAuxRsq(auxRegs)
  vif <- getVIF(auxRsq)

  cat("Drum roll please! \n")
  cat("\n")
  cat("And your R_j Squareds are (auxiliary Rsq)\n")
  print(auxRsq)
  cat("The Corresponding VIF, 1/(1-R_j^2)\n")
  print(vif)
  cat("Bivariate Correlations for design matrix \n")
  mm <- model.matrix(model)[ , -1] ## data, omit intercept
  print(round(cor(mm[,]), 2))
  invisible(auxRegs)
}
NULL

#TODO 1. Find out if there is a more numerically stable/recommended method
#TODO 2. figure smarter way to get rid of intercept variable than -1
##' Calculates partial correlation coefficients after retrieving data matrix froma fitted regression model
##'
##' The input is a fitted regression model, from which the design
##' matrix is retrieved, along with the dependent variable. The
##' partial correlation is calculated using matrix algebra that
##' has not been closely inspected for numerical precision.  That is to
##' say, it is in the stats book style, rather than the numerically
##' optimized calculating format that functions like \code{lm()} have adopted.
##'
##' I often criticize partial correlations because they change in a
##' very unstable way as terms are added or removed in regression
##' models. Nevertheless, I teach with books that endorse them, and in
##' order to have something to criticize, I need to have a function
##' like this. There are other packages that offer partial correlation
##' calculations, but they are either 1) not easy to install from CRAN
##' because of dependencies or 2) do not directly calculate the values
##' we want to see.
##'
##' To students. 1) This gives the same result as the function
##' \code{cov2pcor} in \code{gRbase}, so far as I can tell. Why use
##' this?  Simply for convenenience.  We have found that installing
##' gRbase is a problem because it depends on packages in
##' Bioconductor.  2) By default, I show only one column of output,
##' the partial correlations involving the dependent variable as
##' something being explained. The other columns that would depict the
##' dependent variable as a predictor of the independent variables
##' have been omitted. You can let me know if you think that's wrong.
##'
##' Please note I have not gone out of my way to make this calculation
##' "numerically stable."  It does not use any orthogonal matrix
##' calculations; it is using the same textbook theoretical stats
##' formula that is used by cov2pcor in gRbase and in every other
##' package or online source I could find. I prepared a little
##' WorkingExample file matrix-partial-correlations-1.R that
##' discusses this, in case you are interested (http://pj.freefaculty.org/R).
##'
##' @param model A fitted regression model, such as output from
##' lm(). Any object that has methods model.matrix and model.frame
##' will be sufficient.
##' @param dvonly Default = TRUE. Only show first column of the full
##' partial correlation matrix. That corresponds to the partial correlation of each predictor with y. I mean, r[yx].[others]
##' @return A column or matrix of partial correlation coefficients
##' @export
##' @author Paul E. Johnson <pauljohn@@ku.edu>
getPartialCor <- function(model, dvonly = TRUE){
    x <- model.matrix(model)
    modelterms <- terms(model)
    hasIntercept <- attr(modelterms, "intercept")
    if (hasIntercept) x <- x[ , -1] # removes intercept. assumes intercept in column 1
    yname <- as.character(modelterms[[2]])
    y <- model.response(model.frame(model))
    x <- cbind(y, x)
    colnames(x)[1] <- yname
    R <- cor(x)
    Rinv <- solve(R)
    Id <- diag(1/sqrt(diag(Rinv)))
    p <- - Id %*% Rinv %*% Id
    dimnames(p) <- list(colnames(x), colnames(x))
    if (dvonly) p[ ,1, drop = FALSE] else p
}

## There's got to be a better way.

## Center X, then proceed via QR decomposition?
## S-1 * R'R * S-1

## Solving that:
## S *  inv(R) inv(R') S




