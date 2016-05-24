#' A Comparison of Regression Modelling Strategies
#'
#' Compare two nested modelling strategies and return measures
#' of their relative predictive performance
#'
#' This is the core function in the \pkg{apricom} package. The *compare* function can be used
#' to compare the performance of two prediction model building approaches for either simulated
#' or user-specified data. For further details, see the \pkg{apricom} user manual.
#'
#' The following strategies are currently supported: heuristic shrinkage ("heuristic"),
#' split-sample-derived shrinkage ("split"), cross-validation-derived shrinkage ("kcv"),
#' leave-one-out cross-validation-derived shrinkage ("loocv"), bootstrap-derived shrinkage ("boot")
#'  and penalized logistic regression using Firth's penalty ("pml.firth"). Furthermore, models built
#'  using these methods may be compared with raw models fitted by ordinary least squares
#'  estimation ("lsq") or maximum likelihood estimation ("ml").
#'
#'   Strategies should be specified within the "strat1" and "strat2" arguments in the form of a list,
#'   starting with the strategy name (as listed above in parentheses), followed by relevant
#'   parameters for each respective method. For further details see individual help files for each
#'   strategy, and the examples below. Note that in the *compare* function call, the dataset should
#'   not be specified within the "strat1" or "strat2" arguments, and instead should only be called
#'   within the "data" argument.
#'
#' @import utils
#' @import shrink
#'
#' @param model     the type of regression model. Either "linear" or "logistic".
#' @param Ntrials   number of simulation trials.
#' @param strat1    a list containing the strategy name and strategy-specific parameter values.
#'                  This modelling strategy is taken as the reference for comparison.
#' @param strat2    a list containing the strategy name and strategy-specific parameter values.
#'                  This modelling strategy is compared with the reference strategy, strat1.
#' @param data      a list describing the dataset in which the selected modelling strategies
#'                  will be compared. If the first object in the list is "norm" or "unif",
#'                  the user may submit parameters for generating multivariable simulated
#'                  datasets (see details below. Users may specify their own dataset using
#'                  the format \code{list("dataset", user_data)}, where the second object in
#'                  the list is their own dataset.
#' @param Nrows     the number of rows of observations in simulated datasets of type
#'                  "norm" or "unif".
#' @param Ncomp     the number of rows of observations in the comparison set. This dataset is
#'                  taken to represent the overall population, from which the training set is
#'                  sampled. When \code{data} is of type \code{"dataset"}, if Ncomp is not
#'                  specified, the original data will be used as the comparison dataset.
#' @param int       logical. If \code{int == TRUE} an intercept will be included in the
#'                  regression model.
#' @param int.adj   logical. If \code{int.adj == TRUE} the intercept will be re-estimated
#'                  after shrinkage is applied to the regression coefficients.
#' @param trim      logical. If \code{trim == TRUE} a "trimmed" comparison distribution
#'                  will be returned, along with a victory rate and median precision
#'                  ratio derived using the trimmed distribution. The trimmed distribution
#'                  only contains precision ratios within a range of plus or minus two times
#'                  the interquartile range around the median precision ratio.
#' @param output    logical. If \code{output == TRUE} the function will return two graphical
#'                  representations of the comparison distribution.
#' @return \code{compare} returns a list containing the following:
#' @return \item{VR}{the victory rate of strategy 2 over strategy 1.}
#' @return \item{MPR}{the median precision ratio over Ntrials comparison trials.}
#' @return \item{PR.IQR}{the precision ratio interquartile range over \code{Ntrials} comparison trials.}
#' @return \item{VR.trim}{if \code{trim == TRUE} the "trimmed" victory rate of
#'               strategy 2 over strategy 1 is returned.}
#' @return \item{MPR.trim}{if \code{trim == TRUE} the "trimmed" median precision ratio over
#'                        Ntrials comparison trials is returned.}
#' @return \item{distribution}{the comparison distribution of strategy 2 vs. strategy 1. This is
#'                            the distribution of precision ratios generated from \code{Ntrials} comparison trials}
#' @return \item{distribution.trim}{if \code{trim == TRUE} the "trimmed" comparison distribution is returned.}
#' @return \item{N.rejected}{the number of trials excluded from the comparison distribution by trimming}
#' @return \item{strat1}{modelling strategy 1}
#' @return \item{strat2}{modelling strategy 2}
#' @return \item{shrinkage1}{If strategy 1 is a shrinkage-after-estimation technique, a vector or matrix
#'               containing the shrinkage factor estimated in each trial is returned}
#' @return \item{shrinkage1}{If strategy 1 is a shrinkage-after-estimation technique, a vector or matrix
#'               containing the shrinkage factor estimated in each trial is returned}
#'
#' @note When using \code{compare} it is strongly recommended that ideally 10000
#'       comparison trials are used, to give stable estimates. Comparisons with
#'       logistic regression modelling model adjustment strategies are \emph{considerably}
#'       slower than with linear regression, and 1000-5000 trials may be preferred. The
#'       examples provided in this documentation use considerably fewer comparison trials
#'       and yield highly unstable estimates.
#'
#' @examples
#'## Example 1: Comparison of heuristic formula-derived shrinkage against
#'## a raw least squares model. Data is simulated multivariable random
#'## normally distributed.The comparison set will have 2000 rows. Here only
#'## 10 trial replicates are used, but at least 1000 should be used in practice
#'
#'  mu <- c(rep(0, 21))
#'  rho <- 0.5
#'  comp1 <- compare(model = "linear", Ntrials = 10, strat1 = list("lsq"),
#'           strat2 = list("heuristic", DF = 8),
#'           data = list("norm", mu, rho), Nrows = 200, Ncomp = 2000,
#'           int = TRUE, int.adj = FALSE, trim = FALSE, output = TRUE)
#'
#'
#'## Example 2: A truncated comparison of 10-rep, 10-fold
#'## cross-validation-derived shrinkage against leave-one-out cross-validation.
#'## Data is simulated multivariable random uniformly distributed
#'## (50 rows; 5 predictors with mean=0 ; r^2 = 0.7)
#'## The comparison set will contain 1000 observations.
#'
#' mu <- c(rep(0, 6))
#' rho <- 0.7
#' comp2 <- compare(model = "linear", Ntrials = 10, strat1 = list("loocv"),
#'           strat2 = list("kcv", k = 10, nreps = 10),data = list("unif", mu, rho),
#'           Nrows = 50, Ncomp = 1000, trim = TRUE)
#'
#'
#'## Example 3:  Comparison of penalized logistic regression with
#'## Firth's penalty against raw logistic regression model using
#'## maximum likelihood estimation.
#'## Note that the logistf package is required for pml.firth.
#'
#' library(shrink)
#' data(deepvein)
#' dv.data <- datashape(deepvein, y = 3, x = 4:11)
#' set.seed(123)
#' comp4 <- compare(model = "logistic", Ntrials = 10,
#'          strat1 = list("ml"), strat2 = list("pml.firth"),
#'          data = list("dataset", dv.data), int = TRUE,
#'          int.adj = TRUE, trim = FALSE, output = TRUE)
#'
#'@references Pestman W., Groenwold R. H. H., Teerenstra. S, \emph{"Comparison of
#'            strategies when building linear prediction models."}
#'            Numerical Linear Algebra with Applications (2013)

compare<- function(model, Ntrials, strat1, strat2, data, Nrows, Ncomp,
                   int = TRUE, int.adj, trim = FALSE, output = TRUE) {

  if (nargs() < 6) stop("model, Ntrials, strat1, strat2 and data must be specified")
  if (Ntrials < 1000) warning("Victory rates based on less than 1000
                              trials are unstable")

  if (model == "linear") {
    sse2 <- c(rep(NA, Ntrials))
    sse1 <- c(rep(NA, Ntrials))
  } else if (model == "logistic") {
    neg2LL2 <- c(rep(NA, Ntrials))
    neg2LL1 <- c(rep(NA, Ntrials))
  }

  if (strat1[1] == "heuristic") save.lambda1 <- TRUE
  else if (strat1[1] == "boot") save.lambda1 <- TRUE
  else if (strat1[1] == "kcv") save.lambda1 <- TRUE
  else if (strat1[1] == "loocv") save.lambda1 <- TRUE
  else if (strat1[1] == "split") save.lambda1 <- TRUE
  else save.lambda1 = FALSE
  if (strat2[1] == "heuristic") save.lambda2 <- TRUE
  else if (strat2[1] == "boot") save.lambda2 <- TRUE
  else if (strat2[1] == "kcv") save.lambda2 <- TRUE
  else if (strat2[1] == "loocv") save.lambda2 <- TRUE
  else if (strat2[1] == "split") save.lambda2 <- TRUE
  else save.lambda2 = FALSE

  strat1.lambda <- c(rep(0, Ntrials))
  strat2.lambda <- c(rep(0, Ntrials))
  total <- Ntrials

  pb <- txtProgressBar(min = 0, max = total, style = 3)

### SIMULATIONS WITH MVNORM DATASETS: LINEAR REGRESSION ONLY
  if (data[1] == "norm") {
    if (missing(Nrows) | missing(Ncomp)) stop("Number of rows in dataset (Nrows)
                                      and comparison set (Ncomp) must be specified")
    mu <- data[[2]]
    Cov <- data[[3]]
## If C is rho not a cov matrix
    if (length(Cov) == 1) {
      lmu<- length(mu) - 1
      C0 <- diag(lmu + 1)
      C0[lmu + 1, 1] <- sqrt(Cov)
      C0[1, lmu + 1] <- sqrt(Cov)
      Cov <- C0
    }

    for (k in 1 : Ntrials) {
      setTxtProgressBar(pb, k)
      d.train <- randnor(Nrows, mu, Cov)
      d.comp <- randnor(Ncomp, mu, Cov)

      if (int == TRUE) {
        d.train <- cbind(1, d.train)
        d.comp <- cbind(1, d.comp)
      }

      strat2.out <- strategy(d.train, strat2, int, int.adj, model)
      if (save.lambda2 == FALSE) strat2.coeff <- strat2.out
      else {
        strat2.coeff <- strat2.out$shrunk.coeff
        strat2.lambda[k] <- strat2.out$lambda
      }
      sse2[k] <- sse(strat2.coeff, d.comp)

      strat1.out <- strategy(d.train, strat1, int, int.adj, model)
      if (save.lambda1 == FALSE) strat1.coeff <- strat1.out
      else {
        strat1.coeff <- strat1.out$shrunk.coeff
        strat1.lambda[k] <- strat1.out$lambda
      }
      sse1[k] <- sse(strat1.coeff,d.comp)

      Sys.sleep(0.0000001)
    }
    close(pb)
  }

### SIMULATIONS WITH MVUNIF DATASETS: LINEAR REGRESSION ONLY
  else if (data[1] == "unif") {
    if (missing(Nrows) | missing(Ncomp)) stop("Number of rows in dataset (Nrows)
                                    and comparison set (Ncomp) must be specified")
    mu <- data[[2]]
    Cov <- data[[3]]
    if (length(Cov) == 1) {
      lmu<- length(mu) - 1
      C0 <- diag(lmu + 1)
      C0[lmu + 1, 1] <- sqrt(Cov)
      C0[1, lmu + 1] <- sqrt(Cov)
      Cov <- C0
    }
    Q <- diag(length(mu))

    for (k in 1 : Ntrials) {
      setTxtProgressBar(pb, k)
      d.train <- randunif(Nrows, mu, Cov, Q)
      d.comp <- randunif(Ncomp, mu, Cov, Q)
      if (int == TRUE) {
        d.train <- cbind(1, d.train)
        d.comp <- cbind(1, d.comp)
      }

      strat2.out <- strategy(d.train, strat2, int, int.adj, model)
      if (save.lambda2 == FALSE) strat2.coeff <- strat2.out
      else {
        strat2.coeff <- strat2.out$shrunk.coeff
        strat2.lambda[k] <- strat2.out$lambda
      }
      sse2[k] <- sse(strat2.coeff, d.comp)

      strat1.out <- strategy(d.train, strat1, int, int.adj, model)
      if (save.lambda1 == FALSE) strat1.coeff <- strat1.out
      else {
        strat1.coeff <- strat1.out$shrunk.coeff
        strat1.lambda[k] <- strat1.out$lambda
      }
      sse1[k] <- sse(strat1.coeff,d.comp)

      Sys.sleep(0.0000001)
    }
    close(pb)
  }

### COMPARISONS WITH REAL (USER SUPPLIED) DATASETS
  else if (data[1] == "dataset") {
    data0 <- as.matrix(data.frame(data[2]))
    Nrows <- dim(data0)[1]

    for (k in 1:Ntrials) {
      setTxtProgressBar(pb, k)
      d.train <- randboot(data0, Nrows)
      if (missing(Ncomp)) d.comp <- data0
      else d.comp <- randboot(data0, Ncomp)

      if (int == TRUE) {
        d.train <- cbind(1, d.train)
        d.comp <- cbind(1, d.comp)
      }

      if (model == "linear") {

        strat2.out <- strategy(d.train, strat2, int, int.adj, model)
        if (save.lambda2 == FALSE) strat2.coeff <- strat2.out
        else {
          strat2.coeff <- strat2.out$shrunk.coeff
          strat2.lambda[k] <- strat2.out$lambda
        }
        sse2[k] <- sse(strat2.coeff, d.comp)

        strat1.out <- strategy(d.train, strat1, int, int.adj, model)
        if (save.lambda1 == FALSE) strat1.coeff <- strat1.out
        else {
          strat1.coeff <- strat1.out$shrunk.coeff
          strat1.lambda[k] <- strat1.out$lambda
        }
        sse1[k] <- sse(strat1.coeff,d.comp)

        Sys.sleep(0.0000001)

        } else if (model == "logistic") {

          strat2.out <- strategy(d.train, strat2, int, int.adj, model)
          if (save.lambda2 == FALSE) strat2.coeff <- strat2.out
          else {
            strat2.coeff <- strat2.out$shrunk.coeff
            strat2.lambda[k] <- strat2.out$lambda
          }
          neg2LL2[k] <- loglikelihood(strat2.coeff, d.comp)

          strat1.out <- strategy(d.train, strat1, int, int.adj, model)
          if (save.lambda1 == FALSE) strat1.coeff <- strat1.out
          else {
            strat1.coeff <- strat1.out$shrunk.coeff
            strat1.lambda[k] <- strat1.out$lambda
          }
          neg2LL1[k] <- loglikelihood(strat1.coeff,d.comp)

          Sys.sleep(0.0000001)
        }
      else stop("Model must be of type linear or logistic")
    }
    close(pb)
  }

  else stop("Invalid data source. data must be of type norm, unif or dataset")

  if (model == "linear") {

    if (save.lambda1==TRUE & save.lambda2==TRUE) comparison <- compdist(sse1, sse2,
      model, output, strat1.lambda, strat2.lambda, trim=trim, strat1=strat1, strat2=strat2)
    else if(save.lambda1==TRUE & save.lambda2==FALSE) comparison <- compdist(sse1, sse2,
      model, output, lambda1 = strat1.lambda, trim=trim, strat1=strat1, strat2=strat2)
    else if(save.lambda1==FALSE & save.lambda2==TRUE) comparison <- compdist(sse1, sse2,
      model, output, lambda2 = strat2.lambda, trim=trim, strat1=strat1, strat2=strat2)
    else comparison <- compdist(sse1, sse2, model, output, trim=trim, strat1=strat1,
      strat2=strat2)

    } else {

      if (save.lambda1==TRUE & save.lambda2==TRUE) comparison <- compdist(neg2LL1,
        neg2LL2, model, output, strat1.lambda, strat2.lambda, trim=trim,
        strat1=strat1, strat2=strat2)
      else if(save.lambda1==TRUE & save.lambda2==FALSE) comparison <- compdist(neg2LL1,
        neg2LL2, model, output, lambda1 = strat1.lambda, trim=trim, strat1=strat1,
        strat2=strat2)
      else if(save.lambda1==FALSE & save.lambda2==TRUE) comparison <- compdist(neg2LL1,
        neg2LL2, model, output, lambda2 = strat2.lambda, trim=trim, strat1=strat1,
        strat2=strat2)
      else comparison <- compdist(neg2LL1, neg2LL2, model, output, trim=trim,
        strat1=strat1, strat2=strat2)
    }
}
