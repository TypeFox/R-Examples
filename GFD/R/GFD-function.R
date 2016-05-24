#' Tests for General Factorial Designs
#' 
#' The GFD function calculates the Wald-type statistic (WTS), the ANOVA-type 
#' statistic (ATS) as well as a permutation version of the WTS for general 
#' factorial designs.
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side
#'   contains the response variable and the right hand side contains the factor
#'   variables of interest. An interaction term must be specified.
#' @param data A data.frame, list or environment containing the variables in 
#'   \code{formula}. The default option is \code{NULL}.
#' @param nperm The number of permutations used for calculating the permuted 
#'   Wald-type statistic. The default option is 10000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#'   
#' @details The package provides the Wald-type statistic, a permuted version
#'   thereof as well as the ANOVA-type statistic for general factorial designs,
#'   even with non-normal error terms and/or heteroscedastic variances. It is
#'   implemented for both crossed and hierarchically nested designs and allows
#'   for an arbitrary number of factor combinations as well as different sample
#'   sizes in the crossed design.
#'   The \code{GFD} function returns three p-values: One for the ATS based on an F-quantile and
#'   two for the WTS, one based on the \eqn{\chi^2}
#'    distribution and one based on the permutation procedure.
#'   Since the ATS is only an approximation and the WTS based on the \eqn{\chi^2}
#'   distribution is known 
#'   to be very liberal for small sample sizes, we recommend to use the WTPS in these situations.
#'
#'   
#' @return A \code{GFD} object containing the following components:
#' \item{Descriptive}{Some descriptive statistics of the data for all factor
#'   level combinations. Displayed are the number of individuals per factor
#'   level combination, the mean, variance and 100*(1-alpha)\% confidence
#'   intervals.}
#'  \item{WTS}{The value of the WTS along with degrees of freedom of the central chi-square distribution and p-value, as well as the p-value of the permutation procedure.}
#'  \item{ATS}{The value of the ATS, degrees of freedom of the central F distribution and the corresponding p-value.}
#' 
#' @examples
#' GFD(weightgain ~ source * type, data = HSAUR::weightgain)
#' 
#' data(startup)
#' model <- GFD(Costs ~ company, data = startup)
#' summary(model)
#' 
#' @references Friedrich, S., Konietschke, F., Pauly, M.(2015). GFD - An R-package
#' for the Analysis of General Factorial Designs. Submitted to Journal of Statistical Software.
#' 
#' Pauly, M., Brunner, E., Konietschke, F.(2015). Asymptotic Permutation Tests in General Factorial Designs. Journal of the Royal Statistical Society - Series B 77, 461-473.
#' 
#' @importFrom graphics axis legend par plot title
#' @importFrom stats ecdf formula model.frame pchisq pf qt terms var
#' @importFrom utils read.table
#' @importFrom methods hasArg
#' 
#' @export

GFD <- function(formula, data = NULL, nperm = 10000,
                alpha = 0.05){
  
  input_list <- list(formula = formula, data = data,
                     nperm = nperm, alpha = alpha)
  dat <- model.frame(formula, data)
  subject <- 1:nrow(dat)
  dat2 <- data.frame(dat, subject = subject)
  nf <- ncol(dat) - 1
  nadat <- names(dat)
  nadat2 <- nadat[-1]
  fl <- NA
  for (aa in 1:nf) {
    fl[aa] <- nlevels(as.factor(dat[ ,aa + 1]))
  }
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[ ,jj + 1]))
  }
  lev_names <- expand.grid(levels)
  
  if (nf == 1) {
    # one-way layout
    dat2 <- dat2[order(dat2[, 2]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    # contrast matrix
    hypo <- diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl)
    WTS_out <- matrix(NA, ncol = 3, nrow = 1)
    ATS_out <- matrix(NA, ncol = 4, nrow = 1)
    WTPS_out <- rep(NA, 1)
    rownames(WTS_out) <- fac_names
    rownames(ATS_out) <- fac_names
    names(WTPS_out) <- fac_names
    results <- Stat(data = response, n = n, hypo, nperm = nperm, alpha)
    WTS_out <- results$WTS
    ATS_out <- results$ATS
    WTPS_out <- results$WTPS
    mean_out <- results$Mean
    Var_out <- results$Cov
    CI <- results$CI
    colnames(CI) <- c("CIl", "CIu")
    descriptive <- cbind(lev_names, n, mean_out, Var_out, CI)
    colnames(descriptive) <- c(nadat2, "n", "Means", "Variances",
                               paste("Lower", 100 * (1 - alpha), "%", "CI"),
                               paste("Upper", 100 * (1 - alpha), "%", "CI"))
    WTS_output <- c(WTS_out, WTPS_out)
    names(WTS_output) <- cbind ("Test statistic", "df",
                                "p-value", "p-value WTPS")
    names(ATS_out) <- cbind("Test statistic", "df1", "df2", "p-value")
    output <- list()
    output$input <- input_list
    output$Descriptive <- descriptive
    output$WTS <- WTS_output
    output$ATS <- ATS_out
    output$plotting <- list(levels, fac_names, nf)
    names(output$plotting) <- c("levels", "fac_names", "nf")
    # end one-way layout ------------------------------------------------------
  } else {
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]), ]
    # sorting data according to factors
    dat2 <- dat2[do.call(order, dat2[, 2:(nf + 1)]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    fac_names_original <- fac_names
    perm_names <- t(attr(terms(formula), "factors")[-1, ])
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    if (length(fac_names) == nf) {
      # delete factorcombinations which don't exist
      n <- n[n != 0]
    }
    # mixture of nested and crossed designs is not possible
    if (length(fac_names) != nf && 2 %in% nr_hypo) {
      stop("A model involving both nested and crossed factors is
           not impemented!")
    }
    # only 3-way nested designs are possible
    if (length(fac_names) == nf && nf >= 4) {
      stop("Four- and higher way nested designs are
           not implemented!")
    }
    # no factor combinations with less than 2 observations
    if (0 %in% n || 1 %in% n) {
      stop("There is at least one factor-level combination
           with less than 2 observations!")
    }
    # correct labeling of factors in nested design
    if (length(fac_names) == nf) {
      if (nf == 2) {
        if (all(levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][1]]))
                == levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][2]])))) {
          stop("The levels of the nested factor must be
               named without repetitions!")
        }
      } else if (nf == 3) {
        if (all(levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][1]]))
                == levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][2]]))) ||
              all(levels(as.factor(dat2[, 4][dat2[, 3] == levels[[2]][1]]))
                  == levels(as.factor(dat2[, 4][dat2[, 3] == levels[[2]][fl[2] / fl[1] + 1]])))) {
          stop("The levels of the nested factor must be
                 named without repetitions!")
        }
      }
    }
    if (length(fac_names) == nf) {
      # nested
      TYPE <- "nested"
      hypo_matrices <- HN(fl)
      # create correct level combinations
      blev <- list()
      lev_names <- list()
      for (ii in 1:length(levels[[1]])) {
        blev[[ii]] <- levels(as.factor(dat[, 3][dat[, 2] == levels[[1]][ii]]))
        lev_names[[ii]] <- rep(levels[[1]][ii], length(blev[[ii]]))
      }
      if (nf == 2) {
        lev_names <- as.factor(unlist(lev_names))
        blev <- as.factor(unlist(blev))
        lev_names <- cbind.data.frame(lev_names, blev)
      } else {
        lev_names <- lapply(lev_names, rep,
                            length(levels[[3]]) / length(levels[[2]]))
        lev_names <- lapply(lev_names, sort)
        lev_names <- as.factor(unlist(lev_names))
        blev <- lapply(blev, rep, length(levels[[3]]) / length(levels[[2]]))
        blev <- lapply(blev, sort)
        blev <- as.factor(unlist(blev))
        lev_names <- cbind.data.frame(lev_names, blev, as.factor(levels[[3]]))
      }
    } else {
      # crossed
      TYPE <- "crossed"
      hypo_matrices <- HC(fl, perm_names, fac_names)[[1]]
      fac_names <- HC(fl, perm_names, fac_names)[[2]]
    }
    if (length(fac_names) != length(hypo_matrices)) {
      stop("Something is wrong: Perhaps a missing interaction term in formula?")
    }
    WTS_out <- matrix(NA, ncol = 3, nrow = length(hypo_matrices))
    ATS_out <- matrix(NA, ncol = 4, nrow = length(hypo_matrices))
    WTPS_out <- rep(NA, length(hypo_matrices))
    rownames(WTS_out) <- fac_names
    rownames(ATS_out) <- fac_names
    names(WTPS_out) <- fac_names
    colnames(ATS_out) <- c("Test statistic", "df1", "df2", "p-value")
    # calculate results
    for (i in 1:length(hypo_matrices)) {
      results <- Stat(data = response, n = n, hypo_matrices[[i]],
                      nperm = nperm, alpha)
      WTS_out[i, ] <- results$WTS
      ATS_out[i, ] <- results$ATS
      WTPS_out[i] <- results$WTPS
    }
    mean_out <- results$Mean
    Var_out <- results$Cov
    CI <- results$CI
    colnames(CI) <- c("CIl", "CIu")
    descriptive <- cbind(lev_names, n, mean_out, Var_out, CI)
    colnames(descriptive) <- c(nadat2, "n", "Means", "Variances",
                               paste("Lower", 100 * (1 - alpha),"%", "CI"),
                               paste("Upper", 100 * (1 - alpha),"%", "CI"))
    
    # calculate group means, variances and CIs ----------------------------
    mu <- list()
    sigma <- list()
    n_groups <- list()
    lower <- list()
    upper <- list()
    for (i in 1:nf) {
      mu[[i]] <- c(by(dat2[, 1], dat2[, i + 1], mean))
      sigma[[i]] <- c(by(dat2[, 1], dat2[, i + 1], var))
      n_groups[[i]] <- c(by(dat2[, 1], dat2[, i + 1], length))
      lower[[i]] <- mu[[i]] - sqrt(sigma[[i]] / n_groups[[i]]) *
        qt(1 - alpha / 2, df = n_groups[[i]])
      upper[[i]] <- mu[[i]] + sqrt(sigma[[i]] / n_groups[[i]]) *
        qt(1 - alpha / 2, df = n_groups[[i]])
    }
       
    # Output ------------------------------------------------------
    WTS_output <- cbind(WTS_out, WTPS_out)
    colnames(WTS_output) <- cbind ("Test statistic", "df", "p-value",
                                   "p-value WTPS")
    output <- list()
    output$input <- input_list
    output$Descriptive <- descriptive
    output$WTS <- WTS_output
    output$ATS <- ATS_out
    output$plotting <- list(levels, fac_names, nf, TYPE, mu, lower, upper, fac_names_original, dat2, fl, alpha, nadat2, lev_names)
    names(output$plotting) <- c("levels", "fac_names", "nf", "Type", "mu", "lower", "upper", "fac_names_original", "dat2", "fl", "alpha", "nadat2", "lev_names")
  }
  class(output) <- "GFD"
  return(output)
}
