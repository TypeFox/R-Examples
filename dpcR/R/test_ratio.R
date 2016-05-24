#' Rate ratio test
#' 
#' The test for comparing ratio of two Poisson means: \eqn{r =
#' \frac{\lambda_1}{\lambda_2}}{r = \lambda_1/\lambda_2}.
#' 
#' Objects \code{dpcr1} and \code{dpcr2} can be: \enumerate{ \item numeric
#' vectors of length 2. The first element is assumed to be number of positive
#' partitions and the second one to be the total number of partitions. \item
#' numeric vectors of length greater than 2. The length of vector is assumed to
#' represent total number of partitions. Every element of the vector with value
#' bigger than 0 is assumed to be a positive partitions.  \item
#' \code{\linkS4class{adpcr}} objects with type \code{tnp} (total number of
#' positive wells in panel) or \code{nm} (number of molecules per partition).
#' \code{\linkS4class{ddpcr}} objects with type \code{tnp} (total number of
#' positive droplets) or \code{nm} (number of molecules per droplet).  } Both
#' \code{dpcr1} and \code{dpcr2} must have the same class. See Examples.
#' 
#' The \code{ratio_test} is a wrapper around
#' \link[rateratio.test]{rateratio.test} function with custom input and output
#' tailored specifically for digital PCR experiments.
#' 
#' @aliases test_ratio test_ratio,numeric test_ratio,adpcr test_ratio,ddpcr
#' test_ratio,numeric-method test_ratio,adpcr-method test_ratio,ddpcr-method
#' test_ratio,numeric,numeric-method test_ratio,adpcr,adpcr-method
#' test_ratio,ddpcr,ddpcr-method
#' @param dpcr1 a (non-empty) numeric vector of data values of length 2 or more
#' or an object of class \code{\linkS4class{adpcr}} or
#' \code{\linkS4class{ddpcr}}. See Details.
#' @param dpcr2 a (non-empty) numeric vector of data values of length 2 or more
#' or an object of class \code{\linkS4class{adpcr}} or
#' \code{\linkS4class{ddpcr}}. See Details.
#' @param alternative alternative hypothesis, must be one of: \code{two.sided},
#' \code{greater} or \code{less}.
#' @param conf.level confidence level for the returned confidence interval.
#' @return An object of class `htest' containing the following components:
#' \item{p.value}{the p-value of the test} \item{estimate}{a vector with the
#' means (lambdas) of both experiments and their ratio}
#' \item{conf.int}{confidence interval for ratio between two experiments.}
#' \item{alternative}{type of alternative hypothesis} \item{method}{description
#' of method} \item{data.name}{description of data}
#' @author Michael Fay, Michal Burdukiewicz, Stefan Roediger
#' @seealso See also \code{\link[stats]{poisson.test}}.
#' @references Fay M.P. \emph{Two-sided exact tests and matching confidence
#' intervals for discrete data} R Journal 2 (1), 2010.
#' @keywords compare mean poisson
#' @examples
#' 
#' # Input values are numeric vectors representing dPCR experiments
#' x1 <- rpois(765, 1.1)
#' x2 <- rpois(765, 1.1)
#' test_ratio(x1, x2)
#' 
#' # Input values represent only number of positive partitions and total 
#' # partitions
#' x3 <- sum(rpois(765, 1.1) > 0)
#' x4 <- sum(rpois(765, 1.1) > 0)
#' test_ratio(c(x3, 765), c(x4, 765))
#' 
#' # It is possible to mix different types of input as long as they have 
#' # the same class
#' test_ratio(c(x3, 765), x1)
#' 
#' # The same is true for adpcr and dpcr objects.
#' x5 <- sim_adpcr(400, 1600, 100, pos_sums = TRUE, n_panels = 1)
#' x6 <- sim_adpcr(400, 1600, 100, pos_sums = FALSE, n_panels = 1)
#' test_ratio(x5, x6)
#' 
#' x7 <- sim_ddpcr(400, 1600, 100, pos_sums = TRUE, n_exp = 1)
#' x8 <- sim_ddpcr(400, 1600, 100, pos_sums = FALSE, n_exp = 1)
#' test_ratio(x7, x8)
#' 
#' @export test_ratio
test_ratio <- function(dpcr1, dpcr2, 
                       alternative = c("two.sided", "less", "greater"), 
                       conf.level = 0.95) {
  if (length(dpcr1) != 2) {
    n_x <- length(dpcr1)
    k_x <- sum(dpcr1 > 0)
  } else {
    n_x <- dpcr1[2]
    k_x <- dpcr1[1]
  }
  if (length(dpcr2) != 2) {
    n_y <- length(dpcr2)
    k_y <- sum(dpcr2 > 0)
  } else {
    n_y <- dpcr2[2]
    k_y <- dpcr2[1]
  }
  
  test_res <- rateratio.test(c(k_x, k_y), c(n_x, n_y), RR = 1, alternative = alternative, 
                             conf.level = conf.level)
  test_res[["data.name"]] <- paste0("dPCR 1: positive partitions: ", k_x, "; total partitions: ",
                                    n_x, ".\n       dPCR 2: positive partitions: ", k_y, 
                                    "; total partitions: ", n_y, ".")
  test_res[["estimate"]] <- c(calc_lambda(k_x, n_x)[1,2], calc_lambda(k_y, n_y)[1,2],
                              test_res[["estimate"]][1])
  test_res[["null.value"]] <- NULL
  names(test_res[["estimate"]]) <- c("Lambda1", "Lambda2", "Lambda1/Lambda2")
  test_res
}

setMethod("test_ratio", 
          signature(dpcr1 = "adpcr", dpcr2 = "adpcr"), 
          function(dpcr1, dpcr2, alternative = c("two.sided", "less", "greater"), 
                   conf.level = 0.95) {
            if(ncol(dpcr1) != 1 || ncol(dpcr2) != 1)
              stop("Both 'dpcr1' and 'dpcr2' must contain only one experiment.", 
                   call. = TRUE)
            if(!all(c(slot(dpcr1, "type"), slot(dpcr2, "type")) %in% c("nm", "tnp")))
              stop("Both 'dpcr1' and 'dpcr2' must have type 'nm' or 'tnp'", 
                   call. = TRUE)
            n_x <- slot(dpcr1, "n")
            n_y <- slot(dpcr2, "n")
            k_x <- ifelse(slot(dpcr1, "type") == "nm", sum(dpcr1 > 0), dpcr1[1])
            k_y <- ifelse(slot(dpcr2, "type") == "nm", sum(dpcr2 > 0), dpcr2[1])
            test_ratio(c(k_x, n_x), c(k_y, n_y), alternative = alternative, 
                       conf.level = conf.level)  
          })

setMethod("test_ratio", 
          signature(dpcr1 = "ddpcr", dpcr2 = "ddpcr"), 
          function(dpcr1, dpcr2, alternative = c("two.sided", "less", "greater"), 
                   conf.level = 0.95) {
            if(ncol(dpcr1) != 1 || ncol(dpcr2) != 1)
              stop("Both 'dpcr1' and 'dpcr2' must contain only one experiment.", 
                   call. = TRUE)
            if(!all(c(slot(dpcr1, "type"), slot(dpcr2, "type")) %in% c("nm", "tnp")))
              stop("Both 'dpcr1' and 'dpcr2' must have type 'nm' or 'tnp'", 
                   call. = TRUE)
            n_x <- slot(dpcr1, "n")
            n_y <- slot(dpcr2, "n")
            k_x <- ifelse(slot(dpcr1, "type") == "nm", sum(dpcr1 > 0), dpcr1[1])
            k_y <- ifelse(slot(dpcr2, "type") == "nm", sum(dpcr2 > 0), dpcr2[1])
            test_ratio(c(k_x, n_x), c(k_y, n_y), alternative = alternative, 
                       conf.level = conf.level)  
          })





#a wrapper around rateratio.test
# test_experiments <- function(dpcr_experiments) {
#   m_dat <- melt(dpcr_experiments)
#   colnames(m_dat)[1L:2] <- c("partition", "experiment")
#   m_dat[["experiment"]] <- factor(m_dat[["experiment"]])
#   
#   glm_fit <- glm(value ~ experiment, data = m_dat, family = quasipoisson)
#   
#   tuk <- glht(glm_fit, linfct = mcp(experiment = "Tukey"))
#   cld(tuk)
# }
