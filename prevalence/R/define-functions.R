###=========================================================================#
### 'DEFINE' FUNCTIONS
###=========================================================================#

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- define_x ........................ define observed test results
###-- define_prior .................... define prior for cond prob scheme
###-- define_prior2 ................... define prior for covariance scheme

## -------------------------------------------------------------------------#
## Define observed test results --------------------------------------------#

define_x <-
function(h) {
  ## check if h is defined
  if (missing(h))
    stop("The number of tests 'h' is not defined")

  ## check if h is defined correctly
  if (is.null(h) || is.na(h) || is.infinite(h))
    stop("The number of tests 'h' is not defined")
  if (is.character(h))
    stop("The number of tests 'h' must be a numeric value")
  checkInput(h, "h", class = "integer", min = 1)

  ## print title
  test <- ifelse(h == 1, "test:", "tests:")
  cat("Definition of the apparent test results, 'x', for", h, test)

  ## define test status for all APs
  status <- matrix(nrow = 2 ^ h, ncol = h)
  for (i in seq(h)) {
    status[, i] <- rep(c("+", "-"), each = 2 ^ (h - i), times = 2 ^ (i - 1))
  }

  ## paste output
  T <- character(2 ^ h)
  for (i in seq(2 ^ h))
    T[i] <-
      paste0("T", seq(h), status[i, ], ",", collapse = "")

  out <-
    paste0("\nx[", seq(2^h), "] : ",
           sapply(T, substr, start = 1, stop = nchar(T) - 1))

  cat(out, "\n")
}


## -------------------------------------------------------------------------#
## Define prior for conditional probability scheme -------------------------#

define_prior <-
function(h) {
  ## check if h is defined
  if (missing(h))
    stop("The number of tests 'h' is not defined")

  ## check if h is defined correctly
  if (is.null(h) || is.na(h) || is.infinite(h))
    stop("The number of tests 'h' is not defined")
  if (is.character(h))
    stop("The number of tests 'h' must be a numeric value")
  checkInput(h, "h", class = "integer", min = 1)

  ## print title
  test <- ifelse(h == 1, "test:", "tests:")
  cat("Conditional probability scheme\n")
  cat("Definition of the prior, 'theta', for", h, test, "\n")

  ## print theta[1-3]
  cat("theta[1] : P(D+) = TP\n")
  cat("theta[2] : P(T1+|D+) = SE1\n")
  cat("theta[3] : P(T1-|D-) = SP1\n")

  ## print remaining thetas, if needed
  if (h > 1) {
    t <- 4
    for (i in 2:h) {
      N <- 2 ^ i  # number of thetas for test i
      for (k in seq(N)) {
        D <- ifelse(k <= (N/2), "+", "-")  # true disease status
        T <- ifelse(k <= (N/2), "+", "-")  # current test status
        out <- paste0("P(T", i, T, "|D", D)

        for (j in seq(i-1)) {
          select <- rep(c("+", "-"), each = (2^(i-j-1)), times = (2^(j-1)))
          select <- c(select, rev(select))
          out <- paste0(out, ",T", j, select[k])
        }

        cat(paste("theta[", t, "] : ", out, ")\n", sep = ""))
        t <- t + 1
      }
    }
  }
}


## -------------------------------------------------------------------------#
## Define prior for covariance scheme --------------------------------------#

define_prior2 <-
function(h) {
  ## check if h is defined
  if (missing(h))
    stop("The number of tests 'h' is not defined")

  ## check if h is defined correctly
  if (is.null(h) || is.na(h) || is.infinite(h))
    stop("The number of tests 'h' is not defined")
  if (is.character(h))
    stop("The number of tests 'h' must be a numeric value")
  checkInput(h, "h", class = "integer", min = 1)

  ## print title
  test <- ifelse(h == 1, "test:", "tests:")
  cat("Covariance scheme\n")
  cat("Definition of the prior for", h, test, "\n")

  ## define node labels
  nodes <- get_nodes(h)

  ## define node names
  n_test <- rev(seq(h, 2))
  n_comb <- choose(h, n_test)

  cov <- character()
  for (i in seq_along(n_test))
    cov <- c(cov, apply(combn(h, n_test[i]), 2, paste, collapse = ",T"))

  names <-
    c("   True Prevalence",
      paste0("Sensitity T", seq(h)),
      paste0("Specificity T", seq(h)),
      paste0(" Covariance(", paste0("T", cov), "|D+)"),
      paste0(" Covariance(", paste0("T", cov), "|D-)"))

  ## print labels and names
  for (i in seq_along(nodes)) {
    cat(nodes[i], ":", names[i], "\n")
  }

}