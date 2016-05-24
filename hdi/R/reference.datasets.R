rXb <-
  function(n, p, s0,
           xtype = c("toeplitz", "exp.decay", "equi.corr"),
           btype = "U[-2,2]", permuted = FALSE,
           iteration = NA, do2S = TRUE,
           x.par = switch(xtype,
                          "toeplitz"  = 0.9,
                          "equi.corr" = 0.8,
                          "exp.decay" = c(0.4, 5)),
           verbose = TRUE)
{
  ## Purpose: generate the reference datasets used in the paper accompanying
  ##          this R package
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Oct 2015; 'do2S', 'x.par'
  ## and tweaks: Martin Maechler

  ## Checking arguments
  xtype <- match.arg(xtype)

  stopifnot(is.character(btype), length(btype) == 1,
            n == as.integer(n), length(n) == 1, n >= 1,
            p == as.integer(p), length(p) == 1, p >= 1,
            s0== as.integer(s0),length(s0)== 1, 0 <= s0, s0 <= p)
  do.seed <- is.numeric(iteration) && !is.na(iteration)

  if(do.seed) {
    if(verbose) {
      cat("A value for the argument iteration has been specified:\n")
      cat("The seed will be set for reproducibility, the old RNGstate will be restored after the data generation.\n")
    }

    ## Based on the example of
    ## stats:::simulate.lm
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))

    ## Set seeds in such a way that iteration 1:50
    ## correspond to the setups used in the paper
    seed <- iteration + 2
    set.seed(seed)
  }

  x <- rX(n = n, p = p, xtype = xtype,
          permuted = permuted, do2S = do2S, par = x.par)
  if(do.seed && s0 > 0 && !grepl("^bfix", btype)) # only set if using RNG
    set.seed(seed)
  beta <- rb(p = p, s0 = s0, btype = btype)
  list(x=x, beta=beta)
}

rX <-
  function(n, p, xtype, permuted, do2S = TRUE,
           par = switch(xtype,
                        "toeplitz"  = 0.9,
                        "equi.corr" = 0.8,
                        "exp.decay" = c(0.4, 5)))
{
  ## Purpose: generate the reference design matrix used in the paper accompanying
  ##          this R package
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Oct 2015; 'par': Martin Maechler, 17 Dec 2015

  xtype <- tolower(xtype)## in case capitalization was used
  stopifnot(is.numeric(par), is.finite(par))
  Sigma <- switch(xtype,
                  "toeplitz" = {
                    indices <- toeplitz(0:(p-1))
                    stopifnot(length(par) == 1)
                    cov <- par ^ abs(indices)
                    if(do2S) ## This looks really stupid but the tiny numerical differences
                      ## make it identical to the datasets used in the paper
                      solve(solve(cov))
                    else cov
                  },
                  "equi.corr" = {
                    stopifnot(length(par) == 1)
                    cov <- matrix(par, p,p)
                    diag(cov) <- 1
                    if(do2S) solve(solve(cov)) else cov
                  },
                  "exp.decay" = {
                    stopifnot(length(par) == 2)
                    indices <- toeplitz(0:(p-1))
                    solve(par[1]^ (abs(indices)/par[2]))
                  },

                  stop("Invalid 'xtype': Must be one of 'toeplitz', 'equi.corr' or 'exp.decay'"))

  x <- mvrnorm(n, rep(0,p), Sigma)

  if(permuted)
    x[, sample.int(ncol(x))]
  else
    x
}

rb <- function(p, s0, btype) {
  ## Purpose: generate the reference coefficient vector used in the paper
  ##          accompanying this R package
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Ruben Dezeure, Date: 6 Oct 2015;  tweaks by Martin Maechler

  stopifnot(s0 <= p, p >= 0, length(s0) == 1, length(p) == 1,
            is.character(btype), length(btype) == 1)
  invalid.btype.txt <- "Invalid btype: Please provide a btype of the form 'bfix*' for a fixed value or 'U[*,*]', where * are two numbers, the lower and upper bounds."

  if(grepl("^U\\[",btype) &&
     grepl("\\]$", btype) &&
     sapply(regmatches(btype, gregexpr(",", btype)), length) == 1) {
    ## there should be one and only one comma in btype
    ## generate random coefficients from a uniform distribution
    split.btype <- strsplit(sub("^U\\[", '', sub("\\]$", '', btype)),
                            ",")[[1]]

    ## extract the lower and upper bounds:
    lower <- as.numeric(split.btype[1])
    upper <- as.numeric(split.btype[2])
    if(is.na(lower) || is.na(upper) || length(lower) != 1 || length(upper) != 1)
      ## if there were any errors in the extraction
      stop(invalid.btype.txt)

    b <- runif(s0, lower, upper)

  } else if(grepl("^bfix", btype)) {

    b <- as.numeric(sub("^bfix", "", btype)) ## extract the bfix number
    if(is.na(b)) ## if there were any errors in the extraction
      stop(invalid.btype.txt)

    b <- rep(b, s0)
  }
  else
    stop(invalid.btype.txt)

  c(b, rep(0, p - s0))
}
