truePrev <-
function(x, n, SE = 1, SP = 1, prior = c(1, 1),
         nchains = 2, burnin = 10000, update = 10000,
         verbose = FALSE) {

  ## check x and n
  if (missing(x))  stop("'x' is missing")
  if (missing(n))  stop("'n' is missing")
  checkInput(x, "x", class = "integer", min = 0)
  checkInput(n, "n", class = "integer", minEq = 0)
  binom <- length(x) == 1
  if (!binom & sum(x) != n) stop("'x' does not sum to 'n'")
  if (!binom & length(n) == 1) n <- rep(n, length(x))
  if (any(x > n)) stop("'x' cannot be larger than 'n'")

  ## check SE and SP
  checkInput(SE, "SE", class = c("formula", "list", "numeric"))
  checkInput(SP, "SP", class = c("formula", "list", "numeric"))
  Se <- checkBinPrior(SE, "SE")
  Sp <- checkBinPrior(SP, "SP")

  ## check prior
  checkInput(prior, "prior", class = "numeric", length = 2, minEq = 0)

  ## check nchains, burnin & update
  checkInput(nchains, "nchains", class = "integer", min = 2)
  checkInput(burnin, "burnin", class = "integer", min = 1)
  checkInput(update, "update", class = "integer", min = 1)

  ## check options
  checkInput(verbose, "verbose", class = "logical")

  ## get output
  out <- truePrevBinom(x, n, Se, Sp, prior,
                       nchains, burnin, update, verbose)

  ## return output
  return(out)
}