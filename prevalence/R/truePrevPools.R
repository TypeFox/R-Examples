truePrevPools <-
function(x, n, SE = 1, SP = 1, prior = c(1, 1), 
         nchains = 2, burnin = 10000, update = 10000,
         verbose = FALSE) {

  ## check x and n
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  checkInput(x, "x", class = "integer", value = c(0, 1))
  checkInput(n, "n", class = "integer", minEq = 0)
  if (length(x) > 1 & length(n) == 1)  n <- rep(n, length(x))
  if (length(x) != length(n))  stop("'x' and 'n' must be of same length")
  if (length(x) == 1)  stop("\"truePrevPools\" requires at least 2 pools")

  ## check SE & SP
  checkInput(SE, "SE", class = c("formula", "list", "numeric"))
  checkInput(SP, "SP", class = c("formula", "list", "numeric"))
  Se <- checkBinPrior(SE)
  Sp <- checkBinPrior(SP)

  ## check prior
  checkInput(prior, "prior", class = "numeric", length = 2, minEq = 0)

  ## check nchains, burnin & update
  checkInput(nchains, "nchains", class = "integer", min = 2)
  checkInput(burnin, "burnin", class = "integer", min = 1)
  checkInput(update, "update", class = "integer", min = 1)

  ## check options
  checkInput(verbose, "verbose", class = "logical")

  ## create model
  model <- character()
  model[1] <- "model {"
  model[2] <- "for (i in 1:N) {"
  model[3] <- "x[i] ~ dbern(AP[i])"
  model[4] <- paste("AP[i] <- SEpool[i] * (1 - pow(1 - TP, n[i])) +",
                    "(1 - SPpool[i]) * pow(1 - TP, n[i])")
  model[5] <- paste("SEpool[i] <- 1 - (pow(1 - SE, n[i] * TP) *",
                    "pow(SP, n[i] * (1 - TP)))")
  model[6] <- "SPpool[i] <- pow(SP, n[i])"
  model[7] <- "}"

  model <- c(model, writeSeSp("SE", Se))
  model <- c(model, writeSeSp("SP", Sp))

  model <- c(model, paste0("TP ~ dbeta(", prior[1], ", ", prior[2], ")"))

  model <- c(model, "}")

  class(model) <- "prevModel"

  ## create data
  data <- list(x = x, n = n, N = length(n))

  ## create inits
  inits <- NULL

  ## get results!
  if (verbose) cat("JAGS progress:\n\n")

  JAGSout <- R2JAGS(model = model, data = data, inits = inits,
                    nchains = nchains, burnin = burnin, update = update,
                    nodes = c("SE", "SP", "TP"), verbose = verbose)

  mcmc.list <- JAGSout$mcmc.list
  class(mcmc.list) <- c("list", "mcmc.list")

  nodes <- colnames(mcmc.list[[1]])                      # extract node names
  mcmc.list_list <- list()                               # initiate list
  for (i in seq_along(nodes))                            # assign nodes
    mcmc.list_list[[i]] <- mcmc.list[, i]
  names(mcmc.list_list) <- nodes                         # assign node names
  mcmc.list_list <- mcmc.list_list[c("TP", "SE", "SP")]  # reorder elements

  ## define diagnostics
  # deviance information criterion
  DIC <- JAGSout$dic

  # brooks-gelman-rubin diagnostic
  # exclude fixed nodes
  exclude <- which(apply(mcmc.list[[1]], 2, sd) == 0)
  if (length(exclude) > 0) {
    BGR <- gelman.diag(mcmc.list[, -exclude])
  } else {
    BGR <- gelman.diag(mcmc.list)
  }

  ## get output
  out <- new("prev",
             par = list(x = x, n = n, SE = Se, SP = Sp, prior = prior,
                        nchains = nchains, burnin = burnin, update = update,
                        inits = inits),
             model = model,
             mcmc = mcmc.list_list,
             diagnostics = list(DIC = DIC,
                                BGR = BGR))

  ## return output
  return(out)
}
