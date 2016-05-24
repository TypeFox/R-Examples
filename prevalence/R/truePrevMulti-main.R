###=========================================================================#
### TRUE PREVALENCE FROM MULTIPLE TESTS / main functions
###=========================================================================#

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- truePrevMulti ................... user interface for cond prob scheme
###-- truePrevMulti2 .................. user interface for covariance scheme
###-- truePrevMultinom_conditional .... create model for cond prob scheme
###-- truePrevMultinom_covariance ..... create model for covariance scheme


## -------------------------------------------------------------------------#
## User interface for conditional probability scheme -----------------------#

truePrevMulti <-
function(x, n, prior, nchains = 2, burnin = 10000, update = 10000,
         verbose = FALSE) {

  ## check x and n
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  checkInput(x, "x", class = "integer", min = 0)
  checkInput(n, "n", class = "integer", minEq = 0)
  if (sum(x) != n) stop("'x' does not sum to 'n'")
  if ((log(length(x), 2) %% 1 != 0) | length(x) < 4) {
    stop("'x' is not correctly specified; see ?define_x")
  }

  ## check prior
  if (missing(prior)) stop("'prior' is missing")
  prior <- checkMultiPrior_conditional(substitute(prior))

  ## check nchains, burnin & update
  checkInput(nchains, "nchains", class = "integer", min = 2)
  checkInput(burnin, "burnin", class = "integer", min = 1)
  checkInput(update, "update", class = "integer", min = 1)

  ## check options
  checkInput(verbose, "verbose", class = "logical")

  ## get output
  out <-
    truePrevMultinom_conditional(x, n, prior, nchains, burnin, update, verbose)

  ## return output
  return(out)
}


## -------------------------------------------------------------------------#
## User interface for covariance scheme ------------------------------------#

truePrevMulti2 <-
function(x, n, prior, nchains = 2, burnin = 10000, update = 10000,
         verbose = FALSE) {

  ## check x and n
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  checkInput(x, "x", class = "integer", min = 0)
  checkInput(n, "n", class = "integer", minEq = 0)
  if (sum(x) != n) stop("'x' does not sum to 'n'")
  if ((log(length(x), 2) %% 1 != 0) | length(x) < 4) {
    stop("'x' is not correctly specified; see ?define_x")
  }

  ## check prior
  if (missing(prior)) stop("'prior' is missing")
  prior <-
    checkMultiPrior_covariance(substitute(prior), log2(length(x)))

  ## check nchains, burnin & update
  checkInput(nchains, "nchains", class = "integer", min = 2)
  checkInput(burnin, "burnin", class = "integer", min = 1)
  checkInput(update, "update", class = "integer", min = 1)

  ## check options
  checkInput(verbose, "verbose", class = "logical")

  ## get output
  out <-
    truePrevMultinom_covariance(x, n, prior, nchains, burnin, update, verbose)

  ## return output
  return(out)
}


## -------------------------------------------------------------------------#
## Create model for conditional probability scheme -------------------------#

truePrevMultinom_conditional <-
function(x, n, prior, nchains, burnin, update, verbose) {

  ## create model
  h <- log(length(x), 2)     # number of tests
  ntheta <- 2 ^ (h + 1) - 1  # number of thetas

  model <- character()

  ## write model initiation
  model[1] <- "model {"
  model[2] <- paste0("x[1:", 2 ^ h,
                     "] ~ dmulti(AP[1:", 2 ^ h, "], n)")

  ## write AP[] definitions in terms of theta[]
  s <- multiModel_select(h)  # define theta construct for SE/SP
  p <- multiModel_probs(s)   # define AP[.] in terms of theta[.]
  model <- c(model, "", p, "")

  ## write theta[] prior
  for (i in seq(ntheta))
    model <- c(model,
      writeSeSp(paste0("theta[", i, "]"), prior[[i]]))

  ## write bayesP definition
  bayesP <-
    c(paste0("x2[1:", (2^h), "] ~ dmulti(AP[1:", (2^h), "], n)"),
      paste0("for (i in 1:", (2^h), ") {"),
      "d1[i] <- x[i] * log(max(x[i],1) / (AP[i]*n))",
      "d2[i] <- x2[i] * log(max(x2[i],1) / (AP[i]*n))",
      "}",
      "G0 <- 2 * sum(d1[])",
      "Gt <- 2 * sum(d2[])",
      "bayesP <- step(G0 - Gt)")
  model <- c(model, "", bayesP)

  ## write SE[]/SP[] definition
  model <- c(model, "", multiModel_SeSp(h))

  ## close model
  model <- c(model, "}")

  ## define model class
  class(model) <- "prevModel"

  ## create data
  data <- list(x = x, n = n)

  ## generate inits
  inits <- NULL

  ## get results!
  if (verbose) cat("JAGS progress:\n\n")

  nodes <- paste0(c("SE", "SP"), rep(seq(h), each = 2))
  nodes <- c("TP", nodes, "bayesP")

  JAGSout <- R2JAGS(model = model, data = data, inits = inits,
                    nchains = nchains, burnin = burnin, update = update,
                    nodes = nodes, verbose = verbose)

  ## define mcmc samples
  mcmc.list <- JAGSout$mcmc.list
  class(mcmc.list) <- c("list", "mcmc.list")
  names <- colnames(mcmc.list[[1]])
  mcmc.list_list <- list()
  order <- c(length(names) - 1, c(t(cbind(1:h, 1:h+h))), length(names))
  for (i in seq_along(names))
    mcmc.list_list[[i]] <- mcmc.list[, order[i]]
  names(mcmc.list_list) <- names[order]

  ## define diagnostics
  # deviance information criterion
  DIC <- JAGSout$dic

  # bayes-p
  bayesP <- mean(unlist(mcmc.list_list$bayesP))

  # brooks-gelman-rubin diagnostic
  # exclude bayes-p and fixed nodes
  exclude <-
    c(which(colnames(mcmc.list[[1]]) == "bayesP"),
      which(apply(mcmc.list[[1]], 2, sd) == 0))
  BGR <- gelman.diag(mcmc.list[, -exclude])

  ## define output
  out <-
    new("prev",
        par = list(x = x, n = n, prior = prior,
                   nchains = nchains, burnin = burnin, update = update,
                   inits = inits),
        model = model,
        mcmc = mcmc.list_list,
        diagnostics = list(DIC = DIC,
                           BGR = BGR,
                           bayesP = bayesP))

  ## return output
  return(out)
}


## -------------------------------------------------------------------------#
## Create model for covariance scheme --------------------------------------#

truePrevMultinom_covariance <-
function(x, n, prior, nchains, burnin, update, verbose) {

  ## number of tests
  h <- log(length(x), 2)

  ## number of priors
  n_priors <- 1 + (2 * h) + (2 * sum(choose(h, seq(h, 2))))

  ## define model vector
  model <- character()

  ## write model initiation
  model[1] <- "model {"
  model[2] <- paste0("x[1:", 2 ^ h,
                     "] ~ dmulti(AP[1:", 2 ^ h, "], n)")

  ## write prob_se[], prob_sp[]
  s <- multiModel_select(h)[[1]]

  prob_se <- paste0("prob_se[", seq(nrow(s)), "] <-")

  for (i in seq(nrow(s))) {
    ## first element
    prob_se[i] <-
      paste(prob_se[i],
            paste(ifelse(s[i, ] == 1,
                         paste0("SE[", seq(ncol(s)), "]"),
                         paste0("(1 - SE[", seq(ncol(s)), "])")),
                  collapse = " * "))

    ## define index for 'a'
    a <- c(1, 0)
    
    ## h - 2 elements
    if (h > 2) {
      for (k in seq((h - 2), 1)) {
        se <-
          apply(t(apply(combn(h, k), 1, rev)),
                2,
                function(x) {
                  paste(ifelse(s[i, x] == 1,
                               paste0("SE[", x, "]"),
                               paste0("(1 - SE[", x, "])")),
                        collapse = " * ")
                })

        sign <-
          apply(t(apply(combn(h, k), 1, rev)),
                2,
                function(x) prod(2 * s[i, -x] - 1))  # convert (1,0) to (1,-1)

        a[2] <- a[2] + choose(h, k)
        prob_se[i] <-
          paste0(prob_se[i],
                 paste(
                   paste0(ifelse(sign == 1,
                                 " + ", " - "),
                          "a[", seq(a[1], a[2]), "] * ", se),
                   collapse = ""))
        a[1] <- a[2] + 1
      }
    }

    ## final element
    prob_se[i] <-
      paste(prob_se[i],
            paste0(ifelse(prod(2 * s[i, ] - 1) == 1,
                   "+ ","- "),
                   "a[", a[1], "]"))
  }

  prob_sp <- gsub("SE", "SP", prob_se)
  prob_sp <- gsub("prob_se", "prob_sp", prob_sp)
  prob_sp <- gsub("a", "b", prob_sp)
  for (i in seq_along(prob_sp)) {
    prob_sp[i] <-
      gsub(paste0("prob_sp[", i ,"]"),
           paste0("prob_sp[", 1 + length(prob_sp) - i ,"]"),
           fixed = TRUE,
           prob_sp[i])
  }

  model <-
    c(model,
      "",
      prob_se,
      "",
      prob_sp, "")

  ## write definition of AP and constraints
  model <-
    c(model,
      paste0("for (i in 1:", 2 ^ h, ") {"),
      "AP[i] <- TP * prob_se[i] + (1 - TP) * prob_sp[i]",
      "",
      write_constraint("AP", 1, 1),
      write_constraint("AP", 2, 2),
      write_constraint("prob_se", 1, 3),
      write_constraint("prob_se", 2, 4),
      write_constraint("prob_sp", 1, 5),
      write_constraint("prob_sp", 2, 6),
      "}",
      "")

  ## write prior
  priors <- get_nodes(h)

  for (i in seq(n_priors)) {
    model <-
      c(model,
        writeSeSp(prior[[i]][[1]], prior[[i]][[2]]))
  }

  ## write Bayes-P definition
  model <- c(model, "", write_bayesP(h))

  ## close model
  model <- c(model, "}")

  ## define model class
  class(model) <- "prevModel"

  ## create data
  data <- list(x = x, x2 = x, n = n,
               O1 = rep(1, 2 ^ h), O2 = rep(0, 2 ^ h),
               O3 = rep(1, 2 ^ h), O4 = rep(0, 2 ^ h),
               O5 = rep(1, 2 ^ h), O6 = rep(0, 2 ^ h))

  ## generate inits
  inits <- NULL

  ## get results!
  if (verbose) cat("JAGS progress:\n\n")

  nodes <- c("TP", "SE", "SP", "a", "b", "bayesP")

  JAGSout <- R2JAGS(model = model, data = data, inits = inits,
                    nchains = nchains, burnin = burnin, update = update,
                    nodes = nodes, verbose = verbose)

  ## define mcmc samples
  mcmc.list <- JAGSout$mcmc.list
  class(mcmc.list) <- c("list", "mcmc.list")
  names <- colnames(mcmc.list[[1]])
  if (h == 2) {
    names[which(names == "a")] <- "a[1]"
    names[which(names == "b")] <- "b[1]"
  }
  mcmc.list_list <- list()
  order <- match(c(priors, "bayesP"), names)
  for (i in seq_along(names))
    mcmc.list_list[[i]] <- mcmc.list[, order[i]]
  names(mcmc.list_list) <- names[order]

  ## define diagnostics
  # deviance information criterion
  DIC <- JAGSout$dic

  # bayes-p
  bayesP <- mean(unlist(mcmc.list_list$bayesP))

  # brooks-gelman-rubin diagnostic
  # exclude bayes-p and fixed nodes
  exclude <-
    c(which(colnames(mcmc.list[[1]]) == "bayesP"),
      which(apply(mcmc.list[[1]], 2, sd) == 0))
  BGR <- gelman.diag(mcmc.list[, -exclude])

  ## define output
  out <-
    new("prev",
        par = list(x = x, n = n, prior = prior, nchains = nchains,
                   burnin = burnin, update = update, inits = inits),
        model = model,
        mcmc = mcmc.list_list,
        diagnostics = list(DIC = DIC,
                           BGR = BGR,
                           bayesP = bayesP))

  ## return output
  return(out)
}
