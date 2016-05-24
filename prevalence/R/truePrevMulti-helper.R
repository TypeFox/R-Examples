###=========================================================================#
### TRUE PREVALENCE FROM MULTIPLE TESTS / helper functions
###=========================================================================#

###=========================================================================#
###== FUNCTIONS ============================================================#
###-- checkMultiPrior_conditional ..... check prior for cond prob scheme
###-- checkMultiPrior_covariance ...... check prior for covariance scheme
###---| explode ....................... main explode function
###----| explode_theta ................ explode theta (cond prob)
###----| explode_nodes ................ explode nodes (covariance)
###----| explode_operator ............. explode operator
###----| explode_dist ................. explode distribution
###---| get_nodes ..................... define nodes of conditional model

###-- multiModel_select ............... define theta construct for SE/SP
###-- multiModel_probs ................ define AP[] in terms of theta[]
###-- multiModel_SeSp ................. expressions for TP, SE[] and SP[]
###---| multiModel_build ..............
###---| multiModel_collapse ...........

###-- write_bayesP .................... write definition of Bayes-P
###-- write_bayesP0 ................... write definition of Bayes-P
###-- write_constraint ................ write constraint on prob_se, prob_sp


## -------------------------------------------------------------------------#
## Check prior for conditional probability scheme --------------------------#

checkMultiPrior_conditional <-
function(prior) {
  ## evaluate whether prior is defined correctly
  first_element <- as.character(prior)[1]
  if (!any(c("{", "list") == first_element))
    stop("'prior' is not defined correctly")

  ## if prior is defined as a list
  ## note: list element names currently not taken in account!
  if (first_element == "list") {
    n <- length(prior) - 1
    priors_list <- vector("list", n)
    for (i in seq(n))
      priors_list[[i]] <-
        checkSeSp(eval(parse(text = prior)[[i + 1]]),
                  type = "prob")
  }

  ## if prior is defined as a function
  if (first_element == "{") {
    n <- length(prior) - 1
    priors_list0 <- vector("list", n)
    for (i in seq(n)) {
      priors_list0[[i]] <-
        explode(as.character(prior[[i + 1]]), "conditional")
    }

    ## get indices from priors_list0
    index <- sapply(priors_list0, function(x) as.numeric(x[[1]]))

    ## check if all indices exist
    if (length(index) != max(index)) {
      stop("The indices of 'theta[.]' are not correctly specified.\n",
           "See ?theta for more info.")
    }
    if (!all(unique(index) == index)) {
      stop("The indices of 'theta[.]' are not correctly specified.\n",
           "See ?theta for more info.")
    }
    if (length(index) == 1 | log(length(index) + 1, 2)%%1 != 0) {
      stop("The number of specified theta values is incorrect.\n",
           "See ?theta for more info.")
    }

    ## re-arrange list elements if needed
    order <- order(index)
    priors_list <- vector("list", n)
    for (i in seq(n))
      priors_list[[i]] <- priors_list0[[order[i]]][[2]]
    
  }

  ## return prior in list format
  return(priors_list)
}


## -------------------------------------------------------------------------#
## Check prior for covariance scheme ---------------------------------------#

checkMultiPrior_covariance <-
function(prior, h) {
  ## evaluate whether prior is defined correctly
  first_element <- as.character(prior)[1]
  if (!any(c("{", "list") == first_element))
    stop("'prior' is not defined correctly")

  ## based on h tests, these priors are expected
  priors <- get_nodes(h)

  ## if prior is defined as a list
  ## note: list element names currently not taken in account!
  if (first_element == "list") {
    n <- length(prior) - 1

    # check if length is as expected
    if (n != length(priors))
      stop("'prior' is not defined correctly; ",
           "expected ", length(priors), " elements, ",
           "got ", n)

    # check if names are as expected
    prior_names <- names(eval(prior))
    priors_names <- gsub("\\]", "", gsub("\\[", "", priors))
    if (!all(prior_names %in% priors_names))
      stop("'prior' is not defined correctly; ",
           "expected priors are ", paste(priors_names, collapse = ", "))
    if (!all(priors_names %in% prior_names))
      stop("'prior' is not defined correctly; ",
           "expected priors are ", paste(priors_names, collapse = ", "))

    # check priors and put in list
    priors_list <- vector("list", n)
    for (i in seq(n)) {
      xi <-
        suppressWarnings(
          as.numeric(substr(prior_names[i], 2, nchar(prior_names[i]))))
      type <-
        ifelse(substr(prior_names[i], 1, 1) %in% c("a", "b"),
               cov_depth(h, xi),
               "prob")
      priors_list[[i]] <-
        list(prior_names[i],
             checkSeSp(eval(parse(text = prior)[[i + 1]]),
                       type = type))
    }

    ## re-arrange list elements if needed
    order <- match(prior_names, priors_names)
    priors_list <- priors_list[order]

    ## rename priors (add brackets)
    for (i in seq(n))
      priors_list[[i]][[1]] <- priors[i]

  ## if prior is defined as a function
  } else if (first_element == "{") {
    n <- length(prior) - 1
    priors_list <- vector("list", n)
    for (i in seq(n)) {
      priors_list[[i]] <-
        explode(as.character(prior[[i + 1]]), "covariance", h)
    }

    ## check if all priors are defined
    priors_nodes <- sapply(priors_list, function(x) x[[1]])
    if (!all(priors_nodes %in% priors))
      stop("'prior' is not defined correctly; ",
           "expected priors are ", paste(priors_names, collapse = ", "))
    if (!all(priors %in% priors_nodes))
      stop("'prior' is not defined correctly; ",
           "expected priors are ", paste(priors_names, collapse = ", "))

    ## re-arrange list elements if needed
    order <- match(priors_nodes, priors)
    priors_list <- priors_list[order]
  }

  ## return prior in list format
  return(priors_list)
}


## -------------------------------------------------------------------------#
## Check depth of covariance parameter -------------------------------------#

cov_depth <-
function(h, x) {
  n_test <- rev(seq(h, 2))
  n_comb <- choose(h, n_test)
  rep(n_test, n_comb)[x]
}


## -------------------------------------------------------------------------#
## Main explode function ---------------------------------------------------#

explode <-
function(x, method, h = NULL) {
  ## create list of 2 (node & dist)
  priors <- vector("list", 2)

  ## extract node
  priors[[1]] <-
    switch(method,
           conditional = explode_theta(x[2]),
           covariance = explode_nodes(x[2]))

  ## check if operator is correctly specified
  explode_operator(x[1])

  ## define type
  xi <-
    suppressWarnings(
      as.numeric(substr(priors[[1]], 3, nchar(priors[[1]])-1)))
  type <-
    ifelse(substr(priors[[1]], 1, 1) %in% c("a", "b"),
           cov_depth(h, xi),
           "prob")

  ## extract distribution
  priors[[2]] <- explode_dist(x[3], type)
  return(priors)
}


## -------------------------------------------------------------------------#
## explode theta (conditional probability scheme) --------------------------#

explode_theta <-
function(x) {
  ## find 'theta[]'
  if (length(grep("theta", x, fixed = TRUE)) != 1)
    stop("Priors must be defined as vector 'theta'")
  if (length(grep("[", x, fixed = TRUE)) != 1 |
      length(grep("]", x, fixed = TRUE)) != 1)
    stop("The different values of theta must be defined as 'theta[.]'")

  ## extract '.' in 'theta[.]'
  x <- strsplit(x, "theta[", fixed = TRUE)[[1]][2]
  theta <- strsplit(x, "]", fixed = TRUE)[[1]][1]

  ## theta should be an integer
  if (!is.numeric(theta) && as.numeric(theta) %% 1 != 0)
    stop("'theta[.]' not specified correctly")

  return(theta)
}


## -------------------------------------------------------------------------#
## explode nodes (covariance scheme) ---------------------------------------#

explode_nodes <-
function (x) {
  if (!any(c(grepl("TP", x, fixed = TRUE),
             grepl("SE", x, fixed = TRUE),
             grepl("SP", x, fixed = TRUE),
             grepl("a", x, fixed = TRUE),
             grepl("b", x, fixed = TRUE)))) {
      stop("Priors must be named 'TP', 'SE', 'SP', 'a' or 'b'")
  }

  if (!(strsplit(x, "[", fixed = T)[[1]][1] %in%
        c("TP", "SE", "SP", "a", "b"))) {
      stop("Priors must be named 'TP', 'SE', 'SP', 'a' or 'b'")
  }

  if (x != "TP" &&
      !all(c(grepl("[", x, fixed = TRUE),
             grepl("]", x, fixed = TRUE)))) {
      stop("Priors must be defined as vectors")
  }

  if (x != "TP") {
    rhs <- strsplit(x, "[", fixed = TRUE)[[1]][2]
    i <- strsplit(rhs, "]", fixed = TRUE)[[1]][1]
    if (!grepl("^[[:digit:]]+$", i) || i < 1) {
      stop("Prior '", strsplit(x, "[", fixed = T)[[1]][1],
           "' not correctly indexed")
    }
  }

  return(x)
}


## -------------------------------------------------------------------------#
## explode operator --------------------------------------------------------#

explode_operator <-
function(operator) {
  ## operator should be '~' or '<-' or '='
  if (!any(c("~", "<-", "=") == operator))
    stop("Operator should be either '~', '<-' or '='")
}


## -------------------------------------------------------------------------#
## explode distribution ----------------------------------------------------#

explode_dist <-
function(x, type) {
  d <- dist2list(x, type)
  return(d)
}


## -------------------------------------------------------------------------#
## Define nodes of conditional model ---------------------------------------#

get_nodes <-
function(h) {
  nodes <-
    c("TP",
      paste0("SE[", seq(h), "]"),
      paste0("SP[", seq(h), "]"),
      paste0("a[", seq(sum(choose(h, seq(h, 2)))), "]"),
      paste0("b[", seq(sum(choose(h, seq(h, 2)))), "]"))
  return(nodes)
}


## -------------------------------------------------------------------------#
## Define theta construct for SE/SP ----------------------------------------#

## out01 defines construct of SE/SP: 0=(1-theta[.]), 1=theta[.]
## outSE and outSP define which thetas in expression for SE and SP

multiModel_select <-
function(n){
  out01 <- array(dim = c(2^n, n))
  outSE <- array(dim = c(2^n, n))
  outSP <- array(dim = c(2^n, n))
  for (i in seq(n)){
    out01[, i] <- rev(rep(c(0, 1), each = 2^(n-i), times = 2^(i-1)))
    outSE[, i] <- rev(rep(c((2^i+2^(i-1)-1):(2^i)), each = 2^(n-i+1)))
    outSP[, i] <- rev(rep(c((2^i+2^(i-1)):(2^(i+1)-1)), each = 2^(n-i+1)))
  }
  return(list(out01, outSE, outSP))
}


## -------------------------------------------------------------------------#
## Define AP[] in terms of theta[] -----------------------------------------#

multiModel_probs <-
function(s) {
  p <- character(dim(s[[1]])[1])
  for (i in seq(dim(s[[1]])[1])) {
    p[i] <- paste0("AP[", i, "] <- theta[1]")
    for (j in seq(dim(s[[1]])[2])) {
      p[i] <-
        paste0(p[i],
               ifelse(s[[1]][i,j] == 1, "*theta[", "*(1-theta["),
               s[[2]][i,j],
               ifelse(s[[1]][i,j] == 1, "]", "])"))
    }
    p[i] <- paste0(p[i], " + (1-theta[1])")
    for (j in seq(dim(s[[1]])[2])) {
      p[i] <-
        paste0(p[i],
               ifelse(s[[1]][i,j] == 0, "*theta[", "*(1-theta["),
               s[[3]][i,j],
               ifelse(s[[1]][i,j] == 0, "]", "])"))
    }
  }

  return(p)
}


## -------------------------------------------------------------------------#
## Expressions for TP, SE[] and SP[] ---------------------------------------#

multiModel_SeSp <-
function(n){
  TPSESP <- character(1 + 2 * n)
  TPSESP[1] <- "TP <- theta[1]"
  TPSESP[2] <- "SE1 <- theta[2]"
  TPSESP[3] <- "SP1 <- theta[3]"

  if (n > 1){
    for (i in 2:n){
      buildSE <- multiModel_build(i, SP = FALSE)
      buildSP <- multiModel_build(i, SP = TRUE)
      if (i > 2){
        for (j in i:3){
          buildSE <- multiModel_collapse(buildSE, j, SP = FALSE)
          buildSP <- multiModel_collapse(buildSP, j, SP = TRUE)
        }
      }
      TPSESP[(2*i)]   <- paste0("SE", i, " <- ", buildSE)
      TPSESP[(2*i)+1] <- paste0("SP", i, " <- ", buildSP)
    }
  }
  return(TPSESP)
}

multiModel_build <-
function(n, SP){
  N <- 2^(n-1)
  Next <- c(2*N, 2*N+1)
  out <- character(N/2)
  for (i in seq(N/2)){
    out[i] <-
      paste0("theta[", N+(i-1)+(SP*N/2),
             "] * theta[", Next[1]+(SP*N),
             "] + (1-theta[", N+(i-1)+(SP*N/2), 
             "]) * theta[", Next[2]+(SP*N),
             "]")
    Next <- Next + 2
  }
  return(out)
}


multiModel_collapse <-
function(build, n, SP){
  N <- length(build) / 2
  Next <- 2^(n-2)
  out <- character(N)
  for (i in seq(N)){
    ii <- (2*i)-1
    out[i] <-
      paste0("theta[", Next+(i-1)+(SP*Next/2),
             "] * (", build[ii],
             ") + (1-theta[", Next+(i-1)+(SP*Next/2),
             "]) * (", build[ii+1],
             ")")
  }
  return(out)
}


## -------------------------------------------------------------------------#
## Write definition of Bayes-P ---------------------------------------------#

write_bayesP <-
function(h) {
  bayesP <-
    c(paste0("x2[1:", (2^h), "] ~ dmulti(AP[1:", (2^h), "], n)"),
      #"d1 <- pow(x[] - AP[] * n, 2) / (n * AP[] * (1 - AP[]))",
      #"d2 <- pow(x2[] - AP[] * n, 2) / (n * AP[] * (1 - AP[]))",
      paste0("for (i in 1:", (2^h), ") {"),
      "d1[i] <- x[i] * log(max(x[i],1) / (AP[i]*n))",
      "d2[i] <- x2[i] * log(max(x2[i],1) / (AP[i]*n))",
      "}",
      "G0 <- sum(d1[])",
      "Gt <- sum(d2[])",
      "bayesP <- step(G0 - Gt)")
  return(bayesP)
}

write_bayesP0 <-
function(h) {
  bayesP <-
    c(paste0("x2[1:", (2^h), "] ~ dmulti(AP[1:", (2^h), "], n)"),
      paste0("for (i in 1:", (2^h), ") {"),
      "d1[i] <- x[i] * log(max(x[i],1) / (AP[i]*n))",
      "d2[i] <- x2[i] * log(max(x2[i],1) / (AP[i]*n))",
      "}",
      "G0 <- 2 * sum(d1[])",
      "Gt <- 2 * sum(d2[])",
      "bayesP <- step(G0 - Gt)")
  return(bayesP)
}


## -------------------------------------------------------------------------#
## Write constraint on 'prob_se', 'prob_sp' --------------------------------#

write_constraint <-
function(node, constraint, i) {
  add <- ifelse (constraint == 2, " - 1", "")
  constr <-
    c(paste0("constraint", i, "[i] <- step(", node, "[i]", add, ")"),
      paste0("O", i, "[i] ~ dbern(constraint", i, "[i])"))
  return(constr)
}