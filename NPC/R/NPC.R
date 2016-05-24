#' @export
NPC <- function (data, ## Data frame with treatment, response, and other vars.
                 keep = TRUE, ## Subset of observations (default keeps all)
                 tr.var, ## Name of treatment variable (character)
                 tr.label, ## Level of 'tr.var' indicating treated units (char.)
                 y.vars, ## Names of response variables (character)
                 comb.fun = "ProductCF", ## Combining function (char. or fun.)
                 n.perms = 1000,   ## Number of permutations
                 block.var = NULL, ## Variable defining blocks within which to
                 ## restrict permutations (character)
                 clust.var = NULL, ## Variable defining clusters of observations
                 ## assigned to treatment en bloc (character)
                 event.var = NULL, ## Logical variable indicating whether duration
                 ## variables were observed rather than censored
                 ## (character)
                 alternative = "greater",  ## direction(s) of alternative
                 seed = 1, ## Random seed
                 na.rm = TRUE,     ## Delete observations with missing values?
                 FWE.adj = TRUE,   ## Calculate FWE-adjusted p-values
                 step.down = identical(comb.fun, "MinimumCF"), ## Adjust p-values
                 ## w/ stepdown MinP
                 test.statistic = "StudentsT", ## Test statistic(s)
                 return.matrix = FALSE, ## Return the permutation distribution of
                 ## test statistics and p-values?
                 print.steps = TRUE, ## Print progress of the function?
                 adapt.test = logical(length(y.vars)) ## Use the adaptive test
                 ## of Hogg et al. (1975)
                 ) {
  n.tests <- length(y.vars)
  ## (1) ERRORS
  stopifnot(length(alternative) == n.tests || length(alternative) == 1)
  stopifnot(all(alternative %in% c("two.sided", "less", "greater")))
  stopifnot(!FWE.adj || step.down || length(y.vars) < 15)
  ## (2) PREPARE FUNCTION ARGUMENTS
  ## Treatment variable(s)
  if (length(tr.label) == 1) tr.label <- rep(tr.label, n.tests)
  stopifnot(length(tr.label) == n.tests)
  ## Subset data
  keep.vars <- c(tr.var, y.vars, block.var, clust.var, event.var)
  if (na.rm) {
    newdata <- stats::na.omit(data[keep, keep.vars])
    n.drops <- nrow(data) - nrow(newdata)
    cat('\nNumber of observations deleted due to missingness:', n.drops, '\n')
    data <- newdata
  } else {
    data <- data[keep, keep.vars]
  }
  ## Variables
  Y <- subset(data,, y.vars)
  Tr <- factor(data[, tr.var])
  if (!is.null(event.var)) event.var <- data[, event.var]
  if (!is.null(block.var)) {
    block.var <- data[, block.var]
  } else {
    block.var <- gl(1, nrow(data))
  }
  if (!is.null(clust.var)) {
    clust.var <- data[, clust.var]
  } else {
    clust.var <- gl(nrow(data), 1)
  }
  clust.var <- droplevels(clust.var)
  if (is.character(comb.fun)) {
    CF <- eval(parse(text=comb.fun))
    stopifnot(is.function(CF))
  }
  ## Test statistic
  if (!is.list(test.statistic)) test.statistic <- as.list(test.statistic)
  if (length(test.statistic) == 1) {
    test.statistic <- rep(test.statistic, n.tests)
  }
  for (j in seq_along(test.statistic)) {
    if (is.character(test.statistic[[j]])) {
      test.statistic[[j]] <- eval(parse(text=test.statistic[[j]]))
    }
    if (adapt.test[j]) {
      if (length(unique(Y[, j])) < 0.1 * length(Y[, j])) {
        cat('Overriding choice of adaptive test because fewer than',
            '10% of observations are unique.\n', sep='\n')
      } else {
        test.statistic[[j]] <- HoggAdapt(y=Y[, j])
      }
    }
  }
  ## CALCULATE P-VALUES
  ## Distinguish clusters by block so permutations are independent across blocks
  blcl <- interaction(block.var, clust.var, drop = TRUE)
  cl.obs <- match(levels(blcl), blcl) ## first obs. in each cluster
  set.seed(seed)
  ## Permute clusters within blocks
  hw <- permute::how(blocks = block.var[cl.obs], maxperm = n.perms)
  cl.obs.shuffle <- permute::shuffleSet(n = length(cl.obs),
                                        nset = n.perms, control = hw)
  ## Add observed permutation
  cl.obs.shuffle <- rbind(seq_along(cl.obs), cl.obs.shuffle)
  O <- apply(cl.obs.shuffle, 1, function (i) Tr[cl.obs][i])
  Omega <- apply(O, 2, function (x) x[match(blcl, blcl[cl.obs])])
  actual.n.perms <- ncol(Omega)
  T0 <- vector("numeric", n.tests)
  TP <- matrix(NA, ncol = n.tests, nrow = actual.n.perms)
  for (j in 1:n.tests) {
    if (print.steps) {
      cat('\n', y.vars[j], '\n')
      print(test.statistic[[j]])
    }
    ## observed value of test statistic(s)
    stat.j <- test.statistic[[j]]
    T0[j] <- stat.j(y = Y[, j], tr = Tr, tl = tr.label[j], event = event.var,
                    block = block.var, clust = blcl)
    ## permutation distribution of test statistics
    TP[, j] <- apply(Omega, 2, function (z) {
      stat.j(y = Y[, j], tr = z, tl = tr.label[j], event = event.var,
             block = block.var, clust = blcl)
    })
  }
  if (!is.null(colnames(Y))) names(T0) <- colnames(TP) <- colnames(Y)
  if (print.steps) {
    cat("\nObserved test statistics:\n")
    print(T0)
  }
  T0P <- rbind(T0, TP)
  ## transform according to alternative
  alternative <- rep(alternative, length.out = n.tests)
  for (j in 1:n.tests) {
    if (alternative[j] == "two.sided") {
      T0P[, j] <- abs(T0P[, j])
    }
    if (alternative[j] == "less") {
      T0P[, j] <- -(T0P[, j])
    }    
  } 
  ## Partial p-values
  p.raw <- as.matrix(t2p(T0P))
  rownames(T0P) <- c("Observed", paste("Perm", 1:actual.n.perms))
  rownames(p.raw) <- c("Observed", paste("Perm", 1:actual.n.perms))
  if (is.null(colnames(Y))) {
    colnames(p.raw) <- paste0("p", 1:n.tests)
  } else {
    colnames(p.raw) <- names(T0)
  }
  ## NPC
  p.joint <- t2p(apply(p.raw, 1, CF, B=actual.n.perms))[1]
  ## FWE correction
  if (FWE.adj) {
    p.adj <- FWE(p.raw, step.down, cfun=CF)
    PV <- c(p.raw[1, ], p.adj, p.joint)
    names(PV) <- c(colnames(p.raw), paste(colnames(p.raw), "(adj.)"), "NPC")
  } else {
    PV <- c(p.raw[1, ], p.joint)
    names(PV) <- c(colnames(p.raw), "NPC")
  }
  out <- list(p.values=PV)
  if (return.matrix) out <- list(p.values=PV, p.matrix=p.raw, T.matrix=T0P)
  return(out)
}
