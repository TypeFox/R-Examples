compliance.score <- function(
  D, Z, W,
  weights = NULL,
  link = qnorm,
  inv.link = pnorm,
  genoud = TRUE,
  num.iter = ifelse(genoud, 200, 1e4),
  one.sided = FALSE
  ) {

  if(!all.equal(sort(unique(D)), 0:1))
    stop("Treatment must be binary.")
  if(!all.equal(sort(unique(Z)), 0:1))
    stop("Assignment must be binary.")

  W <- as.matrix(W)
  if(length(D) != length(Z) || length(D) != nrow(W))
    stop("Number of observations must be equal for all variables.")

  n.obs <- length(D)
  if (is.null(weights))
    weights <- rep(1, n.obs)

  if (min(apply(W, 2, var)) == 0)
    stop("One or more predictors are constants.")

  # normalize all inputs
  W <- apply(W, 2, function(x) (x - mean(x)) / sd(x))
  W1 <- cbind(1, W)
  df <- data.frame(W = W)

  if (one.sided) {
    if (max(D[Z == 0]) == 1) {
      if (min(D[Z == 1]) == 0)
        stop("Non-compliance is not one-sided.")
      else {
        warning(paste("One-sided non-compliance not in expected direction.",
                      "Reversing coding of Z."))
        Z <- ifelse(Z == 1, 0, 1)
      }
    }

    fit <- glm(
      D ~ W,
      df, subset = (Z == 1),
      family = binomial(link="probit")
      )
    C.score <- predict(fit, newdata = df, type = "response")

    return(list(C.score = C.score, A.score = rep(0, n.obs), opt.log = NA))

  } else { # two-sided non-compliance
    if (genoud) {
      if(!requireNamespace("rgenoud", quietly = TRUE))
        stop("rgenoud package required to use genoud = TRUE option.")
    } else {
      if(!requireNamespace("minqa", quietly = TRUE))
        stop("minqa package required to use genoud = FALSE option.")
    }

    suppressWarnings({
      glm.CA <- glm(
        D ~ W, df, subset = (Z == 1),
        family = binomial(link="probit")
        )
      glm.A <- glm(
        D ~ W, df, subset = (Z == 0),
        family = binomial(link="probit")
        )
    })
    CA.score.cheap <- predict(glm.CA, df, type = "response")
    A.score.cheap <- predict(glm.A, df, type = "response")
    p.A.given.AC <- pmin(pmax(A.score.cheap / CA.score.cheap, 0.01), .99)
    suppressWarnings(
      lm.A.given.AC <- lm(link(p.A.given.AC) ~ W, weights = weights)
      )
    starting <- c(glm.CA$coefficients, lm.A.given.AC$coefficients)
    starting <- pmin(4.5, pmax(-4.5, starting))
    starting[is.nan(starting)] <- 0

    complier.neg.logLik.fnc <- function (par) {
      theta <- par[1:(length(par)/2)]
      theta.A <- par[(length(par)/2+1):(length(par))]

      p.AC <- inv.link(rowSums(t(t(W1) * theta)))
      p.A <- inv.link(rowSums(t(t(W1) * theta.A)))

      p.D <- p.AC * (1 - p.A) * Z + p.AC * p.A
      p.D <- pmax(pmin(p.D, 1 - 1e-5), 1e-5)

      return(-sum(log((p.D)^D * (1 - p.D)^(1 - D) * weights)))
    }

    # optimize
    if (genoud) {
      opt.log <- rgenoud::genoud(
        fn = complier.neg.logLik.fnc,
        nvars = (ncol(W1) * 2),
        starting.values = starting,
        pop.size = 100,
        Domains = 10,
        max.generations = num.iter,
        hessian = TRUE,
        print.level = 0,
        optim.method = "Nelder-Mead",
        hard.generation.limit = TRUE,
        gradient.check = FALSE
        )

    } else {
      opt.log <- minqa::bobyqa(
        starting, complier.neg.logLik.fnc,
        control = list(maxit = num.iter, trace = 0)
        )
      if (opt.log$convergence > 0) warning("Nelder-Mead did not converge!")
    }

    CA.score <- inv.link(
      rowSums(t(t(W1) * opt.log$par[1:(length(opt.log$par)/2)]))
      )
    A.score <- inv.link(
      rowSums(t(
        t(W1) * opt.log$par[(length(opt.log$par)/2+1):(length(opt.log$par))]))
      ) * CA.score
    C.score <- CA.score - A.score
    return(list(C.score = C.score, A.score = A.score, opt.log = opt.log))
  }
}
