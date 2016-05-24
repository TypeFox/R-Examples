cloglog.sample.size <- function(p.alt, n = NULL, p = 0.5,
                                power = 0.80, alpha = 0.05,
                                alternative = c("two.sided", "greater", "less"),
                                exact.n = FALSE,
                                recompute.power = FALSE, phi = 1) {
  compute.power <- function(Zalpha, p, p2, n, delta, alternative, phi, ...) {
    gamma1 <- log(-log(p))
    gamma2 <- log(-log(p2))
    sd1 <- sqrt(phi * var.cloglog(p, n))
    sd2 <- sqrt(phi * var.cloglog(p2, n))
    pz1 <- pnorm((- sd1 * Zalpha - gamma1 + gamma2)/sd2) # decrease in p
    pz2 <- pnorm((- sd1 * Zalpha + gamma1 - gamma2)/sd2) # increase in p
    power <- switch(alternative,
                    less      = ifelse(delta == 0, 1 - pnorm(Zalpha), pz1),
                    greater   = ifelse(delta == 0, 1 - pnorm(Zalpha), pz2),
                    two.sided = ifelse(delta == 0, 2 * (1 - pnorm(Zalpha)), pz1 + pz2))
    power
  }
  if(is.null(version$language))
    assign("compute.power", compute.power)
  compute.delta <- function(Zalpha, Zpower, p, n, alternative, phi, ...) {
    find.delta <- function(delta, Zalpha, power, p, n, alternative, phi)
      compute.power(Zalpha, p, p + delta, n, delta, alternative, phi) - power
    power <- pnorm(Zpower)
    lower <- if(alternative == "less") - p + 1e-008 else rep(0, length(Zalpha))
    upper <- if(alternative == "less") rep(0,length(Zalpha)) else 1 - p - 1e-008
    delta <- rep(NA, length(Zalpha))
    for(i in seq(along = Zalpha))
      delta[i] <- uniroot(f = find.delta,
                          lower = lower[i], upper = upper[i],
                          Zalpha = Zalpha[i],
                          power = power[i],
                          p = p[i], n = n[i],
                          alternative = alternative,
                          phi = phi[i])$root
    delta
  }
  compute.sample.size <- function(Zalpha, Zpower, p, p2, delta, exact.n, phi, ...) {
    gamma1 <- log(-log(p))
    gamma2 <- log(-log(p2))
    sd1 <- sqrt(phi * var.cloglog(p))
    sd2 <- sqrt(phi * var.cloglog(p2))
    n <- ((sd1 * Zalpha + sd2 * Zpower)/(gamma2 - gamma1))^2
    if(!exact.n) n <- ceiling(n)
    else n <- pmax(n, 1)
    n
  }
  if(!missing(alternative)) {
    alt.expanded <- match.arg(alternative)
  } else {
    alt.expanded <- alternative[1]
  }
  ###
  ### Determine what to compute and set that arg to NA for passing to 
  ### build.power.table.  The as.numeric(NA) is to prevent as.data.frame
  ### from converting to a factor (a bug).
  ###
  if(is.null(n)) {
    compute.what <- "sample.size"
    compute.function <- "compute.sample.size"
    n <- as.numeric(NA)
  } else if((missing(p.alt) || is.null(p.alt))) {
    compute.what <- "delta"
    compute.function <- "compute.delta"
    delta <- as.numeric(NA)
    p2 <- as.numeric(NA)
  } else {
    compute.what <- "power"
    compute.function <- "compute.power"
    power <- as.numeric(NA)
  }
  if(compute.what != "delta") p2 <- p.alt
  arg.names <- c("p", "p2", "delta", "alpha", "power", "n", "phi")
  table.names <- c("p.null", "p.alt", "delta", "alpha", "power", "n", "phi")
  power.table <- expand.grid(p = p, p2 = p2, alpha = alpha,
                             power = power, n = n, phi = phi)
  power.table$delta <- power.table$p2 - power.table$p
  power.table <- power.table[, arg.names]
  if(alt.expanded == "two.sided") {
    Zalpha <- qnorm(1 - power.table$alpha/2)
  } else {
    Zalpha <- qnorm(1 - power.table$alpha)
  }
  if(!all(is.na(power.table$power))) {
    Zpower <- qnorm(power.table$power)
  } else {
    Zpower <- as.numeric(NA)
  }
  arglist <- list(Zalpha = Zalpha, Zpower = Zpower,
                  exact.n = exact.n, alternative = alt.expanded)
  arglist <- c(power.table, arglist)
  if(compute.what == "sample.size") compute.what <- "n"
  power.table[, compute.what] <- do.call(compute.function, arglist)
  if(recompute.power && !exact.n && compute.function == "compute.sample.size") {
    arglist[[compute.what]] <- power.table[, compute.what]
    power.table[, "power"] <- do.call("compute.power", arglist)
  }
  if(compute.function == "compute.delta")
    power.table$p2 <- power.table$p + power.table$delta
  names(power.table) <- table.names
  power.table
}
