LogrankA <-
function(surv, group, weight) {
  # are all necessary arguments specified?
  if (missing(surv)) {
    stop("No survival argument specified")
  }
  if (missing(group)) {
    stop("No group argument specified")
  }
  # is the survival argument a survival object?
  if (is.Surv(surv) == FALSE) {
    stop("Survival argument not a survival object")
  }
  
  # extract time and status from survival object
  time <- surv[ , 1]
  status <- surv[ , 2]

  # if weight is not specified, default to a dummy weight vector of 1s
  if (missing(weight)) {
    weight <- rep(1, length(time))
    # otherwise, check if weight is numeric...
  } else if (!is.numeric(weight)) {
    stop("Non-numeric weight value found")
    # ...and non-negative...
  } else if (any(na.omit(weight) < 0)) {
    stop("Negative weight value found")
    # ...and integer (check for existence fractional part)
  } else if (any(na.omit(weight) %% 1 > 0)) {
    stop("Non-integer weight value found")
  }
  # is the censoring of type "right"?
  if (attr(surv, "type") != "right") {
    stop("Censorings are not of type right")
  }
  # are all arguments of equal length?
  if (any(c(length(time), length(status),
            length(group), length(weight)) != length(time))) {
    stop("Arguments differ in length")
  }
  # is more than one group specified?
  if (length(unique(group)) == 1) {
    stop("Only a single group specified")
  }
  
  # drop observations containing missing values
  survdata <- data.frame(time, status, group, weight)
  if (any(is.na(survdata))) {
    na.dropped <- na.omit(survdata);
    time <- na.dropped$time;
    status <- na.dropped$status;
    group <- na.dropped$group;
    weight <- na.dropped$weight;
    n.dropped <- sum(survdata$weight) - sum(na.dropped$weight)
  } else {
    n.dropped <- 0
  }
  
  # risk population at start of process time
  r.t0 <- sum(weight)
  # number of observations (events and censored cases) for each time interval i
  w <- sort(unique(weight))
    # each count of weight gets multiplicated with corresponding
    # value for weight. result is the number of occurrences
  obs.i <- colSums(table(weight, time) * w)
  # risk population at the beginning of each time interval i
  RevCumsumRev <- function (x) {rev(cumsum(rev(x)))}
  r.i <- RevCumsumRev(obs.i)
  # number of group-independent events in each time interval i
    # censored cases excluded from the total number of observations by
    # multiplying the status (0 for censored cases) with the weight
  sw <- status * weight
  d <- sort(unique(sw))
  d.i <- colSums(table(sw, time) * d)
  
  # number of observations for each time interval i by group k
  obs.ki <- colSums(table(weight, time, group) * w)
  # risk population of group k at the beginning of each time interval i
  r.ki <- apply(obs.ki, 2, RevCumsumRev)
  # expected events of group k in each time interval i
  e.ki <- d.i * r.ki / r.i
  # total observed events group k
  O.k <- colSums(table(sw, group) * d)
  # total expected events group k
  E.k <- colSums(e.ki)
  
  # actual logrank test statistic
  LR.k <- (O.k - E.k)^2 / E.k
  LR <- sum(LR.k)
  # degrees of freedom for chi^2 test of LR
  df <- length(O.k) - 1
  # chi^2 test of LR
  p.chi2 <- pchisq(LR, df = df, lower.tail = FALSE)
  
  # output
  r.kt0 <- t(t(r.ki[1, ]))
  logrank.parameter <- cbind(r.kt0, t(t(O.k)), t(t(E.k)), LR.k)
  colnames(logrank.parameter) <- c("N", "Obs. events", "Exp. events",
                                   "(O-E)^2/E")
  cat("Valid observations: ", r.t0, "\n",
      "Dropped observations: ", n.dropped, "\n\n",
      "Logrank test statistic = ", LR, " on ", df,
      " degrees of freedom,\np = ", p.chi2, "\n\n", sep = "")
  list(p.chi2 = p.chi2, df = df, LR = LR,
       lr.parameter = logrank.parameter)
}
