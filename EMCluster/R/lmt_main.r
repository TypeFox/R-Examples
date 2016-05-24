# For LMT.

lmt <- function(emobj.0, emobj.a, x, tau = 0.5, n.mc.E.delta = 1000,
    n.mc.E.chi2 = 1000, verbose = FALSE){
  if(class(emobj.0) != "emret" || class(emobj.a) != "emret"){
    stop("emobj.0 and emobj.a should be both in \"emret\" class.")
  }
  if(emobj.0$nclass == emobj.a$nclass){
    stop("emobj.0 and emobj.a should have different numbers of clusters.")
  }

  # K
  k.0 <- emobj.0$nclass
  k.a <- emobj.a$nclass

  # logL
  ll.0 <- emobj.0$llhdval
  ll.a <- emobj.a$llhdval

  # Likelihood ratio statistics
  delta.hat <- emobj.a$llhdval - emobj.0$llhdval
  E.delta <- get.E.delta(x, emobj.0, emobj.a, tau = tau, n.mc = n.mc.E.delta)

  # Chi-squared statistics
  E.chi2.0 <- get.E.chi2(x, emobj.0, emobj.a, "0", tau = tau,
                         n.mc = n.mc.E.chi2, verbose = verbose)
  E.chi2.a <- get.E.chi2(x, emobj.0, emobj.a, "a", tau = tau,
                         n.mc = n.mc.E.chi2, verbose = verbose)

  # Testing statistics.
  T <- 2 * (delta.hat - E.delta)

  # p-values.
  pv.0 <- pchisq.my(T, E.chi2.0[1], E.chi2.0[2], lower.tail = FALSE)
  pv.a <- pchisq.my(T, E.chi2.a[1], E.chi2.a[2], lower.tail = FALSE)
  pv <- pv.0 * tau + pv.a * (1 - tau)

  # Return.
  ret <- list(K.0 = k.0, K.a = k.a, ll.0 = ll.0, ll.a = ll.a,
              delta.hat = delta.hat,
              E.delta = E.delta, E.chi2.0 = E.chi2.0, E.chi2.a = E.chi2.a,
              T = T, pv.0 = pv.0, pv.a = pv.a, pv = pv)
  class(ret) <- "lmt"
  ret
} # End of lmt().

print.lmt <- function(x, digits = max(4, getOption("digits") - 3), ...){
  K.0 <- x$K.0
  K.a <- x$K.a
  ll.0 <- x$ll.0
  ll.a <- x$ll.a
  delta.hat <- x$delta.hat
  E.delta <- x$E.delta
  E.chi2.0 <- x$E.chi2.0
  E.chi2.a <- x$E.chi2.a
  T <- x$T
  pv.0 <- x$pv.0
  pv.a <- x$pv.a
  pv <- x$pv

  cat("- H.0: K = ", K.0, " vs H.a: K = ", K.a, "\n",
      "    ll.0 = ", ll.0, ", ll.a = ", ll.a, "\n",
      "    df.0 = ", E.chi2.0[1], ", nc.0 = ", E.chi2.0[2],
      ", df.a = ", E.chi2.a[1], ", nc.a = ", E.chi2.a[2], "\n",
      "    delta.hat = ", delta.hat, ", E.delta = ", E.delta,
      ", T = ", T, "\n",
      "    pv.0 = ", pv.0, ", pv.a = ", pv.a, " pv = ", pv,
      "\n", sep = "")

  invisible()
} # End of print.lmt().
