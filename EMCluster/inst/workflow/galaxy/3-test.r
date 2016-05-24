rm(list = ls())

library(EMCluster)
# source("./R/information.r")
load("./data/galaxy.RndEM.rda")

tau <- 0.5
n.mc.base <- 1000

ret.test <- NULL
for(k0 in 1:5){
  for(ka in (k0+1):6){
    seed <- 1234 + k0 + ka
    set.seed(seed)
    n.mc <- n.mc.base * max(c(k0, ka))

    emobj0 <- ret.save[[k0]]
    emobja <- ret.save[[ka]]

    ll.0 <- emobj0$llhdval
    ll.a <- emobja$llhdval
    delta.hat <- emobja$llhdval - emobj0$llhdval
    E.delta <- get.E.delta(x, emobj0, emobja, tau = tau, n.mc = n.mc)
    E.chi2.0 <- get.E.chi2(x, emobj0, emobja, "0", tau = tau, n.mc = n.mc,
                           verbose = FALSE)
    E.chi2.a <- get.E.chi2(x, emobj0, emobja, "a", tau = tau, n.mc = n.mc,
                           verbose = FALSE)
    T <- 2 * (delta.hat - E.delta)

    pv0 <- pchisq.my(T, E.chi2.0[1], E.chi2.0[2], lower.tail = FALSE)
    pva <- pchisq.my(T, E.chi2.a[1], E.chi2.a[2], lower.tail = FALSE)
    pv <- pv0 * tau + pva * (1 - tau)

    cat("  H0: K = ", k0, " vs Ha: K = ", ka, "\n",
        "    ll0 = ", ll.0, ", lla = ", ll.a, "\n",
        "    df0 = ", E.chi2.0[1], ", nc0 = ", E.chi2.0[2],
        ", dfa = ", E.chi2.a[1], ", nca = ", E.chi2.a[2], "\n",
        "    delta.hat = ", delta.hat, ", E.delta = ", E.delta,
        ", T = ", T, "\n",
        "    pv0 = ", pv0, ", pva = ", pva, " pv = ", pv,
        "\n", sep = "")

    ret.test <- rbind(ret.test,
                      c(k0, ka, ll.0, ll.a, E.chi2.0, E.chi2.a, delta.hat,
                        E.delta, T, pv0, pva, pv, seed))
  }
}

colnames(ret.test) <- c("K0", "Ka", "ll0", "lla",
                        "df0", "nc0", "dfa", "nca",
                        "delta.hat", "E.delta", "T",
                        "pv0", "pva", "pv", "seed")

save(list = c("ret.test"), file = "./data/galaxy.test.rda")

print(ret.test)
