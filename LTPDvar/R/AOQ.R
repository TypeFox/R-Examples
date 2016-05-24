AOQ <- function (p, n, k, N, type = c("exact", "napprox","ewmaSK","ewma2"),lam=1)
(1 - n/N) * p * OC(p, n, k, type,lam)
