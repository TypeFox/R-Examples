psyfun.boot <- function(obj, N = 100){
   n <- obj$prior.weights
   f <- fitted(obj)
   resp.bt <- matrix(rbinom(N * length(n), n, f), ncol = N)
   bt.res <- sapply(seq_len(N), function(x) {
      r <- resp.bt[, x]
      res.bt <- glm(cbind(r, n - r) ~ model.matrix(obj) - 1,
         binomial(obj$family$link))
      cc <- coef(res.bt)
      names(cc) <- names(coef(obj))
   cc
    })
    bt.res
}
