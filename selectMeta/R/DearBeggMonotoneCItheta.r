DearBeggMonotoneCItheta <- function(res, lam = 2, conf.level = 0.95, maxiter = 500){

print("Now computing confidence interval. May take a moment (up to a few minutes).")

## now profile over theta
CIlow <- uniroot(f = DearBeggProfileLL, interval = c(-20, res$theta), res0 = res, lam = lam, conf.level = conf.level, maxiter = maxiter)$root
CIup <- uniroot(f = DearBeggProfileLL, interval = c(res$theta, 20), res0 = res, lam = lam, conf.level = conf.level, maxiter = maxiter)$root

res <- list("ci.theta" = c(CIlow, CIup))
return(res)
}
