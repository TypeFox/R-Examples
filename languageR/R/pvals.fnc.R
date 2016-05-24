pvals.fnc <- function(object, ...) {
  stop("MCMC sampling is no longer supported by lme4.
       For p-values, use the lmerTest package, which provides
       functions summary() and anova() which give p-values of
       various kinds.")
}
