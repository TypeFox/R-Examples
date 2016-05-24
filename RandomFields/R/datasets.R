
CHOICE <- c("low", "medium", "high")
TAILCHOICE <- c("compact", "exponential", "power")
UNIVARIATE <- 1

## make sure that the results can be reobtained in the future
## -> use explicite the method for simiulation, not RPgauss

## get peoples suggestion on nu(x) and also other
## non-stationary cov fcts

## can covariates be involved?

## non-Gaussian marginal distribution: empirically based; Box-Cox based?!

## Gneiting, whittle, Cauchy, ex(?!) multivariate verallg. von cauchy?

xRFdatasets <- function(nonstationarity=CHOICE,
                       trendnonstationarity=CHOICE,
                       anisotropy=CHOICE,
                       differentiability=CHOICE,
                       tail = TAILCHOICE,
                       grid= c(FALSE, TRUE),
                       
                       locations = CHOICE, ## or number
                       spacedimension=1:3,
                       multivariate = UNIVARIATE,
                       time=c(FALSE, TRUE),

                       holdout_points,

                       covariates,

                       trend,
                       
                       n = 1,                       
                       seed = 0,
                       marginal= "Gaussian"
                      ) {
  res <- array(dim=c(length(nonstationarity),
                 length(trendnonstationarity),
                 length(anisotropy),
                 length(differentiability),
                 length(tail),
                 length(grid),
                 n
                 ))



 # Bsplines

}
                       
