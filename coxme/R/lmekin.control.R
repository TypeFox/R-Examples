#
# Gather all of the control parameters for lmekin into one spot
#
lmekin.control <- function(optpar=list(method='BFGS', 
                                      control=list(reltol=1e-8)),
                           varinit=c(.02, .1, .8, 1.5)^2,
                           corinit = c(0, .3)) {
    if (length(varinit)<1 || !is.numeric(varinit))
        stop("varinit must be a vector of numeric values")
    if (length(corinit)<1 || !is.numeric(corinit))
        stop("corinit must be a vector of numeric values")
    if (any(varinit <=0)) stop ("varinit values must be >0")
    if (any(corinit <0))  stop("corinit values must be >=0")

    # The sparsity issues don't apply to lmekin models.
    # Set the option so that all sparsity is assumed: any
    #  factor with 1 or more levels, and level of the factor
    #  that accounts for 100% or less of the level.
    # Since coxme and lmekin share the routine that sets up
    #  the random effects design matrices, we have to set this.
    list(optpar=optpar, varinit=varinit, corinit=corinit,
         sparse=c(1,1))
    }
