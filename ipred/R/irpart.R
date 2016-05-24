#
#  use update to fit multiple trees to bootstrap samples
#
irpart <- function(formula, data=NULL, weights, subset,
		   na.action=na.rpart, method, model=FALSE, x=FALSE, y=TRUE,
		   parms, control, cost, bcontrol, ...)
{

    mc <- match.call()
    mc$bcontrol <- NULL
    mc[[1]] <- as.name("rpart")

    m <- match.call(expand.dots=FALSE)
    m$model <- m$method <- m$control <- m$bcontrol <- NULL
    m$x <- m$y <- m$parms <- m$... <- NULL
    m$cost <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame.default")
    m <- eval(m, parent.frame())

    init_tree <- eval(mc, parent.frame())
    nobs <- length(init_tree$where)
    if (missing(weights)) { 
        weights <- rep(1.0, nobs)
    } else {
        warning("weights argument ignored in irpart")
    }

    yclasses <- c(class = "sclass", exp = "ssurv", anova = "sreg", poisson = "sreg")

    # 
    # Bagging: repeat this several times!
    #

    if (is.null(bcontrol)) stop("bcontrol not given")
    mod <- vector(mode="list", length=bcontrol$nbagg)

    for (b in 1:bcontrol$nbagg) {
        if (bcontrol$nbagg > 1)
            bindx <- sample(1:nobs, bcontrol$ns, replace=bcontrol$replace)
        else
            bindx <- 1:nobs
        tab <- tabulate(bindx, nbins = nobs)

        mc$data <- m[bindx,,drop = FALSE] ### tab * weights
        ans <- eval(mc, parent.frame())

        # return the appropriate class
        this <- list(bindx = bindx, btree = ans)

        class(this) <- yclasses[init_tree$method]

        mod[[b]] <- this
        }
    mod
}

