bootmex <-
    # Bootstrap inference for a conditional multivaratiate extremes model.
function (x, R = 100, nPass = 3, trace = 10) {
    theCall <- match.call()
    if (class(x) != "mex"){
      stop("object must be of type 'mex'")
    }

# Construct the object to be returned.
    ans <- list()
    ans$call <- theCall

    getTran <- function(i, x, data, mod, th, qu, margins) {
        x <- c(x[, i])
        param <- mod[[i]]$coefficients
        th <- th[i]
        qu <- qu[i]
        data <- c(data[, i])
        if (margins == "gumbel"){
			     res <- revTransform(exp(-exp(-x)), data = data, th = th,
            					qu = qu, sigma = exp(param[1]), xi = param[2])
		    } else {
			    y <- x
			    y[x < 0] <- exp(x[x < 0]) / 2
			    y[x >= 0] <- 1 - exp(-x[x >= 0]) / 2
			    res <- revTransform(y, data = data, th = th,
					     			qu = qu, sigma = exp(param[1]), xi = param[2])
		    }
        res
    }

    mar <- x$margins
    dep <- x$dependence
    which <- dep$which
    constrain <- dep$constrain
    v <- dep$v
    dqu <- dep$dqu
    dth <- dep$dth
    margins <- dep$margins
    penalty <- mar$penalty
    priorParameters <- mar$priorParameters
    start <- 0.75* coef(x)$dependence[1:2,] # scale back towards zero in case point est on edge of original parameter space and falls off edge of constrained space for bootstrap sample

    n <- dim(mar$transformed)[[1]]
    d <- dim(mar$transformed)[[2]]
    dqu <- rep(dqu, d)
    dependent <- (1:d)[-which]

    ans$simpleDep <- dep$coefficients
    ans$dqu <- dqu
    ans$which <- which
    ans$R <- R
    ans$simpleMar <- mar
    ans$margins <- margins
    ans$constrain <- constrain

    innerFun <- function(i, x, which, dth, dqu, margins, penalty, priorParameters, constrain, v=v, start=start,
        pass = 1, trace = trace, n=n, d=d, getTran=getTran, dependent=dependent) {

        g <- sample(1:(dim(mar$transformed)[[1]]), size = n, replace = TRUE)
        g <- mar$transformed[g, ]
        ok <- FALSE

        while (!ok) {
          for (j in 1:(dim(g)[[2]])){
            u <- runif(nrow(g))
            if (margins == "gumbel"){
              g[order(g[, j]), j] <- sort(-log(-log(u)))
            } else {
              g[order(g[, j]), j] <- sort(sign(u - .5) * log(1 - 2*abs(u - .5)))
            }
          }
          if (sum(g[, which] > dth) > 1  &   all(g[g[,which] > dth , which] > 0)){ ok <- TRUE }
        }

        g <- sapply(1:d, getTran, x = g, data = mar$data, margins=margins,
                    mod = mar$models, th = mar$mth, qu = mar$mqu)

        dimnames(g)[[2]] <- names(mar$models)

        ggpd <- migpd(g, mth = mar$mth,
					  penalty = penalty, priorParameters = priorParameters)

        gd <- mexDependence(ggpd, dqu = dqu, which = which, margins=margins, constrain=constrain, v=v, start=start)
        res <- list(GPD = coef(ggpd)[3:4, ],
                    dependence = gd$dependence$coefficients,
                    Z = gd$dependence$Z,
                    Y = g)

        if (pass == 1) {
            if (i%%trace == 0) {
                cat(paste(i, "replicates done\n"))
            }
        }
        res
    } # Close innerFun

    res <- lapply(1:R, innerFun, x = x, which = which, dth = dth, margins=margins,
        dqu = dqu, penalty = penalty, priorParameters = priorParameters, constrain=constrain, v=v, start=start,
        pass = 1, trace = trace, getTran=getTran, n=n, d=d, dependent=dependent)

    # Sometimes samples contain no extreme values. Need to have another pass or two
    if (nPass > 1) {
        for (pass in 2:nPass) {
            rerun <- sapply(res, function(x) any(sapply(x, function(x) any(is.na(x)))))
            wh <- !unlist(lapply(res, function(x) dim(x$Z)[[1]] > 0))
            rerun <- apply(cbind(rerun, wh), 1, any)
            if (sum(rerun) > 0) {
                cat("Pass", pass, ":", sum(rerun), "samples to rerun.\n")
                rerun <- (1:R)[rerun]
                res[rerun] <- lapply((1:R)[rerun], innerFun,
                  x = x, which = which, dth = dth, dqu = dqu, margins=margins,
                  penalty = penalty, priorParameters = priorParameters, constrain=constrain, v=v, start=start,
                  pass = pass, trace = trace, getTran=getTran, n=n, d=d, dependent=dependent)
            }
        }
    }

    ans$boot <- res
    oldClass(ans) <- c("bootmex", "mex")
    ans
}

