`weightedCor` <- function(x, ...) {
    UseMethod("weightedCor")
}

`weightedCor.default` <- function(x, env, fossil, method  = c("rda","cca"),
                                  test = TRUE, type = c("simulate","permute"),
                                  sim = 999, verbose = TRUE, ...) {
    ## convert to data matrices
    x <- data.matrix(x)
    fossil <- data.matrix(fossil)
    ## find the FUN
    method <- match.arg(method)
    FUN <- match.fun(method)
    ## match type
    type <- match.arg(type)
    ## intersection of species in x and fossil
    nams <- intersect(colnames(x), colnames(fossil))
    x <- x[, nams]
    fossil <- fossil[, nams]
    ## WA model and optima
    Args <- head(formals(wa.default), -1)
    dots <- list(x = x, y = env, env = NULL, ...)
    Args <- modifyList(Args, dots)
    Args <- lapply(Args,
                   function(x) if(typeof(x) == "language") {eval(x)[1]} else {x})
    mod <- do.call(waFit, Args)
    opt <- mod$wa.optima
    ## and predictions for fossil samples
    if(Args$tol.dw) {
        pArgs <- list(X = fossil, optima = opt, tol = mod$model.tol,
                      nr = nrow(fossil), nc = ncol(fossil))
        pFun <- WATpred
    } else {
        pArgs <- list(X = fossil, optima = opt)
        pFun <- WApred
    }
    pred <- drop(deshrinkPred(do.call(pFun, pArgs), mod$coefficients,
                              Args$deshrink))
    ## ordination
    ord <- rdaFit(fossil, pred)
    ## axis 1 species scores
    scrs <- drop(scores(ord, display = "species", scaling = 0,
                        choices = 1))
    ## mean abundance
    abund <- colMeans(x)

    ## weighted correlation
    wtdCorrel <- abs(cov.wt(X <- cbind(opt, scrs), wt = sqrt(abund),
                         cor = TRUE)$cor[1,2])
    Correl <- abs(cor(X)[1,2])

    ## test, either simulate or permute
    if(test) {
        if(type == "permute")
            warning("`type = \"permute\"` not currently implemented")
        if(type == "simulate") {
            if(verbose) {
                writeLines(paste("Simulating", sim, "Weighted Correlations:"))
                pb <- txtProgressBar(min = 0, max = sim, style = 3)
                on.exit(close(pb))
                on.exit(cat("\n"), add = TRUE)
            }
            wtdCorDist <- corDist <- numeric(length = sim)
            unifs <- matrix(runif(sim * nrow(x)), nrow = sim, byrow = TRUE)
            for(i in seq_len(sim)) {
                if(verbose)
                    setTxtProgressBar(pb, i)
                Args$y <- unifs[i,]
                rmod <- do.call(waFit, Args)
                ropt <- rmod$wa.optima
                pArgs$optima <- ropt
                rpred <- drop(deshrinkPred(do.call(pFun, pArgs),
                                           rmod$coefficients,
                                           Args$deshrink))
                rord <- rdaFit(fossil, rpred)
                rscrs <- drop(scores(rord, display = "species",
                                     scaling = 0,
                                     choices = 1))
                wtdCorDist[i] <-
                    abs(cov.wt(rX <- cbind(ropt, rscrs), wt = sqrt(abund),
                               cor = TRUE)$cor[1,2])
                corDist[i] <- abs(cor(rX)[1,2])
            }
        }
        ndist <- list(wtdCorDist = wtdCorDist, corDist = corDist)
    } else {
        ndist <- corDist <- wtdCorDist <- NULL
    }

    ## function call
    .call <- match.call()
    .call[[1]] <- as.name("weightedCor")

    ## return object
    retval <- list(wtdCorrel = wtdCorrel, Correl = Correl,
                   data = data.frame(Optima = opt, SppScores = scrs,
                   meanAbund = abund),
                   ord = ord, model = mod, method = class(ord),
                   ndist = ndist, sim = sim, type = type,
                   env = deparse(substitute(env)),
                   call = .call)
    class(retval) <- "weightedCor"
    retval
}

`plot.weightedCor` <- function(x, type = c("bubble", "null"),
                               weighted = TRUE,
                               size = 0.25,
                               xlab = paste(x$env, "WA Optima"),
                               ylab = "Axis 1 Score",
                               xlim,
                               main = "",
                               sub = NULL,
                               border = "gray75",
                               col = "gray75",
                               obscol = "red",
                               fg = "black", ...) {
    type <- match.arg(type)
    if(missing(sub))
        sub <- with(x, bquote(rho[w] == .(round(wtdCorrel, 3)) ~~~~
                              rho == .(round(Correl, 3))))
    if(type == "bubble") {
        with(x$data, symbols(Optima, SppScores, circles = meanAbund,
                             inches = size,
                             xlab = xlab, ylab = ylab, fg = fg, sub = sub, ...))
    }
    if(type == "null") {
        if(weighted) {
            ndist <- x$ndist$wtdCorDist
            obs <- x$wtdCorrel
        } else {
            ndist <- x$ndist$corDist
            obs <- x$Correl
        }
        if(missing(xlim))
            xlim <- range(0, ndist, obs)
        Dens <- density(ndist)
        Hist <- hist(ndist, plot = FALSE, ...)
        ylim <- range(0, Hist$density, Dens$y)
        plot(Hist, freq = FALSE, border = border, col = col, xlab = xlab,
             ylab = "Density", main = main, xlim = xlim, ylim = ylim,
             sub = sub, ...)
        abline(h = 0, col = col)
        lines(Dens)
        rug(obs, col = obscol, lwd = 2)
        rug(obs, col = obscol, lwd = 2, side = 3)
        box()
    }
    invisible(x)
}


`print.weightedCor` <- function(x, digits = 3, ...) {
    ## compute things we need
    wtdMsg <- paste("Weighted Correlation:", round(x$wtdCorrel, digits = digits))
    corMsg <- paste("Correlation         :", round(x$Correl, digits = digits))
    if(TEST <- !is.null(x$ndist)) {
        wtdMsg <- paste(wtdMsg, " (p = ",
                        format.pval(sum(x$ndist$wtdCorDist >= x$wtdCorrel) /
                                    (x$sim + 1)),
                        ")", sep = "")
        corMsg <- paste(corMsg, " (p = ",
                        format.pval(sum(x$ndist$CorDist >= x$Correl) / (x$sim +1)),
                        ")", sep = "")
        if(x$type == "simulate")
            type <- "Simulation"
        else
            type <- "Permutation"
        testMsg <- paste("Test Type           :", type)
    }

    writeLines(strwrap("Weighted correlation of WA Transfer Function",
                       prefix = "\n\t"))
    writeLines("\nCall:")
    print(x$call)
    writeLines("")
    if(TEST) {
        writeLines(testMsg)
    }
    writeLines(wtdMsg)
    writeLines(corMsg)
    writeLines("")
}
