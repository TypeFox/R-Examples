## User function to compute weighted average optima for taxa

`optima` <- function(x, ...)
    UseMethod("optima")

`optima.default` <- function(x, env, boot = FALSE, nboot = 1000,
                             alpha = 0.05, ...) {
    x <- data.matrix(x)
    opt <- colSums(x * env) / colSums(x)
    ## bootstrap?
    if (boot) {
        bootOpt <- function(x, env) {
            nr <- NROW(x)
            take <- sample(nr, nr, replace = TRUE)
            x <- x[take, ]           ## subset
            x <- x[, colSums(x) > 0] ## species with no occurences now
            env <- env[take]
            colSums(x * env) / colSums(x)
        }
        opt.star <- replicate(nboot, bootOpt(x, env))
        opt.bar <- rowMeans(opt.star, na.rm = TRUE)
        opt.sd <- apply(opt.star, 1, sd, na.rm = TRUE)
        ## transpose the quantile as apply returns probs *rows* and we
        ## want to arrange these as /cols/ in a matrix
        opt.pci <- t(apply(opt.star, 1, quantile,
                           probs = c(alpha / 2, 1 - (alpha / 2))))
        opt <- cbind(Optima = opt,
                     optBoot = opt.bar,
                     optSD = opt.sd,
                     opt.pci)
    }
    class(opt) <- "optima"
    if (is.matrix(opt)) {
        rownames(opt) <- colnames(x)
        class(opt) <- c(class(opt), "matrix")
    } else {
        names(opt) <- colnames(x)
    }
    attr(opt, "env") <- deparse(substitute(env))
    opt
}

`print.optima` <- function(x, ...) {
    cat("\n")
    msg <- paste("Weighted Average Optima For:", attr(x, "env"))
    writeLines(strwrap(msg, prefix = "\t"),
               sep = "\n\n")
    attr(x, "env") <- NULL
    print(unclass(x), ...)
}
