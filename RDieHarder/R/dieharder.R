## $Date$
## $Id$

dieharderGenerators <- function() {
    val <- .Call("dieharderGenerators", PACKAGE="RDieHarder")
    return(data.frame(names=val[[1]], id=val[[2]]))
}

dieharderTests <- function() {
    val <- .Call("dieharderTests", PACKAGE="RDieHarder")
    return(data.frame(names=val[[1]], id=val[[2]]))
}

dieharder <- function(rng="mt19937",
                      test=1,
                      psamples=100,
                      seed=0,
                      verbose=FALSE,
                      inputfile="",
                      ntuple=5) {
    UseMethod("dieharder")
}

dieharder.default <- function(rng="mt19937",
                              test="diehard_runs",
                              psamples=100,
                              seed=0,
                              verbose=FALSE,
                              inputfile="",
                              ntuple=5) {

    if (length(rng) > 1) {
        warning("Only one rng argument supported in dieharder")
        return(NULL)
    }

    if (is.character(rng)) {
        genpos <- charmatch(rng, as.character(.dieharder.generators$names))
    } else {
        genpos = which(.dieharder.generators[,"id"] == rng)
    }
    if (length(genpos)==0 || is.na(genpos)) {
        warning("rng argument ", rng, " unknown")
        return(NULL)
    }
    if (genpos == 0 && rng != 1) {
        warning("rng argument ", rng, " ambiguous")
        return(NULL)
    }
    gen <- .dieharder.generators$id[genpos]
    ##cat("Genpos: ", genpos, " Gen: ", gen, "\n", sep="")

    if (length(test) > 1) {
        warning("Only one test argument supported in dieharder")
        return(NULL)
    }
    if (is.character(test)) {
        runtestpos <- charmatch(test, as.character(.dieharder.tests$names))
        if (is.na(runtestpos)) {
            warning("test argumement ", test, " unknown")
            return(NULL)
        }
        if (runtestpos == 0) {
            warning("test argumment ", test, " ambiguous")
            return(NULL)
        }
        runtest <- .dieharder.tests$id[runtestpos]
    } else {
        runtest <- test
    }

    val <- .Call("dieharder",
                 as.integer(gen),
                 as.integer(runtest),
                 as.integer(seed),
                 as.integer(psamples),
                 #as.integer(rngdraws),
                 as.integer(verbose),
                 as.character(inputfile),
                 as.integer(ntuple),
                 PACKAGE="RDieHarder")
    obj <- list(p.value=val[[1]][1],
                data=val[[2]],      ## not used by htest methods
                method=val[[3]],
                data.name=paste("Created by RNG `",
                                .dieharder.generators[genpos,"names"], "' with seed=",
                                as.integer(seed), ", sample of size ",
                                as.integer(psamples), sep=""),
                generator=as.character(.dieharder.generators[genpos,"names"]),
                nvalues=val[[4]],
                p.values=val[[1]]
                )
    class(obj) <- c("dieharder", "htest")
    return(obj)
}

plot.dieharder <- function(x, ...) {

    local.par <- par(mfrow=c(2,1), mar=c(2,4,3,1), oma=c(0,0,3.5,0))

    vec <- x$data

    hist(vec, probability=TRUE, main="Histogram and Density estimate",
         xlab="", ylab="density")
    lines(density(x$data, from=0, to=1), col='darkgray')

    ##qqplot(vec, seq(0, 1, length.out=length(vec)),
    ##       main="QQ-Plot", ylab="Uniform sequence", xlab="")

    plot(ecdf(x$data), main="ECDF", ylab="", xlab="",
         verticals= TRUE, do.p = FALSE)
    segments(0,0,1,1, col='darkgray', lty="dotted")

    mtext(text = x$method, outer = TRUE, cex = 1.2, font = 2, line = 2)
    mtext(text = x$data.name, outer = TRUE, cex = 1.0, font = 1, line = 1)
    pksk <- round(x$p.value, 4)
    pks <- round(ks.test(x$data, "punif", 0, 1, exact=TRUE)$p.value, 4)
    pw <- round(wilcox.test(x$data, mu=0.5)$p.value, 4)
    mtext(text = paste("Test p-values: ",
          pksk, " (Kuiper-K-S),  ",
          pks, " (K-S),  ",
          pw, "  (Wilcoxon)"),
          outer = TRUE, cex = 1.0, font = 1, line = 0)

    par(local.par)
    invisible(x)
}

print.dieharder <- function(x, ...) {
    z <- x
    class(z) <- "htest"
    print(z)
    invisible(x)
}

summary.dieharder <- function(object, ...) {
    print(object, ...)
    cat("\nSummary for test data\n")
    print(summary(object$data))
    if (any(!is.na(object$data))) {
        cat("\n\nStem and leaf plot for test data\n")
        print(stem(object$data))
        print(ks.test(object$data, "punif", 0, 1, exact=TRUE))
        print(wilcox.test(object$data, mu=0.5))
    }
    invisible(object)
}

##

#dieharderVec <- function(gen=c(22,69), test=1, verbose=FALSE) {
#
#  val <- .Call("dieharderVectorised",
#                as.list(as.integer(gen)),
#                as.integer(test),
#                as.integer(verbose),
#                PACKAGE="RDieHarder")
#   return(val)
#}


