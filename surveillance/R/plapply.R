################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Parallelized lapply (wrapping around mclapply and parLapply)
### taking care of the random seed and printing progress information
###
### Copyright (C) 2015 Sebastian Meyer
### $Revision: 1273 $
### $Date: 2015-03-21 17:39:13 +0100 (Sam, 21. Mär 2015) $
################################################################################


plapply <- function (X, FUN, ...,
                     .parallel = 1, .seed = NULL, .verbose = TRUE)
{
    if (!(useCluster <- inherits(.parallel, "cluster"))) {
        stopifnot(isScalar(.parallel), .parallel >= 1)
        .parallel <- as.vector(.parallel, mode = "integer")
        if (.Platform$OS.type == "windows" && .parallel > 1L) {
            useCluster <- TRUE
            .parallel <- parallel::makeCluster(.parallel)
            on.exit(parallel::stopCluster(.parallel))
        }
    }
    FUN <- match.fun(FUN)
    .FUN <- if (useCluster || is.primitive(FUN)) {
        FUN  # no support for reporting to the master || add.on.exit
    } else { # be verbose on.exit of FUN
        verboseExpr <- if (isTRUE(.verbose)) {
            ## progress bar or dots
            if (.parallel == 1L && interactive()) {
                env <- new.env(hash = FALSE, parent = environment(FUN))
                environment(FUN) <- env  # where the progress bar lives
                env$pb <- txtProgressBar(min = 0, max = length(X), initial = 0, style = 3)
                on.exit(close(env$pb), add = TRUE)
                quote(setTxtProgressBar(pb, pb$getVal() + 1L))
            } else {
                on.exit(cat("\n"), add = TRUE)
                quote(cat("."))
            }
        } else if (is.call(.verbose) || is.expression(.verbose)) {
            ## custom call or expression
            .verbose
        } else if (is.character(.verbose)) {
            ## custom progress symbol
            on.exit(cat("\n"), add = TRUE)
            substitute(cat(.verbose))
        } # else NULL (no output)
        ## add on.exit(verboseExpr) to body(FUN)
        do.call(add.on.exit, list(FUN, verboseExpr))
    }
    
    ## set random seed for reproducibility
    if (!is.null(.seed)) {
        if (useCluster) {
            parallel::clusterSetRNGStream(cl = .parallel, iseed = .seed)
        } else {
            if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                set.seed(NULL)  # initialize
            }
            .orig.seed <- get(".Random.seed", envir = .GlobalEnv)
            on.exit(assign(".Random.seed", .orig.seed, envir = .GlobalEnv),
                    add = TRUE)
            if (.parallel == 1L) {
                set.seed(seed = .seed)
            } else {
                stopifnot(requireNamespace("parallel", quietly = TRUE))
                ## Note @ R 3.1.3: this loading of package "parallel"
                ## before set.seed() is crucial; otherwise, the first run of
                ## plapply() would not be reproducible !!!
                set.seed(seed = .seed, kind = "L'Ecuyer-CMRG")
                parallel::mc.reset.stream()
            }
        }
    }

    ## rock'n'roll
    if (useCluster) {
        parallel::parLapply(cl = .parallel, X = X, fun = .FUN, ...)
    } else if (.parallel == 1L) {
        lapply(X = X, FUN = .FUN, ...)
    } else { # use forking
        parallel::mclapply(X = X, FUN = .FUN, ...,
                           mc.preschedule = TRUE, mc.set.seed = TRUE,
                           mc.silent = FALSE, mc.cores = .parallel)
    }
}


## add an on.exit() statement at the beginning of a function
add.on.exit <- function (FUN, expr)
{
    FUN <- match.fun(FUN)
    if (is.null(expr <- substitute(expr))) {
        return(FUN)
    }
    if (is.primitive(FUN)) { # body(FUN) is NULL
        stop("not implemented for primitive functions")
    }
    onexitexpr <- substitute(on.exit(expr))
    obody <- body(FUN)
    body(FUN) <- if (is.call(obody) && identical(as.name("{"), obody[[1L]])) {
        ## body(FUN) is a braced expression (usual case)
        ## and we insert on.exit(expr) directly after "{"
        as.call(append(x = as.list(obody), values = onexitexpr, after = 1L))
    } else {
        ## body(FUN) is a symbol or a single call like UseMethod("print")
        as.call(c(as.name("{"), onexitexpr, obody))
    }
    FUN
}
