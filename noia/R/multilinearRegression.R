multilinearRegression <-
function (phen, gen = NULL, genZ = NULL, reference = "noia", 
    max.level = NULL, max.dom = NULL, fast = FALSE, e.unique = FALSE, 
    start.algo = "linear", start.values = NULL, robust = FALSE, 
    bilinear.steps = 1, ...) 
{
    "matrix2list" <- function(mat) {
        ans <- list()
        for (c in 1:ncol(mat)) {
            label <- colnames(mat)[c]
            ans[[label]] <- c(mat[, c])
        }
        return(ans)
    }
    if (is.null(max.level) || max.level > 2) {
        max.level <- 2
    }
    if (is.null(max.dom) || max.dom > max.level) {
        max.dom <- max.level
    }
    prep <- prepareRegression(phen = phen, gen = gen, genZ = genZ, 
        reference = reference, max.level = max.level, max.dom = max.dom, 
        fast = fast)
    if (is.null(start.values)) {
        start.values <- startingValues(prep$phen, prep$genZ, 
            reference = reference, max.level = max.level, max.dom = max.dom, 
            fast = fast, e.unique = e.unique, start.algo = start.algo, 
            bilinear.steps, ...)
    }
    ans <- NULL
    if (robust) {
        cat("Starting robust fitting procedure\n")
        cat("\tTrying starting algorithm: linear...")
        try(ans <- multilinearRegression(phen = phen, gen = gen, 
            reference = reference, genZ = genZ, max.level = max.level, 
            max.dom = max.dom, fast = fast, e.unique = e.unique, 
            start.algo = "linear", control = nls.control(maxiter = 500), 
            ...), silent = TRUE)
        if (is.null(ans)) {
            cat("FAILED\n")
            cat("\tTrying starting algorithm: multilinear...")
            try(ans <- multilinearRegression(phen = phen, gen = gen, 
                reference = reference, genZ = genZ, max.level = max.level, 
                max.dom = max.dom, fast = fast, e.unique = e.unique, 
                start.algo = "multilinear", control = nls.control(maxiter = 500), 
                ...), silent = TRUE)
            if (is.null(ans)) {
                cat("FAILED\n")
                cat("\tTrying starting algorithm: subset...")
                try(ans <- multilinearRegression(phen = phen, 
                  gen = gen, reference = reference, genZ = genZ, 
                  max.level = max.level, max.dom = max.dom, fast = fast, 
                  e.unique = e.unique, start.algo = "subset", 
                  control = nls.control(maxiter = 500), ...), 
                  silent = TRUE)
                if (is.null(ans)) {
                  cat("FAILED\n")
                  cat("\tTrying starting algorithm: bilinear...")
                  try(ans <- multilinearRegression(phen = phen, 
                    gen = gen, reference = reference, genZ = genZ, 
                    max.level = max.level, max.dom = max.dom, 
                    fast = fast, e.unique = e.unique, start.algo = "bilinear", 
                    bilinear.steps = 1, control = nls.control(maxiter = 500), 
                    ...), silent = TRUE)
                  if (is.null(ans)) {
                    cat("FAILED\n")
                    cat("\tTrying starting algorithm: bilinear2...")
                    try(ans <- multilinearRegression(phen = phen, 
                      gen = gen, reference = reference, genZ = genZ, 
                      max.level = max.level, max.dom = max.dom, 
                      fast = fast, e.unique = e.unique, start.algo = "bilinear", 
                      bilinear.steps = 2, control = nls.control(maxiter = 500), 
                      ...), silent = TRUE)
                    if (is.null(ans)) {
                      cat("FAILED\n")
                      cat("\tTrying starting algorithm: bilinear3...")
                      try(ans <- multilinearRegression(phen = phen, 
                        gen = gen, reference = reference, genZ = genZ, 
                        max.level = max.level, max.dom = max.dom, 
                        fast = fast, e.unique = e.unique, start.algo = "bilinear", 
                        bilinear.steps = 3, control = nls.control(maxiter = 500), 
                        ...), silent = TRUE)
                      if (is.null(ans)) {
                        cat("FAILED\n")
                      }
                      else {
                        cat("OK\n")
                      }
                    }
                    else {
                      cat("OK\n")
                    }
                  }
                  else {
                    cat("OK\n")
                  }
                }
                else {
                  cat("OK\n")
                }
            }
            else {
                cat("OK\n")
            }
        }
        else {
            cat("OK\n")
        }
    }
    else {
        nn <- effectsNamesMultilinear(prep$nloc, max.level, max.dom)
        form <- formulaMultilinear(prep$nloc, max.level, max.dom, 
            e.unique)
        X <- matrix2list(prep$x)
        phen <- prep$phen
        regression <- nls(formula = as.formula(form), start = start.values, 
            ...)
        ans <- prep
        ans$E <- coef(regression)
        ans$std.err <- summary(regression)$coef[, 2]
        ans$pvalues <- summary(regression)$coef[, 4]
        if (e.unique) {
            ee <- ans$E["ee"]
            ss <- ans$std.err["ee"]
            pp <- ans$pvalues["ee"]
            ans$E <- ans$E[names(ans$E) != "ee"]
            ans$std.err <- ans$std.err[names(ans$std.err) != 
                "ee"]
            ans$pvalues <- ans$pvalues[names(ans$pvalues) != 
                "ee"]
            for (i in 1:(prep$nloc - 1)) {
                for (j in (i + 1):prep$nloc) {
                  ans$E[paste("e", i, j, sep = "")] <- ee
                  ans$std.err[paste("e", i, j, sep = "")] <- ss
                  ans$pvalues[paste("e", i, j, sep = "")] <- pp
                }
            }
        }
        names(ans$E) <- nn
        names(ans$std.err) <- nn
        ans$resvar <- var(residuals(regression)) * (length(residuals(regression)) - 
            1)/length(residuals(regression))
        names(ans$pvalues) <- nn
        ans$variances <- rep(NA, length(nn))
        names(ans$variances) <- nn
        ans$regression <- regression
        class(ans) <- "noia.multilinear"
    }
    return(ans)
}
