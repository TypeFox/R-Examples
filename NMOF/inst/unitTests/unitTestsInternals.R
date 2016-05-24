## -*- truncate-lines: t; -*-
## internals

test.checkList <- function() {
    ## error for all opt-functions: NA in list
    OF <- function(x) NULL
    algo <- list(printBar = NA)
    checkException(DEopt(OF = OF, algo = algo), silent = TRUE)
    checkException(GAopt(OF = OF, algo = algo), silent = TRUE)
    checkException(LSopt(OF = OF, algo = algo), silent = TRUE)
    checkException(PSopt(OF = OF, algo = algo), silent = TRUE)
    checkException(TAopt(OF = OF, algo = algo), silent = TRUE)

    ## warning changed to error
    op <- options("warn")
    on.exit(options(op))
    options(warn=2)
    algo <- list(test = 42, nB = 100)
    checkException(DEopt(OF = OF, algo = algo), silent = TRUE)
    checkException(GAopt(OF = OF, algo = algo), silent = TRUE)
    checkException(LSopt(OF = OF, algo = algo), silent = TRUE)
    checkException(PSopt(OF = OF, algo = algo), silent = TRUE)
    checkException(TAopt(OF = OF, algo = algo), silent = TRUE)
}
