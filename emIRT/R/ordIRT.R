ordIRT <- function(.rc,
                    .starts = NULL,
                    .priors = NULL,
                    .D = 1L,
                    .control = NULL
                    ) {
    cl <- match.call()

    ## Input change: Instead of rollcall() object, .rc is now a sparse matrix

    divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

    ## Default Control
    default <- list(threads = 1L,
                    verbose = FALSE,
                    maxit = 300,
                    thresh = 1e-6,
                    checkfreq = 50
                    )
    cat("\n")
    cat(divider)
    cat("ordIRT: Ordinal IRT via Expectation Maximization\n\n")

    ## Main Call to Computation
    ret <- .Call('ordIRT_estimate',
                 PACKAGE = 'emIRT',
                 .starts$tau,
                 .starts$DD,
                 .starts$beta,
                 .starts$x,
                 .rc,
                 .priors$x$mu,
                 .priors$x$sigma,
                 .priors$beta$mu,
                 .priors$beta$sigma,
                 1L,
                 ifelse(!is.null(.control$threads), .control$threads, default$threads),
                 ifelse(!is.null(.control$verbose), .control$verbose, default$verbose),
                 ifelse(!is.null(.control$maxit), .control$maxit, default$maxit),
                 ifelse(!is.null(.control$thresh), .control$thresh, default$thresh),
                 ifelse(!is.null(.control$checkfreq), .control$checkfreq, default$checkfreq)
                 )

    cat(paste("\t",
              "Done in ",
              ret$runtime$iters,
              " iterations, using ",
              ret$runtime$threads,
              " threads.",
              "\n",
              sep = ""
              )
        )

	rownames(ret$means$x) <- rownames(.rc)
	rownames(ret$vars$x) <- rownames(.rc)

	rownames(ret$means$beta) <- colnames(.rc)
	rownames(ret$vars$beta) <- colnames(.rc)

	rownames(ret$means$tau) <- colnames(.rc)
	rownames(ret$vars$tau) <- colnames(.rc)

	rownames(ret$means$Delta_sq) <- colnames(.rc)
	rownames(ret$means$Delta) <- colnames(.rc)

    cat(divider)
    ret$call <- cl
    class(ret) <- c("ordIRT", "emIRT")
    return(ret)
}
