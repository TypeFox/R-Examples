dynIRT <- function(.data,
                    .starts = NULL,
                    .priors = NULL,
                    .control = NULL
                    ) {
    cl <- match.call()

    divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

    ## Default Control
    default <- list(threads = 1L,
                    verbose = FALSE,
                    maxit = 500,
                    thresh = 1e-6,
                    checkfreq = 50
                    )
    cat("\n")
    cat(divider)
    cat("dynIRT: Dynamic IRT via Variational Inference\n\n")

    ## Main Call to Computation
    ret <- .Call('dynIRT_estimate',
                 PACKAGE = 'emIRT',
                 .starts$alpha,
                 .starts$beta,
                 .starts$x,
                 .data$rc,
                 .data$startlegis,
                 .data$endlegis,
                 .data$bill.session,
                 .data$T,
                 .priors$x.mu0,
                 .priors$x.sigma0,
                 .priors$beta.mu,
                 .priors$beta.sigma,
                 .priors$omega2,
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

	rownames(ret$means$x) <- rownames(.data$rc)
	rownames(ret$vars$x) <- rownames(.data$rc)

	rownames(ret$means$alpha) <- colnames(.data$rc)
	rownames(ret$means$beta) <- colnames(.data$rc)

	rownames(ret$vars$alpha) <- colnames(.data$rc)
	rownames(ret$vars$beta) <- colnames(.data$rc)

    cat(divider)
    ret$call <- cl
    class(ret) <- "dynIRT"
    return(ret)
}
