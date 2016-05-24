poisIRT <- function(.rc,
					i = 0:(nrow(.rc)-1),
					NI = nrow(.rc),
                    .starts = NULL,
                    .priors = NULL,
                    .control = NULL
                    ) {
    cl <- match.call()

    divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

    ## Default Control
    default <- list(threads = 1L,
                    verbose = FALSE,
                    maxit = 1500,
                    thresh = 1e-6,
                    checkfreq = 50
                    )
    cat("\n")
    cat(divider)
    cat("poisIRT: approximate Bayesian inference for textual data\n\n")

    ## Main Call to Computation
    ret <- .Call('poisIRT_estimate',
                 PACKAGE = 'emIRT',
                 .starts$alpha,
                 .starts$psi,
                 .starts$beta,
                 .starts$x,
                 .rc,
                 matrix(i,ncol=1),
                 NI,
                 .priors$psi$mu,
                 .priors$psi$sigma2,
                 .priors$alpha$mu,
                 .priors$alpha$sigma2,
                 .priors$beta$mu,
                 .priors$beta$sigma2,
                 .priors$x$mu,
                 .priors$x$sigma2,
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

    ## Labelling of Output
	rownames(ret$means$psi) <- rownames(.rc)
	rownames(ret$means$beta) <- rownames(.rc)
	rownames(ret$vars$beta) <- rownames(.rc)
	rownames(ret$means$alpha) <- colnames(.rc)

	#Would only work if NI=NK, so this is omitted
	#rownames(ret$means$x) <- colnames(.rc)
	#rownames(ret$vars$x) <- colnames(.rc)

    cat(divider)
    ret$call <- cl
    class(ret) <- c("poisIRT","emIRT")

    return(ret)
}
