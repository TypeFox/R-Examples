hierIRT <- function(.data = NULL,
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
    cat("hierIRT: approximate Bayesian inference for roll call data\n\n")

    ## Main Call to Computation
    ret <- .Call('hierIRT_estimate',
                 PACKAGE = 'emIRT',
                 .starts$alpha,
                 .starts$beta,
                 .starts$gamma,
                 .starts$sigma,
                 .starts$eta,
                 .data$y,
                 .data$z,
                 .data$g,
                 .data$i,
                 .data$j,
                 .priors$gamma.mu,
                 .priors$gamma.sigma,
                 .priors$beta.mu,
                 .priors$beta.sigma,
                 .priors$sigma.v,
                 .priors$sigma.s,
                 .data$ND,
                 .data$NG,
                 .data$NI,
                 .data$NJ,
                 .data$NL,
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

    cat(divider)
    ret$call <- cl
    class(ret) <- c("hierIRT", "emIRT")
    return(ret)
}
