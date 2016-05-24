networkIRT <- function(.y,
                       .starts = NULL,
                       .priors = NULL,
                       .control = NULL,
                       .anchor_subject = NULL,
                       .anchor_item = NULL
                       ) {
    cl <- match.call()

    divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

    ## Default Control
    default <- list(threads = 1L,
                    verbose = FALSE,
                    maxit = 500,
                    thresh = 1e-6,
                    checkfreq = 50,
                    convtype = 1
                    )
    cat("\n")
    cat(divider)
    cat("networkIRT: Network IRT via Variational Inference\n\n")

    ret <- list()

    ## Main Call to Computation
    ret <- .Call('endorseIRT_estimate',
                 PACKAGE = 'emIRT',
                 .starts$alpha,
                 .starts$beta,
                 .starts$w,
                 .starts$theta,
                 .starts$gamma,
                 .y,
                 .priors$alpha$mu,
                 .priors$alpha$sigma,
                 .priors$beta$mu,
                 .priors$beta$sigma,
                 .priors$w$mu,
                 .priors$w$sigma,
                 .priors$theta$mu,
                 .priors$theta$sigma,
                 .priors$gamma$mu,
                 .priors$gamma$sigma,
                 ifelse(!is.null(.control$threads), .control$threads, default$threads),
                 ifelse(!is.null(.control$verbose), .control$verbose, default$verbose),
                 ifelse(!is.null(.control$maxit), .control$maxit, default$maxit),
                 ifelse(!is.null(.control$thresh), .control$thresh, default$thresh),
                 ifelse(!is.null(.control$checkfreq), .control$checkfreq, default$checkfreq),
                 ifelse(!is.null(.control$convtype), .control$convtype, default$convtype)
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


    ## ########
    ## Rotation
    ## ########

    if (!is.null(.anchor_subject)) {
        cat(paste("\n\tPolarity anchoring will use manually selected subject:",
                  "\n\t\t", .anchor_subject,
                  "\n",
                  sep = ""
                  )
            )
        cPivot <- mean(ret$means$theta)
        cAdj <- ifelse(ret$means$theta[.anchor_subject] < cPivot,
                       -1,
                       1
                       )

        ret$means$theta <- ret$means$theta * cAdj
        ret$means$w <- ret$means$w * cAdj

    } else if (!is.null(.anchor_item)) {
        cat(paste("\n\tPolarity anchoring will use manually selected item:",
                  "\n\t\t", .anchor_item,
                  "\n",
                  sep = ""
            )
            )

        cPivot <- mean(ret$means$w)
        cAdj <- ifelse(ret$means$w[.anchor_item] < cPivot,
                       -1,
                       1
                       )

        ret$means$w <- ret$means$w * cAdj
        ret$means$theta <- ret$means$theta * cAdj

    }

    cat(divider)

    ret$call <- cl

    class(ret) <- c("networkIRT", "emIRT")

    return(ret)
}
