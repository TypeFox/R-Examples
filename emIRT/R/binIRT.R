binIRT <- function(.rc,
                   .starts = NULL,
                   .priors = NULL,
                   .D = 1L,
                   .control = NULL,
                   .anchor_subject = NULL,
                   .anchor_outcomes = FALSE
                   ) {
    cl <- match.call()

    divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

    ## Default Control
    default <- list(threads = 1L,
                    verbose = FALSE,
                    maxit = 500,
                    convtype = 1,
                    thresh = 1e-6,
                    checkfreq = 50,
                    withLB = FALSE,
                    withProbs = FALSE,
                    asEM = FALSE
                    )
    cat("\n")
    cat(divider)
    cat("binIRT: Binary IRT via Expectation Maximization\n\n")

    ## Currently Force 1 Dimension
    if (.D != 1) stop("Only 1 dimension is currently supported. Use '.D = 1.'")

    ## Main Call to Computation
    ret <- .Call('FastEst_estimate',
                 PACKAGE = 'emIRT',
                 .starts$alpha,
                 .starts$beta,
                 .starts$x,
                 .rc$votes,
                 .priors$x$mu,
                 .priors$x$sigma,
                 .priors$beta$mu,
                 .priors$beta$sigma,
                 .D,
                 ifelse(!is.null(.control$threads), .control$threads, default$threads),
                 ifelse(!is.null(.control$verbose), .control$verbose, default$verbose),
                 ifelse(!is.null(.control$maxit), .control$maxit, default$maxit),
                 ifelse(!is.null(.control$convtype), .control$convtype, default$convtype),
                 ifelse(!is.null(.control$thresh), .control$thresh, default$thresh),
                 ifelse(!is.null(.control$checkfreq), .control$checkfreq, default$checkfreq),
                 ifelse(!is.null(.control$withLB), .control$withLB, default$withLB),
                 ifelse(!is.null(.control$withProbs), .control$withProbs, default$withProbs),
                 ifelse(!is.null(.control$asEM), .control$asEM, default$asEM)
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


    ret$means$beta <- cbind(ret$means$a, ret$means$b)
    ret$means$b <- NULL
    ret$means$a <- NULL

    .anchor_outcomes <- .anchor_outcomes[1]
    ## Rotation for 1 Dimension
    if (.D > 1L) {
        if (!is.null(.anchor_subject) | .anchor_outcomes) {
            cat(paste("\n\t",
                      "Polarity anchors: specified for more than 1 dimensions.",
                      "\n\t\tNot supported -- ignoring anchors.",
                      "\n",
                      sep = ""
                      )
                )
        }
    }

    if (!is.null(.anchor_subject) & .anchor_outcomes) {
        cat(paste("\n\t",
                  "Polarity anchors: both manual and outcome-based anchors specified.",
                  "\n\t\tAnchoring will use manually selected subjects.",
                  "\n",
                  sep = ""
                  )
            )

        idx <- .anchor_subject[1]
        thisx <- ret$means$x[idx, 1]
        cMean <- mean(ret$means$x[, 1])

        if (thisx < cMean) {
            ret$means$x[, 1] <- ret$means$x[, 1] * -1
            ret$means$beta[, 2] <- ret$means$beta[, 2] * -1
        }

    }

    if (!is.null(.anchor_subject) & !.anchor_outcomes) {
        cat(paste("\n\t",
                  "Polarity anchors: manually selected subjects specified.",
                  "\n",
                  sep = ""
                  )
            )

        idx <- .anchor_subject[1]
        thisx <- ret$means$x[idx, 1]
        cMean <- mean(ret$means$x[, 1])

        if (thisx < cMean) {
            ret$means$x[, 1] <- ret$means$x[, 1] * -1
            ret$means$beta[, 2] <- ret$means$beta[, 2] * -1
        }

    }

    if (is.null(.anchor_subject) & .anchor_outcomes) {
        cat(paste("\n\t",
                  "Polarity anchors: outcomes-based anchors specified.",
                  "\n",
                  sep = ""
                  )
            )

        nSucc <- apply(.rc$votes, 2, function (X) sum(X == 1))
        nFail <- apply(.rc$votes, 2, function (X) sum(X == -1))
                                        # number of successes and failures for
                                        # each item e.g., number of votes for,
                                        # number of correct responses

        ## print(nSucc)
        ## print(nFail)

        propFail <- (nFail) / (nFail + nSucc)
        cDiff <- propFail
                                        # model-free "difficulty"

        ## print(propFail)

        getScore <- function (X) {
            idxP <- which(X == 1)
            idxN <- which(X == -1)
            sum(propFail[idxP]) / sum(propFail[c(idxP, idxN)])
        }

        vScores <- apply(.rc$votes, 1, getScore)

        names(vScores) <- NULL
        ## print(vScores)

        cMax <- max(vScores, na.rm = TRUE)
        ## print(cMax)
        idx <- which(vScores == cMax)[1]

        ## print(idx)

        thisx <- ret$means$x[idx, 1]
        cMean <- mean(ret$means$x[, 1])

        if (thisx < cMean) {
            ret$means$x[, 1] <- ret$means$x[, 1] * -1
            ret$means$beta[, 2] <- ret$means$beta[, 2] * -1
        }


    }


    ## Labelling of Output
    dlx <- paste("d", 1:.D, sep = "")
    dla <- "d0"
    dlb <- paste("d", 1:.D, sep = "")

    rownames(ret$means$x) <- rownames(.rc$legis.data)
    colnames(ret$means$x) <- dlx
    colnames(ret$means$beta) <- c(dla, dlb)

    colnames(ret$vars$x) <- dlx
    rownames(ret$vars$x) <- dlx

    colnames(ret$vars$beta) <- c(dla, dlb)
    rownames(ret$vars$beta) <- c(dla, dlb)

    cat(divider)

    ret$call <- cl

    class(ret) <- c("binIRT", "emIRT")
    return(ret)
}
