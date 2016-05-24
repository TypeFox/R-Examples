getStarts <- function(.N,
                      .J,
                      .D,
                      .type = "zeros"
                      ) {

    if (.type == "zeros"){
        starts <- list(alpha = {matrix(0,
                                       nrow = .J,
                                       ncol = 1
                                       )
                            },
                       beta = {matrix(0,
                                      nrow = .J,
                                      ncol = .D
                                      )
                           },
                       x = {matrix(rnorm(.N * .D),
                                   nrow = .N,
                                   ncol = .D
                                   )
                        }
                       )
    } else if (.type == "random") {
        starts <- list(alpha = {matrix(rnorm(.J * 1) * 10,
                                       nrow = .J,
                                       ncol = 1
                                       )
                            },
                       beta = {matrix(rnorm(.J * .D) * 10,
                                      nrow = .J,
                                      ncol = .D
                                      )
                               },
                       x = {matrix(rnorm(.N * .D) * 1,
                                   nrow = .N,
                                   ncol = .D
                                   )
                        }
                       )
    } else {
        stop("Unknown type.")
    }


    return(starts)
}



getStartsE <- function(.N,
                       .J,
                       .y
                       ) {


    vAlpha <- qnorm(apply(.y, 2, mean))
    vBeta <- qnorm(apply(.y, 1, mean))

    vAlpha[vAlpha == -Inf] <- -2
    vAlpha[vAlpha == Inf] <- 2
    vBeta[vBeta == -Inf] <- -2
    vBeta[vBeta == Inf] <- 2

    mS1 <- matrix(vAlpha,
                  nrow = .N,
                  ncol = .J,
                  byrow = TRUE
                  )

    mS2 <- matrix(vBeta,
                  nrow = .N,
                  ncol = .J,
                  byrow = FALSE
                  )

    mS <- mS1 + mS2
    mSP1 <- pnorm(mS)
    mSP0 <- 1 - mSP1

    ## .y2 <- .y

    .y2 <- (.y * 2 - 1)

    q1 <- (1 - mSP1)
                                        # unexp prob of choosing 1

    q2 <- (1 - mSP0)
                                        # unexp prob of choosing 0

    q <- q1 * (.y == 1) - q2 * (.y == 0)

    q <- q2 * (.y == 0)

    ## plot(.y2, q)
    ## plot(lData$dist, q)

    ## .y2 <- lData$mu - mS




    ## dN <- dist(.y2, "euclidean")
    ## dJ <- dist(t(.y2), "euclidean")

    ## vTheta <- scale(cmdscale(dN, 1))
    ## vW <- scale(cmdscale(dJ, 1))


    ## plot(lData$alpha, vAlpha)
    ## plot(lData$beta, vBeta)

    ## plot(lData$theta, vTheta)
    ## cor(lData$theta, vTheta)

    ## plot(lData$w, vW)
    ## cor(lData$w, vW)

    starts <- list(alpha = {matrix(vAlpha,
                                   nrow = .J,
                                   ncol = 1
                                   )
                        },
                   beta = {matrix(vBeta,
                                  nrow = .N,
                                  ncol = 1
                                  )
                       },
                   w = {matrix(rnorm(.J, sd = .1),
                               nrow = .J,
                               ncol = 1
                               )
                    },
                   theta = {matrix(rnorm(.N, sd = .1),
                                   nrow = .N,
                                   ncol = 1
                                   )
                        },
                   gamma = {matrix(1,
                                   nrow = 1,
                                   ncol = 1
                                   )
                        }
                   )

    return(starts)
}
