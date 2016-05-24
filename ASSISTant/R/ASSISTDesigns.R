#' A class to encapsulate the adaptive clinical trial design of Lai, Lavori and Liao
#'
#' @description \code{ASSISTDesign} objects are used to design, simulate and analyze
#' adaptive group sequential clinical trial with three stages.
#'
#' @docType class
#' @seealso \code{LLL.SETTINGS} for an explanation of trial parameters
#' @importFrom R6 R6Class
#' @importFrom mvtnorm pmvnorm Miwa
#' @importFrom stats uniroot rnorm pnorm qnorm
#' @section Methods:
#'
#' \describe{
#'   \item{\code{ASSISTDesign$new(designParameters, trialParameters, generateData)}}{Create a new
#'         \code{ASSISTDesign} instance object using the parameters specified}
#'   \item{\code{getDesignParameters},\code{getTrialParameters},
#'         \code{getBoundaries}}{Accessor methods for (obvious) object fields}
#'   \item{\code{print()}}{Print the object in a human readable form}
#'   \item{\code{computeCriticalValues()}}{Compute the critical boundary values \eqn{\tilde{b}},
#'         \eqn{b} and \eqn{c} for futility, efficacy and final efficacy decisions; saved in field
#'         \code{boundaries}}
#'   \item{\code{explore(numberOfSimulations = 5000, rngSeed = 12345, effectiveParameters = self$getDesignParameters(), showProgress = TRUE)}}{Explore the
#'         design using the specified number of simulations and random number seed. \code{trueParameters} is by default the same
#'         as \code{designParameters} as would be the case for a Type I error calculation. If changed, would yield power.
#'         Show progress if so desired. Returns a data frame of results}
#'   \item{\code{analyze(trialHistory)}}{Analyze
#'         the design given the \code{trialHistory} which is the result of a call to \code{explore} to
#'         simulate the design. Return a list of summary quantities}
#'   \item{\code{summary(analysis)}}{Print the operating characteristics of the design, using the analysis
#'         result from the \code{analyze} call}
#' }
#'
#' @references Adaptive Choice of Patient Subgroup for Comparing Two Treatments
#' by Tze Leung Lai and Philip W. Lavori and Olivia Yueh-Wen Liao. Contemporary Clinical Trials,
#' Vol. 39, No. 2, pp 191-200 (2014). \url{http://www.sciencedirect.com/science/article/pii/S1551714414001311}
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @examples
#' \dontrun{
#' data(LLL.SETTINGS)
#' prevalence <- LLL.SETTINGS$prevalences$table1
#' scenario <- LLL.SETTINGS$scenarios$S0
#' designParameters <- list(prevalence = prevalence,
#'                        mean = scenario$mean,
#'                        sd = scenario$sd)
#' designA <- ASSISTDesign$new(trialParameters = LLL.SETTINGS$trialParameters,
#'                             designParameters = designParameters)
#' print(designA)
#' ## A realistic design uses 5000 simulations or more!
#' result <- designA$explore(showProgress = interactive())
#' analysis <- designA$analyze(result)
#' designA$summary(analysis)
#' }
#' ## For full examples, try:
#' ## browseURL(system.file("full_doc/ASSISTant.html", package="ASSISTant"))
#'
ASSISTDesign <- R6Class("ASSISTDesign",
                        private = list(
                            designParameters = NA,
                            trialParameters = NA,
                            boundaries = NA,
                            checkParameters = function(designParameters, trialParameters, generateData) {
                                if (length(trialParameters$N) != 3) {
                                    stop("Improper sample size vector; this design assumes three interim looks")
                                }
                                if (!integerInRange(trialParameters$N, low = 1)) {
                                    stop("Improper values for sample sizes")
                                }
                                if (!identical(order(trialParameters$N), seq_along(trialParameters$N))) {
                                    stop("Improper values for sample sizes; need increasing sequence")
                                }
                                if (!scalarInRange(designParameters$J, low = 3, high = 10)) {
                                    stop("Improper number of subgroups; need at least 3; max 10")
                                }

                                if (!scalarInRange(trialParameters$type1Error, low=0.0001, 0.2)) {
                                    stop("Improper type 1 error")
                                }
                                if (!scalarInRange(trialParameters$type2Error, low=0.0001, 0.2)) {
                                    stop("Improper type 2 error")
                                }
                                if (!scalarInRange(trialParameters$eps, low=1e-5, high = 1 - 1e-5)) {
                                    stop("Improper epsilon specified")
                                }
                                if (any(designParameters$prevalence <= 0)) {
                                    stop("Improper prevalence specified")
                                }
                                if (!all.equal(dim(designParameters$mean), c(2, designParameters$J))) {
                                    stop("Mean dimension does not match number of groups")
                                }
                                if (!all.equal(dim(designParameters$sd), c(2, designParameters$J))) {
                                    stop("Mean dimension does not match number of groups")
                                }
                                if (!all(designParameters$sd > 0)) {
                                    stop("SDs are not all positive")
                                }
                                if (!is.null(generateData)) {
                                    if (!is.function(generateData)) {
                                        stop("The generateData argument must be a function")
                                    }
                                    argNames <- sort(names(formals(self$generateData)))
                                    if (any(is.na(match(argNames, names(formals(generateData)))))) {
                                        stop("Arguments for data generator are incorrect")
                                    }
                                }
                                TRUE
                            },
                            wilcoxon = function(x, y, theta = 0) {
                                ## R function wilcox.test return statistic=sum(rank of first sample)-m(m+1)/2
                                ## for a first sample of size m
                                ## The standardized Wilcoxon statistics we want is sum(rank of first sample)-m(m+n+1)/2
                                ## =[wilcox.test(x,y)-mn/2]/sqrt(mn(m+n+1)/12)
                                nx <- length(x)
                                ny <- length(y)
                                (wilcox.test(x, y, exact = FALSE)$statistic - nx * ny * (1/2 + theta)) /
                                    sqrt(nx * ny * (nx + ny + 1) / 12)
                            },
                            selectSubpg = function (data) {
                                data <- data[order(data$subGroup), ]
                                ## ASSUME EACH GROUP is represented!!
                                ## FIX needed if you use names for groups instead of 0, 1, 2, .., J
                                counts <- cumsum(table(data$subGroup))
                                wcx <- sapply(counts[-length(counts)],
                                              function(n) {
                                                  d <- data[seq_len(n), ]
                                                  score <- split(d$score, d$trt)
                                                  private$wilcoxon(score$`1`, score$`0`)
                                              })
                                which.max(wcx)
                            },
                            performInterimLook = function (data, stage) {
                                d <- split(data$score, data$trt)
                                wcx <- private$wilcoxon(d$`1`, d$`0`)
                                bdy <- private$boundaries
                                if (stage < 3) {
                                    if (wcx >= bdy["b"]) { ## Reject
                                        decision <- 1
                                    } else {
                                        effectSize <- private$trialParameters$effectSize
                                        wcx.fut <- private$wilcoxon(d$`1`, d$`0`, theta = effectSize)
                                        if (wcx.fut < bdy["btilde"]) { ## Futility, so accept
                                            decision <- -1
                                        } else {
                                            decision <- 0 ## continue
                                        }
                                    }
                                } else {
                                    if (wcx >= bdy["c"]) { ## Final boundary
                                        decision <- 1 ## reject
                                    } else {
                                        decision <- -1 ## accept
                                    }
                                }
                                list(decision = decision, wcx = wcx)
                            },
                            den.vs = function(v, i, mu.prime, Sigma.prime, fut ) {
                                ## Density function used in integration
                                mu.prime <- mu.prime * v
                                pmvnorm(upper = c(rep(v, nrow(mu.prime) - 1), fut), mean = mu.prime[, i],
                                        sigma = Sigma.prime[[i]], algorithm = Miwa()) * dnorm(v)
                            },
                            mHP.btilde = function (beta = 0.1, cov.J) {
                                ## Function for computing futility boundary btilde
                                numLooks <- length(private$trialParameters$N)
                                crossingProb <- function(btilde = -3) {
                                    if (numLooks == 1) {
                                        pnorm(btilde , mean = 0) - beta
                                    } else {
                                        1 - pmvnorm(lower = rep(btilde, numLooks - 1),
                                                    upper = rep(Inf, numLooks - 1),
                                                    mean = rep(0, numLooks - 1),
                                                    sigma = cov.J[-numLooks, -numLooks],
                                                    algorithm = Miwa()) - beta
                                    }
                                }
                                btilde <- uniroot(f = crossingProb, lower = qnorm(beta) - 1,
                                                  upper = qnorm(beta^(1 / (numLooks - 1))) + 1)

                                btilde$root
                            },
                            mHP.b = function (cov.J, mu.prime, Sigma.prime, alpha, btilde) {
                                ## Function for computing efficary boundary b
                                trialParameters <- private$trialParameters
                                designParameters <- private$designParameters
                                J <- designParameters$J
                                N <- trialParameters$N
                                numLooks <- length(N)
                                effectSize <- trialParameters$effectSize
                                q <- cumsum(designParameters$prevalence)
                                crossingProb <- function(b) {
                                    f <- function( stage, time.accept.J, i) {
                                        ##i=sub-population selected
                                        ##time.accept.J = time when H_J is accepted (leading to subgroup selection),
                                        ## stage = time when H_J is rejected
                                        btilde <- btilde + effectSize * sqrt(3 * N[time.accept.J])
                                        ssi <- replace(N, time.accept.J, N[time.accept.J] * q[i])
                                        if (stage == time.accept.J) {
                                            return(integrate(function(x) {
                                                sapply(x,
                                                       function(x) private$den.vs(x, i, mu.prime,
                                                                                  Sigma.prime, btilde))}, b,
                                                Inf)$value)
                                        } else if (stage == (time.accept.J + 1) ) {
                                            sigma <- sqrt(ssi[stage - 1] / ssi[stage])
                                            integrand <- function(u) {
                                                private$den.vs(u, i, mu.prime, Sigma.prime, btilde) *
                                                    pnorm(b, mean = u * sigma, sd = sqrt(1 - sigma^2),
                                                          lower.tail = FALSE)
                                            }
                                            return(integrate(function(u) { sapply(u, integrand) },
                                                             -Inf, b)$value)
                                        }
                                    }
                                    ##type I error at interims=P(accept H_J at stage 1, reject H_I at stage 1)+
                                    ##                         P(accept H_J at stage 1, reject H_I at stage 2)+
                                    ##                         P(accept H_J at stage 2, reject H_I at stage 2)+
                                    ##                         P(reject H_J at stage 1 or 2)
                                    sum(sapply(seq_len(J - 1),
                                               function(i)  f(1, 1, i) + f(2, 1, i) + f(2, 2, i))) +
                                        1 - pmvnorm(lower = -Inf, upper = b, mean = rep(0, numLooks - 1),
                                                    sigma = cov.J[-numLooks, -numLooks], algorithm = Miwa()) - alpha
                                }
                                uniroot(f = crossingProb, lower = 1, upper = 4, maxiter = 20)$root
                            },
                            mHP.c = function (cov.J, mu.prime, Sigma.prime, alpha, btilde, b) {
                                ## Function for computing final boundary c
                                trialParameters <- private$trialParameters
                                designParameters <- private$designParameters
                                J <- designParameters$J
                                N <- trialParameters$N
                                numLooks <- length(N)
                                q <- cumsum(trialParameters$prevalence)
                                effectSize <- trialParameters$effectSize
                                q <- cumsum(designParameters$prevalence)
                                crossingProb <- function(c) {
                                    f <- function(time.accept.J, i) {#i=sub-populatino selected
                                        btilde <- btilde + effectSize * sqrt(3 * N[time.accept.J])
                                        ssi <- replace(N, time.accept.J, N[time.accept.J] * q[i])
                                        if (time.accept.J == 3 ) {
                                            return(integrate(function(x) { sapply(x, function(x)
                                                private$den.vs(x, i, mu.prime, Sigma.prime, btilde))}, c, Inf)$value)
                                        } else if (time.accept.J == 2 ) {
                                            sigma <- sqrt(ssi[2] / ssi[3])
                                            integrand <- function(u) {
                                                private$den.vs(u, i, mu.prime, Sigma.prime, btilde) *
                                                    pnorm(c, mean = u * sigma, sd = sqrt(1 - sigma^2), lower.tail = FALSE)
                                            }
                                            return(integrate(function(u) {sapply(u, integrand)}, -Inf, b)$value)
                                        } else {
                                            v23 <- c(sqrt(ssi[1] / ssi[2]), sqrt(ssi[1] / ssi[3]))
                                            sigma23 <- matrix(sqrt(ssi[2] / ssi[3]), 2, 2)
                                            diag(sigma23) <- 1

                                            sigma <- sigma23 - v23 %*% t(v23)
                                            integrand <- function(u) {
                                                private$den.vs(u, i, mu.prime, Sigma.prime, btilde) *
                                                    pmvnorm(lower = c(-Inf, c), upper = c(b, Inf),
                                                            mean = u * v23, sigma = sigma)
                                            }
                                            return(integrate(function(u) { sapply(u, integrand) }, -Inf, b)$value)
                                        }
                                    }
                                    ##type I error at final=P(accept H_J at stage 1, reject H_I at stage 3)+
                                    ##                      P(accept H_J at stage 2, reject H_I at stage 3)+
                                    ##                      P(accept H_J at stage 3, reject H_I at stage 3)+
                                    ##                      P(reject H_J at stage 3)
                                    sum(sapply(seq_len(J-1), function(i) f(1, i) + f(2, i) + f(3, i))) +
                                        pmvnorm(lower = c(rep(-Inf, numLooks - 1), c),
                                                upper = c(rep(b, numLooks - 1), Inf), sigma = cov.J,
                                                mean = rep(0, numLooks), algorithm = Miwa()) - alpha
                                }
                                uniroot(f = crossingProb, lower= b - 1, upper = b + 1, maxiter = 20)$root
                            }
                        ),
                        public = list(
                            generateData = function(prevalence = rep(1/6, 6), N,
                                                    mean = matrix(0, 2, 6),
                                                    sd = matrix(1, 2, 6)) {
                                if (N == 0) {
                                    data.frame(subGroup = integer(0), trt = integer(0),
                                               score = numeric(0))
                                } else {
                                    subGroup <- sample(seq_along(prevalence), N, replace = TRUE,
                                                       prob = prevalence)
                                    trt <- sample(c(0L, 1L), N, replace = TRUE)
                                    rankin <- unlist(
                                        Map(function(i, j)
                                            rnorm(n=1, mean = mean[i, j], sd = sd[i, j]),
                                            trt + 1, subGroup))
                                    data.frame(subGroup = subGroup, trt = trt, score = rankin)
                                }
                            },
                            initialize = function(designParameters, trialParameters, generateData = NULL) {
                                designParameters$J <- length(designParameters$prevalence)
                                private$checkParameters(designParameters, trialParameters, generateData)
                                trialParameters$effectSize <- (qnorm(1 - trialParameters$type1Error) +
                                                               qnorm(1 - trialParameters$type2Error)) /
                                    sqrt(3 * trialParameters$N[3])
                                designParameters$prevalence <- designParameters$prevalence / sum(designParameters$prevalence)
                                private$designParameters <- designParameters
                                private$trialParameters <- trialParameters
                                if (!is.null(generateData)) {
                                    self$generateData <- generateData
                                }
                                private$boundaries <- self$computeCriticalValues()
                            },
                            getDesignParameters = function() private$designParameters,
                            getTrialParameters = function() private$trialParameters,
                            getBoundaries  = function() private$boundaries,
                            print = function() {
                                cat("Design Parameters:\n")
                                str(private$designParameters)
                                cat("Trial Parameters:\n")
                                str(private$trialParameters)
                                cat("Boundaries:\n")
                                str(private$boundaries)
                                cat("Data Generating function:\n")
                                print(self$generateData)
                            },
                            computeCriticalValues = function() {
                                ## Find eff boundary for subgp
                                ## Sigma = covariance matrix between subgroup,
                                ## which is roughly stage independent
                                trialParameters <- private$trialParameters
                                designParameters <- private$designParameters
                                J <- designParameters$J
                                N <- trialParameters$N
                                alpha <- trialParameters$type1Error
                                beta <- trialParameters$type2Error
                                eps <- trialParameters$eps
                                effectSize <- trialParameters$effectSize
                                q <- cumsum(designParameters$prevalence)
                                Sigma <- matrix(0, J, J)
                                for (i in seq_len(J - 1)) {
                                    for (j in (i + 1):J) {
                                        Sigma[i, j] <- sqrt(q[i] / q[j])
                                    }
                                }
                                Sigma <- Sigma + t(Sigma)
                                diag(Sigma) <- 1
                                mu.prime <- matrix(0, (J - 1), J)
                                Sigma.prime <- vector("list", J)
                                for (i in seq_len(J)) {
                                    mu.prime[, i] <- Sigma[-i, i]
                                    Sigma.prime[[i]] <- Sigma[-i, -i] - Sigma[-i, i] %*% t(Sigma[i, -i])
                                }
                                numLooks <- length(private$trialParameters$N)
                                cov.J <- matrix(0, numLooks, numLooks)
                                for (i in seq_len(numLooks - 1)) {
                                    for (j in (i + 1):numLooks) {
                                        cov.J[i, j] <- sqrt(N[i]^2 * (N[j] + 1) / N[j]^2 / (N[i] + 1) )
                                    }
                                }
                                cov.J <- cov.J + t(cov.J)
                                diag(cov.J) <- 1
                                btilde <- private$mHP.btilde(beta * eps, cov.J )
                                b <- private$mHP.b(cov.J, mu.prime, Sigma.prime, alpha * eps, btilde)
                                ##b.I=2.379879
                                c <- private$mHP.c(cov.J, mu.prime, Sigma.prime, alpha * (1 - eps), btilde, b)
                                c(btilde = btilde, b = b, c = c)
                            },
                            explore = function (numberOfSimulations = 5000, rngSeed = 12345,
                                                trueParameters = self$getDesignParameters(),
                                                showProgress = TRUE) {
                                ## Save rng state
                                oldRngState <- if (exists(".Random.seed", envir = .GlobalEnv)) {
                                    get(x = ".Random.seed", envir=.GlobalEnv)
                                } else {
                                    NULL
                                }
                                ## set our seed
                                set.seed(seed = rngSeed, normal.kind = NULL)

                                naVec <- rep(NA, numberOfSimulations)
                                zeroVec <- integer(numberOfSimulations)
                                trialHistory <- data.frame(stage = naVec, decision.ITT = naVec,
                                                           decision.subgp = naVec, select = naVec,
                                                           statistic = naVec, lost = zeroVec,
                                                           ITTfutStage = zeroVec,
                                                           ITTfutSubgp = zeroVec)

                                if (showProgress) {
                                    pb <- txtProgressBar(min = 0, max = numberOfSimulations, style = 3)
                                }
                                trialParameters <- private$trialParameters
                                if (is.null(trueParameters$J)) {
                                    trueParameters$J <- length(trueParameters$prevalence)
                                }

                                glrBoundary <- private$boundaries
                                prevalence <- trueParameters$prevalence

                                for (i in seq_len(numberOfSimulations)) {
                                    dataSoFar <- data.frame(subGroup = integer(0), trt = integer(0),
                                                            score = numeric(0))
                                    subgp <- trueParameters$J
                                    previousN <- 0
                                    j <- 1
                                    interim <- NULL
                                    while (j <= 3) {
                                        dataSoFar <- rbind(dataSoFar[dataSoFar$subGroup <= subgp, ],
                                                           self$generateData(prevalence = prevalence[1:subgp],
                                                                             N = trialParameters$N[j] - previousN,
                                                                             mean = trueParameters$mean[, 1:subgp,
                                                                                                             drop = FALSE],
                                                                             sd = trueParameters$sd[, 1:subgp,
                                                                                                         drop = FALSE]))
                                        interim <- private$performInterimLook(dataSoFar, stage = j)
                                        if (interim$decision == 1) {
                                            break #interim$decision=1 if reject, -1 if accept
                                        } else if (interim$decision == 0) { #interim$decision=0 if continue
                                            previousN <- nrow(dataSoFar)
                                            j <- j + 1
                                        } else { #interim$decision=-1 if accept
                                            if (subgp == trueParameters$J) {
                                                subgp <- private$selectSubpg(dataSoFar)
                                                previousN <- nrow(dataSoFar)
                                                trialHistory$ITTfutStage[i] <- j
                                                trialHistory$ITTfutSubgp[i] <- subgp
                                                trialHistory$lost[i] <- trialParameters$N[j] - sum(dataSoFar$subGroup <= subgp)

                                            } else {
                                                break
                                            }
                                        }
                                    }
                                    if (subgp == trueParameters$J) {
                                        trialHistory[i, 1:5] <- c(stage = j, decision.ITT = interim$decision,
                                                                  decision.subgp = NA, select = subgp,
                                                                  statistic = interim$wcx)
                                    } else {
                                        trialHistory[i, 1:5] <- c(stage = j, decision.ITT = -1,
                                                                  decision.subgp = interim$decision,
                                                                  select = subgp, statistic = interim$wcx)
                                    }
                                    if (showProgress) {
                                        setTxtProgressBar(pb, i)
                                    }
                                }
                                if (showProgress) {
                                    close(pb)
                                }
                                ## Restore rng state
                                if (is.null(oldRngState)) {
                                    rm(".Random.seed", envir = .GlobalEnv)
                                } else {
                                    assign(x = ".Random.seed", value = oldRngState, envir = .GlobalEnv)
                                }
                                trialHistory
                            },
                            analyze = function (trialHistory) {
                                numberOfSimulations <- nrow(trialHistory)
                                reject.ITT <- (trialHistory$decision.ITT == 1)
                                reject.subgp <- (trialHistory$decision.subgp == 1)
                                reject.subgp[is.na(reject.subgp)] <- FALSE
                                reject <- !(reject.ITT + reject.subgp == 0)
                                ##  numReject <- sum(Reject)
                                earlyStop <- (trialHistory$stage < 3)
                                earlyStopEff <- (reject & earlyStop)
                                earlyStopFut <- (!reject & earlyStop)
                                popReject <- table(trialHistory$select[reject]) / numberOfSimulations
                                exitRandSS <- private$trialParameters$N[trialHistory$stage]
                                exitAnalyzeSS <- exitRandSS - trialHistory$lost
                                stageAtExitProportion <- table(trialHistory$stage) / numberOfSimulations

                                futilityTable <- table(trialHistory$ITTfutStage, trialHistory$ITTfutSubgp)
                                meanLossFutility <- tapply(trialHistory$lost,
                                                           list(trialHistory$ITTfutStage,trialHistory$ITTfutSubgp),
                                                           mean)
                                sdLossFutility <- tapply(trialHistory$lost,
                                                         list(trialHistory$ITTfutStage,trialHistory$ITTfutSubgp),
                                                         sd)

                                list(numberOfSimulations = numberOfSimulations,
                                     reject.ITT = reject.ITT,
                                     reject.subgp = reject.subgp,
                                     reject = reject,
                                     earlyStopEff = earlyStopEff,
                                     earlyStopFut = earlyStopFut,
                                     popReject = popReject,
                                     exitRandSS = exitRandSS,
                                     exitAnalyzeSS = exitAnalyzeSS,
                                     lost = trialHistory$lost,
                                     stageAtExitProportion = stageAtExitProportion,
                                     futilityTable = futilityTable,
                                     meanLossFutility = meanLossFutility,
                                     sdLossFutility = sdLossFutility
                                     )
                            },
                            summary = function(analysis) {
                                cat(sprintf("P(Reject H0_ITT) = %f; P(Reject H0_subgp) = %f; P(Reject H0) = %f\n",
                                            mean(analysis$reject.ITT), mean(analysis$reject.subgp),
                                            mean(analysis$reject)))
                                cat(sprintf("P(Early stop for efficacy [futility]) = %f [%f]\n",
                                            mean(analysis$earlyStopEff), mean(analysis$earlyStopFut)))
                                cat(sprintf("Mean [SD] Randomized N = %f [%f]\n",
                                            mean(analysis$exitRandSS), sd(analysis$exitRandSS)))
                                cat("\nStage at exit (proportion)\n")
                                print(analysis$stageAtExitProportion)
                                cat(sprintf("\nMean [SD] Lost N = %f [%f]\n",
                                            mean(analysis$lost), sd(analysis$lost)))
                                cat(sprintf("Mean [SD] Analyzed N = %f [%f]\n",
                                            mean(analysis$exitAnalyzeSS), sd(analysis$exitAnalyzeSS)))
                                cat("\nChance of each subpopulation rejected\n")
                                print(analysis$popReject)
                                cat("\nCounts by futility stage and subgroup choice\n")
                                print(analysis$futilityTable)
                                cat("\nMean loss by futility stage and subgroup\n")
                                print(analysis$meanLossFutility)
                                cat("\nSD loss by futility stage and subgroup\n")
                                print(analysis$sdLossFutility)
                            }
                        ))


#' A fixed sample design to compare against the adaptive clinical trial design of Lai, Lavori and Liao.
#'
#'
#' @description \code{ASSISTDesignB} objects are used to design a trial with certain
#' characteristics provided in the object instantiation method. This design differs from
#' \code{ASSISTDesign} in only how it computes the critical boundaries, how it performs the interim look,
#' and what quantities are computed in a trial run.
#'
#' @docType class
#' @seealso \code{ASSISTDesign} which is a superclass of this object
#' @importFrom R6 R6Class
#' @importFrom mvtnorm pmvnorm Miwa
#' @importFrom stats uniroot rnorm pnorm qnorm
#' @section Methods:
#'
#' \describe{
#'   \item{\code{ASSISTDesignB$new(designParameters, trialParameters, generateData)}}{Create a new \code{ASSISTDesign}
#'         instance object using the parameters specified. }
#'   \item{\code{getDesignParameters},\code{getTrialParameters},
#'         \code{getBoundaries}}{Accessor methods for (obvious) object slots}
#'   \item{\code{print()}}{Print the object in a human readable form}
#'   \item{\code{computeCriticalValues()}}{Compute the critical boundary value \eqn{c_\alpha}}
#'   \item{\code{explore(numberOfSimulations = 5000, rngSeed = 12345, trueParameters = self$getDesignParameters(), showProgress = TRUE)}}{Explore the design
#'         using the specified number of simulations and random number seed.  \code{trueParameters} is by default the same
#'         as \code{designParameters} as would be the case for a Type I error calculation. If changed, would yield power.
#'         Show progress if so desired. Returns a data frame of results}
#'   \item{\code{analyze(trialHistory)}}{Analyze
#'         the design given the \code{trialHistory} which is the result of a call to \code{explore} to
#'         simulate the design. Return a list of summary quantities}
#'   \item{\code{summary(analysis)}}{Print the operating characteristics of the design, using the analysis
#'         result from the \code{analyze} call}
#' }
#'
#' @references Adaptive Choice of Patient Subgroup for Comparing Two Treatments
#' by Tze Leung Lai and Philip W. Lavori and Olivia Yueh-Wen Liao. Contemporary Clinical Trials,
#' Vol. 39, No. 2, pp 191-200 (2014). \url{http://www.sciencedirect.com/science/article/pii/S1551714414001311}
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @examples
#' \dontrun{
#' data(LLL.SETTINGS)
#' prevalence <- LLL.SETTINGS$prevalences$table1
#' scenario <- LLL.SETTINGS$scenarios$S0
#' designParameters <- list(prevalence = prevalence,
#'                        mean = scenario$mean,
#'                        sd = scenario$sd)
#' designB <- ASSISTDesignB$new(trialParameters = LLL.SETTINGS$trialParameters,
#'                             designParameters = designParameters)
#' print(designB)
#' ## A realistic design uses 5000 simulations or more!
#' result <- designB$explore(showProgress = interactive())
#' analysis <- designB$analyze(result)
#' designB$summary(analysis)
#' }
#' ## For full examples, try:
#' ## browseURL(system.file("full_doc/ASSISTant.html", package="ASSISTant"))
#'
ASSISTDesignB <- R6Class("ASSISTDesignB",
                         inherit = ASSISTDesign,
                         private = list(
                             performInterimLook = function (data) {
                                 d <- split(data$score, data$trt)
                                 wcx <- private$wilcoxon(d$`1`, d$`0`)
                                 bdy <- private$boundaries
                                 if (wcx >= bdy["cAlpha"]) { ## Final boundary
                                     decision <- 1 ## reject
                                 } else {
                                     decision <- -1 ## accept
                                 }

                                 list(decision = decision, wcx = wcx)
                             },
                             ## Function for computing futility boundary btilde
                             mHP.ITT = function (mu.prime, Sigma.prime, alpha) {
                                 J <- private$designParameters$J
                                 ## Derive interim eff boundary b.I for subgp
                                 crossingProb <- function(c) {
                                     f <- function(i) { #i=sub-population selected
                                         integrate(function(x) {
                                             sapply(x, function(x)
                                                 private$den.vs(x, i, mu.prime, Sigma.prime, c))}, c, Inf)$value
                                     }
                                     sum(sapply(seq_len(J - 1), function(i) f(i))) +
                                         pnorm(c, lower.tail = FALSE) - alpha
                                 }
                                 uniroot(f = crossingProb, lower = 1, upper = 4, maxiter = 20)$root
                             }
                             ),
                         public = list(
                             computeCriticalValues = function() {
                                 ## Find eff boundary for subgp
                                 ## Sigma = covariance matrix between subgroup,
                                 ## which is roughly stage independent
                                 trialParameters <- private$trialParameters
                                 designParameters <- private$designParameters
                                 J <- designParameters$J
                                 alpha <- trialParameters$type1Error
                                 q <- cumsum(private$designParameters$prevalence)
                                 Sigma <- matrix(0, J, J)
                                 for (i in seq_len(J - 1)) {
                                     for (j in (i + 1):J) {
                                         Sigma[i, j] <- sqrt(q[i] / q[j])
                                     }
                                 }
                                 Sigma <- Sigma + t(Sigma)
                                 diag(Sigma) <- 1
                                 mu.prime <- matrix(0, (J - 1), J)
                                 Sigma.prime <- vector("list", J)
                                 for (i in seq_len(J)) {
                                     mu.prime[, i] <- Sigma[-i, i]
                                     Sigma.prime[[i]] <- Sigma[-i, -i] - Sigma[-i, i] %*% t(Sigma[i, -i])
                                 }
                                 cAlpha <- private$mHP.ITT(mu.prime, Sigma.prime, alpha)
                                 list(cAlpha = cAlpha)
                             },
                             explore = function (numberOfSimulations = 100, rngSeed = 12345,
                                                 trueParameters = self$getDesignParameters(),
                                                 showProgress = TRUE) {
                                 ## Save rng state
                                 oldRngState <- if (exists(".Random.seed", envir = .GlobalEnv)) {
                                     get(x = ".Random.seed", envir=.GlobalEnv)
                                 } else {
                                     NULL
                                 }
                                 ## set our seed
                                 set.seed(seed = rngSeed, normal.kind = NULL)

                                 trialParameters <- private$trialParameters

                                 if (is.null(trueParameters$J)) {
                                     trueParameters$J <- length(trueParameters$prevalence)
                                 }
                                 J <- trueParameters$J

                                 glrBoundary <- private$boundaries

                                 naVec <- rep(NA, numberOfSimulations)
                                 zeroVec <- integer(numberOfSimulations)
                                 trialHistory <- data.frame(decision = naVec, select = naVec,
                                                            statistic = naVec,
                                                            matrix(0, numberOfSimulations, J))

                                 if (showProgress) {
                                     pb <- txtProgressBar(min = 0, max = numberOfSimulations, style = 3)
                                 }

                                 for (i in seq_len(numberOfSimulations)) {
                                     dataSoFar <- self$generateData(prevalence = trueParameters$prevalence,
                                                                    N = trialParameters$N[3],
                                                                    mean = trueParameters$mean,
                                                                    sd = trueParameters$sd)
                                     interim <- private$performInterimLook(dataSoFar)
                                     subgp <- J ## Last group
                                     if (interim$decision == -1) { ## continue
                                         subgp <- private$selectSubpg(dataSoFar)
                                         interim <- private$performInterimLook(dataSoFar[dataSoFar$subGroup <= subgp, ])
                                     }
                                     trialHistory[i, ] <- c(decision = interim$decision,
                                                            select = subgp,
                                                            statistic = interim$wcx,
                                                            table(dataSoFar$subGroup))
                                     if (showProgress) {
                                         setTxtProgressBar(pb, i)
                                     }
                                 }
                                 if (showProgress) {
                                     close(pb)
                                 }
                                 ## Restore rng state
                                 if (is.null(oldRngState)) {
                                     rm(".Random.seed", envir = .GlobalEnv)
                                 } else {
                                     assign(x = ".Random.seed", value = oldRngState, envir = .GlobalEnv)
                                 }
                                 names(trialHistory) <- c("decision", "select", "statistic",
                                                          sapply(seq_len(J), function(i) paste0("G", i)))
                                 trialHistory
                             },
                             analyze = function (trialHistory) {
                                 J <- private$designParameters$J
                                 numberOfSimulations <- nrow(trialHistory)
                                 reject <- (trialHistory$decision == 1)
                                 rejectGroupTable <- table(trialHistory$select[reject])

                                 list(reject = reject, rejectGroupTable = rejectGroupTable,
                                      rejectSubgroup = sum(rejectGroupTable[-J]) / numberOfSimulations)
                             },
                             summary = function(analysis) {
                                 numberOfSimulations <- length(analysis$reject)
                                 cat(sprintf("P(Reject H0) = %f\n",
                                             mean(analysis$reject)))
                                 cat(sprintf("P(Reject H0_ITT) = %f\n",
                                             mean(analysis$reject) - analysis$rejectSubgroup))
                                 cat(sprintf("P(Reject H0_subgp) = %f\n",
                                             analysis$rejectSubgroup))
                                 cat("\nChance of each subpopulation rejected\n")
                                 print(analysis$rejectGroupTable / numberOfSimulations)
                             }
                         ))

#' A fixed sample RCT design to compare against the adaptive clinical trial design of Lai, Lavori and Liao.
#'
#'
#' @description \code{ASSISTDesignC} objects are used to design a trial with certain
#' characteristics provided in the object instantiation method. This design differs from
#' \code{ASSISTDesign} in only how it computes the critical boundaries, how it performs the interim look,
#' and what quantities are computed in a trial run.
#'
#' @docType class
#' @seealso \code{ASSISTDesignB} which is a superclass of this object
#' @importFrom R6 R6Class
#' @importFrom mvtnorm pmvnorm Miwa
#' @importFrom stats uniroot rnorm pnorm qnorm
#' @section Methods:
#'
#' \describe{
#'   \item{\code{ASSISTDesignC$new(designParameters, trialParameters, generateData)}}{Create a new
#'         \code{ASSISTDesign} instance object using the parameters specified. }
#'   \item{\code{getDesignameters},\code{getTrialParameters},
#'         \code{getBoundaries}}{Accessor methods for (obvious) object slots}
#'   \item{\code{print()}}{Print the object in a human readable form}
#'   \item{\code{computeCriticalValues()}}{Compute the critical boundary value \eqn{c_\alpha}}
#'   \item{\code{explore(numberOfSimulations = 5000, rngSeed = 12345, trueParameters = self$getDesignParameters(), showProgress = TRUE)}}{Explore the design
#'         using the specified number of simulations and random number seed.  \code{trueParameters} is by default the same
#'         as \code{designParameters} as would be the case for a Type I error calculation. If changed, would yield power.
#'         Show progress if so desired. Returns a data frame of results}
#'   \item{\code{analyze(trialHistory)}}{Analyze
#'         the design given the \code{trialHistory} which is the result of a call to \code{explore} to
#'         simulate the design. Return a list of summary quantities}
#' \item{\code{summary(analysis)}}{Print the operating characteristics of the design, using the analysis
#'         result from the \code{analyze} call}
#' }
#'
#' @references Adaptive Choice of Patient Subgroup for Comparing Two Treatments
#' by Tze Leung Lai and Philip W. Lavori and Olivia Yueh-Wen Liao. Contemporary Clinical Trials,
#' Vol. 39, No. 2, pp 191-200 (2014). \url{http://www.sciencedirect.com/science/article/pii/S1551714414001311}
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @examples
#' data(LLL.SETTINGS)
#' prevalence <- LLL.SETTINGS$prevalences$table1
#' scenario <- LLL.SETTINGS$scenarios$S0
#' designParameters <- list(prevalence = prevalence,
#'                        mean = scenario$mean,
#'                        sd = scenario$sd)
#' ## A realistic design uses 5000 simulations or more!
#' designC <- ASSISTDesignC$new(trialParameters = LLL.SETTINGS$trialParameters,
#'                             designParameters = designParameters)
#' print(designC)
#' result <- designC$explore(numberOfSimulations = 100, showProgress = interactive())
#' analysis <- designC$analyze(result)
#' designC$summary(analysis)
#' ## For full examples, try:
#' ## browseURL(system.file("full_doc/ASSISTant.html", package="ASSISTant"))
#'
ASSISTDesignC <- R6Class("ASSISTDesignC",
                         inherit = ASSISTDesignB,
                         public = list(
                             computeCriticalValues = function() {
                                 ## Find eff boundary for subgp
                                 ## Sigma = covariance matrix between subgroup,
                                 ## which is roughly stage independent
                                 trialParameters <- private$trialParameters
                                 alpha <- trialParameters$type1Error
                                 list(cAlpha = qnorm(1 - alpha))
                             },
                             explore = function (numberOfSimulations = 5000, rngSeed = 12345,
                                                 trueParameters = self$getDesignParameters(),
                                                 showProgress = TRUE) {
                                 ## Save rng state
                                 oldRngState <- if (exists(".Random.seed", envir = .GlobalEnv)) {
                                     get(x = ".Random.seed", envir=.GlobalEnv)
                                 } else {
                                     NULL
                                 }
                                 ## set our seed
                                 set.seed(seed = rngSeed, normal.kind = NULL)
                                 trialParameters <- private$trialParameters

                                 if (is.null(trueParameters$J)) {
                                     trueParameters$J <- length(trueParameters$prevalence)
                                 }
                                 J <- trueParameters$J
                                 glrBoundary <- private$boundaries

                                 naVec <- rep(NA, numberOfSimulations)
                                 zeroVec <- integer(numberOfSimulations)
                                 trialHistory <- data.frame(decision = naVec,
                                                            statistic = naVec,
                                                            matrix(0, numberOfSimulations, J))

                                 if (showProgress) {
                                     pb <- txtProgressBar(min = 0, max = numberOfSimulations, style = 3)
                                 }

                                 for (i in seq_len(numberOfSimulations)) {
                                     dataSoFar <- self$generateData(prevalence = trueParameters$prevalence,
                                                                    N = trialParameters$N[3],
                                                                    mean = trueParameters$mean,
                                                                    sd = trueParameters$sd)
                                     interim <- private$performInterimLook(dataSoFar)
                                     trialHistory[i, ] <- c(decision = interim$decision,
                                                            statistic = interim$wcx,
                                                            table(dataSoFar$subGroup))
                                     if (showProgress) {
                                         setTxtProgressBar(pb, i)
                                     }
                                 }
                                 if (showProgress) {
                                     close(pb)
                                 }
                                 ## Restore rng state
                                 if (is.null(oldRngState)) {
                                     rm(".Random.seed", envir = .GlobalEnv)
                                 } else {
                                     assign(x = ".Random.seed", value = oldRngState, envir = .GlobalEnv)
                                 }
                                 names(trialHistory) <- c("decision", "statistic",
                                                          sapply(seq_len(J), function(i) paste0("G", i)))
                                 trialHistory
                             },
                             analyze = function (trialHistory) {
                                 J <- private$designParameters$J
                                 numberOfSimulations <- nrow(trialHistory)
                                 reject <- (trialHistory$decision == 1)
                                 list(reject = reject)
                             },
                             summary = function(analysis) {
                                 numberOfSimulations <- length(analysis$reject)
                                 cat(sprintf("P(Reject H0) = %f\n",
                                             mean(analysis$reject)))
                             }
                         ))


#' The DEFUSE3 design
#'
#'
#' @description \code{DEFUSE3Design} is a slight variant of the the adaptive
#' clinical trial design of Lai, Lavori and Liao. Simulation is used to compute
#' the expected maximum sample size and the boundary for early futility is adjusted to
#' account as well.
#'
#' @docType class
#' @seealso \code{ASSISTDesign} which is a superclass of this object
#' @importFrom R6 R6Class
#' @importFrom mvtnorm pmvnorm Miwa
#' @importFrom stats uniroot rnorm pnorm qnorm
#' @section Methods:
#'
#' \describe{
#'   \item{\code{DEFUSE3Design$new(designParameters, trialParameters, generateData, numberOfSimulations = 5000, rngSeed = 54321, showProgress = TRUE)}}{Create
#'         a new \code{ASSISTDesign} instance object using the parameters specified. }
#'   \item{\code{getDesignParameters},\code{getTrialParameters},
#'         \code{getBoundaries}}{Accessor methods for (obvious) object slots}
#'   \item{\code{print()}}{Print the object in a human readable form}
#'   \item{\code{adjustCriticalValues(numberOfSimulations, rngSeed, showProgress)}}{Adjust the critical values
#'         by performing simulations using the parameters provided}
#'   \item{\code{computeCriticalValues()}}{Compute the critical boundary value \eqn{c_\alpha}}
#'   \item{\code{explore(numberOfSimulations = 5000, rngSeed = 12345, trueParameters = self$getDesignParameters(). showProgress = TRUE)}}{Explore the design
#'         using the specified number of simulations and random number seed.  \code{trueParameters} is by default the same
#'         as \code{designParameters} as would be the case for a Type I error calculation. If changed, would yield power.
#'         Show progress if so desired. Returns a data frame of results}
#'   \item{\code{analyze(trialHistory)}}{Analyze
#'         the design given the \code{trialHistory} which is the result of a call to \code{explore} to
#'         simulate the design. Return a list of summary quantities}
#'   \item{\code{summary(analysis)}}{Print the operating characteristics of the design, using the analysis
#'         result from the \code{analyze} call}
#' }
#'
#' @references Adaptive design of confirmatory trials: Advances and challenges,
#' \url{http://www.sciencedirect.com/science/article/pii/S1551714415300239} by
#' Tze Leung Lai and Philip W. Lavori and Ka Wai Tsang. Contemporary Clinical Trials, Vol. 45, Part A,
#' pp 93-102 (2015).
#'
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @examples
#' trialParameters <- list(N = c(200, 340, 476), type1Error = 0.025,
#'                         eps = 1/2, type2Error = 0.1)
#' designParameters <- list(
#'    nul0 = list(prevalence = rep(1/6, 6), mean = matrix(0, 2, 6),
#'                sd = matrix(1, 2, 6)),
#'    alt1 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
#'                c(0.5, 0.4, 0.3, 0, 0, 0)),
#'                sd = matrix(1, 2, 6)),
#'    alt2 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
#'                c(0.5, 0.5, 0, 0, 0, 0)),
#'                sd = matrix(1,2, 6)),
#'    alt3 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(0.36, 6)),
#'                sd = matrix(1,2, 6)),
#'    alt4 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6), rep(0.30, 6)),
#'                sd = matrix(1,2, 6)),
#'    alt5 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
#'                c(0.4, 0.3, 0.2, 0, 0, 0)),
#'                sd = matrix(1,2, 6)),
#'    alt6 = list(prevalence = rep(1/6, 6), mean = rbind(rep(0, 6),
#'                c(0.5, 0.5, 0.3, 0.3, 0.1, 0.1)),
#'                sd = matrix(1,2, 6)))
#'
#'\dontrun{
#' ## A realistic design uses 5000 simulations or more!
#' defuse3 <- DEFUSE3Design$new(trialParameters = trialParameters,
#'                              numberOfSimulations = 25,
#'                              designParameters = designParameters$nul0,
#'                              showProgress = FALSE)
#' print(defuse3)
#' result <- defuse3$explore(showProgress = interactive())
#' analysis <- defuse3$analyze(result)
#' print(defuse3$summary(analysis))
#' }
#' ## For full examples, try:
#' ## browseURL(system.file("full_doc/defuse3.html", package="ASSISTant"))
#'
DEFUSE3Design <- R6Class("DEFUSE3Design",
                         inherit = ASSISTDesign,
                         public = list(
                             initialize = function(designParameters, trialParameters, generateData = NULL,
                                                   numberOfSimulations = 5000, rngSeed = 54321,
                                                   showProgress = TRUE,
                                                   trueParameters = NULL) {
                                 super$initialize(designParameters, trialParameters, generateData)
                                 self$adjustCriticalValues(numberOfSimulations, rngSeed, showProgress)
                             },
                             adjustCriticalValues = function(numberOfSimulations, rngSeed, showProgress) {
                                 designParameters <- private$designParameters
                                 trialParameters <- private$trialParameters
                                 ## Run simulations to estimate expect max sample sizes
                                 result <- self$explore(numberOfSimulations = numberOfSimulations,
                                                        rngSeed = rngSeed,
                                                        showProgress = showProgress)
                                 simDN <- matrix(NA, numberOfSimulations, 3)
                                 q <- cumsum(designParameters$prevalence)
                                 for (i in seq_len(numberOfSimulations)) {
                                     j <- result$stage[i]
                                     simDN[i, seq_len(j)] <- diff(c(0, trialParameters$N[seq_len(j)]))
                                     jfut <- result$ITTfutStage[i]
                                     if (jfut > 0) {
                                         subgp <- result$select[i]
                                         simDN[i, seq_len(jfut)] <- simDN[i, seq_len(jfut)] * q[subgp]
                                     }
                                 }

                                 ## End Simulation
                                 ## Adjust bTilde and effectSize
                                 DEM <- apply(simDN, 2, mean, na.rm = TRUE)
                                 EM <- cumsum(DEM) ## Expected N actually
                                 J <- designParameters$J
                                 beta <- trialParameters$type2Error
                                 eps <- trialParameters$eps
                                 numLooks <- length(trialParameters$N)
                                 cov.J <- matrix(0, numLooks, numLooks)
                                 for (i in seq_len(numLooks - 1)) {
                                     for (j in (i + 1):numLooks) {
                                         cov.J[i, j] <- sqrt(EM[i]^2 * (EM[j] + 1) / EM[j]^2 / (EM[i] + 1) )
                                     }
                                 }
                                 cov.J <- cov.J + t(cov.J)
                                 diag(cov.J) <- 1
                                 ## Recompute btilde
                                 private$boundaries["btilde"] <- private$mHP.btilde(beta * eps, cov.J )
                                 ## Recompute Effect Size
                                 private$trialParameters$effectSize <- (qnorm(1 - trialParameters$type1Error) +
                                                                        qnorm(1 - trialParameters$type2Error)) /
                                     sqrt(3 * EM[3])
                             }
                         ))

