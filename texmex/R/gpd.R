gpd <- texmexFamily(name = 'GPD',
                   log.lik = function(data, th, ...) {
                               y <- data$y
                               X.phi <- data$D$phi
                               X.xi <- data$D$xi
                               n.phi <- ncol(X.phi)
                               n.end <- n.phi + ncol(X.xi)
                               function(param) {
                                 stopifnot(length(param) == n.end)
                                 phi <- X.phi %*% param[1:n.phi]
                                 xi <- X.xi %*% param[(1 + n.phi):n.end]
                                 sum(dgpd(y, exp(phi), xi, u=th, log.d=TRUE))
                                }
                   }, # Close log.lik
                   param = c('phi', 'xi'),
                   info = gpd.info,
                   start = function(data){
                             y <- data$y
                             X.phi <- data$D[[1]]
                             X.xi <- data$D[[2]]
                             c(log(mean(y)), rep(1e-05, -1 + ncol(X.phi) + ncol(X.xi)))
                    }, # Close start

                    resid = function(o){
                              p <- texmexMakeParams(coef(o), o$data$D)
                              delta <- (o$data$y - o$threshold) / exp(p[,1])
                              .log1prel(delta * p[,2]) * delta # Standard exponential
                    }, # Close resid

                    endpoint = function(param, model){
                                 model$threshold - exp(param[, 1]) / param[, 2]
                    },
                    rl = function(m, param, model){
                      ## write in terms of qgpd; let's not reinvent the wheel
                      qgpd(1/(m * model$rate),
                           exp(param[,1]), param[,2], u=model$threshold,
                           lower.tail=FALSE)
                    },
                    delta = function(param, m, model){
                              # This is not exact if a prior (penalty) function is used, but
                              # the CI is approximate anyway.
                             param <- c(model$rate, param)
                             out <- matrix(0, nrow=2, ncol=length(m))
                             if (param[3] == 0){ # exponential case
                               out[1,] <- exp(param[2]) * log(m * param[1])
                             } else {
                             out[1,] <- exp(param[2]) / param[3] * ((m * param[1])^param[3] - 1)
                             out[2,] <- -exp(param[2]) / (param[3] * param[3]) * ( (m * param[1] )^param[3] - 1 ) +
                             exp(param[2]) / param[3] * (m * param[1])^param[3] * log(m * param[1])
                             }
                             out
                    }, # Close delta
                    density = function(n, param, model){
                                dgpd(n, exp(c(param[, 1])), c(param[, 2]), u=model$threshold)
                    },

                    rng = function(n, param, model){
                            rgpd(n, exp(c(param[, 1])), c(param[, 2]), u=model$threshold)
                    },
                    prob = function(n, param, model){
                             pgpd(n, exp(c(param[, 1])), c(param[, 2]), u=model$threshold)
                    },
                    quant = function(n, param, model){
                              qgpd(n, exp(c(param[, 1])), c(param[, 2]), u=model$threshold)
                    }
)


