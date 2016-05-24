gev <- texmexFamily(name = 'GEV',
                    param = c('mu', 'phi', 'xi'),
                    log.lik = function(data, ...) {
                                y <- data$y
                                X.mu <- data$D$mu
                                X.phi <- data$D$phi
                                X.xi <- data$D$xi

                                n.mu <- ncol(X.mu)
                                n.phi <- n.mu + ncol(X.phi)
                                n.end <- n.phi + ncol(X.xi)

                                function(param) {
                                  stopifnot(length(param) == n.end)
                                  mu <- X.mu %*% param[1:n.mu]
                                  phi <- X.phi %*% param[(1 + n.mu):n.phi]
                                  xi <- X.xi %*% param[(1 + n.phi):n.end]
                                  if (any(1 + xi/exp(phi)*(y-mu) <= 0)) -Inf
                                  else sum(dgev(y, mu, exp(phi), xi, log.d=TRUE))
                                }
                    }, # Close log.lik
                    info = NULL, # will mean that numerical approx gets used
                    delta = function(param, m, model){ # model not used but required by a calling function
                              y <- -log(1 - 1/m)
                              out <- rep(1, 3)

                              out[2] <- -exp(param[2])/param[3] * (1 - y^(-param[3])) # Coles p.56
                              out[3] <- exp(param[2]) * param[3]^(-2) * (1 - y^(-param[3])) - # change of variable from sigma to phi gives exp(param[2]) in out[2]
                              exp(param[2])/param[3] * y^(-param[3]) * log(y)
                              out
                    }, # Close delta
                    start = function(data){
                              y <- data$y
                              X.mu <- data$D[[1]]
                              X.phi <- data$D[[2]]
                              X.xi <- data$D[[3]]

                              c(mean(y), rep(0, ncol(X.mu)-1), log(IQR(y)/2),
                              rep(.001, -1 + ncol(X.phi) + ncol(X.xi)))
                    }, # Close start
                    endpoint = function(param, model){
                                 param[, 1] - exp(param[, 2]) / param[, 3]
                    },
                    rng = function(n, param, model){
                            rgev(n, c(param[, 1]), exp(c(param[, 2])), c(param[, 3]))
                    },
                    density = function(x, param, model){
                                dgev(x, c(param[, 1]), exp(c(param[, 2])), c(param[, 3]))
                    },
                    prob = function(x, param, model){
                             pgev(x, c(param[, 1]), exp(c(param[, 2])), c(param[, 3]))
                    },
                    quant = function(p, param, model){
                              qgev(p, c(param[, 1]), exp(c(param[, 2])), c(param[, 3]))
                    },
                    resid = function(o) {
                      p <- texmexMakeParams(coef(o), o$data$D)
                      shift <- (o$data$y - p[,1]) / exp(p[,2])
                      .log1prel(shift * p[,3]) * shift # standard Gumbel see Coles p.110 eq (6.6)
                    }, # Close resid

                    rl = function(m, param, model){
                      qgev(1/m, param[,1], exp(param[,2]), param[,3], lower.tail=FALSE)
                    } # Close lp
)

