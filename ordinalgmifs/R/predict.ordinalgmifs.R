predict.ordinalgmifs <-
function(object, neww=NULL, newx=NULL, model.select="AIC", ...) {
    if (model.select == "AIC") {
        model.select = object$model.select
    }
    else if (model.select == "BIC") {
        model.select = which.min(object$BIC)
    }
    y <- object$y
    w <- object$w
    x <- object$x
    link <- object$link
    if (!is.null(x)) {
        alpha <- object$alpha[model.select, ]
        if (dim(w)[2] != 0) {
            if (is.null(dim(object$zeta))) {
                zeta <- object$zeta[model.select]
            }
            else {
                zeta <- object$zeta[model.select, ]
            }
        }
        beta <- object$beta[model.select, ]
        if (object$probability.model == "Stereotype") {
            if (is.null(dim(object$phi))) {
                phi <- object$phi[model.select]
            }
            else {
                phi <- object$phi[model.select, ]
            }
        }
    }
    else {
        alpha <- object$alpha
        if (dim(w)[2] != 0) {
            zeta <- object$zeta
        }
        if (object$probability.model == "Stereotype") {
            phi <- object$phi
        }
    }
    k <- length(unique(y))
    if (!is.null(neww)) 
        neww <- as.matrix(neww)
    if (!is.null(newx)) 
        newx <- as.matrix(newx)
    if (is.null(neww) & is.null(newx)) {
        neww <- object$w
        newx <- object$x
    }
    n <- max(dim(neww)[1], dim(newx)[1])
    if (!is.null(newx) & identical(newx, x)) {
        if (object$scale) {
            sd <- apply(newx, 2, sd)
            for (i in 1:dim(newx)[2]) {
                if (sd[i] == 0) {
                  newx[, i] <- scale(newx[, i], center = TRUE, 
                    scale = FALSE)
                }
                else {
                  newx[, i] <- scale(newx[, i], center = TRUE, 
                    scale = TRUE)
                }
            }
        }
    }
    else if (!is.null(newx) && object$scale) {
        newx <- rbind(x, newx)
        sd <- apply(newx, 2, sd)
        for (i in 1:dim(newx)[2]) {
            if (sd[i] == 0) {
                newx[, i] <- scale(newx[, i], center = TRUE, 
                  scale = FALSE)
            }
            else {
                newx[, i] <- scale(newx[, i], center = TRUE, 
                  scale = TRUE)
            }
        }
        newx <- matrix(newx[-(1:dim(x)[1]), ], ncol = dim(x)[2])
    }
    levels <- sort(unique(y))
    if (dim(w)[2] != 0) {
        if (is.null(x)) {
            Xb <- neww %*% zeta
        }
        else if (!is.null(x)) {
            Xb <- neww %*% zeta + newx %*% beta
        }
    }
    else if (!is.null(x)) {
        Xb <- newx %*% beta
    }
    else {
        Xb <- 0
    }
    if (object$probability.model == "Cumulative") {
        z <- matrix(ncol = k - 1, nrow = n)
        for (i in 1:(k - 1)) {
            z[, i] <- alpha[i] + Xb
        }
        pi <- matrix(ncol = k, nrow = n)
        if (link == "logit") {
            for (i in 1:k) {
                if (i == 1) {
                  pi[, i] <- exp(z[, i])/(1 + exp(z[, i]))
                }
                else if (i <= k - 1) {
                  pi[, i] <- exp(z[, i])/(1 + exp(z[, i])) - 
                    exp(z[, i - 1])/(1 + exp(z[, i - 1]))
                }
                else if (i == k) {
                  pi[, i] <- 1 - exp(z[, i - 1])/(1 + exp(z[, 
                    i - 1]))
                }
            }
        }
        else if (link == "probit") {
            for (i in 1:k) {
                if (i == 1) {
                  pi[, i] <- pnorm(z[, i])
                }
                else if (i <= k - 1) {
                  pi[, i] <- pnorm(z[, i]) - pnorm(z[, i - 1])
                }
                else if (i == k) {
                  pi[, i] <- 1 - pnorm(z[, i - 1])
                }
            }
        }
        else if (link == "cloglog") {
            for (i in 1:k) {
                if (i == 1) {
                  pi[, i] <- 1 - exp(-exp(z[, i]))
                }
                else if (i <= k - 1) {
                  pi[, i] <- exp(-exp(z[, i - 1])) - exp(-exp(z[, 
                    i]))
                }
                else if (i == k) {
                  pi[, i] <- exp(-exp(z[, i - 1]))
                }
            }
        }
    }
    else if (object$probability.model == "AdjCategory") {
        eta <- matrix(0, ncol = k - 1, nrow = n)
        for (i in 1:(k - 1)) {
            eta[, i] <- alpha[i] + Xb
        }
        if (n==1) {
        	eta.cumsum <- matrix(cumsum(eta), nrow = nrow(eta), byrow = T)
        } else {
        	eta.cumsum <- matrix(apply(eta, 1, cumsum), nrow = nrow(eta), byrow = T)        	
        }
        numer <- rep(0, dim(eta.cumsum)[1])
        for (i in 1:dim(eta.cumsum)[2]) {
            numer <- numer + exp(eta.cumsum[, i])
        }
        pi <- matrix(0, ncol = k, nrow = n)
        pi[, 1] <- 1 - numer/(1 + numer)
        for (i in 2:k) {
            pi[, i] <- exp(eta.cumsum[, i - 1] + log(pi[, 1]))
        }
    }
    else if (object$probability.model == "ForwardCR") {
        pi <- matrix(0, nrow = n, ncol = k)
        if (link == "logit" | link == "cloglog") {
            pi[, 1] <- G(alpha[1] + Xb, link)
            pi[, 2] <- G(alpha[2] + Xb, link) * (1 - pi[, 1])
            if (k > 3) {
                for (i in 3:(k - 1)) {
                  if (n==1) {
                  	pi[, i] <- G(alpha[i] + Xb, link) * (1 - sum(pi[,1:(i - 1)]))
                  } else if (n>1) {
                    pi[, i] <- G(alpha[i] + Xb, link) * (1 - rowSums(pi[,1:(i - 1)]))
                  }
                }
            }
        }
        else if (link == "probit") {
            pi[, 1] <- pnorm(alpha[1] + Xb)
            pi[, 2] <- pnorm(alpha[2] + Xb) * (1 - pi[, 1])
            if (k > 3) {
                for (i in 3:(k - 1)) {
                  if (n==1) {
                  	pi[, i] <- pnorm(alpha[i] + Xb) * (1 - sum(pi[,1:(i - 1)]))
                  } else {
                  	pi[, i] <- pnorm(alpha[i] + Xb) * (1 - rowSums(pi[,1:(i - 1)]))                 	
                  }
                }
            }
        }
        if (n==1) {
            pi[, k] <- 1 - sum(pi[, 1:(k - 1)])
        } else {
            pi[, k] <- 1 - rowSums(pi[, 1:(k - 1)])
		}
    }
    else if (object$probability.model == "BackwardCR") {
        pi <- matrix(0, nrow = n, ncol = k)
        if (link == "logit" | link == "cloglog") {
            pi[, k] <- G(alpha[k - 1] + Xb, link)
            pi[, k - 1] <- G(alpha[k - 2] + Xb, link) * (1 - 
                pi[, k])
            if (k > 3) {
                for (i in (k - 2):2) {
                  if (n==1) {
                  	pi[, i] <- G(alpha[i - 1] + Xb, link) * (1 - sum(pi[, k:(i + 1)]))
                  } else {
                   	pi[, i] <- G(alpha[i - 1] + Xb, link) * (1 - rowSums(pi[, k:(i + 1)]))                 	
                  }
                }
            }
        }
        else if (link == "probit") {
            pi[, k] <- pnorm(alpha[k - 1] + Xb)
            pi[, k - 1] <- pnorm(alpha[k - 2] + Xb) * (1 - pi[, k])
            if (k > 3) {
                for (i in (k - 2):2) {
                  if (n==1) {
                  	pi[, i] <- pnorm(alpha[i - 1] + Xb) * (1 - sum(pi[, k:(i + 1)]))                  	
                  } else {
                  	pi[, i] <- pnorm(alpha[i - 1] + Xb) * (1 - rowSums(pi[, k:(i + 1)]))
                  }
                }
            }
        }
        if (n==1) {
        	pi[, 1] <- 1 - sum(pi[, k:2])
        } else {
        	pi[, 1] <- 1 - rowSums(pi[, k:2])        	
        }
    }
    else if (object$probability.model == "Stereotype") {
        eta <- matrix(0, ncol = k - 1, nrow = n)
        for (i in 1:(k - 1)) {
            eta[, i] <- exp(alpha[i] + phi[i] * Xb)
        }
        if (n==1) {
        	pik <- 1 - sum(eta)/(1 + sum(eta))
        } else {
        	pik <- 1 - rowSums(eta)/(1 + rowSums(eta))       	
        }
        pi <- matrix(0, ncol = k, nrow = n)
        pi[, k] <- pik
        pi[, 1:(k - 1)] <- eta * pik
    }
    class <- levels[apply(pi, 1, which.max)]
    list(predicted = pi, class = class)
}
