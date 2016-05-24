gpd.info <-
    # Compute the observed information matrix from a gpd object.
    # The expressions are given in Appendix A of Davison & Smith 1990.
    # Note we are using a simpler parameterisation in which phi = log(sigma)
    # and xi are both linear in their covariates. Xi is -k used in Davison and Smith.
	# If penalization is used, the calculation accounts for this, but the resulting
	# estimates of variance will be too low and bias might dominate MSE
function(o, method="observed"){
    if (class(o) != "evmOpt"){ stop("object must be of class 'evmOpt'") }
    if (method != "observed"){ stop("only 'observed' information is implemented") }

    x <- o$data$D$phi; z <- o$data$D$xi
    ns <- ncol(x); nk <- ncol(z)
    phi <- coef(o)[1:ns]
    xi <- coef(o)[(ns+1):(ns + nk)]

    phi.i <- colSums(phi * t(x))
    xi.i <- colSums(xi * t(z))
    w.i <- (o$data$y - o$threshold) / exp(phi.i)

    if (any(xi.i < -.50)){ warning("Fitted values of xi < -0.5") }

    # Second derivatives of penalties
    p <- matrix(0, nrow=ns+nk, ncol=ns+nk)
    if (o$penalty %in% c("gaussian", "quadratic")){ # note if Lasso penalty used then 2nd deriv is zero hence no term for this
        Si <- solve(o$priorParameters[[2]])
        for (i in 1:(ns+nk)){
            for (j in 1:(ns + nk)){
            p[i,j] <- 2*Si[i,j]
         }
        }
    }

    # Second and mixed derivatives of log-lik wrt coefficients of linear predictors

    d2li.dphi2 <- -(1 + 1/xi.i) * xi.i * w.i / (1 + xi.i*w.i)^2
    d2li.dphidxi <- 1/xi.i^2 * (1/(1 + xi.i*w.i) - 1) + (1+1/xi.i)*w.i/(1 + xi.i*w.i)^2
    d2li.dxi2 <- -2/xi.i^3 * log(1 + xi.i*w.i) + 2*w.i/(xi.i^2 * (1 + xi.i*w.i)) + (1 + 1/xi.i)*w.i^2/(1 + xi.i*w.i)^2

    # Matrix has 4 blocks, 2 of which are transposes of each other. Need block for phi parameters,
    # block for xi parameters and block for the cross of them.

    Ip <- matrix(0, ncol=ns, nrow=ns)
    for (u in 1:ns){
        for (v in 1:ns){
            Ip[u,v] <- -sum(x[,u] * x[,v] * d2li.dphi2)
        }
    }

    Ix <- matrix(0, ncol=nk, nrow=nk)
    for (s in 1:nk){
        for (t in 1:nk){
            Ix[s,t] <- -sum(z[,s] * z[,t] * d2li.dxi2)
        }
    }

    Ipx <- matrix(0, ncol=nk, nrow=ns)
    for (u in 1:ns){
        for (s in 1:nk){
            Ipx[u,s] <- -sum(z[,s] * x[,u] * d2li.dphidxi )
        }
    }

    i <- rbind( cbind(Ip, Ipx), cbind(t(Ipx), Ix))

    # return observed Information matrix.   Note that an estimate of the covariance matrix is given by the inverse of this matrix.
    i - p
}

