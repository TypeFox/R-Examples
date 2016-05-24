# qml.R
#
# created: Feb/04/2015, NU
# last mod: Oct/26/2015, NU

#--------------- main functions ---------------

qml <- function(model, data, start, max.iter=150, 
                optimizer=c("nlminb", "optim"), neg.hessian=TRUE, ...) {

    if (model$info$num.eta > 1) stop("QML is not implemented for more than one eta (yet).")

    if (is.matrix(data)) {
        data <- data
    } else if (is.data.frame(data)) {
        data <- as.matrix(data)
    } else {
        stop("data need to be a matrix or a data frame.")
    }

    suppressWarnings(
    est <- mstep_qml(model=model, data=data, parameters=start,
                       neg.hessian=neg.hessian, optimizer=optimizer,
                       max.iter=max.iter, ...)
    )

    names(est$par) <- model$info$par.names

    if (sum(est$par - start) == 0) {
      stop("NA/NaN function evaluation. Please try different set of starting parameters.")
    }
    
    out <- list(model.class=class(model), coefficients=est$par,
                objective=-est$objective,
                convergence=est$convergence,
                neg.hessian=est$hessian,
                iterations=est$iterations,
                info=model$info[c("num.xi", "num.eta", "num.x", "num.y",
                "constraints", "num.classes")])

    class(out) <- "qmlEst"
    out
}

mu_qml <- function(model, data) {

    m <- model$matrices$class1 

    # extract x and y from data frame
    x <- data[, 1:model$info$num.x]
    y <- data[, (model$info$num.x + 1):dim(data)[2]]
    
    if (model$info$num.y > 1) {
        # transformation of y
        beta <- m$Lambda.y[-1,] 
        R <- cbind(-beta, diag(length(beta)))
        u <- y %*% t(R)
    } else {
        u <- 0
    }
    
    # Eqs 15, 16
    Sigma1 <- m$Phi - m$Phi %*% t(m$Lambda.x) %*% solve(m$Lambda.x %*%
              m$Phi %*% t(m$Lambda.x) + m$Theta.d) %*% m$Lambda.x %*% m$Phi 

    L1 <- m$Phi %*% t(m$Lambda.x) %*% solve(m$Lambda.x %*% 
          m$Phi %*% t(m$Lambda.x) + m$Theta.d)

    L2 <- -m$Theta.e[1,1] %*% t(beta) %*% solve(R %*% m$Theta.e %*% t(R))
    
    # m.x is unconditional
    # m.u mean vector for R %*% epsilon (0)
    # Eq 14: but mu.x <- 0, i.e., without means for xi
    mu.x <- m$nu.x + m$Lambda.x %*% m$tau
    #mu.u <- as.matrix(rep(0, length(beta)))
    mu.u <- R %*% m$nu.y
    
    N <- nrow(data)

    mtau.m <- matrix(rep(m$tau, N), nrow=model$info$num.xi, ncol=N, byrow=FALSE)
    mux.m <- matrix(rep(mu.x, N), nrow=N, ncol=model$info$num.x, byrow=TRUE)
    malpha.m <- matrix(rep(m$alpha, N), nrow=model$info$num.eta, ncol=N, byrow=TRUE)
    mnuy1.m <- matrix(rep(m$nu.y[1], N), nrow=model$info$num.eta, ncol=N, byrow=TRUE)
    muu.m <- matrix(rep(mu.u, N), nrow=N, ncol=model$info$num.y-model$info$num.eta, byrow=TRUE)
 
    # m.y1 is conditional given x and u
    # Eq 12 

    mu.y1 <- mnuy1.m + malpha.m + m$Gamma %*% (mtau.m + L1 %*% t(x - mux.m)) + 
           diag(t((mtau.m + L1 %*% t(x - mux.m))) %*% m$Omega %*% (mtau.m + L1 %*%
           t(x - mux.m))) + L2 %*% t(u - muu.m) + sum(diag(m$Omega %*% Sigma1))

    #mu.y1 <- sum(diag(m$Omega %*% Sigma1)) + c(m$alpha) + m$Gamma %*% 
    #         L1 %*% t(x) + diag(x %*% t(L1) %*% m$Omega %*% L1 %*% t(x)) +
    #         L2 %*% t(u)
    
    mu.xu <- c(mu.x, mu.u)    # mean for f2
    mu <- list(mu.xu, mu.y1)  

    mu
}


sigma_qml <- function(model, data) {

    m <- model$matrices$class1

    # extract x and y from data frame
    x <- data[, 1:model$info$num.x]
    y <- data[, (model$info$num.x + 1):dim(data)[2]]
    
    if (model$info$num.y > 1) {
        # transformation of y
        beta <- m$Lambda.y[-1,] 
        R <- cbind(-beta, diag(length(beta)))
        u <- y %*% t(R)
    } else {
        u <- 0
    }
    
    # Eqs 15, 16
    Sigma1 <- m$Phi - m$Phi %*% t(m$Lambda.x) %*% solve(m$Lambda.x %*%
              m$Phi %*% t(m$Lambda.x) + m$Theta.d) %*% m$Lambda.x %*% m$Phi 
    Sigma2 <- m$Psi + m$Theta.e[1,1] - m$Theta.e[1,1]^2 %*% t(beta) %*%
              solve(R %*% m$Theta.e %*% t(R)) %*% beta
    L1 <- m$Phi %*% t(m$Lambda.x) %*% solve(m$Lambda.x %*% m$Phi %*%
          t(m$Lambda.x) + m$Theta.d)
    L2 <- -m$Theta.e[1,1] %*% t(beta) %*% solve(R %*% m$Theta.e %*% t(R))

    mu.x <- m$nu.x + m$Lambda.x %*% m$tau

    N <- nrow(data)

    mux.m <- matrix(rep(mu.x, N), nrow=N, ncol=model$info$num.x, byrow=TRUE)
    mtau.m <- matrix(rep(m$tau, N), nrow=model$info$num.xi, ncol=N, byrow=FALSE)
    mGamma.m <- matrix(rep(m$Gamma, N), nrow=N, ncol=model$info$num.xi, byrow=TRUE)

    # Eq 18
    Sigma3 <- var.z(m$Omega, Sigma1)

    Sigma4 <- diag((((mGamma.m + t(mtau.m + L1 %*% t(x - mux.m))) %*%
              (m$Omega + t(m$Omega)))) %*% Sigma1 %*% t(((mGamma.m +
              t(mtau.m + L1 %*% t(x - mux.m))) %*% (m$Omega +
              t(m$Omega))))) + Sigma3
    
    # Eq 14
    # Cov(x), Cov(u)
    sigma.x <- m$Lambda.x %*% m$Phi %*% t(m$Lambda.x) + m$Theta.d
    sigma.u <- R %*% m$Theta.e %*% t(R)

    # Cov(x,u)
    p1 <- length(x[1,])
    p2 <- length(u[1,])
    sigma.xu <- matrix(0, p1+p2, p1+p2)
    sigma.xu[1:p1, 1:p1] <- sigma.x
    sigma.xu[(p1+1):(p1+p2), (p1+1):(p1+p2)] <- sigma.u
    
    # sigma.y1 is conditional given x and u
    # Eq 17
    # sigma.y1 <- diag((c(m$Gamma) + 2*x %*% t(L1) %*% m$Omega) %*% Sigma1
    #             %*% t(c(m$Gamma) + 2*x %*% t(L1) %*% m$Omega)) + Sigma2 + 
    #             Sigma3
    # --> original equation from paper: faulty!

    #sigma.y1 <- diag((c(m$Gamma) + x %*% t(L1) %*% (m$Omega + t(m$Omega))) %*% Sigma1
    #            %*% t(c(m$Gamma) + x %*% t(L1) %*% (m$Omega + t(m$Omega)))) + Sigma2 +
    #            Sigma3
 
    sigma.y1 <- Sigma2 + Sigma4

    sigma.xy <- list(sigma.xu, sigma.y1)

    sigma.xy
}

loglikelihood_qml <- function(parameters, model, data) {
    
    mod.filled <- fill_model(model = model, parameters = parameters)

    # extract x and y from data frame
    x <- data[, 1:model$info$num.x]
    y <- data[, (model$info$num.x + 1):dim(data)[2]]
       
    res <- 0
    
    mean.qml <- mu_qml(model = mod.filled, data=data)
    sigma.qml  <- sigma_qml(model = mod.filled, data=data)
    
    if (model$info$num.y > 1) {
        # transformation of y
        beta <- mod.filled$matrices$class1$Lambda.y[-1,] 
        R <- cbind(-beta, diag(length(beta)))
        u <- y %*% t(R)
    } else {
        u <- 0
    }

    # Eq 10: densities
    f2 <- dmvnorm(cbind(x, u), mean = mean.qml[[1]], sigma = sigma.qml[[1]])
    # original implementation: produces NaN when sds are negative
    f3 <- dnorm(y[,1], mean = mean.qml[[2]], sd = sqrt(sigma.qml[[2]]))

    lls <- sum(log(f2*f3))
    res <- res + lls
    
  return(-res)
}

mstep_qml <- function(model, parameters, data, neg.hessian=FALSE,
                      optimizer=c("nlminb", "optim"), max.iter=1,
                      control=list(), ...) {

    # optimizer
    optimizer <- match.arg(optimizer)

    if (optimizer == "nlminb") {

        if (is.null(control$iter.max)) {
            control$iter.max <- max.iter
        } else warning("iter.max is set for nlminb. max.iter will be ignored.")

        est <- nlminb(start=parameters, objective=loglikelihood_qml,
                      data=data, model=model,
                      upper=model$info$bounds$upper,
                      lower=model$info$bounds$lower, control=control, ...)
        names(est$par) <- model$info$par.names

        if (neg.hessian == TRUE){
            est$hessian <- fdHess(pars=est$par, fun=loglikelihood_qml,
                                        model=model, data=data)$Hessian
        }

    } else {

        if (is.null(control$maxit)){
            control$maxit <- max.iter
        } else warning("maxit is set for optim. max.iter will be ignored.")

        est <- optim(par=parameters, fn=loglikelihood_qml, data=data,
                     model=model, upper=model$info$bounds$upper,
                     lower=model$info$bounds$lower, method="L-BFGS-B",
                     control=control, ...) 
        # fit est to nlminb output
        names(est) <- gsub("value", "objective", names(est))
        #names(est) <- gsub("counts", "iterations", names(est))
        names(est$par) <- model$info$par.names

        if (neg.hessian == TRUE){
            est$hessian <- fdHess(pars=est$par, fun=loglikelihood_qml,
                                        model=model, data=data)$Hessian
  
        }
    }

    est
}


#--------------- helper functions ---------------

var.z <- function(Omega, Sigma1){

  ds <- dim(Sigma1)[1]
  varz <- 0
  
  # Eq 18
  for(i in 1:ds){
    for(j in 1:ds){
      for(k in 1:ds){
        for(l in 1:ds){
          varzij <- Omega[i,j]*Omega[k,l]*(Sigma1[i,k]*Sigma1[j,l]+Sigma1[i,l]*Sigma1[j,k])
          varz <- varz+varzij
        }
      }
    }
  }

  varz
}



