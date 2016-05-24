logdet <- function(A)
  {
    return(2 * sum(log(diag(chol(A)))))
  }

make.D <- function(x.1, x.2)
  {
    ##AFL, 31.07.13
    ##returns the distance matrix used by a Gaussian Process
    n.1 <- dim(x.1)[1]
    n.2 <- dim(x.2)[1]
    d <- dim(x.1)[2]
    D <- matrix(0, n.1, n.2)
    for(j in 1:n.2)
          {
            D[,j] <- sqrt(colSums( (t(x.1) - x.2[j,])^2))
          }
    return(D)
  }

gev.z.p <- function(p,mu, sigma, xi)
  {
    if(abs(xi[1]) < 1e-5)
      {
        return(mu - sigma * log( -log(1 - p)))
      }else{
        return(mu - sigma/xi * (1 - (-log(1 - p))^(-xi)))
      }
  }


gev.like <- function(Y,mu,kappa,xi)
  {
    h <- 1 + xi * kappa * (Y - mu)
    l <- log(kappa) - (xi + 1)/xi * log(h) - h^(-1/xi)
    return(l)
  }

dmvnorm <- function(x, mu, Sigma)
  {
    p <- length(x)
    dd <- 2 * sum(log(diag(chol(Sigma))))
    K <- solve(Sigma)
    return( -1/2 * t(x - mu) %*% K %*% (x-mu) - dd/2)
  }

##----------------- Generic Functions for the linear model parts ----------------
gev.update.M <- function(Y,X,M,alpha,lambda,D, beta.0, Omega.0)
  {
    p <- dim(X)[2]
    M.curr <- M
    ww <- sample(2:p,1)
    if(any(ww == M.curr))
      {
        M.new <- M.curr[-which(M.curr == ww)]
      }else{
        M.new <- sort(c(M.curr, ww))
      }

    if(length(M.new) == 0)
      {
        return(M)
      }

    ##----------- Shared Objects ------------
    C <- 1/alpha * exp(-D/lambda)
    diag(C) <- diag(C) + 1e-5
    C.inv <- solve(C)
    ##---------------------------------------
    
    ##--------- Get statistics for current model -------
    p.M.curr <- length(M.curr)
    X.M.curr <- X[,M.curr,drop=FALSE]
    Omega.0.M <- Omega.0[M.curr,M.curr,drop=FALSE]
    beta.0.M <- beta.0[M.curr]
    Xi <- t(X.M.curr) %*% C.inv %*% X.M.curr + Omega.0.M
    Xi.inv <- solve(Xi)
    theta.hat <- Xi.inv %*% (Omega.0.M %*% beta.0.M + t(X.M.curr) %*% C.inv %*% Y)
    score.curr <- 0.5 * t(theta.hat) %*% Xi %*% theta.hat - 0.5 * logdet(Xi)
    ##----------------------------------------------------

    ##--------- Get statistics for new model -------
    p.M.new <- length(M.new)
    X.M.new <- X[,M.new,drop=FALSE]
    Omega.0.M <- Omega.0[M.new, M.new, drop=FALSE]
    beta.0.M <- beta.0[M.new]
    Xi <- t(X.M.new) %*% C.inv %*% X.M.new + Omega.0.M
    Xi.inv <- solve(Xi)
    theta.hat <- Xi.inv %*% (Omega.0.M %*% beta.0.M + t(X.M.new) %*% C.inv %*% Y)
    score.new <- 0.5 * t(theta.hat) %*% Xi %*% theta.hat - 0.5 * logdet(Xi)
    ##----------------------------------------------------

    mh <- score.new - score.curr
    if(log(runif(1)) < mh)
      {
       M.curr <- M.new
      }
    
    return(M.curr)
  }

gev.update.theta <- function(Y,X,M,alpha,lambda,D, beta.0, Omega.0)
  {
    p <- dim(X)[2]
    p.M <- length(M)
    X.M <- X[,M,drop=FALSE]

    C <- 1/alpha * exp(-1/lambda * D)
    diag(C) <- diag(C) + 1e-5
    C.inv <- solve(C)
    Omega.0.M <- Omega.0[M,M,drop=FALSE]
    beta.0.M <- beta.0[M]
    Xi <- t(X.M) %*% C.inv %*% X.M + Omega.0.M
    Xi.inv <- solve(Xi)
    theta.hat <- Xi.inv %*% (Omega.0.M %*% beta.0.M + t(X.M) %*% C.inv %*% Y)
    theta <- rep(0, p)
    
    theta[M] <- t(chol(Xi.inv)) %*% rnorm(p.M) + theta.hat

    return(theta)
    
  }
##----------- End generic linear model functions -------------------------

##---------- Updating Mu random effects ----------------------------------
f.prime <- function(tau,tau.hat,varsigma,xi,kappa,R)
  {
    h <- 1 + xi  * kappa * (R-tau)
    if(any(h < 0))return(-Inf)
    L <- (xi + 1) * kappa * h^(-1) - kappa * h^(-1/xi - 1)
    res <- sum(L) - (tau - tau.hat)/varsigma
    return(res)
  }

f.double.prime <- function(tau, tau.hat, varsigma, xi, kappa, R)
  {
    h <- 1 + xi * kappa * (R-tau)
    if(any(h < 0))return(-Inf)
    L <- xi * (xi + 1) * kappa^2 * h^(-2) - (xi + 1) * kappa^2 * h^(-1/xi - 2)
    res <- sum(L) - 1/varsigma
    return(res)
  }

gev.update.tau.mu <- function(G)
  {
    n.s <- G$n.s
    theta.mu <- G$theta.mu

    C <- exp(-1/G$lambda.mu * G$D) / G$alpha.mu
    diag(C) <- diag(C) + 1e-5
    C.inv <- solve(C)
    tau <- G$tau.mu
    G$accept.tau.mu <- rep(0, n.s)
    for(s in 1:n.s)
      {
##        print(s)
        Y.s <- G$Y.list[[s]]
        fit.s <- sum(theta.mu * G$X[s,])
        R.s <- Y.s - fit.s
        kappa.s <- sum(G$theta.kappa * G$X[s,]) + G$tau.kappa[s]
        xi.s <- sum(G$theta.xi * G$X[s,]) + G$tau.xi[s]
        tau.hat.s <- sum(-C.inv[s,-s]/C.inv[s,s] * tau[-s])
        varsigma.s <- 1/C.inv[s,s] ##precision matrix stuff

        f.s <- f.prime(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)
        ff.s <- f.double.prime(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)
        if(is.finite(f.s) && (ff.s < 0))
          {
            b <- f.s - ff.s * tau[s]
            d <- -ff.s
            tau.new <- rnorm(1, b/d, sd =sqrt(1/d))
            
            f.s <- f.prime(tau.new, tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)
            ff.s <- f.double.prime(tau.new, tau.hat.s, varsigma.s, xi.s, kappa.s, R.s)
            if(is.finite(f.s) && (ff.s < 0))
              {
                b.new <- f.s - ff.s * tau.new
                d.new <- -ff.s
                
                mu.curr <- sum(G$theta.mu * G$X[s,]) + tau[s]
                mu.new <- sum(G$theta.mu * G$X[s,]) + tau.new
                
                L.curr <- sum(gev.like(G$Y.list[[s]], mu.curr, kappa.s, xi.s))
                L.new <- sum(gev.like(G$Y.list[[s]], mu.new, kappa.s, xi.s))
                prior.curr <- dnorm(tau[s], tau.hat.s, sd=sqrt(varsigma.s), log=TRUE)
                prior.new <- dnorm(tau.new, tau.hat.s, sd=sqrt(varsigma.s), log=TRUE)
                prop.curr <- dnorm(tau.new, b/d, sqrt(1/d), log=TRUE)
                prop.new <- dnorm(tau[s], b.new/d.new, sqrt(1/d.new), log=TRUE)
                
                alpha <- L.new - L.curr + prior.new - prior.curr + prop.new - prop.curr
                if(log(runif(1)) < alpha)
                  {
                    G$accept.tau.mu[s] <- 1
                    tau[s] <- tau.new
                  }
              }
          }
      }
    G$tau.mu <-tau
    return(G)
  }


##------------------ End Updating Mu Random Effects ----------------------------

##---------  Updating Kappa Random Effects -------------------------------------
g.prime <- function(tau,tau.hat,varsigma,xi,kappa.hat,eps)
  {
    h <- 1 + xi * (kappa.hat + tau) * eps
    if(any(h < 0))return(-Inf)
    L <- 1/(kappa.hat + tau) - (xi + 1) * eps * h^(-1) + eps * h^(-1/xi - 1)
    res <- sum(L) - (tau - tau.hat)/varsigma
    return(res)
  }

g.double.prime <- function(tau, tau.hat, varsigma, xi, kappa.hat, eps)
  {
    h <- 1 + xi * (kappa.hat + tau) * eps
    if(any(h <0))return(-Inf)
    L <- -(kappa.hat + tau)^(-2) + (xi + 1) * xi * eps^2 * h^(-2) - eps^2 * (xi + 1) * h^(-xi^(-1) - 2)
    res <- sum(L) - 1/varsigma
    return(res)
  }

gev.update.tau.kappa <- function(G)
  {

    n.s <- G$n.s
    theta.kappa <- G$theta.kappa

    C <- 1/G$alpha.kappa * exp(-1/G$lambda.kappa * G$D)
    diag(C) <- diag(C) + 1e-5
    C.inv <- solve(C)
    tau <- G$tau.kappa
    accept <- 0
    G$accept.tau.kappa <- rep(0,n.s)
    for(s in 1:n.s)
      {
        Y.s <- G$Y.list[[s]]
        mu.s <- sum(G$theta.mu * G$X[s,]) + G$tau.mu[s]
        eps.s <- Y.s - mu.s
        kappa.hat <- sum(G$theta.kappa * G$X[s,])
        xi.s <- sum(G$theta.xi * G$X[s,]) + G$tau.xi[s]

        tau.hat.s <- sum(-C.inv[s,-s]/C.inv[s,s] * tau[-s])
        varsigma.s <- 1/C.inv[s,s] ##precision matrix stuff

        g.s <- g.prime(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.hat, eps.s)
        gg.s <- g.double.prime(tau[s], tau.hat.s, varsigma.s, xi.s, kappa.hat, eps.s)
        if(is.finite(g.s) && (gg.s < 0))
           {
             b <- g.s - gg.s * tau[s]
             d <- -gg.s
             tau.new <- rnorm(1, b/d, sd =sqrt(1/d))
             
             g.s <- g.prime(tau.new, tau.hat.s, varsigma.s, xi.s, kappa.hat, eps.s)
             gg.s <- g.double.prime(tau.new, tau.hat.s, varsigma.s, xi.s, kappa.hat, eps.s)
             if(is.finite(g.s) && (gg.s < 0))
               {
                 b.new <- g.s - gg.s * tau.new
                 d.new <- -gg.s
                 
                 kappa.curr <- sum(G$theta.kappa * G$X[s,]) + tau[s]
                 kappa.new <- sum(G$theta.kappa * G$X[s,]) + tau.new
                 if(kappa.new > 0)
                   {
                     L.curr <- sum(gev.like(G$Y.list[[s]], mu.s, kappa.curr, xi.s))
                     L.new <- sum(gev.like(G$Y.list[[s]], mu.s, kappa.new, xi.s))
                     prior.curr <- dnorm(tau[s], tau.hat.s, sd=sqrt(varsigma.s), log=TRUE)
                     prior.new <- dnorm(tau.new, tau.hat.s, sd=sqrt(varsigma.s), log=TRUE)
                     prop.curr <- dnorm(tau.new, b/d, sqrt(1/d), log=TRUE)
                     prop.new <- dnorm(tau[s], b.new/d.new, sqrt(1/d.new), log=TRUE)
                     
                     alpha <- L.new - L.curr + prior.new - prior.curr + prop.new - prop.curr
                     if(log(runif(1)) < alpha)
                       {
                         G$accept.tau.kappa[s] <- 1
                         accept <- accept + 1
                         tau[s] <- tau.new
                       }
                   }
               }
           }
      }
    G$tau.kappa <-tau
    return(G)
    
  }
##----------------- End Updating Kappa Random Effects -------------------------------------

##-------------- Updating Xi Random Effects -----------------------------------------------
j.prime <- function(tau, tau.hat, varsigma, kappa, xi.hat, eps)
  {
    h <- 1 + kappa * eps * (xi.hat + tau)
    if(any(h < 0))return(-Inf)
    f.1 <- (xi.hat + tau + 1)/(xi.hat + tau) * log(h)
    f.1.dot <- -log(h)/(xi.hat + tau)^2 + (xi.hat + tau + 1)/(xi.hat + tau) * h^(-1) * eps * kappa
    f.2 <-  exp(-(xi.hat + tau)^(-1) * log(h))
    f.2.dot <- f.2 * (log(h)/(xi.hat + tau)^2 - h^(-1) * eps * kappa / (xi.hat + tau) )
    res <- -sum(f.1.dot) - sum(f.2.dot) - (tau - tau.hat)/varsigma
    return(res)
  }

j.double.prime <- function(tau, tau.hat, varsigma, kappa, xi.hat, eps)
  {
    ## What an ugly function
    h <- 1 + kappa * eps * (xi.hat + tau)
    if(any(h < 0))return(-Inf)
    f.1 <- (xi.hat + tau + 1)/(xi.hat + tau) * log(h)
    f.2 <-  exp(-(xi.hat + tau)^(-1) * log(h))
    
    f.1.dot <- -log(h)/(xi.hat + tau)^2 + (xi.hat + tau + 1)/(xi.hat + tau) * h^(-1) * eps * kappa
    f.2.dot <- f.2 * ( log(h)/(xi.hat + tau)^2 - h^(-1) * eps * kappa / (xi.hat + tau) )
    g.1.dot <- -2 * (xi.hat + tau)^(-3) * log(h) + (xi.hat + tau)^(-2) * h^(-1) * kappa * eps
    g.2.dot <- -h^(-1) * eps * kappa * (xi.hat + tau)^(-2) - (xi.hat + tau + 1)/(xi.hat + tau) * h^(-2) * eps^2 * kappa^2
    g.3.dot.1 <- f.2.dot * (log(h) * (xi.hat + tau)^(-2) )
    g.3.dot.2 <- f.2 * ( -2 * log(h) * (xi.hat + tau)^(-3) + h^(-1) * kappa * eps * (xi.hat + tau)^(-2) )
    g.3.dot <- g.3.dot.1 + g.3.dot.2
    g.4.dot.1 <- f.2.dot * (h^(-1) * kappa * eps * (xi.hat + tau)^(-1))
    g.4.dot.2 <- -f.2 * eps * kappa * ( h^(-1) * (xi.hat + tau)^(-2) + h^(-2) * eps * kappa * (xi.hat + tau)^(-1) )
    g.4.dot <- g.4.dot.1 + g.4.dot.2
    res <- sum(g.1.dot) - sum(g.2.dot) - sum(g.3.dot) + sum(g.4.dot) - 1/varsigma
    return(res)
  }

gev.update.tau.xi <- function(G)
  {

    n.s <- G$n.s
    theta.xi <- G$theta.xi

    C <- 1/G$alpha.xi * exp(-1/G$lambda.xi * G$D)
    diag(C) <- diag(C) + 1e-5
    C.inv <- solve(C)
    tau <- G$tau.xi
    accept <- 0
    G$accept.tau.xi <- rep(0,n.s)
    for(s in 1:n.s)
      {
        Y.s <- G$Y.list[[s]]
        mu.s <- sum(G$theta.mu * G$X[s,]) + G$tau.mu[s]
        eps.s <- Y.s - mu.s
        kappa.s <- sum(G$theta.kappa * G$X[s,]) + G$tau.kappa[s]
        xi.hat <- sum(G$theta.xi * G$X[s,])

        tau.hat.s <- sum(-C.inv[s,-s]/C.inv[s,s] * tau[-s])
        varsigma.s <- 1/C.inv[s,s] ##precision matrix stuff

        j.s <- j.prime(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
        jj.s <- j.double.prime(tau[s], tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
        if(is.finite(j.s) && (jj.s < 0))
          {
            b <- j.s - jj.s * tau[s]
            d <- -jj.s
            tau.new <- rnorm(1, b/d, sd =sqrt(1/d))
            temp <- 1 + kappa.s * (xi.hat + tau.new) * eps.s
            if(all(temp > 0))
              {
                j.s <- j.prime(tau.new, tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
                jj.s <- j.double.prime(tau.new, tau.hat.s, varsigma.s, kappa.s, xi.hat, eps.s)
                if(is.finite(j.s) && (jj.s < 0))
                  {
                    b.new <- j.s - jj.s * tau.new
                    d.new <- -jj.s
                    
                    xi.curr <- xi.hat + tau[s]
                    xi.new <- xi.hat + tau.new
                    
                    L.curr <- sum(gev.like(G$Y.list[[s]], mu.s, kappa.s, xi.curr))
                    L.new <- sum(gev.like(G$Y.list[[s]], mu.s, kappa.s, xi.new))
                    prior.curr <- dnorm(tau[s], tau.hat.s, sd=sqrt(varsigma.s), log=TRUE)
                    prior.new <- dnorm(tau.new, tau.hat.s, sd=sqrt(varsigma.s), log=TRUE)
                    prop.curr <- dnorm(tau.new, b/d, sqrt(1/d), log=TRUE)
                    prop.new <- dnorm(tau[s], b.new/d.new, sqrt(1/d.new), log=TRUE)
                    
                    alpha <- L.new - L.curr + prior.new - prior.curr + prop.new - prop.curr
                    if(log(runif(1)) < alpha)
                      {
                        G$accept.tau.xi[s] <- 1
                        accept <- accept + 1
                        tau[s] <- tau.new
                      }
                  }
              }
          }
      }
    G$tau.xi <-tau
    return(G)
  }
##--------------- End Updating Xi Random Effects -----------------------------------------

##--------------  Updating Gaussian Process Hyper Parameters -----------------------------
l.prime <- function(tau, alpha, lambda,D,a,b)
  {
    E.l <- exp(-1/lambda * D)
    diag(E.l) <- diag(E.l) + 1e-5
    E.inv <- solve(E.l)
    F.l <- 1/lambda^2 * D * E.l
    M.l <- E.inv %*% (-F.l) %*% E.inv

    res <- -0.5 * sum(diag(E.inv %*% F.l)) - 0.5 * alpha * t(tau) %*% M.l %*% tau - b + (a - 1)/lambda
    return(res[1])
  }

l.double.prime <- function(tau, alpha, lambda, D,a,b)
  {
    E.l <- exp(-1/lambda * D)
    diag(E.l) <- diag(E.l) + 1e-5
    E.inv <- solve(E.l)
    F.l <- 1/lambda^2 * D * E.l
    M.l <- E.inv %*% (-F.l) %*% E.inv
    G.l <- -2/lambda^3 * (D * E.l) + 1/lambda^2 * (D * F.l)
    L.l <- M.l %*% F.l + E.inv %*% G.l
    N.l <- M.l %*% (-F.l) %*% E.inv + E.inv %*% (-G.l) %*% E.inv + E.inv %*% (-F.l) %*% M.l
    res <- -0.5 * sum(diag(L.l)) - 0.5 * alpha * t(tau) %*% N.l %*% tau - (a - 1)*lambda^(-2)
    return(res[1])
  }
 
gev.update.lambda <- function(tau, alpha, lambda, D, a, b)
  {

    l.curr <- l.prime(tau, alpha, lambda, D, a, b)
    l.double.curr <- l.double.prime(tau, alpha, lambda, D, a, b)
    b.curr <- l.curr - l.double.curr * lambda
    d.curr <- -l.double.curr
    if(d.curr > 0)
      {
        lambda.new <- rnorm(1, b.curr/d.curr, sd = sqrt(1/d.curr))
        if(lambda.new > 0)
          {
            l.new <- l.prime(tau, alpha, lambda.new, D, a, b)
            l.double.new <- l.double.prime(tau, alpha, lambda.new, D, a, b)
            b.new <- l.new - l.double.new * lambda.new
            d.new <- -l.double.new
            if(d.new > 0)
              {
                E.l.curr <- exp(-1/lambda * D)
                diag(E.l.curr) <- diag(E.l.curr) + 1e-5
                E.l.curr.inv <- solve(E.l.curr)
                E.l.new <- exp(-1/lambda.new * D)
                diag(E.l.new) <- diag(E.l.new) + 1e-5
                E.l.new.inv <- solve(E.l.new)
                
                L.curr <- -0.5 * alpha * t(tau) %*% E.l.curr.inv %*% tau - 0.5 * logdet(E.l.curr)
                L.new <- -0.5 * alpha * t(tau) %*% E.l.new.inv %*% tau - 0.5 * logdet(E.l.new)
                prior.curr <- dgamma(lambda, a,b, log=TRUE)
                prior.new <- dgamma(lambda.new, a,b, log=TRUE)
                prop.curr <- dnorm(lambda.new, b.curr/d.curr, sd = sqrt(1/d.curr), log=TRUE)
                prop.new <- dnorm(lambda, b.new/d.new, sd = sqrt(1/d.new), log = TRUE)
                
                mh <- L.new - L.curr + prior.new - prior.curr + prop.new - prop.curr
                if(log(runif(1)) < mh)
                  {
                    lambda <- lambda.new
                  }
              }
          }
      }
    return(lambda)
  }

gev.update.hyper <- function(G)
  {
    ##--------- Update Mu -------------
    tau <- G$tau.mu
    lambda <- G$lambda.mu
    E.lambda <- exp(-G$D/lambda)
    diag(E.lambda) <- diag(E.lambda) + 1e-5
    n.s <- G$n.s
    G$alpha.mu <- rgamma(1,n.s/2 + G$prior$mu$alpha.a, (t(tau) %*% solve(E.lambda) %*% tau)/2 + G$prior$mu$alpha.b)
    G$lambda.mu <- gev.update.lambda(tau, G$alpha.mu, lambda, G$D, G$prior$mu$lambda.a, G$prior$mu$lambda.b)
    ##-----------------------------------

    ##--------- Update Kappa -------------
    tau <- G$tau.kappa
    lambda <- G$lambda.kappa
    E.lambda <- exp(-G$D/lambda)
    diag(E.lambda) <- diag(E.lambda) + 1e-5
    n.s <- G$n.s
    G$alpha.kappa <- rgamma(1,n.s/2 + G$prior$kappa$alpha.a, (t(tau) %*% solve(E.lambda) %*% tau)/2 + G$prior$kappa$alpha.b)
    G$lambda.kappa <- gev.update.lambda(tau, G$alpha.kappa, lambda, G$D, G$prior$kappa$lambda.a, G$prior$kappa$lambda.b)
    ##-----------------------------------

    ##--------- Update Xi -------------
    if(!G$fixed.xi)
      {
        tau <- G$tau.xi
        lambda <- G$lambda.xi
        E.lambda <- exp(-G$D/lambda)
        diag(E.lambda) <- diag(E.lambda) + 1e-5
        n.s <- G$n.s
        G$alpha.xi <- rgamma(1,n.s/2 + G$prior$xi$alpha.a, (t(tau) %*% solve(E.lambda) %*% tau)/2 + G$prior$xi$alpha.b)
        G$lambda.xi <- gev.update.lambda(tau, G$alpha.xi, lambda, G$D, G$prior$xi$lambda.a, G$prior$xi$lambda.b)
      }
    ##-----------------------------------

    return(G)

  }
##----------------- End Update Hyper parameters ----------------------------------------

##----------------- Initialize Object and helper functions -----------------------------
gp.like.lambda <- function(lambda, alpha, tau, D)
  {
    if(alpha < 0) return(-Inf)
    if(lambda < 0) return(-Inf)
    C <- exp(-1/lambda * D)/alpha
    diag(C) <- diag(C) + 1e-5
    C.inv <- solve(C)
    l <- -1/2 * t(tau) %*% C.inv %*% tau - 0.5 * logdet(C) + dgamma(lambda,1,1,log=TRUE)
    return(-l)
  }

gev.init <- function(Y.list, X.all,S, prior.user, full, fixed.xi)
  {
    ## I literally have no idea how this function got so long.
    G <- NULL
    n.s <- length(Y.list)
    G$n.s <- n.s
    D <- make.D(S,S)
    p.max <- dim(X.all)[2]
    p <- p.max
    G$p <- p
    G$Y.list <- Y.list
    G$X <- X.all
    G$Y.vec <- NULL
    G$I.vec <- NULL
    G$X.long <- NULL
    G$full <- full
    for(i in 1:length(Y.list))
      {
        n.i <- length(Y.list[[i]])
        G$Y.vec <- c(G$Y.vec, Y.list[[i]])
        G$I.vec <- c(G$I.vec, rep(i, n.i))
        G$X.long <- rbind(G$X.long, matrix(X.all[i,], n.i, p, byrow=TRUE))
      }
    G$D <- D

    ML <- matrix(unlist(lapply(Y.list,"gevmle")),ncol=3, byrow=TRUE)
    G$prior <- NULL
    
    ##----------- Set Mu model -------------------------
    mu.temp <- ML[,1]
    if(G$full){
      G$M.mu <- 1:p.max
    }else{
      G$M.mu <- sort(unique(c(1,sort(sample(1:p.max,sum(rbinom(p.max,1,.5)), replace=FALSE)))))
    }
    p.M <- length(G$M.mu)
    X.M <- G$X[,G$M.mu,drop=FALSE]
    Xi.mu <- t(X.M) %*% X.M + diag(p.M)
    Xi.mu.inv <- solve(Xi.mu)
    theta.hat <- Xi.mu.inv %*% t(X.M) %*% mu.temp
    G$theta.mu <- rep(0, p.max)
    G$theta.mu[G$M.mu] <- theta.hat
    tau <- mu.temp - X.M %*% theta.hat
    G$tau.mu <- tau
    G$alpha.mu <- 1/var(tau)[1]
    temp <- optimize(gp.like.lambda,interval=c(0,1e4),G$alpha.mu, tau, D)
    G$lambda.mu <- temp[[1]]

    G$prior$mu <- NULL
    if(is.null(prior.user$mu$beta.0))
      {
        G$prior$mu$beta.0 <- rep(0, p)
      }else{
        G$prior$mu$beta.0 <- prior.user$mu$beta.0
      }
    if(is.null(prior.user$mu$Omega.0))
      {
        G$prior$mu$Omega.0 <- diag(p)
      }else{
        G$prior$mu$Omega.0 <- prior.user$mu$Omega.0
      }
    if(is.null(prior.user$mu$alpha.a))
      {
        G$prior$mu$alpha.a <- 1
      }else{
        G$prior$mu$alpha.a <- prior.user$mu$alpha.a
      }
    if(is.null(prior.user$mu$alpha.b))
      {
        G$prior$mu$alpha.b <- 1
      }else{
        G$prior$mu$alpha.b <- prior.user$mu$alpha.b
      }
    if(is.null(prior.user$mu$lambda.a))
      {
        G$prior$mu$lambda.a <- 1
      }else{
        G$prior$mu$lambda.a <- prior.user$mu$lambda.a
      }
    if(is.null(prior.user$mu$lambda.b))
      {
        G$prior$mu$lambda.b <- 1
      }else{
        G$prior$mu$lambda.b <- prior.user$mu$lambda.b
      }
    ##--------------------------------------------------

    ##----------- Set Kappa model -------------------------
    kappa.temp <- 1/ML[,2]
    if(G$full)
      {
        G$M.kappa <- 1:p.max
      }else{
        G$M.kappa <- sort(unique(c(1,sort(sample(1:p.max,sum(rbinom(p.max,1,.5)), replace=FALSE)))))
      }
    p.M <- length(G$M.kappa)
    X.M <- G$X[,G$M.kappa,drop=FALSE]
    Xi.kappa <- t(X.M) %*% X.M + diag(p.M)
    Xi.kappa.inv <- solve(Xi.kappa)
    theta.hat <- Xi.kappa.inv %*% t(X.M) %*% kappa.temp
    G$theta.kappa <- rep(0, p.max)
    G$theta.kappa[G$M.kappa] <- theta.hat
    tau <- kappa.temp - X.M %*% theta.hat
    G$tau.kappa <- tau
    G$alpha.kappa <- 1/var(tau)[1]
    temp <- optimize(gp.like.lambda,interval=c(0,1e4),G$alpha.kappa, tau, D)
    G$lambda.kappa <- temp[[1]]
    if(is.null(prior.user$kappa$beta.0))
      {
        G$prior$kappa$beta.0 <- rep(0, p.max)
      }else{
        G$prior$kappa$beta.0 <- prior.user$kappa$beta.0
      }
    if(is.null(prior.user$kappa$Omega.0))
      {
        G$prior$kappa$Omega.0 <- diag(p.max)
      }else{
        G$prior$kappa$Omega.0 <- prior.user$kappa$Omega.0
      }
    if(is.null(prior.user$kappa$alpha.a))
      {
        G$prior$kappa$alpha.a <- 1
      }else{
        G$prior$kappa$alpha.a <- prior.user$kappa$alpha.a
      }
    if(is.null(prior.user$kappa$alpha.b))
      {
        G$prior$kappa$alpha.b <- 1
      }else{
        G$prior$kappa$alpha.b <- prior.user$kappa$alpha.b
      }
    if(is.null(prior.user$kappa$lambda.a))
      {
        G$prior$kappa$lambda.a <- 1
      }else{
        G$prior$kappa$lambda.a <- prior.user$kappa$lambda.a
      }
    if(is.null(prior.user$kappa$lambda.b))
      {
        G$prior$kappa$lambda.b <- 1
      }else{
        G$prior$kappa$lambda.b <- prior.user$kappa$lambda.b
      }
    ##--------------------------------------------------

    ##---------- Set Xi model --------------------------
    if(!is.null(fixed.xi))
      {
        G$fixed.xi <- TRUE
        G$M.xi <- 1
        G$theta.xi <- c(fixed.xi, rep(0, p.max - 1))
        G$tau.xi <- rep(0,n.s)
        G$alpha.xi <- 1
        G$lambda.xi <- 1
      }else{
        G$fixed.xi <- FALSE
        xi.temp <- ML[,3]/5
        if(G$full)
          {
            G$M.xi <- 1:p.max
          }else{
            G$M.xi <- sort(unique(c(1,sort(sample(1:p.max,sum(rbinom(p.max,1,.5)), replace=FALSE)))))
          }
        p.M <- length(G$M.xi)
        X.M <- G$X[,G$M.xi,drop=FALSE]
        Xi.xi <- t(X.M) %*% X.M + diag(p.M)
        Xi.xi.inv <- solve(Xi.xi)
        theta.hat <- Xi.xi.inv %*% t(X.M) %*% xi.temp
        G$theta.xi <- rep(0, p.max)
        G$theta.xi[G$M.xi] <- theta.hat
        tau <- xi.temp - X.M %*% theta.hat
        G$tau.xi <- tau
        G$alpha.xi <- 1/var(tau)[1]
        temp <- optimize(gp.like.lambda,interval=c(0,1e4),G$alpha.xi, tau, D)
        G$lambda.xi <- temp[[1]]
        if(is.null(prior.user$xi$beta.0))
          {
            G$prior$xi$beta.0 <- rep(0, p.max)
          }else{
            G$prior$xi$beta.0 <- prior.user$xi$beta.0
          }
        if(is.null(prior.user$xi$Omega.0))
          {
            G$prior$xi$Omega.0 <- diag(p.max)
          }else{
            G$prior$xi$Omega.0 <- prior.user$xi$Omega.0
          }
        if(is.null(prior.user$xi$alpha.a))
          {
            G$prior$xi$alpha.a <- 1
          }else{
            G$prior$xi$alpha.a <- prior.user$xi$alpha.a
          }
        if(is.null(prior.user$xi$alpha.b))
          {
            G$prior$xi$alpha.b <- 1
          }else{
            G$prior$xi$alpha.b <- prior.user$xi$alpha.b
          }
        if(is.null(prior.user$xi$lambda.a))
          {
            G$prior$xi$lambda.a <- 1
          }else{
            G$prior$xi$lambda.a <- prior.user$xi$lambda.a
          }
        if(is.null(prior.user$xi$lambda.b))
          {
            G$prior$xi$lambda.b <- 1
          }else{
            G$prior$xi$lambda.b <- prior.user$xi$lambda.b
          }
      }
    ##---------------------------------------------------

    return(G)
  }
##--------------------------------------------------------------

##--------------- Routines -------------------------------------
gev.update <- function(G)
  {

    ##----------- Mu Model -------------------------
    Y <- G$X %*% G$theta.mu + G$tau.mu
    if(!G$full) {G$M.mu <- gev.update.M(Y,G$X,G$M.mu,G$alpha.mu,G$lambda.mu,G$D, G$prior$mu$beta.0, G$prior$mu$Omega.0)}
    G$theta.mu <- gev.update.theta(Y,G$X,G$M.mu,G$alpha.mu,G$lambda.mu,G$D, G$prior$mu$beta.0, G$prior$mu$Omega.0) 
    G$tau.mu <- Y - G$X %*% G$theta.mu
    G <- gev.update.tau.mu(G)
    ##----------------------------------------------

    ##---------- Kappa Model -----------------------
    Y <- G$X %*% G$theta.kappa + G$tau.kappa
    if(!G$full) {G$M.kappa <- gev.update.M(Y, G$X, G$M.kappa, G$alpha.kappa, G$lambda.kappa, G$D, G$prior$kappa$beta.0, G$prior$kappa$Omega.0)}
    G$theta.kappa <- gev.update.theta(Y, G$X, G$M.kappa, G$alpha.kappa, G$lambda.kappa, G$D, G$prior$kappa$beta.0, G$prior$kappa$Omega.0)
    G$tau.kappa <- Y - G$X %*% G$theta.kappa
    G <- gev.update.tau.kappa(G)
    ##----------------------------------------------
    
    ##----------- Xi Model -------------------------
    if(G$fixed.xi)
      {
        G$accept.tau.xi <- rep(1,G$n.s) ## A simple hack to make everything downstream play nice
      }else{
        Y <- G$X %*% G$theta.xi + G$tau.xi
        if(!G$full) {G$M.xi <- gev.update.M(Y, G$X, G$M.xi, G$alpha.xi, G$lambda.xi, G$D, G$prior$xi$beta.0, G$prior$xi$Omega.0)}
        G$theta.xi <- gev.update.theta(Y, G$X, G$M.kappa, G$alpha.kappa, G$lambda.kappa, G$D, G$prior$xi$beta.0, G$prior$xi$Omega.0)
        G$tau.xi <- Y - G$X %*% G$theta.xi
        G <- gev.update.tau.xi(G)
      }
    ##----------------------------------------------
    
    G <- gev.update.hyper(G)

    return(G)

  }

gev.results.init <- function(n.s, p, n.reps)
  {
    R <- NULL
    R$THETA <- array(NA, dim = c(n.reps, p, 3))
    R$M <- array(0, dim = c(n.reps, p, 3))
    R$TAU <- array(0, dim = c(n.reps,n.s,3))
    R$LAMBDA <- matrix(0, n.reps,3)
    R$ALPHA <- matrix(0, n.reps,3)
    R$ACCEPT.TAU <- array(0,dim=c(n.reps,n.s,3))
    return(R)
  }

gev.process.results <- function(R, burn=1e2)
  {
    output <- NULL
    S <- dim(R$TAU)[1]
    l <- (burn + 1):S
    tbl.mu <- cbind(colMeans(R$M[l,,1]), colMeans(R$THETA[l,,1]), apply(R$THETA[l,,1],2,"sd"), apply(R$THETA[l,,1],2,"quantile",.025),apply(R$THETA[l,,1],2,"quantile",.975),apply(R$THETA[l,,1],2,"effectiveSize"))
    colnames(tbl.mu) <- c("Prob","Mean","SD","2.5%","97.5%","ESS")
    rownames(tbl.mu) <- paste("V",1:dim(R$THETA)[2],sep="")
    tbl.kappa <- cbind(colMeans(R$M[l,,2]), colMeans(R$THETA[l,,2]), apply(R$THETA[l,,2],2,"sd"), apply(R$THETA[l,,2],2,"quantile",.025),apply(R$THETA[l,,2],2,"quantile",.975),apply(R$THETA[l,,2],2,"effectiveSize"))
    colnames(tbl.kappa) <- c("Prob","Mean","SD","2.5%","97.5%","ESS")
    rownames(tbl.kappa) <- paste("V",1:dim(R$THETA)[2],sep="")
    
    tbl.xi <- cbind(colMeans(R$M[l,,3]), colMeans(R$THETA[l,,3]), apply(R$THETA[l,,3],2,"sd"), apply(R$THETA[l,,3],2,"quantile",.025),apply(R$THETA[l,,3],2,"quantile",.975),apply(R$THETA[l,,3],2,"effectiveSize"))
    colnames(tbl.xi) <- c("Prob","Mean","SD","2.5%","97.5%","ESS")
    rownames(tbl.xi) <- paste("V",1:dim(R$THETA)[2],sep="")

    EF <- t(apply(R$TAU,2,"effectiveSize"))
    tbl.re <- cbind(colMeans(R$LAMBDA),apply(R$LAMBDA,2,"sd"),apply(R$LAMBDA,2,"effectiveSize"),colMeans(R$ALPHA), apply(R$ALPHA,2,"sd"),apply(R$ALPHA,2,"effectiveSize"), apply(EF,2,"min"),colMeans(EF), apply(EF,2,"max"))
    rownames(tbl.re) <- c("mu","kappa","xi")
    colnames(tbl.re) <- c("AlphaMean","AlphaSD","AlphaESS","LambdaMean","LambdaSD","LambdaESS","MinESS","MeanESS","MaxESS")
    output$tbl.mu <- tbl.mu
    output$tbl.kappa <- tbl.kappa
    output$tbl.xi <- tbl.xi
    output$tbl.re <- tbl.re

    D.TAU <- apply(R$ACCEPT.TAU,c(2,3),"mean")
    D.LAMBDA <- R$LAMBDA[l[2]:l[length(l)],] - R$LAMBDA[l[1]:l[length(l) - 1],]
    A.LAMBDA <- (abs(D.LAMBDA) > 1e-5) * 1
    tbl.accept <- cbind(colMeans(A.LAMBDA),apply(D.TAU,2,"min"),apply(D.TAU,2,"mean"),apply(D.TAU,2,"max"))
    colnames(tbl.accept) <- c("Lambda","MinTau","MeanTau","MaxTau")
    rownames(tbl.accept) <- c("mu","kappa","xi")
    output$tbl.accept <- tbl.accept
    return(output)

  }

gev.crps <- function(Y.obs,Y.samp)
  {
    n <- length(Y.obs)
    s <- NULL
    n.samp <- length(Y.samp)
    ss <- sample(1:n.samp,n.samp/2)
    C.2 <- mean( abs(Y.samp[ss] - Y.samp[-ss]))
    for(i in 1:n)
      {
        s[i] <- mean( abs(Y.obs[i] - Y.samp)) - 0.5 * C.2
      }
    return(mean(s))
  }

gev.logscore <- function(Y.obs, Y.samp)
  {
    den <- density(Y.samp)
    den.fun <- splinefun(den$x,den$y)
    ls <- mean(-log(den.fun(Y.obs)))
    return(ls)

  }

gev.impute <- function(R,X.drop, S.drop, burn = NULL, n.each = NULL)
  {
    reps <- dim(R$THETA)[1]
    if(is.null(burn))burn <- round(reps/10)
    if(is.null(n.each))
      {
        n.each <- round(1e6/reps) ##let's get about a million draws
      }
    I <- (burn + 1):reps
    Y <- matrix(0,length(I),n.each)
    S.all <- rbind(S.drop, R$S)
    D.all <- make.D(S.all,S.all)
    for(i in 1:length(I))
      {
        it <- I[i]
        ##------------ Get mu_s --------------
        alpha <- R$ALPHA[it,1]
        lambda <- R$LAMBDA[it,1]
        C <- 1/alpha * exp(-D.all/lambda)
        diag(C) <- diag(C) + 1e-5
        C.inv <- solve(C)
        tau.hat <- sum(-C.inv[1,-1]/C.inv[1,1] * R$TAU[it,,1])
        varsigma <- 1/C.inv[1,1]
        tau.new <- rnorm(1,tau.hat,sd=sqrt(varsigma))
        mu.s <- sum(R$THETA[it,,1] * X.drop) + tau.new
        ##---------------------------------------

        ##------------ Get kappa_s --------------
        alpha <- R$ALPHA[it,2]
        lambda <- R$LAMBDA[it,2]
        C <- 1/alpha * exp(-D.all/lambda)
        diag(C) <- diag(C) + 1e-5
        C.inv <- solve(C)
        tau.hat <- sum(-C.inv[1,-1]/C.inv[1,1] * R$TAU[it,,2])
        varsigma <- 1/C.inv[1,1]
        kappa.hat <- sum(R$THETA[it,,2] * X.drop)
        kappa.s <- rtnorm(1,kappa.hat + tau.hat, sd=sqrt(varsigma),lower=0)
        ##---------------------------------------

        ##------------ Get xi_s --------------
        if(is.null(R$fixed.xi))
          {
            alpha <- R$ALPHA[it,3]
            lambda <- R$LAMBDA[it,3]
            C <- 1/alpha * exp(-D.all/lambda)
            diag(C) <- diag(C) + 1e-5
            C.inv <- solve(C)
            tau.hat <- sum(-C.inv[1,-1]/C.inv[1,1] * R$TAU[it,,3])
            varsigma <- 1/C.inv[1,1]
            tau.new <- rnorm(1,tau.hat,sd=sqrt(varsigma))
            xi.s <- sum(R$THETA[it,,3] * X.drop) + tau.new
          }else{
            xi.s <- R$fixed.xi
          }
        ##---------------------------------------

        Y[i,] <- rgev(n.each,mu.s,1/kappa.s,xi.s)
      }
    return(as.vector(Y))
  }

spatial.gev.bma <- function(Y.list, X.all,S,n.reps,prior.user= NULL, full = FALSE, fixed.xi = NULL, print.every=0)
  {
    G <- gev.init(Y.list,X.all,S, prior.user,full,fixed.xi)
    R <- gev.results.init(length(Y.list), dim(X.all)[2], n.reps)
    R$S <- S
    R$fixed.xi <- fixed.xi
    for(i in 1:n.reps)
      {
        if( (print.every >0) && (i %% print.every == 0))print(paste("On Interation", i))
        G <- gev.update(G)
        R$THETA[i,,1] <- G$theta.mu
        R$THETA[i,,2] <- G$theta.kappa
        R$THETA[i,,3] <- G$theta.xi
        R$M[i,G$M.mu,1] <- 1
        R$M[i,G$M.kappa,2] <- 1
        R$M[i,G$M.xi,3] <- 1
        R$TAU[i,,1] <- G$tau.mu
        R$TAU[i,,2] <- G$tau.kappa
        R$TAU[i,,3] <- G$tau.xi
        R$LAMBDA[i,1] <- G$lambda.mu
        R$LAMBDA[i,2] <- G$lambda.kappa
        R$LAMBDA[i,3] <- G$lambda.xi
        R$ALPHA[i,1] <- G$alpha.mu
        R$ALPHA[i,2] <- G$alpha.kappa
        R$ALPHA[i,3] <- G$alpha.xi
        R$ACCEPT.TAU[i,,] <- cbind(G$accept.tau.mu,G$accept.tau.kappa,G$accept.tau.xi)
      }
    return(R)

  }
