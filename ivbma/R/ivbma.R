get.Y.tilde <- function(Y,Sigma,eta)
  {
    dd <- dim(Sigma)[2]
    not.y <- (1:dd)[-1]
    Y.tilde <- Y - eta %*% solve(Sigma[not.y,not.y,drop=FALSE]) %*% Sigma[not.y,1]
    return(Y.tilde)
  }

get.xi <- function(Sigma)
  {
    K <- solve(Sigma)
    return(1/K[1,1])
  }

get.YX.star <- function(Y,X,W,rho,Sigma,eta,r,j)
  {

    gamma <- rho[-(1:r)]
    if(r == 1)
      {
        Y.star <- (Y - W %*% gamma) / rho[1]
        X.star <- X
      }else{
        both <- c(1, j + 1)
        beta <- rho[1:r]
        Psi <- Sigma[both,-both] %*% solve(Sigma[-both,-both])
        Y.star <- (Y- W %*% gamma - X[,-j,drop=FALSE] %*% beta[-j] - eta[,-j,drop=FALSE] %*% Psi[1,]) / beta[j]
        X.star <- X[,j] - eta[,-j,drop=FALSE] %*% Psi[2,]
      }
    return(cbind(Y.star,X.star))
  }

get.Sigma.star <- function(Sigma,r,j)
  {
    both <- c(1, j + 1)
    Sigma.star <- Sigma[both,both] - Sigma[both,-both] %*% solve(Sigma[-both,-both]) %*% Sigma[-both,both]
    return(Sigma.star)
  }

get.Omega <- function(Sigma, beta)
  {
    Omega <- matrix(c(Sigma[2,2] + Sigma[1,1] / beta^2 + 2 * Sigma[2,1] / beta,
                      Sigma[2,2] + (1 / beta) * Sigma[1,2],
                      Sigma[2,2] + (1 / beta) * Sigma[1,2],
                      Sigma[2,2]),2,2)
    return(Omega)
  }

get.X.star <- function(Sigma,X,eps,eta,r,j)
  {
    if(r == 1)
      {
        X.star <- X - Sigma[2,1] / Sigma[1,1] * eps
      }else{
        X.star <- X[,j] - cbind(eps, eta[,-j]) %*% solve(Sigma[-(j + 1),-(j + 1)]) %*% Sigma[-(j + 1),(j + 1)]  
      }
    return(X.star)
  }

get.omega <- function(Sigma,j)
  {
    K <- solve(Sigma)
    return(1/K[j + 1, j + 1])
  }

update.L <- function(Y.tilde,xi, L, V)
{

  p.V <- dim(V)[2]
  
  ##----- Flip a variable -----
  w <- sample(1:p.V,1)
  L.new <- L
  L.new[w] <- 1 - L.new[w]
  if(sum(L.new) == 0)
    {
      p.1 <- sum(L)
      V.1 <- V[, (1:p.V)[L==1],drop=FALSE]
      Xi.1 <- diag(p.1) + (t(V.1)%*%V.1) / xi
      rho.hat.1 <- xi^(-1) %*% t(Y.tilde) %*% V.1 %*% solve(Xi.1)
      rho <- rep(0,p.V)
      rho[(1:p.V)[L==1]] <- rmvnorm.precision(rho.hat.1, Xi.1)
      return(list(L = L, rho = rho))
    }
  ##-----------------------------

  ##---- Stats for old model ----
  p.1 <- sum(L)
  V.1 <- V[, (1:p.V)[L==1],drop=FALSE]
  Xi.1 <- diag(p.1) + (t(V.1)%*%V.1) / xi
  rho.hat.1 <- xi^(-1) %*% t(Y.tilde) %*% V.1 %*% solve(Xi.1)
  score.1 <- 0.5 * rho.hat.1 %*% Xi.1 %*% t(rho.hat.1) - 0.5 * mydet(Xi.1)
  ##------------------------------

  ##---- Stats for new model -----
  p.2 <- sum(L.new)
  V.2 <- V[, (1:p.V)[L.new==1],drop=FALSE]
  Xi.2 <- diag(p.2) + (t(V.2)%*%V.2) / xi
  rho.hat.2 <- xi^(-1) %*% t(Y.tilde) %*% V.2 %*% solve(Xi.2)
  score.2 <- 0.5 * rho.hat.2 %*% Xi.2 %*% t(rho.hat.2) - 0.5 * mydet(Xi.2)
  ##------------------------------

  ##----- Make a Decision --------
  alpha <- score.2 - score.1
  if(log(runif(1)) < alpha)
    {
      L <- L.new
      rho.hat.1 <- rho.hat.2
      Xi.1 <- Xi.2
    }
  ##-----------------------------

  rho <- rep(0,p.V)
  rho[(1:p.V)[L==1]] <- rmvnorm.precision(rho.hat.1, Xi.1)

  return(list(L = L, rho = rho))

}

update.rho <- function(Y.tilde,xi,L,V)
  {
    p.V <- dim(V)[2]
    rho <- rep(0, p.V)
    V.L <- V[, (1:p.V)[L==1],drop=FALSE]
    Xi.L <- diag(sum(L)) + (t(V.L)%*%V.L) / xi
    rho.hat.L <- xi^(-1) %*% t(Y.tilde) %*% V.L %*% solve(Xi.L)
    rho[(1:p.V)[L==1]] <- rmvnorm.precision(rho.hat.L, Xi.L)
    return(rho)
  }

update.M <- function(YX.star, U, M, Omega)
  {
    ##------ Keep book --------
    n <- dim(YX.star)[1]
    p.U <- dim(U)[2]
    ##------------------------

    ##----- Flip a variable -----
    w <- sample(1:p.U,1)
    M.new <- M
    M.new[w] <- 1 - M.new[w]
    if(sum(M.new) == 0)
      {
        p.1 <- sum(M)
        Phi.inv <- solve(chol(Omega))
        S.row <- YX.star %*% Phi.inv
        S <- as.vector(S.row)
        U.1 <- U[,(1:p.U)[M == 1],drop=FALSE]
        U.1.row <- cbind(as.vector(U.1),as.vector(U.1)) %*% Phi.inv
        T.1 <- rbind(matrix(U.1.row[,1],n,p.1),matrix(U.1.row[,2],n,p.1))
        Omega.1 <- diag(p.1) + t(T.1) %*% T.1
        lambda.1 <- t(S) %*% T.1 %*% solve(Omega.1)
        lambda <- rep(0,p.U)
        lambda[ (1:p.U)[M==1]] <- rmvnorm.precision(lambda.1,Omega.1)
        return(list(M = M, lambda = lambda))
      }
    ##-----------------------------
    
    ##------ Dependent Var ---
    Phi.inv <- solve(chol(Omega))
    S.row <- YX.star %*% Phi.inv
    S <- as.vector(S.row)
    ##-------------------------

    ##---- New Score ------------
    p.2 <- sum(M.new)
    U.2 <- U[,(1:p.U)[M.new == 1],drop=FALSE]
    U.2.row <- cbind(as.vector(U.2),as.vector(U.2)) %*% Phi.inv
    T.2 <- rbind(matrix(U.2.row[,1],n,p.2),matrix(U.2.row[,2],n,p.2))
    Omega.2 <- diag(p.2) + t(T.2) %*% T.2
    lambda.2 <- t(S) %*% T.2 %*% solve(Omega.2)
    score.2 <- 0.5 * lambda.2 %*% Omega.2 %*% t(lambda.2) - 0.5 * mydet(Omega.2)
    ##---------------------------

    ##---- Old Score ------------
    p.1 <- sum(M)
    U.1 <- U[,(1:p.U)[M == 1],drop=FALSE]
    U.1.row <- cbind(as.vector(U.1),as.vector(U.1)) %*% Phi.inv
    T.1 <- rbind(matrix(U.1.row[,1],n,p.1),matrix(U.1.row[,2],n,p.1))
    Omega.1 <- diag(p.1) + t(T.1) %*% T.1
    lambda.1 <- t(S) %*% T.1 %*% solve(Omega.1)
    score.1 <- 0.5 * lambda.1 %*% Omega.1 %*% t(lambda.1) - 0.5 * mydet(Omega.1)
    ##---------------------------

    ##---- Decide ---------------
    alpha <- score.2 - score.1
    if(log(runif(1)) < alpha)
      {
        M <- M.new
        Omega.1 <- Omega.2
        lambda.1 <- lambda.2
      }
    ##---------------------------

    lambda <- rep(0,p.U)
    lambda[ (1:p.U)[M==1]] <- rmvnorm.precision(lambda.1,Omega.1)

    return(list(M = M, lambda = lambda))
  }

update.lambda <- function(YX.star,U,M,Omega)
  {
    ##------ Keep book --------
    n <- dim(YX.star)[1]
    p.U <- dim(U)[2]
    p.M <- sum(M)
    ##------------------------

    ##------ Dependent Var ---
    Phi.inv <- solve(chol(Omega))
    S.row <- YX.star %*% Phi.inv
    S <- as.vector(S.row)
    ##-------------------------

    ##---- New Score -----------
    U.M <- U[,(1:p.U)[M == 1],drop=FALSE]
    U.M.row <- cbind(as.vector(U.M),as.vector(U.M)) %*% Phi.inv
    T.M <- rbind(matrix(U.M.row[,1],n,p.M),matrix(U.M.row[,2],n,p.M))
    Omega.M <- diag(p.M) + t(T.M) %*% T.M
    lambda.M <- t(S) %*% T.M %*% solve(Omega.M)
    ##---------------------------

    lambda <- rep(0,p.U)
    lambda[ (1:p.U)[M==1]] <- rmvnorm.precision(lambda.M,Omega.M)

    return(lambda)
  }


update.M.SUR <-function(X.star, U, omega, M)
  {

    n <- length(X.star)
    p.U <- dim(U)[2]
    p.1 <- sum(M)

    ##----- Propose ---------
    w <- sample(1:p.U,1)
    M.new <- M
    M.new[w] <- 1 - M.new[w]
    p.2 <- sum(M.new)
    if(p.2 == 0)
      {
        U.1 <- U[,(1:p.U)[M == 1],drop=FALSE]
        Omega.1 <- diag(p.1) + (t(U.1) %*% U.1) / omega
        lambda.1 <- omega^(-1) * t(X.star) %*% U.1 %*% solve(Omega.1)
        lambda <- rep(0,p.U)
        lambda[ (1:p.U)[M==1]] <- rmvnorm.precision(lambda.1,Omega.1)
        
        return(list(M = M, lambda = lambda))
      }
    ##------------------------

    ##----- Score Old --------
    U.1 <- U[,(1:p.U)[M == 1],drop=FALSE]
    Omega.1 <- diag(p.1) + (t(U.1) %*% U.1) / omega
    lambda.1 <- omega^(-1) * t(X.star) %*% U.1 %*% solve(Omega.1)
    score.1 <- 0.5 * lambda.1 %*% Omega.1 %*% t(lambda.1) - 0.5 * mydet(Omega.1)
    ##----------------------------

    ##------ Score New -----------
    U.2 <- U[,(1:p.U)[M.new == 1],drop=FALSE]
    Omega.2 <- diag(p.2) + (t(U.2) %*% U.2) / omega
    lambda.2 <- omega^(-1) * t(X.star) %*% U.2 %*% solve(Omega.2)
    score.2 <- 0.5 * lambda.2 %*% Omega.2 %*% t(lambda.2) - 0.5 * mydet(Omega.2)
    ##----------------------------

    ##------- Decide -------------
    alpha <- score.2 - score.1
    if(log(runif(1)) < alpha)
      {
        M <- M.new
        Omega.1 <- Omega.2
        lambda.1 <- lambda.2
      }
    ##---------------------------

    lambda <- rep(0,p.U)
    lambda[ (1:p.U)[M==1]] <- rmvnorm.precision(lambda.1,Omega.1)

    return(list(M = M, lambda = lambda))
    
  }

update.lambda.SUR <- function(X.star, U, omega, M)
{
  p.U <- dim(U)[2]
  p.M <- sum(M)
  U.M <- U[,(1:p.U)[M == 1],drop=FALSE]
  Omega.M <- diag(p.M) + (t(U.M) %*% U.M) / omega
  lambda.M <- omega^(-1) * t(X.star) %*% U.M %*% solve(Omega.M)
  lambda <- rep(0, p.U)
  lambda[(1:p.U)[M == 1]] <- rmvnorm.precision(lambda.M, Omega.M)
  return(lambda)
}

update.Sigma <- function(eps,eta,r)
  {
    E <- cbind(eps,eta)
    UU <- t(E) %*% E
    n <- dim(E)[1]
    
    Sigma <- rinvwish(3 + n,diag(r +1) + UU)

  }

ivbma.sample.theta <- function(theta,D,full)
  {

    ##---- Update Outcome Parameters ---------
    Y.tilde <- get.Y.tilde(D$Y,theta$Sigma,theta$eta)
    xi <- get.xi(theta$Sigma)
    if(!full)
      {
        l <- update.L(Y.tilde,xi,theta$L,D$V)
        theta$L <- l$L
        theta$rho <- l$rho
      }else{
        theta$rho <- update.rho(Y.tilde,xi,theta$L,D$V)
      }
    theta$eps <- D$Y - D$V %*% theta$rho
    ##----------------------------------------

    ##----- Update Instrument Parameters -----
    for(i in 1:D$r)
      {
        if(theta$L[i])
          {
            YX.star <- get.YX.star(D$Y,D$X,D$W,theta$rho,theta$Sigma,theta$eta,D$r,i)
            if(D$r > 1)
              {
                Sigma.star <- get.Sigma.star(theta$Sigma,D$r,i)
                Omega <- get.Omega(Sigma.star,theta$rho[i])
              }else{
                Omega <- get.Omega(theta$Sigma,theta$rho[i])
              }
            if(!full)
              {
                l <- update.M(YX.star,D$U,theta$M[,i],Omega)
                theta$M[,i] <- l$M
                theta$lambda[,i] <- l$lambda
              }else{
                theta$lambda[,i] <- update.lambda(YX.star, D$U, theta$M[,i],Omega)
              }
          }else{
            X.star <- get.X.star(theta$Sigma,D$X,theta$eps,theta$eta,D$r,i)
            omega <- get.omega(theta$Sigma,i)
            if(!full)
              {
                l <- update.M.SUR(X.star,D$U,omega, theta$M[,i])
                theta$M[,i] <- l$M
                theta$lambda[,i] <- l$lambda
              }else{
                theta$lambda[,i] <- update.lambda.SUR(X.star,D$U, omega, theta$M[,i])
              }
          }
        theta$eta[,i] <- D$X[,i] - D$U %*% theta$lambda[,i]
      }
    ##------------------------------------------

    ##------- Update Sigma ---------------------
    theta$Sigma <- update.Sigma(theta$eps,theta$eta,D$r)
    ##------------------------------------------
    
    return(theta)
  }

ivbma.diagnostics <- function(theta,D)
  {
    ##EVERYTHING HERE IS COMPLETELY EXPERIMENTAL
    ##It does NOT constitute peer reviewed methodology
    n <- dim(D$Z)[1]
    q <- dim(D$Z)[2]
    r <- D$r
    if(r == 1)
      {
        keep <- (1:q)[theta$M[1:q,1] > 0.5]
      }else{
        keep <- (1:q)[rowSums(theta$M[1:q,]) > 0.5]
      }
    k.2 <- length(keep)
    if(k.2 == 0)return(1)
    Hat <- D$Z %*% solve(t(D$Z) %*% D$Z) %*% t(D$Z) %*% theta$eps
    r2 <- sum( (Hat - theta$eps) ^2)/sum(theta$eps^2)
    Sargan <- 1 - pchisq(n * r2, k.2 - 1)

    ##----------------------------------------------
    ##The Following is a proto-type of a Bayesian Sargan Test
    ## "A Bayesian Diagnostic of Instrument Validity via Conditional Bayes Factors"
    ## T. Kaiser and A. Lenkoski, in progress
    ## If you're going to use this don't be a jerk and not cite us.
    K.tilde <- (t(D$Z[,keep,drop=FALSE]) %*% D$Z[,keep,drop=FALSE] + diag(k.2)) / theta$Sigma[1,1]
    b.tilde <- solve(K.tilde) %*% t(D$Z[,keep,drop=FALSE]) %*% theta$eps
    score.top <- 0.5 * t(b.tilde) %*% K.tilde %*% b.tilde - 0.5 * mydet(K.tilde)
    Bayesian.Sargan <- 1 - 1 / (1 + exp(-score.top))
    ## Seems to be right, at least.
    ##--------------------------------------------------------

    return(c(Sargan, Bayesian.Sargan))
  }

ivbma.init <- function(D,full)
  {

    n <- length(D$Y)
    q <- dim(D$Z)[2]
    p <- dim(D$W)[2]
    r <- dim(D$X)[2]
    d <- 3  ##Prior degrees of freedom
    R <- diag(r + 1) ##Prior covariance

    ##---------- Initialization --------------------
    p.V <- p + r
    p.U <- p + q
    rho <- NULL
    lambda <- matrix(NA,p + q,r) 
    eta <- matrix(0,n,r)
    M <- matrix(0,p + q,r)
    L <- NULL
    for(i in 1:r)
      {
        lambda[,i] <- solve(t(D$U) %*% D$U + diag(p.U)) %*% t(D$U) %*% D$X[,i]
        if(full)
          {
            M[,i] <- 1
          }else{
            M[,i] <- rbinom(p + q,1,.5)
          }
        lambda[ (1:p.U)[M[,i] == 0] ,i] <- 0
        eta[,i] <- D$X[,i] - D$U %*% lambda[,i]
      }
    eps <- D$Y - D$V %*% (solve(t( D$V) %*% D$V + diag(p.V)) %*% t(D$V) %*% D$Y)
    Sigma <- rinvwish(d + n,R + (n-1) * cov(cbind(eps, eta)))
    if(full)
      {
        L <- rep(1, p.V)
      }else{
        L <- rbinom(p.V,1,.5)
      }
    ##-----------------------------------------------

    ##--------- Load Up -----------------------------
    theta <- NULL
    theta$eps <- eps
    theta$eta <- eta
    theta$Sigma <- Sigma
    theta$lambda <- lambda
    theta$rho <- rho
    theta$M <- M
    theta$L <- L
    ##-----------------------------------------------

    return(theta)
  }



ivbma.results.init <- function(D,odens, run.diagnostics)
{
  p.V <- dim(D$V)[2]
  p.U <- dim(D$U)[2]
  r <- D$r
  ##-------- Information to be returned ----------
  results <- NULL
  results$rho  <- matrix(0, odens, p.V)
  results$rho.bar <- rep(0,p.V)
  results$lambda <- array(0, dim = c(p.U,r, odens))
  results$lambda.bar <- rep(0, p.U)
  results$Sigma <- array(dim=c(r + 1,r + 1,odens))
  results$Sigma.bar <- matrix(0,r + 1, r + 1)
  results$M <- array(dim = c(p.U,r, odens))
  results$M.bar <- matrix(0, p.U, r)
  results$L <- matrix(0,odens,p.V)
  results$L.bar <- rep(0, p.V)
  results$run.diagnostics <- FALSE
  if(run.diagnostics)
    {
      results$run.diagnostics <- TRUE
      results$Sargan <- 0
      results$Bayesian.Sargan <- 0
    }
  ##----------------------------------------------

  return(results)
}

ivbma <- function(Y,X,Z,W,s=1e3,b = round(s/10),
                  full = FALSE,odens = min(c(5e3,s-b)),
                  print.every = round(s/10),
                  run.diagnostics = FALSE)
  {

    ##----------- Fill Data ---------
    D <- NULL
    D$Y <- Y
    D$U <- cbind(as.matrix(Z),as.matrix(W))
    D$Z <- as.matrix(Z)
    D$V <- cbind(as.matrix(X),as.matrix(W))
    D$W <- as.matrix(W)
    D$X <- as.matrix(X)
    D$r <- dim(D$X)[2]
    ##-------------------------------

    ##--------- Initialize ----------
    theta <- ivbma.init(D,full)
    results <- ivbma.results.init(D,odens, run.diagnostics)
    which.save <- round(seq(b + 1, s, length = odens))
    save.loc <- 1
    next.save <- which.save[save.loc]
    ##--------------------------------
    print(paste("Running IVBMA for",s,"iterations",Sys.time()))
    for(i in 1:s)
      {
        if(i %% print.every == 0)print(paste("On Iteration", i,Sys.time()))

        theta <- ivbma.sample.theta(theta,D,full)

        
        ##---------- Record ----------------------------
        if(i == next.save)
          {
            ##--------------- Save -----------------------------------
            results$rho[save.loc,] <- theta$rho
            results$lambda[,,save.loc] <- theta$lambda
            results$Sigma[,,save.loc] <- theta$Sigma
            results$L[save.loc,] <- theta$L
            results$M[,,save.loc] <- theta$M
            save.loc <- save.loc + 1
            next.save <- which.save[save.loc]
            ##--------------------------------------------------------
            if(run.diagnostics)
              {
                dd <- ivbma.diagnostics(theta,D)
                results$Sargan <- results$Sargan + dd[1] / odens
##                results$Bayesian.Sargan[i] <- results$Bayesian.Sargan + dd[2] / odens
                results$Bayesian.Sargan[i] <- dd[2]
              }
          }
        if(i > b)
          {
            results$rho.bar <- results$rho.bar + theta$rho/ (s - b)
            results$lambda.bar <- results$lambda.bar + theta$lambda / (s - b)
            results$L.bar <- results$L.bar + theta$L / (s - b)
            results$M.bar <- results$M.bar + theta$M / (s - b)
            results$Sigma.bar <- results$Sigma.bar + theta$Sigma / (s - b)
          }
        ##------------------------------------------------------
      }
    class(results) <- c("ivbma",class(results))
    return(results)
  }

