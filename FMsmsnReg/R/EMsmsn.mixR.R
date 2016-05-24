EMsmsn.mixR <- function(y, x1, Abetas = NULL, medj= NULL, sigma2 = NULL, shape = NULL, pii = NULL, g = NULL, get.init = TRUE, criteria = TRUE, group = FALSE,
                        family = "Skew.t", error = 0.00001, iter.max = 500, obs.prob= FALSE, kmeans.param = NULL)
{
  if(get.init == TRUE)
  {
    if(length(g) == 0) stop("g is not specified correctly.\n")

    k.iter.max <- 50
    n.start    <- 1
    algorithm  <- "Hartigan-Wong"
    if(length(kmeans.param) > 0)
    {
      if(length(kmeans.param$iter.max)  > 0) {k.iter.max <- kmeans.param$iter.max }
      if(length(kmeans.param$n.start)   > 0) {n.start    <- kmeans.param$n.start  }
      if(length(kmeans.param$algorithm) > 0) {algorithm  <- kmeans.param$algorithm}
    }

    if(g > 1)
    {
      Abetas <- solve(t(x1)%*%x1)%*%t(x1)%*%y
      yr     <- y-x1%*%Abetas
      init   <- kmeans(yr,g,k.iter.max,n.start,algorithm)
      pii    <- init$size/length(yr)
      medj   <- as.vector(init$centers)
      sigma2 <- init$withinss/init$size
      shape  <- c()
      if(family == "Skew.cn")
      {
        nu    <- c(0.1,0.1)
      }else{
        nu    <- 3
      }
      for (j in 1:g)
      {
        m3       <- (1/init$size[j])*sum( (yr[init$cluster == j] - medj[j])^3)
        shape[j] <- sign(m3/sigma2[j]^(3/2))
      }
    }else{
      Abetas <- solve(t(x1)%*%x1)%*%t(x1)%*%y
      yr     <- y-x1%*%Abetas
      medj   <- mean(yr)
      sigma2 <- var(yr)
      pii    <- 1
      m3     <- (1/length(yr))*sum( (yr - medj)^3 )
      shape  <- sign(m3/sigma2^(3/2))
      if(family == "Skew.cn")
      {
        nu    <- c(0.1,0.1)
      }else{
        nu     <- 3
      }
    }
  }

  ################################################################################
  #                               Begin Family
  ################################################################################

  ################################################################################
  #                              (1)   Skew-t
  ################################################################################

  if (family == "Skew.t")
  {
    start.time  <- Sys.time() #Begin Time
    q           <- ncol(x1)
    tetam       <- matrix(data=NA,nrow=4*g + q + 1,ncol=iter.max)
    #Parameters for "optimize" function
    TOLERANCIA  <- 1e-6
    MAX_NU      <- 20
    MIN_NU      <- 2.01

    k1          <- sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
    b           <- -sqrt(2/pi)*k1
    n           <- length(y)
    p           <- ncol(x1)-1
    x           <- as.matrix(x1[,2:(p+1)])

    beta0       <- Abetas[1]                  ## slope
    betas       <- as.matrix(Abetas[2:(p+1)]) ## parameters of regression

    delta       <- Delta <- Gama <- varphi<-rep(0,g)
    mu          <- matrix(0,n,g)

    for (k in 1:g)
    {
      delta[k]   <- shape[k]/(sqrt(1 + shape[k]^2))
      Delta[k]   <- sqrt(sigma2[k])*delta[k]
      Gama[k]    <- sigma2[k] - Delta[k]^2
      varphi[k]  <- beta0+medj[k]
      mu[,k]     <- x%*%betas+varphi[k]+b*Delta[k]
    }

    teta        <- c(Abetas, medj, Delta, Gama, pii, nu)

    criterio    <- 1
    count       <- 0

    lk          <- sum(log(d.mixedST(y, pii, mu , sigma2, shape, nu)))## log-likelihood

    while((criterio > error) && (count <= iter.max))
    {#begin while
      count      <- count + 1
      tal        <- matrix(0, n, g)
      S1         <- matrix(0, n, g)
      S2         <- matrix(0, n, g)
      S3         <- matrix(0, n, g)
      soma1      <- matrix(0, p, p)
      soma2      <- matrix(0, p, 1)

      for (j in 1:g)
      {
        ### E-step: computing ui, tui, tui2 ###
        ### mu[,j]<-x%*%betas+b*Delta[j]+varphi[j]
        dj         <- ((y - mu[,j])/sqrt(sigma2[j]))^2
        Mtij2      <- 1/(1 + (Delta[j]^2)*(Gama[j]^(-1)))
        Mtij       <- sqrt(Mtij2)
        mutij      <- Mtij2*Delta[j]*(Gama[j]^(-1))*(y - mu[,j])+b
        A          <- (mutij-b) / Mtij

        E          <- (2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2[j])*dt.ls(y, mu[,j], sigma2[j],shape[j] ,nu))
        u          <- ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2[j])*dt.ls(y, mu[,j], sigma2[j],shape[j] ,nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)

        d1         <- dt.ls(y, mu[,j], sigma2[j], shape[j], nu)
        if(length(which(d1 == 0)) > 0) {d1[which(d1 == 0)] <- .Machine$double.xmin}
        d2         <- d.mixedST(y, pii, mu, sigma2, shape, nu)
        if(length(which(d2 == 0)) > 0) {d2[which(d2 == 0)] <- .Machine$double.xmin}

        tal[,j]    <- d1*pii[j] / d2
        S1[,j]     <- tal[,j]*u
        S2[,j]     <- tal[,j]*(mutij*u + Mtij*E)
        S3[,j]     <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*(mutij+b)*E)

        ### M-step: atualizar mu, Delta, Gama, sigma2, shape ###
        pii[j]     <- (1/n)*sum(tal[,j])
        Delta[j]   <- sum(S2[,j]*(y - mu[,j] + b*Delta[j]))/sum(S3[,j])
        Gama[j]    <- sum(S1[,j]*(y - mu[,j] + b*Delta[j])^2 - 2*(y - mu[,j] + b*Delta[j])*Delta[j]*S2[,j] + Delta[j]^2*S3[,j]) / sum(tal[,j])
        sigma2[j]  <- Gama[j] + Delta[j]^2
        shape[j]   <- sigma2[j]^(-1/2)*Delta[j]/sqrt(1 - Delta[j]^2*sigma2[j]^(-1))
        varphi[j]  <- sum(S1[,j]*(y - x%*%betas) - Delta[j]*S2[,j])/sum(S1[,j])
        soma1      <- soma1 + t(x)%*%diag(S1[,j]/Gama[j])%*%x
        soma2      <- soma2 + t(x)%*%diag(tal[,j]/Gama[j])%*%(u*(y - varphi[j]) - (mutij*u + Mtij*E)*Delta[j])
      }

      betas       <- solve(soma1)%*%soma2
      pii[g]      <- 1 - (sum(pii) - pii[g])

      zero.pos    <- NULL
      zero.pos    <- which(pii == 0)
      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }

      beta0       <- sum(pii*varphi)
      Abetas      <- c(beta0,betas)

      for (k in 1:g)
      {
        mu[,k]    <- x%*%betas+varphi[k]+b*Delta[k]
        medj[k]   <- varphi[k]-beta0
      }

      if (pii[2] < 0.5 && g==2)
      {
        mu       <- cbind(mu[,2],mu[,1])
        medj     <- as.vector(c(medj[2],medj[1]))
        varphi   <- as.vector(c(varphi[2], varphi[1]))
        pii      <- as.vector(c(pii[2], pii[1]))
        sigma2   <- as.vector(c(sigma2[2], sigma2[1]))
      }

      #Here begining the update for estimating the value of nu
      logvero.ST  <- function(nu) sum(log(d.mixedST(y, pii, mu, sigma2, shape, nu)))
      nu          <- optimize(f=logvero.ST, interval=c(MIN_NU,MAX_NU),maximum=TRUE,tol=TOLERANCIA)$maximum

      teta        <- c(Abetas, medj, Delta, Gama, pii, nu)
      tetam[,count] <- c(Abetas, medj, sigma2, shape, pii, nu)
      auxlog      <- d.mixedST(y, pii, mu, sigma2,shape,nu)

      if(length(which(auxlog == 0)) > 0) auxlog[which(auxlog == 0)] <- .Machine$double.xmin
      lk1         <- sum(log(auxlog))
      criterio    <- abs(lk1/lk-1)
      lk          <- lk1
    } # fim while

    if(criteria == TRUE)
    {
      cl  <- apply(tal, 1, which.max)
      icl <- 0
      for (j in 1:g) {icl <- icl+sum(log(pii[j]*dt.ls(y, mu[,j], sigma2[j], shape[j], nu)))}
    }

    tetam           <- tetam[,1:count]
    EPim            <- im.FMsmsnReg(y, x1, g, model = NULL, initial.model=FALSE,family,Abetas, medj, sigma2=sigma2, shape = shape, pii,nu)
    EP              <- EPim$EP
    end.time        <- Sys.time()
    time.taken      <- end.time - start.time
  }# end skew.t

  ################################################################################
  ##                              (2)     Skew-cn
  ################################################################################

  if (family == "Skew.cn"){
    start.time  <- Sys.time()
    q           <- ncol(x1)
    tetam       <- matrix(data=NA,nrow=4*g + q + 2,ncol=iter.max)
    if(length(nu) != 2) stop("The Skew.cn need a vector of two components in nu both between (0,1).")
    if((nu[1] <= 0) || (nu[1] >= 1) || (nu[2] <= 0) || (nu[2] >= 1)) stop("nu component(s) out of range.")

    k1           <- nu[1]/nu[2]^(1/2)+1-nu[1]
    b            <- -sqrt(2/pi)*k1
    n            <- length(y)
    p            <- ncol(x1)-1
    x            <- as.matrix(x1[,2:(p+1)])

    beta0        <- Abetas[1]                  ## slope
    betas        <- as.matrix(Abetas[2:(p+1)]) ## parameters of regression

    delta        <- Delta <- Gama <- varphi<-rep(0,g)
    mu           <- matrix(0,n,g)

    for (k in 1:g)
    {
      delta[k]    <- shape[k]/(sqrt(1 + shape[k]^2))
      Delta[k]    <- sqrt(sigma2[k])*delta[k]
      Gama[k]     <- sigma2[k] - Delta[k]^2
      varphi[k]   <- beta0+medj[k]
      mu[,k]      <- x%*%betas+varphi[k]+b*Delta[k]
    }

    teta         <- c(Abetas, medj, Delta, Gama, pii, nu)

    criterio     <- 1
    count        <- 0
    lk           <- sum(log(d.mixedSNC(y, pii, mu , sigma2, shape, nu)))## log-likelihood

    while((criterio > error) && (count <= iter.max))
    { #begin while
      count <- count + 1
      tal        <- matrix(0, n, g)
      S1         <- matrix(0, n, g)
      S2         <- matrix(0, n, g)
      S3         <- matrix(0, n, g)
      soma1      <- matrix(0, p, p)
      soma2      <- matrix(0, p, 1)

      for (j in 1:g)
      {
        ### E-step: computing ui, tui, tui2 ###
        ### mu[,j]<-x%*%betas+b*Delta[j]+varphi[j]
        dj       <- ((y - mu[,j])/sqrt(sigma2[j]))^2
        Mtij2    <- 1/(1 + (Delta[j]^2)*(Gama[j]^(-1)))
        Mtij     <- sqrt(Mtij2)
        mutij    <- Mtij2*Delta[j]*(Gama[j]^(-1))*(y - mu[,j]) + b
        A        <- (mutij - b) / Mtij

        u        <- (2/dSNC(y,mu[,j],sigma2[j],shape[j],nu))*(nu[1]*nu[2]*dnorm(y,mu[,j],sqrt(sigma2[j]/nu[2]))*pnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu[,j],sqrt(sigma2[j]))*pnorm(A,0,1))
        E        <- (2/dSNC(y,mu[,j],sigma2[j],shape[j],nu))*(nu[1]*sqrt(nu[2])*dnorm(y,mu[,j],sqrt(sigma2[j]/nu[2]))*dnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu[,j],sqrt(sigma2[j]))*dnorm(A,0,1))

        d1       <- dSNC(y, mu[,j], sigma2[j], shape[j], nu)
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
        d2       <- d.mixedSNC(y, pii, mu, sigma2, shape, nu)
        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

        tal[,j]  <- d1*pii[j]/d2
        S1[,j]   <- tal[,j]*u
        S2[,j]   <- tal[,j]*(mutij*u + Mtij*E)
        S3[,j]   <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*(mutij+b)*E)

        ### M-step: atualizar mu, Delta, Gama, sigma2 ###
        pii[j]    <- (1/n)*sum(tal[,j])
        Delta[j]  <- sum(S2[,j]*(y - mu[,j] + b*Delta[j]))/sum(S3[,j])
        Gama[j]   <- sum(S1[,j]*(y - mu[,j] + b*Delta[j])^2 - 2*(y - mu[,j] + b*Delta[j])*Delta[j]*S2[,j] + Delta[j]^2*S3[,j]) / sum(tal[,j])
        sigma2[j] <- Gama[j] + Delta[j]^2
        shape[j]  <- sigma2[j]^(-1/2)*Delta[j]/sqrt(1 - Delta[j]^2*sigma2[j]^(-1))
        varphi[j] <- sum(S1[,j]*(y - x%*%betas) - Delta[j]*S2[,j])/sum(S1[,j])
        soma1     <- soma1 + t(x)%*%diag(S1[,j]/Gama[j])%*%x
        soma2     <- soma2 + t(x)%*%diag(tal[,j]/Gama[j])%*%(u*(y - varphi[j]) - (mutij*u + Mtij*E)*Delta[j])
      }

      betas      <- solve(soma1)%*%soma2
      pii[g]     <- 1 - (sum(pii) - pii[g])

      zero.pos   <- NULL
      zero.pos   <- which(pii == 0)
      if(length(zero.pos) != 0){
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }

      beta0      <- sum(pii*varphi)
      Abetas     <- c(beta0,betas)

      for (k in 1:g)
      {
        mu[,k]   <- x%*%betas+varphi[k]+b*Delta[k]
        medj[k]  <- varphi[k]-beta0
      }

      if (pii[2] < 0.5  && g==2)
      {
        nu         <- cbind(nu[2],nu[1])
        mu         <- cbind(mu[,2],mu[,1])
        medj       <- as.vector(c(medj[2],medj[1]))
        varphi     <- as.vector(c(varphi[2], varphi[1]))
        pii        <- as.vector(c(pii[2], pii[1]))
        sigma2     <- as.vector(c(sigma2[2], sigma2[1]))
      }

      logvero.SNC  <- function(nu) sum(log( d.mixedSNC(y, pii, mu, sigma2, shape, nu) ))
      nu           <- optim(nu, logvero.SNC, control = list(fnscale = -1), method = "L-BFGS-B", lower = rep(0.01, 2), upper = rep(0.99,2))$par

      teta         <- c(Abetas, medj, Delta, Gama, pii, nu)
      tetam[,count] <-  c(Abetas, medj, sigma2, shape, pii, nu)
      auxlog       <- d.mixedSNC(y, pii, mu, sigma2,shape,nu)

      if(length(which(auxlog == 0)) > 0) auxlog[which(auxlog == 0)] <- .Machine$double.xmin
      lk1          <- sum(log(auxlog))
      criterio     <- abs(lk1/lk-1)
      lk           <- lk1
    } # fim while

    if (criteria == TRUE)
    {
      cl  <- apply(tal, 1, which.max)
      icl <- 0
      for (j in 1:g) icl <- icl+sum(log(pii[j]*dSNC(y[cl==j], mu[j], sigma2[j], shape[j], nu)))
    }
    tetam           <- tetam[,1:count]
    EPim            <- im.FMsmsnReg(y, x1, g, model = NULL, initial.model=FALSE,family,Abetas, medj, sigma2=sigma2, shape = shape, pii,nu)
    EP              <- EPim$EP
    end.time        <- Sys.time()
    time.taken      <- end.time - start.time
  } #End skew.cn

  ################################################################################
  ##                               (3)    Skew-slash
  ################################################################################
  if (family == "Skew.slash")
  {
    start.time   <- Sys.time()
    q            <- ncol(x1)
    tetam        <- matrix(data=NA,nrow=4*g + q + 1,ncol=iter.max)
    cat("\nThe Skew.slash can take a long time to run.\n")
    #Parameters for "optimize" function
    TOLERANCIA   <- 1e-8
    MAX_NU       <- 6.5
    MIN_NU       <- 1.01

    k1            <- 2*nu/(2*nu-1)
    b             <- -sqrt(2/pi)*k1
    n             <- length(y)
    p             <- ncol(x1)-1
    x             <- as.matrix(x1[,2:(p+1)])

    beta0         <- Abetas[1]                  ## slope
    betas         <- as.matrix(Abetas[2:(p+1)]) ## parameters of regression

    delta         <- Delta <- Gama <- varphi<-rep(0,g)
    mu            <- matrix(0,n,g)

    for (k in 1:g)
    {
      delta[k]     <- shape[k]/(sqrt(1 + shape[k]^2))
      Delta[k]     <- sqrt(sigma2[k])*delta[k]
      Gama[k]      <- sigma2[k] - Delta[k]^2
      varphi[k]    <- beta0+medj[k]
      mu[,k]       <- x%*%betas+varphi[k]+b*Delta[k]
    }

    teta          <- c(Abetas, medj, Delta, Gama, pii, nu)

    criterio      <- 1
    count         <- 0

    lk            <- sum(log(d.mixedSS(y, pii, mu , sigma2, shape, nu)))## log-likelihood

    while((criterio > error) && (count <= iter.max))
    {#begin while
      count       <- count + 1#;print(count)
      tal         <- matrix(0, n, g)
      S1          <- matrix(0, n, g)
      S2          <- matrix(0, n, g)
      S3          <- matrix(0, n, g)
      soma1       <- matrix(0, p, p)
      soma2       <- matrix(0, p, 1)

      for (j in 1:g)
      {
        ### E-step: computing ui, tui, tui2 ###
        ### mu[,j]<-x%*%betas+b*Delta[j]+varphi[j]
        dj         <- ((y - mu[,j])/sqrt(sigma2[j]))^2
        Mtij2      <- 1/(1 + (Delta[j]^2)*(Gama[j]^(-1)))
        Mtij       <- sqrt(Mtij2)
        mutij      <- Mtij2*Delta[j]*(Gama[j]^(-1))*(y - mu[,j])+b
        A          <- (mutij-b) / Mtij

        u <- vector(mode="numeric",length=n)
        E <- vector(mode="numeric",length=n)

        for(i in 1:n){
          E[i]  <- (((2^(nu + 1))*nu*gamma(nu + 1))/(dSS(y[i], mu[i,j], sigma2[j], shape[j], nu)*pi*sqrt(sigma2[j])))* ((dj[i]+A[i]^2)^(-nu-1))*pgamma(1,nu+1,(dj[i]+A[i]^2)/2)
          faux  <- function(u) u^(nu+0.5)*exp(-u*dj[i]/2)*pnorm(u^(1/2)*A[i])
          aux22 <- integrate(faux,0,1)$value
          u[i]  <- ((sqrt(2)*nu) / (dSS(y[i], mu[i,j], sigma2[j], shape[j], nu)*sqrt(pi)*sqrt(sigma2[j])))*aux22
        }

        d1         <- dSS(y, mu[,j], sigma2[j], shape[j], nu)
        if(length(which(d1 == 0)) > 0) {d1[which(d1 == 0)] <- .Machine$double.xmin}
        d2         <- d.mixedSS(y, pii, mu, sigma2, shape, nu)
        if(length(which(d2 == 0)) > 0) {d2[which(d2 == 0)] <- .Machine$double.xmin}

        tal[,j]    <- d1*pii[j] / d2
        S1[,j]     <- tal[,j]*u
        S2[,j]     <- tal[,j]*(mutij*u + Mtij*E)
        S3[,j]     <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*(mutij+b)*E)

        ### M-step: atualizar mu, Delta, Gama, sigma2, shape ###
        pii[j]     <- (1/n)*sum(tal[,j])
        Delta[j]   <- sum(S2[,j]*(y - mu[,j] + b*Delta[j]))/sum(S3[,j])
        Gama[j]    <- sum(S1[,j]*(y - mu[,j] + b*Delta[j])^2 - 2*(y - mu[,j] + b*Delta[j])*Delta[j]*S2[,j] + Delta[j]^2*S3[,j]) / sum(tal[,j])
        sigma2[j]  <- Gama[j] + Delta[j]^2
        shape[j]   <- sigma2[j]^(-1/2)*Delta[j]/sqrt(1 - Delta[j]^2*sigma2[j]^(-1))
        varphi[j]  <- sum(S1[,j]*(y - x%*%betas) - Delta[j]*S2[,j])/sum(S1[,j])
        soma1      <- soma1 + t(x)%*%diag(S1[,j]/Gama[j])%*%x
        soma2      <- soma2 + t(x)%*%diag(tal[,j]/Gama[j])%*%(u*(y - varphi[j]) - (mutij*u + Mtij*E)*Delta[j])
      }

      betas       <- solve(soma1)%*%soma2
      pii[g]      <- 1 - (sum(pii) - pii[g])

      zero.pos    <- NULL
      zero.pos    <- which(pii == 0)
      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }

      beta0       <- sum(pii*varphi)
      Abetas      <- c(beta0,betas)

      for (k in 1:g)
      {
        mu[,k]    <- x%*%betas+varphi[k]+b*Delta[k]
        medj[k]   <- varphi[k]-beta0
      }

      if (pii[2] < 0.5  && g==2)
      {
        mu          <- cbind(mu[,2],mu[,1])
        medj        <- as.vector(c(medj[2],medj[1]))
        varphi      <- as.vector(c(varphi[2], varphi[1]))
        pii         <- as.vector(c(pii[2], pii[1]))
        sigma2      <- as.vector(c(sigma2[2], sigma2[1]))
      }

      #Here begining the update for estimating the value of nu
      logvero.SS  <- function(nu){ sum(log(d.mixedSS(y, pii, mu, sigma2, shape, nu))) }
      nu          <- optimize(f=logvero.SS, interval=c(MIN_NU,MAX_NU), maximum=TRUE,tol=TOLERANCIA)$maximum
      teta        <- c(Abetas, medj, Delta, Gama, pii, nu)
      tetam[,count] <-  c(Abetas, medj, sigma2, shape, pii, nu)
      auxlog      <- logvero.SS(nu)

      if(length(which(auxlog == 0)) > 0) auxlog[which(auxlog == 0)] <- .Machine$double.xmin
      lk1         <- auxlog
      criterio    <- abs(lk1/lk-1)
      lk          <- lk1
    } # fim while

    if(criteria == TRUE)
    {
      cl           <- apply(tal, 1, which.max)
      icl          <- 0
      for (j in 1:g) {icl <- icl+sum(log(pii[j]*dSS(y, mu[,j], sigma2[j], shape[j], nu)))}
    }
    tetam           <- tetam[,1:count]
    EPim            <- im.FMsmsnReg(y, x1, g, model = NULL, initial.model=FALSE,family,Abetas, medj, sigma2=sigma2, shape = shape, pii,nu)
    EP              <- EPim$EP
    end.time        <- Sys.time()
    time.taken      <- end.time - start.time
  }# end skew.slash


  ################################################################################
  ###                              (4)  Skew.normal
  ################################################################################

  if (family == "Skew.normal")
  {
    start.time   <- Sys.time()
    q            <- ncol(x1)
    tetam        <- matrix(data=NA,nrow=4*g + q,ncol=iter.max)
    k1           <- 1
    b            <- -sqrt(2/pi)*k1

    n            <- length(y)
    p            <- ncol(x1)-1
    x            <- as.matrix(x1[,2:(p+1)])

    beta0        <- Abetas[1]                   ## slope
    betas        <- as.matrix(Abetas[2:(p+1)])  ##  parameters of regression

    delta        <- Delta <- Gama <- varphi<-rep(0,g)
    mu           <- matrix(0,n,g)

    for (k in 1:g)
    {
      delta[k]   <- shape[k] / (sqrt(1 + shape[k]^2))
      Delta[k]   <- sqrt(sigma2[k])*delta[k]
      Gama[k]    <- sigma2[k] - Delta[k]^2
      varphi[k]  <- beta0+medj[k]
      mu[,k]     <- x%*%betas+varphi[k]+b*Delta[k]
    }

    teta         <- c(Abetas, medj, Delta, Gama, pii)

    criterio     <- 1
    count        <- 0

    lk           <- sum(log(d.mixedSN(y, pii, mu, sigma2, shape))) ## log-likelihood

    while((criterio > error) && (count <= iter.max))
    {#begin while
      count      <- count + 1
      tal        <- matrix(0, n, g)
      S1         <- matrix(0, n, g)
      S2         <- matrix(0, n, g)
      S3         <- matrix(0, n, g)
      soma1      <- matrix(0, p,p)
      soma2      <- matrix(0, p,1)

      for (j in 1:g)
      {
        ### E-step: calculando ui, tui, tui2 ###
        ### mu[,j]<-x%*%betas+b*Delta[j]+varphi[j]
        Mtij2     <- 1/(1 + (Delta[j]^2)*(Gama[j]^(-1)))
        Mtij      <- sqrt(Mtij2)
        mutij     <- Mtij2*Delta[j]*(Gama[j]^(-1))*(y - mu[,j])+b

        prob      <- pnorm((mutij-b)/Mtij)
        if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin

        E         <- dnorm((mutij-b)/Mtij)/prob
        u         <- rep(1, n)

        d1        <- dSN(y, mu[,j], sigma2[j], shape[j])
        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
        d2        <- d.mixedSN(y, pii, mu, sigma2, shape)
        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

        tal[,j]   <- d1*pii[j]/d2
        S1[,j]    <- tal[,j]*u
        S2[,j]    <- tal[,j]*(mutij*u + Mtij*E)
        S3[,j]    <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*(mutij+b)*E)

        ### M-step: atualizar mu, Delta, Gama, sigma2 ###
        pii[j]    <- (1/n)*sum(tal[,j])
        Delta[j]  <- sum(S2[,j]*(y - mu[,j] + b*Delta[j]))/sum(S3[,j])
        Gama[j]   <- sum(S1[,j]*(y - mu[,j] + b*Delta[j])^2 - 2*(y - mu[,j] + b*Delta[j])*Delta[j]*S2[,j] + Delta[j]^2*S3[,j]) / sum(tal[,j])
        sigma2[j] <- Gama[j] + Delta[j]^2
        shape[j]  <- sigma2[j]^(-1/2)*Delta[j]/sqrt(1 - Delta[j]^2*sigma2[j]^(-1))
        varphi[j] <- sum(S1[,j]*(y - x%*%betas) - Delta[j]*S2[,j])/sum(S1[,j])
        soma1     <- soma1 + t(x)%*%diag(S1[,j]/Gama[j])%*%x
        soma2     <- soma2 + t(x)%*%diag(tal[,j]/Gama[j])%*%(u*(y - varphi[j]) - (mutij*u + Mtij*E)*Delta[j])
      }

      betas      <- solve(soma1)%*%soma2
      pii[g]     <- 1 - (sum(pii) - pii[g])

      zero.pos   <- NULL
      zero.pos   <- which(pii == 0)
      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }

      beta0      <- sum(pii*varphi)
      Abetas     <- c(beta0,betas)

      for (k in 1:g)
      {
        mu[,k]   <- x%*%betas+varphi[k]+b*Delta[k]
        medj[k]  <- varphi[k]-beta0
      }

      if (pii[2] < 0.5 && g==2)
      {
        mu         <- cbind(mu[,2],mu[,1])
        medj       <- as.vector(c(medj[2],medj[1]))
        varphi     <- as.vector(c(varphi[2], varphi[1]))
        pii        <- as.vector(c(pii[2], pii[1]))
        sigma2     <- as.vector(c(sigma2[2], sigma2[1]))
      }
      teta       <- c(Abetas, medj, Delta, Gama, pii)
      tetam[,count] <-  c(Abetas, medj, sigma2, shape, pii)
      auxlog     <- d.mixedSN(y, pii, mu, sigma2, shape)

      if(length(which(auxlog == 0)) > 0) auxlog[which(auxlog == 0)] <- .Machine$double.xmin
      lk1        <- sum(log(auxlog))
      criterio   <- abs(lk1/lk-1)
      lk         <- lk1
    } # fim while

    if (criteria == TRUE)
    {
      cl           <- apply(tal, 1, which.max)
      icl          <- 0
      for (j in 1:g) {icl<-icl+sum(log(pii[j]*dSN(y, mu[,j], sigma2[j], shape[j])))}
    }
    tetam           <- tetam[,1:count]
    EPim            <- im.FMsmsnReg(y, x1, g, model = NULL, initial.model=FALSE,family,Abetas, medj, sigma2=sigma2, shape = shape, pii)
    EP              <- EPim$EP
    end.time        <- Sys.time()
    time.taken      <- end.time - start.time
  } #end skew.normal

  #Begin of criteria
  if (criteria == TRUE){
    if((family == "t") | (family == "Normal")) d <- p+1+ g*2 + (g-1) #Abeta + varphi, Sigma + pi
    else d <- p+1+ g*3 + (g-1) #Abeta + varphi+ Sigma + lambda+pi
    aic <- -2*lk + 2*d
    bic <- -2*lk + log(n)*d
    edc <- -2*lk + 0.2*sqrt(n)*d
    icl <- -2*icl + log(n)*d
    if(family=="Skew.normal")
    {
      ttable           <- data.frame(cbind(c(Abetas,medj,sigma2, shape, pii[1:(g-1)]),c(EP)))
      #Row names of ttable
      namesrowAbetas   <- c(); for(i in 1:length(Abetas)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
      namesrowMedj     <- c(); for(i in 1:g){namesrowMedj[i]   <- paste("mu",i,sep="")}
      namesrowSigmas   <- c(); for(i in 1:g){namesrowSigmas[i] <- paste("sigma",i,sep="")}
      namesrowShape    <- c(); for(i in 1:g){namesrowShape[i]  <- paste("shape",i,sep="")}
      namesrowPii      <- c(); for(i in 1:(g-1)){namesrowPii[i]  <- paste("pii",i,sep="")}
      rownames(ttable) <- c(namesrowAbetas,namesrowMedj,namesrowSigmas,namesrowShape,namesrowPii)
      colnames(ttable) <- c("Estimate","SE")
      obj.out <- list(tetam = tetam, ttable = ttable, EP=EP, criterio = criterio, time = time.taken, loglik=lk, Abetas = Abetas, medj=medj, sigma2 = sigma2, shape = shape, pii = pii, aic = aic, bic = bic, edc = edc, icl= icl, iter = count, n = length(y), group = cl)
    }

    if(family=="Skew.t" || family=="Skew.slash")
    {
      ttable           <- data.frame(cbind(c(Abetas,medj,sigma2, shape, pii[1:(g-1)], nu),c(EP,NA)))
      #Row names of ttable
      namesrowAbetas   <- c(); for(i in 1:length(Abetas)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
      namesrowMedj     <- c(); for(i in 1:g){namesrowMedj[i]     <- paste("mu",i,sep="")}
      namesrowSigmas   <- c(); for(i in 1:g){namesrowSigmas[i]   <- paste("sigma",i,sep="")}
      namesrowShape    <- c(); for(i in 1:g){namesrowShape[i]    <- paste("shape",i,sep="")}
      namesrowPii      <- c(); for(i in 1:(g-1)){namesrowPii[i]  <- paste("pii",i,sep="")}
      rownames(ttable) <- c(namesrowAbetas,namesrowMedj,namesrowSigmas,namesrowShape,namesrowPii,"nu")
      colnames(ttable) <- c("Estimate","SE")
      obj.out <- list(tetam = tetam, ttable = ttable, EP = EP, criterio = criterio, time = time.taken, loglik=lk, Abetas = Abetas, medj=medj, sigma2 = sigma2, shape = shape, pii = pii, nu = nu, aic = aic, bic = bic, edc = edc, icl= icl, iter = count, n = length(y), group = cl)
    }

    if(family=="Skew.cn")
    {
      ttable           <- data.frame(cbind(c(Abetas,medj,sigma2, shape, pii[1:(g-1)], nu[1],nu[2]),c(EP,NA,NA)))
      #Row names of ttable
      namesrowAbetas   <- c(); for(i in 1:length(Abetas)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
      namesrowMedj     <- c(); for(i in 1:g){namesrowMedj[i]   <- paste("mu",i,sep="")}
      namesrowSigmas   <- c(); for(i in 1:g){namesrowSigmas[i] <- paste("sigma",i,sep="")}
      namesrowShape    <- c(); for(i in 1:g){namesrowShape[i]  <- paste("shape",i,sep="")}
      namesrowPii      <- c(); for(i in 1:(g-1)){namesrowPii[i]  <- paste("pii",i,sep="")}
      rownames(ttable) <- c(namesrowAbetas,namesrowMedj,namesrowSigmas,namesrowShape,namesrowPii,"nu","gamma")
      colnames(ttable) <- c("Estimate","SE")
      obj.out <- list(tetam = tetam, ttable = ttable, EP = EP, criterio = criterio, time = time.taken, loglik=lk, Abetas = Abetas, medj=medj, sigma2 = sigma2, shape = shape, pii = pii, nu = nu, aic = aic, bic = bic, edc = edc, icl= icl, iter = count, n = length(y), group = cl)
    }


  }else{
    if(family=="Skew.normal")
    {
      ttable            <- data.frame(cbind(c(Abetas,medj,sigma2, shape, pii[1:(g-1)]),c(EP)))
      #Row names of ttable
      namesrowAbetas   <- c(); for(i in 1:length(Abetas)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
      namesrowMedj     <- c(); for(i in 1:g){namesrowMedj[i]   <- paste("mu",i,sep="")}
      namesrowSigmas   <- c(); for(i in 1:g){namesrowSigmas[i] <- paste("sigma",i,sep="")}
      namesrowShape    <- c(); for(i in 1:g){namesrowShape[i]  <- paste("shape",i,sep="")}
      rownames(ttable) <- c(namesrowAbetas,namesrowMedj,namesrowSigmas,namesrowShape,"pii1")
      colnames(ttable) <- c("Estimate","SE")
      obj.out <- list(tetam = tetam, ttable = ttable, EP = EP, criterio = criterio, time = time.taken, loglik=lk, Abetas = Abetas, medj=medj, sigma2 = sigma2, shape = shape, pii = pii, iter = count, n = length(y), group = apply(tal, 1, which.max))
    }

    if(family=="Skew.t" || family=="Skew.slash")
    {
      ttable           <- data.frame(cbind(c(Abetas,medj,sigma2, shape, pii[1:(g-1)], nu),c(EP,NA)))
      #Row names of ttable
      namesrowAbetas   <- c(); for(i in 1:length(Abetas)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
      namesrowMedj     <- c(); for(i in 1:g){namesrowMedj[i]   <- paste("mu",i,sep="")}
      namesrowSigmas   <- c(); for(i in 1:g){namesrowSigmas[i] <- paste("sigma",i,sep="")}
      namesrowShape    <- c(); for(i in 1:g){namesrowShape[i]  <- paste("shape",i,sep="")}
      namesrowPii      <- c(); for(i in 1:(g-1)){namesrowPii[i]  <- paste("pii",i,sep="")}
      rownames(ttable) <- c(namesrowAbetas,namesrowMedj,namesrowSigmas,namesrowShape,namesrowPii,"nu")
      colnames(ttable) <- c("Estimate")
      obj.out <- list(tetam = tetam, ttable = ttable, EP = EP, criterio = criterio, time = time.taken, loglik=lk, Abetas = Abetas, medj=medj, sigma2 = sigma2, shape = shape, pii = pii, nu = nu, iter = count, n = length(y), group = apply(tal, 1, which.max))
    }

    if(family=="Skew.cn")
    {
      ttable           <- data.frame(cbind(c(Abetas,medj,sigma2, shape, pii[1:(g-1)], nu[1],nu[2]),c(EP,NA,NA)))
      #Row names of ttable
      namesrowAbetas   <- c(); for(i in 1:length(Abetas)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
      namesrowMedj     <- c(); for(i in 1:g){namesrowMedj[i]   <- paste("mu",i,sep="")}
      namesrowSigmas   <- c(); for(i in 1:g){namesrowSigmas[i] <- paste("sigma",i,sep="")}
      namesrowShape    <- c(); for(i in 1:g){namesrowShape[i]  <- paste("shape",i,sep="")}
      namesrowPii      <- c(); for(i in 1:(g-1)){namesrowPii[i]  <- paste("pii",i,sep="")}
      rownames(ttable) <- c(namesrowAbetas,namesrowMedj,namesrowSigmas,namesrowShape,namesrowPii,"nu","gamma")
      colnames(ttable) <- c("Estimate","SE")
      obj.out <- list(tetam = tetam, ttable = ttable, EP = EP, criterio = criterio, time = time.taken, loglik=lk, Abetas = Abetas, medj=medj, sigma2 = sigma2, shape = shape, pii = pii, nu = nu, iter = count, n = length(y), group = apply(tal, 1, which.max))
    }
  }
  if (group == FALSE) obj.out <- obj.out[-(length(obj.out))]
  if (obs.prob == TRUE)
  {
    nam   <- c()
    for (i in 1:ncol(tal)) nam <- c(nam,paste("Group ",i,sep=""))
    if(ncol(tal) == 1){ dimnames(tal)[[2]] <- list(nam) }
    if(ncol(tal) > 1) { dimnames(tal)[[2]] <- nam }
    obj.out$obs.prob <- tal
    if((ncol(tal) - 1) > 1) obj.out$obs.prob[,ncol(tal)] <- 1 - rowSums(obj.out$obs.prob[,1:(ncol(tal)-1)])
    else obj.out$obs.prob[,ncol(tal)] <- 1 - obj.out$obs.prob[,1]
    obj.out$obs.prob[which(obj.out$obs.prob[,ncol(tal)] <0),ncol(tal)] <- 0.0000000000
    obj.out$obs.prob <- round(obj.out$obs.prob,10)
  }

  class(obj.out) <- family
  return(obj.out)
}


#g = 2; get.init = TRUE; criteria = TRUE; group = FALSE; family = "Skew.slash"; error = 10^-8; iter.max = 2000; obs.prob= FALSE; kmeans.param = NULL
