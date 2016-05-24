#######################################################################
######## Aplicacao do Artigo MIX-SNI Classico - dados BMI  ############
#######################################################################
##library(mvtnorm)
matrix.sqrt <- function(A) {
  sva <- svd(A)
  if (min(sva$d)>=0)
    Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
  else
    stop("Matrix square root is not defined.\n")
  return(Asqrt)
}
#####################################################################
#########         Matriz de Info - Univariado         ###############

im.smsn <- function(y, model){
  if((class(model) != "t") && (class(model) != "Skew.t") && (class(model) != "Skew.cn") && (class(model) != "Skew.slash") && (class(model) != "Skew.normal") && (class(model) != "Normal")) stop(paste("Family",class(model),"not recognized.\n",sep=" "))
  y <- as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  g <- length(model$pii)

  Sipi <- Simu <- Silambda <- c()
  Ssigma11 <- c()

  nam.nu <- NULL
  nam.gr <- NULL

  mu <- model$mu
  Sigma <- model$sigma2
  lambda <- model$shape
  pii <- model$pii
  nu <- model$nu

  if (class(model) == "t"){
    soma <- soma2 <- 0

    I.Phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric((( 2^w*nu^(nu/2)*gamma(w + nu/2))/(gamma(nu/2)*(nu + di)^(nu/2 + w)))*pt( ((Ai)/(di + nu)^(0.5))*sqrt(2*w + nu), 2*w + nu))
    I.phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric(((2^w*nu^(nu/2))/(sqrt(2*pi)*gamma(nu/2)))*(1/(di + Ai^2 + nu))^((nu + 2*w)/2)*gamma((nu + 2*w)/2))

    for (i in 1:n){
      S <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu <- 0
      for (j in 1:g){
        yi <- matrix(y[i,], 1,1)
        Ai <- as.numeric(t(lambda[j])%*%solve(matrix.sqrt(Sigma[j]))%*%(y[i,] - mu[j]))
        di <- as.numeric(mahalanobis(yi, mu[j], Sigma[j]))

        #derivadinhas
        Dr <- as.numeric(matrix.sqrt(Sigma[[j]]))
        Dadj <- 1/(2*Dr)
        dir.dmu <- (-2*(y[i,] - mu[j]))/Sigma[j]
        dAir.dmu <- -lambda[j]/Dr

        #para os elementos de sigma
        ddet.ds11 <- -(1/Dr^2)*Dadj
        dir.ds11 <-  - (y[i,] - mu[j])^2/Dr^4
        dAir.ds11 <- -(0.5* lambda[j]*(y[i,] - mu[j]))/Dr^3

        dPsi.dmu <- 2/sqrt(2*pi*Sigma[j])*( dAir.dmu * I.phi((p+1)/2, Ai, di, nu) - (1/2)*dir.dmu*I.Phi((p/2)+1, Ai, di, nu) )
        dPsi.dsigma11 <- (2/(2*pi)^(p/2))*( ddet.ds11*I.Phi(p/2, Ai, di, nu) - (1/2)*dir.ds11/sqrt(Sigma[j])*I.Phi(p/2+1, Ai, di, nu) + dAir.ds11/sqrt(Sigma[j])*I.phi((p+1)/2, Ai, di, nu) )

        ui <- rgamma(10000, shape = nu/2, rate = nu/2)
        resto <- mean(ui^(p/2)*log(ui)*exp(-ui*di/2)*pnorm(ui^(1/2)*Ai))
        dPsi.dnu <- dPsi.dnu + pii[j]*((det(as.matrix(Sigma[j]))^(-1/2))/(2*pi)^(p/2))*((log(nu/2)+1-digamma(nu/2))*I.Phi(p/2, Ai, di, nu) - I.Phi((p+2)/2, Ai, di, nu) + resto)

        Simu <- as.vector((pii[j]/ d.mixedST(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dmu)
        Ssigma11 <- (pii[j]/ d.mixedST(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dsigma11
        Sipi <- (1/d.mixedST(yi, pii, mu, Sigma, lambda, nu))*( dt.ls(yi, mu[j], Sigma[j], lambda[j], nu) - dt.ls(yi, mu[g], Sigma[g], lambda[g], nu))

        S <- c(S, Simu, Ssigma11, Sipi)

      }
      S <- S[-length(S)]
      Sinu <- (1/d.mixedST(yi, pii, mu, Sigma, lambda, nu))*dPsi.dnu
      S <- c(S, Sinu)

      soma <- soma + S%*%t(S)
##      soma2 <- soma2 + S
    }
      for (i in 1:g) nam.gr <- c(nam.gr,paste("mu",i,sep=""),paste("sigma",i,sep=""),paste("p",i,sep=""))
      nam.gr <- nam.gr[-length(nam.gr)]
      nam.nu <- c("nu")
  }

  if (class(model) == "Skew.t"){
    soma <- soma2 <- 0

    I.Phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric((( 2^w*nu^(nu/2)*gamma(w + nu/2))/(gamma(nu/2)*(nu + di)^(nu/2 + w)))*pt( ((Ai)/(di + nu)^(0.5))*sqrt(2*w + nu), 2*w + nu))
    I.phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric(((2^w*nu^(nu/2))/(sqrt(2*pi)*gamma(nu/2)))*(1/(di + Ai^2 + nu))^((nu + 2*w)/2)*gamma((nu + 2*w)/2))

    for (i in 1:n){
      S <- c() #vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu <- 0
      for (j in 1:g){
        yi <- matrix(y[i,], 1,1)
        Ai <- as.numeric(t(lambda[j])%*%solve(matrix.sqrt(Sigma[j]))%*%(y[i,] - mu[j]))
        di <- as.numeric(mahalanobis(yi, mu[j], Sigma[j]))

        #derivadinhas
        Dr <- as.numeric(matrix.sqrt(Sigma[[j]]))
        Dadj <- 1/(2*Dr)
        dir.dmu <- (-2*(y[i,] - mu[j]))/Sigma[j]
        dAir.dmu <- -lambda[j]/Dr
        dAir.dlambda <- (y[i,] - mu[j])/Dr

        #para os elementos de sigma
        ddet.ds11 <- -(1/Dr^2)*Dadj
        dir.ds11 <-  - (y[i,] - mu[j])^2/Dr^4
        dAir.ds11 <- -(0.5* lambda[j]*(y[i,] - mu[j]))/Dr^3

        dPsi.dmu <- 2/sqrt(2*pi*Sigma[j])*( dAir.dmu * I.phi((p+1)/2, Ai, di, nu) - (1/2)*dir.dmu*I.Phi((p/2)+1, Ai, di, nu) )
        dPsi.dlambda <- (2/sqrt(2*pi*Sigma[j]))*dAir.dlambda*I.phi((p+1)/2, Ai, di, nu)
        dPsi.dsigma11 <- (2/(2*pi)^(p/2))*( ddet.ds11*I.Phi(p/2, Ai, di, nu) - (1/2)*dir.ds11/sqrt(Sigma[j])*I.Phi(p/2+1, Ai, di, nu) + dAir.ds11/sqrt(Sigma[j])*I.phi((p+1)/2, Ai, di, nu) )

        ui <- rgamma(10000, shape = nu/2, rate = nu/2)
        resto <- mean(ui^(p/2)*log(ui)*exp(-ui*di/2)*pnorm(ui^(1/2)*Ai))
        dPsi.dnu <- dPsi.dnu + pii[j]*((det(as.matrix(Sigma[j]))^(-1/2))/(2*pi)^(p/2))*((log(nu/2)+1-digamma(nu/2))*I.Phi(p/2, Ai, di, nu) - I.Phi((p+2)/2, Ai, di, nu) + resto)

        Simu <- as.vector((pii[j]/ d.mixedST(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dmu)
        Silambda <- as.vector((pii[j]/ d.mixedST(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dlambda )
        Ssigma11 <- (pii[j]/ d.mixedST(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dsigma11
        Sipi <- (1/d.mixedST(yi, pii, mu, Sigma, lambda, nu))*( dt.ls(yi, mu[j], Sigma[j], lambda[j], nu) - dt.ls(yi, mu[g], Sigma[g], lambda[g], nu))

        S <- c(S, Simu, Ssigma11, Silambda, Sipi)

      }
      S <- S[-length(S)]
      Sinu <- (1/d.mixedST(yi, pii, mu, Sigma, lambda, nu))*dPsi.dnu
      S <- c(S, Sinu)

      soma <- soma + S%*%t(S)
##      soma2 <- soma2 + S
    }
      for (i in 1:g) nam.gr <- c(nam.gr,paste("mu",i,sep=""),paste("sigma",i,sep=""),paste("shape",i,sep=""),paste("p",i,sep=""))
      nam.gr <- nam.gr[-length(nam.gr)]
      nam.nu <- c("nu")
  }

  if (class(model) == "Skew.cn"){
    soma <- soma2 <- 0

    I.Phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric( sqrt(2*pi)*(nu[1]*nu[2]^(w -0.5)*dnorm(sqrt(di), 0, sqrt(1/nu[2]))*pnorm(nu[2]^(1/2)*Ai) + (1 - nu[1])*(dnorm(sqrt(di), 0,1)*pnorm(Ai)) )   )
    I.phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric( nu[1]*nu[2]^(w - 0.5)*dnorm(sqrt(di + Ai^2), 0, sqrt(1/nu[2])) + (1 - nu[1])*dnorm(sqrt(di + Ai^2))   )

    for (i in 1:n){
      S <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu1 <- dPsi.dnu2 <- 0
      for (j in 1:g){
        yi <- matrix(y[i,], 1,1)
        Ai <- as.numeric(t(lambda[j])%*%solve(matrix.sqrt(Sigma[j]))%*%(y[i,] - mu[j]))
        di <- as.numeric(mahalanobis(yi, mu[j], Sigma[j]))

        #derivadinhas
        Dr <- as.numeric(matrix.sqrt(Sigma[[j]]))
        Dadj <- 1/(2*Dr)
        dir.dmu <- (-2*(y[i,] - mu[j]))/Sigma[j]
        dAir.dmu <- -lambda[j]/Dr
        dAir.dlambda <- (y[i,] - mu[j])/Dr

        #para os elementos de sigma
        ddet.ds11 <- -(1/Dr^2)*Dadj
        dir.ds11 <-  - (y[i,] - mu[j])^2/Dr^4
        dAir.ds11 <- -(0.5* lambda[j]*(y[i,] - mu[j]))/Dr^3

        dPsi.dmu <- 2/sqrt(2*pi*Sigma[j])*( dAir.dmu * I.phi((p+1)/2, Ai, di, nu) - (1/2)*dir.dmu*I.Phi((p/2)+1, Ai, di, nu) )
        dPsi.dlambda <- (2/sqrt(2*pi*Sigma[j]))*dAir.dlambda*I.phi((p+1)/2, Ai, di, nu)
        dPsi.dsigma11 <- (2/(2*pi)^(p/2))*( ddet.ds11*I.Phi(p/2, Ai, di, nu) - (1/2)*dir.ds11/sqrt(Sigma[j])*I.Phi(p/2+1, Ai, di, nu) + dAir.ds11/sqrt(Sigma[j])*I.phi((p+1)/2, Ai, di, nu) )
        dPsi.dnu1 <- dPsi.dnu1 + pii[j]*2*(dmvnorm(yi, mu[j], nu[2]^(-1)*as.matrix(Sigma[j]))*pnorm(nu[2]^(1/2)*Ai) - dmvnorm(yi, mu[j], as.matrix(Sigma[j]))*pnorm(Ai) )
        dPsi.dnu2 <- dPsi.dnu1 + pii[j]*((nu[1]*det(as.matrix(Sigma[j]))^(-1/2)*nu[2]^(p/2))/(2*pi)^(p/2))*exp(-nu[2]*di/2)*(p*nu[2]^(-1)*pnorm(nu[2]^(1/2)*Ai) + dnorm(nu[2]^(1/2)*Ai)*Ai*nu[2]^(-1/2) - pnorm(nu[2]^(1/2)*Ai)*di )

        Simu <- as.vector((pii[j]/ d.mixedSNC(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dmu)
        Silambda <- as.vector((pii[j]/ d.mixedSNC(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dlambda )

        Ssigma11 <- (pii[j]/ d.mixedSNC(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dsigma11
        Sipi <- (1/d.mixedSNC(yi, pii, mu, Sigma, lambda, nu))*( dSNC(yi, mu[j], Sigma[j], lambda[j], nu) - dSNC(yi, mu[g], Sigma[g], lambda[g], nu))

        S <- c(S, Simu, Ssigma11, Silambda, Sipi)
      }
      S <- S[-length(S)]
      Sinu1 <- (1/d.mixedSNC(yi, pii, mu, Sigma, lambda, nu))*dPsi.dnu1
      Sinu2 <- (1/d.mixedSNC(yi, pii, mu, Sigma, lambda, nu))*dPsi.dnu2
      S <- c(S, Sinu1, Sinu2)
      soma <- soma + S%*%t(S)
##      soma2 <- soma2 + S
    }
      for (i in 1:g) nam.gr <- c(nam.gr,paste("mu",i,sep=""),paste("sigma",i,sep=""),paste("shape",i,sep=""),paste("p",i,sep=""))
      nam.gr <- nam.gr[-length(nam.gr)]
      nam.nu <- c("nu1","nu2") 

  }

    if (class(model) == "Skew.slash"){
    soma <- soma2 <- 0

    I.Phi <- function(w=0, Ai=0, di=0, nu=0) {
      f <- function(u) pnorm(u^(0.5)*Ai)*dgamma(u,nu+w,di/2)
      resp <- integrate(f,0,1)$value
      res1 <- (nu*(2^(w + nu)*gamma(w + nu))/(di^(w + nu)))*resp
      return(res1)
    }

    I.phi <- function(w=0, Ai=0, di=0, nu=0){
      res2 <- ((nu*2^(w + nu)*gamma(w + nu))/(sqrt(2*pi)*(di + Ai^2)^(w + nu)))*pgamma(1, w + nu, (di + Ai^2)/2)
      return(res2)
    }

    for (i in 1:n){
      S <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      dPsi.dnu <- 0
      for (j in 1:g){
        yi <- matrix(y[i,], 1,1)
        Ai <- as.numeric(t(lambda[j])%*%solve(matrix.sqrt(Sigma[j]))%*%(y[i,] - mu[j]))
        di <- as.numeric(mahalanobis(yi, mu[j], Sigma[j]))

        #derivadinhas
        Dr <- as.numeric(matrix.sqrt(Sigma[[j]]))
        Dadj <- 1/(2*Dr)
        dir.dmu <- (-2*(y[i,] - mu[j]))/Sigma[j]
        dAir.dmu <- -lambda[j]/Dr
        dAir.dlambda <- (y[i,] - mu[j])/Dr

        #para os elementos de sigma
        ddet.ds11 <- -(1/Dr^2)*Dadj
        dir.ds11 <-  - (y[i,] - mu[j])^2/Dr^4
        dAir.ds11 <- -(0.5* lambda[j]*(y[i,] - mu[j]))/Dr^3

        dPsi.dmu <- 2/sqrt(2*pi*Sigma[j])*( dAir.dmu * I.phi((p+1)/2, Ai, di, nu) - (1/2)*dir.dmu*I.Phi((p/2)+1, Ai, di, nu) )
        dPsi.dlambda <- (2/sqrt(2*pi*Sigma[j]))*dAir.dlambda*I.phi((p+1)/2, Ai, di, nu)
        dPsi.dsigma11 <- (2/(2*pi)^(p/2))*( ddet.ds11*I.Phi(p/2, Ai, di, nu) - (1/2)*dir.ds11/sqrt(Sigma[j])*I.Phi(p/2+1, Ai, di, nu) + dAir.ds11/sqrt(Sigma[j])*I.phi((p+1)/2, Ai, di, nu) )

        #f <- function(u) 2*u^(nu - 1)*(1 + nu*log(u))*dnorm(yi,mu[j],sqrt(Sigma[j]/u))*pnorm(u^(1/2)*lambda[j]*(Sigma[j]^(-1/2))*(yi-mu[j]))
        u <- runif(5000)
        dPsi.dnu <- dPsi.dnu + pii[j]*mean(2*u^(nu - 1)*(1 + nu*log(u))*dnorm(yi,mu[j],sqrt(Sigma[j]/u))*pnorm(u^(1/2)*lambda[j]*(Sigma[j]^(-1/2))*(yi-mu[j])))#integrate(f,0,1)$value

        Simu <- as.vector((pii[j]/ d.mixedSS(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dmu)
        Silambda <- as.vector((pii[j]/ d.mixedSS(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dlambda )
        Ssigma11 <- (pii[j]/ d.mixedSS(yi, pii, mu, Sigma, lambda, nu) )*dPsi.dsigma11
        Sipi <- (1/d.mixedSS(yi, pii, mu, Sigma, lambda, nu))*( dSS(yi, mu[j], Sigma[j], lambda[j], nu) - dSS(yi, mu[g], Sigma[g], lambda[g], nu))

        S <- c(S, Simu, Ssigma11, Silambda, Sipi)
      }
      S <- S[-length(S)]
      Sinu <- (1/d.mixedSS(yi, pii, mu, Sigma, lambda, nu))*dPsi.dnu
      S <- c(S, Sinu)
      soma <- soma + S%*%t(S)
##      soma2 <- soma2 + S
    }
      for (i in 1:g) nam.gr <- c(nam.gr,paste("mu",i,sep=""),paste("sigma",i,sep=""),paste("shape",i,sep=""),paste("p",i,sep=""))
      nam.gr <- nam.gr[-length(nam.gr)]
      nam.nu <- c("nu") 
  }

  if (class(model) == "Skew.normal"){
    soma <- soma2 <- 0

    I.Phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric( exp(-di/2)*pnorm(Ai) )
    I.phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric( exp(-di/2)*dnorm(Ai) )
    
    for (i in 1:n){
      S <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      for (j in 1:g){
        yi <- matrix(y[i,], 1,1)
        Ai <- as.numeric(t(lambda[j])%*%solve(matrix.sqrt(Sigma[j]))%*%(y[i,] - mu[j]))
        di <- as.numeric(mahalanobis(yi, mu[j], Sigma[j]))

        #derivadinhas
        Dr <- as.numeric(matrix.sqrt(Sigma[[j]]))
        Dadj <- 1/(2*Dr)
        dir.dmu <- (-2*(y[i,] - mu[j]))/Sigma[j]
        dAir.dmu <- -lambda[j]/Dr
        dAir.dlambda <- (y[i,] - mu[j])/Dr

        #para os elementos de sigma
        ddet.ds11 <- -(1/Dr^2)*Dadj
        #ddet.ds11 <- -0.5*(1/det(Dr)^6)*Dadj
        #ddet.ds11 <- -0.5*(1/det(Dr)^3)
        #dir.ds11 <- - t(y[i,] - mu[j])%*%solve(Dr)%*%(D11%*%solve(Dr) + solve(Dr)%*%D11)%*%solve(Dr)%*%(y[i,] - mu[j])
        dir.ds11 <-  - (y[i,] - mu[j])^2/Dr^4
        #dAir.ds11 <- - t(lambda[j])%*%solve(Dr)%*%D11%*%solve(Dr)%*%(y[i,] - mu[j])
        dAir.ds11 <- -(0.5* lambda[j]*(y[i,] - mu[j]))/Dr^3

        dPsi.dmu <- 2/sqrt(2*pi*Sigma[j])*( dAir.dmu * I.phi(Ai=Ai, di=di) - (1/2)*dir.dmu*I.Phi(Ai=Ai, di=di) )
        dPsi.dlambda <- (2/sqrt(2*pi*Sigma[j]))*dAir.dlambda*I.phi(Ai=Ai, di=di)
        dPsi.dsigma11 <- (2/(2*pi)^(p/2))*( ddet.ds11*I.Phi(Ai=Ai, di=di) - (1/2)*dir.ds11/sqrt(Sigma[j])*I.Phi(Ai=Ai, di=di) + dAir.ds11/sqrt(Sigma[j])*I.phi(Ai=Ai, di=di) )
        #dPsi.dsigma11 <- (2/(2*pi)^(p/2))*( ddet.ds11*I.Phi(Ai=Ai, di=di) - dir.ds11*I.Phi(Ai=Ai, di=di) + 2*dAir.ds11*I.phi(Ai=Ai, di=di) )

        Simu <- as.vector((pii[j]/ d.mixedSN(yi, pii, mu, Sigma, lambda) )*dPsi.dmu)
        Silambda <- as.vector((pii[j]/ d.mixedSN(yi, pii, mu, Sigma, lambda) )*dPsi.dlambda )
        Ssigma11 <- (pii[j]/ d.mixedSN(yi, pii, mu, Sigma, lambda)) * dPsi.dsigma11
        Sipi <- (1/d.mixedSN(yi, pii, mu, Sigma, lambda))*( dSN(yi, mu[j], Sigma[j], lambda[j]) - dSN(yi, mu[g], Sigma[g], lambda[g]))

        S <- c(S, Simu, Ssigma11, Silambda, Sipi)
      }
      S <- S[-length(S)]
      soma <- soma + S%*%t(S)
##      soma2 <- soma2 + S

    }
      for (i in 1:g) nam.gr <- c(nam.gr,paste("mu",i,sep=""),paste("sigma",i,sep=""),paste("shape",i,sep=""),paste("p",i,sep=""))
      nam.gr <- nam.gr[-length(nam.gr)]
  }

  if (class(model) == "Normal"){
    soma <- soma2 <- 0

    I.Phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric( exp(-di/2)*1/2 )
    I.phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric( exp(-di/2)*dnorm(0) )


    for (i in 1:n){
      S <- c() # vetor com todas as derivadas em relacao a cada parametro desconhecido do modelo
      for (j in 1:g){
        yi <- matrix(y[i,], 1,1)
        Ai <- as.numeric(t(lambda[j])%*%solve(matrix.sqrt(as.matrix(Sigma[j])))%*%(y[i,] - mu[j]))
        di <- as.numeric(mahalanobis(yi, mu[j], Sigma[j]))

        #derivadinhas
        Dr <- as.numeric(matrix.sqrt(Sigma[[j]]))
        Dadj <- 1/(2*Dr)
        dir.dmu <- (-2*(y[i,] - mu[j]))/Sigma[j]
        dAir.dmu <- -lambda[j]/Dr
        #dAir.dlambda <- (y[i,] - mu[j])/Dr

        #para os elementos de sigma
        ddet.ds11 <- -(1/Dr^2)*Dadj
        dir.ds11 <-  - (y[i,] - mu[j])^2/Dr^4
        dAir.ds11 <- -(0.5* lambda[j]*(y[i,] - mu[j]))/Dr^3

        dPsi.dmu <- 2/sqrt(2*pi*Sigma[j])*( dAir.dmu * I.phi(di=di) - (1/2)*dir.dmu*I.Phi(di=di) )
        dPsi.dsigma11 <- (2/(2*pi)^(p/2))*( ddet.ds11*I.Phi(di=di) - (1/2)*dir.ds11/sqrt(Sigma[j])*I.Phi(di=di) + dAir.ds11/sqrt(Sigma[j])*I.phi(di=di) )
        
        Simu <- as.vector((pii[j]/ d.mixedSN(yi, pii, mu, Sigma, lambda) )*dPsi.dmu)
        Ssigma11 <- (pii[j]/ d.mixedSN(yi, pii, mu, Sigma, lambda) )*dPsi.dsigma11
        Sipi <- (1/d.mixedSN(yi, pii, mu, Sigma, lambda))*( dSN(yi, mu[j], Sigma[j], lambda[j]) - dSN(yi, mu[g], Sigma[g], lambda[g]))

        S <- c(S, Simu, Ssigma11, Sipi)

      }
      S <- S[-length(S)]
      soma <- soma + S%*%t(S)
##      soma2 <- soma2 + S
      }
      for (i in 1:g) nam.gr <- c(nam.gr,paste("mu",i,sep=""),paste("sigma",i,,sep=""),paste("p",i,sep=""))
      nam.gr <- nam.gr[-length(nam.gr)]

  }

  dimnames(soma)[[1]] <- c(nam.gr,nam.nu)
  dimnames(soma)[[2]] <- c(nam.gr,nam.nu)
  if(class(model) == "t"){
     
  }
  ##return(list(soma, soma2))
  return(list(IM=soma))
}

#####################################################################
#########         Procura o melhor fitting         ###############
smsn.search <- function(y, nu, g.min = 1, g.max = 3, family = "Skew.normal", criteria = "bic", error = 0.0001, iter.max = 100, calc.im = FALSE, uni.Gama = FALSE, kmeans.param = NULL)
{
 y <- as.matrix(y)
 if( g.min > g.max) stop("The number g.min should be less than g.max.\n")
 if( (criteria != "aic") && (criteria != "bic") && (criteria != "edc") && (criteria != "icl")) stop ("Criterion should be \"aic\" or \"bic\" or \"edc\".\n")

 g.size <- paste("g=",g.min,sep="")
 num <- g.max - g.min
  if (num > 0){             
     for (i in 1:num) {      
       k <- paste("g=",g.min+i,sep="")  
       g.size <- c(g.size,k) 
     }                     
  }                      
  g <- seq(g.min,g.max,by=1)

output <- list()
crit <- rep(0,length(g))

if(ncol(y) == 1){
  for(i in 1:length(g)){
     output[[g.size[i]]] <- smsn.mix(y=y, nu=nu, g = g[i], get.init = TRUE, criteria = TRUE, group = TRUE, family = family, error = error, iter.max = iter.max, calc.im = calc.im, kmeans.param)
     crit[i] <- output[[g.size[i]]][[criteria]]
  }
}
else{
    for(i in 1:length(g)){
       output[[g.size[i]]] <- smsn.mmix(y=y, nu=nu, g = g[i], get.init = TRUE, criteria = TRUE, group = TRUE, family = family, error = error, iter.max = iter.max, uni.Gama = uni.Gama, calc.im = calc.im, kmeans.param)
       crit[i] <- output[[g.size[i]]][[criteria]]
    }
}
best <- which(crit == min(crit))
names(crit) <- g.size

out <- list(criteria = crit, best.model = output[[g.size[best]]])

out
}

