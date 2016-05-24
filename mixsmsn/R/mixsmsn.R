#####################################################################
##########       Algoritmo EM para mistura            ###############
#                   alterado em 26/08/08                            #

smsn.mix <- function(y, nu, mu = NULL, sigma2 = NULL, shape = NULL, pii = NULL, g = NULL, get.init = TRUE, criteria = TRUE, group = FALSE, 
                     family = "Skew.normal", error = 0.00001, iter.max = 100, calc.im = TRUE, obs.prob= FALSE, kmeans.param = NULL){
  #y: eh o vetor de dados (amostra) de tamanho n
  #mu, sigma2, shape, pii: sao os valores iniciais para o algoritmo EM. Cada um deles deve ser um vetor de tamanho g
  #                       (o algoritmo entende o numero de componentes a ajustar baseado no tamanho desses vetores)
  #nu: valor inicial para o nu (no caso da SNC deve ser um vetor bidimensional com valores entre 0 e 1)
  #loglik: TRUE ou FALSE - Caso queira que seja calculada a log-verossimilhanca do modelo ajustado
  #cluster: TRUE ou FALSE - Caso queira um vetor de tamanho n indicando a qual componente a i-esima observacao pertence
  #type: c("ST","SNC","SS","SN","Normal") - ST: Ajusta por misturas de Skew-T
  #                                       - SNC: Ajusta por misturas de Skew Normal Contaminada
  #                                       - SS: Ajusta por misturas de Skew-Slash (ATENCAO: AINDA COM PROBLEMAS)
  #                                       - SN: Ajusta por misturas de Skew Normal
  #                                       - Normal: Misturas de Normais
  #error: define o criterio de parada do algoritmo.
  if(ncol(as.matrix(y)) > 1) stop("This function is only for univariate response y!") 
  if((family != "t") && (family != "Skew.t") && (family != "Skew.cn") && (family != "Skew.slash") && (family != "Skew.normal") && (family != "Normal")) stop(paste("Family",family,"not recognized.\n",sep=" "))
  if((length(g) == 0) && ((length(mu)==0) || (length(sigma2)==0) || (length(shape)==0) || (length(pii)==0)))  stop("The model is not specified correctly.\n")
  if(get.init == FALSE){
       g <- length(mu)
       if((length(mu) != length(sigma2)) || (length(mu) != length(pii))) stop("The size of the initial values are not compatibles.\n")
       if((family == "Skew.t" || family == "Skew.cn" || family == "Skew.slash" || family == "Skew.normal") & (length(mu) != length(shape))) stop("The size of the initial values are not compatibles.\n")
       if(sum(pii) != 1) stop("probability of pii does not sum to 1.\n")
  }  
  if((length(g)!= 0) && (g < 1)) stop("g must be greater than 0.\n")

  if (get.init == TRUE){
    if(length(g) == 0) stop("g is not specified correctly.\n")

    k.iter.max <- 50
    n.start <- 1
    algorithm <- "Hartigan-Wong"
    if(length(kmeans.param) > 0){
       if(length(kmeans.param$iter.max) > 0 ) k.iter.max <- kmeans.param$iter.max
       if(length(kmeans.param$n.start) > 0 ) n.start <- kmeans.param$n.start
       if(length(kmeans.param$algorithm) > 0 ) algorithm <- kmeans.param$algorithm
    }

    if(g > 1){
      init <- kmeans(y,g,k.iter.max,n.start,algorithm)
      pii <- init$size/length(y)
      mu <- as.vector(init$centers)
      sigma2 <- init$withinss/init$size
      shape <- c()
      for (j in 1:g){
        m3 <- (1/init$size[j])*sum( (y[init$cluster == j] - mu[j])^3 )
        shape[j] <- sign(m3/sigma2[j]^(3/2))
      }
    }
   else{
     mu <- mean(y)
     sigma2 <- var(y)
     pii <- 1
     m3 <- (1/length(y))*sum( (y - mu)^3 )
     shape <- sign(m3/sigma2^(3/2))
   }
  }

  if (family == "t"){
      shape <- rep(0,g)
      lk <- sum(log(d.mixedST(y, pii, mu, sigma2, shape, nu) ))
      n <- length(y)
      delta <- Delta <- Gama <- rep(0,g)
      for (k in 1:g){
        delta[k] <- 0 ##shape[k] / (sqrt(1 + shape[k]^2))
        Delta[k] <- 0 ##sqrt(sigma2[k])*delta[k]
        Gama[k] <- sigma2[k] - Delta[k]^2
      }

      teta <- c(mu, Delta, Gama, pii, nu)
      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama

      criterio <- 1
      count <- 0

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
 ##       print(count)
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu[j])/sqrt(sigma2[j]))^2
          Mtij2 <- 1/(1 + (Delta[j]^2)*(Gama[j]^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta[j]*(Gama[j]^(-1))*(y - mu[j])
          A <- mutij / Mtij

          E=(2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2[j])*dt.ls(y, mu[j], sigma2[j],shape[j] ,nu))
          u= ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2[j])*dt.ls(y, mu[j], sigma2[j],shape[j] ,nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)

          d1 <- dt.ls(y, mu[j], sigma2[j], shape[j], nu)
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <-d.mixedST(y, pii, mu, sigma2, shape, nu)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2
          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)


          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])
          mu[j] <- sum(S1[,j]*y - Delta.old[j]*S2[,j]) / sum(S1[,j])
          Delta[j] <- 0
          Gama[j] <- sum(S1[,j]*(y - mu[j])^2 - 2*(y - mu[j])*Delta[j]*S2[,j] + Delta[j]^2*S3[,j]) / sum(tal[,j])       
          sigma2[j] <- Gama[j] + Delta[j]^2
          shape[j] <- 0

        }
        #aqui comecam as alteracoes para estimar o valor de nu
        logvero.ST <- function(nu) sum(log(d.mixedST(y, pii, mu, sigma2, shape, nu)))
        nu <- optimize(logvero.ST, c(0,100), tol = 0.000001, maximum = TRUE)$maximum
        lk1 <- sum(log( d.mixedST(y, pii, mu, sigma2, shape, nu) ))
        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

        param <- teta
        teta <- c(mu, Delta, Gama, pii, nu)
        #criterio <- sqrt((teta-param)%*%(teta-param))
        criterio <- abs(lk1/lk-1)

        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama
        lk<-lk1
        
      }

    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl <- icl+sum(log(pii[j]*dt.ls(y[cl==j], mu[j], sigma2[j], shape[j], nu)))
       #icl=-2*icl+4*g*log(n)
    }

      #### grafico ajuste
      ##hist(y, breaks = 40,probability=T,col="grey")
      ##xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      ##lines(xx,d.mixedST(xx, pii, mu, sigma2, shape, nu),col="red")
      ######################
  }


  if (family == "Skew.t"){
      lk <- sum(log(d.mixedST(y, pii, mu, sigma2, shape, nu) ))
      n <- length(y)
      delta <- Delta <- Gama <- rep(0,g)
      for (k in 1:g){
        delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
        Delta[k] <- sqrt(sigma2[k])*delta[k]
        Gama[k] <- sigma2[k] - Delta[k]^2
      }

      teta <- c(mu, Delta, Gama, pii, nu)
      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama

      criterio <- 1
      count <- 0

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
 ##       print(count)
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu[j])/sqrt(sigma2[j]))^2
          Mtij2 <- 1/(1 + (Delta[j]^2)*(Gama[j]^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta[j]*(Gama[j]^(-1))*(y - mu[j])
          A <- mutij / Mtij

          E=(2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2[j])*dt.ls(y, mu[j], sigma2[j],shape[j] ,nu))
          u= ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2[j])*dt.ls(y, mu[j], sigma2[j],shape[j] ,nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)

          d1 <- dt.ls(y, mu[j], sigma2[j], shape[j], nu)
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <-d.mixedST(y, pii, mu, sigma2, shape, nu)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2
          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)


          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])
          mu[j] <- sum(S1[,j]*y - Delta.old[j]*S2[,j]) / sum(S1[,j])
          Delta[j] <- sum(S2[,j]*(y - mu[j])) / sum(S3[,j])
          Gama[j] <- sum(S1[,j]*(y - mu[j])^2 - 2*(y - mu[j])*Delta[j]*S2[,j] + Delta[j]^2*S3[,j]) / sum(tal[,j])       
          sigma2[j] <- Gama[j] + Delta[j]^2
          shape[j] <- ((sigma2[j]^(-1/2))*Delta[j] )/(sqrt(1 - (Delta[j]^2)*(sigma2[j]^(-1))))

        }
        #aqui comecam as alteracoes para estimar o valor de nu
        logvero.ST <- function(nu) sum(log(d.mixedST(y, pii, mu, sigma2, shape, nu)))
        nu <- optimize(logvero.ST, c(0,100), tol = 0.000001, maximum = TRUE)$maximum
        lk1 <- sum(log( d.mixedST(y, pii, mu, sigma2, shape, nu) ))
        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

        param <- teta
        teta <- c(mu, Delta, Gama, pii, nu)
        #criterio <- sqrt((teta-param)%*%(teta-param))
        criterio <- abs(lk1/lk-1)

        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama
        lk<-lk1
        
      }

    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl <- icl+sum(log(pii[j]*dt.ls(y[cl==j], mu[j], sigma2[j], shape[j], nu)))
       #icl=-2*icl+4*g*log(n)
    }

      #### grafico ajuste
      ##hist(y, breaks = 40,probability=T,col="grey")
      ##xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      ##lines(xx,d.mixedST(xx, pii, mu, sigma2, shape, nu),col="red")
      ######################
  }

  if (family == "Skew.cn"){
      if(length(nu) != 2) stop("The Skew.cn need a vector of two components in nu both between (0,1).")
      if((nu[1] <= 0) || (nu[1] >= 1) || (nu[2] <= 0) || (nu[2] >= 1)) stop("nu component(s) out of range.")

      lk <- sum(log( d.mixedSNC(y, pii, mu, sigma2, shape, nu)))
      n <- length(y)
      delta <- Delta <- Gama <- rep(0,g)
      for (k in 1:g){
        delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
        Delta[k] <- sqrt(sigma2[k])*delta[k]
        Gama[k] <- sigma2[k] - Delta[k]^2
      }


      teta <- c(mu, Delta, Gama, pii, nu)
      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama
      nu.old <- nu

      criterio <- 1
      count <- 0

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
##        print(count)
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu[j])/sqrt(sigma2[j]))^2
          Mtij2 <- 1/(1 + (Delta[j]^2)*(Gama[j]^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta[j]*(Gama[j]^(-1))*(y - mu[j])
          A <- mutij / Mtij

          u=(2/dSNC(y,mu[j],sigma2[j],shape[j],nu))*(nu[1]*nu[2]*dnorm(y,mu[j],sqrt(sigma2[j]/nu[2]))*pnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu[j],sqrt(sigma2[j]))*pnorm(A,0,1))
          E=(2/dSNC(y,mu[j],sigma2[j],shape[j],nu))*(nu[1]*sqrt(nu[2])*dnorm(y,mu[j],sqrt(sigma2[j]/nu[2]))*dnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu[j],sqrt(sigma2[j]))*dnorm(A,0,1))

          d1 <-dSNC(y, mu[j], sigma2[j], shape[j], nu)
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <- d.mixedSNC(y, pii, mu, sigma2, shape, nu)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2
          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)


          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])
          mu[j] <- sum(S1[,j]*y - Delta.old[j]*S2[,j]) / sum(S1[,j])
          Delta[j] <- sum(S2[,j]*(y - mu[j])) / sum(S3[,j])
          Gama[j] <- sum(S1[,j]*(y - mu[j])^2 - 2*(y - mu[j])*Delta[j]*S2[,j] + Delta[j]^2*S3[,j]) / sum(tal[,j])
          sigma2[j] <- Gama[j] + Delta[j]^2
          shape[j] <- ((sigma2[j]^(-1/2))*Delta[j] )/(sqrt(1 - (Delta[j]^2)*(sigma2[j]^(-1))))
        }

        #aqui comecam as alteracoes para estimar o valor de nu
        logvero.SNC <- function(nu) sum(log( d.mixedSNC(y, pii, mu, sigma2, shape, nu) ))
        nu <- optim(nu.old, logvero.SNC, control = list(fnscale = -1), method = "L-BFGS-B", lower = rep(0.01, 2), upper = rep(0.99,2))$par

        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

        param <- teta
        teta <- c(mu, Delta, Gama, pii, nu)
        lk1 <- sum(log( d.mixedSNC(y, pii, mu, sigma2, shape, nu) ))
        #criterio <- sqrt((teta-param)%*%(teta-param))
        criterio <- abs(lk1/lk-1)
        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama
        nu.old <- nu
        lk<-lk1
      }

    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl <- icl+sum(log(pii[j]*dSNC(y[cl==j], mu[j], sigma2[j], shape[j], nu)))
       #icl <- -2*icl+(4*g+1)*log(n)
    }
     
      #### grafico ajuste
      ##hist(y, breaks = 40,probability=T,col="grey")
      ##xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      ##lines(xx,d.mixedSNC(xx, pii, mu, sigma2, shape, nu),col="red")
      ######################
      
  }

  if (family == "Skew.slash"){
      cat("\nThe Skew.slash can take a long time to run.\n")
      lk <- sum(log( d.mixedSS(y, pii, mu, sigma2, shape, nu) ))
      n <- length(y)
      delta <- Delta <- Gama <- rep(0,g)
      for (k in 1:g){
        delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
        Delta[k] <- sqrt(sigma2[k])*delta[k]
        Gama[k] <- sigma2[k] - Delta[k]^2
      }

      teta <- c(mu, Delta, Gama, pii, nu)
      #teta <- c(mu, Delta, Gama, pii)
      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama

      criterio <- 1
      count <- 0

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
 ##       print(count)
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          u <- vector(mode="numeric",length=n)
          E <- vector(mode="numeric",length=n)

          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu[j])/sqrt(sigma2[j]))^2
          Mtij2 <- 1/(1 + (Delta[j]^2)*(Gama[j]^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta[j]*(Gama[j]^(-1))*(y - mu[j])
          A <- mutij / Mtij

          for(i in 1:n){
            #U <- runif(2500)
            #V <- pgamma(1,(2*nu + 3)/2, dj[i]/2)*U
            #S <- qgamma(V,(2*nu + 3)/2, dj[i]/2)
            #u[i] <- ( ((2^(nu + 2))*nu*gamma((2*nu + 3)/2)) / (dSS(y[i], mu[j], sigma2[j], shape[j], nu)*sqrt(pi)*sqrt(sigma2[j]))) * pgamma(1,(2*nu + 3)/2,dj[i]/2)*(dj[i]^(-(2*nu + 3)/2))*mean(pnorm(S^(1/2)*A[i]))
            E[i] <- (((2^(nu + 1))*nu*gamma(nu + 1))/(dSS(y[i], mu[j], sigma2[j], shape[j], nu)*pi*sqrt(sigma2[j])))* ((dj[i]+A[i]^2)^(-nu-1))*pgamma(1,nu+1,(dj[i]+A[i]^2)/2)
            faux <- function(u) u^(nu+0.5)*exp(-u*dj[i]/2)*pnorm(u^(1/2)*A[i])
            aux22 <- integrate(faux,0,1)$value
            u[i] <- ((sqrt(2)*nu) / (dSS(y[i], mu[j], sigma2[j], shape[j], nu)*sqrt(pi)*sqrt(sigma2[j])))*aux22  
          }
          #E <- (((2^(nu + 1))*nu*gamma(nu + 1))/(dSS(y, mu[j], sigma2[j], shape[j], nu)*pi*sqrt(sigma2[j])))*((dj+A^2)^(-nu-1))*pgamma(1,nu+1,(dj+A^2)/2)
          d1 <- dSS(y, mu[j], sigma2[j], shape[j], nu)
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <- d.mixedSS(y, pii, mu, sigma2, shape, nu)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2

          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)


          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])
          mu[j] <- sum(S1[,j]*y - Delta.old[j]*S2[,j]) / sum(S1[,j])
          Delta[j] <- sum(S2[,j]*(y - mu[j])) / sum(S3[,j])
          Gama[j] <- sum(S1[,j]*(y - mu[j])^2 - 2*(y - mu[j])*Delta[j]*S2[,j] + Delta[j]^2*S3[,j]) / sum(tal[,j])
          sigma2[j] <- Gama[j] + Delta[j]^2
          shape[j] <- ((sigma2[j]^(-1/2))*Delta[j] )/(sqrt(1 - (Delta[j]^2)*(sigma2[j]^(-1))))
        }
        #aqui comecam as alteracoes para estimar o valor de nu
        logvero.SS <- function(nu) sum(log( d.mixedSS(y, pii, mu, sigma2, shape, nu) ))
        nu <- optimize(logvero.SS, c(0,100), tol = 0.000001, maximum = TRUE)$maximum
        lk1 <- sum(log( d.mixedSS(y, pii, mu, sigma2, shape, nu) ))
        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

        param <- teta
        teta <- c(mu, Delta, Gama, pii, nu)
        #teta <- c(mu, Delta, Gama, pii)
        #criterio <- sqrt((teta-param)%*%(teta-param))
        criterio<-abs(lk1/lk-1)

        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama
        lk<-lk1
        #print(criterio)
      }

    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl <- icl+sum(log(pii[j]*dSS(y[cl==j], mu[j], sigma2[j], shape[j], nu)))
       #icl <- -2*icl+4*g*log(n)
    }
      
      #### grafico ajuste
      ##hist(y, breaks = 40,probability=T,col="grey")
      ##xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      ##lines(xx,d.mixedSS(xx, pii, mu, sigma2, shape, nu),col="red")
      ######################
      
  }
  
  if (family == "Skew.normal"){
          lk <- sum(log( d.mixedSN(y, pii, mu, sigma2, shape) ))
      n <- length(y)
      delta <- Delta <- Gama <- rep(0,g)
      for (k in 1:g){
        delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
        Delta[k] <- sqrt(sigma2[k])*delta[k]
        Gama[k] <- sigma2[k] - Delta[k]^2
      }

      teta <- c(mu, Delta, Gama, pii)
      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama

      criterio <- 1
      count <- 0

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
   ##     print(count)
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          ### E-step: calculando ui, tui, tui2 ###
          
          Mtij2 <- 1/(1 + (Delta[j]^2)*(Gama[j]^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta[j]*(Gama[j]^(-1))*(y - mu[j])

          prob <- pnorm(mutij/Mtij)
          if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin
          E = dnorm(mutij/Mtij) / prob
          u = rep(1, n)

          
          d1 <- dSN(y, mu[j], sigma2[j], shape[j])
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <- d.mixedSN(y, pii, mu, sigma2, shape)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2

          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])
          mu[j] <- sum(S1[,j]*y - Delta.old[j]*S2[,j]) / sum(S1[,j])
          Delta[j] <- sum(S2[,j]*(y - mu[j])) / sum(S3[,j])
          Gama[j] <- sum(S1[,j]*(y - mu[j])^2 - 2*(y - mu[j])*Delta[j]*S2[,j] + Delta[j]^2*S3[,j]) / sum(tal[,j])
          sigma2[j] <- Gama[j] + Delta[j]^2
          shape[j] <- ((sigma2[j]^(-1/2))*Delta[j] )/(sqrt(1 - (Delta[j]^2)*(sigma2[j]^(-1))))

        }
        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

        param <- teta
        teta <- c(mu, Delta, Gama, pii)
        lk1 <- sum(log( d.mixedSN(y, pii, mu, sigma2, shape) ))
       # criterio <- sqrt((teta-param)%*%(teta-param))
        criterio <- abs(lk1/lk-1)

        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama
        lk<-lk1

      }

    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl<-icl+sum(log(pii[j]*dSN(y[cl==j], mu[j], sigma2[j], shape[j])))
       #icl <- -2*icl+(4*g-1)*log(n)
    }

      #### grafico ajuste
      ##hist(y, breaks = 40,probability=T,col="grey")
      ##xx=seq(min(y),max(y),(max(y)-min(y))/1000)
      ##lines(xx,d.mixedSN(xx, pii, mu, sigma2, shape),col="red")
      ######################
  }
  
  if (family == "Normal"){
      shape <- rep(0,g) 
      lk <- sum(log( d.mixedSN(y, pii, mu, sigma2, shape) ))
      n <- length(y)
      delta <- Delta <- Gama <- rep(0,g)

      for (k in 1:g){
        delta[k] <- 0 ##shape[k] / (sqrt(1 + shape[k]^2))
        Delta[k] <- 0 ##sqrt(sigma2[k])*delta[k]
        Gama[k] <- sigma2[k] - Delta[k]^2
      }

      teta <- c(mu, Delta, Gama, pii)
      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama

      criterio <- 1
      count <- 0

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          ### E-step: calculando ui, tui, tui2 ###
          Mtij2 <- 1/(1 + (Delta[j]^2)*(Gama[j]^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta[j]*(Gama[j]^(-1))*(y - mu[j])

          prob <- pnorm(mutij/Mtij)
          if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin
          E = dnorm(mutij/Mtij) / prob
          u = rep(1, n)
 
          d1 <- dSN(y, mu[j], sigma2[j], shape[j])
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <- d.mixedSN(y, pii, mu, sigma2, shape)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2
          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])
          mu[j] <- sum(S1[,j]*y - Delta.old[j]*S2[,j]) / sum(S1[,j])
          Delta[j] <- 0          
          Gama[j] <- sum(S1[,j]*(y - mu[j])^2 - 2*(y - mu[j])*Delta[j]*S2[,j] + Delta[j]^2*S3[,j]) / sum(tal[,j])
          sigma2[j] <- Gama[j] + Delta[j]^2
          shape[j] <- 0
        }
        
        pii[g] <- 1 - (sum(pii) - pii[g])
        
        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }
             
        param <- teta
        teta <- c(mu, Delta, Gama, pii)
        lk1 <- sum(log(d.mixedSN(y, pii, mu, sigma2, shape)))
        #criterio <- sqrt((teta-param)%*%(teta-param))
        criterio <- abs(lk1/lk-1)
        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama
        lk<-lk1
      }       

    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl <- icl+sum(log(pii[j]*dSN(y[cl==j], mu[j], sigma2[j], shape[j])))
       #icl <- -2*icl+(3*g-1)*log(n)
    }

    #### grafico ajuste
    ##hist(y, breaks = 40,probability=T,col="grey")
    ##xx=seq(min(y),max(y),(max(y)-min(y))/1000)
    ##lines(xx,d.mixedSN(xx, pii, mu, sigma2, shape),col="red")
    ######################

  }


  if (criteria == TRUE){       
    if((family == "t") | (family == "Normal")) d <- g*2 + (g-1) #mu + Sigma + pi
    else d <- g*3 + (g-1) #mu + Sigma + pi
    aic <- -2*lk + 2*d
    bic <- -2*lk + log(n)*d
    edc <- -2*lk + 0.2*sqrt(n)*d
    icl <- -2*icl + log(n)*d
    obj.out <- list(mu = mu, sigma2 = sigma2, shape = shape, pii = pii, nu = nu, aic = aic, bic = bic, edc = edc, icl= icl, iter = count, n = length(y), group = cl)#, im.sdev=NULL)
  }
 else obj.out <- list(mu = mu, sigma2 = sigma2, shape = shape, pii = pii, nu = nu, iter = count, n = length(y), group = apply(tal, 1, which.max))#, im.sdev=NULL)
# if (group == FALSE) obj.out <- obj.out[-(length(obj.out)-1)]
 if (group == FALSE) obj.out <- obj.out[-(length(obj.out))]
 if (obs.prob == TRUE){
     nam <- c()
     for (i in 1:ncol(tal)) nam <- c(nam,paste("Group ",i,sep=""))
     if(ncol(tal) == 1) dimnames(tal)[[2]] <- list(nam)
     if(ncol(tal) > 1) dimnames(tal)[[2]] <- nam
     obj.out$obs.prob <- tal
     if((ncol(tal) - 1) > 1) obj.out$obs.prob[,ncol(tal)] <- 1 - rowSums(obj.out$obs.prob[,1:(ncol(tal)-1)])
     else obj.out$obs.prob[,ncol(tal)] <- 1 - obj.out$obs.prob[,1]
     obj.out$obs.prob[which(obj.out$obs.prob[,ncol(tal)] <0),ncol(tal)] <- 0.0000000000
     obj.out$obs.prob <- round(obj.out$obs.prob,10)
 }
 class(obj.out) <- family

 if(calc.im){
   IM <-  im.smsn(as.matrix(y), obj.out)
   sdev <- sqrt(diag(solve(IM[[1]])))
   obj.out$im.sdev = sdev
  }
 
# else obj.out <- obj.out[-length(obj.out)] 
 
 class(obj.out) <- family
 obj.out
}


##########     FIM DO ALGORITMO EM PARA MISTURAS      ###############
#####################################################################

