##############################################################################################
###########     ALGORITMO EM MULTIVARIADO PARA MISTURAS DE DISTRI SMSN        ################
###########     Alteracoes feitas em 03/05/2010 por Celso Romulo              ################
###########     Esta eh a versao utilizada no artigo submetido ao CSDA         ################
smsn.mmix <- function(y, nu=1, mu = NULL, Sigma = NULL, shape = NULL, pii = NULL, g = NULL, get.init = TRUE, criteria = TRUE,
                      group = FALSE, family = "Skew.normal", error = 0.0001, iter.max = 100, uni.Gama = FALSE, calc.im=FALSE,
                      obs.prob= FALSE, kmeans.param = NULL){
  #mu, Sigma, shape devem ser do tipo list(). O numero de entradas no list eh o numero g de componentes de misturas
  #cada entrada do list deve ser de tamanho igual ao numero de colunas da matriz de dados y
  y <- as.matrix(y)
  dimnames(y) <- NULL
  if(!is.matrix(y)) stop("The response is not in a matrix format\n")
  if(ncol(y) <= 1) stop("For the univariate case use the smsn.mix function\n")
  if((family != "t") && (family != "Skew.t") && (family != "Skew.cn") && (family != "Skew.slash") && (family != "Skew.normal") && (family != "Normal")) stop(paste("Family",family,"not recognized.",sep=" "))
  if((length(g) == 0) && ((length(mu)==0) || (length(Sigma)==0) || (length(shape)==0) || (length(pii)==0)))  stop("The model is not specified correctly.\n")
  if(get.init == FALSE){
       g <- length(mu)
       if((length(mu) != length(Sigma)) || (length(mu) != length(pii))) stop("The size of the initial values are not compatibles.\n")
       if((family == "Skew.t" || family == "Skew.cn" || family == "Skew.slash" || family == "Skew.normal") & (length(mu) != length(shape))) stop("The size of the initial values are not compatibles.\n")

       if(sum(pii) != 1) stop("probability of pii does not sum to 1.\n")
       for (j in 1:g){
        dimnames(Sigma[[j]]) <- NULL
        names(mu[[j]]) <- NULL
        names(shape[[j]]) <- NULL
      }
  }
  if((length(g)!= 0) && (g < 1)) stop("g must be greater than 0.\n")

  p <- ncol(y)
  n <- nrow(y)

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
      mu <- shape <- Sigma <- list()
      for (j in 1:g){
        mu[[j]] <- init$centers[j,]
        shape[[j]] <- sign( apply( (y[init$cluster == j, ] - matrix(rep(mu[[j]], nrow(y[init$cluster == j, ])), 
                      nrow = nrow(y[init$cluster == j, ]), ncol=p, byrow = TRUE))^3, 2, sum))
        Sigma[[j]] <- var(y[init$cluster == j,])
        dimnames(Sigma[[j]]) <- NULL
        names(mu[[j]]) <- NULL
        names(shape[[j]]) <- NULL
      }
    }
    
    else{
      pii <- 1
      mu[[1]] <- colMeans(y)
      Sigma[[1]] <- var(y)
      shape[[1]] <- sign( apply( (y - matrix(rep(mu[[1]], nrow(y)), nrow(y), p, byrow = TRUE) )^3, 2, sum ))
      dimnames(Sigma[[1]]) <- NULL
      names(mu[[1]]) <- NULL
      names(shape[[1]]) <- NULL
    }

  }

  if (family == "t"){
      delta <- Delta <- Gama <- list()
#      teta <- c()
      for (k in 1:g){
        Delta[[k]] <- shape[[k]] <- rep(0,p)
        Gama[[k]] <- Sigma[[k]] - Delta[[k]]%*%t(Delta[[k]])
#        teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)], Delta[[k]])
      }
#      teta <- c(teta, pii, nu) #teta para nu desconhecido
      if(uni.Gama){
         Gama.uni <- Gama[[1]]
         if(g > 1) for(k in 2:g) Gama.uni <- Gama.uni + Gama[[k]]
         Gama.uni <- Gama.uni / g
         for(k in 1:g) Gama[[k]] <- Gama.uni
      }
      
      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama

      criterio <- 1
      count <- 0
      lkante <- 1

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          ### E-step: calculando ui, tui, tui2 ###
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Mtij2 <- as.numeric( 1/(1 + t(Delta[[j]])%*%solve(Gama[[j]])%*%Delta[[j]]) )
          Mtij <- sqrt(Mtij2)
          mutij <- apply(  matrix(rep( Mtij2*t(Delta[[j]])%*%solve(Gama[[j]]),n), n, p, byrow = TRUE ) * Dif, 1, sum)
          A <- mutij / Mtij
          dj <- mahalanobis(y, mu[[j]], Sigma[[j]])

          E = (2*(nu)^(nu/2)*gamma((p+nu+1)/2)*((dj + nu + A^2))^(-(p+nu+1)/2)) / (gamma(nu/2)*(sqrt(pi))^(p+1)*sqrt(det(Sigma[[j]]))*dmvt.ls(y, mu[[j]], Sigma[[j]], shape[[j]] ,nu))
          u = ((4*(nu)^(nu/2)*gamma((p+nu+2)/2)*(dj + nu)^(-(p+nu+2)/2)) / (gamma(nu/2)*sqrt(pi^p)*sqrt(det(Sigma[[j]]))*dmvt.ls(y, mu[[j]], Sigma[[j]], shape[[j]] ,nu)) )*pt(sqrt((p+nu+2)/(dj+nu))*A,p+nu+2)

          d1 <- dmvt.ls(y, mu[[j]], Sigma[[j]], shape[[j]], nu)
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <- d.mixedmvST(y, pii, mu, Sigma, shape, nu)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2
          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])

          mu[[j]] <- apply(S1[,j]*y - S2[,j] * matrix(rep(Delta.old[[j]], n), n, p, byrow = TRUE), 2, sum) / sum(S1[,j])
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Delta[[j]] <- rep(0,p)
          
          sum2 <- matrix(0,p,p)
          for (i in 1:n) sum2 <- sum2 + (S1[i,j]*(y[i,] - mu[[j]])%*%t(y[i,] - mu[[j]]) - S2[i,j]*Delta[[j]]%*%t(y[i,] - mu[[j]]) - S2[i,j]*(y[i,] - mu[[j]])%*%t(Delta[[j]])  + S3[i,j]*Delta[[j]]%*%t(Delta[[j]]))

          Gama[[j]] <- sum2 / sum(tal[,j])
          
          if(!uni.Gama){
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- rep(0,p)
          }
        }
        if(uni.Gama){
           GS <- 0
           for (j in 1:g) GS <- GS+kronecker(tal[,j],Gama[[j]])
           Gama.uni <- t(rowSums(array(t(GS),dim=c(p,p,n)),dims=2))/n
           for (j in 1:g){
             Gama[[j]] <- Gama.uni
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- rep(0,p)
           }
        }
        #aqui comecam as alteracoes para estimar o valor de nu
        logvero.ST <- function(nu) sum(log( d.mixedmvST(y, pii, mu, Sigma, shape, nu) ))
        nu <- optimize(logvero.ST, c(0,100), tol = 0.000001, maximum = TRUE)$maximum

        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

#        param <- teta
#        teta <- c()
#        for (k in 1:g) teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)], Delta[[k]])
#        teta <- c(teta, pii, nu) #teta para nu desconhecido
        lk <- sum(log( d.mixedmvST(y, pii, mu, Sigma, shape, nu) ))
        criterio <- abs((lk/lkante)-1)

        lkante <- lk
        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama
      }
    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl <- icl+sum(log(pii[j]*dmvt.ls(y[cl==j,], mu[[j]], Sigma[[j]], shape[[j]], nu)))
       #icl=-2*icl+4*g*log(n)
    }

  }


  if (family == "Skew.t"){
      delta <- Delta <- Gama <- list()
#      teta <- c()
      for (k in 1:g){
        delta[[k]] <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
        Delta[[k]] <- as.vector(matrix.sqrt(Sigma[[k]])%*%delta[[k]])
        Gama[[k]] <- Sigma[[k]] - Delta[[k]]%*%t(Delta[[k]])
#        teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)], Delta[[k]])
      }
#      teta <- c(teta, pii, nu) #teta para nu desconhecido
      if(uni.Gama){
         Gama.uni <- Gama[[1]]
         if(g > 1) for(k in 2:g) Gama.uni <- Gama.uni + Gama[[k]]
         Gama.uni <- Gama.uni / g
         for(k in 1:g) Gama[[k]] <- Gama.uni
      }
      
      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama

      criterio <- 1
      count <- 0
      lkante <- 1

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          ### E-step: calculando ui, tui, tui2 ###
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Mtij2 <- as.numeric( 1/(1 + t(Delta[[j]])%*%solve(Gama[[j]])%*%Delta[[j]]) )
          Mtij <- sqrt(Mtij2)
          mutij <- apply(  matrix(rep( Mtij2*t(Delta[[j]])%*%solve(Gama[[j]]),n), n, p, byrow = TRUE ) * Dif, 1, sum)
          A <- mutij / Mtij
          dj <- mahalanobis(y, mu[[j]], Sigma[[j]])

          E = (2*(nu)^(nu/2)*gamma((p+nu+1)/2)*((dj + nu + A^2))^(-(p+nu+1)/2)) / (gamma(nu/2)*(sqrt(pi))^(p+1)*sqrt(det(Sigma[[j]]))*dmvt.ls(y, mu[[j]], Sigma[[j]], shape[[j]] ,nu))
          u = ((4*(nu)^(nu/2)*gamma((p+nu+2)/2)*(dj + nu)^(-(p+nu+2)/2)) / (gamma(nu/2)*sqrt(pi^p)*sqrt(det(Sigma[[j]]))*dmvt.ls(y, mu[[j]], Sigma[[j]], shape[[j]] ,nu)) )*pt(sqrt((p+nu+2)/(dj+nu))*A,p+nu+2)

          d1 <- dmvt.ls(y, mu[[j]], Sigma[[j]], shape[[j]], nu)
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <- d.mixedmvST(y, pii, mu, Sigma, shape, nu)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2
          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])

          mu[[j]] <- apply(S1[,j]*y - S2[,j] * matrix(rep(Delta.old[[j]], n), n, p, byrow = TRUE), 2, sum) / sum(S1[,j])
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Delta[[j]] <- apply(S2[,j]*Dif, 2, sum) / sum(S3[,j])
          
          sum2 <- matrix(0,p,p)
          for (i in 1:n) sum2 <- sum2 + (S1[i,j]*(y[i,] - mu[[j]])%*%t(y[i,] - mu[[j]]) - S2[i,j]*Delta[[j]]%*%t(y[i,] - mu[[j]]) - S2[i,j]*(y[i,] - mu[[j]])%*%t(Delta[[j]])  + S3[i,j]*Delta[[j]]%*%t(Delta[[j]]))

          Gama[[j]] <- sum2 / sum(tal[,j])
          
          if(!uni.Gama){
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- (solve(matrix.sqrt(Sigma[[j]]))%*%Delta[[j]]) / as.numeric( (1 - t(Delta[[j]])%*%solve(Sigma[[j]])%*%Delta[[j]]) )^(1/2)
          }
        }
        if(uni.Gama){
           GS <- 0
           for (j in 1:g) GS <- GS+kronecker(tal[,j],Gama[[j]])
           Gama.uni <- t(rowSums(array(t(GS),dim=c(p,p,n)),dims=2))/n
           for (j in 1:g){
             Gama[[j]] <- Gama.uni
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- (solve(matrix.sqrt(Sigma[[j]]))%*%Delta[[j]]) / as.numeric( (1 - t(Delta[[j]])%*%solve(Sigma[[j]])%*%Delta[[j]]) )^(1/2)
           }
        }
        #aqui comecam as alteracoes para estimar o valor de nu
        logvero.ST <- function(nu) sum(log( d.mixedmvST(y, pii, mu, Sigma, shape, nu) ))
        nu <- optimize(logvero.ST, c(0,100), tol = 0.000001, maximum = TRUE)$maximum

        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

#        param <- teta
#        teta <- c()
#        for (k in 1:g) teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)], Delta[[k]])
#        teta <- c(teta, pii, nu) #teta para nu desconhecido
        lk <- sum(log( d.mixedmvST(y, pii, mu, Sigma, shape, nu) ))
        criterio <- abs((lk/lkante)-1)

        lkante <- lk
        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama
      }
    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl <- icl+sum(log(pii[j]*dmvt.ls(y[cl==j,], mu[[j]], Sigma[[j]], shape[[j]], nu)))
       #icl=-2*icl+4*g*log(n)
    }

  }

  if (family == "Skew.cn"){
      if(length(nu) != 2) stop("The Skew.cn need a vector of two components in nu both between (0,1).")
      delta <- Delta <- Gama <- list()
#      teta <- c()
      for (k in 1:g){
        delta[[k]] <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
        Delta[[k]] <- as.vector(matrix.sqrt(Sigma[[k]])%*%delta[[k]])
        Gama[[k]] <- Sigma[[k]] - Delta[[k]]%*%t(Delta[[k]])    
##        teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)], Delta[[k]])
      }
##      teta <- c(teta, pii, nu) #teta para nu desconhecido
      if(uni.Gama){
         Gama.uni <- Gama[[1]]
         if(g > 1) for(k in 2:g) Gama.uni <- Gama.uni + Gama[[k]]
         Gama.uni <- Gama.uni / g
         for(k in 1:g) Gama[[k]] <- Gama.uni
      }
      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama
      nu.old <- nu

      criterio <- 1
      count <- 0
      lkante <- 1

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          ### E-step: calculando ui, tui, tui2 ###
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Mtij2 <- as.numeric( 1/(1 + t(Delta[[j]])%*%solve(Gama[[j]])%*%Delta[[j]]) )
          Mtij <- sqrt(Mtij2)
          mutij <- apply(  matrix(rep( Mtij2*t(Delta[[j]])%*%solve(Gama[[j]]),n), n, p, byrow = TRUE ) * Dif, 1, sum)
          A <- mutij / Mtij
          dj <- mahalanobis(y, mu[[j]], Sigma[[j]])

          u = (2/dmvSNC(y,mu[[j]],Sigma[[j]],shape[[j]],nu))*(nu[1]*nu[2]*dmvnorm(y,mu[[j]],Sigma[[j]]/nu[2])*pnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dmvnorm(y,mu[[j]],Sigma[[j]])*pnorm(A,0,1))
          E = (2/dmvSNC(y,mu[[j]],Sigma[[j]],shape[[j]],nu))*(nu[1]*sqrt(nu[2])*dmvnorm(y,mu[[j]],Sigma[[j]]/nu[2])*dnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dmvnorm(y,mu[[j]],Sigma[[j]])*dnorm(A,0,1))

          d1 <-dmvSNC(y, mu[[j]], Sigma[[j]], shape[[j]], nu)
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <- d.mixedmvSNC(y, pii, mu, Sigma, shape, nu)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2
          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])
          
          mu[[j]] <- apply(S1[,j]*y - S2[,j] * matrix(rep(Delta.old[[j]], n), n, p, byrow = TRUE), 2, sum) / sum(S1[,j])
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Delta[[j]] <- apply(S2[,j]*Dif, 2, sum) / sum(S3[,j])
          
          sum2 <- matrix(0,p,p)
          for (i in 1:n) sum2 <- sum2 + (S1[i,j]*(y[i,] - mu[[j]])%*%t(y[i,] - mu[[j]]) - S2[i,j]*Delta[[j]]%*%t(y[i,] - mu[[j]]) - S2[i,j]*(y[i,] - mu[[j]])%*%t(Delta[[j]])  + S3[i,j]*Delta[[j]]%*%t(Delta[[j]]))

          Gama[[j]] <- sum2 / sum(tal[,j])
          
          if(!uni.Gama){
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- (solve(matrix.sqrt(Sigma[[j]]))%*%Delta[[j]]) / as.numeric( (1 - t(Delta[[j]])%*%solve(Sigma[[j]])%*%Delta[[j]]) )^(1/2)
          }
        }
        if(uni.Gama){
           GS <- 0
           for (j in 1:g) GS <- GS+kronecker(tal[,j],Gama[[j]])
           Gama.uni <- t(rowSums(array(t(GS),dim=c(p,p,n)),dims=2))/n
           for (j in 1:g){
             Gama[[j]] <- Gama.uni
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- (solve(matrix.sqrt(Sigma[[j]]))%*%Delta[[j]]) / as.numeric( (1 - t(Delta[[j]])%*%solve(Sigma[[j]])%*%Delta[[j]]) )^(1/2)
           }
        }
        #aqui comecam as alteracoes para estimar o valor de nu
        logvero.SNC <- function(nu) sum(log( d.mixedmvSNC(y, pii, mu, Sigma, shape, nu) ))
        nu <- optim(nu.old, logvero.SNC, control = list(fnscale = -1), method = "L-BFGS-B", lower = rep(0.01, 2), upper = rep(0.99,2))$par

        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

##        param <- teta
##        teta <- c()
##        for (k in 1:g) teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)], Delta[[k]])
##        teta <- c(teta, pii, nu) #teta para nu desconhecido
        lk <- sum(log( d.mixedmvSNC(y, pii, mu, Sigma, shape, nu) ))
        criterio <- abs((lk/lkante)-1)

        lkante <- lk
        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama
        nu.old <- nu
      }
    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl <- icl+sum(log(pii[j]*dmvSNC(y[cl==j,], mu[[j]], Sigma[[j]], shape[[j]], nu)))
       #icl <- -2*icl+(4*g+1)*log(n)
    }

  }


    if (family == "Skew.slash"){
      cat("\nThe Skew.slash takes a long time to run.\n")
      delta <- Delta <- Gama <- list()
#      teta <- c()
      for (k in 1:g){
        delta[[k]] <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
        Delta[[k]] <- as.vector(matrix.sqrt(Sigma[[k]])%*%delta[[k]])
        Gama[[k]] <- Sigma[[k]] - Delta[[k]]%*%t(Delta[[k]])
#        teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)], Delta[[k]])
      }
      #teta <- c(teta, pii, nu) #teta para nu desconhecido
#      teta <- c(teta, pii) #teta para nu desconhecido
      if(uni.Gama){
         Gama.uni <- Gama[[1]]
         if(g > 1) for(k in 2:g) Gama.uni <- Gama.uni + Gama[[k]]
         Gama.uni <- Gama.uni / g
         for(k in 1:g) Gama[[k]] <- Gama.uni
      }

      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama

      criterio <- 1
      count <- 0
      lkante <- 1

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          ### E-step: calculando ui, tui, tui2 ###
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Mtij2 <- as.numeric( 1/(1 + t(Delta[[j]])%*%solve(Gama[[j]])%*%Delta[[j]]) )
          Mtij <- sqrt(Mtij2)
          mutij <- apply(  matrix(rep( Mtij2*t(Delta[[j]])%*%solve(Gama[[j]]),n), n, p, byrow = TRUE ) * Dif, 1, sum)
          A <- mutij / Mtij
          dj <- mahalanobis(y, mu[[j]], Sigma[[j]])

          #int <- c()
          #for (k in 1:n){ #melhorar se der
          #  f <- function(u) pnorm(u^(1/2)*A[k])*pgamma(u,(2*nu+p+2)/2,dj[k]/2)
          #  int[k] <- integrate(f,0,1)$value
          #}

          #u = ((2^(nu+2)*nu*gamma((2*nu+p+2)/(2)))/(dmvSS(y, mu[[j]], Sigma[[j]], shape[[j]], nu)*sqrt(pi^p)*det(Sigma[[j]])^(1/2)))*pgamma(1,(2*nu+p+2)/2,dj/2)*(dj^(-(2*nu+p+2)/2))*int
          #E = ((2^(nu+1)*nu*gamma((2*nu+p+1)/(2)))/(dmvSS(y, mu[[j]], Sigma[[j]], shape[[j]], nu)*sqrt(pi)^(p+1)*det(Sigma[[j]])^(1/2)))*(dj + A^2)^(-(2*nu+p+1)/(2))*pgamma(1,(2*nu+p+1)/(2),(dj+A^2)/2)

          u <- E <- c()
#===========================================================================================================
# Usando integracao de Monte Carlo para calcular kappa_1
#          for(i in 1:n){
#            U <- runif(2500)
#            V <- pgamma(1,(2*nu + 3)/2, dj[i]/2)*U
#            S <- qgamma(V,(2*nu + 3)/2, dj[i]/2)
#            u[i] = ((2^(nu+2)*nu*gamma((2*nu+p+2)/(2)))/
#                   (dmvSS(y[i,], mu[[j]], Sigma[[j]], shape[[j]], nu)*sqrt(pi^p)*det(Sigma[[j]])^(1/2)))*
#                   pgamma(1,(2*nu+p+2)/2,dj[i]/2)*(dj[i]^(-(2*nu+p+2)/2))*mean(pnorm(S^(1/2)*A[i]))
#          }
#============================================================================================================
#===========================================================================================================
# Usando integracao numerica para calcular kappa_1 (Introduzido por Celso Romulo)
          for(i in 1:n){
            faux <- function(u) u^(nu+p/2)*exp(-u*dj[i]/2)*pnorm(u^(1/2)*A[i])
            aux22 <- integrate(faux,0,1)$value
            u[i] = nu*2^(1-p/2)*pi^(-0.5*p)*det(Sigma[[j]])^(-1/2)*aux22/dmvSS(y[i,], mu[[j]], Sigma[[j]], shape[[j]], nu)
          }
#            faux <- function(u) u^(nu+p/2)*exp(-u*dj/2)*pnorm(u^(1/2)*A)
#            aux22 <- integrate(faux,0,1)$value
#            u = nu*2^(1-p/2)*pi^(-0.5*p)*det(Sigma[[j]])^(-1/2)*aux22/dmvSS(y, mu[[j]], Sigma[[j]], shape[[j]], nu)
#============================================================================================================
          E = ((2^(nu+1)*nu*gamma((2*nu+p+1)/(2)))/(dmvSS(y, mu[[j]], Sigma[[j]], shape[[j]], nu)*sqrt(pi)^(p+1)*det(Sigma[[j]])^(1/2)))*
              (dj + A^2)^(-(2*nu+p+1)/(2))*pgamma(1,(2*nu+p+1)/(2),(dj+A^2)/2)

          d1 <- dmvSS(y, mu[[j]], Sigma[[j]], shape[[j]], nu)
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <- d.mixedmvSS(y, pii, mu, Sigma, shape, nu)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2
          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])
 
          mu[[j]] <- apply(S1[,j]*y - S2[,j] * matrix(rep(Delta.old[[j]], n), n, p, byrow = TRUE), 2, sum) / sum(S1[,j])
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Delta[[j]] <- apply(S2[,j]*Dif, 2, sum) / sum(S3[,j])
          
          sum2 <- matrix(0,p,p)
          for (i in 1:n) sum2 <- sum2 + (S1[i,j]*(y[i,] - mu[[j]])%*%t(y[i,] - mu[[j]]) - S2[i,j]*Delta[[j]]%*%t(y[i,] - mu[[j]]) - S2[i,j]*(y[i,] - mu[[j]])%*%t(Delta[[j]])  + S3[i,j]*Delta[[j]]%*%t(Delta[[j]]))

          Gama[[j]] <- sum2 / sum(tal[,j])

#          Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
######## DISCUTIR ISSO AQUI
#          if (any(!is.finite(Sigma[[j]])))  {
#             print("infinite or missing values in Sigma")
             #print(Sigma[[j]])
#             obj.out <- list(Sigma=Sigma[[j]],iter = count,control=1)
#             return(obj.out)
#          }
          if(!uni.Gama){
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- (solve(matrix.sqrt(Sigma[[j]]))%*%Delta[[j]]) / as.numeric( (1 - t(Delta[[j]])%*%solve(Sigma[[j]])%*%Delta[[j]]) )^(1/2)
          }
        }
        if(uni.Gama){
           GS <- 0
           for (j in 1:g) GS <- GS+kronecker(tal[,j],Gama[[j]])
           Gama.uni <- t(rowSums(array(t(GS),dim=c(p,p,n)),dims=2))/n
           for (j in 1:g){
             Gama[[j]] <- Gama.uni
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- (solve(matrix.sqrt(Sigma[[j]]))%*%Delta[[j]]) / as.numeric( (1 - t(Delta[[j]])%*%solve(Sigma[[j]])%*%Delta[[j]]) )^(1/2)
           }
        }
        #aqui comecam as alteracoes para estimar o valor de nu
#===============================================
# Ative se desejar estimar nu
        logvero.SS <- function(nu) sum(log( d.mixedmvSS(y, pii, mu, Sigma, shape, nu) ))
        nu <- optimize(logvero.SS, c(0,100), tol = 0.000001, maximum = TRUE)$maximum

        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

#        param <- teta
#        teta <- c()
#        for (k in 1:g) teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)], Delta[[k]])
        #teta <- c(teta, pii, nu) #teta para nu desconhecido
#        teta <- c(teta, pii)
#        criterio <- sqrt((teta-param)%*%(teta-param))
        lk <- sum(log( d.mixedmvSS(y, pii, mu, Sigma, shape, nu) ))
        criterio <- abs((lk/lkante)-1)

        lkante <- lk
        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama

      }
    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl <- icl+sum(log(pii[j]*dmvSS(y[cl==j,], mu[[j]], Sigma[[j]], shape[[j]], nu)))
       #icl <- -2*icl+4*g*log(n)
    }

  }

  if (family == "Skew.normal"){
      delta <- Delta <- Gama <- list()
#      teta <- c()
      for (k in 1:g){
        delta[[k]] <- shape[[k]] / as.numeric(sqrt(1 + t(shape[[k]])%*%shape[[k]]))
        Delta[[k]] <- as.vector(matrix.sqrt(Sigma[[k]])%*%delta[[k]])
        Gama[[k]] <- Sigma[[k]] - Delta[[k]]%*%t(Delta[[k]])
#       teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)], Delta[[k]])
      }
#      teta <- c(teta, pii) #teta para nu conhecido
      if(uni.Gama){
         Gama.uni <- Gama[[1]]
         if(g > 1) for(k in 2:g) Gama.uni <- Gama.uni + Gama[[k]]
         Gama.uni <- Gama.uni / g
         for(k in 1:g) Gama[[k]] <- Gama.uni
      }

      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama

      criterio <- 1
      count <- 0
      lkante <- 1

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          ### E-step: calculando ui, tui, tui2 ###
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Mtij2 <- as.numeric( 1/(1 + t(Delta[[j]])%*%solve(Gama[[j]])%*%Delta[[j]]) )
          Mtij <- sqrt(Mtij2)
          mutij <- apply(  matrix(rep( Mtij2*t(Delta[[j]])%*%solve(Gama[[j]]),n), n, p, byrow = TRUE ) * Dif, 1, sum)
          A <- mutij / Mtij

          prob <- pnorm(A)
          if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin

          E = dnorm(A) / prob
          u = rep(1, n)

          #print(Sigma)
          #print(d.mixedmvSN(y, pii, mu, Sigma, shape))
          #print("foi")

          d1 <- dmvSN(y, mu[[j]], Sigma[[j]], shape[[j]])
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <- d.mixedmvSN(y, pii, mu, Sigma, shape)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2
          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])

          mu[[j]] <- apply(S1[,j]*y - S2[,j] * matrix(rep(Delta.old[[j]], n), n, p, byrow = TRUE), 2, sum) / sum(S1[,j])
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Delta[[j]] <- apply(S2[,j]*Dif, 2, sum) / sum(S3[,j])
          
          sum2 <- matrix(0,p,p)
          for (i in 1:n) sum2 <- sum2 + (S1[i,j]*(y[i,] - mu[[j]])%*%t(y[i,] - mu[[j]]) - S2[i,j]*Delta[[j]]%*%t(y[i,] - mu[[j]]) - S2[i,j]*(y[i,] - mu[[j]])%*%t(Delta[[j]])  + S3[i,j]*Delta[[j]]%*%t(Delta[[j]]))

          Gama[[j]] <- sum2 / sum(tal[,j])

          if(!uni.Gama){
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- (solve(matrix.sqrt(Sigma[[j]]))%*%Delta[[j]]) / as.numeric( (1 - t(Delta[[j]])%*%solve(Sigma[[j]])%*%Delta[[j]]) )^(1/2)
          }
        }
        if(uni.Gama){
           GS <- 0
           for (j in 1:g) GS <- GS+kronecker(tal[,j],Gama[[j]])
           Gama.uni <- t(rowSums(array(t(GS),dim=c(p,p,n)),dims=2))/n
           for (j in 1:g){
             Gama[[j]] <- Gama.uni
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- (solve(matrix.sqrt(Sigma[[j]]))%*%Delta[[j]]) / as.numeric( (1 - t(Delta[[j]])%*%solve(Sigma[[j]])%*%Delta[[j]]) )^(1/2)
           }
        }
        #print("saiu do for g")
        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

#        param <- teta
#        teta <- c()
#        for (k in 1:g) teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)], Delta[[k]])
#        teta <- c(teta, pii) #teta para nu conhecido
        lk <- sum(log( d.mixedmvSN(y, pii, mu, Sigma, shape) ))
        criterio <- abs((lk/lkante)-1)

        lkante <- lk
        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama
      }

    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl <- icl+sum(log(pii[j]*dmvSN(y[cl==j,], mu[[j]], Sigma[[j]], shape[[j]])))
       #icl <- -2*icl+(4*g-1)*log(n)
    }

  }

  if (family == "Normal"){
      delta <- Delta <- Gama <- list()
#      teta <- c()
      for (k in 1:g){
        Delta[[k]] <- shape[[k]] <- rep(0,p)
        Gama[[k]] <- Sigma[[k]]
#        teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)])
      }
#      teta <- c(teta, pii) #teta para nu conhecido
      if(uni.Gama){
         Gama.uni <- Gama[[1]]
         if(g > 1) for(k in 2:g) Gama.uni <- Gama.uni + Gama[[k]]
         Gama.uni <- Gama.uni / g
         for(k in 1:g) Gama[[k]] <- Gama.uni
      }

      mu.old <- mu
      Delta.old <- Delta
      Gama.old <- Gama

      criterio <- 1
      count <- 0
      lkante <- 1

      while((criterio > error) && (count <= iter.max)){
        count <- count + 1
        tal <- matrix(0, n, g)
        S1 <- matrix(0, n, g)
        S2 <- matrix(0, n, g)
        S3 <- matrix(0, n, g)
        for (j in 1:g){
          ### E-step: calculando ui, tui, tui2 ###
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Mtij2 <- as.numeric( 1/(1 + t(Delta[[j]])%*%solve(Gama[[j]])%*%Delta[[j]]) )
          Mtij <- sqrt(Mtij2)
          mutij <- apply(  matrix(rep( Mtij2*t(Delta[[j]])%*%solve(Gama[[j]]),n), n, p, byrow = TRUE ) * Dif, 1, sum)
          A <- mutij / Mtij

          prob <- pnorm(A)
          if(length(which(prob == 0)) > 0) prob[which(prob == 0)] <- .Machine$double.xmin
          E = dnorm(A) / prob
          u = rep(1, n)

          d1 <- dmvSN(y, mu[[j]], Sigma[[j]], shape[[j]])
          if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin
          d2 <- d.mixedmvSN(y, pii, mu, Sigma, shape)
          if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

          tal[,j] <- d1*pii[j] / d2
          S1[,j] <- tal[,j]*u
          S2[,j] <- tal[,j]*(mutij*u + Mtij*E)
          S3[,j] <- tal[,j]*(mutij^2*u + Mtij2 + Mtij*mutij*E)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
          pii[j] <- (1/n)*sum(tal[,j])
          
          mu[[j]] <- apply(S1[,j]*y - S2[,j] * matrix(rep(Delta.old[[j]], n), n, p, byrow = TRUE), 2, sum) / sum(S1[,j])
          Dif <- y - matrix(rep(mu[[j]], n), n, p, byrow = TRUE)
          Delta[[j]] <- rep(0,p)
          
          sum2 <- matrix(0,p,p)
          for (i in 1:n) sum2 <- sum2 + (S1[i,j]*(y[i,] - mu[[j]])%*%t(y[i,] - mu[[j]]))

 
          Gama[[j]] <- sum2 / sum(tal[,j])
          
          if(!uni.Gama){
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- rep(0,p)
          }
        }
        if(uni.Gama){
           GS <- 0
           for (j in 1:g) GS <- GS+kronecker(tal[,j],Gama[[j]])
           Gama.uni <- t(rowSums(array(t(GS),dim=c(p,p,n)),dims=2))/n
           for (j in 1:g){
             Gama[[j]] <- Gama.uni
             Sigma[[j]] <- Gama[[j]] + Delta[[j]]%*%t(Delta[[j]])
             shape[[j]] <- rep(0,p)
           }
        }

        pii[g] <- 1 - (sum(pii) - pii[g])

        zero.pos <- NULL
        zero.pos <- which(pii == 0)
        if(length(zero.pos) != 0){
          pii[zero.pos] <- 1e-10
          pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
        }

#        param <- teta
#        teta <- c()
#        for (k in 1:g) teta <- c(teta, mu[[k]], Gama[[k]][upper.tri(Gama[[k]], diag = TRUE)])
#        teta <- c(teta, pii) #teta para nu conhecido
        lk <- sum(log( d.mixedmvSN(y, pii, mu, Sigma, shape) ))
        criterio <- abs((lk/lkante)-1)

        lkante <- lk
        mu.old <- mu
        Delta.old <- Delta
        Gama.old <- Gama
        #print(count)
        #print(mu[[1]])
      }

    if (criteria == TRUE){
       cl <- apply(tal, 1, which.max)
       icl <- 0
       for (j in 1:g) icl <- icl+sum(log(pii[j]*dmvSN(y[cl==j,], mu[[j]], Sigma[[j]], shape[[j]])))
       #icl <- -2*icl+(3*g-1)*log(n)

    }

  }

     if(criteria == TRUE){
        if(uni.Gama){
           if((family == "t") | (family == "Normal")) d <- g*p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) + (g-1) #mu + Sigma + pi       
           else d <- g*2*p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) + (g-1) #mu + shape + Sigma + pi       
        }
        else{
          if((family == "t") | (family == "Normal")) d <- g*(p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) ) + (g-1) #mu + shape + Sigma + pi
          else d <- g*(2*p + length(Sigma[[1]][upper.tri(Sigma[[1]], diag = TRUE)]) ) + (g-1) #mu + shape + Sigma + pi
        }
        aic <- -2*lk + 2*d
        bic <- -2*lk + log(n)*d
        edc <- -2*lk + 0.2*sqrt(n)*d
        icl <- -2*icl + log(n)*d
        obj.out <- list(mu = mu, Sigma = Sigma, shape = shape, pii = pii, nu = nu, logLik = lk, aic = aic, bic = bic, edc = edc, icl=icl, iter = count, n = nrow(y), group = cl)
     }
     if(criteria == FALSE) obj.out <- list(mu = mu, Sigma = Sigma, shape = shape, pii = pii, nu = nu, iter = count, n = nrow(y), group =apply(tal, 1, which.max))

     if (group == FALSE) obj.out <- obj.out[-length(obj.out)]
     obj.out$uni.Gama <- uni.Gama    
     if (obs.prob == TRUE){
     nam <- c()
     for (i in 1:ncol(tal)) nam <- c(nam,paste("Group ",i,sep=""))
     dimnames(tal)[[2]] <- nam
     obj.out$obs.prob <- tal
     if((ncol(tal) - 1) > 1) obj.out$obs.prob[,ncol(tal)] <- 1 - rowSums(obj.out$obs.prob[,1:(ncol(tal)-1)])
     else obj.out$obs.prob[,ncol(tal)] <- 1 - obj.out$obs.prob[,1]
     obj.out$obs.prob[which(obj.out$obs.prob[,ncol(tal)] <0),ncol(tal)] <- 0.0000000000
     obj.out$obs.prob <- round(obj.out$obs.prob,10)
 }

     class(obj.out) <- family
     for(i in 1:length(obj.out$Sigma)) obj.out$Sigma[[i]] <- matrix.sqrt(obj.out$Sigma[[i]])    

     if(calc.im){
        IM <-  imm.smsn(as.matrix(y), obj.out)
        sdev <- sqrt(diag(solve(IM$IM)))
        obj.out$imm.sdev = sdev
     }
     class(obj.out) <- family
 
 obj.out
} #aqui  termina
