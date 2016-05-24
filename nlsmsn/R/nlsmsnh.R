################################################################################
##### Funçao Não Linear para X e Betas......Calculo da 1ra e 2da derivada  #####
################################################################################

smsn.nl.intern <- function(y, x, betas, sigma2, shape, nu, nlf, criteria, family, error, iter.max){
  # y: é o vetor de dados (amostra) de tamanho n
  # mu, sigma2, shape, pii: são os valores iniciais para o algorítmo EM. Cada um deles deve ser um vetor de tamanho g
  #                       (o algorítmo entende o número de componentes a ajustar baseado no tamanho desses vetores)
  #nu: valor inicial para o nu (no caso da Skew.cn deve ser um vetor bidimensional com valores entre 0 e 1)
  #criteria: TRUE ou FALSE - Caso queira que seja calculada a log-verossimilhança do modelo ajustado
  #cluster: TRUE ou FALSE - Caso queira um vetor de tamanho n indicando a qual componente a í-esima observação pertence
  #family: c("Skew.t","Skew.cn","Skew.slash","Skew.normal","Normal") - Skew.t: Ajusta por misturas de Skew-T
  #                                       - Skew.cn: Ajusta por misturas de Skew Normal Contaminada
  #                                       - Skew.slash: Ajusta por misturas de Skew-Slash (ATENÇÃO: AINDA COM PROBLEMAS)
  #                                       - Skew.normal: Ajusta por misturas de Skew Normal
  #                                       - Normal: Misturas de Normais
  #error: define o critério de parada do algorítmo.
  
  
  if (family == "t"){
  shape <- 0
      n <- length(y)
      p <- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
      k2<-(nu/2)*gamma((nu-2)/2)/gamma(nu/2)
        delta <- shape / (sqrt(1 + shape^2))
        Delta <- sqrt(sigma2)*delta
        Gama <- sigma2 - Delta^2
#      teta <- c(betas, Delta, Gama, nu)
      teta <- c(betas, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0
      while(criterio > error & count<=iter.max){
      count <- count + 1
     # print(count)
#        tal <- matrix(0, n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

              b<- -sqrt(2/pi)*k1
          mu.1<-nlf(x,betas)
          mu<-mu.1+b*Delta
          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu)/sqrt(sigma2))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          E=(2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2)*dt.ls(y, mu, sigma2,shape ,nu))
          u= ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2)*dt.ls(y, mu, sigma2,shape ,nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)
          S1<- u
          S2 <- (mutij_b*u + Mtij*E)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
#
#         betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),algorithm = "port",trace=TRUE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),trace=FALSE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
 #         Delta <- sum(S2*(y - mu.1)) / sum(S3)
          Delta <-0 
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
 #       }

        nu <-nu   #nu fixo

        param <- teta
#        teta <- c(betas, Delta, Gama,nu)
        teta <- c(betas, Delta, Gama)

        criterio <- sqrt((teta-param)%*%(teta-param))

#    if( (is.na(shape)|(abs(shape)>19))){
#          criterio=error/10
#          shape=20
#          betas<-betas.old
#          Delta<-Delta.old
#          Gama<-Gama.old
#           }
        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
      }
      vari<-sigma2*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos<-(y-mu.1)/sqrt(vari)
      if (criteria == TRUE){
        lk <- sum(log( dt.ls(y, mu, sigma2,shape ,nu) ))
        AIC <- (-2)*lk + 2*(p+2)     # Param = betas + sigma2 + nu 
        BIC <- (-2)*lk + (p+2)*log(n)
        EDC <- (-2)*lk + (p+2)*0.2*sqrt(n)

#        obj.out <- list(betas = betas, sigma2 = sigma2, shape = shape, nu = nu, loglik = lk, iter = count, n = length(y))
          obj.out <- list(u = u, betas = betas, sigma2 = sigma2, shape = shape, nu = nu, loglik = lk, AIC=AIC, BIC=BIC, EDC=EDC, iter = count,res=residuos, n = length(y))
      } else {
#        obj.out <- list(betas = betas, sigma2 = sigma2, shape = shape, nu = nu, iter = count, n = length(y))
         obj.out <- list(u = u, betas = betas, sigma2 = sigma2, shape = shape, nu = nu, iter = count,res=residuos, n = length(y))

      }
#    if (cluster == TRUE) obj.out[[length(obj.out)+1]] <- apply(tal, 1, which.max)
      class(obj.out) <- family
      return(obj.out)
  }

  if (family == "Skew.t"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
      k2<-(nu/2)*gamma((nu-2)/2)/gamma(nu/2)
        delta <- shape / (sqrt(1 + shape^2))
        Delta <- sqrt(sigma2)*delta
        Gama <- sigma2 - Delta^2

      teta <- c(betas, Delta, Gama, nu)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error & count<=iter.max){
      count <- count + 1
     # print(count)
#        tal <- matrix(0, n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

              b<- -sqrt(2/pi)*k1
          mu.1<-nlf(x,betas)
          mu<-mu.1+b*Delta
          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu)/sqrt(sigma2))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          E=(2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2)*dt.ls(y, mu, sigma2,shape ,nu))
          u= ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2)*dt.ls(y, mu, sigma2,shape ,nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)


          S1<- u
          S2 <- (mutij_b*u + Mtij*E)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
#
#         betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),algorithm = "port",trace=TRUE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),trace=FALSE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <- sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))

 #       }

        nu <-nu   #nu fixo

        param <- teta
        teta <- c(betas, Delta, Gama,nu)
        criterio <- sqrt((teta-param)%*%(teta-param))

#    if( (is.na(shape)|(abs(shape)>16))){
#          criterio=error/10
#          shape=20
#          betas<-betas.old
#          Delta<-Delta.old
#          Gama<-Gama.old
#           }
        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
      }
      vari<-sigma2*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos<-(y-mu.1)/sqrt(vari)
      residuos=1
      if (criteria == TRUE){
        lk <- sum(log( dt.ls(y, mu, sigma2,shape ,nu) ))
        AIC <- (-2)*lk + 2*(p+3)              # Param = betas + sigma2 + shape + nu 
        BIC <- (-2)*lk + (p+3)*log(n)
        EDC <- (-2)*lk + (p+3)*0.2*sqrt(n)
#        obj.out <- list(betas = betas, sigma2 = sigma2, shape = shape, nu = nu, loglik = lk, iter = count, n = length(y))
          obj.out <- list(u = u, betas = betas, sigma2 = sigma2, shape = shape, nu = nu, loglik = lk, AIC=AIC, BIC=BIC, EDC=EDC, iter = count,res=residuos, n = length(y))
      } else {
#        obj.out <- list(betas = betas, sigma2 = sigma2, shape = shape, nu = nu, iter = count, n = length(y))
         obj.out <- list(u = u, betas = betas, sigma2 = sigma2, shape = shape, nu = nu, iter = count,res=residuos, n = length(y))

      }
#    if (cluster == TRUE) obj.out[[length(obj.out)+1]] <- apply(tal, 1, which.max)
      class(obj.out) <- family
      return(obj.out)
  }

  if (family == "Skew.cn"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-nu[1]/nu[2]^(1/2)+1-nu[1]
      k2<-nu[1]/nu[2]+1+nu[1]
        delta <- shape / (sqrt(1 + shape^2))
        Delta <- sqrt(sigma2)*delta
        Gama <- sigma2 - Delta^2

      teta <- c(betas, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error & count<=iter.max){
      count <- count + 1
      #print(count)

        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-nlf(x,betas)
          mu<-mu.1+b*Delta
          ### E-step: calculando ui, tui, tui2 ###

          dj <- ((y - mu)/sqrt(sigma2))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          u=(2/dSNC(y,mu,sigma2,shape,nu))*(nu[1]*nu[2]*dnorm(y,mu,sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu,sqrt(sigma2))*pnorm(A,0,1))
#         E=(2/dSNC(y,mu,sigma2,shape,nu))*(nu[1]*sqrt(nu[2])*dnorm(y,mu,sqrt(sigma2/nu[2]))*dnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu,sqrt(sigma2))*dnorm(A,0,1))
          E=(2/dSNC(y,mu,sigma2,shape,nu))*(nu[1]*sqrt(nu[2])*dnorm(y,mu,sqrt(sigma2/nu[2]))*dnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu,sqrt(sigma2))*dnorm(A,0,1))
          S1<- u
          S2 <- (mutij_b*u + Mtij*E)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
#
#          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),algorithm = "port",trace=TRUE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
         betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),trace=FALSE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <- sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))

 #       }

        nu <-nu   #nu fixo

        param <- teta
        teta <- c(betas, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

#        if( (is.na(shape)|(abs(shape)>19))){
#          criterio=error/10
#          shape=20
#          betas<-betas.old
#          Delta<-Delta.old
#          Gama<-Gama.old
#           }

        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
      }
      vari<-sigma2*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos<-(y-mu.1)/sqrt(vari)

      if (criteria == TRUE){
        lk <- sum(log(dSNC (y, mu, sigma2, shape, nu)))
       AIC <- (-2)*lk + 2*(p+4)                       # Param = betas + sigma2 + shape + nu1 + nu2 
       BIC <- (-2)*lk + (p+4)*log(n)
       EDC <- (-2)*lk + (p+4)*0.2*sqrt(n)
#       obj.out <- list(betas =betas, sigma2 = sigma2, shape = shape, nu = nu, loglik = lk, iter = count, n = length(y))
        obj.out <- list(u = u, betas = betas, sigma2 = sigma2, shape = shape, nu = nu, loglik = lk, AIC=AIC, BIC=BIC,EDC=EDC, iter = count,res=residuos ,n = length(y))
      } else {
#       obj.out <- list(betas = betas, sigma2 = sigma2, shape = shape, nu = nu, iter = count, n = length(y))
        obj.out <- list(u = u, betas = betas, sigma2 = sigma2, shape = shape, nu = nu, iter = count, res=residuos,n = length(y))
      }
#      if (cluster == TRUE) obj.out[[length(obj.out)+1]] <- apply(tal, 1, which.max)
      class(obj.out) <- family
      return(obj.out)
  }

  if (family == "Skew.slash"){
#      print("modelo Slash....")
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-2*nu/(2*nu-1)
      k2<-2*nu/(2*nu-2)
        delta <- shape / (sqrt(1 + shape^2))
        Delta <- sqrt(sigma2)*delta
        Gama <- sigma2 - Delta^2

      teta <- c(betas, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

     while(criterio > error & count<=iter.max){
      count <- count + 1
      #print(count)
          u <- vector(mode="numeric",length=n)
          E <- vector(mode="numeric",length=n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-nlf(x,betas)
          mu<-mu.1+b*Delta
          ### E-step: calculando ui, tui, tui2 ###

          dj <- ((y - mu)/sqrt(sigma2))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          for(i in 1:n){
            E[i] <- (((2^(nu + 1))*nu*gamma(nu + 1))/(dSS(y[i], mu[i], sigma2, shape, nu)*pi*sqrt(sigma2)))* ((dj[i]+A[i]^2)^(-nu-1))*pgamma(1,nu+1,(dj[i]+A[i]^2)/2)
            faux <- function(u) u^(nu+0.5)*exp(-u*dj[i]/2)*pnorm(u^(1/2)*A[i])
            aux22 <- integrate(faux,0,1)$value
            u[i] <- ((sqrt(2)*nu) / (dSS(y[i], mu[i], sigma2, shape, nu)*sqrt(pi)*sqrt(sigma2)))*aux22
          }

          S1<- u
          S2 <- (mutij_b*u + Mtij*E)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
         betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),algorithm = "port",trace=TRUE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
#          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),trace=TRUE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <- sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))

 #       }

          nu <-nu   #nu fixo
       param <- teta
        teta <- c(betas, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

#          if( (is.na(shape)|(abs(shape)>16))){
#          criterio=error/10
#          shape=20
#          betas<-betas.old
#          Delta<-Delta.old
#          Gama<-Gama.old
#           }

        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
      }
      vari<-sigma2*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos <-(y-mu.1)/sqrt(vari)
      if (criteria == TRUE){
        lk <- sum(log(dSS(y, mu, sigma2, shape, nu)))
       AIC <- -2*lk + 2*(p+3)                          # Param = betas + sigma2 + shape + nu 
       BIC <- (-2)*lk + (p+3)*log(n)
       EDC <- (-2)*lk + (p+3)*0.2*sqrt(n)
#        obj.out <- list(betas = betas, sigma2 = sigma2, shape = shape, nu = nu,loglik = lk, iter = count, n = length(y))
         obj.out <- list(u = u, betas = betas, sigma2 = sigma2, shape = shape, nu = nu, loglik = lk,  AIC=AIC, BIC=BIC, EDC=EDC, iter = count,res=residuos,n = length(y))
      } else {
#        obj.out <- list(betas = betas, sigma2 = sigma2, shape = shape, nu = nu,iter = count, n = length(y))
        obj.out <- list(u = u, betas = betas, sigma2 = sigma2, shape = shape, nu = nu, iter = count,res=residuos, n = length(y))
      }
#      if (cluster == TRUE) obj.out[[length(obj.out)+1]] <- apply(tal, 1, which.max)
      class(obj.out) <- family
      return(obj.out)
  }

################################################################################

  if (family == "Skew.normal"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-1
      k2=1
        delta <- shape / (sqrt(1 + shape^2))
        Delta <- sqrt(sigma2)*delta
        Gama <- sigma2 - Delta^2

      teta <- c(betas, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error & count<=iter.max){
      count <- count + 1
      #print(count)
#        tal <- matrix(0, n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-nlf(x,betas.old)
          mu<-mu.1+b*Delta.old
          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu)/sqrt(sigma2))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          E = dnorm(mutij/Mtij) / pnorm(mutij/Mtij)
          u = rep(1, n)


          S1<- u
          S2 <- (mutij_b*u + Mtij*E)
          S3 <- ((mutij_b)^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
#
          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),trace=FALSE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <- sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))


 #       }


        param <- teta
        teta <- c(betas, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))
#          if( (is.na(shape)|(abs(shape)>16))){
#          criterio=error/10
#          shape=20
#          betas<-betas.old
#          Delta<-Delta.old
#         Gama<-Gama.old
#           }

        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
      }
      vari<-sigma2*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos <-(y-mu.1)/sqrt(vari)
      if (criteria == TRUE){
        lk <- sum(log( dSN (y, mu, sigma2, shape) ))
       AIC <- -2*lk + 2*(p+2)                          # Param = betas + sigma2 + shape  
       BIC <- (-2)*lk + (p+2)*log(n)
       EDC <- (-2)*lk + (p+2)*0.2*sqrt(n)
#        obj.out <- list(betas = betas, sigma2 = sigma2, shape = shape, loglik = lk, iter = count, n = length(y))
         obj.out <- list(u = u, betas = betas, sigma2 = sigma2, shape = shape, nu = nu, loglik = lk,  AIC=AIC, BIC=BIC, EDC=EDC, iter = count, res=residuos, n = length(y))
      } else {
#        obj.out <- list(betas = betas, sigma2 = sigma2, shape = shape, iter = count, n = length(y))
         obj.out <- list(u = u, betas = betas, sigma2 = sigma2, shape = shape, nu = nu, iter = count,res=residuos ,n = length(y))
      }
#    if (cluster == TRUE) obj.out[[length(obj.out)+1]] <- apply(tal, 1, which.max)
      class(obj.out) <- family
      return(obj.out)
  }

  if (family == "Normal"){
      shape <- 0
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-1
      k2<-1
        delta <- shape / (sqrt(1 + shape^2))
        Delta <- sqrt(sigma2)*delta
        Gama <- sigma2 - Delta^2

      teta <- c(betas, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error & count<=iter.max){
      count <- count + 1
#      print(count)
#        tal <- matrix(0, n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-nlf(x,betas.old)
          mu<-mu.1+b*Delta.old
          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu)/sqrt(sigma2))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          E = dnorm(mutij/Mtij) / pnorm(mutij/Mtij)
          u = rep(1, n)


          S1<- u
          S2 <- (mutij_b*u + Mtij*E)
          S3 <- ((mutij_b)^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
#


#          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),algorithm = "port",trace=TRUE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old), trace=FALSE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <- 0 #sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))

 #       }


        param <- teta
        teta <- c(betas, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))
#        if( (is.na(shape)|(shape>6))){
#          criterio=error/10
#          shape=20
#          betas<-betas.old
#          Delta<-Delta.old
#          Gama<-Gama.old
#           }
        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
      }
      vari<-sigma2*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos<-(y-mu.1)/sqrt(vari)
      if (criteria == TRUE){
        lk <- sum(log( dSN (y, mu, sigma2, shape) ))
       AIC <- (-2)*lk + 2*(p+1)                                  # Param = betas + sigma2 
       BIC <- (-2)*lk + (p+1)*log(n)
       EDC <- (-2)*lk + (p+1)*0.2*sqrt(n)
        obj.out <- list(betas = betas, sigma2 = sigma2, shape = shape, loglik = lk,  AIC=AIC,  BIC=BIC, EDC=EDC, iter = count,res=residuos, n = length(y))
      } else {
        obj.out <- list(betas = betas, sigma2 = sigma2, shape = shape, iter = count,res=residuos, n = length(y))
      }
      #if (cluster == TRUE) obj.out[[length(obj.out)+1]] <- apply(tal, 1, which.max)
      class(obj.out) <- family
      return(obj.out)
  }
}



####################################################################
##########       Algorítmo EM para SNI            ###############
#                   alterado em 05/02/2010                         ##

smsn.nlh.intern <- function(y, x, z, betas, sigma2, shape, rho, nu, nlf, rho.func, criteria, family , error, iter.max){
  # y: é o vetor de dados (amostra) de tamanho n
  # rho: parametro do Heter
  # z matrix planeja do heter.
  # mu, sigma2, shape, pii: são os valores iniciais para o algorítmo EM. Cada um deles deve ser um vetor de tamanho g
  #                       (o algorítmo entende o número de componentes a ajustar baseado no tamanho desses vetores)
  #nu: valor inicial para o nu (no caso da Skew.cn deve ser um vetor bidimensional com valores entre 0 e 1)
  #loglik: TRUE ou FALSE - Caso queira que seja calculada a log-verossimilhança do modelo ajustado
  #cluster: TRUE ou FALSE - Caso queira um vetor de tamanho n indicando a qual componente a í-esima observação pertence
  #family: c("Skew.t","Skew.cn","Skew.slash","Skew.normal","Normal") - Skew.t: Ajusta por misturas de Skew-T
  #                                       - Skew.slashL: Ajusta por misturas de Skew-Slash (ATENÇÃO: AINDA COM PROBLEMAS)
  #                                       - Skew.normal: Ajusta por misturas de Skew Normal
  #                                       - Normal: Misturas de Normais
  #error: define o critério de parada do algorítmo.


  if (family == "t"){
      shape <- 0
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
      k2<-(nu/2)*gamma((nu-2)/2)/gamma(nu/2)
      delta <- shape / (sqrt(1 + shape^2))
      Delta <- sqrt(sigma2)*delta
      Gama <- sigma2 - Delta^2
      wi<-nlwi(z,rho,rho.func)
      teta <- c(betas,rho, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old <- rho
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error & count<=iter.max){
      count <- count + 1
      #print(count)
#        tal <- matrix(0, n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-nlf(x,betas)
          mu<-mu.1+b*Delta*sqrt(wi)
          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu)/sqrt(sigma2*wi))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b<-mutij+b
          A <- mutij / (Mtij)

          E=(2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2*wi)*dt.ls(y, mu, sigma2*wi,shape ,nu))
          u= ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2*wi)*dt.ls(y, mu, sigma2*wi,shape ,nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)

          S1<- u/wi
          S2 <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
#

          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),trace=FALSE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <-0
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.ST <- function(rho) sum(log(dtH.ls(rho, y, z, mu , sigma2 ,shape, nu, rho.func)))
          rho <- optimize(logvero.ST, c(-5,-0.0001), tol = 0.0000001, maximum = TRUE)$maximum
 #       }

        nu <-nu   #nu fixo
 
       param <- teta
        teta <- c(betas,rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        rho.old <- rho
        Delta.old <- Delta
        Gama.old <- Gama
        wi<-nlwi(z,rho,rho.func)
      #  if (count==1000) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos<-(y-mu.1)/sqrt(vari)
      if (criteria == TRUE){
        lk <- sum(log(dtH.ls(rho, y, z, mu , sigma2 ,shape, nu, rho.func)))
       AIC <- (-2)*lk + 2*(p+3)                     # Param = betas + sigma2 + rho + nu  
       BIC <- (-2)*lk + (p+3)*log(n)
       EDC <- (-2)*lk + (p+3)*0.2*sqrt(n)
        obj.out <- list(betas = betas, rho=rho,sigma2 = sigma2, shape = shape, nu = nu, loglik = lk,  AIC=AIC, BIC=BIC, EDC=EDC, iter = count, res=residuos,n = length(y),mah=dj)
      } else {
        obj.out <- list(betas = betas,rho=rho, sigma2 = sigma2, shape = shape, nu = nu, iter = count, res=residuos,n = length(y),mah=dj)
      }
      class(obj.out) <- family
      return(obj.out)
  }

  if (family == "Skew.t"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2)
      k2<-(nu/2)*gamma((nu-2)/2)/gamma(nu/2)
      delta <- shape / (sqrt(1 + shape^2))
      Delta <- sqrt(sigma2)*delta
      Gama <- sigma2 - Delta^2
      wi<-nlwi(z,rho,rho.func)
      teta <- c(betas,rho, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old <- rho
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error & count<=iter.max){
      count <- count + 1
      #print(count)
#        tal <- matrix(0, n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-nlf(x,betas)
          mu<-mu.1+b*Delta*sqrt(wi)
          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu)/sqrt(sigma2*wi))^2           
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b<-mutij+b
          A <- mutij / (Mtij)

          E=(2*(nu)^(nu/2)*gamma((2+nu)/2)*((dj + nu + A^2))^(-(2+nu)/2)) / (gamma(nu/2)*pi*sqrt(sigma2*wi)*dt.ls(y, mu, sigma2*wi,shape ,nu))
          u= ((4*(nu)^(nu/2)*gamma((3+nu)/2)*(dj + nu)^(-(nu+3)/2)) / (gamma(nu/2)*sqrt(pi)*sqrt(sigma2*wi)*dt.ls(y, mu, sigma2*wi,shape ,nu)) )*pt(sqrt((3+nu)/(dj+nu))*A,3+nu)


          S1<- u/wi
          S2 <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
#

          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),trace=FALSE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <-sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.ST <- function(rho) sum(log(dtH.ls(rho, y, z, mu , sigma2 ,shape, nu, rho.func)))
          rho <- optimize(logvero.ST, c(-5,-0.0001), tol = 0.0000001, maximum = TRUE)$maximum
 #       }

        nu <-nu   #nu fixo

        param <- teta
        teta <- c(betas,rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        rho.old <- rho
        Delta.old <- Delta
        Gama.old <- Gama
        wi<-nlwi(z,rho,rho.func)
#        if (count==1000) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos<-(y-mu.1)/sqrt(vari)
      if (criteria == TRUE){
        lk <- sum(log(dtH.ls(rho, y, z, mu , sigma2 ,shape, nu, rho.func)))
       AIC <- (-2)*lk + 2*(p+4)                            # Param = betas + sigma2 + shape + rho + nu  
       BIC <- (-2)*lk + (p+4)*log(n)       
        EDC <- (-2)*lk + (p+4)*0.2*sqrt(n)
        obj.out <- list(betas = betas, rho=rho,sigma2 = sigma2, shape = shape, nu = nu, loglik = lk,  AIC=AIC, BIC=BIC, EDC=EDC, iter = count, res=residuos,n = length(y),mah=dj)
      } else {
        obj.out <- list(betas = betas,rho=rho, sigma2 = sigma2, shape = shape, nu = nu, iter = count, res=residuos,n = length(y),mah=dj)
      }
      class(obj.out) <- family
      return(obj.out)
  }

  if (family == "Skew.cn"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-nu[1]/nu[2]^(1/2)+1-nu[1]
      k2<-nu[1]/nu[2]+1+nu[1]


        delta <- shape / (sqrt(1 + shape^2))
        Delta <- sqrt(sigma2)*delta
        Gama <- sigma2 - Delta^2
      wi<-nlwi(z,rho,rho.func)
      teta <- c(betas,rho, Delta, Gama)


      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old<-rho
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0


      while(criterio > error & count<=iter.max )
      {
      count <- count + 1
      #print(count)

        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-nlf(x,betas)
          mu<-mu.1+b*Delta*sqrt(wi)
          ### E-step: calculando ui, tui, tui2 ###

          dj <- ((y - mu)/sqrt(wi*sigma2))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          u=(2/dSNC(y,mu,sigma2*wi,shape,nu))*(nu[1]*nu[2]*dnorm(y,mu,sqrt(wi*sigma2/nu[2]))*pnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu,sqrt(sigma2*wi))*pnorm(A,0,1))
          E=(2/dSNC(y,mu,sigma2*wi,shape,nu))*(nu[1]*sqrt(nu[2])*dnorm(y,mu,sqrt(sigma2*wi/nu[2]))*dnorm(sqrt(nu[2])*A,0,1)+(1-nu[1])*dnorm(y,mu,sqrt(sigma2*wi))*dnorm(A,0,1))

          S1<- u/wi
          S2 <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)



          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
#
          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),trace=FALSE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <- sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.SCN <- function(rho) sum(log(dSNCH(rho, y, z, mu , sigma2 ,shape, nu, rho.func)))
          rho <- optimize(logvero.SCN, c(-5,-0.00001), tol = 0.0000001, maximum = TRUE)$maximum
 #       }

        nu <-nu   #nu fixo

        param <- teta
        teta <- c(betas, rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
        rho.old <- rho
        wi<-nlwi(z,rho,rho.func)
#        if (count==1000) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos<-(y-mu.1)/sqrt(vari)
      if (criteria == TRUE){
        lk <- sum(log(dSNCH(rho, y, z, mu , sigma2 ,shape, nu, rho.func)))
       AIC <- (-2)*lk + 2*(p+5)                            # Param = betas + sigma2 + shape + rho + nu1 + nu2  
       BIC <- (-2)*lk + (p+5)*log(n)
       EDC <- (-2)*lk + (p+5)*0.2*sqrt(n)
        obj.out <- list(betas =betas, rho=rho, sigma2 = sigma2, shape = shape, nu = nu, loglik = lk,  AIC=AIC, BIC=BIC, EDC=EDC, iter = count, res=residuos,n = length(y),mah=dj)
      } else {
        obj.out <- list(betas = betas, rho=rho, sigma2 = sigma2, shape = shape, nu = nu, iter = count, res=residuos,n = length(y),mah=dj)
      }
      class(obj.out) <- family
      return(obj.out)
  }

  if (family == "Skew.slash"){
      print("modelo Slash....")
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-2*nu/(2*nu-1)
      k2<-2*nu/(2*nu-2)
        delta <- shape / (sqrt(1 + shape^2))
        Delta <- sqrt(sigma2)*delta
        Gama <- sigma2 - Delta^2

      teta <- c(betas, rho, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old<-rho
      wi<-nlwi(z,rho,rho.func)
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error & count<=iter.max){
      count <- count + 1
      #print(count)
          u <- vector(mode="numeric",length=n)
          E <- vector(mode="numeric",length=n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-nlf(x,betas)
          mu<-mu.1+b*Delta*sqrt(wi)
          ### E-step: calculando ui, tui, tui2 ###

          dj <- ((y - mu)/sqrt(wi*sigma2))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          for(i in 1:n){
            E[i] <- (((2^(nu + 1))*nu*gamma(nu + 1))/(dSS(y[i], mu[i], sigma2*wi[i], shape, nu)*pi*sqrt(sigma2*wi[i])))* ((dj[i]+A[i]^2)^(-nu-1))*pgamma(1,nu+1,(dj[i]+A[i]^2)/2)
            faux <- function(u) u^(nu+0.5)*exp(-u*dj[i]/2)*pnorm(u^(1/2)*A[i])
            aux22 <- integrate(faux,0,1)$value
            u[i] <- ((sqrt(2)*nu) / (dSS(y[i], mu[i], sigma2*wi[i], shape, nu)*sqrt(pi)*sqrt(sigma2*wi[i])))*aux22
          }

          S1<- u/wi
          S2 <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S3 <- (mutij_b^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)


          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
#


          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),trace=FALSE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <- sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.SS <- function(rho) sum(log(dSSH(rho, y, z, mu , sigma2 ,shape, nu, rho.func)))
          rho <- optimize(logvero.SS, c(-5,-0.0001), tol = 0.0000001, maximum = TRUE)$maximum
 #       }

        nu <-nu   #nu fixo

        param <- teta
        teta <- c(betas, rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
        rho.old <- rho
        wi<-nlwi(z,rho,rho.func)
  #      if (count==200) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos <-(y-mu.1)/sqrt(vari)
      if (criteria == TRUE){
        lk <- sum(log(dSSH(rho, y, z, mu , sigma2 ,shape, nu, rho.func)))
       AIC <- (-2)*lk + 2*(p+4)                                     # Param = betas + sigma2 + shape + rho + nu  
       BIC <- (-2)*lk + (p+4)*log(n)
       EDC <- (-2)*lk + (p+4)*0.2*sqrt(n)
        obj.out <- list(betas = betas, rho=rho, sigma2 = sigma2, shape = shape, nu = nu,loglik = lk,  AIC=AIC, BIC=BIC, EDC=EDC, iter = count, res=residuos,n = length(y),mah=dj)
      } else {
        obj.out <- list(betas = betas, rho=rho, sigma2 = sigma2, shape = shape, nu = nu,iter = count, res=residuos, n = length(y),mah=dj)
      }

      class(obj.out) <- family
      return(obj.out)
  }

  if (family == "Skew.normal"){
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-1
      k2=1
        delta <- shape / (sqrt(1 + shape^2))
        Delta <- sqrt(sigma2)*delta
        Gama <- sigma2 - Delta^2
        wi<-nlwi(z,rho,rho.func)
      teta <- c(betas, rho, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old<-rho
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error & count<=iter.max){
      count <- count + 1
      #print(count)
#        tal <- matrix(0, n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-nlf(x,betas.old)
          mu<-mu.1++b*Delta*sqrt(wi)
          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu)/sqrt(sigma2*wi))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          E = dnorm(mutij/Mtij) / pnorm(mutij/Mtij)
          u = rep(1, n)

           S1<- u/wi
          S2 <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S3 <- ((mutij_b)^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
#


          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),trace=FALSE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <- sum(S2*(y - mu.1)) / sum(S3)
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.SN <- function(rho) sum(log(dSNH(rho, y, z, mu , sigma2 ,shape, rho.func)))
          rho <- optimize(logvero.SN, c(-5,-0.0001), tol = 0.0000001, maximum = TRUE)$maximum
 #       }

        param <- teta
        teta <- c(betas, rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
        rho.old<-rho
        wi<-nlwi(z,rho,rho.func)
 #   if (count==1000) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos <-(y-mu.1)/sqrt(vari)
      if (criteria == TRUE){
        lk <- sum(log( dSNH (rho, y, z, mu , sigma2 ,shape, rho.func) ))
       AIC <- (-2)*lk+ 2*(p+3)                                    # Param = betas + sigma2 + shape + rho   
       BIC <- (-2)*lk + (p+3)*log(n)
       EDC <- (-2)*lk + (p+3)*0.2*sqrt(n)
        obj.out <- list(betas = betas, rho=rho, sigma2 = sigma2, shape = shape, loglik = lk, AIC=AIC, BIC=BIC, EDC=EDC, iter = count, res=residuos,n = length(y),mah=dj)
      } else {
        obj.out <- list(betas = betas, rho=rho, sigma2 = sigma2, shape = shape, iter = count,res=residuos,n = length(y),mah=dj)
      }

      class(obj.out) <- family
      return(obj.out)
  }

 if (family == "Normal"){
      shape <- 0
      n <- length(y)
      p<- length(betas)
      delta <- Delta <- Gama <- c(0)
      k1<-1
      k2=1
      delta <- shape / (sqrt(1 + shape^2))
      Delta <- sqrt(sigma2)*delta
      Gama <- sigma2 - Delta^2
      wi<-nlwi(z,rho,rho.func)
      teta <- c(betas, rho, Delta, Gama)
      betas.old <- betas
      Delta.old <- Delta
      Gama.old <- Gama
      rho.old<-rho
      b<- -sqrt(2/pi)*k1
      criterio <- 1
      count <- 0

      while(criterio > error & count<=iter.max){
      count <- count + 1
      #print(count)
#        tal <- matrix(0, n)
        S1 <- matrix(0, n)
        S2 <- matrix(0, n)
        S3 <- matrix(0, n)

          mu.1<-nlf(x,betas.old)
          mu<-mu.1++b*Delta*sqrt(wi)
          ### E-step: calculando ui, tui, tui2 ###
          dj <- ((y - mu)/sqrt(sigma2*wi))^2
          Mtij2 <- 1/(1 + (Delta^2)*(Gama^(-1)))
          Mtij <- sqrt(Mtij2)
          mutij <- Mtij2*Delta*(Gama^(-1))*(y - mu)/sqrt(wi)
          mutij_b<-mutij+b
          A <- mutij / Mtij

          E = dnorm(mutij/Mtij) / pnorm(mutij/Mtij)
          u = rep(1, n)


          S1<- u/wi
          S2 <- (mutij_b*u + Mtij*E)/sqrt(wi)
          S3 <- ((mutij_b)^2*u + Mtij2 + Mtij*(mutij_b+b)*E)
          ymod<-(y-Delta.old*S2/S1)

          ### M-step: atualizar mu, Delta, Gama, sigma2 ###
#


          betas <-as.vector(coef(nls(ymod ~ nlf(x,betas),start=list(betas=betas.old),trace=FALSE,weights=S1))) #sum(S1*y - Delta.old*S2) / sum(S1)
          #print(betas)
          Gama <- sum(S1*(y - mu.1)^2 - 2*(y - mu.1)*Delta.old*S2 + Delta.old^2*S3) / n
          Delta <- 0
          sigma2 <- Gama + Delta^2
          shape <- ((sigma2^(-1/2))*Delta )/(sqrt(1 - (Delta^2)*(sigma2^(-1))))
          logvero.SN <- function(rho) sum(log(dSNH(rho, y, z, mu , sigma2 ,shape, rho.func)))
          rho <- optimize(logvero.SN, c(-5,-0.00001), tol = 0.0000001, maximum = TRUE)$maximum
 #       }


        param <- teta
        teta <- c(betas, rho, Delta, Gama)
        criterio <- sqrt((teta-param)%*%(teta-param))

        betas.old <- betas
        Delta.old <- Delta
        Gama.old <- Gama
        rho.old<-rho
        wi<-nlwi(z,rho,rho.func)
        if (count==200) {criterio= error}
      }
      vari<-sigma2*wi*(k2-2/pi*k1^2*shape^2/(1+shape^2))
      residuos <-(y-mu.1)/sqrt(vari)
      if (criteria == TRUE){
        lk <- sum(log( dSNH (rho, y, z, mu , sigma2 ,shape, rho.func) ))
       AIC <- (-2)*lk + 2*(p+2)                               # Param = betas + sigma2 + rho 
       BIC <- (-2)*lk + (p+2)*log(n)
       EDC <- (-2)*lk + (p+2)*0.2*sqrt(n)
        obj.out <- list(betas = betas, rho=rho, sigma2 = sigma2, shape = shape, loglik = lk,  AIC=AIC, BIC=BIC, EDC=EDC, iter = count, res=residuos,n = length(y),mah=dj)
      } else {
        obj.out <- list(betas = betas, rho=rho, sigma2 = sigma2, shape = shape, iter = count,res=residuos,n = length(y),mah=dj)
      }

      class(obj.out) <- family
      return(obj.out)
  }
}












