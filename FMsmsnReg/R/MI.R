###############################################################################
#Teste
#Se rodar a seguinte linha de embaixo ou seja a "source" vai carregar diretamente o exemplo com do modelo parST
#source("~/Dropbox/Tesis Doutorado IME-USP Luis Benites/IdeaMisturasCensuras/paper 3/MisturasErro_SMSN/R code/Application/Application2.R")

###############################################################################
im.FMsmsnReg     <- function(y, x1, g=NULL, model=NULL,initial.model=TRUE,family = NULL,Abetas = NULL, medj= NULL, sigma2 = NULL, shape = NULL, pii = NULL, nu=NULL)
{
  #Definiendo os parametros
  n              <- length(y) #Tamanho de amostra
  y              <- as.numeric(y)
  if(initial.model==TRUE) #isto es usado cuando es dado o model como argumento mais no caso em que nao seja dado es necesario o argumento
  {                       #initial.model
    if((class(model) != "Skew.t") && (class(model) != "Skew.cn") && (class(model) != "Skew.slash") && (class(model) != "Skew.normal")) stop(paste("Family",class(model),"not recognized.\n",sep=" "))
    g                  <- length(model$pii)
    Abetas             <- model$Abetas
    medj               <- model$medj
    sigma2             <- model$sigma2
    shape              <- model$shape
    pii                <- model$pii
    if(class(model)!="Skew.normal")  nu  <- model$nu
  }else{
    model               <- 0
    class(model)        <- family
  }


  #########################################
  #Parametros estimador pelo Algoritmo EM
  beta0          <- Abetas[1]
  betas          <- as.matrix(Abetas[2:(length(Abetas))])   # parameters of regression dimension "p"
  x              <- as.matrix(x1[,2:(length(Abetas))])
  delta          <- shape/(sqrt(1 + shape^2))
  Delta          <- sqrt(sigma2)*delta
  ##########################################Definiendo os parametros
  n              <- length(y) #Tamanho de amostra

  if(class(model)=="Skew.t")
  {

  k1             <- sqrt(nu/2)*gamma((nu-1)/2)/gamma(nu/2) #Skew.t
  b              <- -sqrt(2/pi)*k1 #Skew.t

  mu=mu1         <- matrix(0,n,g)
  varnu          <- rep(0,g) # beta0 + mu_j #Parametrizacao
  for (k in 1:g)
  {
    varnu[k]     <- beta0+medj[k]
    mu1[,k]      <- x%*%betas
    mu[,k]       <- mu1[,k]+varnu[k] + b*Delta[k]
  }

  #Verosimilhanca da mixtura da Skew.T
  lk             <- d.mixedST(y, pii, mu, sigma2, shape, nu)

  #Funcao do Skew.T
  I.Phi          <- function(w=0, Ai=0, di=0, nu=0) as.numeric(((2^w*nu^(nu/2)*gamma(w + nu/2))/(gamma(nu/2)*(nu + di)^(nu/2 + w)))*pt((Ai/(di + nu)^(0.5))*sqrt(2*w + nu), 2*w + nu))
  I.phi          <- function(w=0, Ai=0, di=0, nu=0) as.numeric((2^w*nu^(nu/2)/(sqrt(2*pi)*gamma(nu/2)))*(1/(di + Ai^2 + nu))^((nu + 2*w)/2)*gamma((nu + 2*w)/2))

  ################################
  #Definiendo os objetos Si
  Sibetas        <- matrix(0,length(betas),n)
  Sivarnu        <- matrix(0,g,n)
  Sisigma2       <- matrix(0,g,n)
  Sishape        <- matrix(0,g,n)
  Sipi           <- matrix(0,g-1,n)
  ################################

  #Nesta parte estou obtenido Sibetas onde considero "beta1", "beta2", pois na derivada respecto de "varnu" ja esta incluso o "beta0"
  #pelo que a derivada respecto a "beta" seria sem considera "beta0"
  #Lembrar que para cada "i" teremos a suma das "g" componentes
  for(i in 1:n)
  {
    Dr           <- sqrt(sigma2)
    B            <- y - mu
    di           <- B[i,]^2/Dr^2
    Ai           <- shape*B[i,]/Dr

    #Para betas
    dmu1         <- B[i,]/Dr^3
    dmu2         <- shape/Dr^2
    for(j in 1:g) Sibetas[,i]    <- Sibetas[,i] + (2/sqrt(2*pi))*pii[j]*(dmu1[j]*I.Phi(3/2, Ai[j], di[j], nu) - dmu2[j]*I.phi(1, Ai[j], di[j], nu))*x[i,]/lk[i]
    for(j in 1:(g-1)) Sipi[j,i]  <- (dt.ls(y[i], mu[i,j], sigma2[j], shape[j], nu) - dt.ls(y[i], mu[i,g], sigma2[g], shape[g], nu))/lk[i]
  }

  #Nesta parte vamos obter os demais Si, onde tambem es considerado o Sibetas (de tres componentes, pois temos 3 variaveis que incluye o intercepto)
  soma           <- 0
  for(i in 1:n)
  {
    Dr           <- sqrt(sigma2)
    B            <- y - mu
    di           <- B[i,]^2/Dr^2
    Ai           <- shape*B[i,]/Dr
    ################################################################
    #Para varnu
    dmu1         <- B[i,]/Dr^3
    dmu2         <- shape/Dr^2

    #Para sigma
    dsigma1      <- 1/Dr^3
    dsigma2      <- B[i,]*(B[i,] + b*Delta)/Dr^5
    dsigma3      <- shape*(B[i,] + b*Delta)/Dr^4

    #Para os elementos de lambda
    dlambda1     <- (B[i,]*b)/(Dr^2*(1+shape^2)^(3/2))
    dlambda2     <- (1/Dr^2)*(B[i,] - b*Delta/(1+shape^2))
    ###############################################################

    Sivarnu[,i]  <- (2/sqrt(2*pi))*as.matrix((pii/lk[i])*( dmu1     * I.Phi(3/2, Ai, di, nu)     - dmu2 * I.phi(1, Ai, di, nu)))
    Sishape[,i]  <- (2/sqrt(2*pi))*as.matrix((pii/lk[i])*( dlambda1 * I.Phi(3/2, Ai, di, nu) + dlambda2 * I.phi(1 , Ai, di, nu)))
    Sisigma2[,i] <- (1/sqrt(2*pi))*as.matrix((pii/lk[i])*(-dsigma1  * I.Phi(1/2, Ai, di, nu) + dsigma2  * I.Phi(3/2, Ai, di, nu) -  dsigma3*I.phi(1, Ai, di, nu)))

    S            <- c(Sibetas[,i],Sipi[,i],Sivarnu[,i],Sisigma2[,i],Sishape[,i])
    soma         <- soma + S%*%t(S)
  }

  #Jacobiano
  p      <- length(betas)
  j11    <- diag(rep(1,p))
  j12    <- matrix(0, nrow=p, ncol=g-1)
  j13    <- matrix(0, nrow=p, ncol=1)
  j14    <- matrix(0, nrow=p, ncol=g)
  j15    <- matrix(0, nrow=p, ncol=g)
  j16    <- matrix(0, nrow=p, ncol=g)
  J1     <- cbind(j11,j12,j13,j14,j15,j16)

  j21    <- matrix(0, nrow=g-1, ncol=p)
  j22    <- diag(rep(1,g-1))
  j23    <- matrix(0, nrow=g-1, ncol=1)
  j24    <- -as.matrix((medj[1:(g-1)] - medj[g])^(-1))%*%t(as.matrix(pii))
  j24    <- matrix(j24,ncol=g,nrow=g-1)
  j25    <- matrix(0, nrow=g-1, ncol=g)
  j26    <- matrix(0, nrow=g-1, ncol=g)
  J2     <- cbind(j21,j22,j23,j24,j25,j26)

  j31    <- matrix(0, nrow=g, ncol=p)
  j32    <- matrix(0, nrow=g, ncol=g-1)
  j33    <- matrix(1, nrow=g, ncol=1)
  j34    <- diag(rep(1,g))
  j35    <- matrix(0, nrow=g, ncol=g)
  j36    <- matrix(0, nrow=g, ncol=g)
  J3     <- cbind(j31,j32,j33,j34,j35,j36)

  j41    <- matrix(0, nrow=g, ncol=p)
  j42    <- matrix(0, nrow=g, ncol=g-1)
  j43    <- matrix(0, nrow=g, ncol=1)
  j44    <- matrix(0, nrow=g, ncol=g)
  j45    <- diag(rep(1,g))
  j46    <- matrix(0, nrow=g, ncol=g)
  J4     <- cbind(j41,j42,j43,j44,j45,j46)

  j51    <- matrix(0, nrow=g, ncol=p)
  j52    <- matrix(0, nrow=g, ncol=g-1)
  j53    <- matrix(0, nrow=g, ncol=1)
  j54    <- matrix(0, nrow=g, ncol=g)
  j55    <- matrix(0, nrow=g, ncol=g)
  j56    <- diag(rep(1,g))
  J5     <- cbind(j51,j52,j53,j54,j55,j56)

  Jacobian       <- rbind(J1,J2,J3,J4,J5)
  IM             <- t(Jacobian)%*%solve(soma)%*%Jacobian
  EP             <- as.matrix(sqrt(diag(t(Jacobian)%*%solve(soma)%*%Jacobian)))

  namesrowBetas  <- c(); for(i in 1:length(betas)){namesrowBetas[i]  <- paste("beta",i,sep="")}
  namesrowMedj   <- c(); for(i in 1:g)            {namesrowMedj[i]   <- paste("mu",i,sep="")}
  namesrowSigmas <- c(); for(i in 1:g)            {namesrowSigmas[i] <- paste("sigma",i,sep="")}
  namesrowShape  <- c(); for(i in 1:g)            {namesrowShape[i]  <- paste("shape",i,sep="")}
  namesrowPii    <- c(); for(i in 1:(g-1))        {namesrowPii[i]    <- paste("pii",i,sep="")}
  rownames(EP)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
  colnames(EP)   <- c("SE")

  colnames(IM)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
  rownames(IM)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
  }

  if(class(model)=="Skew.cn")
  {

    k1                 <- nu[1]/nu[2]^(1/2)+1-nu[1] #Skew.cn
    b                  <- -sqrt(2/pi)*k1 #Skew.cn

    mu=mu1         <- matrix(0,n,g)
    varnu          <- rep(0,g) # beta0 + mu_j #Parametrizacao
    for (k in 1:g)
    {
      varnu[k]     <- beta0+medj[k]
      mu1[,k]      <- x%*%betas
      mu[,k]       <- mu1[,k]+varnu[k] + b*Delta[k]
    }

    #Verosimilhanca da mixtura da Skew.T
    lk             <- d.mixedSNC(y, pii, mu, sigma2, shape, nu)

    #Funcao do Skew.cn
    I.Phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric( sqrt(2*pi)*(nu[1]*nu[2]^(w -0.5)*dnorm(sqrt(di), 0, sqrt(1/nu[2]))*pnorm(nu[2]^(1/2)*Ai) + (1 - nu[1])*(dnorm(sqrt(di), 0,1)*pnorm(Ai)) )   )
    I.phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric( nu[1]*nu[2]^(w - 0.5)*dnorm(sqrt(di + Ai^2), 0, sqrt(1/nu[2])) + (1 - nu[1])*dnorm(sqrt(di + Ai^2))   )

    ################################
    #Definiendo os objetos Si
    Sibetas        <- matrix(0,length(betas),n)
    Sivarnu        <- matrix(0,g,n)
    Sisigma2       <- matrix(0,g,n)
    Sishape        <- matrix(0,g,n)
    Sipi           <- matrix(0,g-1,n)
    ################################

    #Nesta parte estou obtenido Sibetas onde considero "beta1", "beta2", pois na derivada respecto de "varnu" ja esta incluso o "beta0"
    #pelo que a derivada respecto a "beta" seria sem considera "beta0"
    #Lembrar que para cada "i" teremos a suma das "g" componentes
    for(i in 1:n)
    {
      Dr           <- sqrt(sigma2)
      B            <- y - mu
      di           <- B[i,]^2/Dr^2
      Ai           <- shape*B[i,]/Dr

      #Para betas
      dmu1         <- B[i,]/Dr^3
      dmu2         <- shape/Dr^2
      for(j in 1:g) Sibetas[,i]    <- Sibetas[,i] + (2/sqrt(2*pi))*pii[j]*(dmu1[j]*I.Phi(3/2, Ai[j], di[j], nu) - dmu2[j]*I.phi(1, Ai[j], di[j], nu))*x[i,]/lk[i]
      for(j in 1:(g-1)) Sipi[j,i]  <- (dSNC(y[i], mu[i,j], sigma2[j], shape[j], nu) - dSNC(y[i], mu[i,g], sigma2[g], shape[g], nu))/lk[i]
    }

    #Nesta parte vamos obter os demais Si, onde tambem es considerado o Sibetas (de tres componentes, pois temos 3 variaveis que incluye o intercepto)
    soma           <- 0
    for(i in 1:n)
    {
      Dr           <- sqrt(sigma2)
      B            <- y - mu
      di           <- B[i,]^2/Dr^2
      Ai           <- shape*B[i,]/Dr
      ################################################################
      #Para varnu
      dmu1         <- B[i,]/Dr^3
      dmu2         <- shape/Dr^2

      #Para sigma
      dsigma1      <- 1/Dr^3
      dsigma2      <- B[i,]*(B[i,] + b*Delta)/Dr^5
      dsigma3      <- shape*(B[i,] + b*Delta)/Dr^4

      #Para os elementos de lambda
      dlambda1     <- (B[i,]*b)/(Dr^2*(1+shape^2)^(3/2))
      dlambda2     <- (1/Dr^2)*(B[i,] - b*Delta/(1+shape^2))
      ###############################################################
      for(j in 1:g)
      {
       Sivarnu[j,i]  <- (2/sqrt(2*pi))*as.matrix((pii[j]/lk[i])*( dmu1[j]     * I.Phi(3/2, Ai[j], di[j], nu)     - dmu2[j] * I.phi(1, Ai[j], di[j], nu)))
       Sishape[j,i]  <- (2/sqrt(2*pi))*as.matrix((pii[j]/lk[i])*( dlambda1[j] * I.Phi(3/2, Ai[j], di[j], nu) + dlambda2[j] * I.phi(1 , Ai[j], di[j], nu)))
       Sisigma2[j,i] <- (1/sqrt(2*pi))*as.matrix((pii[j]/lk[i])*(-dsigma1[j]  * I.Phi(1/2, Ai[j], di[j], nu) + dsigma2[j]  * I.Phi(3/2, Ai[j], di[j], nu) -  dsigma3[j]*I.phi(1, Ai[j], di[j], nu)))
      }
      S            <- c(Sibetas[,i],Sipi[,i],Sivarnu[,i],Sisigma2[,i],Sishape[,i])
      soma         <- soma + S%*%t(S)
    }

    #Jacobiano
    p      <- length(betas)
    j11    <- diag(rep(1,p))
    j12    <- matrix(0, nrow=p, ncol=g-1)
    j13    <- matrix(0, nrow=p, ncol=1)
    j14    <- matrix(0, nrow=p, ncol=g)
    j15    <- matrix(0, nrow=p, ncol=g)
    j16    <- matrix(0, nrow=p, ncol=g)
    J1     <- cbind(j11,j12,j13,j14,j15,j16)

    j21    <- matrix(0, nrow=g-1, ncol=p)
    j22    <- diag(rep(1,g-1))
    j23    <- matrix(0, nrow=g-1, ncol=1)
    j24    <- -as.matrix((medj[1:(g-1)] - medj[g])^(-1))%*%t(as.matrix(pii))
    j24    <- matrix(j24,ncol=g,nrow=g-1)
    j25    <- matrix(0, nrow=g-1, ncol=g)
    j26    <- matrix(0, nrow=g-1, ncol=g)
    J2     <- cbind(j21,j22,j23,j24,j25,j26)

    j31    <- matrix(0, nrow=g, ncol=p)
    j32    <- matrix(0, nrow=g, ncol=g-1)
    j33    <- matrix(1, nrow=g, ncol=1)
    j34    <- diag(rep(1,g))
    j35    <- matrix(0, nrow=g, ncol=g)
    j36    <- matrix(0, nrow=g, ncol=g)
    J3     <- cbind(j31,j32,j33,j34,j35,j36)

    j41    <- matrix(0, nrow=g, ncol=p)
    j42    <- matrix(0, nrow=g, ncol=g-1)
    j43    <- matrix(0, nrow=g, ncol=1)
    j44    <- matrix(0, nrow=g, ncol=g)
    j45    <- diag(rep(1,g))
    j46    <- matrix(0, nrow=g, ncol=g)
    J4     <- cbind(j41,j42,j43,j44,j45,j46)

    j51    <- matrix(0, nrow=g, ncol=p)
    j52    <- matrix(0, nrow=g, ncol=g-1)
    j53    <- matrix(0, nrow=g, ncol=1)
    j54    <- matrix(0, nrow=g, ncol=g)
    j55    <- matrix(0, nrow=g, ncol=g)
    j56    <- diag(rep(1,g))
    J5     <- cbind(j51,j52,j53,j54,j55,j56)

    Jacobian        <- rbind(J1,J2,J3,J4,J5)
    IM              <- t(Jacobian)%*%solve(soma)%*%Jacobian
    EP             <- as.matrix(sqrt(diag(t(Jacobian)%*%solve(soma)%*%Jacobian)))

    namesrowBetas  <- c(); for(i in 1:length(betas)){namesrowBetas[i]  <- paste("beta",i,sep="")}
    namesrowMedj   <- c(); for(i in 1:g)            {namesrowMedj[i]   <- paste("mu",i,sep="")}
    namesrowSigmas <- c(); for(i in 1:g)            {namesrowSigmas[i] <- paste("sigma",i,sep="")}
    namesrowShape  <- c(); for(i in 1:g)            {namesrowShape[i]  <- paste("shape",i,sep="")}
    namesrowPii    <- c(); for(i in 1:(g-1))        {namesrowPii[i]    <- paste("pii",i,sep="")}
    rownames(EP)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
    colnames(EP)   <- c("SE")

    colnames(IM)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
    rownames(IM)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
  }

  if(class(model)=="Skew.slash")
  {

    k1                <- 2*nu/(2*nu-1) #Skew.slash
    b                 <- -sqrt(2/pi)*k1 #Skew.slash

    mu=mu1         <- matrix(0,n,g)
    varnu          <- rep(0,g) # beta0 + mu_j #Parametrizacao
    for (k in 1:g)
    {
      varnu[k]     <- beta0+medj[k]
      mu1[,k]      <- x%*%betas
      mu[,k]       <- mu1[,k]+varnu[k] + b*Delta[k]
    }

    #Verosimilhanca da mixtura da Skew.T
    lk             <- d.mixedSS(y, pii, mu, sigma2, shape, nu)

    #Funcao do Skew.slash
    I.Phi  <- function(w=0, Ai=0, di=0, nu=0)
    {
      f    <- function(u) pnorm(u^(0.5)*Ai)*dgamma(u,nu+w,di/2)
      resp <- integrate(f,0,1)$value
      res1 <- (nu*(2^(w + nu)*gamma(w + nu))/(di^(w + nu)))*resp
      return(res1)
    }

    I.phi  <- function(w=0, Ai=0, di=0, nu=0)
    {
      res2 <- ((nu*2^(w + nu)*gamma(w + nu))/(sqrt(2*pi)*(di + Ai^2)^(w + nu)))*pgamma(1, w + nu, (di + Ai^2)/2)
      return(res2)
    }

    ################################
    #Definiendo os objetos Si
    Sibetas        <- matrix(0,length(betas),n)
    Sivarnu        <- matrix(0,g,n)
    Sisigma2       <- matrix(0,g,n)
    Sishape        <- matrix(0,g,n)
    Sipi           <- matrix(0,g-1,n)
    ################################

    #Nesta parte estou obtenido Sibetas onde considero "beta1", "beta2", pois na derivada respecto de "varnu" ja esta incluso o "beta0"
    #pelo que a derivada respecto a "beta" seria sem considera "beta0"
    #Lembrar que para cada "i" teremos a suma das "g" componentes
    for(i in 1:n)
    {
      Dr           <- sqrt(sigma2)
      B            <- y - mu
      di           <- B[i,]^2/Dr^2
      Ai           <- shape*B[i,]/Dr

      #Para betas
      dmu1         <- B[i,]/Dr^3
      dmu2         <- shape/Dr^2
      for(j in 1:g) Sibetas[,i]    <- Sibetas[,i] + (2/sqrt(2*pi))*pii[j]*(dmu1[j]*I.Phi(3/2, Ai[j], di[j], nu) - dmu2[j]*I.phi(1, Ai[j], di[j], nu))*x[i,]/lk[i]
      for(j in 1:(g-1)) Sipi[j,i]  <- (dSS(y[i], mu[i,j], sigma2[j], shape[j], nu) - dSS(y[i], mu[i,g], sigma2[g], shape[g], nu))/lk[i]
    }

    #Nesta parte vamos obter os demais Si, onde tambem es considerado o Sibetas (de tres componentes, pois temos 3 variaveis que incluye o intercepto)
    soma           <- 0
    for(i in 1:n)
    {
      Dr           <- sqrt(sigma2)
      B            <- y - mu
      di           <- B[i,]^2/Dr^2
      Ai           <- shape*B[i,]/Dr
      ################################################################
      #Para varnu
      dmu1         <- B[i,]/Dr^3
      dmu2         <- shape/Dr^2

      #Para sigma
      dsigma1      <- 1/Dr^3
      dsigma2      <- B[i,]*(B[i,] + b*Delta)/Dr^5
      dsigma3      <- shape*(B[i,] + b*Delta)/Dr^4

      #Para os elementos de lambda
      dlambda1     <- (B[i,]*b)/(Dr^2*(1+shape^2)^(3/2))
      dlambda2     <- (1/Dr^2)*(B[i,] - b*Delta/(1+shape^2))
      ###############################################################
      for(j in 1:g)
      {
       Sivarnu[j,i]  <- (2/sqrt(2*pi))*as.matrix((pii[j]/lk[i])*( dmu1[j]     * I.Phi(3/2, Ai[j], di[j], nu)     - dmu2[j] * I.phi(1, Ai[j], di[j], nu)))
       Sishape[j,i]  <- (2/sqrt(2*pi))*as.matrix((pii[j]/lk[i])*( dlambda1[j] * I.Phi(3/2, Ai[j], di[j], nu) + dlambda2[j] * I.phi(1 , Ai[j], di[j], nu)))
       Sisigma2[j,i] <- (1/sqrt(2*pi))*as.matrix((pii[j]/lk[i])*(-dsigma1[j]  * I.Phi(1/2, Ai[j], di[j], nu) + dsigma2[j]  * I.Phi(3/2, Ai[j], di[j], nu) -  dsigma3[j]*I.phi(1, Ai[j], di[j], nu)))
      }

      S            <- c(Sibetas[,i],Sipi[,i],Sivarnu[,i],Sisigma2[,i],Sishape[,i])
      soma         <- soma + S%*%t(S)
    }

    #Jacobiano
    p      <- length(betas)
    j11    <- diag(rep(1,p))
    j12    <- matrix(0, nrow=p, ncol=g-1)
    j13    <- matrix(0, nrow=p, ncol=1)
    j14    <- matrix(0, nrow=p, ncol=g)
    j15    <- matrix(0, nrow=p, ncol=g)
    j16    <- matrix(0, nrow=p, ncol=g)
    J1     <- cbind(j11,j12,j13,j14,j15,j16)

    j21    <- matrix(0, nrow=g-1, ncol=p)
    j22    <- diag(rep(1,g-1))
    j23    <- matrix(0, nrow=g-1, ncol=1)
    j24    <- -as.matrix((medj[1:(g-1)] - medj[g])^(-1))%*%t(as.matrix(pii))
    j24    <- matrix(j24,ncol=g,nrow=g-1)
    j25    <- matrix(0, nrow=g-1, ncol=g)
    j26    <- matrix(0, nrow=g-1, ncol=g)
    J2     <- cbind(j21,j22,j23,j24,j25,j26)

    j31    <- matrix(0, nrow=g, ncol=p)
    j32    <- matrix(0, nrow=g, ncol=g-1)
    j33    <- matrix(1, nrow=g, ncol=1)
    j34    <- diag(rep(1,g))
    j35    <- matrix(0, nrow=g, ncol=g)
    j36    <- matrix(0, nrow=g, ncol=g)
    J3     <- cbind(j31,j32,j33,j34,j35,j36)

    j41    <- matrix(0, nrow=g, ncol=p)
    j42    <- matrix(0, nrow=g, ncol=g-1)
    j43    <- matrix(0, nrow=g, ncol=1)
    j44    <- matrix(0, nrow=g, ncol=g)
    j45    <- diag(rep(1,g))
    j46    <- matrix(0, nrow=g, ncol=g)
    J4     <- cbind(j41,j42,j43,j44,j45,j46)

    j51    <- matrix(0, nrow=g, ncol=p)
    j52    <- matrix(0, nrow=g, ncol=g-1)
    j53    <- matrix(0, nrow=g, ncol=1)
    j54    <- matrix(0, nrow=g, ncol=g)
    j55    <- matrix(0, nrow=g, ncol=g)
    j56    <- diag(rep(1,g))
    J5     <- cbind(j51,j52,j53,j54,j55,j56)

    Jacobian       <- rbind(J1,J2,J3,J4,J5)
    IM             <- t(Jacobian)%*%solve(soma)%*%Jacobian
    EP             <- as.matrix(sqrt(diag(t(Jacobian)%*%solve(soma)%*%Jacobian)))

    namesrowBetas  <- c(); for(i in 1:length(betas)){namesrowBetas[i]  <- paste("beta",i,sep="")}
    namesrowMedj   <- c(); for(i in 1:g)            {namesrowMedj[i]   <- paste("mu",i,sep="")}
    namesrowSigmas <- c(); for(i in 1:g)            {namesrowSigmas[i] <- paste("sigma",i,sep="")}
    namesrowShape  <- c(); for(i in 1:g)            {namesrowShape[i]  <- paste("shape",i,sep="")}
    namesrowPii    <- c(); for(i in 1:(g-1))        {namesrowPii[i]    <- paste("pii",i,sep="")}
    rownames(EP)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
    colnames(EP)   <- c("SE")

    colnames(IM)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
    rownames(IM)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
  }

  if(class(model)=="Skew.normal")
  {

    k1             <- 1  #Skew.normal
    b              <- -sqrt(2/pi)*k1 #Skew.normal

    mu=mu1         <- matrix(0,n,g)
    varnu          <- rep(0,g) # beta0 + mu_j #Parametrizacao
    for (k in 1:g)
    {
      varnu[k]     <- beta0+medj[k]
      mu1[,k]      <- x%*%betas
      mu[,k]       <- mu1[,k]+varnu[k] + b*Delta[k]
    }

    #Verosimilhanca da mixtura da Skew.T
    lk             <- d.mixedSN(y, pii, mu, sigma2, shape)

    #Funcao do Skew.normal
    I.Phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric( exp(-di/2)*pnorm(Ai) )
    I.phi <- function(w=0, Ai=0, di=0, nu=0) as.numeric( exp(-di/2)*dnorm(Ai) )

    ################################
    #Definiendo os objetos Si
    Sibetas        <- matrix(0,length(betas),n)
    Sivarnu        <- matrix(0,g,n)
    Sisigma2       <- matrix(0,g,n)
    Sishape        <- matrix(0,g,n)
    Sipi           <- matrix(0,g-1,n)
    ################################

    #Nesta parte estou obtenido Sibetas onde considero "beta1", "beta2", pois na derivada respecto de "varnu" ja esta incluso o "beta0"
    #pelo que a derivada respecto a "beta" seria sem considera "beta0"
    #Lembrar que para cada "i" teremos a suma das "g" componentes
    for(i in 1:n)
    {
      Dr           <- sqrt(sigma2)
      B            <- y - mu
      di           <- B[i,]^2/Dr^2
      Ai           <- shape*B[i,]/Dr

      #Para betas
      dmu1         <- B[i,]/Dr^3
      dmu2         <- shape/Dr^2
      for(j in 1:g) Sibetas[,i]    <- Sibetas[,i] + (2/sqrt(2*pi))*pii[j]*(dmu1[j]*I.Phi(Ai[j], di[j]) - dmu2[j]*I.phi(Ai[j], di[j]))*x[i,]/lk[i]
      for(j in 1:(g-1)) Sipi[j,i]  <- (dSN(y[i], mu[i,j], sigma2[j], shape[j]) - dSN(y[i], mu[i,g], sigma2[g], shape[g]))/lk[i]
    }

    #Nesta parte vamos obter os demais Si, onde tambem es considerado o Sibetas (de tres componentes, pois temos 3 variaveis que incluye o intercepto)
    soma           <- 0
    for(i in 1:n)
    {
      Dr           <- sqrt(sigma2)
      B            <- y - mu
      di           <- B[i,]^2/Dr^2
      Ai           <- shape*B[i,]/Dr
      ################################################################
      #Para varnu
      dmu1         <- B[i,]/Dr^3
      dmu2         <- shape/Dr^2

      #Para sigma
      dsigma1      <- 1/Dr^3
      dsigma2      <- B[i,]*(B[i,] + b*Delta)/Dr^5
      dsigma3      <- shape*(B[i,] + b*Delta)/Dr^4

      #Para os elementos de lambda
      dlambda1     <- (B[i,]*b)/(Dr^2*(1+shape^2)^(3/2))
      dlambda2     <- (1/Dr^2)*(B[i,] - b*Delta/(1+shape^2))
      ###############################################################

      Sivarnu[,i]  <- (2/sqrt(2*pi))*as.matrix((pii/lk[i])*( dmu1     * I.Phi(Ai, di)     - dmu2 * I.phi(Ai, di)))
      Sishape[,i]  <- (2/sqrt(2*pi))*as.matrix((pii/lk[i])*( dlambda1 * I.Phi(Ai, di) + dlambda2 * I.phi(Ai, di)))
      Sisigma2[,i] <- (1/sqrt(2*pi))*as.matrix((pii/lk[i])*(-dsigma1  * I.Phi(Ai, di) + dsigma2  * I.Phi(Ai, di) -  dsigma3*I.phi(Ai, di)))

      S            <- c(Sibetas[,i],Sipi[,i],Sivarnu[,i],Sisigma2[,i],Sishape[,i])
      soma         <- soma + S%*%t(S)
    }

    #Jacobiano
    p      <- length(betas)
    j11    <- diag(rep(1,p))
    j12    <- matrix(0, nrow=p, ncol=g-1)
    j13    <- matrix(0, nrow=p, ncol=1)
    j14    <- matrix(0, nrow=p, ncol=g)
    j15    <- matrix(0, nrow=p, ncol=g)
    j16    <- matrix(0, nrow=p, ncol=g)
    J1     <- cbind(j11,j12,j13,j14,j15,j16)

    j21    <- matrix(0, nrow=g-1, ncol=p)
    j22    <- diag(rep(1,g-1))
    j23    <- matrix(0, nrow=g-1, ncol=1)
    j24    <- -as.matrix((medj[1:(g-1)] - medj[g])^(-1))%*%t(as.matrix(pii))
    j24    <- matrix(j24,ncol=g,nrow=g-1)
    j25    <- matrix(0, nrow=g-1, ncol=g)
    j26    <- matrix(0, nrow=g-1, ncol=g)
    J2     <- cbind(j21,j22,j23,j24,j25,j26)

    j31    <- matrix(0, nrow=g, ncol=p)
    j32    <- matrix(0, nrow=g, ncol=g-1)
    j33    <- matrix(1, nrow=g, ncol=1)
    j34    <- diag(rep(1,g))
    j35    <- matrix(0, nrow=g, ncol=g)
    j36    <- matrix(0, nrow=g, ncol=g)
    J3     <- cbind(j31,j32,j33,j34,j35,j36)

    j41    <- matrix(0, nrow=g, ncol=p)
    j42    <- matrix(0, nrow=g, ncol=g-1)
    j43    <- matrix(0, nrow=g, ncol=1)
    j44    <- matrix(0, nrow=g, ncol=g)
    j45    <- diag(rep(1,g))
    j46    <- matrix(0, nrow=g, ncol=g)
    J4     <- cbind(j41,j42,j43,j44,j45,j46)

    j51    <- matrix(0, nrow=g, ncol=p)
    j52    <- matrix(0, nrow=g, ncol=g-1)
    j53    <- matrix(0, nrow=g, ncol=1)
    j54    <- matrix(0, nrow=g, ncol=g)
    j55    <- matrix(0, nrow=g, ncol=g)
    j56    <- diag(rep(1,g))
    J5     <- cbind(j51,j52,j53,j54,j55,j56)

    Jacobian       <- rbind(J1,J2,J3,J4,J5)
    IM             <- t(Jacobian)%*%solve(soma)%*%Jacobian
    EP             <- as.matrix(sqrt(diag(t(Jacobian)%*%solve(soma)%*%Jacobian)))

    namesrowBetas  <- c(); for(i in 1:length(betas)){namesrowBetas[i]  <- paste("beta",i,sep="")}
    namesrowMedj   <- c(); for(i in 1:g)            {namesrowMedj[i]   <- paste("mu",i,sep="")}
    namesrowSigmas <- c(); for(i in 1:g)            {namesrowSigmas[i] <- paste("sigma",i,sep="")}
    namesrowShape  <- c(); for(i in 1:g)            {namesrowShape[i]  <- paste("shape",i,sep="")}
    namesrowPii    <- c(); for(i in 1:(g-1))        {namesrowPii[i]    <- paste("pii",i,sep="")}
    rownames(EP)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
    colnames(EP)   <- c("SE")

    colnames(IM)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
    rownames(IM)   <- c(namesrowBetas,namesrowPii,"beta0",namesrowMedj,namesrowSigmas,namesrowShape)
  }

  return(list(IM=IM,EP=EP))
}

#im.FMsmsnReg2(y, x1, model=parSS)
#model=TrySMSN
#Jacobiano



