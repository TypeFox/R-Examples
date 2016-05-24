##################################################################
##########     Funcoes para gerar SNI e mistura       ############


  rmix <- function(n, pii, family, arg, cluster=FALSE) {
    #Funcao para gerar misturas de g populacoes
    #n: numero de amostras geradas
    #p: vetor de proporcoes das misturas (tamanho g)
    #arg: deve ser um tipo list com cada entrada contendo um vetor de tamanho g de agrumentos a ser passado para rF1
 if((family != "t") && (family != "Skew.t") && (family != "Skew.cn") && (family != "Skew.slash") && (family != "Skew.normal") && (family != "Normal")) stop(paste("Family",family,"not recognized.",sep=" "))

    if((family == "Normal") || (family == "Skew.normal") ) {
                             rF1 <- gen.Skew.normal
                             for (i in 1:length(arg)) if(length(arg[[i]]) != 4 && length(arg[[i]]) != 3) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                             if(family == "Normal") for (i in 1:length(arg)) arg[[i]][3] <- 0                              }
                                                           
    if ((family == "t") || (family == "Skew.t") ){
                             rF1 <- gen.Skew.t
                             for (i in 1:length(arg)) if(length(arg[[i]]) != 4) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                             if(family == "t") for (i in 1:length(arg)) arg[[i]][3] <- 0
                           }
    if (family == "Skew.cn"){
                              rF1 <- gen.Skew.cn
                              for (i in 1:length(arg)) if(length(arg[[i]]) != 5) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                            }
    if (family == "Skew.slash") {
                                 rF1 <- gen.Skew.slash
                                 for (i in 1:length(arg)) if(length(arg[[i]]) != 4) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                                }

    x1 <- vector(mode = "numeric", length = n)
    clu <- vector(mode = "numeric", length = n)
    g <- length(pii)
    interval <- c(0)
    for (j in 1:g-1) interval <- cbind(interval, interval[j] + pii[j])
    interval <- cbind(interval, 1)
    for(i in 1:n) {
      u <- runif(1)
      clu[i] <- findInterval(u, interval)
      x1[i] <- do.call("rF1", c(list(1), arg[[clu[i]]]))
      
    }
    if(cluster) return(list(y=x1, cluster=clu))
    else return(x1)
  }


  gen.Skew.normal <- function(n, mu, sigma2, shape, nu=NULL){
    #Funcao para gerar valores aleatorios de uma Skew-Normal
    #n: qtd de valores a ser gerado
    #mu, sigma2, shape: locacao, escala e assimetria, respectivamente
    delta <- shape / sqrt(1 + shape^2)
    y <- mu*rep(1,n) + sqrt(sigma2)*(delta*abs(rnorm(n)) + (1 - delta^2)^(1/2)*rnorm(n))
    return(y)
  }

  gen.Skew.t <- function(n, mu, sigma2, shape, nu ){
    #Funcao para gerar Skew-t
    #n: qtd de valores a ser gerado
    #mu, sigma2, shape: locacao, escala e assimetria, respectivamente
    y <- mu + (rgamma(n, nu/2, nu/2))^(-1/2)*gen.Skew.normal(n, 0, sigma2, shape)
  }


  gen.Skew.cn <- function(n, mu, sigma2, shape, nu1, nu2){
    #Funcao para gerar Skew Normal Contaminada
    #n: qtd de valores a ser gerado
    #mu, sigma2, shape: locacao, escala e assimetria, respectivamente
    #rmix(n, nu[1], gen.Skew.normal, list(c(mu,sigma2/nu[2],shape), c(mu,sigma2,shape)))
    rmix(n, c(nu1,1-nu1), "Skew.normal", list(c(mu,sigma2/nu2,shape), c(mu,sigma2,shape)))
  }


  gen.Skew.slash <- function(n, mu, sigma2, shape, nu){
    # Funcao para gerar Skew Slash
    #n: qtd de valores a ser gerado
    #mu, sigma2, shape: locacao, escala e assimetria, respectivamente
    u1 <- runif(n)
    u2 <- u1^(1/(nu))   # formula 10 do artigo e metodo da inversao
    ys <- mu + (u2)^(-1/2)*gen.Skew.normal(n, 0, sigma2, shape)
    return(ys)
  }


##########  FIM   Funcoes para gerar SNI e mistura     ###########
##################################################################

################################################################
##########      Funcoes para numeros aleatorios     ############

  rmmix <- function(n, pii, family, arg, cluster=FALSE) {
  ##require(mvtnorm)
    #Funcao para gerar misturas de g populacoes
    #n: numero de amostras geradas
    #p: vetor de proporcoes das misturas (tamanho g)
    #arg: deve ser um tipo list com cada entrada contendo um vetor de tamanho g de agrumentos a ser passado para rF1
    if((family != "t") && (family != "Skew.t") && (family != "Skew.cn") && (family != "Skew.slash") && (family != "Skew.normal") && (family != "Normal")) stop(paste("Family",family,"not recognized.",sep=" "))

    if (family == "Normal") {
                                rF1 <- gen.SN.multi
                                for (i in 1:length(arg)){
                                   arg[[i]]$shape <- rep(0,length(arg[[i]]$shape))     
                                   if(length(arg[[i]]) != 4 && length(arg[[i]]) != 3) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                                }
                             }
    if (family == "Skew.normal"){
                                   rF1 <- gen.SN.multi
                                   for (i in 1:length(arg)) if(length(arg[[i]]) != 4 && length(arg[[i]]) != 3) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                                 }
    if ((family == "t") || (family == "Skew.t")){
                              rF1 <- gen.ST.multi
                              for (i in 1:length(arg)) if(length(arg[[i]]) != 4) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                              if(family == "t") for (i in 1:length(arg)) arg[[i]]$shape <- rep(0,length(arg[[i]]$shape))
                            }
    if (family == "Skew.cn"){
                               rF1 <- gen.SCN.multi
                               for (i in 1:length(arg)) if(length(arg[[i]]) != 4) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                             }
    if (family == "Skew.slash"){
                                  rF1 <- gen.SS.multi
                                  for (i in 1:length(arg)) if(length(arg[[i]]) != 4) stop(paste("Number of arguments is not comformidable for argument ",i,".\n",sep=" "))
                                }

    x1 <- matrix(data = NA, ncol = length(arg[[1]]$mu),nrow = n)
    clu <- vector(mode = "numeric", length = n)
    g <- length(pii)
    interval <- c(0)
    for (j in 1:g-1) interval <- cbind(interval, interval[j] + pii[j])
    interval <- cbind(interval, 1)
    for(i in 1:n) {
      u <- runif(1)
      clu[i] <- findInterval(u, interval)
      x1[i,] <- do.call("rF1", c(list(1), arg[[findInterval(u, interval)]]))
    }
    if(cluster) return(list(y=x1, cluster=clu))
    else return(x1)
  }

  gen.SN.multi <- function(n, mu, Sigma, shape, nu=NULL){
    p <- length(mu)
    delta <- shape / (sqrt(1 + t(shape)%*%shape))
    y <- matrix(0,n,p)
    for (i in 1:n) y[i,] <- mu + matrix.sqrt(Sigma)%*%(delta*abs(rnorm(1)) + matrix.sqrt(diag(p) - delta%*%t(delta))%*%as.vector(rmvnorm(1, mean = rep(0,p), sigma = diag(p))))
    return(y)
  }

  gen.ST.multi <- function(n, mu, Sigma, shape, nu){
    y <- matrix(rep(mu, n), n, length(mu), byrow = T) + (rgamma(n, nu/2, nu/2))^(-1/2)*gen.SN.multi(n, rep(0, length(mu)), Sigma, shape)
    return(y)
  }


  gen.SCN.multi  <- function(n, mu, Sigma, shape, nu){
    x1 <- matrix(0,n,length(mu))
    for (i in 1:n){
      u <- runif(1)
      if (u < nu[1]) x1[i,] <- gen.SN.multi(1, mu, Sigma/nu[2], shape)
      if (u > nu[1]) x1[i,] <- gen.SN.multi(1, mu, Sigma, shape)
    }
    return(x1)
  }

  gen.SS.multi  <- function(n, mu, Sigma, shape, nu){
    u1 <- runif(n)
    u2 <- u1^(1/(nu))   # formula 10 do artigo e metodo da inversao
    ys <- mu + (u2)^(-1/2)*gen.SN.multi(n, c(0,0), Sigma, shape)
    return(ys)
  }
