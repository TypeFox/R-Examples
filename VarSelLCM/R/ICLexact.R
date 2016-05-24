########################################################################################################################
## Fonction calculant la vraisemblance completee integree en zMAP (ICL exact)
########################################################################################################################

## fonction outils
logdinvgamma <- function(x, alpha, beta)
  alpha*log(beta) - lgamma(alpha) - (alpha+1)*log(x) - beta/x

## Integrale sur une variable continue pour 1 classe ou non discriminante
IntegreOneVariableContinuous <- function(x, priors){
  nu <- priors[1]
  s0 <- priors[2]
  mu <- priors[3]
  n0 <- priors[4]
  n <- length(x)
  n1 <- n + n0
  s1 <- sqrt( sum((x-mean(x))**2) + s0**2 + ((mu - mean(x))**2)/(1/n0 + 1/n) )
  # theta1 <- (n0*mu + n*mean(x))/n1
  integre <- -log(sqrt(pi))*n + lgamma((n + nu)/2) - lgamma(nu/2) +   nu * log(s0/s1) - n*log(s1) + log(sqrt(n0 / n1) )
  return(integre)
}


## ICL exact dans le cas de variables continues
ICLcontinuous <- function(obj){
  ICLexact <- lgamma(obj@model@g/2) - obj@model@g*lgamma(1/2) + 
    sum(lgamma(table(c(1:obj@model@g, obj@partitions@zMAP)) - 1/2)) - lgamma(length(obj@partitions@zMAP) + obj@model@g/2)
  for (j in 1:obj@data@d){
    if (obj@model@omega[j]==0){
      ICLexact <- ICLexact  + IntegreOneVariableContinuous(obj@data@data[which(obj@data@notNA[,j]==1),j], obj@data@priors[j,])
    }else{
      for (k in unique(obj@partitions@zMAP))
        ICLexact <- ICLexact  + IntegreOneVariableContinuous(obj@data@data[which( (obj@partitions@zMAP==k)*obj@data@notNA[,j] ==1),j], obj@data@priors[j,])
    }
  }
  names(ICLexact) <- NULL
  return(ICLexact)
}