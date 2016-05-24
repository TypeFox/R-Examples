
statBO <- function(data,formula,family,groupRef,groupInd,indice){
  
  call <- match.call()
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data)) 
    data <- environment(formula)
  
  groupInd<-groupInd[indice]
  data<-data[indice,]
  dataA <- subset(data,groupInd==1)
  nA <- length(dataA[,1])
  
  
  dataB <- subset(data,groupInd==0)
  nB <- length(dataB[,1])
  
  
  
  # glm model fore each group
  glmA <- stats::glm(formula, family, data=dataA, x=TRUE)
  betaA <- glmA$coefficients
  Xa <- glmA$x
  EbetaAXa <- mean(glmA$fitted.values)
  
  
  glmB <- stats::glm(formula, family, data=dataB,x=TRUE)
  betaB <- glmB$coefficients
  Xb <- glmB$x
  EbetaBXb <- mean(glmB$fitted.values)
  
  
  regoutput <- list(GroupA=glmA, GroupB=glmB)
  
  # function that compute approximation of the expected mean
  expect <- function(X,beta) mean(stats::make.link(family$link)$linkinv((X)%*%beta))
  EbetaAXb <- expect(Xb,betaA)
  EbetaBXa <- expect(Xa,betaB)
  
  # twofold decomposition
  
  if (groupRef=="A"){
    char= EbetaAXa- EbetaAXb
    coef= EbetaAXb-EbetaBXb
  }
  if (groupRef=="B"){
    char= EbetaBXa- EbetaBXb
    coef= EbetaAXa-EbetaBXa
  }
  diff=char+coef
  
  
  # threefold decomposition
  
  char3=EbetaBXa-EbetaBXb
  coef3=EbetaAXb-EbetaBXb
  int=EbetaAXa-EbetaBXa+EbetaAXb-EbetaBXb
  
  
  result<-c(char,coef,diff,char3,coef3,int)
  #out<-list(char=char,coef=coef,diff=diff,char3=char3,coef3=coef3,int=int)
  return(result)
  
  
}
