# Blinder-Oaxaca decomposition for generalized linear model
#
# Date: 22/07/15
# Auteur: Aurelien Nicosia and Simon Baillargeon-Ladouceur
#
###############################################################################

GeneralOaxaca <- function(formula,  family= stats::gaussian, data, groupInd, groupRef="A", B=1000, control=list()){

  


  
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
  
  # subset of each group
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
      if (requireNamespace("boot", quietly = TRUE)) {
        Bot<-boot::boot(data=data, statistic=statBO, R=B, 
                  formula=formula,
                  family=family, groupInd = groupInd, groupRef="A"
        )
      }
      
      
    }
    if (groupRef=="B"){
      char= EbetaBXa- EbetaBXb
      coef= EbetaAXa-EbetaBXa
      if (requireNamespace("boot", quietly = TRUE)) {
        Bot<-boot::boot(data=data, statistic=statBO, R=B, 
                        formula=formula,
                        family=family, groupInd = groupInd, groupRef="B"
        )}
    }
    
    # bootstrap estimation of s.e.
    SD<-apply(Bot$t,2,stats::sd)
    sechar <-SD[1]
    secoef <- SD[2]
    sediff <- SD[3]
    
    
    diff= char+coef
    prop = c(char/diff,coef/diff,diff/diff)*100
    C<-c(char,coef,diff)
    SEC<-c(sechar,secoef,sediff)
    ZVAL<-c(char,coef,diff)/c(sechar,secoef,sediff)
    twofold<- cbind(C,prop,SEC,ZVAL,round(1-stats::pnorm(abs(ZVAL)),4))
    
    ic2<-NULL
    for (i in (1:3)){
      if (requireNamespace("boot", quietly = TRUE)) {
      ic2<-rbind(ic2,boot::boot.ci(Bot,type= "norm",index=i)$normal[c(2,3)])}
      
    }
    twofold<-cbind(twofold,ic2)
    
    colnames(twofold) <- c("value","prop (%)","s.e","z value","P(>|z|)","CI_l","CI_u")
    rownames(twofold) <- c("char","coeff","diff tot")
  
  # threefold decomposition
    char3=EbetaBXa-EbetaBXb
    coef3=EbetaAXb-EbetaBXb
    int=EbetaAXa-EbetaBXa+EbetaAXb-EbetaBXb
    
    sechar3 <- SD[4]
      secoef3 <- SD[5]
      seint <-  SD[6]
  
      diff= char3+coef3+int
      prop = c(char3/diff,coef3/diff,int/diff)*100
      C<-c(char3,coef3,int)
      SEC<-c(sechar3,secoef3,seint)
      ZVAL<-c(char3,coef3,int)/c(sechar3,secoef3,seint)
      threefold<- cbind(C,prop,SEC,ZVAL,round(1-stats::pnorm(abs(ZVAL)),4))
      ic3<-NULL
      for (i in (4:6)){
        if (requireNamespace("boot", quietly = TRUE)) {
        ic3<-rbind(ic3,boot::boot.ci(Bot,type= "norm",index=i)$normal[c(2,3)])}
        
      }
      threefold<-cbind(threefold,ic3)
      colnames(threefold) <- c("value","prop (%)","s.e","z value","P(>|z|)","CI_l","CI_u")
      rownames(threefold) <- c("char","coeff","int")
      
      # other usefull informations
      n<-list(nA=nA,nB=nB)
      summaryStat<-list(summaryA=summary(glmA$y),summaryB=summary(glmB$y),meandiff=mean(glmA$y)-mean(glmB$y))
      
      
  out <- list(regoutput=regoutput, twofold=twofold,  threefold=threefold,n=n,summaryStat=summaryStat)
  class(out) <- "Blinder-Oaxaca"
  out
}

  