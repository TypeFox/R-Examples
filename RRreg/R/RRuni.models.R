# RRuni: estimate pi and var(pi)

########################################
# Warner model (Warner, 1965)
RRuni.Warner<-function(response,p){
  n <- length(response)
  lambda <- sum(response)/n
  pi <- (lambda+p-1)/(2*p-1)
  piSE <- sqrt((lambda*(1-lambda))/(n*(2*p-1)^2))  
  res <- list(model="Warner",call=paste("Warner Model with p =",p),pi=pi,piSE=piSE,n=n)
  return(res)
} 

########################################
# Crosswise model (Yu, Tian, Tang, 2007)
RRuni.Crosswise<-function(response,p){
  n <- length(response)
  lambda <- sum(response)/n
  pi <- (lambda+p-1)/(2*p-1)
  piSE <- sqrt((lambda*(1-lambda))/(n*(2*p-1)^2))  
  res <- list(model="Crosswise",call=paste("Crosswise Model with p =",p),pi=pi,piSE=piSE,n=n)
  return(res)
} 

########################################
# Unrelated Question Technique (UQTknown) ()
RRuni.UQTknown<-function(response,p){
  n <- length(response)
  lambda <- sum(response)/n
  pi <- (lambda-(1-p[1])*p[2])/p[1]
  piSE <- sqrt((lambda*(1-lambda))/(n*p[1]^2))  
  pstring <- paste("'p1' =",p[1],"(probability of answering sensitive question) and 'p2' =",p[2],"(probability of answer 'Yes' in unrelated question")
  res <- list(model="UQTknown",call=paste0("Unrelated Question Technique with known prevalence of unrelated question. ",pstring),pi=pi,piSE=piSE,n=n)
  return(res)
} 

########################################
# Unrelated Question Technique (UQTunknown) (Greenberg et al. 1969)
RRuni.UQTunknown<-function(response,p,group){
  s <- group ==1
  n1 <- length(response[s])
  n2 <- length(response[!s])
  lambda1 <- sum(response[s])/n1
  lambda2 <- sum(response[!s])/n2
  pi <- (lambda1*(1-p[2])-lambda2*(1-p[1])) / (p[1]-p[2])
  piSE <- sqrt( (lambda1*(1-lambda1)*(1-p[2])^2/n1 + lambda2*(1-lambda2)*(1-p[1])^2/n2) / (p[1]-p[2])^2 ) 
  piUQ <- (p[2]*lambda1-p[1]*lambda2) / ( p[2]-p[1])
  piUQSE <- sqrt ( (p[2]^2*lambda1*(1-lambda1)/n1 +p[1]^2*lambda2*(1-lambda2)/n2)  / (p[1]-p[2])^2 )
  pstring <- paste("'p1' =",p[1]," and 'p2' =",p[2],"(probability of answering the sensitive question in group 1 and group 2 respectively)")
  res <- list(model="UQTunknown",call=paste0("Unrelated Question Technique with unknown prevalence of unrelated question.",pstring), pi=pi, piSE=piSE, piUQ=piUQ, piUQSE=piUQSE, n=n1+n2)
  return(res)
} 

########################################
# Mangat's model (Mangat, 1994)
RRuni.Mangat<-function(response,p){
  n <- length(response)
  lambda <- sum(response)/(n)
  pi <- (lambda-1+p)/p
  piSE <- sqrt(lambda*(1-lambda)/( (n-1)*p^2 )  )  
  res <- list(model="Mangat",call=paste0("Mangat's Model with p =",p, " (probability for noncarriers to respond with 'Yes'"),pi=pi,piSE=piSE,n=n)
  return(res)
}

########################################
# Kuk's model (Kuk, 1990) 
RRuni.Kuk<-function(response,p){
  # number of repetitions: how often the procedure was applied 
  rep <- max(response)
  n <- length(response)
  lambda <- sum(response)/(n*rep)
  pi <- (lambda-p[2])/(p[1]-p[2])
  phi <- pi*p[1]+(1-pi)*p[2]
  piSE <- sqrt( phi*(1-phi)/(rep*n*(p[1]-p[2])^2) + (1-1/rep)*pi*(1-pi)/n)  
  res <- list(model="Kuk",call=paste("Kuk's Model with ",rep," repetitions and p1 =",p[1],", p2 =",p[2]),pi=pi,piSE=piSE,n=n,Kukrep=rep)
  return(res)
} 

########################################
# Forced Response Model (formulas from Himmelfarb, 2008)
RRuni.FR<-function(response,p){
  numCat=length(p)
    
  # parameter for mixture distribution: proportion of true answers
  pTrue <- 1-sum(p);  
  n <- length(response)
  pi <- vector(length=numCat)
  piSE <- vector(length=numCat)
  lambda <- vector(length=numCat)    
  callText <- ""
  
  for (i in 1:numCat){
    lambda[i] <- table(factor(response, levels=0:(numCat-1)))[i]/n
    pi[i] <- (lambda[i]-p[i])/pTrue
    piSE[i] <- sqrt((lambda[i]*(1-lambda[i]))/(pTrue^2*n))
    callText <- paste(callText,"p",i-1,"=",p[i]," , ",sep="")
  }
  res <- list(model="FR",call=paste("Forced Response (FR) Model with ",numCat,
                                    " response categories and randomization probabilities ",callText),pi=pi,piSE=piSE,n=n)
  return(res)
} 



########################################
# stochastic lie detector (Moshagen, Musch, Erdfelder,2012)
RRuni.SLD<-function(response,p,group){
  group <- as.numeric(group)
  p1 <- p[1]
  p2 <- p[2]
  n1 <- length(response[group==1])
  n2 <- length(response[group==2])
  lambda1 <- sum(response[group==1])/n1
  lambda2 <- sum(response[group==2])/n2
  
  pi <- ((lambda2-lambda1)+(p2-p1))/(p2-p1)
  piSE <- sqrt(
    (lambda1*(1-lambda1))/(n1*(p1-p2)^2)+
      (lambda2*(1-lambda2))/(n2*(p2-p1)^2))
  t <- (lambda2*(1-p1)-lambda1*(1-p2))/((lambda2-lambda1)+(p2-p1))
  tSE <- sqrt(
    (lambda1*(1-lambda1))/n1*(((p1-p2)*(p2+lambda2-1))/(p1-p2+lambda1-lambda2)^2)^2+
      (lambda2*(1-lambda2))/n2*(((p1-p2)*(p1+lambda1-1))/(p1-p2+lambda1-lambda2)^2)^2  )
  
  res <- list(model="SLD",call=paste("Stochastic Lie Detector with p1 =",
                                     round(p1,4),", p2 =",round(p2,4)),
              pi=pi,piSE=piSE,t=t,tSE=tSE,n=c(n1,n2))
  return(res)
} 

########################################
# Cheater Detection Model (Clark, Desharnais, 1998)
RRuni.CDM<-function(response,p,group){
  group <- as.numeric(group)
  p1 <- p[1]
  p2 <- p[2]
  n1 <- sum(group==1)
  n2 <- sum(group==2)
  lambda1 <- sum(response[group==1])/n1
  lambda2 <- sum(response[group==2])/n2
  
  pi <- (p2*lambda1-p1*lambda2)/(p2-p1)
  piSE <- sqrt(( p1^2*lambda2*(1-lambda2)/n2+
                p2^2*lambda1*(1-lambda1)/n1)/(p1-p2)^2)
  beta <- (lambda2-lambda1)/(p2-p1)
  betaSE <- sqrt(
    (lambda2*(1-lambda2)/n2+lambda1*(1-lambda1)/n1)/(p1-p2)^2)
  gamma <- 1 + ((p2-1)*lambda1+(1-p1)*lambda2)/(p1-p2)
  gammaSE <- sqrt(((p2-1)^2*lambda1*(1-lambda1)/n1+(1-p1)^2*lambda2*(1-lambda2)/n2)/((p1-p2)^2))
  
  res=list(model="CDM",call=paste("Cheater Detection Model with p1 =",
                                  round(p1,4),", p2 =",round(p2,4)),
           pi=pi,piSE=piSE,beta=beta,betaSE=betaSE,gamma=gamma,gammaSE=gammaSE,n=c(n1,n2))
  return(res)
} 

########################################
# Cheater Detection Model symetric (Ostapczuk et al., 2009)
RRuni.CDMsym<-function(response,p,group){
  group <- as.numeric(group)
  n1 <- sum(group==1)
  n2 <- sum(group==2) 
  lambda1 <- sum(response[group==1])/n1
  lambda2 <- sum(response[group==2])/n2
  pi.denom <- (p[3]*(1-p[2])-p[1]*(1-p[4]))
  pi <- (p[3]*lambda1-p[1]*lambda2) / pi.denom
  gamma.denom <- p[1]*(1-p[4])-p[3]*(1-p[2])
  gamma <- ( (p[1]-lambda1)*(1-p[3]-p[4])-(p[3]-lambda2)*(1-p[1]-p[2])) / gamma.denom
  beta <- 1-pi-gamma
  l1.var <- lambda1*(1-lambda1)/n1
  l2.var <- lambda2*(1-lambda2)/n2
  piSE <- sqrt( p[3]^2 *l1.var +p[1]^2*l2.var ) / abs(pi.denom)
  gammaSE <- sqrt( (1-p[3]-p[4])^2*l1.var + (1-p[1]-p[2])^2*l2.var ) / abs(gamma.denom)
  res=list(model="CDMsym",call=paste0("Cheater Detection Model (symmetric); for group 1: p1 = ",p[1]," , p2 = ",p[2],"; for group 2: p3 = ",p[3],", p4 = ",p[4]," (probabilities of directed 'Yes'/'No' answers)"),
           pi=pi,piSE=piSE,beta=beta, betaSE=NA, gamma=gamma,gammaSE=gammaSE,n=c(n1,n2))
  return(res)
} 

RRuni.CDMsym.ll <- function(param,x,p,group){
  pi <- param[1]
  gg <-param[2]   # greek gamma =: gg
  s <- group==1   # group selection vector
  vec1 <- x[s]*log(pi*(1-p[2])+(1-pi-gg)*p[1])+(1-x[s])*log(pi*p[2]+(1-pi-gg)*(1-p[1])+gg)
  vec2 <- x[!s]*log(pi*(1-p[4])+(1-pi-gg)*p[3])+(1-x[!s])*log(pi*p[4]+(1-pi-gg)*(1-p[3])+gg)
  ll <- sum(vec1)+sum(vec2)
  return(ll)
}

RRuni.CDMsym.llgrad <- function(param,x,p,group){
  pi <- param[1]
  gg <-param[2]   # greek gamma =: gg
  grad = rep(0,2)
  s <- group==1   # group selection vector
  # pi
  vec1 <- (x[s]*(1-p[2]-p[1])) / (pi*(1-p[2])+(1-pi-gg)*p[1])  +
            ( (1-x[s])* (p[2]-1+p[1]) ) / (pi*p[2]+(1-pi-gg)*(1-p[1])+gg)
  vec2 <- (x[!s]*(1-p[4]-p[3])) / (pi*(1-p[4])+(1-pi-gg)*p[3])  +
            ( (1-x[!s])*(p[4]-1+p[3]) ) / (pi*p[4]+(1-pi-gg)*(1-p[3])+gg)
  grad[1] <- sum(vec1)+sum(vec2)
  # gamma
  vec1 <- - (x[s]*p[1]) / (pi*(1-p[2])+(1-pi-gg)*p[1]) +
            ((1-x[s])*p[1]) / (pi*p[2]+(1-pi-gg)*(1-p[1])+gg)
  vec2 <- - (x[!s]*p[3]) / (pi*(1-p[4])+(1-pi-gg)*p[3]) +
            ((1-x[!s])*p[3]) / (pi*p[4]+(1-pi-gg)*(1-p[3])+gg)
  grad[2] <- sum(vec1)+sum(vec2)
  return(grad)
}
# pp <- c(.2,0,.8,0)
# par <- c(0.3017339, 0.1960341)
# RRuni.CDMsym.ll(c(.3,0.2),datS$response,pp,datS$group)
# RRuni.CDMsym(datS$response,pp,datS$group)
# print(RRuni.CDMsym.llgrad(par,datS$response,pp,datS$group))
# print(grad(func=RRuni.CDMsym.ll,x=par,xxx=datS$response,p=pp,group=datS$group))

########################################
# Continuous RR: Greenberg kuebler abernathy horvitz 1971 (p. 247)
RRuni.mix.norm<-function(response,p){
  mean.obs <- mean(response)
  pi <- (mean.obs - (1-pt) * p[2]) / pt
  n <- length(response)
  piSE <- sd(response)/(sqrt(n)*p)
  pstring <- paste0("'p[1]'=",round(p[1],4)," (probability of answering sensitive question), 'p[2]'=",p[2]," (p[3]=",p[3],"): mean (SD) of masking distribution")
  res <- list(model="mix.norm",call=paste0("Continuous mixture RR design with probability of truthful responding ",pstring),pi=pi,piSE=piSE,n=n)
  return(res)
} 

# Continuous RR
RRuni.mix.exp<-function(response,p){
  mean.obs <- mean(response)
  pi <- (mean.obs - (1-pt) * p[2]) / pt
  n <- length(response)
  piSE <- sd(response)/(sqrt(n)*p)  
  pstring <- paste0("'p[1]'=",round(p[1],4)," (probability of answering sensitive question), 'p[2]'=",p[2],": mean and SD of exponential masking distribution")
  res <- list(model="mix.norm",call=paste0("Continuous mixture RR design with probability of truthful responding ",pstring),pi=pi,piSE=piSE,n=n)
  return(res)
} 

# two questions with two unknown distributions (two-group model)
RRuni.mix.unknown<-function(response,p,group){
  mean.obs1 <- mean(response[group==1])
  mean.obs2 <- mean(response[group==2])
  pi <- ((1-p[2])*mean.obs1 - (1-p[1])*mean.obs2) /(p[1]-p[2])
  piUQ <- (p[2]*mean.obs1 - p[1]*mean.obs2)/(p[2]-p[1])
  
  v.obs1 <- var(response[group==1])/sum(group==1)
  v.obs2 <- var(response[group==2])/sum(group==2)
  piSE <- sqrt((1-p[2])^2*v.obs1 + (1-p[1])^2*v.obs2)/abs(p[1]-p[2])
  piUQSE <- sqrt(p[2]^2*v.obs1 + p[1]^2*v.obs2)/abs(p[2]-p[1])
  
  # for RRcor: calculate sigma^2(sensitive) ... not SE! (seems to work fine)
  num <- p[2]*(p[1]*pi+(1-p[1])*piUQ)^2-p[1]*(p[2]*pi+(1-p[2]*piUQ))^2
  y.var <- var(response[group==1]) - piUQ^2 + num/(p[2]-p[1])
  num <- (1-p[2])*(p[1]*pi+(1-p[1])*piUQ)^2-(1-p[1])*(p[2]*pi+(1-p[2]*piUQ))^2
  x.var<- var(response[group==2]) - pi^2 + num/(p[1]-p[2])
  
  pstring <- paste0("'p[1]='",round(p[1],4)," and p[2]=",round(p[2],4),": probability of responding to sensitive question in group 1 and 2, respectively")
  res <- list(model="mix.unknown",call=paste0("Continuous mixture RR design with probability of truthful responding ",pstring),pi=pi,piSE=piSE,piUQ=piUQ, piUQSE=piUQSE,n=length(response),y.var=y.var, x.var=x.var)
#   print(x.var)
#   print(y.var)
  return(res)
} 


RRuni.custom <- function(response, p, group){
  freq <- table(factor(response, levels=0:(ncol(p)-1)))
  n <- sum(freq)
  lambda <- freq/n
  Pinv <- solve(p)
  pi <- c(Pinv %*% lambda)
  # van den hout, van der heijden 2002
  piCov <- (Pinv %*% ( diag(lambda) - lambda %*% t(lambda))%*% t(Pinv) )/(n-1)
  piSE <- sqrt(diag(piCov))
  res <- list(model="custom",call=paste0("RR as misclassification"),
              pi=pi,piSE=piSE, n=n) #,y.var=y.var, x.var=x.var)
  return(res)
}