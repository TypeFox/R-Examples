#### generazione dati discreti
###                         el 06/02/2013
####################################################
#rm(list=ls())
#library(polycor)
#require(MASS)
## default
#n <- 10
#Q <- NULL
#R <- matrix(c(1,.6,.6,1),2,2)
#th <- NULL
#probs <- NULL

#th <- list(c(-Inf,qchisq(pbinom(0:(Q[1]-2),(length(1:Q[1])-1),.5),1),Inf),
#           c(-Inf,qnorm(pbinom(0:(Q[1]-2),(length(1:Q[1])-1),.5)),Inf))

#th <- list(c(-Inf,qchisq(pbinom(0:3,4,.5),1),Inf),
#           c(-Inf,qnorm(pbinom(0:2,3,.5)),Inf))

#for (j in 1:length(th)) {
#  l <- length(th[[j]])
#  probs <- c(probs,list(pnorm(th[[j]])[2:l]-pnorm(th[[j]])[1:(l-1)]))
#}

#probs <- list(c(.125,.375,.375,.125),c(.125,.375,.375,.125))
#th <- NULL

## se le soglie o le probs non sono nulle deve controllare che Q sia coerente.
#n <- 10
#R <- matrix(c(1,.4,.4,1),2,2)
#Q <- th <- NULL
#probs <- list(c(.125,.375,.375,.125),c(.125,.375,.375,.125))
#Dx <- rdatagen(n=10000,R=R,probs=probs)$data

### DA FARE  29/10/2014
### gestione liste nel caso in cui siano tutte uguali
# rdatagen(10,probs=c(.2,.2,.2,.2,.2),Q=5)
# Q dovrebbe essere lungo 1 oppure ncol(R)
# probs deve essere una lista ed essere probabilita se non lo e si normalzza con warning
#n <- 10; R <- diag(1,2); Q <- 4; probs <- list(c(.2,.2,.2),c(0,.5,.5)); th <- NULL
#probs <- c(2,5,0,3)
#rdatagen(10,Q=4,probs=probs)

#########################################################
rdatagen <- function(n=100,R=diag(1,2),Q=NULL,th=NULL,probs=NULL) {
  
  ## controlli probs  aggiunto 29/10/2014
  if (!is.null(probs)) {
    if (!is.list(probs)) {
      warning("probs set to list") ## inglese da controllare
      probs <- list(probs)
    }
    for (j in 1:length(probs)) {
      if (sum(probs[[j]])!=1) {
        warning("probs normalized") ## inglese da controllare
        probs[[j]] <- probs[[j]]/sum(probs[[j]])
      }
    }
  }
  
  if ((is.null(Q))&(is.null(th))&(is.null(probs))) {
    Q <- rep(1,ncol(R))
  } else {
    if ((is.null(Q))&(!is.null(th))) for (j in 1:length(th)) Q <- c(Q,length(th[[j]])-1)
    if (is.null(Q)&(is.null(th))&(!is.null(probs))) for (j in 1:length(probs)) Q <- c(Q,length(probs[[j]]))
  }
  #Q
  
  if (length(Q)==1) { # Q dovrebbe essere lungo 1 oppure ncol(R)
    Q <- rep(Q,ncol(R))
  }
  #Q
  
  if ((is.null(th))&(is.null(probs))) {
    for (j in 1:ncol(R)) {
      if (Q[j]>1) {
        th <- c(th,list(c(-Inf,qnorm(pbinom(0:(Q[j]-2),(length(1:Q[j])-1),.5)),Inf)))    
      } else {
        th <- c(th,list(c(-Inf,Inf)))
      }
    }
  }
  #th
  
  if ((is.null(probs))&(!(is.null(th)))) {
    for (j in 1:length(th)) {
      if (length(th[[j]])!=(Q[j]+1)) {
        cat(paste("found thresholds:", paste(th[[j]],collapse=","),"and Q:",Q[j]),"\n")
        stop("number of thresholds must be equal to Q+1") 
      }
      l <- length(th[[j]])
      probs <- c(probs,list(pnorm(th[[j]])[2:l]-pnorm(th[[j]])[1:(l-1)]))
    }
  }
  probs
  
  if ((is.null(th))&(!is.null(probs))) {
    
    ## se le probs sono tutte uguali per ciascun item 
    if (length(probs)==1) {
      PROBS <- probs
      for (j in 1:(ncol(R)-1)) probs <- c(probs,PROBS)
    }
    probs
    
    for (j in 1:length(probs)) {
      if (length(probs[[j]])!=Q[j]) {
        cat(paste("found probs:", paste(probs[[j]],collapse=","),"and Q:",Q[j]),"\n")
        stop("number of probabilities must be equal to Q")
      }
      th <- c(th,list(c(-Inf,qnorm(cumsum(probs[[j]])))))
    }
  }

  ## dati continui
  X <- mvrnorm(n,rep(0,ncol(R)),R)

  ### discretizzazione
  Dx <- X
  colonne <- which(Q>1)
  #for (j in colonne) Dx[,j] <- cut(Dx[,j],th[[j]])
  for (j in colonne) Dx[,j] <-.bincode(Dx[, j], th[[j]]) ### 
  
  return(list(data=Dx,R=R,thresholds=th,probs=probs))
}

#X <- rdatagen(1000,Q=5,probs=c(.1,8,0,0,5))
#par(mfrow=c(1,2))
#for (j in 1:ncol(X$data)) barplot(table(X$data[,j]))


#Q <- c(5,0)
#Q <- 2
#th <- list(c(-Inf,qchisq(pbinom(0:(Q-2),(length(1:Q)-1),.5),1),Inf),
#           c(-Inf,Inf))
           #c(-Inf,qnorm(pbinom(0:(Q-2),(length(1:Q)-1),.5)),Inf))
#probs <- list(c(.125,.375,.375,.125),c(.125,.375,.375,.125))
#X <- simdata(n=1000,Q=Q,R=R) #,probs=probs)
#X <- simdata()
#X <- simdata(n=10000,R=R,Q=c(1,3))
#X <- simdata(n=10000,R=R,Q=c(2,3),th=th)
#X <- simdata(n=10000,R=R,probs=probs)

#(th <- X$thresholds)
#(probs <- X$probs)
#Dx <- X$data

#summary(Dx)
#Dx <- as.data.frame(Dx)

#par(mfrow=c(1,2))
#for (j in 1:ncol(Dx)) {
#  if (length(probs[[j]])>1) {
#    barplot(table(Dx[,j])/sum(table(Dx[,j])))  
#    Dx[,j] <- ordered(Dx[,j])
#  } else {
#    hist(Dx[,j],main="")
#  }
#}
#hetcor(Dx)$correlations 

#round(probs[[1]],4)
#round(probs[[2]],4)
