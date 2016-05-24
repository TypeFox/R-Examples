library("mvtnorm")

#' Hypothesis: bootstrapping data causes estimated pcor network to be more strongly connected

#' Dependent variable: average absolute edge weight
#' Replicate 100 times: sample independent data, estimate pcor, bootstrap, estimate pcor

nVar <- 20
nObs <- 100
nRep <- 100

Res <- do.call(rbind,replicate(nRep,{
  Data <- matrix(rnorm(nVar*nObs),nObs,nVar)
  net <- -cov2cor(solve(cov(Data)))
  diag(net) <- 0
  samp <- mean(abs(net[upper.tri(net,diag=FALSE)]))
  
  # nonparametric bootstrap:
  newData <- Data[sample(1:nrow(Data),nrow(Data),TRUE),]    
  net2 <- -cov2cor(solve(cov(newData)))
  nonPar <- mean(abs(net2[upper.tri(net2,diag=FALSE)]))
  
  # parametrics bootstrap:
  newData <- rmvnorm(nrow(Data),sigma = solve(diag(1,nrow(net)) - net))
  net3 <- -cov2cor(solve(cov(newData)))
  Par <- mean(abs(net3[upper.tri(net3,diag=FALSE)]))
  
  return(data.frame(sample=samp,nonparametric=nonPar,parametric=Par))
},simplify=FALSE))


#' Non-parametric:
t.test(Res$sample,Res$parametric,paired=FALSE)


#' Parametric:
t.test(Res$sample,Res$parametric,paired=TRUE)
