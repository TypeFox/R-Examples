.AddDiag <-
function(mat,Vec){
  return(mat+diag(Vec))
}
.ClassIBSUnphased <-
function(ObsSnp){
  Obs <- paste(as.character(ObsSnp[1]),as.character(ObsSnp[2]))
  
  res <- switch(Obs,
                "0 0" = 0,
                "0 1" = 1,
                "0 2" = 2,
                "1 0" = 3,
                "1 1" = 4,
                "1 2" = 5,
                "2 0" = 6,
                "2 1" = 7,
                "2 2" = 8)
  return(res)
}
.ConditionIdentifiabilityCouple <-
function(Crossing,ParentPop){
if ((length(unique(c(Crossing)))<4)||(length(unique(ParentPop))>2)){
return(1)
}else{
return(NA)
}
}
.EM <-
function(ProbCond , DeltaInit , Prec){
ProbCondVec <- c(t(ProbCond))
.C("BoucleEMacc" , as.double(ProbCondVec) , as.integer(dim(ProbCond)[2]) , as.integer(dim(ProbCond)[1]) , as.double(Prec) , as.double(DeltaInit) , as.double(rep(0,dim(ProbCond)[1]*dim(ProbCond)[2])))
}
.LogVrai <-
function(CondVrais){
  res <- sum(log(rowSums(CondVrais)))
  return(res)
}
.MatTriSup <-
function(mat,Vecteur){
  mat[lower.tri(mat,diag=F)] <- Vecteur
  return(t(mat))
}
.RelatednessCouple <-
function(HybridGenom,Freq,Crossing,LinePop,Phased=FALSE,NbInit=5,Prec=10^(-4)){
NbSnp <- dim(HybridGenom)[1]
LinePopChar <- paste((1:4)[factor(LinePop,levels=unique(LinePop))],collapse=" ")
VecCrossChar <- paste((1:4)[factor(c(Crossing),levels=unique(c(Crossing)))],collapse=" ")


  
FreqPop <- Freq[,LinePop]
if (Phased==TRUE){
ImpossibleRelatednessCrossing <- switch(VecCrossChar,
                                          "1 1 1 1" = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),
                                          "1 1 1 2" = c(1,2,3,4,5,6,7,8,9,10,12,13,14),
                                          "1 1 2 1" = c(1,2,3,4,5,6,7,8,9,10,11,13,14),
                                          "1 1 2 2" = c(1,2,3,4,5,6,7,9,10,11,12,13,14),
                                          "1 1 2 3" = c(1,3,4,5,6,7,9,10,13,14),
                                          "1 2 1 1" = c(1,2,3,4,5,6,7,8,9,10,11,12,14),
                                          "1 2 1 2" = c(1,2,3,4,5,6,7,8,10,11,12,13,14),
                                          "1 2 1 3" = c(1,2,3,5,6,7,8,10,12,14),
                                          "1 2 2 1" = c(1,2,3,4,5,6,7,8,9,11,12,13,14),
                                          "1 2 2 2" = c(1,2,3,4,5,6,7,8,9,10,11,12,13),
                                          "1 2 2 3" = c(1,2,3,4,5,7,8,9,12,13),
                                          "1 2 3 1" = c(1,2,3,4,6,7,8,9,11,14),
                                          "1 2 3 2" = c(1,2,3,4,5,6,8,9,12,13),
                                          "1 2 3 3" = c(1,2,4,5,6,7,9,10,11,12),
                                          "1 2 3 4" = c())
  
ImpossibleRelatednessPop <- switch(LinePopChar,
                                     "1 1 1 1" = c(),
                                     "1 1 1 2" = c(3,5,7,8,9,10,12,13,14,15),
                                     "1 1 2 1" = c(3,4,6,8,9,10,11,13,14,15),
                                     "1 1 2 2" = c(4,5,6,7,9,10,11,12,13,14,15),
                                     "1 1 2 3" = c(3,4,5,6,7,8,9,10,11,12,13,14,15),
                                     "1 2 1 1" = c(2,6,7,8,9,10,11,12,14,15),
                                     "1 2 1 2" = c(2,3,5,6,8,10,11,12,13,14,15),
                                     "1 2 1 3" = c(2,3,5,6,7,8,9,10,11,12,13,14,15),
                                     "1 2 2 1" = c(2,3,4,7,8,9,11,12,13,14,15),
                                     "1 2 2 2" = c(2,4,5,8,9,10,11,12,13,15),
                                     "1 2 2 3" = c(2,3,4,5,7,8,9,10,11,12,13,14,15),
                                     "1 2 3 1" = c(2,3,4,6,7,8,9,10,11,12,13,14,15),
                                     "1 2 3 2" = c(2,3,4,5,6,8,9,10,11,12,13,14,15),
                                     "1 2 3 3" = c(2,4,5,6,7,8,9,10,11,12,13,14,15),
                                     "1 2 3 4" = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15))
  
ImpossibleRelatedness <- c(ImpossibleRelatednessCrossing,ImpossibleRelatednessPop)


GenoChar <- paste(HybridGenom[,1],HybridGenom[,2],HybridGenom[,3],HybridGenom[,4],sep='')
IBS <- strtoi(GenoChar, base = 2)
ListProbCond <- lapply(0:15 , function(x) .VecProbCondPhased(FreqPop,x,(IBS==x)))

DeltaInit <- matrix(runif(15*NbInit,min=0,max=1),NbInit,15)
DeltaInit[,ImpossibleRelatedness] <- 0
DeltaInit <- DeltaInit/rowSums(DeltaInit)
NbConfIBD <- 15
Class <- c(15,8,11,12,2,13,14,3,9,10,4,5,6,7,1)
}else{


  
ImpossibleRelatednessCrossing <- switch(VecCrossChar,
                                          "1 1 1 1" = c(1,2,3,4,5,6,7,8),
                                          "1 1 1 2" = c(1,2,3,4,5,6,8),
                                          "1 1 2 1" = c(1,2,3,4,5,6,8),
                                          "1 1 2 2" = c(1,2,3,4,6,7,8),
                                          "1 1 2 3" = c(1,3,4,6,8),
                                          "1 2 1 1" = c(1,2,3,4,5,6,7),
                                          "1 2 1 2" = c(1,2,3,4,5,7,8),
                                          "1 2 1 3" = c(1,2,3,5),
                                          "1 2 2 1" = c(1,2,3,4,5,7,8),
                                          "1 2 2 2" = c(1,2,3,4,5,6,7),
                                          "1 2 2 3" = c(1,2,3,5),
                                          "1 2 3 1" = c(1,2,3,5),
                                          "1 2 3 2" = c(1,2,3,5),
                                          "1 2 3 3" = c(1,2,4,8),
                                          "1 2 3 4" = c())
  
ImpossibleRelatednessPop <- switch(LinePopChar,
                                     "1 1 1 1" = c(),
                                     "1 1 1 2" = c(3,5,6,8,9),
                                     "1 1 2 1" = c(3,5,6,8,9),
                                     "1 1 2 2" = c(4,6,7,8,9),
                                     "1 1 2 3" = c(3,4,5,6,7,8,9),
                                     "1 2 1 1" = c(2,5,6,7,9),
                                     "1 2 1 2" = c(2,3,5,7,8,9),
                                     "1 2 1 3" = c(2,3,5,6,7,8,9),
                                     "1 2 2 1" = c(2,3,5,7,8,9),
                                     "1 2 2 2" = c(2,5,6,7,9),
                                     "1 2 2 3" = c(2,3,5,6,7,8,9),
                                     "1 2 3 1" = c(2,3,5,6,7,8,9),
                                     "1 2 3 2" = c(2,3,5,6,7,8,9),
                                     "1 2 3 3" = c(2,4,5,6,7,8,9),
                                     "1 2 3 4" = c(2,3,4,5,6,7,8,9))
  
ImpossibleRelatedness <- c(ImpossibleRelatednessCrossing,ImpossibleRelatednessPop)


HybGen <- matrix(NA,NbSnp,2)
HybGen[,1] <- rowSums(HybridGenom[,c(1,2)])
HybGen[,2] <- rowSums(HybridGenom[,c(3,4)])
IBS <- sapply(1:NbSnp , function(x) .ClassIBSUnphased(HybGen[x,]))
ListProbCond <- lapply(0:8 , function(x) .VecProbCondUnphased(FreqPop,x,(IBS==x),LinePop))

DeltaInit <- matrix(runif(9*NbInit,min=0,max=1),NbInit,9)
DeltaInit[,ImpossibleRelatedness] <- 0
DeltaInit <- DeltaInit/rowSums(DeltaInit)
NbConfIBD <- 9
Class <- c(9,5,7,2,8,3,6,4,1)
}
ProbCond <- Reduce("+",ListProbCond)

List <- sapply(1:NbInit, function(x) .EM(ProbCond,DeltaInit[x,],Prec))
LogV <- sapply(1:NbInit, function(x) .LogVrai(matrix(List[6,][[x]],NbSnp,NbConfIBD,byrow=TRUE)))
DeltaEst <- List[5,][[which.max(LogV)]][Class]
LogVEst <- max(LogV)

return(list(Delta=DeltaEst,LogV=LogVEst))
}
.RelatednessLineCouple <-
function(Genotype,Freq,LinePop,NbInit=5,Prec=10^(-4)){
NbSnp <- dim(Genotype)[1]
GenoChar <- paste(Genotype[,1],Genotype[,2],sep="")
IBS <- strtoi(GenoChar,base=2)
ListProbCond <- lapply(0:3 , function(x) .VecProbCondLine(Freq,x,(IBS==x)))
DeltaInit <- matrix(runif(2*NbInit,0,1),NbInit,2)
DeltaInit <- DeltaInit/rowSums(DeltaInit)
ProbCond <- Reduce("+",ListProbCond)
List <- sapply(1:NbInit, function(x) .EM(ProbCond,DeltaInit[x,],Prec))
LogV <- sapply(1:NbInit, function(x) .LogVrai(matrix(List[6,][[x]],NbSnp,2,byrow=TRUE)))
DeltaEst <- List[5,][[which.max(LogV)]]
LogVEst <- max(LogV)
return(list(Delta=DeltaEst,LogV=LogVEst))
}
.TheFifteenDeltaGraph <-
function(){
dev.new()
layout(matrix(seq(1,15),byrow=TRUE,3,5),1,1)
plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta1',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(0,1),col='red')
lines(c(0,0),c(0,1),col='red')
lines(c(0,1),c(1,1),col='red')
lines(c(0,1),c(0,0),col='red')
lines(c(0,1),c(1,0),col='red')
lines(c(1,1),c(0,1),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta2',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(1,1),col='red')
lines(c(0,1),c(0,0),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta3',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(0,1),col='red')
lines(c(0,0),c(0,1),col='red')
lines(c(0,1),c(1,1),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta4',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(1,0),col='red')
lines(c(1,1),c(0,1),col='red')
lines(c(0,1),c(1,1),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta5',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(1,1),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta6',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(1,0),col='red')
lines(c(0,0),c(0,1),col='red')
lines(c(0,1),c(0,0),col='red')
  
plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta7',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(0,1),col='red')
lines(c(1,1),c(0,1),col='red')
lines(c(0,1),c(0,0),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta8',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(0,0),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta9',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,0),c(0,1),col='red')
lines(c(1,1),c(0,1),col='red')
  
plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta10',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(0,1),col='red')
lines(c(0,1),c(1,0),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta11',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,0),c(0,1),col='red')
  
plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta12',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(1,0),col='red')
  
plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta13',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(0,1),col='red')
  
plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta14',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(1,1),c(0,1),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta15',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
}
.TheNineDeltaGraph <-
function(){  
dev.new()
layout(matrix(seq(1,9),byrow=TRUE,3,3),1,1)
plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta1',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(0,0),col='red')
lines(c(0,1),c(1,1),col='red')
lines(c(0,0),c(0,1),col='red')
lines(c(1,1),c(0,1),col='red')
lines(c(0,1),c(0,1),col='red')
lines(c(0,1),c(1,0),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta2',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(0,0),col='red')
lines(c(0,1),c(1,1),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta3',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(1,1),col='red')
lines(c(0,1),c(0,1),col='red')
lines(c(0,0),c(0,1),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta4',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(1,1),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta5',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(0,0),col='red')
lines(c(0,1),c(1,0),col='red')
lines(c(0,0),c(0,1),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta6',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,1),c(0,0),col='red')


plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta7',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,0),c(0,1),col='red')
lines(c(1,1),c(0,1),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta8',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
lines(c(0,0),c(0,1),col='red')

plot(c(0,1,0,1),c(0,0,1,1),type='p',main='Delta9',xlab='Gamete',ylab='Hybrids',yaxt='n',xaxt='n')
}
.VecProbCondLine <-
function(Freq,IBS,Position){
FreqIBS <- Freq*Position
q <- Position-FreqIBS
p <- FreqIBS
IBS <- as.character(IBS)
res <- switch(IBS,
"0" = cbind(q^2,q),
"1" = cbind(q*p,0),
"2" = cbind(q*p,0),
"3" = cbind(p^2,p))
return(res)
}
.VecProbCondPhased <-
function(Freq , IBS , Position){
FreqIBS <- Freq*Position
q1 <- Position-FreqIBS[,1]
p1 <- FreqIBS[,1]
q2 <- Position-FreqIBS[,2]
p2 <- FreqIBS[,2]
q3 <- Position-FreqIBS[,3]
p3 <- FreqIBS[,3]
q4 <- Position-FreqIBS[,4]
p4 <- FreqIBS[,4]
IBS <- as.character(IBS)
  
res <- switch(IBS,
                "0" = cbind(q1*q2*q3*q4 , q1*q3*q4 , q1*q2*q3 , q1*q2*q4 , q1*q2*q3 , q1*q2*q4 , q1*q2*q3 , q1*q3 , q1*q2 , q1*q2 , q1*q4 , q1*q3 , q1*q2 , q1*q2 , q1),
                "1" = cbind(q1*q2*q3*p4 , q1*q3*p4 , 0 , q1*q2*p4 , 0 , q1*q2*p4 , 0 , 0 , 0 , 0 , q1*p4 , 0 , 0 , 0 , 0),
                "2" = cbind(q1*q2*p3*q4 , q1*p3*q4 , 0 , 0 , q1*q2*p3 , 0 , q1*q2*p3 , 0 , 0 , 0 , 0 , q1*p3 , 0 , 0 , 0),
                "3" = cbind(q1*q2*p3*p4 , q1*p3*p4 , q1*q2*p3 , 0 , 0 , 0 , 0 , q1*p3 , 0 , 0 , 0 , 0 , 0 , 0 , 0),
                "4" = cbind(q1*p2*q3*q4 , 0 , q1*p2*q3 , q1*p2*q4 , q1*p2*q3 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , q1*p2 , 0 , 0),
                "5" = cbind(q1*p2*q3*p4 , 0 , 0 , q1*p2*p4 , 0 , 0 , q1*p2*q3 , 0 , q1*p2 , 0 , 0 , 0 , 0 , 0 , 0),
                "6" = cbind(q1*p2*p3*q4 , 0 , 0 , 0 , q1*p2*p3 , q1*p2*q4 , 0 , 0 , 0 , q1*p2 , 0 , 0 , 0 , 0 , 0),
                "7" = cbind(q1*p2*p3*p4 , 0 , q1*p2*p3 , 0 , 0 , q1*p2*p4 , q1*p2*p3 , 0 , 0 , 0 , 0 , 0 , 0 , q1*p2 , 0),
                "8" = cbind(p1*q2*q3*q4 , 0 , p1*q2*q3 , 0 , 0 , p1*q2*q4 , p1*q2*q3 , 0 , 0 , 0 , 0 , 0 , 0 , p1*q2 , 0),
                "9" = cbind(p1*q2*q3*p4 , 0 , 0 , 0 , p1*q2*q3 , p1*q2*p4 , 0 , 0 , 0 , p1*q2 , 0 , 0 , 0 , 0 , 0),
                "10" = cbind(p1*q2*p3*q4 , 0 , 0 , p1*q2*q4 , 0 , 0 , p1*q2*p3 , 0 , p1*q2 , 0 , 0 , 0 , 0 , 0 , 0),
                "11" = cbind(p1*q2*p3*p4 , 0 , p1*q2*p3 , p1*q2*p4 , p1*q2*p3 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , p1*q2 , 0 , 0),
                "12" = cbind(p1*p2*q3*q4 , p1*q3*q4 , p1*p2*q3 , 0 , 0 , 0 , 0 , p1*q3 , 0 , 0 , 0 , 0 , 0 , 0 , 0),
                "13" = cbind(p1*p2*q3*p4 , p1*q3*p4 , 0 , 0 , p1*p2*q3 , 0 , p1*p2*q3 , 0 , 0 , 0 , 0 , p1*q3 , 0 , 0 , 0),
                "14" = cbind(p1*p2*p3*q4 , p1*p3*q4 , 0 , p1*p2*q4 , 0 , p1*p2*q4 , 0 , 0 , 0 , 0 , p1*q4 , 0 , 0 , 0 , 0),
                "15" = cbind(p1*p2*p3*p4 , p1*p3*p4 , p1*p2*p3 , p1*p2*p4 , p1*p2*p3 , p1*p2*p4 , p1*p2*p3 , p1*p3 , p1*p2 , p1*p2 , p1*p4 , p1*p3 , p1*p2 , p1*p2 , p1))
 
return(res)
}
.VecProbCondUnphased <-
function(Freq , IBS , Position , LinePop){
FreqIBS <- Freq*Position
q1 <- Position-FreqIBS[,1]
p1 <- FreqIBS[,1]
q2 <- Position-FreqIBS[,2]
p2 <- FreqIBS[,2]
q3 <- Position-FreqIBS[,3]
p3 <- FreqIBS[,3]
q4 <- Position-FreqIBS[,4]
p4 <- FreqIBS[,4]
IBS <- as.character(IBS)
  
res <- switch(IBS,
"0" = cbind(q1*q2*q3*q4 , q1*q3*q4 , q1*q2*q3 , q1*q2*q4*(LinePop[1]==LinePop[3]) + q1*q2*q3*(LinePop[1]==LinePop[4]) + q1*q2*q4*(LinePop[2]==LinePop[3]) + q1*q2*q3*(LinePop[2]==LinePop[4]) , q1*q3 , q1*q2*((LinePop[1]==LinePop[3])*(LinePop[2]==LinePop[4]) + (LinePop[1]==LinePop[4])*(LinePop[2]==LinePop[3])) , q1*q4*(LinePop[1]==LinePop[3]) + q1*q3*(LinePop[1]==LinePop[4]) , q1*q2*(LinePop[2]==LinePop[3]) + q1*q2*(LinePop[1]==LinePop[3]) , q1),
                "1" = cbind(q1*q2*q3*p4 + q1*q2*p3*q4 , q1*q3*p4 + q1*p3*q4 , 0 , q1*q2*p4*(LinePop[1]==LinePop[3]) + q1*q2*p3*(LinePop[1]==LinePop[4]) + q1*q2*p4*(LinePop[2]==LinePop[3]) + q1*q2*p3*(LinePop[2]==LinePop[4]) , 0 , 0 , q1*p4*(LinePop[1]==LinePop[3]) + q1*p3*(LinePop[1]==LinePop[4]) , 0 , 0),
                "2" = cbind(q1*q2*p3*p4 , q1*p3*p4 , q1*q2*p3 , 0 , q1*p3 , 0 , 0 , 0 , 0),
                "3" = cbind(p1*q2*q3*q4 + q1*p2*q3*q4 , 0 , q1*p2*q3 + p1*q2*q3 , q1*p2*q4*(LinePop[1]==LinePop[3]) + q1*p2*q3*(LinePop[1]==LinePop[4]) + p1*q2*q4*(LinePop[2]==LinePop[3]) + p1*q2*q3*(LinePop[2]==LinePop[4]) , 0 , 0 , 0 , q1*p2*(LinePop[1]==LinePop[3]) + p1*q2*(LinePop[2]==LinePop[3]) , 0),
                "4" = cbind(q1*p2*q3*p4 + q1*p2*p3*q4 + p1*q2*q3*p4 + p1*q2*p3*q4 , 0 , 0 , (q1*p2*p4 + p1*q2*q4)*(LinePop[1]==LinePop[3]) + (q1*p2*p3 + p1*q2*q3)*(LinePop[1]==LinePop[4]) + (q1*p2*q4 + p1*q2*p4)*(LinePop[2]==LinePop[3]) + (q1*p2*q3 + p1*q2*p3)*(LinePop[2]==LinePop[4]) , 0 , (q1*p2 + p1*q2)*(LinePop[1]==LinePop[3])*(LinePop[2]==LinePop[4]) + (q1*p2 + p1*q2)*(LinePop[1]==LinePop[4])*(LinePop[2]==LinePop[3]) , 0 , 0 , 0),
                "5" = cbind(q1*p2*p3*p4 + p1*q2*p3*p4 , 0 , q1*p2*p3 + p1*q2*p3 , p1*q2*p4*(LinePop[1]==LinePop[3]) + p1*q2*p3*(LinePop[1]==LinePop[4]) + q1*p2*p4*(LinePop[2]==LinePop[3]) + q1*p2*p3*(LinePop[2]==LinePop[4]) , 0 , 0 , 0 , p1*q2*(LinePop[1]==LinePop[3]) + q1*p2*(LinePop[2]==LinePop[3]) , 0),
                "6" = cbind(p1*p2*q3*q4 , p1*q3*q4 , p1*p2*q3 , 0 , p1*q3 , 0 , 0 , 0 , 0),
                "7" = cbind(p1*p2*q3*p4 + p1*p2*p3*q4 , p1*q3*p4 + p1*p3*q4 , 0 , p1*p2*q4*(LinePop[1]==LinePop[3]) + p1*p2*q3*(LinePop[1]==LinePop[4]) + p1*p2*q4*(LinePop[2]==LinePop[3]) + p1*p2*q3*(LinePop[2]==LinePop[4]) , 0 , 0 , p1*q4*(LinePop[1]==LinePop[3]) + p1*q3*(LinePop[1]==LinePop[4]) , 0 , 0),
                "8" = cbind(p1*p2*p3*p4 , p1*p3*p4 , p1*p2*p3 , p1*p2*p4*(LinePop[1]==LinePop[3]) + p1*p2*p3*(LinePop[1]==LinePop[4]) + p1*p2*p4*(LinePop[2]==LinePop[3]) + p1*p2*p3*(LinePop[2]==LinePop[4]) , p1*p3 , p1*p2*((LinePop[1]==LinePop[3])*(LinePop[2]==LinePop[4]) + (LinePop[1]==LinePop[4])*(LinePop[2]==LinePop[3])) , p1*p4*(LinePop[1]==LinePop[3]) + p1*p3*(LinePop[1]==LinePop[4]) , p1*p2*(LinePop[2]==LinePop[3]) + p1*p2*(LinePop[1]==LinePop[3]) , p1))


res[,4] <- res[,4]/max((LinePop[1]==LinePop[3])+(LinePop[1]==LinePop[4])+(LinePop[2]==LinePop[3])+(LinePop[2]==LinePop[4]),1)
res[,6] <- res[,6]/max((LinePop[1]==LinePop[3])*(LinePop[2]==LinePop[4])+(LinePop[2]==LinePop[3])*(LinePop[1]==LinePop[4]),1)
res[,7] <- res[,7]/max((LinePop[1]==LinePop[3])+(LinePop[1]==LinePop[4]),1)
res[,8] <- res[,8]/max((LinePop[1]==LinePop[3])+(LinePop[2]==LinePop[3]),1)
return(res)
}
