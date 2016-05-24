simul_data_YX <- function(totdim,ncomp) {
dimok <- FALSE
varsR <- c(10,8,6,4,2,1/2) 
varepsilon <- .01 
varsF <- c(.25,.125,.05,.0125,.005,.00125) 

if(totdim==1){stop("'totdim' must be > 1")}
if(ncomp==1){stop("'ncomp' must be > 1")}

if(totdim==2){
dimok <- TRUE
if(!(ncomp %in% 1:totdim)){
cat(paste("ncomp must be <= ",totdim,"\n"))
cat(paste("ncomp was set to ",totdim,"\n"))
ncomp <- totdim
}
ksi1 <- c(1,1)/sqrt(2)
ksi2 <- c(1,-1)/sqrt(2)
ksi <- cbind(ksi1,ksi2)[1:totdim,1:ncomp]
}

if(totdim==3){
dimok <- TRUE
if(!(ncomp %in% 1:totdim)){
cat(paste("ncomp must be <= ",totdim,"\n"))
cat(paste("ncomp was set to ",totdim,"\n"))
ncomp <- totdim
}
ksi1 <- c(1,1,1)/sqrt(3)
ksi2 <- c(-1/2,-1/2,1)/sqrt(3/2)
ksi3 <- c(-1,1,0)/sqrt(2)
ksi <- cbind(ksi1,ksi2,ksi3)[1:totdim,1:ncomp]
}

if(totdim==4){
dimok <- TRUE
if(!(ncomp %in% 1:totdim)){
cat(paste("ncomp must be <= ",totdim,"\n"))
cat(paste("ncomp was set to ",totdim,"\n"))
ncomp <- totdim
}
ksi1 <- c(1,1,1,1)/2
ksi2 <- c(1,-1,1,-1)/2
ksi3 <- c(1,1,-1,-1)/2
ksi4 <- c(1,-1,-1,1)/2
ksi <- cbind(ksi1,ksi2,ksi3,ksi4)[1:totdim,1:ncomp]
}

if(totdim==5){
dimok <- TRUE
if(!(ncomp %in% 1:totdim)){
cat(paste("ncomp must be <= ",totdim,"\n"))
cat(paste("ncomp was set to ",totdim,"\n"))
ncomp <- totdim
}
ksi1 <- c(1,1,1,1,1)/sqrt(5)
ksi2 <- c(1,-1,1,-1,0)/2
ksi3 <- c(1,1,-1,-1,0)/2
ksi4 <- c(1,-1,-1,1,0)/2
ksi5 <- c(1/4,1/4,1/4,1/4,-1)/sqrt(5)*2
ksi <- cbind(ksi1,ksi2,ksi3,ksi4,ksi5)[1:totdim,1:ncomp]
}

if((totdim %% 6 == 0)&(totdim>=6)){
dimok <- TRUE
if(!(ncomp %in% 1:6)){
cat(paste("ncomp must be <= ",6,"\n"))
cat(paste("ncomp was set to ",6,"\n"))
ncomp <- 6
}
ksi1 <- rep(c(1,1,1,1,1,1),totdim/6)/sqrt(totdim)
ksi2 <- rep(c(-1/2,-1/2,1,-1/2,-1/2,1),totdim/6)/sqrt(totdim/2)
ksi3 <- rep(c(-1,1,0,-1,1,0),totdim/6)/sqrt(2*totdim/3)
ksi4 <- rep(c(1,1,1,-1,-1,-1),totdim/6)/sqrt(totdim)
ksi5 <- rep(c(-1/2,-1/2,1,1/2,1/2,-1),totdim/6)/sqrt(totdim/2)
ksi6 <- rep(c(-1,1,0,1,-1,0),totdim/6)/sqrt(2*totdim/3)
ksi <- cbind(ksi1,ksi2,ksi3,ksi4,ksi5,ksi6)[1:totdim,1:ncomp]
}

if((totdim %% 6 == 1)&(totdim>=6)){
dimok <- TRUE
if(!(ncomp %in% 1:6)){
cat(paste("ncomp must be <= ",6,"\n"))
cat(paste("ncomp was set to ",6,"\n"))
ncomp <- 6
}
ksi1 <- c(1,1,1,1,1,1,1,rep(c(1,1,1,1,1,1),(totdim-1)/6-1))/sqrt(totdim)
ksi2 <- c(-1/2,-1/2,1,1,-1,1,-1,rep(c(-1/2,-1/2,1,-1/2,-1/2,1),(totdim-1)/6-1))/sqrt(3*((totdim-1)/6-1)+11/2)
ksi3 <- c(-1,1,0,1,1,-1,-1,rep(c(-1,1,0,-1,1,0),(totdim-1)/6-1))/sqrt(4*((totdim-1)/6-1)+6)
ksi4 <- c(-8/11,-8/11,16/11,-6/11,6/11,-6/11,6/11,rep(c(1,1,1,-1,-1,-1),(totdim-1)/6-1))/sqrt(6*((totdim-1)/6-1)+4*sqrt(33)/11)
ksi5 <- c(2/3,-2/3,0,1/3,1/3,-1/3,-1/3,rep(c(-1/2,-1/2,1,1/2,1/2,-1),(totdim-1)/6-1))/sqrt(3*((totdim-1)/6-1)+2*sqrt(3)/3)
ksi6 <- c(0,0,0,1,-1,-1,1,rep(c(-1,1,0,1,-1,0),(totdim-1)/6-1))/sqrt(4*((totdim-1)/6-1)+4)
ksi <- cbind(ksi1,ksi2,ksi3,ksi4,ksi5,ksi6)[1:totdim,1:ncomp]
}

if((totdim %% 6 == 2)&(totdim>=6)){
dimok <- TRUE
if(!(ncomp %in% 1:6)){
cat(paste("ncomp must be <= ",6,"\n"))
cat(paste("ncomp was set to ",6,"\n"))
ncomp <- 6
}
ksi1 <- c(1,1,1,1,1,1,1,1,rep(c(1,1,1,1,1,1),(totdim-2)/6-1))/sqrt(totdim)
ksi2 <- c(1,-1,1,-1,1,-1,1,-1,rep(c(-1/2,-1/2,1,-1/2,-1/2,1),(totdim-2)/6-1))/sqrt(3*((totdim-2)/6-1)+8)
ksi3 <- c(1,1,-1,-1,1,1,-1,-1,rep(c(-1,1,0,-1,1,0),(totdim-2)/6-1))/sqrt(4*((totdim-2)/6-1)+8)
ksi4 <- c(1,-1,-1,1,1,-1,-1,1,rep(c(1,1,1,-1,-1,-1),(totdim-2)/6-1))/sqrt(6*((totdim-2)/6-1)+8)
ksi5 <- c(1,1,1,1,-1,-1,-1,-1,rep(c(-1/2,-1/2,1,1/2,1/2,-1),(totdim-2)/6-1))/sqrt(3*((totdim-2)/6-1)+8)
ksi6 <- c(1,-1,1,-1,-1,1,-1,1,rep(c(-1,1,0,1,-1,0),(totdim-2)/6-1))/sqrt(4*((totdim-2)/6-1)+8)
ksi <- cbind(ksi1,ksi2,ksi3,ksi4,ksi5,ksi6)[1:totdim,1:ncomp]
}

if((totdim %% 6 == 3)&(totdim>=6)){
dimok <- TRUE
if(!(ncomp %in% 1:6)){
cat(paste("ncomp must be <= ",6,"\n"))
cat(paste("ncomp was set to ",6,"\n"))
ncomp <- 6
}
ksi1 <- c(1,1,1,1,1,1,1,1,1,rep(c(1,1,1,1,1,1),(totdim-3)/6-1))/sqrt(totdim)
ksi2 <- c(-1/2,-1/2,1,-1/2,-1/2,1,-1/2,-1/2,1,rep(c(-1/2,-1/2,1,-1/2,-1/2,1),(totdim-3)/6-1))/sqrt(3*((totdim-3)/6-1)+9/2)
ksi3 <- c(-1,1,0,-1,1,0,-1,1,0,rep(c(-1,1,0,-1,1,0),(totdim-3)/6-1))/sqrt(4*((totdim-3)/6-1)+6)
ksi4 <- c(-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,1,1,1,rep(c(1,1,1,-1,-1,-1),(totdim-3)/6-1))/sqrt(6*((totdim-3)/6-1)+9/2)
ksi5 <- c(1/4,1/4,-1/2,1/4,1/4,-1/2,-1/2,-1/2,1,rep(c(-1/2,-1/2,1,1/2,1/2,-1),(totdim-3)/6-1))/sqrt(3*((totdim-3)/6-1)+9/4)
ksi6 <- c(1/2,-1/2,0,1/2,-1/2,0,-1,1,0,rep(c(-1,1,0,1,-1,0),(totdim-3)/6-1))/sqrt(4*((totdim-3)/6-1)+3)
ksi <- cbind(ksi1,ksi2,ksi3,ksi4,ksi5,ksi6)[1:totdim,1:ncomp]
}

if((totdim %% 6 == 4)&(totdim>=6)){
dimok <- TRUE
if(!(ncomp %in% 1:6)){
cat(paste("ncomp must be <= ",6,"\n"))
cat(paste("ncomp was set to ",6,"\n"))
ncomp <- 6
}
ksi1 <- c(1,1,1,1,1,1,1,1,1,1,rep(c(1,1,1,1,1,1),(totdim-4)/6-1))/sqrt(totdim)
ksi2 <- c(1,-1,1,-1,0,1,-1,1,-1,0,rep(c(-1/2,-1/2,1,-1/2,-1/2,1),(totdim-4)/6-1))/sqrt(3*((totdim-4)/6-1)+8)
ksi3 <- c(1,1,-1,-1,0,1,1,-1,-1,0,rep(c(-1,1,0,-1,1,0),(totdim-4)/6-1))/sqrt(4*((totdim-4)/6-1)+8)
ksi4 <- c(1,-1,-1,1,0,1,-1,-1,1,0,rep(c(1,1,1,-1,-1,-1),(totdim-4)/6-1))/sqrt(6*((totdim-4)/6-1)+8)
ksi5 <- c(1/4,1/4,1/4,1/4,-1,1/4,1/4,1/4,1/4,-1,rep(c(-1/2,-1/2,1,1/2,1/2,-1),(totdim-4)/6-1))/sqrt(3*((totdim-4)/6-1)+5/2)
ksi6 <- c(1,1,1,1,1,-1,-1,-1,-1,-1,rep(c(-1,1,0,1,-1,0),(totdim-4)/6-1))/sqrt(4*((totdim-4)/6-1)+10)
ksi <- cbind(ksi1,ksi2,ksi3,ksi4,ksi5,ksi6)[1:totdim,1:ncomp]
}

if((totdim %% 6 == 5)&(totdim>=6)){
dimok <- TRUE
if(!(ncomp %in% 1:6)){
cat(paste("ncomp must be <= ",6,"\n"))
cat(paste("ncomp was set to ",6,"\n"))
ncomp <- 6
}
ksi1 <- c(1,1,1,1,1,1,1,1,1,1,1,rep(c(1,1,1,1,1,1),(totdim-5)/6-1))/sqrt(totdim)
ksi2 <- c(1,-1,1,-1,0,1,-1,1,-1,0,0,rep(c(-1/2,-1/2,1,-1/2,-1/2,1),(totdim-5)/6-1))/sqrt(3*((totdim-5)/6-1)+8)
ksi3 <- c(1,1,-1,-1,0,1,1,-1,-1,0,0,rep(c(-1,1,0,-1,1,0),(totdim-5)/6-1))/sqrt(4*((totdim-5)/6-1)+8)
ksi4 <- c(10/11,-12/11,-12/11,10/11,-1/11,10/11,-12/11,-12/11,10/11,-1/11,10/11,rep(c(1,1,1,-1,-1,-1),(totdim-5)/6-1))/sqrt(6*((totdim-5)/6-1)+98/11)
ksi5 <- c(13/196,53/196,53/196,13/196,-53/49,13/196,53/196,53/196,13/196,-53/49,40/49,rep(c(-1/2,-1/2,1,1/2,1/2,-1),(totdim-5)/6-1))/sqrt(3*((totdim-5)/6-1)+325/98)
ksi6 <- c(4/5,62/65,62/65,4/5,77/65,-6/5,-68/65,-68/65,-6/5,-53/65,8/13,rep(c(-1,1,0,1,-1,0),(totdim-5)/6-1))/sqrt(4*((totdim-5)/6-1)+138/13)
ksi <- cbind(ksi1,ksi2,ksi3,ksi4,ksi5,ksi6)[1:totdim,1:ncomp]
}

if(!dimok) {stop("Incorrect value for 'totdim'. 'totdim' must be > 1")}

epsilon <- stats::rnorm(totdim,mean=rep(0,totdim),sd=varepsilon)

r <- stats::rnorm(ncomp,mean=rep(0,ncomp),sd=varsR[1:ncomp])

simX <- r%*%t(ksi)+epsilon


if(ncomp==2) {
HH <- 3
eta21 <- c(1,2,1)/sqrt(6)
eta22 <- c(0,1,-2)/sqrt(5)
eta <- cbind(eta21,eta22)
}

if(ncomp==3) {
HH <- 3
eta31 <- c(1,2,1)/sqrt(6)
eta32 <- c(0,1,-2)/sqrt(5)
eta33 <- c(-5,2,1)/sqrt(30)
eta <- cbind(eta31,eta32,eta33)
}

if(ncomp==4) {
HH <- 4
eta41 <- c(1,2,1,0)/sqrt(6)
eta42 <- c(0,1,-2,1)/sqrt(5)
eta43 <- c(0,1,-2,-5)/sqrt(30)
eta44 <- c(-5,2,1,0)/sqrt(30)
eta <- cbind(eta41,eta42,eta43,eta44)
}

if(ncomp==5) {
HH <- 4
eta51 <- c(1,1,1,1)/2
eta52 <- c(1,1,1,1)/2
eta53 <- c(1,1,1,1)/2
eta54 <- c(1,1,1,1)/2
eta55 <- c(1,1,1,1)/2
eta <- cbind(eta51,eta52,eta53,eta54,eta55)
}

if(ncomp==6) {
HH <- 4
eta61 <- c(1,1,1,1)/2
eta62 <- c(1,1,1,1)/2
eta63 <- c(1,1,1,1)/2
eta64 <- c(1,1,1,1)/2
eta65 <- c(1,1,1,1)/2
eta66 <- c(1,1,1,1)/2
eta <- cbind(eta61,eta62,eta63,eta64,eta65,eta66)
}

f <- stats::rnorm(ncomp,mean=rep(0,ncomp),sd=varsF[1:ncomp])
z <- f+r


sigmaScarre <- .001
lambda <- .6

sigmaPsi <- sigmaScarre*((1-lambda)*diag(rep(1,HH))+lambda*rep(1,HH)%*%t(rep(1,HH)))
Psi <- mvtnorm::rmvnorm(1,mean=rep(0,HH),sigma=sigmaPsi)


Y <- z%*%t(eta)+Psi

res <- c(Y,simX)
names(res) <- c(paste("Y",1:HH),paste("X",1:totdim,sep=""))

return(res)
}
