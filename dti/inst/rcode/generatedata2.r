#
#
#   create temporary file containing the data 
#
#
#scalefs0 <- 20
#
#  Scalefactor for images to avoid discritisation problems
#
bvec <- t(bvec)
btb <- matrix(0,6,dim(bvec)[2])
btb[1,] <- bvec[1,]*bvec[1,]*bvalue/1000
btb[4,] <- bvec[2,]*bvec[2,]*bvalue/1000
btb[6,] <- bvec[3,]*bvec[3,]*bvalue/1000
btb[2,] <- 2*bvec[1,]*bvec[2,]*bvalue/1000
btb[3,] <- 2*bvec[1,]*bvec[3,]*bvalue/1000
btb[5,] <- 2*bvec[2,]*bvec[3,]*bvalue/1000

# a useful function to create a tensor of specified anisotropy
eta <- function(ai){
  aindex <- function(tensor){
    values <- eigen(matrix(tensor[c(1,2,3,2,4,5,3,5,6)],3,3))$values
    sqrt(3/2*sum((values-mean(values))^2)/sum(values^2))
  }
  risk <- function(par,ai) (ai-aindex((1-par)*c(1,0,0,0,0,0)+par*c(1,0,0,1,0,1)))^2
  optimize(f = risk, interval = c(0,1),ai=ai)$minimum
}

#
# create Phantom data
#
etas <- numeric(1001)
for(i in 1:1001) etas[i] <- eta((i-1)/1000)

ind <- array(0,ddim)
ind[seq(3,65,3),,seq(2,26,3)] <- .75
ind[,seq(3,65,3),seq(1,26,3)] <- .75
ind[seq(2,63,3),seq(2,63,3),] <- .75
dtiso <- array(c(1,0,0,1,0,1),dim=c(6,ddim))
etai <- eta(.6)
dtiso[,seq(3,65,6),,seq(2,26,6)] <- (1-etai)*c(0,0,0,1,0,0)+etai*c(1,0,0,1,0,1)
dtiso[,,seq(3,65,6),seq(1,26,6)] <- (1-etai)*c(1,0,0,0,0,0)+etai*c(1,0,0,1,0,1)
dtiso[,seq(2,63,6),seq(2,63,6),] <- (1-etai)*c(0,0,0,0,0,1)+etai*c(1,0,0,1,0,1)
etai <- eta(.7)
dtiso[,seq(6,65,6),,seq(2,26,6)] <- (1-etai)*c(0,0,0,1,0,0)+etai*c(1,0,0,1,0,1)
dtiso[,,seq(6,65,6),seq(1,26,6)] <- (1-etai)*c(1,0,0,0,0,0)+etai*c(1,0,0,1,0,1)
dtiso[,seq(5,63,6),seq(2,63,6),] <- (1-etai)*c(0,0,0,0,0,1)+etai*c(1,0,0,1,0,1)
etai <- eta(.8)
dtiso[,seq(3,65,6),,seq(5,26,6)] <- (1-etai)*c(0,0,0,1,0,0)+etai*c(1,0,0,1,0,1)
dtiso[,,seq(3,65,6),seq(4,26,6)] <- (1-etai)*c(1,0,0,0,0,0)+etai*c(1,0,0,1,0,1)
dtiso[,seq(2,63,6),seq(5,63,6),] <- (1-etai)*c(0,0,0,0,0,1)+etai*c(1,0,0,1,0,1)
etai <- eta(.9)
dtiso[,seq(6,65,6),,seq(5,26,6)] <- (1-etai)*c(0,0,0,1,0,0)+etai*c(1,0,0,1,0,1)
dtiso[,,seq(6,65,6),seq(4,26,6)] <- (1-etai)*c(1,0,0,0,0,0)+etai*c(1,0,0,1,0,1)
dtiso[,seq(5,63,6),seq(5,63,6),] <- (1-etai)*c(0,0,0,0,0,1)+etai*c(1,0,0,1,0,1)
dtiso <- factor*dtiso


# now we want to view the projection of the zylinders onto a plane
# reset S0 image
s0offa <- read.table(system.file("dat/S0ofFA.txt",package="dti"))
s0 <- scalefs0*s0offa[as.integer(as.vector(ind)*500+1),2]
dim(s0) <- dim(ind)

# create noisy data
createdata.dti <- function(file,dtensor,btb,s0,sigma,level=250){
  ngrad <- dim(btb)[2]
  ddim <- dim(s0)
  dim(dtensor)<-c(6,prod(ddim))
  dtensor <- t(dtensor)
  si <- exp(-dtensor%*%btb)*as.vector(s0)
  dim(si)<-c(ddim,ngrad)
  for (j in 1:ngrad) {
    for (i in 1:ddim[3]) {
      si[,,i,j] <- abs(fft(fft(si[,,i,j])+complex(real=rnorm(ddim[1]*ddim[2],0,sigma),imaginary=rnorm(ddim[1]*ddim[2],0,sigma)),inverse=TRUE))/ddim[1]/ddim[2]
    }
  }
  con <- file(file,"wb")
  writeBin(as.integer(si),con,2)
  close(con)
}


#   create phantom - object
tmpfile1 <- tempfile("S_all")
createdata.dti(tmpfile1,dtiso,btb,s0,0.5)


# create noisy data
cat("Creating noisy data with standard deviation ",sigma,"\n")
#set.seed(1)
tmpfile2 <- tempfile("S_noise_all")
createdata.dti(tmpfile2,dtiso,btb,s0,scalefs0*sigma)
