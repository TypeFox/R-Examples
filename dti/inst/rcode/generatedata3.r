#
#
#   create temporary file containing the data 
#
#
#scalefs0 <- 25
#
#  Scalefactor for images to avoid discritisation problems
#
cat("Using data set 3 \n")

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
etas <- numeric(2001)
for(i in 1:2001) etas[i] <- .9
dtiso <- array(c(1,0,0,1,0,1),dim=c(6,ddim))
phi <- (1:1000)*6*pi/1000 
etai <- .2



phi <- (1:2000)*7.6*pi/2000 
rad1 <- 20
rad2 <- 21
h <- 1
for(rad in seq(rad1,rad2,.5)){
  sphi <- sin(phi)
  cphi <- cos(phi)
  x <- as.integer(rad*sphi+32.5)
  y <- as.integer(rad*cphi+32.5)
  z <- as.integer(pmax(1,pmin(26,h*phi+1)))
    for(i in 1:2000) dtiso[,x[i],y[i],z[i]] <-
        (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],2.5*h/rad*cphi[i],
          sphi[i]^2,-2.5*h/rad*sphi[i],(2.5*h/rad)^2)+etai*c(1,0,0,1,0,1)
    for(i in 1:2000) dtiso[,x[i],y[i],z[i]+1] <-
        (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],2.5*h/rad*cphi[i],
          sphi[i]^2,-2.5*h/rad*sphi[i],(2.5*h/rad)^2)+etai*c(1,0,0,1,0,1)
   }
for(rad in seq(rad1,rad2,.5)){
  sphi <- sin(phi)
  cphi <- cos(phi)
  x <- as.integer(rad*sphi+32.5)
  y <- as.integer(rad*cphi+32.5)
  z <- as.integer(pmax(1,pmin(26,h*phi+1.5)))
    for(i in 1:2000) dtiso[,x[i],y[i],z[i]] <-
        (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],2.5*h/rad*cphi[i],
          sphi[i]^2,-2.5*h/rad*sphi[i],(2.5*h/rad)^2)+etai*c(1,0,0,1,0,1)
    for(i in 1:2000) dtiso[,x[i],y[i],z[i]+1] <-
        (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],2.5*h/rad*cphi[i],
          sphi[i]^2,-2.5*h/rad*sphi[i],(2.5*h/rad)^2)+etai*c(1,0,0,1,0,1)
   }

for(rad in seq(rad1,rad2,.5)){
  sphi <- sin(phi)
  cphi <- cos(phi)
  x <- as.integer(rad*sphi+32.5)
  y <- as.integer(rad*cphi+32.5)
  z <- as.integer(pmax(1,pmin(26,h*phi+2)))
    for(i in 1:2000) dtiso[,x[i],y[i],z[i]] <-
        (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],2.5*h/rad*cphi[i],
          sphi[i]^2,-2.5*h/rad*sphi[i],(2.5*h/rad)^2)+etai*c(1,0,0,1,0,1)
    for(i in 1:2000) dtiso[,x[i],y[i],z[i]+1] <-
        (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],2.5*h/rad*cphi[i],
          sphi[i]^2,-2.5*h/rad*sphi[i],(2.5*h/rad)^2)+etai*c(1,0,0,1,0,1)
   }


# now we want to view the projection of the zylinders onto a plane
# reset S0 image
s0offa <- read.table(system.file("dat/S0ofFA.txt",package="dti"))
s0 <- array(s0offa[450,2]*scalefs0+rnorm(prod(ddim)),ddim)

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
createdata.dti(tmpfile1,dtiso,btb,s0,1)


# create noisy data
cat("Creating noisy data with standard deviation ",sigma,"\n")
#set.seed(1)
tmpfile2 <- tempfile("S_noise_all")
createdata.dti(tmpfile2,dtiso,btb,s0,scalefs0*sigma)
