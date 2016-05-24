#
#
#   create temporary file containing the data 
#
#
#scalefs0 <- 25
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
dtiso <- array(c(1,0,0,1,0,1),dim=c(6,ddim))
phi <- (1:1000)*2*pi/1000 
rad1 <- 7
rad2 <- 8
for(rad in seq(rad1,rad2,.25)){
  x <- as.integer(rad*sin(phi)+32.5)
  y <- as.integer(rad*cos(phi)+32.5)
  for(i in 1:125) ind[x[i],y[i],2:25] <- 0
  for(i in 126:250) ind[x[i],y[i],2:25] <- .6
  for(i in 251:375) ind[x[i],y[i],2:25] <- .2
  for(i in 376:500) ind[x[i],y[i],2:25] <- .8
  for(i in 501:625) ind[x[i],y[i],2:25] <- .4
  for(i in 626:750) ind[x[i],y[i],2:25] <- .7
  for(i in 751:875) ind[x[i],y[i],2:25] <- .3
  for(i in 876:1000) ind[x[i],y[i],2:25] <- .9
  for(i in 1:1000){
    etai <- etas[ind[x[i],y[i],12]*1000+1]
    dtiso[,x[i],y[i],] <- (1-etai)*c(0,0,0,0,0,1)+etai*c(1,0,0,1,0,1)
  }
}

rad1 <- 13
rad2 <- 14
for(rad in seq(rad1,rad2,.25)){
  sphi <- sin(phi)
  cphi <- cos(phi)
  x <- as.integer(rad*sphi+32.5)
  y <- as.integer(rad*cphi+32.5)
  for(i in 1:1000) ind[x[i],y[i],2:4] <- .6
  for(i in 1:1000) ind[x[i],y[i],5:7] <- .3
  for(i in 1:1000) ind[x[i],y[i],8:10] <- .9
  for(i in 1:1000) ind[x[i],y[i],11:13] <- 0
  for(i in 1:1000) ind[x[i],y[i],14:16] <- .5
  for(i in 1:1000) ind[x[i],y[i],17:19] <- .2
  for(i in 1:1000) ind[x[i],y[i],20:22] <- .8
  for(i in 1:1000) ind[x[i],y[i],23:25] <- .4 
  for(j in 1:26) {
    etai <- etas[ind[x[i],y[i],j]*1000+1]
    for(i in 1:1000) dtiso[,x[i],y[i],j] <- (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],0,sphi[i]^2,0,0)+etai*c(1,0,0,1,0,1)
  }
}
rad1 <- 19
rad2 <- 20
for(rad in seq(rad1,rad2,.25)){
  sphi <- sin(phi)
  cphi <- cos(phi)
  x <- as.integer(rad*sphi+32.5)
  y <- as.integer(rad*cphi+32.5)
  for(i in 1:1000) ind[x[i],y[i],2:3] <- .6
  for(i in 1:1000) ind[x[i],y[i],4:5] <- .2
  for(i in 1:1000) ind[x[i],y[i],6:7] <- .8
  for(i in 1:1000) ind[x[i],y[i],8:9] <- 0
  for(i in 1:1000) ind[x[i],y[i],10:11] <- .5
  for(i in 1:1000) ind[x[i],y[i],12:13] <- .9
  for(i in 1:1000) ind[x[i],y[i],14:15] <- .2
  for(i in 1:1000) ind[x[i],y[i],16:17] <- .6 
  for(i in 1:1000) ind[x[i],y[i],18:19] <- 0 
  for(i in 1:1000) ind[x[i],y[i],20:21] <- .9 
  for(i in 1:1000) ind[x[i],y[i],22:23] <- .3 
  for(i in 1:1000) ind[x[i],y[i],24:25] <- .7
  for(j in 1:26) {
    etai <- etas[ind[x[i],y[i],j]*1000+1]
    for(i in 1:1000) dtiso[,x[i],y[i],j] <- (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],0,sphi[i]^2,0,0)+etai*c(1,0,0,1,0,1)
  }
}

rad1 <- 25
rad2 <- 26
for(rad in seq(rad1,rad2,.25)){
   x <- as.integer(rad*sin(phi)+32.5)
   y <- as.integer(rad*cos(phi)+32.5)
   sphi <- sin(phi)
   cphi <- cos(phi)
   for(i in 1:1000) ind[x[i],y[i],3:5] <- (1+sphi[i])/3
   for(i in 1:1000) ind[x[i],y[i],6:8] <- (1+cphi[i])/3
   for(i in 1:1000) ind[x[i],y[i],9:11] <- 0
   for(i in 1:1000) ind[x[i],y[i],12:13] <- (1+sphi[i])/2
   for(i in 1:1000) ind[x[i],y[i],14:15] <- (1+cphi[i])/2
   for(i in 1:1000) ind[x[i],y[i],16:18] <- 0
   for(i in 1:1000) ind[x[i],y[i],19:21] <- (2+sphi[i])/3
   for(i in 1:1000) ind[x[i],y[i],22:24] <- (2+cphi[i])/3
   for(j in 1:26) {
      for(i in 1:1000){
         etai <- etas[ind[x[i],y[i],j]*1000+1]
         dtiso[,x[i],y[i],j] <- (1-etai)*c(cphi[i]^2,-sphi[i]*cphi[i],0,sphi[i]^2,0,0)+etai*c(1,0,0,1,0,1)
      }
     cat(".")
   }
}

for( i in 1:64) for (j in 1:64){
  if(max(ind[i,j,])==0&&((i-32.5)^2+(j-32.5)^2>26^2)) {
    dtiso[,i,j,]<-0
  }
}

#  yields a maximum eigenvalue of 2.5 within a cylinder of radius 26, zero tensor outside this cylinder
dtiso <- factor*dtiso


# now we want to view the projection of the zylinders onto a plane
project.cylinder <- function(obj,radius,phi=(1:1000)*2*pi/1000){
  ddim <- dim(obj)
  img <- matrix(0,length(phi),ddim[3])
  x <- as.integer(radius*sin(phi)+32.5)
  y <- as.integer(radius*cos(phi)+32.5)
  for ( i in 1:length(phi)) img[i,] <- obj[x[i],y[i],]
  img
}


# reset S0 image
s0offa <- read.table(system.file("dat/S0ofFA.txt",package="dti"))
s0 <- scalefs0*s0offa[as.integer(as.vector(ind)*500+1),2]
dim(s0) <- dim(ind)
#s0 <- array(as.integer(32000),dim(ind))
for( i in 1:64) for (j in 1:64){
  if(max(ind[i,j,])==0&&((i-32.5)^2+(j-32.5)^2>26^2)) s0[i,j,]<-0
}

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
createdata.dti(tmpfile2,dtiso,btb,s0,sigma*scalefs0)
