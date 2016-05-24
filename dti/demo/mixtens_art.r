cat("Demo for mix-tensor models (single shell, bvalue 1000 \n")
set.seed(1)
source(system.file("rcode/gen_mixtens.r",package="dti"))
ngrad <- readline("Provide number of gradients (default: 136 minimum: 6  maximum: 162):")

ngrad <- if(ngrad!="") as.numeric(ngrad) else 136
if(is.na(ngrad)) ngrad <- 136
ngrad <- max(6,min(162,ngrad))
ns0 <- max(1,(ngrad+10)%/%20) 

cat("Using ",ns0," S0  and ",ngrad,"diffusion weighted images\n")
# read the gradient data, these are 25 gradient directions + one non-zero weighted
# bvec <- read.table(system.file("dat/b-directions.txt",package="dti"))
data(optgradients)
grad <- cbind(matrix(0,3,ns0),optgrad[[ngrad-5]])
# grad <- t(read.table(system.file("dat/hardi-grad.txt",package="dti")))
n <- readline("size of data cube (n x n x n) (default n=6)") 

if (!is.null(n)) n <- as.numeric(n) else n <- 6
if(is.na(n)) n <- 6 else n <- as.integer(min(10,max(1,n)))

snr <- readline("signal-noise-ratio (default 50)") 

snr <- if (snr!="") snr <- as.numeric(snr) else 50
if(is.na(snr)) snr <- 50 
sigma <- 1/snr

scenario <- readline("select scenario \n 1: varying directions (default)
                                      \n 2: varying mixtures
                                      \n 3: varying fa
                                      \n 4: homogeneous")

if (is.null(scenario)) scenario <- 1 else scenario <- as.numeric(scenario)
if(is.na(scenario)) scenario <- 1 else scenario <- as.integer(max(1,min(4,scenario)))

scenario <- max(1,min(4,as.integer(scenario)))
mix <- array(0,c(4,n,n,n))
th <- array(0,c(2,n,n,n))
alpha <- array(0,c(3,n,n,n))
beta <- array(0,c(3,n,n,n))
mix[1,,,] <- 0
if(scenario==1){
th[1,,,] <- 4*1e-3
th[2,,,] <- .8*1e-3
alpha[1,,,] <- 0
beta[1,,,] <- 0
# first component in x -direction
for(i in 1:n) alpha[1,i,,] <- pi/2
for(i in 1:n) beta[1,i,,] <- 0
# second component in x-y -plane
for(i in 1:n) alpha[2,,,i] <- pi/2
for(i in 1:n) beta[2,,,i] <- pi/2*(i-1)/(n-1)
# third component in x-z -plane
for(i in 1:n) alpha[3,,i,] <- pi/2*(1-(i-1)/(n-1))
for(i in 1:n) beta[3,,i,] <- 0
} else { 
angle <- readline("angle between directions (in radiant) default (and maximum) pi/2\n ") 

if (!is.null(angle)) angle <- as.numeric(angle) else angle <- pi/2
if(is.na(angle)) angle <- pi/2 else angle <- min(pi/2,max(0,angle))
# first component in x -direction
for(i in 1:n) alpha[1,i,,] <- pi/2  
for(i in 1:n) beta[1,i,,] <- 0
# second component in x-y -plane
for(i in 1:n) alpha[2,,i,] <- pi/2
for(i in 1:n) beta[2,,i,] <- angle  
# third component in x-z -plane
for(i in 1:n) alpha[3,,,i] <- pi/2-angle
for(i in 1:n) beta[3,,,i] <- 0
}
if(scenario!=2){
mix1 <- readline("first mixture coefficient (default 1/3)\n
                  0 corresponds to mixtures of order 2, 1 to order 1") 

if (!is.null(mix1)) mix1 <- as.numeric(mix1) else mix1 <- 1/3
if(is.na(mix1)) mix1 <- 1/3 else mix1 <- min(1,max(0,mix1))

mix[2,,,] <- mix1
mix[3:4,,,] <- (1-mix1)/2
} else {
for(i in 1:n) mix[2,i,,] <- i/2/(n+1) 
for(i in 1:n) mix[3,,i,] <- i/2/(n+1) 
mix[4,,,] <- 1-mix[2,,,]-mix[3,,,]
}
if(scenario!=3){
fa <- readline("fa (default fa = .8)") 

if (!is.null(fa)) fa <- as.numeric(fa) else fa <- .8
if(is.na(fa)) fa <- .8 else fa <- min(.9,max(.2,fa))
} else {
fa <- seq(.5,.9,length=n^3)
}
cfa <- fa*fa/(1-fa*fa)
l1 <- cfa+sqrt(cfa^2+3*cfa)
l2 <- (3.2/(1+l1))^(1/3)
th[2,,,] <- l2*1e-3
th[1,,,] <- l1*l2*1e-3
bvalue <- c(rep(0,ns0),rep(1000,ngrad))

maxcomp <- readline("maximal order of mix-tensor model (default 3)") 

if (is.null(maxcomp)) maxcomp <- 3 else maxcomp <- as.numeric(maxcomp)
if(is.na(maxcomp)) maxcomp <- 3 else maxcomp <- min(5,max(1,maxcomp))

z0 <- truemixtens(mix,th,alpha,beta,grad,bvalue,sigma,ns0)
z <- tdatamixtens(mix,th,alpha,beta,grad,bvalue,sigma,ns0)
zt <- dtiTensor(z)
zmix <- dwiMixtensor(z,maxcomp=1)
summary(zmix)
if(maxcomp>1) for(m in 2:maxcomp){
zmix0 <- dwiMixtensor(z,maxcomp=m)
zmix <- dwiMtCombine(zmix0,zmix)
summary(zmix)
}
gslexists <- "gsl" %in% .packages(TRUE)
if(gslexists) zqball <- dwiQball(z,order=8,lambda=2e-2)

size <- as.integer(min(.adimpro$xsize/3.2,.adimpro$ysize/2.4))

vodf <- readline("Visualize and compare estimated ODF's (Y/N) :")

if(toupper(vodf)!="N"){
source(system.file("rcode/mousecallbacks.r",package="dti"))
size <- 500
w1 <- show3d(z0,scale=.5,maxobjects=n^3,FOV=1,windowRect = c(1, 1, size, size),odfscale=1)
w2 <- show3d(zt,what="odf",minalpha=1,minfa=0,scale=.5,maxobjects=n^3,FOV=1,windowRect = c(size+11, 1, 2*size+10, size),odfscale=1)
w3 <- show3d(zmix,scale=.5,maxobjects=n^3,FOV=1,windowRect = c(1, size+11 , size, 2*size+10),odfscale=1)
if(gslexists){ 
w4 <- show3d(zqball,scale=.5,maxobjects=n^3,FOV=1,windowRect = c(size+11, size+11, 2*size+10, 2*size+10),odfscale=1)
mouseTrackball(dev=c(w1,w2,w3,w4))
mouseZoom(2,dev=c(w1,w2,w3,w4))
mouseFOV(3,dev=c(w1,w2,w3,w4))
cat("True ODF in device",w1,"\n")
cat("Estimated tensor ODF in device",w2,"\n")
cat("Estimated mixtensor ODF in ",w3,"\n")
cat("Estimated q-Ball ODF in ",w4,"\n")
} else {
mouseTrackball(dev=c(w1,w2,w3))
mouseZoom(2,dev=c(w1,w2,w3))
mouseFOV(3,dev=c(w1,w2,w3))
cat("True ODF in device",w1,"\n")
cat("Estimated tensor ODF in device",w2,"\n")
cat("Estimated mixtensor ODF in ",w3,"\n")
}
}

vfa <- readline("Visualize and compare estimated FA's (Y/N) :")

if(toupper(vfa)!="N"){

X11()
par(mfrow=c(3,n),mar=c(1,1,2,.1),mgp=c(2,1,0))
fa0 <- extract(z0,"fa")$fa
tfa <- extract(zt,"fa")$fa
mfa <- extract(zmix,"fa")$fa
for(i in 1:n) image(fa0[,,i],col=grey((0:255)/255),zlim=c(0,1),main=paste("True fa (Mixture) sl.,",i),xaxt="n",yaxt="n")
for(i in 1:n) image(tfa[,,i],col=grey((0:255)/255),zlim=c(0,1),main=paste("Est. fa (Tensor) sl.,",i),xaxt="n",yaxt="n")
for(i in 1:n) image(mfa[,,i],col=grey((0:255)/255),zlim=c(0,1),main=paste("Est. fa (Mixture) sl.,",i),xaxt="n",yaxt="n")
}


