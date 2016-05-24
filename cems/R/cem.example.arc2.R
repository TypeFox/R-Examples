
cem.example.arc2 <- function(n=150, noise=0.2, type=2, sigmaX= 0.1, sigmaY=0.1,
init = 0, plotEach = 1, file="output.png", testinit=F, nI=100){

#png(file)

### Create arc data set

#create arc
#phi <- runif(n)*pi 
#arc <- cbind(cos(2*phi), sin(phi)) * (1+rnorm(n) * noise)
#arc0 <- cbind(cos(2*phi), sin(phi))

phi<-runif(n)*pi
arc<-cbind(phi,sin(0.5*phi^3)/(phi^2))
nor<-cbind((1.5*(phi^3)*cos(0.5*phi^3)-2*sin(0.5*phi^3))/(phi^3),-phi/phi)
for(i in seq(1,nrow(nor))){
	t<-(nor[i,1]^2+nor[i,2]^2)^0.5;
	nor[i,]<-nor[i,]/t
}
ns=noise
ns=runif(n,-ns,ns)
arc0<-arc
#arc<-(arc+(1+rnorm(n) * noise)*ns*nor)
arc<-(arc+rnorm(n) * noise*nor)

### write data to file
y=arc
zz <- file("Y", "wb")
for(i in seq(1:dim(y)[1]))
{
	writeBin(as.double(y[i,1]), zz)
	writeBin(as.double(y[i,2]), zz)
}
close(zz)





#order by angle for plotting curves    
o0 <- order(phi) 

#test data
#phit <- runif(n*10)*pi 
#arct <- cbind(cos(phit), sin(phit)) * (1+rnorm(10*n) * noise)

phit <- runif(n*10)*pi 
arct<-cbind(phit,sin(0.5*phit^3)/(phit^2))
nort<-cbind((1.5*(phit^3)*cos(0.5*phit^3)-2*sin(0.5*phit^3))/(phit^3),-phit/phit)
for(i in seq(1,nrow(nort))){
	t<-(nort[i,1]^2+nort[i,2]^2)^0.5;
	nort[i,]<-nort[i,]/t
}
ns=noise
ns=runif(10*n,-ns,ns)
#arct<-(arct+(1+rnorm(10*n) * noise)*ns*nort)
arct<-(arct+rnorm(10*n) * noise*nort)





### Initialization

if(init == 0){
 #ground truth initialization
 z <- phi
}else if(init == 1){
 #random initialization
 z <- runif(n)
}
else{
 #close to principal component initialization
 z = arc[,2]
}

#compute initial principal curve
#do not optimize sigmaX after each itertion to show optimization path
pc <-  cem(y=arc, knnX=n, iter=0, optimalSigmaX=F)
#pc <-  cem(y=arc, z=z, knnX=50, knnY=50, iter=0, optimalSigmaX=T)

#set type
pc$type=type

#set intial bandwidth
if(sigmaX){
	pc$sigmaX=rep(sigmaX,length(pc$sigmaX));
} #smoothness of initial curve 
if(sigmaY){
	pc$sigmaY=rep(sigmaY,length(pc$sigmaY));
}
#pc$sigmaY=sigmaY;


#initial prediction (curve)
xi <- predict(pc)
yi <- predict(pc, xi)
oi <- order(xi)

#create sample from initial curve for plotting
r = range(xi)
xpi = seq(r[1], r[2], length.out=500)
ypi = predict(pc, xpi)


if( pc$type > 0 ){
 colName <- "dodgerblue2"
}else{
 colName <- "gold2"
}

tmp <- col2rgb(colName)
col <- rgb(tmp[1], tmp[2], tmp[3], 100, maxColorValue=255)
col2 <- rgb(tmp[1], tmp[2], tmp[3], 255, maxColorValue=255)

lty <- 1
lwd <- 8
pcs <- list()
o <- list()
x <- list()
sigmaX <- list()
sigmaY <- list()
iter <- list()
y <- list()
yt <- list()
mse = NULL
mset = NULL
ortho = NULL
orthot = NULL

lwi = 6
lwg = 6
lwp = 2



### Optimization

#run nIter iterations between plotting
nIter = plotEach


#plot ground truth and initial curve
par(mar=c(5,5,4,2))
plot(arc, xlab=expression(y[1]), ylab=expression(y[2]), col = "#00000010",
       pch=19, asp=1, cex.lab=1.75, cex.axis=1.75, cex=2, bty="n")
lines(arc0[o0,], lwd=lwg, col="black", lty=6)
lines(ypi$y, col="darkgray", lwd=lwi, lty=1)

### write init data to file
tt=ypi$y
zz <- file("Cinit", "wb")
for(i in seq(1:dim(tt)[1]))
{
	writeBin(as.double(tt[i,1]), zz)
	writeBin(as.double(tt[i,2]), zz)
}
close(zz)

if(testinit) break


#cross validation flag (selected curve)
selected = -1


#run a 100 iterations, one at a time 
#for running the whole optimization in one go run either
# pc <-  cem(y=arc, z=z, knnX=n, knnY=50, iter=100)
# pc <- cem.optimize(pc, stepX=1, stepBW=0.1, iter=100)
# here one iterations is run at the time and corss-validation is performed on a
# test set. For minimizing orthogonality cross-validation appears to be not
# necessary.

sigmaX<-list();
sigmaY<-list();

for(k in 1:nI){
  #run one iterations  
  pc <- cem.optimize(pc, stepX=1, stepBW=0.1, iter=1, verbose=2, optimalSigmaX=F)
  
  #store cem values at iteration k
  
  sigmaX[[k]] <- pc$sigmaX
  sigmaY[[k]] <- pc$sigmaY
  

  x[[k]] <-predict(pc)
  
  #train data
  r = range(x[[k]])
  xp = seq(r[1], r[2], length.out=500)
  yp = predict(pc, xp);
  o[[k]] <- order(x[[k]]);
  y[[k]] <- predict(pc, x[[k]])

  #test data
  xt <-predict(pc, arct)
  yt[[k]] <-predict(pc, xt)
  
  #compute mean squared errors
  dt <- (yt[[k]]$y - arct)
  lt <- rowSums(dt^2)
  d <- (y[[k]]$y - arc)
  l <- rowSums(d^2)
  mset[k] <- mean( lt )
  mse[k] <- mean( l )
  
  #compute orthogonaility
  if(pc$type == 3 || pc$type == 0){
    d = d / cbind(sqrt(l), sqrt(l))
    dt = dt / cbind(sqrt(lt), sqrt(lt))
  }
  if(pc$type != 1){
    tl = rowSums( y[[k]]$tangents[[1]]^2)
    tlt = rowSums( yt[[k]]$tangents[[1]]^2)
    y[[k]]$tangents[[1]] = y[[k]]$tangents[[1]] / cbind(sqrt(tl), sqrt(tl))
    yt[[k]]$tangents[[1]] = yt[[k]]$tangents[[1]] / cbind(sqrt(tlt), sqrt(tlt))
  } 
  ortho[k]  = mean( rowSums( y[[k]]$tangents[[1]]  * d  )^2 )
  orthot[k] = mean( rowSums( yt[[k]]$tangents[[1]] * dt )^2 )
 
  #print orthogoanilty and mse values
  print( sprintf("ortho: %f", ortho[k]) ) 
  print( sprintf("mse: %f", mse[k]) )


  #cross-validation 
  if(k>1){
    if(pc$type == 0){
      if(mset[k] > mset[k-1]){
     #   selected = selected+1
      }
    }
    else{
      if(orthot[k] > orthot[k-1]){
      #  selected = selected+1
      }
    }
  }

  #print curve every nIter iteration until selected curve based on cross-validation
  if(selected < 0){
    if(k %% nIter == 0){
      lines(yp$y, col=col, lty=1, lwd=lwp)
    }
  }
  if(selected == 0){
     lines(yp$y, col=col2, lty=1, lwd=lwd)
     #selected=1;
		#break
  }
  
}


#cross-validtion did not select curve - plot curve at end of  optimization
if(selected < 0 ){
  lines(yp$y, col=col2, lty=1, lwd=lwd)
  selected=1;
}

#pretty plot legend
legend("topright", col = c("black", "darkgray", col, col2), lty=c(6, 1, 1, 1), lwd=c(lwg, lwi, lwp, lwd),
legend = c("ground truth", "initialization", "intermediates", "selected"), cex=1.75, bty="n", seg.len=4 )


#dev.off()
#plot test and train error
dev.new()

par(mar=c(5,6,4,2))
if(pc$type==0){
  mset[mset > 5*max(mse)] = 5*max(mse)
  
  plot((1:length(mse)), mse, type="l", lwd=8, cex.lab=1.75, cex.axis=1.75,
       col="darkgray", lty=1, xlab="iteration", ylab=expression(hat(d)(lambda,
       Y)^2), bty="n", ylim = range(mse, mset)) 
  mset = mset -0.005
  lines((1:length(mset)), mset, lwd=8, lty=2, col="black")

  i1 <- which.min(mse)
  i2 <- which.min(mset)
  points(i1, mse[[i1]], col="darkgray", pch=19, cex=3) 
  points(i2, mset[[i2]], col="black", pch=19, cex=3) 
  
  legend("topright", col = c("darkgray", "black"), lty=c(1,2), lwd=c(8,8), legend = c("train", "test"), cex=1.75, bty="n", seg.len=4 )

}
if(pc$type != 0){

  plot((1:length(ortho)), ortho, type="l", lwd=8, cex.lab=1.75, cex.axis=1.75, col="darkgray", lty=1, xlab="iteration", ylab=expression(hat(q)(lambda, Y)^2), bty="n", ylim = range(ortho, orthot))
  lines((1:length(orthot)),orthot, lwd=8, lty=2, col="black")

  i1 <- which.min(ortho)
  i2 <- which.min(orthot)
  points(i1, ortho[[i1]], col="darkgray", pch=19, cex=3) 
  points(i2, orthot[[i2]], col="black", pch=19, cex=3) 

legend("topright", col = c("darkgray", "black"), lty=c(1,2), lwd=c(8,8),
legend = c("train", "test"), cex=1.75, bty="n", seg.len=4 )

}
#
#save output to file
#y=sigmaX[[k-1]]
#zz <- file("sigmaX", "wb")
#for(i in seq(1:length(y)))
#{
#	writeBin(as.double(y[i]), zz)
#}
#close(zz)

#y=sigmaY[[k-1]]
#zz <- file("sigmaY", "wb")
#for(i in seq(1:length(y)))
#{
#	writeBin(as.double(y[i]), zz)
#}
#close(zz)
#
#y=yp$y
#zz <- file("C", "wb")
#for(i in seq(1:dim(y)[1]))
#{
#	writeBin(as.double(y[i,1]), zz)
#	writeBin(as.double(y[i,2]), zz)
#}
#close(zz)


}
  

