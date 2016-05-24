
cem.example.arc <- function(n=150, noise=0.2, risk=2, sigmaX= 0.1, stepX=0.001,
    stepBW=0.01, init = 0, plotEach = 1, noiseInit = 0.5){


### Create arc data set

#create arc
phi <- runif(n)*pi 
arc <- cbind(cos(phi), sin(phi)) * (1+rnorm(n) * noise)
arc0 <- cbind(cos(phi), sin(phi))

#phi <- runif(n)*pi 
#arc <- cbind(phi, sin(phi*phi*phi)/(phi*phi)) * (1+rnorm(n) * noise)
#arc0 <- cbind(phi, sin(phi*phi*phi)/(phi*phi))


#order by angle for plotting curves    
o0 <- order(phi) 

#test data
phit <- runif(n*10)*pi 
arct <- cbind(cos(phit), sin(phit)) * (1+rnorm(10*n) * noise)

#phit <- runif(n*10)*pi 
#arct <- cbind(phit, sin(phit*phit*phit)/(phit*phit)) * (1+rnorm(10*n) * noise)





### Initialization

if(init == 0){
 #ground truth initialization
 z <- phi
}else if(init == 1){
 #random initialization
 z <- runif(n)
}
else if(init == 2){
 #random initialization
 z <- phi+runif(n) * noiseInit
}
else{
 #close to principal component initialization
 z = arc[,2]
}

#compute initial principal curve
#do not optimize sigmaX after each itertion to show optimization path
pc <-  cem(y=arc, x=z, knnX=nrow(arc)/10, iter=0, optimalSigmaX=T)

#set risk
pc$risk=risk

#set intial bandwidth
pc$sigmaX=sigmaX; #smoothness of initial curve 


#initial prediction (curve)
xi <- pc$x
yi <- predict(pc, xi)
oi <- order(xi)

#create sample from initial curve for plotting
r = range(xi)
xpi = seq(r[1], r[2], length.out=500)
ypi = predict(pc, xpi)


if( pc$risk > 0 ){
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


#cross validation flag (selected curve)
selected = -1


#run a 100 iterations, one at a time 
#for running the whole optimization in one go run either
# pc <- cem(y=arc, z=z, knnX=n, knnY=50, iter=100)
# pc <- cem.optimize(pc, stepX=1, stepBW=0.1, iter=100)
# here one iterations is run at the time and corss-validation is performed on a
# test set. For minimizing orthogonality cross-validation appears to be not
# necessary.

for(k in 1:100){
  #run one iterations  
  pc <- cem.optimize(pc, stepX=stepX, stepBW=stepBW, iter=1, verbose=2,
      optimalSigmaX=T, nPoints=100)
  
  #store cem values at iteration k
  
  sigmaX[k] <- pc$sigmaX
  

  x[[k]] <- pc$x
  
  #train data
  r = range(x[[k]])
  xp = seq(r[1], r[2], length.out=500)
  yp = predict(pc, xp);
  o[[k]] <- order(x[[k]]);
  y[[k]] <- predict(pc, x[[k]])

  
  #compute mean squared errors
  d <- (y[[k]]$y - arc)
  l <- rowSums(d^2)
  mse[k] <- mean( l )
  
  #compute orthogonaility
  if(pc$risk == 3 || pc$risk == 0){
    d = d / cbind(sqrt(l), sqrt(l))
  }
  if(pc$risk != 1){
    tl = rowSums( y[[k]]$tangents[[1]]^2)
    y[[k]]$tangents[[1]] = y[[k]]$tangents[[1]] / cbind(sqrt(tl), sqrt(tl))
  } 
  ortho[k]  = mean( rowSums( y[[k]]$tangents[[1]]  * d  )^2 )
 
  #print orthogoanilty and mse values
  print( sprintf("ortho: %f", ortho[k]) ) 
  print( sprintf("mse: %f", mse[k]) )



  #print curve every nIter iteration until selected curve based on cross-validation
#if(selected < 0){
    if(k %% nIter == 0){
      lines(yp$y, col=col, lty=1, lwd=lwp)
    }
# }
  if(selected == 0){
     lines(yp$y, col=col2, lty=1, lwd=lwd)
     selected=1;
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



#plot test and train error
dev.new()

par(mar=c(5,6,4,2))
if(pc$risk==0){
  
  plot((1:length(mse)), mse, type="l", lwd=8, cex.lab=1.75, cex.axis=1.75,
       col="darkgray", lty=1, xlab="iteration", ylab=expression(hat(d)(lambda,
       Y)^2), bty="n", ylim = range(mse)) 

  i1 <- which.min(mse)
  points(i1, mse[[i1]], col="darkgray", pch=19, cex=3) 
  
  legend("topright", col = c("darkgray", "black"), lty=c(1,2), lwd=c(8,8), legend = c("train", "test"), cex=1.75, bty="n", seg.len=4 )

}
if(pc$risk != 0){

  plot((1:length(ortho)), ortho, type="l", lwd=8, cex.lab=1.75, cex.axis=1.75, col="darkgray", lty=1, xlab="iteration", ylab=expression(hat(q)(lambda, Y)^2), bty="n", ylim = range(ortho))

  i1 <- which.min(ortho)
  points(i1, ortho[[i1]], col="darkgray", pch=19, cex=3) 

legend("topright", col = c("darkgray", "black"), lty=c(1,2), lwd=c(8,8),
legend = c("train", "test"), cex=1.75, bty="n", seg.len=4 )

}





}
  

