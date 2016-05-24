cem.example.sr <- function(n =1000, nstd=0.1, init=0, risk=2, stepX=0.1){

  
#Create data 
d     <- swissroll(N=n, nstd = nstd)
dtest <- swissroll(N=n, nstd = nstd)

#init
if(init == 0){
  print("Computing Isomap inital embedding")
  iso <- isomap(dist(d$Xn), k=15)
  z = iso$points[, 1:2] #d$X[, 2:3]
}else{
  z = d[, 3]
}


ps0 <- cem(y = d$Xn, x = z, knnX=100,  iter=0, verbose=2,
    stepX=1,sigmaX = 0.3, stepBW=0.1, risk=risk)


  ps = ps0

  if(ps$risk==0){
    col = "gold2"  
  }else{
    col = "dodgerblue2"
  }

  px <- ps$x;
  perr = Inf


  #plot(iso$points[, 1:2], pch=19, col="darkgray") 
  plot(px, pch=19, col="darkgray") 


  for(i in 1:10){
    ps <- cem.optimize(ps, iter=1, verbose=2, stepX=stepX, stepBW=0.1,
        nPoints=500, optimalSigmaX=T)
    x <- ps$x

    #segments(ps$z[,1], ps$z[, 2], pz[,1], pz[,2], lwd=2, col="#1C86EE55")
    segments(px[,1], px[, 2], x[,1], x[,2], lwd=2, col="#1C86EE55")

    px = x
  } 
   
  coords <- ps$x
  p <- predict(ps, coords)
     




  #build surface model
  range(x[,1])
  range(x[,2])
  r1 <- range(x[,1])
  r2 <- range(x[,2]) #*0.95
  s1 <- seq(r1[1], r1[2], length.out=20)
  s2 <- seq(r2[1], r2[2], length.out=20)
  #s2 <- s2[4:7]
  l1 <- length(s1)
  l2 <- length(s2)

  gx <- matrix(nrow=l1*l2, ncol=2)
  index = 1;
  for(x1 in s1){
    for(x2 in s2){
      gx[index, 1] = x1
      gx[index, 2] = x2
      index = index +1
    }
  }

  gy <- predict(ps, gx)

  indices = c()
  index = 1
  for(i in 1:(l1-1) ){
    for(j in 1:(l2-1) ){
      indices[index] =   (i-1)*l2+1   +j-1
      indices[index+1] = i*l2+1       +j-1
      indices[index+2] = i*l2+1       +j
      indices[index+3] = (i-1)*l2+1   +j
      index = index + 4
    }
  }

  qm <- qmesh3d(vertices=t(gy$y), indices=indices, homogeneous=F)

  #rgl.open()
  #rgl.bg(color="white")
  #plot3d(d$Xn, risk="s", radius=0.35, box=F, alpha=0.25)
  #material3d(color=col, alpha=0.75, ambient=col, depth_mask=F)
  #shade3d(qm)
  #material3d(col="darkgray", ambient="darkgray", lit=F, alpha=0.5, depth_mask=T)
  #wire3d(qm)



#pgx <- predict(ps, dtest$gYn)
  rgl.open()
  rgl.bg(color="white")
  plot3d(gy$y, type="s", col = col,  radius=0.2, box=F, axes=F, alpha=1, xlab="", ylab="", zlab="")
  

  material3d(lwd=2, alpha=0.75)
  wire3d(d$qm)
  #y <- predict(ps, x)
  spheres3d(p$y, col = "black", radius=0.1)

  spheres3d(d$Xn, col = "gray", radius=0.1)


if(F){
plot3d(d$Xn, type="s", col = "gray", radius=0.5, box=F, axes=F, , alpha=0.5, xlab="", ylab="", zlab="")

material3d(lwd=2, alpha=0.75)
wire3d(d$qm)


col = "darkgray"
alpha=0.5

#qm = qm1
#col = "dodgerblue2"
#alpha=0.9

#qm = qm0
#col = "gold2"
#alpha = 0.3

material3d(color=col, alpha=alpha, ambient=col, lit=T)
shade3d(qm)
material3d(col=col, ambient="black", lit=F, alpha=0.5, lwd=3)
wire3d(qm)



}


}	
