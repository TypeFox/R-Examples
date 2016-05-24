swissroll <- function(N=1000, phi = pi, height=4, radius=1, nstd = 0.1, mphi=50, mh =10){
r=radius
h = height
tt = phi*(1+2*runif(N))  
height = height*runif(N)
X = rbind(r*tt * cos(tt), height, r*tt * sin(tt));

Xa = rbind(r*cos(tt)- r*tt * sin(tt), matrix(0, 1, N), r*sin(tt) + r*tt * cos(tt))
Xo = rbind(-Xa[3,], 0, Xa[1, ]);
Xo = Xo / matrix(sqrt(colSums(Xo^2, 1)), ncol=N, nrow=3);
Xn = X;
Xn = X + Xo* t(matrix(nstd*rnorm(N), N, 3));

rownames(X) <- c("x", "y", "z")
rownames(Xn) <- c("x", "y", "z")


s1 <- seq(phi, 3*phi, length.out=mphi)
s2 <- seq(0, h, length.out=mh)
#s2 <- s2[4:7]
l1 <- length(s1)
l2 <- length(s2)

gy <- matrix(nrow=l1*l2, ncol=3)
gyn <- matrix(nrow=l1*l2, ncol=3)

index = 1;
for(x1 in s1){
  for(x2 in s2){
    gy[index, 1] = r*cos(x1)*x1
    gy[index, 2] = x2
    gy[index, 3] = r*sin(x1)*x1

    xa = c(r*cos(x1)- r* x1 * sin(x1), 1, r*sin(x1) + r* x1 * cos(x1))
    xo = c(-xa[3], 0, xa[1]);
    xo = xo / sqrt(sum(xo^2));
    gyn[index, ] = gy[index,] + xo*rnorm(1, 0, nstd)

    index = index +1
  }
}


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

qm <- qmesh3d(vertices=t(gy), indices=indices, homogeneous=FALSE)



obj <- structure(list( Xn = t(Xn), X = t(X), qm = qm, gY = gy, gYn = gyn ) )

obj
}
