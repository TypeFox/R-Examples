sBF <-
function(
    dat, #matrix of data
    depCol = 1, # position of dependent variable
    m = 100, # grid points
    windows = rep(20, ncol(dat)-1), # number of windows. (covariate range)/windows provide the bandwidth for the kernel regression smoother
    bw = NULL, # bandwidth
    method = "gaussian", #kernel method
    mx = 100,  # maximum iterations number
    epsilon = 0.0001, # convergence limit
    PP = NULL, # matrix of joint probabilities
    G = NULL # grid
) {
n <- nrow(dat)
d <- ncol(dat)-1
Y <- dat[,depCol]
c <- mean(Y)
Ybar <- Y-c
X <- as.matrix(dat[, 1:(d+1)!=depCol])
genG <- FALSE
if (is.null(G)) {
	G <- matrix(, m, d)
	genG <- TRUE
} else {
G <- as.matrix(G)
m <- nrow(G)
}
P <- matrix(, m, d)
h <- h2 <- matrix(,d, 1)
if (is.null(h2)) h2 <-  matrix(h2,d, 1)
if (length(windows) != d) windows <- rep(windows[1], d)
for (j in 1:d) {           
  minima <- min(X[,j])
  maxima <- max(X[,j])
  h2[j] <- (maxima-minima)/windows[j]
  if (genG) G[,j] <- seq(from = minima, to = maxima, length = m)
  densit <- density(X[,j], from = minima, to = maxima, n = m)
  P[,j] <- densit$y
  h[j] <- densit$bw
}
if (!is.null(bw)) h2 <- matrix(bw, d, 1) 
Pp <- pP <- matrix(, m^2, choose(d,2))
if(d > 1) {
if (is.null(PP)) {
PP <- matrix(, m^2, choose(d,2))
i <- 1:m %x% rep(1, m)
v <- rep(1:m, m)
j <- rep(1:(d-1), (d-1):1)
k <- matrix(rep(2:d, d-1), d-1, d-1)[lower.tri(matrix(, d-1, d-1), diag=TRUE)]
prJ <- function(z, s) (1/h[j[s]])*(1/h[k[s]])* mean(K((G[i[z],j[s]]-X[,j[s]])/h[j[s]], method = method)*K((G[v[z],k[s]]-X[,k[s]] )/h[k[s]], method = method))
print("Processing joint distributions (matrix PP)", quote=FALSE)
pb <- txtProgressBar(style=3)
for (z in 1:m^2) {
  for (s in 1:choose(d,2))
  	PP[z,s] <- prJ(z, s)
  perc <- z/m^2
  setTxtProgressBar(pb, perc)
}
close(pb)
}
Pr <- NULL
for (j in 1:(d-1)) Pr <- c(Pr,rep(P[,j],d-j))
Pr <- matrix(Pr,m,choose(d,2))
D <- Pr %x% rep(1,m)
pP <- PP/D
Pr <- NULL
for (j in 2:d) Pr <- c(Pr,rep(P[,j],j-1))
Pr <- matrix(Pr,m,choose(d,2))
D  <-  rep(1,m) %x% Pr
u <- matrix(1:(d-1), 1, d-1)
for(rip in 1:(d-2)) u <- rbind(u, u[rip,]+d-rip-1)
perm <- u[upper.tri(u,diag=TRUE)]
PPperm <- PP[,perm]
Pp <- PPperm/D
PpPermRig <-order(rep(1:m, m))
Pp <- as.matrix(Pp[PpPermRig,])    # offset problem when d=2
}
M  <- matrix(, m, d)
for (j in 1:d) {
  for(i in 1:m) {
      M[i,j]  <- 1/h2[j]*(mean(K(((G[i,j]-X[,j])/h2[j]), method = method)*Y)/P[i,j])
  }
}
m0 <- matrix(,d,1)
for (j in 1:d) {
m0[j] <- sum(M[,j]*P[,j])/sum(P[,j])
}
Mcurr <- M
I     <- t(diag(1,m,m) %x% matrix(1,m,1))
e1    <- matrix(1,m,1)
for (r in 1:mx) {
  Mpre <- Mcurr
  for(j in 1:d) {
    A <- B <- matrix(0, m, 1)
    minA <- d*(j-1)-(j-1)*j/2+1
    maxA <- j*d - j*(j+1)/2 
    if (j < d) A <- I %*% (pP[,minA:maxA]*(e1 %x% Mpre[,(j+1):d]) * ((I %*% pP[,minA:maxA]) %x% e1)^(-1) )
    minB <- (j-1)*(j-2)/2+1
    maxB <- j*(j-1)/2
    if (j > 1) B <- I %*% (Pp[,minB:maxB]*(e1 %x% Mcurr[,1:(j-1)]) * ((I %*% Pp[,minB:maxB]) %x% e1)^(-1))
    Mcurr[,j] <- M[,j]  - c -rowSums(cbind(A,B))
    Mcurr[,j] <- Mcurr[,j] - sum(Mcurr[,j]*P[,j])/sum(P[,j])
  }
  if(prod((colSums((Mcurr-Mpre)^2)/(colSums(Mpre^2)+epsilon))<epsilon)==1) break
}
return(
       list(
       mxhat = Mcurr, # estimated univariate functions 
       m0 = c,
       grid = G, # grid
       conv = (colSums((Mcurr-Mpre)^2)/(colSums(Mpre^2)+epsilon))<epsilon, # obtained convergences
       nit = r, # number of iterations
       PP = PP, # matrix of joint probabilities
       bw = h2 # bandwidth
       )
      )
}

