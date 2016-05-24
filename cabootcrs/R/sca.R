
sca <- function(X) { 

# Performs simple CA of a matrix X with rows >= cols

# Output:
# R = row profile matrix
# C = column profile matrix
# Rweights = D_c^-1/2
# Cweights = D_r^-1/2
# Raxes = V
# Caxes = U
# r = rank
# mu = mu

X <- X/sum(X)

rmax <- min(dim(X))-1

rsums <- as.vector(rowSums(X))
csums <- as.vector(colSums(X))

drm1 <- diag( 1/( rsums + (rsums==0) ) * (1-(rsums==0)) )
dcm1 <- diag( 1/( csums + (csums==0) ) * (1-(csums==0)) )
drmh <- sqrt(drm1)
dcmh <- sqrt(dcm1)

Z <- drmh %*% ( X - rsums %*% t(csums) ) %*% dcmh

Y <- svd(Z)

mu <- Y$d[1:rmax]
r <- sum(mu>1e-15)

R <- drm1 %*% X
Rweights <- dcmh
Raxes = Y$v[ ,1:rmax]

C <- dcm1 %*% t(X);
Cweights <- drmh;
Caxes = Y$u[, 1:rmax]

if (r < rmax) { 
  mu[ (r+1):rmax ] <- 0;
  Raxes[ , (r+1):rmax ] <- 0;
  Caxes[ , (r+1):rmax ] <- 0;
}

sca <- new("cabasicresults",
          Rprofile=R,Cprofile=C,Rweights=Rweights,Cweights=Cweights,
          Raxes=Raxes,Caxes=Caxes,r=r,mu=mu)

}


