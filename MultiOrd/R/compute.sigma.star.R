compute.sigma.star <-
function (prop.vec.bin, corr.mat) 
{
    no.bin= length(prop.vec.bin)
 
    sigma = corr.mat
    p = prop.vec.bin
    q = 1 - p

sigmaBB = diag(no.bin)
for (i in 1:no.bin) {
for (j in 1:no.bin) {
if (i != j) {
  sigmaBB[i, j] = phi2poly(sigma[i, j], p[i], p[j]) }
}
}
if(!is.positive.definite(sigmaBB)){
warning( "Tetrachoric correlation matrix is not positive definite)")
sigmaBB=as.matrix(nearPD(sigmaBB, corr = TRUE, 
            keepDiag = TRUE)$mat)
}

sigmaBB = ( sigmaBB+t(sigmaBB) )/2
return(sigmaBB )
}
