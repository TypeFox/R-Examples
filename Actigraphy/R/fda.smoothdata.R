fda.smoothdata <-
function(data, basistype="fourier", nbasis=9, norder=4){
if(missing(data)) 
stop("Missing data")

mat <- data$mat
cov <- data$cov
L <- nrow(mat)

if(tolower(basistype) == "fourier"){
fbase <- create.fourier.basis(rangeval=c(0,L), nbasis)
}else if(tolower(basistype) == "bspline"){
fbase <- create.bspline.basis(rangeval=c(0,L), nbasis, norder)
}else{
stop("basistype must be 'fourier' or 'bspline'.")
}

fpar <- fdPar(fbase) 
fd <- smooth.basis(c(1:L), mat, fpar)
FD <- list(fd=fd, cov=cov)
return(FD)
}
