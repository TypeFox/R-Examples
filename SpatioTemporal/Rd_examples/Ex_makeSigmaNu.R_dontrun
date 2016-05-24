##First create some random locations
x <- rnorm(5)
y <- rnorm(5)

##compute distance matrix
D <- crossDist( cbind(x,y) )

#a vector of locations
I <- c(1,2,3,1,4,4,3,2,1,1)
T <- c(1,1,1,2,2,3,3,3,3,4)

##create a block diagonal matrix consisting of four parts with
##exponential covariance.
sigma.nu <- makeSigmaNu(c(.4,2), D, "exp", nugget=0.1,
                        blocks1 = c(3,2,4,1), ind1 = I)
##and cross covariance
sigma.nu.c <- makeSigmaNu(c(.4,2), D, "exp", nugget=0.1,
                          blocks1 = c(3,2,4,1), ind1 = I, 
                          blocks2 = c(0,0,3,1), ind2 = I[7:10])

##compare the cross-covariance with the relevant part of sigma.nu
range(sigma.nu.c-sigma.nu[,7:10])
\dontshow{
  if( max(abs(sigma.nu.c-sigma.nu[,7:10])) > 1e-13 ){
    stop("make.sigma.nu.cross.cov 1: Results not equal")
  }
  sigma.nu.s <- makeSigmaNu(c(.4,2), D, "exp", nugget=0.1,
                            blocks1 = c(3,2,4,1), ind1 = I,
                            sparse=TRUE)
  sigma.nu.c.s <- makeSigmaNu(c(.4,2), D, "exp", nugget=0.1,
                              blocks1 = c(3,2,4,1), ind1 = I, 
                              blocks2 = c(0,0,3,1), ind2 = I[7:10],
                              sparse=TRUE)
  if( max(abs(sigma.nu.s-sigma.nu)) > 1e-13 ){
    stop("make.sigma.nu sparse not equal")
  }
  if( max(abs(sigma.nu.c.s-sigma.nu.c)) > 1e-13 ){
    stop("make.sigma.nu sparse/cross not equal")
  }
}
##an alternative showing the use of loc.ind2.to.1
sigma.nu.c <- makeSigmaNu(c(.4,2), D[,4:3], "exp", nugget=0.1,
                          blocks1 = c(3,2,4,1), ind1 = I, 
                          blocks2 = c(0,0,2,0), ind2 = 1:2,
                          ind2.to.1=4:3)
##compare the cross-covariance with the relevant part of sigma.nu
range(sigma.nu.c-sigma.nu[,6:7])
\dontshow{
  if( max(abs(sigma.nu.c-sigma.nu[,6:7])) > 1e-13 ){
    stop("make.sigma.nu.cross.cov 2: Results not equal")
  }
}
