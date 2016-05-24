#' Return the Toy Signal.
#'
#' @param ImRange A vector with two components giving the range of the region on
#'               which the Toy Signal is to be computed.
#' @param NPixel Number of pixels of the result in one direction. The resulting
#'               picture will have NPixel x NPixel pixels.           
#' @return A list with components "x", "y" and "z". Here, x and y are the 
#'         coordinates of the grid and z is matrix of dimensions 
#'         c(NPixel,NPixel) giving the Toy Signal.
#' @export
ToySignal = function(ImRange = c(0,10), NPixel = 64){
  
  s = seq(0, 10, length.out = NPixel)
  ds = s[2] - s[1]
  
  #Create single peak test signal.
  single.peak = function(sx,sy,x, y, b) {
    matrix(mvtnorm::dmvnorm(expand.grid(sx,sy),mean=c(x,y),sigma=diag(c(b,b))),
           nrow=length(sy))
  }
  
  #Defining the signal.
  mu1 = 50*single.peak(s, s, 2.6, 5.7, 1.5) + 
    100*single.peak(s, s, 6.8, 7.8, 2.5) +
    50*single.peak(s, s, 7.8, 2.1, 1.3) 
  mu = mu1/max(mu1) * 3
  
  s = seq(from = ImRange[1], to = ImRange[2], length.out = NPixel)
  
  list(x=s,y=s,z=mu)
} 


#' Generate a realization of the Toy Noise 1.
#'
#'@param n The number of realizations to produce.
#' @return A list containing x and y, the coordinates of the grid and
#'        z and array of dimensions c(64,64,n) giving n reallizations of the 
#'        Toy Noise 1.
#' @export
ToyNoise1 = function(n = 1){
  #Parameters of the domain.
  #Pixel size.
  Ns = 64
  #Grid coordinates
  s = seq(0,10,length.out=Ns)
  ds = s[2]-s[1]
  theta = 1
  l1 = 4
  l2 = 1
  tau1 = 25
  z = array(dim=c(64,64,n))
  for(i in 1:n){
    Z = cbind(matrix(rnorm(Ns^2/2/l1^2),Ns/l1,Ns/2/l1)%x% matrix(1,l1,l1),
              matrix(rnorm(Ns^2/2/l2^2),Ns/l2,Ns/2/l2)%x% matrix(1,l2,l2))
    z[,,i] = tau1 * fields::image.smooth(Z,theta=theta,dx=ds,dy=ds)$z
  }
  
  if(n == 1) z = matrix(z,64,64)
  list(x = s, y = s, z = z)
}

#' Generate a realization of the Toy Noise 1 before smoothing.
#'
#'@param n The number of realizations to produce.
#' @return A list containing x and y, the coordinates of the grid and
#'        z and array of dimensions c(64,64,n) giving n reallizations of the 
#'        Toy Noise 1 before smoothing.
#' @export
ToyNoise1Presmooth =function(n = 1){
  #Parameters of the domain.
  #Pixel size.
  Ns = 64
  #Grid coordinates
  s = seq(0,10,length.out=Ns)
  ds = s[2]-s[1]
  l1 = 4
  l2 = 1
  z = array(dim=c(64,64,n))
  for(i in 1:n){
    z[,,i] = cbind(matrix(rnorm(Ns^2/2/l1^2),Ns/l1,Ns/2/l1)%x% matrix(1,l1,l1),
              matrix(rnorm(Ns^2/2/l2^2),Ns/l2,Ns/2/l2)%x% matrix(1,l2,l2))
  }
  
  if(n == 1) z = matrix(z,64,64)
  list(x = s, y = s, z = z)
}


#' Generate a realization of the Toy Noise 2.
#'
#'@param n The number of realizations to produce.
#' @return A list containing x and y, the coordinates of the grid and
#'        z and array of dimensions c(64,64,n) giving n reallizations of the 
#'        Toy Noise 2.
#' @export
ToyNoise2 = function(n = 1){
  #Parameters of the domain.
  #Pixel size.
  Ns = 64
  #Grid coordinates
  s = seq(0,10,length.out=Ns)
  ds = s[2]-s[1]
  theta = 1
  l1 = 4
  l2 = 1
  tau2 = 100
  laplaceker = function(x) fields::double.exp(sqrt(x))
  z = array(dim=c(64,64,n))
  for(i in 1:n){
    Z = cbind(matrix(rnorm(Ns^2/2/l1^2),Ns/l1,Ns/2/l1)%x% matrix(1,l1,l1),
              matrix(rnorm(Ns^2/2/l2^2),Ns/l2,Ns/2/l2)%x% matrix(1,l2,l2))
    z[,,i] = tau2 * fields::image.smooth(Z,theta=theta,dx=ds,dy=ds,
                             kernel.function = laplaceker)$z
  }
  
  if(n == 1) z = matrix(z,64,64)
  list(x = s, y = s, z = z)
}


#' Generate a realization of the Toy Noise 3.
#'
#' @param n The number of realizations to produce.
#' @return A list containing x and y, the coordinates of the grid and
#'        z and array of dimensions c(64,64,n) giving n reallizations of the 
#'        Toy Noise 3.
#' @export
ToyNoise3 = function(n = 1){
  #Parameters of the domain.
  #Pixel size.
  Ns = 64
  #Grid coordinates
  s = seq(0,10,length.out=Ns)
  ds = s[2]-s[1]
  theta = 1
  l1 = 1
  l2 = 1
  tau3 = 50
  rlaplace = function(n,b=1){
    (2*rbinom(n,size=1,prob=0.5)-1)*rexp(n,1/b)
  }
  z = array(dim=c(64,64,n))
  for(i in 1:n){
    Z = cbind(matrix(rlaplace(Ns^2/2/l1^2,b=1),Ns/l1,Ns/2/l1)%x% matrix(1,l1,l1),
              matrix(rt(Ns^2/2/l2^2,df=10),Ns/l2,Ns/2/l2)%x% matrix(1,l2,l2))
    z[,,i] = tau3 * fields::image.smooth(Z,theta=theta,dx=ds,dy=ds)$z
  }
  
  if(n == 1) z = matrix(z,64,64)
  list(x = s, y = s, z = z)
}
