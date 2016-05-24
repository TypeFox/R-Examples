#' Cube preset dataset
#'
#' @title gng.preset.cube
#' @description Generate sample cube dataset
#' @export

#' @param N Number of points
#' @param r Length of the side of cube
#' @param center Center of the plane
#' 
#' @examples
#' X <- gng.preset.cube(100)
#' gng <- GNG(X)
#' 
gng.preset.cube<-function(N, r=0.5, center=c(0.5,0.5,0.5)){
  .gng.box_point<-function(r, center, prob=-1){ 
    point <- c()
    
    if(prob == -1)
      point<-center
    else
      point<-c(center, prob)
    
    point[1:3] = point[1:3] + runif(3, min=-r/2.0, max=r/2.0)   
    
    point
}


  mat<-matrix(0,N,3)
  
  for(i in 1:N){
    mat[i,] = .gng.box_point(r=r, center=center)
  }
  
  mat
}


#' Plane preset dataset
#'
#' @title gng.preset.plane
#' @description Generate sample plane dataset
#' @export

#' @param N Number of points
#' @param side Length of the side of plane
#' @param center Center of the plane
#' 
#' @examples
#' X <- gng.preset.plane(100)
#' gng <- GNG(X)
#' 
gng.preset.plane<-function(N, side=0.5, center=c(0.5,0.5,0.5)){
.gng.plane_point<-function(r,center){
  if(!hasArg(r)) r<-1.0
  if(!hasArg(center)) center<-c(0,0,0)
  
  point<-center
  
  point[1]<-point[1]+r*runif(1.0)
  point[2]<-point[2]+r*runif(1.0)
  point[3]<-point[3]
  
  return(point)
}
  mat<-matrix(0,N,3)
  
  for(i in 1:N){
    mat[i,] = .gng.plane_point(side, center)
    mat[i,3] = mat[i,1]
  }
  
  mat
}


#' Sphere preset dataset
#'
#' @title gng.preset.sphere
#' @description Generate sample sphere dataset
#' @export

#' @param N Number of points
#' @param r Radius of the sphere
#' @param center Center of the sphere
#' 
#' @examples
#' X <- gng.preset.sphere(100)
#' gng <- GNG(X)
gng.preset.sphere<-function(N, r=0.5, center=c(0.5,0.5,0.5)){
.gng.sphere_point<-function(r,center){
  if(!hasArg(r)) r<-1.0
  if(!hasArg(center)) center<-c(0,0,0)
  
  alpha<-runif(1)*2*pi
  beta<-runif(1)*pi
  
  point<-center
  
  point[1]<-point[1]+r*cos(alpha)*sin(beta)
  point[2]<-point[2]+r*sin(alpha)*sin(beta)
  point[3]<-point[3]+r*cos(beta)
  
  return(point)
}



  mat<-matrix(0,N,3)
  
  for(i in 1:N){
    mat[i,] = .gng.sphere_point(r, center)
  }
  
  mat
}


.sigmoid <- function(x){
   1./(1.+exp(-x))
}

gng.preset_potential<-function(N, r=0.5, center=c(0.5,0.5,0.5), prob=-1){
  
  mat <- c()
  if(prob == -1){
    mat <- matrix(rnorm(20,mean=1), N,3)
  }
  else{
    mat <- matrix(rnorm(20,mean=1), N,4)
  }
  
  
  for(j in 1:N){
    t<-rnorm(1,mean=0,sd=1)
    u<-rnorm(1,mean=0,sd=1)
    val<-.sigmoid(t^2+u^2);
    mat[j,1] = t
    mat[j,2] = u
    mat[j,3] = val 	
    if(prob!=-1) mat[j,4] = prob
    
  }
  
  mat
}







