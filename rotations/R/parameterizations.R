#' SO3 class.
#'
#' Class for \eqn{3\times 3}{3-by-3} matrix representation of rotations.
#'
#' @name SO3-class
#' @seealso See the \code{\link{SO3}} functions.
#'
#' @exportClass SO3
setOldClass("SO3")


#' Q4 class.
#'
#' Class for quaterion representation of rotations.
#' 
#' @name Q4-class
#' @seealso See the \code{\link{Q4}} functions.
#'
#' @exportClass Q4
setOldClass("Q4")

#' Quaternions
#' 
#' Creates or tests for objects of class "Q4".
#' 
#' Construct a single or sample of rotations in 3-dimensions in quaternion form.  Several possible inputs for \code{x}
#' are possible and they are differentiated based on their class and dimension.
#' 
#' For \code{x} an n-by-3 matrix or a vector of length 3, the angle-axis representation of rotations is utilized.  More specifically,
#' each quaternion can be interpreted as a rotation of some reference frame about the axis 
#' \eqn{U} (of unit length) through the angle \eqn{\theta}.  For each axis and angle the quaternion is formed through
#' \deqn{q=[cos(\theta/2),sin(\theta/2)U]^\top.}{q=[cos(theta/2),sin(theta/2)U]'.}  The object \code{x} is treated as 
#' if it has rows \eqn{U} and \code{theta} is a vector or angles. If no angle is supplied then the 
#' length of each axis is taken to be the angle of rotation theta.
#' 
#' For \code{x} an n-by-9 matrix of rotation matrices or an object of class \code{"SO3"}, this function will
#' return the quaternion equivalent of \code{x}.  See \code{\link{SO3}} or the vignette "rotations-intro"
#' for more details on rotation matrices.  
#' 
#' For \code{x} an n-by-4 matrix, rows are treated as quaternions; rows that aren't of unit length
#' are made unit length while the rest are returned untouched.  A message is printed if any of the rows are not quaternions.
#' 
#' \code{x} a \code{"data.frame"} it is translated into a matrix of the same dimension and
#' the dimensionality of \code{x} is used to determine the data type: angle-axis, quaternion or rotation (see above).
#' As demonstrated below, \code{is.Q4} may return \code{TRUE} for a data frame, but the functions defined for objects of class 
#' \code{'Q4'} will not be called until \code{as.Q4} has been used.
#' 
#'
#' @export
#' @rdname Q4
#' @param x object to be coerced or tested
#' @param theta vector or single rotation angle; if \code{length(theta)==1}, the same theta is used for all axes
#' @param ... additional arguments.
#' @format \code{id.Q4} is the identity rotation given by the matrix \eqn{[1,0,0,0]^\top}{[1,0,0,0]'}.
#' @return 	\item{as.Q4}{coerces its object into a Q4 type} 
#' 					\item{is.Q4}{returns \code{TRUE} or \code{False} depending on whether its argument satisfies the conditions to be an
#' 					quaternion; namely it must be four-dimensional and of unit length}
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @examples
#' data(drill)                    #Pull off subject 1's wrist measurements
#' Subj1Wrist <- subset(drill, Subject == '1' & Joint == 'Wrist') 
#' 
#'                                ## The measurements are in columns 5:8
#' all(is.Q4(Subj1Wrist[,5:8]))   #TRUE, even though Qs is a data.frame, the rows satisfy the 
#'                                #conditions necessary to be quaternions BUT, 
#'                                #S3 methods (e.g. 'mean' or 'plot') for objects of class 
#'                                #'Q4' will not work until 'as.Q4' is used
#'                 
#' Qs <- as.Q4(Subj1Wrist[,5:8])  #Coerce measurements into 'Q4' type using as.Q4.data.frame
#' all(is.Q4(Qs))                 #TRUE  
#' mean(Qs)                       #Estimate central orientation for subject 1's wrist, see ?mean.Q4
#' Rs <- as.SO3(Qs)               #Coerse a 'Q4' object into rotation matrix format, see ?as.SO3
#' 
#' #Visualize the measuremenets, see ?plot.Q4 for more
#' \dontrun{
#' plot(Qs, col = c(1, 2, 3))} 


as.Q4<-function(x,...){
  UseMethod("as.Q4")
}

#' @rdname Q4
#' @method as.Q4 default
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @export

as.Q4.default <- function(x,theta=NULL,...){  
  
  p<-ncol(x)
  n<-nrow(x)
  
  if(is.null(p)){
    p<-length(x)
    #q<-matrix(q,ncol=3)
    #p<-ncol(q)
    #n<-nrow(q)
  }
  
  if(p==3){
    
    #If input is length 3, x is assumed to be vectors in R^3
    U<-x
    U<-matrix(U,ncol=3)

    ulen<-sqrt(rowSums(U^2))
  
    if(is.null(theta)){ 
      theta<-ulen%%pi
    }
  
    ntheta<-length(theta)
  
    if(nrow(U)!=ntheta){
      if(ntheta==1){
        theta<-rep(theta,n)
      }else{
        stop("Number of angles must match number of axes")
      }
    }
    nonZ<-which(ulen!=0)
    
    U[nonZ,]<-U[nonZ,]/ulen[nonZ]
    #CPP version is causing seg faults, try just doing it in R
    #x <- Q4defaultC(U,theta)
  
    q <- cbind(cos(theta/2), sin(theta/2) * U)
    
  }else if(p==4){
    
    #If input has length divisible by 4, data are normalized and made into class "Q4"
    q<-matrix(x,ncol=4)
    rowLens<-(rowSums(q^2))^0.5
    
    notRots<-which(abs(rowLens-1)>10e-8)
    
    if(length(notRots)>0){
      txt<-"Row(s)"
      for(i in notRots){
        txt<-paste(txt,i)
      }
      txt<-paste(txt,"was(were) not quaternions.")
      message(txt)
    }
    
    q<-q/rowLens
    
  }else if(p==9){
    
    #If input has 9 columns, q is assumed to be rotations
    q<-as.Q4.SO3(x)
    
  }else{
    
    stop("Unknown data type.  Pease see ?as.Q4 for more details.")
    
  }
  class(q)<-"Q4"
  return(q)
}

#' @rdname Q4
#' @method as.Q4 SO3
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @export

as.Q4.SO3 <- function(x,...) {
  
  R<-x
  R<-formatSO3(R)
  theta <- mis.angle(R)
  u <- mis.axis(R)
  x <- as.Q4.default(u,theta)
  
  return(x)
}

#' @rdname Q4
#' @method as.Q4 Q4
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @export

as.Q4.Q4 <- function(x,...) {
  
  return(x)
}

#' @rdname Q4
#' @method as.Q4 data.frame
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @export

as.Q4.data.frame <- function(x,...) {
  #n<-nrow(x)
  #p<-ncol(x)
  q<-data.matrix(x)
  q<-matrix(q,ncol=4)
  return(as.Q4.default(q))
}

#' @rdname Q4
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @export

is.Q4 <- function(x) {
  
  if(length(x)==9){
    FALSE
  }else{
  
	  apply(x,1,function(q){sum(q^2)-1<10e-10 & length(q)==4})
	
  }
    
}

#' @rdname Q4
#' @aliases Q4 is.Q4 id.Q4 as.Q4.default as.Q4.SO3 as.Q4.Q4 as.Q4.data.frame
#' @export

id.Q4 <- as.Q4(c(1,0,0,0))


#' Rotation matrices
#' 
#' Creates or tests for objects of class "SO3".
#' 
#' Construct a single or sample of rotations in 3-dimensions in 3-by-3 matrix form.  Several possible inputs for \code{x}
#' are possible and they are differentiated based on their class and dimension.
#' 
#' For \code{x} an n-by-3 matrix or a vector of length 3, the angle-axis representation of rotations is utilized.  More specifically,
#' each rotation matrix can be interpreted as a rotation of some reference frame about the axis \eqn{U} (of unit length)
#' through the angle \eqn{\theta}.  If a single axis (in matrix or vector format) or matrix of axes are provided for \code{x}, 
#' then for each axis and angle the matrix is formed through
#' \deqn{R=\exp[\Phi(U\theta)]}{R=exp[\Phi(U\theta)]} where \eqn{U} is replace by \code{x}.  If axes are provided but \code{theta}
#'  is not provided then the length of each axis is taken to be the angle of rotation, theta.  
#' 
#' For \code{x} an n-by-4 matrix of quaternions or an object of class \code{"Q4"}, this function will
#' return the rotation matrix equivalent of \code{x}.  See \code{\link{Q4}} or the vignette "rotations-intro"
#' for more details on quaternions.  
#' 
#' For \code{x} an n-by-9 matrix, rows are treated as 3-by-3 matrices; rows that don't form matrices in SO(3)
#' are projected into SO(3) and those that are already in SO(3) are returned untouched.  See \code{\link{project.SO3}} for
#' more on projecting arbitrary matrices into SO(3).  A message is printed if any of the rows are not proper rotations.
#' 
#' \code{x} a \code{"data.frame"} is translated into a matrix of the same dimension and
#' the dimensionality of \code{x} is used to determine the data type: angle-axis, quaternion or rotation.
#' As demonstrated below, \code{is.SO3} may return \code{TRUE} for a data frame, but the functions defined for objects of class 
#' \code{"SO3"} will not be called until \code{as.SO3} has been used. 
#' 
#'
#' @export
#' @rdname SO3
#' @param x object to be coerced or tested; see details for possible forms
#' @param theta vector or single rotation angle; if \code{length(theta)==1} the same theta is used for all axes
#' @param ... additional arguments.
#' @format \code{id.SO3} is the identity rotation given by the the 3-by-3 identity matrix.
#' @return 	\item{as.SO3}{coerces provided data into an SO3 type.} 
#' 					\item{is.SO3}{returns \code{TRUE} or \code{False} depending on whether its argument satisfies the conditions to be an
#' 					rotation matrix.  Namely, has determinant one and its transpose is its inverse.}
#' @aliases SO3 as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @examples
#' data(nickel)                   #Select one location to focus on
#' Loc698 <- subset(nickel, location == 698)
#' 
#' is.SO3(Loc698[,5:13])          #Some of the rows are not rotations due to rounding or entry errors  
#'                                #as.SO3 will project matrices not in SO(3) to SO(3)
#' 
#' Rs <- as.SO3(Loc698[,5:13])    #Translate the Rs data.frame into an object of class 'SO3'
#'                                #Rows 4, 6 and 13 are not in SO(3) so they are projected to SO(3)
#'                                
#' mean(Rs)                       #Estimate the central orientation with the average
#' median(Rs)                     #Re-estimate central orientation robustly
#' Qs <- as.Q4(Rs)                #Coerse into "SO3" format, see ?as.SO3 for more
#' 
#' #Visualize the location, there appears to be two groups
#' \dontrun{
#' plot(Rs, col = c(1, 2, 3))}   


as.SO3 <- function(x,...){
  UseMethod("as.SO3")
}

#' @rdname SO3
#' @method as.SO3 default
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

as.SO3.default <- function(x, theta=NULL,...) {

  p<-ncol(x)
  n<-nrow(x)
  
  if(is.null(p)){
    p<-length(x)
    n<-1
  }  
  
  if(n==3 && p==3 && is.SO3(x) && is.null(theta)){
    
    #If there are 3 rows and columns and the object is already a rotation matrix, the same rotation is returned
    class(x) <- "SO3"
    return(x)
    
  }else if(p==9){
    
    rots<-is.SO3(x)
    notRots<-which(!rots)
    
    if(all(rots)){
      #If there are 9 columns and the data are already rotation matrices then the SO3 class is appeneded and object returned
      
      class(x) <- "SO3"
      return(x)
      
    }else{
      #If there are 9 columns and some are not rotations,
      #those that aren't rotations are projected to SO(3) and the others are left alone

      for(i in notRots){
        x[i,]<-project.SO3(x[i,])
      }
      
      txt<-"Row(s)"
      for(i in notRots){
        txt<-paste(txt,i)
      }
      txt<-paste(txt,"was(were) not in SO(3).")
      message(txt)
      class(x)<-"SO3"
      return(x)
      
    }
  }else if(p==3){
    
  #If there are 3 columns, it's assumed the input R is the matrix of unit axes of rotations and the theta vector are the angles,
    #or the length of the axes is the angle of rotation
    
    U<-matrix(x,n,3)
  
    ulen<-sqrt(rowSums(U^2)) 
    ntheta<-length(theta)
    
    if(is.null(theta)){
      
      theta<-ulen%%(pi)
      
    }else if(ntheta!=n){
      
      if(ntheta==1){
        
        theta<-rep(theta,n)
        
      }else{
        
        stop("Number of angles must match number of axes")
        
      }
      
    }
  
    R<-matrix(NA,n,9)

    for(i in 1:n){
    
      if(ulen[i]!=0)
        U[i,]<-U[i,]/ulen[i]
    
      P <- U[i,] %*% t(U[i,])
    
      R[i,] <- P + (diag(3) - P) * cos(theta[i]) + eskew(U[i,]) * sin(theta[i])
    }
    class(R) <- "SO3"
    return(R)
    
  }else if(p==4){
    
    #If there are 4 columns, it's assumed the input is an n-by-4 matrix with rows corresponding to quaternions 
    R<-as.Q4(x)
    return(as.SO3(R))
    
  }
  
  stop("Unknown data type.  Please see ?SO3 for more details.")

}

# C++ version still isn't working, comment out for now
# SO3.default <- function(U, theta=NULL) {
#   
# 	n<-length(U)/3
# 	
# 	if(n%%1!=0)
# 		stop("Each axis must be in three-dimensions")
# 	
# 	U<-matrix(U,n,3)
# 	ulen<-sqrt(rowSums(U^2)) 
# 	
# 	if(is.null(theta)){ 
# 		theta<-ulen%%pi
# 	}
# 	
# 	ntheta<-length(theta)	
# 	
# 	if(n!=ntheta)
# 		stop("Number of angles must match number of axes")
# 	
# 	if(any(ulen!=1))
# 		U<-U/ulen
# 
# 	R<-SO3defaultC(U,theta)
#  		
#  	class(R) <- "SO3"
#   return(R)
# }


#' @rdname SO3
#' @method as.SO3 Q4
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

as.SO3.Q4<-function(x,...){
  
  q<-formatQ4(x)
  
  if(any((rowSums(q^2)-1)>10e-10)){
    warning("Unit quaternions required.  Input was normalized.")
    nonq<-which((rowSums(q^2)-1)>10e-10)
    q[nonq,]<-as.Q4(q[nonq,]/sqrt(rowSums(q[nonq,]^2)))
  }else{
    class(q)<-"Q4"
  }
  
  theta<-mis.angle(q)
  
  u<-mis.axis(q)
  
  return(as.SO3.default(u, theta)) 
}

#' @rdname SO3
#' @method as.SO3 SO3
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

as.SO3.SO3<-function(x,...){
  return(x)
}

#' @rdname SO3
#' @method as.SO3 data.frame
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

as.SO3.data.frame <- function(x,...) {
  n<-nrow(x)
  p<-ncol(x)
  R<-as.matrix(x,n,p)
  return(as.SO3.default(R))
}

#' @rdname SO3
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

is.SO3 <- function(x) {
	
  Rlen<-length(x)
  
  if(Rlen%%9!=0){
    
    return(FALSE)
    
  }
  R<-x  
  if(class(x)=='data.frame')
    R<-data.matrix(R)
    
  R<-matrix(R,ncol=9)

    
  apply(R,1,
	  function(R){R <- matrix(R, 3, 3)
	  if(any(is.na(R))) return(FALSE)
	  if(abs(det(R)-1)>10e-10) return(FALSE)
	  return(all(abs(t(R) %*% R - diag(1, 3))<10e-5))}) 
    
	
}

#' @rdname SO3
#' @aliases as.SO3 is.SO3 id.SO3 as.SO3.default as.SO3.Q4 as.SO3.SO3 as.SO3.data.frame
#' @export

id.SO3 <- as.SO3(c(1,0,0),0)
