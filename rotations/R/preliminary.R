#' Rotational distance
#'
#' Calculate the extrinsic or intrinsic distance between two rotations.
#'
#' This function will calculate the intrinsic (Riemannian) or extrinsic (Euclidean) distance between two rotations.
#' \code{R2} and \code{Q2} are set to the identity rotations by default.  For rotations \eqn{R_1}{R1} and \eqn{R_2}{R2}
#' both in \eqn{SO(3)}, the Euclidean distance between them is \deqn{||R_1-R_2||_F}{||R1-R2||} where \eqn{||\cdot||_F}{|| ||} is the Frobenius norm.
#' The Riemannian distance is defined as \deqn{||Log(R_1^\top R_2)||_F}{||Log(R1'R2)||} where \eqn{Log} is the matrix logarithm, and it corresponds
#' to the misorientation angle of \eqn{R_1^\top R_2}{R1'R2}.  See the vignette `rotations-intro' for a comparison of these 
#' two distance measures.
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param R2,Q2 the second rotation in the same parameterization as x.
#' @param method string indicating "extrinsic" or "intrinsic" method of distance. 
#' @param p the order of the distance.
#' @param ... additional arguments.
#' @return The rotational distance between each rotation in x and R2 or Q2.
#' @export
#' @examples
#' rs <- rcayley(20, kappa = 10)
#' Rs <- genR(rs, S = id.SO3)
#' dEs <- rot.dist(Rs,id.SO3)
#' dRs <- rot.dist(Rs, id.SO3 , method = "intrinsic")
#' 
#' #The intrinsic distance between the true central orientation and each observation
#' #is the same as the absolute value of observations' respective misorientation angles
#' all.equal(dRs, abs(rs))              #TRUE
#' 
#' #The extrinsic distance is related to the intrinsic distance
#' all.equal(dEs, 2*sqrt(2)*sin(dRs/2)) #TRUE

rot.dist<-function(x,...){
  UseMethod("rot.dist")
}


#' @rdname rot.dist
#' @method rot.dist SO3
#' @export 

rot.dist.SO3 <- function(x, R2=id.SO3, method='extrinsic' , p=1,...) {
  
  R1<-formatSO3(x)
  
  method <- try(match.arg(method,c('projected','extrinsic','intrinsic')),silent=T)
  
  if (class(method)=="try-error")
    stop("method needs to be one of 'projected', 'extrinsic' or 'intrinsic'.")
  
  if(method%in%c('projected','extrinsic')){
    
    n<-nrow(R1)
    R1<-matrix(R1,n,9)
    R2<-matrix(R2,n,9,byrow=TRUE)
    
    so3dist<-sqrt(rowSums((R1-R2)^2))^p
    
  }else if(method=='intrinsic'){
    
    R2<-matrix(R2,3,3)
    
    thetas<-c(rdistSO3C(R1,R2))
    
    so3dist<-thetas^p
    
  }else{
    stop("Incorrect usage of method argument.  Please choose intrinsic or extrinsic")
  }
  
  return(so3dist)
  
}


#' @rdname rot.dist
#' @method rot.dist Q4
#' @export 

rot.dist.Q4 <- function(x, Q2=id.Q4 ,method='extrinsic', p=1,...) {

  Q1<-formatQ4(x)
  Q2<-formatQ4(Q2)
  
  method <- try(match.arg(method,c('projected','extrinsic','intrinsic')),silent=T)
  
  if (class(method)=="try-error")
    stop("method needs to be one of 'projected', 'extrinsic' or 'intrinsic'.")
  
  if(method=='intrinsic'){
    
    q4dist<-c(RdistC(Q1,Q2))^p
    
  }else if(method%in%c('projected','extrinsic')){
    
    q4dist<-c(EdistC(Q1,Q2))^p
    
  }else{
    stop("Incorrect usage of method argument.  Please choose intrinsic or extrinsic")
  }
  
  return(q4dist)
}


#' Misorientation angle
#' 
#' Compute the misorientation angle of a rotation.
#' 
#' Every rotation can be thought of as some reference coordinate system rotated about an axis through an angle.  These quantities
#' are referred to as the misorientation axis and misorientation angle, respectively, in the material sciences literature.
#' This function returns the misorentation angle associated with a rotation assuming the reference coordinate system
#' is the identity.
#'  
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @return Angle of rotation.
#' @seealso \code{\link{mis.axis}}
#' @export
#' @examples
#' rs <- rcayley(20, kappa = 20)
#' Rs <- genR(rs, S = id.SO3)
#' mis.angle(Rs)
#' 
#' #If the central orientation is id.SO3 then mis.angle(Rs) and abs(rs) are equal
#' all.equal(mis.angle(Rs), abs(rs))  #TRUE
#' 
#' #For other reference frames, the data must be centered first
#' S <- genR(pi/2)
#' RsS <- genR(rs, S = S)
#' mis.axis(RsS-S)
#' all.equal(mis.angle(RsS-S),abs(rs)) #TRUE
#' 
#' #If the central orientation is NOT id.SO3 then mis.angle(Rs) and abs(rs) are usual unequal
#' Rs <- genR(rs, S = genR(pi/8))
#' all.equal(mis.angle(Rs), abs(rs))  #Mean relative difference > 0

mis.angle<-function(x){
  UseMethod("mis.angle")
}


#' @rdname mis.angle
#' @method mis.angle SO3
#' @export 

mis.angle.SO3 <- function(x){
	
	Rs<-formatSO3(x)
	theta<-c(rdistSO3C(Rs,diag(1,3,3)))
  return(theta)
}


#' @rdname mis.angle
#' @method mis.angle Q4
#' @export 

mis.angle.Q4 <- function(x){
	
  Qs<-formatQ4(x)
	theta<-2*acos(Qs[,1])
  class(theta)<-"numeric"
	return(theta)
}


#' Misorientation axis
#' 
#' Determine the misorientation axis of a rotation.
#' 
#' Every rotation can be interpreted as some reference coordinate system rotated about an axis through an angle.  These quantities
#' are referred to as the misorientation axis and misorientation angle, respectively, in the material sciences literature.
#' This function returns the misorentation axis associated with a rotation assuming the reference coordinate system
#' is the identity.  The data must be centered before calling \code{mis.axis} if a different coordinate system is required.
#' 
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param ... additional arguments.
#' @return Axis in form of three dimensional vector of length one.
#' @seealso \code{\link{mis.angle}}
#' @export
#' @examples
#' rs <- rcayley(20, kappa = 20)
#' 
#' #If the reference frame is set to id.SO3 then no centering is required
#' Rs <- genR(rs, S = id.SO3)
#' mis.axis(Rs)
#' all.equal(Rs, as.SO3(mis.axis(Rs), mis.angle(Rs)))
#' 
#' #For other reference frames, the data must be centered first
#' S <- genR(pi/2)
#' RsS <- genR(rs, S = S)
#' mis.axis(RsS-S)
#' all.equal(mis.angle(RsS-S),abs(rs)) #TRUE
#' 
#' Qs <- genR(rs, S = id.Q4, space = "Q4")
#' mis.axis(Qs)
#' all.equal(Qs, as.Q4(mis.axis(Qs), mis.angle(Qs)))

mis.axis<-function(x,...){
  UseMethod("mis.axis")
}

#' @rdname mis.axis
#' @method mis.axis SO3
#' @export 

mis.axis.SO3<-function(x,...){
  
	R<-formatSO3(x)
  n<-nrow(R)
	u<-matrix(NA,n,3)
	
	for(i in 1:n){
		Ri<-matrix(R[i,],3,3)
  	X <- Ri - t(Ri)
  	u[i,] <- rev(X[upper.tri(X)])*c(-1,1,-1)
    norm<-sqrt(sum(u[i,]^2))
    
    if(norm!=0){
      u[i,]<-u[i,]/norm
    }

	}

  return(u) # will be trouble, if R is symmetric, i.e. id,  .... 
  
}

#' @rdname mis.axis
#' @method mis.axis Q4
#' @export 

mis.axis.Q4<- function(x,...){
  
  q<-formatQ4(x)
  theta<-mis.angle(q)
  
  u <- q[,2:4]/sin(theta/2)

	if(any(is.infinite(u)|is.nan(u))){
		infs<-which(is.infinite(u)|is.nan(u))
		u[infs]<-0
	}  
  
  u<-matrix(u,ncol=3)
  
  return(u)
}

eskew <- function(U) {
  
  ulen<-sqrt(sum(U^2))
  
  if(ulen!=0){
    U<-U/ulen
  }
  
  u <- U[1]
  v <- U[2]
  w <- U[3]
  
  res <- matrix((-1) * c(0, -w, v, w, 0, -u, -v, u, 0), ncol = 3)
  return(res)
}



#' Generate rotations
#'
#' Generate rotations in matrix format using Rodrigues' formula or quaternions.
#'
#' Given a vector \eqn{U=(u_1,u_2,u_3)^\top\in R^3}{U=(u1,u2,u3)' in R^3} of length one and angle of rotation \eqn{r}, a \eqn{3\times 3}{3-by-3} rotation 
#' matrix is formed using Rodrigues' formula
#' \deqn{\cos(r)I_{3\times 3}+\sin(r)\Phi(U)+(1-\cos(r))UU^\top}{cos(r)I+sin(r)\Phi(U)+(1-cos(r))UU'} 
#' where \eqn{I_{3\times 3}}{I} is the \eqn{3\times 3}{3-by-3} identity matrix, \eqn{\Phi(U)} is a \eqn{3\times 3}{3-by-3} skew-symmetric matrix
#' with upper triangular elements \eqn{-u_3}{-u3}, \eqn{u_2}{u2} and \eqn{-u_1}{-u1} in that order.
#' 
#' For the same vector and angle a quaternion is formed according to \deqn{q=[cos(\theta/2),sin(\theta/2)U]^\top.}{q=[cos(theta/2),sin(theta/2)U]'.}
#'
#' @param r vector of angles.
#' @param S central orientation.
#' @param space indicates the desired representation: rotation matrix "SO3" or quaternions "Q4." 
#' @return A \eqn{n\times p}{n-by-p} matrix where each row is a random rotation matrix (\eqn{p=9}) or quaternion (\eqn{p=4}).
#' @export
#' @examples
#' r <- rvmises(20, kappa = 0.01)
#' Rs <- genR(r, space = "SO3")
#' Qs <- genR(r, space = "Q4")

genR <- function(r, S = NULL, space='SO3') {
  
  if(!(space %in% c("SO3","Q4")))
    stop("Incorrect space argument.  Options are: SO3 and Q4. ")
  
  n<-length(r)

  theta <- acos(runif(n, -1, 1))
  
  # Generate angles phi from a uniform distribution from -pi to pi
  
  phi <- runif(n, -pi, pi)
  u <- matrix(c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)),n,3)
  
  if(space=="SO3"){
  	
  	#For now the C++ code is broken, use R functions
  	#S<-matrix(S,3,3)
  	#o<-SO3defaultC(u,r)
  	#o<-genrC(r,S,1,u)

  	o<-as.SO3.default(x=u,theta=r)
    
  	if(is.null(S)){
  		
  		class(o) <- "SO3"
  		return(o)
  		
  	}else{
      
      if(is.Q4(S)){
        S <- as.SO3(S)
      }
      
  	  S<-formatSO3(S)
      St<-t(matrix(S,3,3))
  	  o<-center.SO3(o,St)
  	  class(o) <- "SO3"
  	  return(o)
  	
  	}
  	
  }else{
  	
  	#S<-matrix(S,1,4)
  	#q<-Q4defaultC(u,r)
  	#q<-genrC(r,S,2,u)
  	
  	q<-matrix(c(cos(r/2),sin(r/2)*u),n,4)
  	
  	if(is.null(S)){
  		
  		class(q)<-"Q4"
  		return(q)
  		
  	}else{
  	
      if(is.SO3(S)){
        S <- as.Q4(S)
      }
      
  		S<-formatQ4(S)
      S<--S
  		q<-center.Q4(q,S)
  	
  		class(q)<-"Q4"
  		return(q)
  	}
  	
  }

}


#' Matrix exponential
#'
#' Compute the matrix exponential for skew-symmetric matrices according to the usual Taylor expansion.
#' The expansion is significantly simplified for skew-symmetric matrices, see \cite{moakher02}.
#' Maps a matrix belonging to the lie algebra \eqn{so(3)} into the lie group \eqn{SO(3)}.
#'
#' @param x single \eqn{3\times 3}{3-by-3} skew-symmetric matrix or \eqn{n\times 9}{n-by-9} sample of skew-symmetric matrices.
#' @return Matrix \eqn{e^{\bm H}}{e^H} in \eqn{SO(3)} .
#' @@cite moakher02
#' @export
#' @examples
#' Rs <- ruars(20, rcayley)
#' lRs <- log(Rs)           #Take the matrix logarithm for rotation matrices
#' Rs2 <- skew.exp(lRs)     #Go back to rotation matrices
#' all.equal(Rs, Rs2)

skew.exp <- function(x) {

  if(length(x)==9){
    
    H<-matrix(x,3,3)
    Hmat<-expskewC(H)
    class(Hmat)<-"SO3"
    return(Hmat)
    
  }else{
    Hmat<-expskewCMulti(x)
    class(Hmat)<-"SO3"
    return(Hmat)
  }
}


#' Rotation logarithm
#'
#' Compute the logarithm of a rotation matrix, which results in a \eqn{3\times 3}{3-by-3} skew-symmetric matrix.  This function maps
#' the lie group \eqn{SO(3)} into its tangent space, which is the space of all \eqn{3\times 3}{3-by-3} skew symmetric matrices,
#' the lie algebra \eqn{so(3)}.  For details see e.g. \cite{moakher02}.
#'
#' @param x \eqn{n\times 9}{n-by-9} matrix where each row corresponds to a random rotation matrix.
#' @param ... additional arguments.
#' @return Skew symmetric matrix \eqn{\log(R)}{log(R)}.
#' @@cite moakher02
#' @export 
#' @method log SO3
#' @examples
#' Rs <- ruars(20, rcayley)
#' 
#' #Here we demonstrate how the logarithm can be used to determine the angle and 
#' #axis corresponding to the provided sample
#' 
#' lRs <- log(Rs)               #Take the logarithm of the sample
#' Ws <- lRs[,c(6, 7, 2)]       #The appropriate diagonal entries are the axis*angle
#' lens <- sqrt(rowSums(Ws^2))  
#' axes <- mis.axis(Rs)
#' angs <- mis.angle(Rs)
#' all.equal(axes, Ws/lens)
#' all.equal(angs, lens)


log.SO3 <- function(x,...) {
  if(length(x)==9){
	  x<-matrix(x,3,3)
	  return(logSO3C(x))
  }else{
    return(logSO3CMulti(x))
  }
  
}

#' Projection into SO(3)
#'
#' Project an arbitrary \eqn{3\times 3}{3-by-3} matrix into \eqn{SO(3)}.
#'
#' This function uses the process detailed in Section 3.1 of \cite{moakher02} to project an arbitrary \eqn{3\times 3}{3-by-3} matrix into \eqn{SO(3)}.
#' More specifically it finds the closest orthogonal 3-by-3 matrix with determinant one to the provided matrix.
#' 
#' @param M \eqn{3\times 3}{3-by-3} matrix to project into \eqn{SO(3)}.
#' @return Projection of \eqn{\bm M}{M} into \eqn{SO(3)}.
#' @seealso \code{\link{mean.SO3}}, \code{\link{median.SO3}}
#' @export
#' @examples
#' #Project an arbitrary 3x3 matrix into SO(3)
#' M<-matrix(rnorm(9), 3, 3)
#' project.SO3(M)
#' 
#' #Project a sample arithmetic mean into SO(3), same as 'mean'
#' Rs <- ruars(20, rcayley)
#' Rbar <- colSums(Rs)/nrow(Rs)
#' project.SO3(Rbar)              #The following is equivalent
#' mean(Rs)         

project.SO3 <- function(M) {
  
	M<-matrix(M,3,3)
  R<-projectSO3C(M)
  return(R)
}


#' Sample distance
#'
#' Compute the sum of the \eqn{p^{th}}{pth} order distances between each row of x and S.
#'
#' @name rotdist.sum
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param S the individual matrix of interest, usually an estimate of the mean.
#' @param method type of distance used method in "extrinsic" or "intrinsic"
#' @param p the order of the distances to compute.
#' @return The sum of the pth order distance between each row of x and S.
#' @seealso \code{\link{rot.dist}}
#' @aliases rotdist.sum.SO3 rotdist.sum.Q4
#' @export
#' @examples
#' Rs <- ruars(20, rvmises, kappa = 10)
#' 
#' SE1 <- median(Rs)                      #Projected median
#' SE2 <- mean(Rs)                        #Projected mean
#' SR2 <- mean(Rs, type = "geometric")    #Geometric mean
#' 
#' #I will use "rotdist.sum" to verify these three estimators minimize the
#' #loss function they are designed to minimize relative to the other esimators.
#' #All of the following statements should evaluate to "TRUE"
#' 
#' #The projected mean minimizes the sum of squared Euclidean distances 
#' rotdist.sum(Rs, S = SE2, p = 2) < rotdist.sum(Rs, S = SE1, p = 2)
#' rotdist.sum(Rs, S = SE2, p = 2) < rotdist.sum(Rs, S = SR2, p = 2)
#' 
#' #The projected median minimizes the sum of first order Euclidean distances 
#' rotdist.sum(Rs, S = SE1, p = 1) < rotdist.sum(Rs, S = SE2, p = 1)
#' rotdist.sum(Rs, S = SE1, p = 1) < rotdist.sum(Rs, S = SR2, p = 1)
#' 
#' #The geometric mean minimizes the sum of squared Riemannian distances 
#' rotdist.sum(Rs, S = SR2, p = 2, method = "intrinsic") < 
#'                  rotdist.sum(Rs, S = SE1, p = 2, method = "intrinsic")
#' rotdist.sum(Rs, S = SR2, p = 2, method = "intrinsic") < 
#'                  rotdist.sum(Rs, S = SE2, p = 2, method = "intrinsic")


rotdist.sum<-function(x, S = genR(0, space=class(x)), method='extrinsic', p=1){
  
  UseMethod( "rotdist.sum" )

}

#' @rdname rotdist.sum
#' @method rotdist.sum SO3
#' @export 

rotdist.sum.SO3 <- function(x, S = id.SO3, method='extrinsic', p=1) {

  return(sum(rot.dist(x,S, method=method, p=p)))
  
}

#' @rdname rotdist.sum
#' @method rotdist.sum Q4
#' @export 

rotdist.sum.Q4 <- function(x, S = id.Q4, method='extrinsic', p=1) {
  
  return(sum(rot.dist(x,S, method=method, p=p)))
  
}

#' Center rotation data
#' 
#' This function will take the sample Rs and return the sample Rs centered at
#' S.  That is, the ith observation of Rs denoted \eqn{R_i}{Ri} is returned as \eqn{S^\top R_i}{S'Ri}.  
#' If S is the true center then the projected mean should be close to the 3-by-3 identity matrix. 
#' 
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param S the rotation or a matrix of \eqn{n\times p}{n-by-p} rotations about which to center each row of x.
#' @return The sample centered about S
#' @export
#' @examples
#' Rs <- ruars(5, rcayley)
#' cRs <- center(Rs, mean(Rs))
#' mean(cRs)                      #Close to identity matrix
#' 
#' all.equal(cRs, Rs - mean(Rs))  #TRUE, center and '-' have the same effect
#'                                #See ?"-.SO3" for more details
#'                                
#' center(Rs,Rs)                  #n-Identity matrices: If the second argument is of the same dimension
#'                                #as Rs then each row is centered around the corresponding
#'                                #row in the first argument

center<-function(x,S){
  
  UseMethod( "center" )
  
}

#' @rdname center
#' @method center SO3
#' @export 

center.SO3<-function(x,S){
	#This takes a set of observations in SO3 and centers them around S
	
	Rs<-formatSO3(x)
  
  if(length(S)==9){
    
	  S<-matrix(formatSO3(S),3,3)
    Rs<-centerCpp(Rs,S)
    
  }else if(nrow(x)==nrow(S)){
    
    for(i in 1:nrow(x)){
      Rs[i,]<-centerCpp(matrix(Rs[i,],1,9),matrix(S[i,],3,3))
    }
    
  }else{
    stop("S must either be a single rotation or have as many rows as x.")
  }
  
  class(Rs)<-"SO3"
	return(Rs)
}


#' @rdname center
#' @method center Q4
#' @export 

center.Q4<-function(x,S){
	#This takes a set of observations in Q4 and centers them around S
	Qs<-formatQ4(x)
	S<-formatQ4(S)
  
  if(length(S)==4){
    S<--S
	
	  for(i in 1:nrow(Qs)){
	  	Qs[i,]<-qMult(S,Qs[i,])
	  }
    
  }else if(nrow(x)==nrow(S)){
    
    for(i in 1:nrow(Qs)){
      Si <- -S[i,]
      Qs[i,]<-qMult(Si,Qs[i,])
    }
    
  }else{
    stop("S must either be a single rotation or have as many rows as x.")
  }
  class(Qs)<-"Q4"
	return(Qs)
}



formatSO3<-function(Rs){
	#This function will take input and format it to work with our functions
	#It also checks that the data is actually SO3 and of appropriate dimension
	
	len<-length(Rs)
	if(len%%9!=0)
		stop("Data needs to have length divisible by 9.")
	
	Rs<-matrix(Rs,len/9,9)
	
	if (!all(is.SO3(Rs))) 
		warning("At least one of the given observations is not in SO(3).  Use result with caution.")
	
	class(Rs)<-"SO3"
	return(Rs)

}

formatQ4<-function(Qs){
	
  #This condition is checked later on
  #if(length(Qs)%%4!=0)
  #  stop("Data needs to have length divisible by 4.")
  
  Qs<-matrix(Qs,length(Qs)/4,4)
  
  if (!all(is.Q4(Qs))) 
  	warning("At least one of the given observations is not a unit quaternion.  Use result with caution.")
  
  
  #if(length(Qs)==4)
  #  return(as.Q4(Qs))
  #else
  class(Qs)<-"Q4"
  return(Qs)
}

pMat<-function(p){
	#Make the matrix P from quaternion p according to 3.1 of Rancourt, Rivest and Asselin (2000)
	#This is one way to multiply quaternions
  p<-as.vector(p)
	Pmat<-matrix(0,4,4)
	Pmat[,1]<-p
	Pmat[,2]<-p[c(2,1,4,3)]*c(-1,1,1,-1)
	Pmat[,3]<-c(-p[3:4],p[1:2])
	Pmat[,4]<-p[4:1]*c(-1,1,-1,1)
	return(Pmat)
}

qMult<-function(q1,q2){
	#Forms quaternion product q1 x q2, i.e., rotate q2 by q1
	#This functions utilizes the 
	q1<-formatQ4(q1)
	q2<-formatQ4(q2)
	q1q2<-pMat(q1)%*%matrix(q2,4,1)
	return(formatQ4(q1q2))
}

proj<-function(u,v){
	#Project the vector v orthogonally onto the line spanned by the vector u
	num<-t(u)%*%v
	denom<-t(u)%*%u
	return(num*u/denom)
}

tLogMat <- function(x, S) {
  tra <- log.SO3(t(S) %*% matrix(x, 3, 3))
  return(as.vector(tra))
}

tLogMat2 <- function(x, S) {
  tra <- log.SO3(matrix(x, 3, 3)%*%t(S))
  return(as.vector(tra))
}

vecNorm <- function(x, S, ...) {
  n <- sqrt(length(x))
  cenX <- x - as.vector(S)
  return(norm(matrix(cenX, n, n), ...))
}
