
ellipseToConicMatrix <- function(saxes=c(1,1), loc=c(0,0), theta = 0){
  if(length(saxes)<2) stop("Arg saxes contains less than 2 elt.\n")
  if(length(loc)<2)  stop("Arg loc contains less than 2 elt.\n")
  if(length(theta)>1){
    theta <- theta[1]
    warning("First element of theta is used\n")
  }
    
  # creation of the matrix in homogeneous coordinates of the ellipse
  CC = matrix(c(1/saxes[1]^2, 0 ,0, 0, 1/saxes[2]^2, 0, 0,0,-1),nrow=3, byrow=TRUE)  
  # Rotation + Translation
  CC <- t(translation(-loc)) %*% t(rotation(theta)) %*% CC %*% rotation(theta) %*% translation(-loc)
  # symmetrization
  CC <- (CC + t(CC))/2
  CC[abs(CC) < .Machine$double.eps^0.95] <- 0
  return(CC)
}

conicMatrixToEllipse <- function(A){
  if(!all.equal(dim(A),c(3,3))) stop("Conic matrix should be of size 3x3\n")
  # if rank 3 and if detB > 0 => ellipse
  detB <- A[1,1]*A[2,2] - 2*A[1,2]
  if(!(qr(A)$rank == 3 && detB > 0)) stop("This is not an ellipse \n")
  b2mac <- (A[1,2]^2-A[1,1]*A[2,2])
  # location
  x0 <- (A[2,2]*A[1,3] - A[1,2]*A[2,3])/b2mac
  y0 <- (A[1,1]*A[2,3] - A[1,2]*A[1,3])/b2mac
  # parameter (semi-axis)
  NUM <- 2*(A[1,1]*A[2,3]^2 + A[2,2]*A[1,3]^2 + A[3,3]*A[1,2]^2 - 2*A[1,2]*A[1,3]*A[2,3]-A[1,1]*A[2,2]*A[3,3])
  sqrtac24b2 <- sqrt((A[1,1]-A[2,2])^2 + 4*A[1,2]^2)
  ac <- A[1,1] + A[2,2]
  a <- sqrt(NUM/(b2mac*(sqrtac24b2 - ac)))
  b <- sqrt(NUM/(b2mac*(-sqrtac24b2 - ac)))
  # theta (orientation)
  if(isTRUE(all.equal(b,0)) && A[1,1] < A[2,2]){
    theta <- 0
  }else if(isTRUE(all.equal(b,0)) && A[1,1] > A[2,2]){
    theta <- 0.5*pi
  }else if(!isTRUE(all.equal(b,0)) && A[1,1] < A[2,2]){
    theta <- 0.5*atan(2*A[1,2]/(A[1,1]-A[2,2]))
  }else if(!isTRUE(all.equal(b,0)) && A[1,1] > A[2,2]){
    theta <- pi/2 + 0.5*atan(2*A[1,2]/(A[1,1]-A[2,2]))
  }
  return(list(loc=c(x0,y0), saxes=c(a,b), theta = theta))
}

quadraticFormToMatrix <- function(v){
  if(length(v)!=6 && !as.vector(v,mode="numeric")) stop("v must have 6 elements and be numeric!\n")
  return(matrix(c(v[1], v[2]/2, v[4]/2, v[2]/2, v[3], v[5]/2, v[4]/2, v[5]/2, v[6]),nrow=3, ncol=3))
}

arcLengthEllipse <- function(p1,p2=NULL,saxes,n=5){
	if( (abs(p1[1]) - saxes[1]) >  .Machine$double.eps^0.75) stop("x-coordinate of point 1 larger/smaller than a!\n")
	L1 <- pEllipticInt(as.numeric(p1[1]),saxes,n=n)
	if(!is.null(p2)){
		if( abs(p2[1]) - saxes[1] > .Machine$double.eps^0.75) stop("x-coordinate of point 2 larger/smaller than a!\n")
		L2 <- pEllipticInt(as.numeric(p2[1]),saxes,n=n)
		# check 1 - m?me quadrant
		if( (length(p1)==1 && length(p2)==1) || sign(p1[2]) == sign(p2[2])){
			#
			LT <- abs(L1 - L2)		
		}else{
			# arc length of a quadrant
			Lqua <- pEllipticInt(as.numeric(saxes[1]),saxes,n=n)
			LT <- 2*Lqua - abs(max(L1,L2)) + abs(min(L1,L2)) 
		}
		return(LT)
	}else{
		return(L1)
	}
}

pEllipticInt <- function(x,saxes,n=5){
	theta <- asin(round(x/saxes[1],12))
	# excentricity
	if(saxes[2]>saxes[1]){
		ex <- sqrt(1-(saxes[1]/saxes[2])^2)
	}else{
		ex <- sqrt(1-(saxes[2]/saxes[1])^2)
	}
	I <- theta
	E <- I
	for(i in 1:n){
		I <- (2*i - 1)/(2*i) * I - sin(theta)^(2*i-1) * cos(theta)/(2*i)
		E <- E - abs(choose(-0.5,i)) * ex^(2*i)/(2*i - 1) *  I
	}
	return(as.numeric(saxes[1]*E))
}

addLine <- function(l,...){
  if(!is.numeric(l) && length(l)<3) stop("arg l is either not numeric or has less than 3 elements")
  if(l[1]==0) abline(h=-l[3]/l[2],...)  	# cas 1. y = y0
  if(l[2]==0) abline(v=-l[3]/l[1],...)		# cas 2. x = x0
  if(l[1]!=0 && l[2]!=0) abline(a=-l[3]/l[2],b=-l[1]/l[2],...)	# cas general
}

join <- function(p,q){
  if(!(is.vector(p, mode="numeric") && is.vector(q, mode="numeric")
       && length(p)==3 && length(q)==3)) 
    stop("p or q are either not a vector or have less than 3 elements\n")
  p <- as.numeric(p)
  q <- as.numeric(q)
  l <- c( p[2]*q[3] - p[3]*q[2],
         -p[1]*q[3] + p[3]*q[1],
          p[1]*q[2] - p[2]*q[1])
  l[abs(l)<.Machine$double.eps^0.75] <- 0
  if(l[3]!=0) l<-l/l[3]
  return(l)
}

meet <- function(l,m){
  if(!(is.vector(l, mode="numeric") && is.vector(m, mode="numeric")
       && length(l)==3 && length(m)==3)) 
    stop("l or m are either not a vector or have less than 3 elements\n")
  p <- c(  l[2]*m[3] - l[3]*m[2],
          -l[1]*m[3] + l[3]*m[1],
           l[1]*m[2] - l[2]*m[1])
  if(p[3]!=0) p <- p/p[3]
  return(p)
}

parallel <- function(p,l){
  if(!(is.vector(p, mode="numeric") && is.vector(l, mode="numeric")
       && length(p)==3 && length(l)==3))
    stop("p or/and l is/are not a vector and/or have less than 3 elements\n")
  join(p,meet(l,c(0,0,1)))
  
}

polar<- function(p,C){
  if(all(dim(C)==c(3,3)) && !isTRUE(all.equal(0,det(C)))){
    ( C %*% p)
  }else{
    if(isTRUE(all.equal(0,det(C)))){
      stop("det(C) == 0\n")
    }else{
      stop("p and C should have dimension [1x3] and [3x3], respectively\n")
    }
  }
}

colinear <- function(p1,p2,p3){
  if(isTRUE(all.equal(0, det(cbind(p1,p2,p3))))){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

rotation <- function(theta, pt = NULL){
  if(!is.numeric(theta)) stop("theta must be numeric!\n")
  rot <- matrix(c(cos(theta), +sin(theta), 0,-sin(theta),cos(theta),0,0,0,1), nrow=3, byrow=TRUE)
  if(is.null(pt)){
    return(rot)
  }else if(is.vector(pt,mode="numeric") && length(pt)>1){
    return(translation(pt)%*%rot%*%translation(-pt))
  }else{
    stop("Error, pt must be numeric and have at least 2 elements")
  }
  
}

translation <- function(v){
  matrix(c(1, 0, v[1], 0 ,1, v[2], 0,0,1), nrow=3, byrow=TRUE)
}

scaling <- function(s){
  matrix(c(s[1], 0, 0, 0 ,s[2], 0, 0,0,1), nrow=3, byrow=TRUE)
}

reflection <- function(alpha){
  rot <- matrix(c(cos(alpha), +sin(alpha), 0,sin(alpha),-cos(alpha),0,0,0,1), nrow=3, byrow=TRUE)
}

minor <- function(A, i, j){
  det( A[-i,-j] )
} 

cofactor <- function(A, i, j){
  (-1)^(i+j) * minor(A,i,j)
}

adjoint <- function(A) {
  n <- nrow(A)
  t(outer(1:n, 1:n, Vectorize(
    function(i,j) cofactor(A,i,j)
  )))
}

skewSymmetricMatrix <- function(p){
  matrix(c(0,p[3], -p[2],-p[3], 0, p[1],p[2], -p[1], 0),byrow=TRUE,ncol=3,nrow=3)
}

conicThrough5Points <- function(p1,p2,p3,p4,p5){
  g1 <- join(p1, p3)
  g2 <- join(p2, p4)
  h1 <- join(p1, p4)
  h2 <- join(p2, p3)
  G <- g1 %*% t(g2)
  H <- h1 %*% t(h2)
  M <- as.numeric(t(p5) %*% H %*% p5) * G - as.numeric(t(p5) %*% G %*% p5) * H
  M + t(M)
}

splitDegenerateConic <- function(C){
  if( !isTRUE(all.equal(0,det(C)))) stop("The conic defined by C is not degenerated\n")
  if(qr(C)$rank==2){
    B <- -adjoint(C)
    i <- which(diag(B)!=0)
    if(length(i)<1) stop("problem 1")
    bet <- sqrt(as.complex(B[i[1],i[1]]))
    p <- B[,i[1]]/bet
    A <- C + skewSymmetricMatrix(p)
  }else if(qr(C)$rank == 1){
    A <- C
  }
  ij <- which(A!=0, arr.ind=TRUE)
  if(nrow(ij) < 1)  stop("problem 2")
  P <- matrix(nrow=3,ncol=2)
  P[,1] <- A[ij[1,1],]
  P[,2] <- A[,ij[1,2]]
  return(P)
}

intersectConicLine <- function(C,l){  
  if(isTRUE(all.equal(c(0,0,0),l))){
	warning("Line undefined (0,0,0), return NULL\n")
	return(NULL)
  }
  M <- skewSymmetricMatrix(l)
  B <- t(M) %*% C %*% M
  l[abs(l) < .Machine$double.eps^0.5]  <-0
  i <- which(l!=0)
  subB_index <- c(1:3)[-i[1]]
  detB <- B[subB_index[1],subB_index[1]]*B[subB_index[2],subB_index[2]] - B[subB_index[1],subB_index[2]]*B[subB_index[2],subB_index[1]]
  if(Re(detB) < 0){
    alpha  <- (1/l[i[1]]) * sqrt(- detB )
    A <- B + alpha*M
    # index for a non-zero value
    idn0 <- which(A!=0, arr.ind=TRUE)
    if(nrow(idn0)>0){
      np <- c(A[idn0[1,1],3]!=0, A[3,idn0[1,2]]!=0)
      P <- matrix(nrow=3,ncol=2)
      P[,1] <- A[idn0[1,1],]  /A[idn0[1,1],3]
      P[,2] <- A[,idn0[1,2]]	/A[3,idn0[1,2]]
      return( P[,np] )
    }
  }
  return(NULL)
}

cubic <- function(p){
  if(abs(p[1])<.Machine$double.eps^0.95){
    stop('Coefficient of highest power must not be zero!\n')
  } 
  if(!(as.vector(p, mode="numeric") && length(p)==4)){
    stop('p is not a numeric or/and has not 4 elements!\n')
  } 
  if(any(is.complex(p))){
    stop("the coefficient must be real")
  }
  a <- numeric(3)
  for(i in 2:4) {a[i-1]=p[i]/p[1]}
  Q <- (a[1]^2 - 3*a[2])/9
  R <- (2*a[1]^3 - 9*a[1]*a[2]+27*a[3])/54
  x <- numeric(3)
  # case 1 - 3 real roots
  if(R^2 < Q^3){
    theta <- acos(R/sqrt(Q^3))
    x[1] <- -2 * sqrt(Q) * cos(theta/3) - a[1]/3
    x[2] <- -2 * sqrt(Q) * cos((theta + 2*pi)/3) - a[1]/3
    x[3] <- -2 * sqrt(Q) * cos((theta - 2*pi)/3) - a[1]/3
  }else{
    A = - sign(R)*( abs(R) + sqrt(R^2-Q^3))^(1/3)
    if(isTRUE(all.equal(0, A))){
      B <- 0
    }else{
      B <- Q/A
    }
    x[1] <- (A + B) - a[1]/3
    x[2] <- -0.5*(A+B) - a[1]/3 + sqrt(3)*complex(real=0,imaginary=1)*(A-B)/2
    x[3] <- -0.5*(A+B) - a[1]/3 - sqrt(3)*complex(real=0,imaginary=1)*(A-B)/2
  }
  return(x)
}

intersectConicConic <- function(C1,C2){
	alp <- det(t(C1))
	bet <- det(rbind(C1[,1], C1[,2],C2[,3])) + det(rbind(C1[,1],C2[,2],C1[,3])) + det(rbind(C2[,1],C1[,2],C1[,3]))
	gam <- det(rbind(C1[,1],C2[,2],C2[,3])) + det(rbind(C2[,1],C1[,2],C2[,3])) + det(rbind(C2[,1],C2[,2],C1[,3]))
	del <- det(t(C2))
	# solve det(lambda*C1 + C2) = det(CC) = 0, CC = degenerate
	# = solve  alp*lambda[2]^3 + bet*lambda[2]^2 * mu + gam*lambda[2]*mu^2 + del*mu^3
	# with mu = 1 	
	lambda <- cubic(c(alp, bet, gam, del))	# fx cubic 
	# select a real value for lambda
	lambdaRe <- Re(lambda[Im(lambda)==0])
	CC <- lambdaRe[1]*C1 + C2  
	CC[abs(CC)<.Machine$double.eps^0.95] <- 0
	# split the degenerate conic CC into two lines
	BB <- -adjoint(CC)
	# indice for a non-zero element of the diagonal
	idn0 <- which(abs(diag(BB))>.Machine$double.eps^0.95)
	if(length(idn0) > 0){
	  Bi <- BB[idn0[1],idn0[1]]
	  if(Re(Bi) <0 ){
		return(NULL)
	  }else{
		beta2 <- sqrt(Bi)
		p <- BB[,idn0[1]]/beta2
		Mp <- matrix(c(0,p[3], -p[2],-p[3], 0, p[1],p[2], -p[1], 0),byrow=TRUE,ncol=3,nrow=3)
		newC <- CC + Mp 
		idn0_2 <- which(abs(newC)>.Machine$double.eps^0.95, arr.ind=TRUE)
		if(nrow(idn0_2)>0){
		   l1 <- newC[idn0_2[1,1],]
		   l2 <- newC[,idn0_2[1,2]]
			myP1 <- intersectConicLine(C1,Re(l1))
			myP2 <- intersectConicLine(C1,Re(l2))
			if(all(myP1!=FALSE) && all(myP2!=FALSE)){
				return(cbind(myP1,myP2))
			}else if(all(myP1==FALSE)){
				return(myP2)
			}else if(all(myP2==FALSE)){
				return(myP1)
			}else{
				return(cbind(myP1,myP2))
			}
		}
	  } 
	  
	}
	return(NULL)
}

ellipse <- function(saxes=c(1,1), loc = c(0,0), theta = 0, n = 201, 
                    method=c("default","angle","distance")){
  method = match.arg(method, c("default","angle","distance"))  
  if(method=="default"){
    phi <- 2*pi*seq(0,1, len = n)
    P <- matrix(nrow=n,ncol=2)
    P[,1] <- saxes[1] * cos(phi)
    P[,2] <- saxes[2] * sin(phi)
  }else if(method=="angle"){
    b <- min(saxes[1],saxes[2])
    a <- max(saxes[1],saxes[2])
    d2 <- (a-b)*(a+b)                   #= a^2 - b^2
    phi <- 2*pi*seq(0,1, len = n)
    sp <- sin(phi)
    cp <- cos(phi)
    r <- a*b / sqrt((saxes[2] * cp)^2 +  (saxes[1] * sp)^2)
    P <- matrix(nrow=n,ncol=2)
    P[,1] <- r * cp
    P[,2] <- r * sp
  }else if(method=="distance"){
    n <- round(n/4)*4
    phi <- 0.5*pi*seq(0,1, by = 1/n)
    P <- matrix(nrow=length(phi),ncol=2)
    P[,1] <- saxes[1] * cos(phi)
    P[,2] <- saxes[2] * sin(phi)
    d2<-c(0,cumsum(sqrt(apply((diff(P))^2,1,sum))))
    phi_new <- approx(d2,phi,xout = seq(0,tail(d2, n=1), length.out=n/4))$y
    phi_new <- phi_new[-length(phi_new)]
    phi_new <- c(phi_new, pi/2,pi - phi_new, phi_new + pi,3*pi/2, 2*pi - phi_new)
    P <- matrix(nrow=length(phi_new),ncol=2)
    P[,1] <- saxes[1] * cos(phi_new)
    P[,2] <- saxes[2] * sin(phi_new)
  }
  if(theta != 0){
    P <- P %*% matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),byrow=TRUE,nrow=2,ncol=2)
  }
  P <- P + matrix(loc[1:2],nrow=nrow(P),ncol=2,byrow=TRUE)
  return(P)
}

