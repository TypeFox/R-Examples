


#######################################################################
## hausdorff(A,B)
## Compute the Hausdorff distances between two point sets, A and B, in
## an arbitrary dimensional vector space.
## A: the rows of this matrix correspond to the first point set
## B: the rows of this matrix corresponds to the second point set
## A and B may have different numbers of rows (different cardinality)
## but must have the same number of columns (i.e. dimensionality)
##
## H(A,B) = max(h(A,B),h(B,A))
## h(A,B) = max_{a in A} min_{b in B} {d(a,b)}
## d is the Euclidean distance
##
#######################################################################

hausdorff = function(A,B){

  if(ncol(A) != ncol(B)){
    cat("WARNING: dimensionality must be the same \n")
    return(NA)
  }

  d1 = max(apply(A,1,compute.dist,B0=B))
  d2 = max(apply(B,1,compute.dist,B0=A))
  result = max(d1,d2)
  result 
}

######################################################################
## compute.dist      Auxiliary function
#######################################################################

compute.dist = function(a0,B0){
  C0 = matrix(rep(a0,each=nrow(B0)),byrow=F,ncol=ncol(B0))
  result = min(apply(C0-B0,1,L2norm))
  result
}

#######################################################################
## L2norm(x)   Auxiliary function
## Computes the L2-norm of vector x
#######################################################################

L2norm = function(x)  sqrt(sum(t(x)*x)) 

#######################################################################
## dhn = function(A,B)
## Hausdorff distance between the point sets A and B
#######################################################################

dhn = function(A,B){
  result = NA
  if(ncol(A) > 2) result = hausdorff(A,B)
  if(ncol(A) == 2) result = distconvhull(A,B)
  result 
}


#######################################################################
## distconvhull(pp1,pp2)
## Hausdorff distance between the convex hulls of the point sets
## pp1 and pp2
#######################################################################

distconvhull = function(pp1,pp2)
  max(max(apply(pp1,1,distShape,pp=pp2)),max(apply(pp2,1,distShape,pp=pp1)))



#######################################################################
## simplex()
##  From the original samples x and y (of dimension d) we obtain
##  Hausdorff distances between the convex hulls of samples
##  x1 and x2 obtained from x and they are saved into dhx
##  Analogously for y and the distances are saved into dhy
## 
##  bootstrap = TRUE x1 and x2 are random selections of size d+1 from x
##                  obtained independently (probably repeating rows of x)
##                  i.e. x1 and x2 are bootstrap samples. In this
##                  case nresamples is the total number of pairs (x1,x2)
##                  used.
##  
##               FALSE then the pairs (x1,x2) are taken from a partition
##                  of x where each subset has d+1 elements. In this case the
##                  total number of distances is limited by the number of
##                  rows of x. 
## 
## 
#######################################################################

simplex = function(x,bootstrap=TRUE,nresamples=10){
  dhx = NULL;
  rx = nrow(x)
  cx = ncol(x)
  ix = 1:rx;
  nn = cx + 1;
  for(i in 1:nresamples){
    if(bootstrap){
      ixa1 = sample(ix,nn)
      ixa2 = sample(setdiff(ix,ixa1),nn);
      xa1 = x[ixa1,]
      xa2 = x[ixa2,]
      dhx = c(dhx,dhn(xa1,xa2));
    }
    else {
      ##Number of samples of x
      nsamplesx = floor(rx/nn);
      for(i in seq(0,(nsamplesx-2),2)){
        xa1 = x[(i*nn+1):((i+1)*nn),]
        xa2 = x[((i+1)*nn+1):((i+2)*nn),]
        dhx = c(dhx,dhn(xa1,xa2))
      }
    }
  }
  dhx
}


#############################################################################
## FUNCTIONS FOR THE EXACT ALGORITHM OF THE 2D SIMPLEX DISPERSION ORDERING 
#############################################################################

#######################################################################
##  rotatePHA
##Given a 2D point w, the point w is rotated to place it in positive part
##of the x-axis. The angle of the rotation is calculated 
##
  
rotatePHA = function(w){
    if(w[1] >0 & w[2] >= 0) tau = atan(abs(w[2] / w[1])) 
    if(w[1] == 0 & w[2] >= 0) tau = pi/2 
    if(w[1] < 0 & w[2] >= 0) tau = pi -  atan(abs(w[2] / w[1]))
    if(w[1] < 0 & w[2] < 0) tau = pi +  atan(abs(w[2] / w[1]))
    if(w[1] == 0 & w[2] < 0) tau = 3*pi/2
    if(w[1] > 0 & w[2] < 0) tau = 2*pi - atan(abs(w[2] / w[1]))
    tau
  }

#######################################################################
##  whichShape
##Given the 3x2 matrix, we evaluate which shape is the
##convex hull of the 2D points corresponding with the rows
##
#######################################################################

whichShape = function(pp){
    C = matrix(data=c(pp[2,1] - pp[1,1],pp[2,2] - pp[1,2],pp[3,1] - pp[1,1],
                 pp[3,2] - pp[1,2]),ncol=2,byrow=T)
    detC = det(C)
    if(detC != 0) typeShape = "triangle"
    if(detC == 0) typeShape = "segment"
    if(pp[1,1] == pp[2,1] & pp[1,1] == pp[3,1] &
       pp[1,2] == pp[2,2] & pp[1,2] == pp[3,2]) typeShape = "point"
    typeShape
  }

#########################################################################
###  moveShapeAndPoint
###  Placing the shape determined by pp and the point z
###  in an appropiate position
#########################################################################

moveShapeAndPoint = function(z,pp){
  typeShape = whichShape(pp)
  result = NA
  if(typeShape == "triangle"){
    cosenos = NULL 
    for(j in 1:3){
      ik = setdiff(1:3,j)
      pp0 = pp[ik,] - matrix(data=rep(pp[j,],2),ncol=2,byrow=T)
      cosenos = c(cosenos,sum(pp0[1,]*pp0[2,])/sqrt(sum(pp0[1,]*pp0[1,]))/
        sqrt(sum(pp0[2,]*pp0[2,])))
    }
    aa = acos(cosenos)
    if(sum(is.na(aa)) == 0){
      pp1 = pp[sort(aa,index.return=T)$ix,]
      ## Translation of the triangle
      pp2 = pp1 - matrix(data=rep(pp1[1,],3),ncol=2,byrow=T)
      z2 = z - pp1[1,]
      ## Angle for the rotation of pp2[2,] to the positive x-axis
      tau = rotatePHA(pp2[2,])
      ## Rotation of pp2 and z2
      AA = rbind(c(cos(tau),sin(tau)),c(-sin(tau),cos(tau)))
      pp3 = t(AA %*% t(pp2))      
      z3 = t(AA %*% t(t(z2)))
      if(pp3[3,2] < 0) {
        pp3[3,2] = abs(pp3[3,2])
        z3[2] = -z3[2]
      }
      result = list(z = z3,pp = pp3,typeShape = typeShape,isna = FALSE)
    }
  }
    if(typeShape == "segment"){
      ##Relabeling the points 
      pp1 = pp
      if(pp[1,1] == pp[2,1] & pp[1,2] == pp[2,2]) pp1 = pp[c(1,3,2),]
      ##Translation
      pp2 = pp1 - matrix(data=rep(pp1[1,],3),ncol=2,byrow=T)
      z2 = z - pp1[1,]
      ## Angle for the rotation of pp2[2,] to the positive x-axis
      tau = rotatePHA(pp2[2,])
      ## Rotation of pp2 and z2
      AA = rbind(c(cos(tau),sin(tau)),c(-sin(tau),cos(tau)))
      pp3 = t(AA %*% t(pp2))      
      z3 = t(AA %*% t(t(z2)))
      result = list(z = z3,pp = pp3,typeShape = typeShape,isna = FALSE)
    }
    if(typeShape == "point"){
      z3 = z
      pp3 = pp
      result = list(z = z3,pp = pp3,typeShape = typeShape,isna = FALSE)
    }
  result
}

#################################################################
####    distTriangle
## It evaluates the distance from the point z to the triangle 
## with vertices pp  
#################################################################
                               
distTriangle = function(z,pp,correctPosition=FALSE){
	if(correctPosition)	{
		z0 = z
		pp0 = pp} else {
		zpp = moveShapeAndPoint(z,pp)
	    z0 = zpp$z
	    pp0 = zpp$pp
	    }
	
	distance = 0
        
   ##Region 1
   if((z0[2] - pp0[3,2]* z0[1]/pp0[3,1] > 0) & 
   (z0[2] + pp0[3,1]* z0[1]/pp0[3,2] > 0) &
   (z0[2]  - pp0[3,2] + pp0[3,1]*(z0[1] - pp0[3,1])/pp0[3,2] <= 0)){
   distance = abs(z0[2] - pp0[3,2]*z0[1]/pp0[3,1])/sqrt(1+(pp0[3,2]/pp0[3,1])^2)
   regionR = 1
   }

   ##Region 2
   if((z0[2] - pp0[3,2] + (pp0[3,1] -pp0[2,1])*(z0[1] - pp0[3,1])/pp0[3,2] > 0) &
   (z0[2] - pp0[3,2] + pp0[3,1]*(z0[1] - pp0[3,1])/pp0[3,2] > 0)){
   distance = dist(rbind(z0,pp0[3,]))
      regionR = 2
   }
   
   ##Region 3
   if((z0[2] - pp0[3,2]*(z0[1] - pp0[2,1])/(pp0[3,1] - pp0[2,1]) >  0 ) &
   (z0[2] + (pp0[3,1] - pp0[2,1])*(z0[1] - pp0[2,1])/pp0[3,2] >  0 ) &
   (z0[2] - pp0[3,2] + (pp0[3,1]-pp0[2,1])*(z0[1] - pp0[3,1])/pp0[3,2] <= 0))
     {
       distance = abs(z0[2] - pp0[3,2] * (z0[1] - pp0[2,1])/
         (pp0[3,1]-pp0[2,1]))/sqrt(1+ (pp0[3,2]/(pp0[3,1]-pp0[2,1]))^2)
               regionR = 3
     }
   
   ##Region 4
   if((z0[2] + (pp0[3,1]-pp0[2,1])*(z0[1] - pp0[2,1])/pp0[3,2] < 0) &
   (z0[1] > pp0[2,1])){ 
   distance = dist(rbind(z0,c(pp0[,1],0)))
      regionR = 4
   }
 
   ##Region 5
   if((z0[1] >= 0) & (pp0[2,1] >= z0[1]) & (z0[2] < 0)){
   distance = abs(z0[2])
      regionR = 5
   }

   ##Region 6
   if((z0[1] <0) & (z0[2] + pp0[3,1] * z0[1]/pp0[3,2] < 0)){ 
   distance = dist(rbind(z0,0))
      regionR = 6
   }

   distance 
}
 
#################################################################
####    distSegment
## It evaluates the distance from the point z to the segment 
## with vertices in pp                                 
#################################################################

distSegment = function(z,pp,correctPosition=FALSE){
		if(correctPosition)	{
		z0 = z
		pp0 = pp} else {
		zpp = moveShapeAndPoint(z,pp)
	    z0 = zpp$z
	    pp0 = zpp$pp
	    }
	rightend = max(pp0[,1])
	if(0 <= z0[1] & z0[1] <= rightend) distance = abs(z0[2])
	if(rightend < z0[1]) distance = dist(rbind(z0,c(rightend,0)))
	if(z0[1] < 0) distance = dist(rbind(z0,c(0,0)))
    distance
}
 
#################################################################
####    distShape
## Given the point z and the points pp, the function gives us the 
## distance from the point z to the convex hull of the points pp
#################################################################
    
distShape = function(z,pp){
  zinpp = 0
  for(i in 1:nrow(pp)){
    if(z[1] == pp[i,1] & z[2] == pp[i,2]){
      zinpp = 1
      distance = 0
    }
  }

  if(zinpp == 0){
    distance = NA
    zpp = moveShapeAndPoint(z,pp)
    if(!is.na(zpp)){
      z0 = zpp$z
      pp0 = zpp$pp
      typeShape0 = zpp$typeShape
      if(typeShape0 == "triangle") 
        distance = distTriangle(z0,pp0,correctPosition=TRUE)
      if(typeShape0 == "segment") 
        distance = distSegment(z0,pp0,correctPosition=TRUE)
      if(typeShape0 == "point") 
        distance = dist(rbind(z,pp[1,]))
    } else {
      distance = 0
    }
  }
  distance
}

  

#############################################################################
## FUNCTIONS FOR THE EXACT ALGORITHM OF THE 2D HAUSDORFF DISPERSION ORDERING ##
####                     hdo08
#############################################################################

#######################################################################
##  rotateVA
##Given a 2D point w, the point w is rotated to place it in 
##the positive y-axis. The function gives us the rotation angle 
##
#######################################################################
  
rotateVA = function(w){
  if(w[1] >0 & w[2] >= 0) tau = pi/2 - atan(abs(w[2] / w[1])) 
  if(w[1] == 0 & w[2] >= 0) tau = 0
  if(w[1] < 0 & w[2] >= 0) tau = 3*pi/2 + atan(abs(w[2] / w[1]))
  if(w[1] < 0 & w[2] < 0) tau = 3*pi/2 -  atan(abs(w[2] / w[1]))
  if(w[1] == 0 & w[2] < 0) tau = pi
  if(w[1] > 0 & w[2] < 0) tau = pi/2 + atan(abs(w[2] / w[1]))
  tau
}

####################################################################
## doRowAux
## Auxiliary function for exactHausdorff
####################################################################

doRowAux = function(xj,xi,r){
  H = NA
  mxj =  L2norm(xj)
  if(mxj <  r){
    H = L2norm(xi) - r
  } else {
    if(xj[2] <= r){
      H = xi[2] - r
    } else {
      if(xj[1] >0) xj[1] = - xj[1]
      thetaj = acos(xj[1]/mxj) - acos(r/mxj)
      Deltaij = r*r + xi[2]*xi[2]*cos(thetaj)*cos(thetaj)
      if(mxj*mxj <= Deltaij){
        H = L2norm(xi - xj)
      } else {
        H = max(abs(xi[2]*sin(thetaj)-r),
          abs(-sqrt(1-(r*r)/(xi[2]*xi[2]))*xj[1]+(r*xj[2])/xi[2]-r))
      }
    }
  }
  H
}


####################################################################
## exactHausdorff(A,prob,r)
## The rows of A are the possible values
## prob are the probabilities of the rows of A
## r value for the Hausdorff dispersion ordering
## Value:
## distance= Hausdorff distances
## probability = the probability distribution over distance
## alldistances = the whole set of distances with repetitions
####################################################################


exactHausdorff = function(A,prob,r){
  prob = prob %*% t(prob)
  n0 = nrow(A)
  D0 = matrix(data=0,nrow=n0,ncol=n0)
  
  ## Relabeling indices
  
  y0 = apply(A,1,L2norm)
  y1 = sort(y0,decreasing=T,index.return=T)
  A1 = A[y1$ix,]
  n = nrow(A1)
  m =  n - sum(y1$x <= r)
  for(i in 1:(m-1)){
    tau = rotateVA(A1[i,])
    BB = rbind(c(cos(tau),-sin(tau)),c(sin(tau),cos(tau)))
    A2 = A1[i:nrow(A1),]
    A3 = t(BB %*% t(A2))
    if(nrow(A3) > 2){
      D0[i,(i+1):n] = apply(A3[2:nrow(A3),],1,doRowAux,xi=A3[1,],r)
    } else {
      D0[i,i+1] =  doRowAux(A3[2,],A3[1,],r)
    }
  }
  D0 = D0 + t(D0)
  D0 = as.vector(D0)
  D0 = as.data.frame(D0)
  prob = as.vector(prob)
  prob = as.data.frame(prob)
  result = aggregate(prob,D0,sum)
  ## result[,1] are the Hausdorff distances
  ## result[,2] are the probabilities
  distance = result[,1]
  probability = result[,2]
  alldistances = rep(distance,round(n0*(n0+1)*probability/2.))
  resultado = list(distance = distance,probability = probability,
    alldistances = alldistances)
  resultado
}



