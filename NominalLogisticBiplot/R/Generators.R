# file NominalLogisticBiplot/R/Generators.R
# copyright (C) 2012-2013 J.C. Hernandez and J.L. Vicente-Villardon
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#Function that calculates an object with the information about the tesselation, centers, real an dummy points
#, hide categories, etc... from the variable of study.
#----------------------Parameters--------------
  #beta: matrix of parameters for the multinomial logistic regression calculated for one specific dependent variable.
Generators <- function(beta) {
	ngrup = nrow(beta) + 1
	a = matrix(0, ngrup, ngrup)
	b = matrix(0, ngrup, ngrup)

	# Calculate the straight lines separating each pair of categories

	for (i in 1:(ngrup - 1)) {
		a[1, (i + 1)] = -1 * beta[i, 1]/beta[i, 3]
		b[1, (i + 1)] = -1 * beta[i, 2]/beta[i, 3]
	}

	if (ngrup > 2) {
		for (i in 1:(ngrup - 2)) for (j in (i + 1):(ngrup - 1)) {
			a[(i + 1), (j + 1)] = (beta[j, 1] - beta[i, 1])/(beta[i, 3] - beta[j, 3])
			b[(i + 1), (j + 1)] = (beta[j, 2] - beta[i, 2])/(beta[i, 3] - beta[j, 3])
		}
	}

	# Calculate the coordinates of the points in which three lines intersect
	# Those are the candidates to form the voronoi diagram

	triplepoints = (ngrup - 2) * (ngrup - 1) * ngrup/6
	xp = matrix(0, triplepoints, 1)
	yp = matrix(0, triplepoints, 1)
	indices = matrix(0, triplepoints, 3)

	if (ngrup > 3) {
		l = 0
		for (i in 1:(ngrup - 2)) for (j in (i + 1):(ngrup - 1)) for (k in (j + 1):(ngrup)) {
			l = l + 1
			xp[l] = (a[i, k] - a[i, j])/(b[i, j] - b[i, k])
			yp[l] = a[i, j] + b[i, j] * xp[l]
			indices[l, ] = c(i, j, k)
		}
	}else if(ngrup == 3){
          xp[1] = (a[1, 3] - a[1, 2])/(b[1, 2] - b[1, 3])
    			yp[1] = a[1, 2] + b[1, 2] * xp[1]
    	    indices[1, ] = c(1,2,3)
        }

	# Decide if the candidates are real or virtual using the probabilities given by the model
	# Probabilities for each category at the candidate points.
  # Real points will be in the diagram, virtual points will be used to determine the dummy points

	xn = cbind(matrix(1, triplepoints, 1), xp, yp)
	probab = matrix(0, triplepoints, ngrup)

	for (i in 1:triplepoints) {
		suma = 1
		for (j in 1:(ngrup - 1)) suma = suma + exp(sum(beta[j, ] * xn[i, ]))
		for (j in 1:(ngrup - 1)) probab[i, (j + 1)] = exp(sum(beta[j, ] * xn[i, ]))/suma
		probab[i, 1] = 1/suma
	}
	# When a point is too far away the probability of one category is 1 and the rest are 0, but sometimes
	# the calculation results in a NaN

  NoNumbers =which(is.nan(probab))
  for (i in NoNumbers) probab[i]=1

	# Calculating the status (real or virtual) for each point. A point is real when one its three
	# common categories has the highest probability, and virtual otherwise
  # the procedure is good when there are more than 3 categories

	status = matrix(1, triplepoints, 1) #1 for real, 0 for virtual

	if (ngrup > 3) {
		l = 0
		nvirtual = 0
		for (i in 1:(ngrup - 2))
     for (j in (i + 1):(ngrup - 1))
      for (k in (j + 1):(ngrup)){
  			l = l + 1
  			if (max(probab[l, ]) > max(c(probab[l, i], probab[l, j], probab[l, k]))) {
  				status[l] = 0
  				nvirtual = nvirtual + 1
  			}
		  }

		# Separating the real and virtual points, and the corresponding indices
		nreal = triplepoints - nvirtual
		coorreal = matrix(0, nreal, 2)
		coorvirt = matrix(0, nvirtual, 2)
		IndReal = matrix(1, nreal, 3)
		IndVirt = matrix(1, nvirtual, 3)

		nr = 0
		nv = 0
		for (i in 1:triplepoints) {
			if (status[i] == 1) {
				nr = nr + 1
				coorreal[nr, ] = matrix(c(xp[i], yp[i]), 2, 1)
				IndReal[nr, ] = indices[i, ]
			}
			if (status[i] == 0) {
				nv = nv + 1
				coorvirt[nv, ] = c(xp[i], yp[i])
				IndVirt[nv, ] = indices[i, ]
			}
		}
	}else if(ngrup == 3){
	       nvirtual = 0
	       nreal = 1
	       coorreal = matrix(0, nreal, 2)
	       IndReal = matrix(1, nreal, 3)
	       nr = 1
	       nv = 0
	       coorreal[nr, ] = matrix(c(xp[1], yp[1]), 2, 1)
				 IndReal[nr, ] = indices[1, ]
	     }


	# Calculating the tesellation #############
	# Each real point should be connected to other three (either real or virtual)
# at least another one has to be real.
# the virtual points are connected just to one real point.

	Joint = matrix(0, nr, nr)
	# This is the incidence matrix, square matrix with the real points in rows and columns
	# the elements are 1 when te points are connected and 0 otherwise

  if((ngrup > 3) & (nr > 1)){
  	for (i in 1:(nr - 1))
     for (j in (i + 1):nr)
      if (length(intersect(IndReal[i, ], IndReal[j, ])) == 2) {
    		Joint[i, j] = 1
    		Joint[j, i] = 1
    	}
	}

	# variables for the coordinates of the dummy points
	dummy.x = NULL
	dummy.y = NULL
	JointDummy = NULL
	IndDummy = NULL

	# n1, n2 & n3 contain the nodes to which each real point is connected
	# when the node is dummy, the number is negative
  n1 = rep(0, nr)
	n2 = rep(0, nr)
	n3 = rep(0, nr)
	ndummy = 0
	for (i in 1:nr) {

		if (sum(Joint[i, ]) == 3) {
			# When a real point is connected to other three real points, search which are the nodes and store them
			neighbors = which(Joint[i, ] == 1)
			n1[i] = neighbors[1]
			n2[i] = neighbors[2]
			n3[i] = neighbors[3]
		}

		if (sum(Joint[i, ]) == 2) {
			# When a real point is connected to other 2 real points, search which are the nodes and store them
			# then we need a dummy point. The dummy point is obtained using the virtual points.
# The virtual point is located on the remaining line (at each real point intersect exactly 3 lines and we have already used 2)
# on the remaining line there should be several virtual points (depending on the number of categories), one of them
# is enough to calculate the dummy point. Calculating the parameters of the straight line passing through the
# real and the virtual points, the  dumme point is located on the opposite direction to the virtual
# point relative to the real one. The dummy point is placed symmetrically, i. e., the distance
# from the real to the virtual and the distance to the dummy are the same

			ndummy = ndummy + 1
			ColJointDummy = matrix(0, nr, 1)
			ColJointDummy[i] = 1
			JointDummy = cbind(JointDummy, ColJointDummy)
			neighbors = which(Joint[i, ] == 1)
			n1[i] = -1 * ndummy
			n2[i] = neighbors[1]
			n3[i] = neighbors[2]
			A = intersect(IndReal[i, ], IndReal[neighbors[1], ])
			B = intersect(IndReal[i, ], IndReal[neighbors[2], ])
			C = setdiff(union(A, B), intersect(A, B))
			IndDummy = rbind(IndDummy, sort(C))
			for (j in 1:nv) if (length(intersect(C, IndVirt[j, ])) == 2) {
				pend = (coorvirt[j, 2] - coorreal[i, 2])/(coorvirt[j, 1] - coorreal[i, 1])
				const = coorreal[i, 2] - pend * coorreal[i, 1]
				cx = coorreal[i, 1] - 3 * (coorvirt[j, 1] - coorreal[i, 1])
				cy = const + pend * cx
				dummy.x = c(dummy.x, cx)
				dummy.y = c(dummy.y, cy)
				break
			}
		}

		if (sum(Joint[i, ]) == 1) {
			# When a real point is connected to just another real point, we need two virtual points on the non used lines
			# The procedure to calculate the dummy points is as before

			n3[i] = which(Joint[i, ] == 1)
			neighbors = which(Joint[i, ] == 1)
			A = intersect(IndReal[i, ], IndReal[neighbors[1], ])
			C = setdiff(IndReal[i, ], A)
			D = c(C, A[1])
			E = c(C, A[2])
			IndDummy = rbind(IndDummy, sort(D))
			IndDummy = rbind(IndDummy, sort(E))

			for (j in 1:nv) if (length(intersect(D, IndVirt[j, ])) == 2) {
				ndummy = ndummy + 1
				ColJointDummy = matrix(0, nr, 1)
				ColJointDummy[i] = 1
				JointDummy = cbind(JointDummy, ColJointDummy)
				n1[i] = -1 * ndummy
				pend = (coorvirt[j, 2] - coorreal[i, 2])/(coorvirt[j, 1] - coorreal[i, 1])
				const = coorreal[i, 2] - pend * coorreal[i, 1]
				cx = coorreal[i, 1] - 3 * (coorvirt[j, 1] - coorreal[i, 1])
				cy = const + pend * cx
				dummy.x = c(dummy.x, cx)
				dummy.y = c(dummy.y, cy)
				break
			}
			for (j in 1:nv) if (length(intersect(E, IndVirt[j, ])) == 2) {
				ndummy = ndummy + 1
				ColJointDummy = matrix(0, nr, 1)
				ColJointDummy[i] = 1
				JointDummy = cbind(JointDummy, ColJointDummy)
				n2[i] = -1 * ndummy
				pend = (coorvirt[j, 2] - coorreal[i, 2])/(coorvirt[j, 1] - coorreal[i, 1])
				const = coorreal[i, 2] - pend * coorreal[i, 1]
				cx = coorreal[i, 1] - 3 * (coorvirt[j, 1] - coorreal[i, 1])
				cy = const + pend * cx
				dummy.x = c(dummy.x, cx)
				dummy.y = c(dummy.y, cy)
				break
			}
		}
		if (sum(Joint[i, ]) == 0) {

			 IndDummy = matrix(0,3,2)
	     IndDummy[1,] = c(IndReal[1,1],IndReal[1,2])
	     IndDummy[2,] = c(IndReal[1,1],IndReal[1,3])
       IndDummy[3,] = c(IndReal[1,2],IndReal[1,3])
			 for(i in 1:3){
			   if(i == 1){n1[1] = -1}
			   if(i == 2){n2[1] = -2}
			   if(i == 3){n3[1] = -3}
			   xnTR = cbind(matrix(1, 1, 1), xn[1,2] + 2, a[IndDummy[i,1],IndDummy[i,2]] + b[IndDummy[i,1],IndDummy[i,2]]*(xn[1,2] + 2))
			   xnTL = cbind(matrix(1, 1, 1), xn[1,2] - 2, a[IndDummy[i,1],IndDummy[i,2]] + b[IndDummy[i,1],IndDummy[i,2]]*(xn[1,2] - 2))
			   
      	 pisubiR = matrix(0, 1, ngrup)
    		 denom = 1
    		 for (j in 1:(ngrup - 1)) denom = denom + exp(sum(beta[j, ] * xnTR[1, ]))
    		 for (j in 1:(ngrup - 1)) pisubiR[1, (j + 1)] = exp(sum(beta[j, ] * xnTR[1, ]))/denom
    		 pisubiR[1, 1] = 1/denom

       	 pisubiL = matrix(0, 1, ngrup)
    		 denom = 1
    		 for (j in 1:(ngrup - 1)){
             denom = denom + exp(sum(beta[j, ] * xnTL[1, ]))
         }
     
    		 for (j in 1:(ngrup - 1)) pisubiL[1, (j + 1)] = exp(sum(beta[j, ] * xnTL[1, ]))/denom               
    		 pisubiL[1, 1] = 1/denom
         if ((pisubiL[1,IndDummy[i,1]] > pisubiR[1,IndDummy[i,1]]) &
                (pisubiL[1,IndDummy[i,2]] > pisubiR[1,IndDummy[i,2]])) {
  				cx = xn[1,2] - 2
  				cy = a[IndDummy[i,1],IndDummy[i,2]] + b[IndDummy[i,1],IndDummy[i,2]]*(xn[1,2] - 2)
  			 }else{
  				cx = xn[1,2] + 2
  				cy = a[IndDummy[i,1],IndDummy[i,2]] + b[IndDummy[i,1],IndDummy[i,2]]*(xn[1,2] + 2)
  			 }
  		   dummy.x = c(dummy.x, cx)
  			 dummy.y = c(dummy.y, cy)
			 }
			 ndummy = 3
		}
	}

  if((ngrup > 3)&(nr > 1)){
  	Borders = matrix(0, sum(sum(Joint))/2, 2)

   	l = 0
  	for (i in 1:(nr - 1)) for (j in (i + 1):nr) if (Joint[i, j] == 1) {
  		l = l + 1
  		Borders[l, ] = sort(intersect(IndReal[i, ], IndReal[j, ]))
  	}

  	Borders = rbind(Borders, IndDummy)
  	nborders = sum(sum(Joint))/2 + ndummy
	}else if((ngrup == 3) | (nr == 1)){
	   Borders = matrix(0,3,2)
	   Borders[1,] = c(IndReal[1,1],IndReal[1,2])
	   Borders[2,] = c(IndReal[1,1],IndReal[1,3])
     Borders[3,] = c(IndReal[1,2],IndReal[1,3])
	   nborders = 3
	 }

  hc = HideCategories(IndReal,ngrup,nreal)
  hc_tras = t(hc)
  nhide = sum(hc_tras)
  ngrupvisible = ngrup - nhide

  equivRegiones = matrix(0,2,ngrup)
  contVisible = 1
  for(i in 1:ngrup){
    equivRegiones[1,i] = i
    if(hc[i,1] == 0){
      equivRegiones[2,i] = contVisible
      contVisible = contVisible + 1
    }
  }
  BordersWH = Borders
  for(i in 1:nborders){
     l = Borders[i, 1]
     m = Borders[i, 2]
     ind_l = which(equivRegiones[1,]==l)
     ind_m = which(equivRegiones[1,]==m)
     BordersWH[i, 1] = equivRegiones[2,ind_l]
     BordersWH[i, 2] = equivRegiones[2,ind_m]
  }

  Coord = InvertTesselation(beta,Borders,BordersWH,nborders,ngrupvisible)

	Centres = matrix(0, ngrupvisible, 2)
	for (i in 1:ngrupvisible) {
		Centres[i, 1] = Coord[2 * i - 1]
		Centres[i, 2] = Coord[2 * i]
	}

 if(nr==1)
  {
    bP12 = (IntersectSegments(Coord[1,1],Coord[2,1],Coord[3,1],Coord[4,1],coorreal[1,1],coorreal[1,2],dummy.x[1],dummy.y[1])
     | IntersectSegments(Coord[1,1],Coord[2,1],Coord[3,1],Coord[4,1],coorreal[1,1],coorreal[1,2],dummy.x[2],dummy.y[2])
     | IntersectSegments(Coord[1,1],Coord[2,1],Coord[3,1],Coord[4,1],coorreal[1,1],coorreal[1,2],dummy.x[3],dummy.y[3]))
    bP13 = (IntersectSegments(Coord[1,1],Coord[2,1],Coord[5,1],Coord[6,1],coorreal[1,1],coorreal[1,2],dummy.x[1],dummy.y[1])
         | IntersectSegments(Coord[1,1],Coord[2,1],Coord[5,1],Coord[6,1],coorreal[1,1],coorreal[1,2],dummy.x[2],dummy.y[2])
         | IntersectSegments(Coord[1,1],Coord[2,1],Coord[5,1],Coord[6,1],coorreal[1,1],coorreal[1,2],dummy.x[3],dummy.y[3]))
    bP23 = (IntersectSegments(Coord[3,1],Coord[4,1],Coord[5,1],Coord[6,1],coorreal[1,1],coorreal[1,2],dummy.x[1],dummy.y[1])
         | IntersectSegments(Coord[3,1],Coord[4,1],Coord[5,1],Coord[6,1],coorreal[1,1],coorreal[1,2],dummy.x[2],dummy.y[2])
         | IntersectSegments(Coord[3,1],Coord[4,1],Coord[5,1],Coord[6,1],coorreal[1,1],coorreal[1,2],dummy.x[3],dummy.y[3]))
    if((bP12 == FALSE) | (bP13 == FALSE) | (bP23 == FALSE)){
        CentresF = matrix(0, ngrupvisible, 2)
       	for (i in 1:ngrupvisible) {
          CentresF[i,1]= 2*coorreal[1,1] - Centres[i, 1]
          CentresF[i,2]= 2*coorreal[1,2] - Centres[i, 2]
        }
       bP12n = (IntersectSegments(CentresF[1,1],CentresF[1,2],CentresF[2,1],CentresF[2,2],coorreal[1,1],coorreal[1,2],dummy.x[1],dummy.y[1])
         | IntersectSegments(CentresF[1,1],CentresF[1,2],CentresF[2,1],CentresF[2,2],coorreal[1,1],coorreal[1,2],dummy.x[2],dummy.y[2])
         | IntersectSegments(CentresF[1,1],CentresF[1,2],CentresF[2,1],CentresF[2,2],coorreal[1,1],coorreal[1,2],dummy.x[3],dummy.y[3]))
        bP13n = (IntersectSegments(CentresF[1,1],CentresF[1,2],CentresF[3,1],CentresF[3,2],coorreal[1,1],coorreal[1,2],dummy.x[1],dummy.y[1])
             | IntersectSegments(CentresF[1,1],CentresF[1,2],CentresF[3,1],CentresF[3,2],coorreal[1,1],coorreal[1,2],dummy.x[2],dummy.y[2])
             | IntersectSegments(CentresF[1,1],CentresF[1,2],CentresF[3,1],CentresF[3,2],coorreal[1,1],coorreal[1,2],dummy.x[3],dummy.y[3]))
        bP23n = (IntersectSegments(CentresF[2,1],CentresF[2,2],CentresF[3,1],CentresF[3,2],coorreal[1,1],coorreal[1,2],dummy.x[1],dummy.y[1])
             | IntersectSegments(CentresF[2,1],CentresF[2,2],CentresF[3,1],CentresF[3,2],coorreal[1,1],coorreal[1,2],dummy.x[2],dummy.y[2])
             | IntersectSegments(CentresF[2,1],CentresF[2,2],CentresF[3,1],CentresF[3,2],coorreal[1,1],coorreal[1,2],dummy.x[3],dummy.y[3]))
       if((bP12n == TRUE) & (bP13n == TRUE) & (bP23n == TRUE)){
           Centres[,1]=CentresF[,1]
           Centres[,2]=CentresF[,2]
        }
    }
    if(sqrt((coorreal[1,1]-Centres[i,1])^2+(coorreal[1,2]-Centres[i,2])^2) < 0.3){  #Si la distancia entre el punto real y los de la inversion es pequeña
      CentresC = matrix(0, ngrupvisible, 2)
      for(i in 1:ngrupvisible){ 
          pte = (Centres[i,2]-coorreal[1,2])/(Centres[i,1]-coorreal[1,1])
          oo = coorreal[1,2] - coorreal[1,1]*pte
          a = pte*pte + 1
          b = -2*coorreal[1,1]+2*pte*(oo-coorreal[1,2])
          c=coorreal[1,1]^2 + oo^2 + coorreal[1,2]^2 -2*coorreal[1,2]*oo - 1
          xs = Eq2gSolve(a,b,c)
          inSegment = ((Centres[i,1] > coorreal[1,1])&(Centres[i,1] < xs[1,1]))
          if(inSegment == TRUE) {
               CentresC[i,1] = xs[1,1]
               CentresC[i,2] = xs[1,1]*pte + oo
          }else{
               CentresC[i,1] = xs[1,2]
               CentresC[i,2] = xs[1,2]*pte + oo
          }
      }
      Centres = CentresC
    }

  }

	result = list()
	result$x = coorreal[, 1]
	result$y = coorreal[, 2]
	result$node = rep(TRUE, nr)
	result$n1 = n1
	result$n2 = n2
	result$n3 = n3
	result$dummy.x = dummy.x
	result$dummy.y = dummy.y
	result$ndummy = ndummy
	result$IndReal = IndReal
	result$Centers = Centres
	result$hideCat = hc_tras
	result$equivRegiones = equivRegiones
	class(result) <- "voronoiprob"
	return(result)
}
