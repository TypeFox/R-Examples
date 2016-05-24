region.P <-
function(region,P,k=100){
			vert0		=cartesian2polar2D(region,P)
			vert0.theta=vert0[,2]

			ind.theta	=sort.list(vert0.theta)
			vert0		=vert0[ind.theta,]
			n.vert0	=nrow(vert0)

## vert0 now contains the n.vert0 polar coordinates of the n.vert0 vertex of "region" (centered on "P")
## (rho, theta) ordered by positive values of theta

			vert		=matrix(0,length(vert0.theta)+1,2)
		vert[1:n.vert0,]=vert0
		vert[n.vert0+1,]=vert0[1,]
			vert.theta	=vert[,2]
			vert.rho	=vert[,1]
			n.vert	=nrow(vert)
			vert.seq	=1:(n.vert-1)
## vert is the cyclic version of vert0; vert.rho and vert.theta are its columns


			vert.gamma	=angle.positive(diff(vert.theta))
## vert.gamma: angles in P of the triangle formed by P and two consecutive vertices
## this generic triangle  T_i has vertices: V_i,P,V_{i+1}
			vert.alpha	=triangle.abgamma.alpha(vert.rho[vert.seq+1],vert.rho[vert.seq],vert.gamma)
## vert.alpha is the angle in  V_i of the triangle T_i
			theta		=angle.positive(seq(0,2*pi,length.out=k+1)+vert.theta[1]-2*pi)
			theta		=theta[2:(k+1)]
## theta is the vector of the k angles in P starting from V_1
			ind1		=findInterval(theta,vert.theta[vert.seq])
			ind1		=ifelse(ind1==0,max(vert.seq),ind1)
## ind1[j]=i if the j-th segment is inside the i-th triangle
			alpha.inf	=vert.alpha[ind1]
			theta.inf	=angle.positive(theta-vert.theta[ind1])
			rho.inf	=vert.rho[ind1]
			rho		=triangle.alphabetac.b(theta.inf,alpha.inf,rho.inf)
## alpha.inf_j,theta.inf_j,rho_j are two angles and the edge of a triangle S_j, with vertex in P and V_{ind_j} and angle in P equal to theta.inf_j
## (triangle S_j and T_{ind_j} has a common edge) 
			return(list(vert.theta=vert.theta,vert.rho=vert.rho,
					vert.gamma=vert.gamma,vert.alpha=vert.alpha,
					theta=theta,vert.seq=vert.seq,theta.inf=theta.inf,
					rho.inf=rho.inf,alpha.inf=alpha.inf,rho=rho,ind1=ind1))
				}
