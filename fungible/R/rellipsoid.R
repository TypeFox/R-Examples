#####################################################
# rellipsoid: R function to sample points uniformly #
# distributed on the surface of an n-dimensionsal   #
# ellipsoid                                         #
#                                                   #
# 9/18/09                                           #
#####################################################


rellipsoid <- function(R, Rsq,Npoints) {

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Arguments~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Arguments:
#	R            Symmetric PD matrix in which b'Rb = Rsq
#	Rsq		     R-squared (squared coefficient of determination)
#  Npoints      Number of desired points in n-space  
#
#  Value:
#  b            Npoints sets of regression vectors 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# compute eigenvalues of R
  Lambda <- eigen(R)$values   
  if(min(Lambda) <=0) stop("\n *** R is not Positive Definite ***\n")

# number of columns of R
  DIM <- ncol(R)              

# compute lengths of the ellipsoid semi-axes
  semi.axes <- sqrt(1/Lambda)	
# Normalize the ellipsoid to have a smallest axis of 1
	ai <- as.matrix(semi.axes/min(semi.axes))
   
# Sample until Npoints accepted
 points.chosen <- 1
 trials <-0
 ellipsoid.points <- matrix(0,DIM,Npoints)
 
 while(points.chosen <=Npoints){
	
	  xyz <- rnorm(DIM)
	       
# Generate uniformally distributed points on an n-spheroid
# Method 2 of Marsaglia 1972
	  sphere.points <- xyz/ sqrt(sum(xyz^2))
	
# Select points via Acceptance/Rejection method 

	if( sum( (sphere.points/ai)^2 )  >= runif(1)^2) {
	# stretch to ellipsoid surface
		ellipsoid.points[,points.chosen]<-sphere.points * semi.axes
		points.chosen <- points.chosen + 1
	}
	
	trials <- trials +1
  } #end while 
		
	acceptance <- Npoints/trials
	cat("\n",Npoints, "points generated on the surface\n", 	"of a", DIM,"dimensional ellipsoid.\n Acceptance rate: ",
	round(acceptance,2),"\n\n")
		
   #back transform from standard position and scale for Rsq
	b <- eigen(R)$vectors %*% ellipsoid.points * sqrt(Rsq)
	list(b=b)
} # End function rellipsoid

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# # R <- matrix(c(1, .56, .56, 1),2,2)
# 
# 
# # Correlations derived from Table A.7 page 224
# # from WAIS-III manual
# # The Psychological Corporation (1997). WAIS-III
# # WMS-III Technical Manual.
# # San Antonio, TX: Author.
# 
# wais3 <- R <- matrix(
# c(1, .76, .58, .43, .75, .75, .42, .54, .41, .57, .64, .54, .50, .53,
# .76, 1, .57, .36, .69, .71, .45, .52, .36, .63, .68, .51, .47, .54,
# .58, .57, 1, .45, .65, .60, .47, .48, .43, .59, .60, .49, .56, .47,
# .43, .36, .45, 1, .37, .40, .60, .30, .32, .34, .35, .28, .35, .29,
# .75, .69, .65, .37, 1, .70, .44, .54, .34, .59, .62, .54, .45, .50,
# .75, .71, .60, .40, .70, 1, .42, .51, .44, .53, .60, .50, .52, .44,
# .42, .45, .47, .60, .44, .42, 1, .46, .49, .47, .43, .27, .50, .42,
# .54, .52, .48, .30, .54, .51, .46, 1, .45, .50, .58, .55, .53, .56,
# .41, .36, .43, .32, .34, .44, .49, .45, 1, .47, .49, .41, .70, .38,
# .57, .63, .59, .34, .59, .53, .47, .50, .47, 1, .63, .62, .58, .66,
# .64, .68, .60, .35, .62, .60, .43, .58, .49, .63, 1, .59, .50, .59,
# .54, .51, .49, .28, .54, .50, .27, .55, .41, .62, .59, 1, .48, .53,
# .50, .47, .56, .35, .45, .52, .50, .53, .70, .58, .50, .48, 1, .51,
# .53, .54, .47, .29, .50, .44, .42, .56, .38, .66, .59, .53, .51, 1),
# nrow=14,ncol=14)
# 	
# 	
# # # generate uniformally distributed regression vectors
# # # on the surface of a 14-dimensional ellipsoid 
# #   N <- 1000000
# #   Rsq <- .21
# #   R <- wais3[1:6,1:6]             #[1:3,1:3]
# #   b<-rellipsoid(R,Rsq, Npoints=N)
# # 
# #  plot(b[1,],b[2,])
# #   
# # #compute validity vectors
# #   r <- R %*% b
# # 	
# # 	Rsq.r <- Rsq.unit <- rep(0,N)
# # 	for(i in 1:N){
# # 		# performance of unit weights
# # 		Rsq.unit[i] <- (t(sign(r[,i])) %*% r[,i])^2 /
# # 		               (t(sign(r[,i])) %*% R %*% sign(r[,i]))
# # 		# performance of correlation weights               
# #      	Rsq.r[i] <- (t(r[,i]) %*% r[,i])^2 /(t(r[,i]) %*% R %*% r[,i])	}
# # 
# #    cat("\nAverage relative performance of unit weights across all criteria:",
# # 	    round(mean(Rsq.unit)/Rsq,3) )     
# # 	cat("\n\nAverage relative performance of r weights across all criteria:",
# # 	    round(mean(Rsq.r)/Rsq,3) ) 
# # 	
# # 	# focus on sector with all positive validity weights    
# # 	all.positive <- apply(sign(r),2,sum)==ncol(R)
# # 	cat("\n\nAverage relative performance of unit weights across all criteria with positive validity weights:",
# # 	    round(mean(Rsq.unit[all.positive])/Rsq,3) )     
# # 	cat("\n\nAverage relative performance of r weights across all criteria with positive validity weights:",
# # 	    round(mean(Rsq.r[all.positive])/Rsq,3) ) 
# # 	
# # 	
