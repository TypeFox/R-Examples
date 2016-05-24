#'@importFrom rgl addNormals
#'@importFrom rgl shade3d
#'@importFrom rgl scale3d
#'@importFrom rgl cylinder3d
#'@importFrom rgl translate3d
#'@importFrom rgl rotate3d
plot3dcylinder <- function(start=c(0,0,0), end=c(0,0,1), radius, mycol="chocolate4", ... ){ 

	# ... cylinder as basis-vector 
	cyl <-  cylinder3d(rbind(c(0,0,0), #start, # start 
							 c(0,0,1)), # end 
							 radius = radius, 
							 e1=cbind(0, 0, 1), 
							 e2=cbind(1, 0, 0), 
							 sides=10 
						) 

	# ... rotate cylinder horizontally and scale it 
	# len <- sqrt(sum(abs(end)*abs(end))) 
	len <- sqrt(sum((end-start)^2))
	cyl <- scale3d(cyl,1,1,len) 

	# Rotate cylinder.
	rotation <- GramSchmidt(end-start+c(1,0,0),end-start+c(0,1,0),end-start,   
		order=c(3,1,2)) 
	cyl <- rotate3d(cyl, matrix=rotation) 

	# Move cylinder:
	cyl <- translate3d(cyl, x=start[1], y=start[2], z=start[3])
	
	# plot cylinder.
	shade3d(addNormals(cyl), col=mycol) 
	
} 
