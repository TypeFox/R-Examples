

#################################################
# load packages or install some packages
# if they are needed.
library_install <- function( pkg , ... ){
	PP <- length(pkg)
	for (vv in 1:PP){
		pp <- pkg[vv]
		ab1 <- do.call( require , list( package=pp , quietly=TRUE ) )
		if ( ! ab1 ){
			install.packages( pp , ...)
				   }
		do.call( require , list( package=pp ) )
					}
            }
################################################			