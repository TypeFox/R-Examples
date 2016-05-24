###################################################
# create delta design matrix
create.deltaDesign <- function(X , main=2 , int=TRUE ){
	D <- ncol(X)
	theta <- X
	deltaD <- matrix( 1 , nrow=nrow(theta) , ncol=1)
#	deltaD <- matrix( 1 , nrow=nrow(theta) , ncol=0)
	for (dd in 1:D){
		for (mm in 1:main){
			deltaD <- cbind( deltaD , theta[,dd]^mm )
							}
						}
	if ( D>1){
		if (int ){
			for (dd1 in 1:(D-1) ){
				for (dd2 in (dd1+1):D){
					deltaD <- cbind( deltaD , theta[,dd1]*theta[,dd2] )			
								}
							}
					}
				}
	return(deltaD)
		}
############################################################