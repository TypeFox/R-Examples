

##################################################
# linear model for clustered data
glm.cluster <- function( data , formula , cluster , ... ){
	mod <- stats::glm( data=data , formula=formula ,  ... )
	if ( length(cluster) > 1 ){
			v1 <- cluster 
				} else {
			v1 <- data[,cluster]
						}	
	dfr <- data.frame( cluster = v1 ) 
	vcov2 <- multiwayvcov::cluster.vcov( model = mod , cluster = dfr)	
	res <- list( "glm_res" = mod , "vcov" = vcov2 )
	class(res) <- "glm.cluster"
	return(res)
			}
###################################################			
coef.glm.cluster <- function( object , ... ){
	coef( object$glm_res)
			}
####################################################			
vcov.glm.cluster <- function( object , ... ){
	 object$vcov
			}
####################################################
summary.glm.cluster <- function( object , ... ){
	smod <- summary( object$glm_res )
	csmod <- smod$coefficients
	csmod[,"Std. Error"] <- sqrt( diag( vcov(object) ))
	csmod[,"z value"] <-  csmod[,"Estimate"] / csmod[,"Std. Error"]
	csmod[,"Pr(>|z|)"] <- stats::pnorm( - abs( csmod[,"z value"] ) )*2
	# R2 <- smod$r.squared
	# cat("R^2 =" , round(R2 , 5),"\n\n" )
	print(csmod)
	invisible(csmod)
			}
#######################################################			