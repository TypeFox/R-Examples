

##################################################
# linear model for clustered data
lm.cluster <- function( data , formula , cluster , ... ){
	# requireNamespace("multiwayvcov")
	mod <- stats::lm( data=data , formula=formula ,  ... )
	if ( length(cluster) > 1 ){
			v1 <- cluster 
				} else {
			v1 <- data[,cluster]
						}	
	dfr <- data.frame( cluster = v1 ) 
	vcov2 <- multiwayvcov::cluster.vcov( model = mod , cluster = dfr)	
	res <- list( "lm_res" = mod , "vcov" = vcov2 )
	class(res) <- "lm.cluster"
	return(res)
			}
###################################################			
coef.lm.cluster <- function( object , ... ){
	coef( object$lm_res)
			}
####################################################			
vcov.lm.cluster <- function( object , ... ){
	 object$vcov
			}
####################################################
summary.lm.cluster <- function( object , ... ){
	smod <- summary( object$lm_res )
	csmod <- smod$coefficients
	csmod[,"Std. Error"] <- sqrt( diag( vcov(object) ))
	csmod[,"t value"] <-  csmod[,"Estimate"] / csmod[,"Std. Error"]
	csmod[,"Pr(>|t|)"] <- stats::pnorm( - abs( csmod[,"t value"] ) )*2
	R2 <- smod$r.squared
	cat("R^2 =" , round(R2 , 5),"\n\n" )
	print(csmod)
	invisible(csmod)
			}
#######################################################			