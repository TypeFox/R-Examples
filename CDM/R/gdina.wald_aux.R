
##############################################
##############################################
# calculate restricted model and output
# distance to estimated GDINA model
calc_dist_restricted_model <- function(pjj , Mj.ii0 , attr.prob.ii , link,
           suffstat_probs.ii, aggr.attr.patt.ii  ){				
		mod0.ii <- stats::lm( pjj ~ 0 + Mj.ii0 , weights = attr.prob.ii )
		fitted.ii <- stats::fitted(mod0.ii)
		if ( link == "logit"){ fitted.ii <- stats::plogis(fitted.ii) }
		if ( link == "log"){ fitted.ii <- exp(fitted.ii) }		
		# weighted distance
		dist.ii <- sum( attr.prob.ii * ( suffstat_probs.ii - fitted.ii)^2 )
		# unweighted distance
		d1 <- ( suffstat_probs.ii - fitted.ii)^2
		d1 <- as.vector(d1[ aggr.attr.patt.ii ])
		d1 <- mean( d1 )
		res <- list( "wgtdist" = dist.ii , "uwgtdist" = d1)
		return(res)
					}
###################################################
###################################################
# creation of constraint matrix
contraint_matrix <- function( delta.ii , rule, Kii, Mj.ii ){
	Dii <- length(delta.ii)
    #***************************
	# DINA
	if ( rule == "DINA"){
			R <- matrix( 0 , nrow=Dii-2 , ncol=Dii)
			for (vv in 2:(Dii-1) ){    
			     R[vv-1,vv] <- 1 
				      }
					  }
    #***************************
	# ACDM
	if (rule == "ACDM"){
			R <- matrix( 0 , nrow=Dii-(Kii+1) , ncol=Dii)
			for (vv in 1:(nrow(R)) ){  
			     vv1 <- vv + ( Kii +1 )
			     R[vv,vv1] <- 1 
							}						  
						}
    #***************************
	# DINO
	if (rule=="DINO"){
	       Dii <- ncol(Mj.ii)
		   R <- matrix( 0 , nrow=Dii-2 , ncol=Dii)	
		   for (vv in 1:(Dii-2)){
				R[vv,] <- Mj.ii[vv+2,] - Mj.ii[vv+1 , ]
					}	
					}
	#**********************************				
		return(R)
			}
####################################################			
