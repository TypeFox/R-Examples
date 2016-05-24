
#######################################################
# mcdina estimate item parameters
mcdina.est.item <- function( n.ik , lr_list , lc_list , delta , I , G , 
		eps , itemstat , itempars , lr_counts ){
		
		for (ii in 1:I){   #		ii <- 2	
			lr.ii <- lr_list[[ii]]
			lr_index.ii <- lr.ii[ , "lr_index"]
			lc.ii <- lc_list[[ii]]
			#******************
			# group specific item parameters
			if ( itempars[ii] == "gr"){	
				for (gg in 1:G){ # gg <- 1
					n.ik.ii.gg <- t( n.ik[ii,,,gg] )
					n.ik.ii.gg_aggr <- as.matrix( rowsum( n.ik.ii.gg , lr_index.ii   ) )
					n.ik.ii.gg_aggr1 <- n.ik.ii.gg_aggr+eps
					cn <- rowSums( n.ik.ii.gg_aggr1 )					
					lr_counts[ii , 1:itemstat$N.lr[ii] ,  gg ] <- cn
					delta.new <- t( n.ik.ii.gg_aggr / 
							matrix( cn  , nrow=nrow(n.ik.ii.gg_aggr) , 
							ncol=ncol(n.ik.ii.gg_aggr) , byrow=FALSE ) )
					delta[ii,,1:itemstat[ii,"N.lr"],gg] <- delta.new
							}
							}
			#******************
			# group specific item parameters
			if ( itempars[ii] != "gr"){	
					n.ik.ii.gg <- t( n.ik[ii,,,1] )

					for (gg in 2:G){ 
					   n.ik.ii.gg <- n.ik.ii.gg + t( n.ik[ii,,,1] )
									}
					n.ik.ii.gg_aggr <- as.matrix( rowsum( n.ik.ii.gg , lr_index.ii   ) )
					n.ik.ii.gg_aggr1 <- n.ik.ii.gg_aggr+eps
					cn <- rowSums( n.ik.ii.gg_aggr1 )
					lr_counts[ii , 1:itemstat$N.lr[ii] ,  1:G ] <- cn
					delta.new <- t( n.ik.ii.gg_aggr / 
							matrix( cn  , nrow=nrow(n.ik.ii.gg_aggr) , 
							ncol=ncol(n.ik.ii.gg_aggr) , byrow=FALSE ) )
					for (gg in 1:G){
						delta[ii,,1:itemstat[ii,"N.lr"],gg] <- delta.new
									}
							}  # end itempars[ii] != "gr"							
				}  # end ii						
		res <- list("delta"= delta , "lr_counts" = lr_counts )
		return(res)
			}
############################################################