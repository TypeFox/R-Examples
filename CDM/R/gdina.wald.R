######################################################
# Wald tests at item level
gdina.wald <- function( object ){
	varmat.delta <- object$varmat.delta
	delta <- object$delta
	rule <- object$control$rule
	#*****
	# checks whether gdina.wald can be applied
	if ( mean(rule=="GDINA") < 1 ){
	   stop("Specify a full GDINA model for performing a Wald test.\n")
								}
	Mj <- object$control$Mj
	q.matrix <- object$q.matrix
	I <- nrow(q.matrix)
	dat <- object$dat
	stats_vars <- c("_X2" , "_df" , "_p", "_sig" ,  "_RMSEA" ,
			"_wgtdist" , "_uwgtdist")
	SV <- length(stats_vars)
	# number of rules
	cdm_rules <- c("DINA" , "DINO" , "ACDM")
	SR <- length(cdm_rules)
	
	
	stats <- matrix( NA , nrow= I  , ncol=SV*SR)
	v1 <- NULL
	for (ss in 1:SR){
		v1 <- c( v1 , paste0( cdm_rules[ss] , stats_vars ) )
			}	
	colnames(stats) <- v1
	
	#-----------------
	# extract relevant output
	suffstat_probs <- object$control$suffstat_probs
	Mj <- object$control$Mj
	Aj <- object$control$Aj
	attr.prob <- object$control$attr.prob
	aggr.attr.patt <- object$control$aggr.attr.patt
	link <- object$link
	
	#******************************************
	# loop over items
	for (ii in 1:I){
		delta.ii <- delta[[ii]]
		var.delta.ii <- varmat.delta[[ii]]
		# number of attributes
		Kii <- sum( q.matrix[ii,] )
		Mj.ii <- Mj[[ii]][[1]]
		suffstat_probs.ii <- suffstat_probs[[ii]]
		pjj <- suffstat_probs.ii
		if ( link == "logit"){ pjj <- stats::qlogis(pjj) }
		if ( link == "log"){ pjj <- log(pjj) }				
		aggr.attr.patt.ii <- aggr.attr.patt[[ii]]
		attr.prob.ii <- rowsum( attr.prob , aggr.attr.patt.ii )
		Aj.ii <- Aj[[ii]]
		
		if (Kii>1){
		
			nobs <- sum( 1 - is.na(dat[,ii] ) )
		    
		    #--------------------------------
			# Tests for different rules
			
			for (rule in cdm_rules){
				# rule <- "DINA"				
				R <- contraint_matrix( delta.ii , rule= rule, Kii, Mj.ii )								
				res <- WaldTest( delta.ii , var.delta.ii , R , nobs )
				stats[ii,paste0(rule,"_X2")] <- res$X2
				stats[ii,paste0(rule,"_df")] <- res$df
				stats[ii,paste0(rule,"_p")] <- res$p
				stats[ii,paste0(rule,"_RMSEA")] <- res$RMSEA
				Mj.ii0 <- .create.Mj( Aj.ii , rule)[[1]]
				res.ii <- calc_dist_restricted_model(pjj , Mj.ii0 , attr.prob.ii, link,
						  suffstat_probs.ii, aggr.attr.patt.ii	)						  
				stats[ii,paste0(rule,"_wgtdist")] <- res.ii$wgtdist
				stats[ii,paste0(rule,"_uwgtdist")] <- res.ii$uwgtdist
								}
			
			
				}  # end if Kii > 1       
		}  # end item

	stats <- data.frame( "item" = colnames(dat) , "NAttr" = rowSums(q.matrix) , stats )

	levels <- c(.01 , .05 )
	labels <- c("**" , "*")	
	for (rule in cdm_rules){
	stats[ ,paste0(rule,"_sig")] <- 
			label_significance_level( stats[,paste0(rule,"_p")] , levels , labels )		
							}
	res <- list("stats"=stats, "cdm_rules" = cdm_rules)
	class(res) <- "gdina.wald"
	return(res)
	}
####################################################	
	

