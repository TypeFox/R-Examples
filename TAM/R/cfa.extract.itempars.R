
######################################################################
# extract item parameters from fitted cfa object (in lavaan)
cfa.extract.itempars <- function( object ){

	if ( object@Options$model.type != "cfa" ){
			stop("Function can only be applied \n if cfa (in lavaan) is used.")
											}								
	# ParTable <- as.data.frame( object@ParTable )
	ParTable <- as.data.frame( lavaan::parameterTable(object) )
	ParTable$parname <- paste0( ParTable$lhs , ParTable$op , ParTable$rhs )	
	labels1 <- paste(ParTable$label)
	ParTable$parname <- ifelse ( labels1 != "" , labels1 , ParTable$parname )	
	ParTable$est <- ParTable$ustart
#	cobj <- coef( object )
#	vars1 <- intersect( names(cobj) , ParTable$parname )
#	NV <- length(vars1)
#	for ( vv in vars1){
#		ParTable[ ParTable$parname %in% vv , "est" ] <- cobj[ vv ]
#					}
	# extract sample statistics
	means <- object@SampleStats@mean[[1]]
	obs.vars <- object@Data@ov.names[[1]]
	covs <- object@SampleStats@cov[[1]]	
	colnames(covs) <- rownames(covs) <- obs.vars
	
	NOV <- length(obs.vars)

	# extract loadings, means and covariance matrix
	part <- ParTable	
	part1 <- part[ paste(part$op) == "=~" , ] 
	lat.vars <- unique( paste(part1$lhs ))
	NLV <- length(lat.vars)
	# loading matrix
	L <- matrix( 0 , nrow=NOV , ncol=NLV)
	rownames(L) <- obs.vars 
	colnames(L) <- lat.vars
	for (ll in lat.vars){
		# ll <- lat.vars[1]
		part1.ll <- part1[ part1$lhs == ll , ]
		L[ paste(part1.ll$rhs) , ll ] <- part1.ll$est
						}

	# covariance matrix of latent variables
	Sigma <- matrix( 0 , nrow=NLV , ncol= NLV )
	rownames(Sigma) <- colnames(Sigma) <- lat.vars
	# vector of latent variable means
	mu <- rep( 0 , NLV )
	names(mu) <- lat.vars
						
	# vector of intercepts
	nu <- rep( 0 , NOV)
	names(nu) <- obs.vars
	names(means) <- obs.vars
	nu[ names(means) ] <- means
	part1 <- part[ paste(part$op) == "~1" , ] 
	if ( nrow( part1) > 0 ){
		part1.ll <- part1[ part1$lhs %in% obs.vars , ]
		nu[ paste( part1.ll$lhs) ] <- part1.ll$est 
		part1.ll <- part1[ part1$lhs %in% lat.vars , ]
		mu[ paste( part1.ll$lhs) ] <- part1.ll$est 
							}
	# extract covariance matrices	
	part1 <- part[ paste(part$op) == "~~" , ] 		
	part1 <- part1[ paste( part1$lhs) %in% lat.vars , ]
    NL <- nrow(part1)
	for (ll in 1:NL){
	   Sigma[ paste(part1[ll,"rhs"]) , paste(part1[ll,"lhs"]) ] <- 
	      Sigma[ paste(part1[ll,"lhs"]) , paste(part1[ll,"rhs"]) ] <- part1[ll,"est"]
					}
	# extract residual variances
	psi <- matrix( 0 , nrow=NOV , ncol=NOV )
	colnames(psi) <- rownames(psi) <- obs.vars
	part1 <- part[ paste(part$op) == "~~" , ] 		
	part1 <- part1[ paste( part1$lhs) %in% obs.vars , ]
    NL <- nrow(part1)
	for (ll in 1:NL){
	   psi[ paste(part1[ll,"rhs"]) , paste(part1[ll,"lhs"]) ] <- 
	      psi[ paste(part1[ll,"lhs"]) , paste(part1[ll,"rhs"]) ] <- part1[ll,"est"]
					}
					
	#**** output				
	res <- list( "L"=L , "nu"= nu , "psi"= psi ,
			"Sigma" = Sigma , "mu" = mu ,
			"obs.means" = means , "obs.cov" = covs ,
			"obs.vars" = obs.vars , "lat.vars" = lat.vars ,
			"NOV" = NOV , "NLV" = NLV )

	return(res)
		}
######################################################