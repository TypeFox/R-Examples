
###################################################################
# Input is a mirt syntax specified with variable names instead
# of numbers
mirt.model.vars <- function( mirtsyn2 , data=NULL , items=NULL ){

	# mirt syntax splitted
	mirtsyn2 <- gsub( ";" , "\n" , mirtsyn2 )
	msyn0 <- strsplit( mirtsyn2 , c(" ") )[[1]]
	syn <- msyn0
	if ( !is.null(data) ){
		items <- colnames(data) 
					}
	if ( is.null(items) ){
		stop("Provide dataset or items.\n")
						}
	# admissible strings
	vecstr <- c( "\n" , "\\(" , "=" , "-" , "," , "\\)" )
	#*******
	# process syntax
	for (vv in vecstr){
    syn <- split_syn_string( syn , vv  )    
            }
    # postprocess syntax                                    
	syn <- syn[ syn != "" ]
	syn[ syn == "\\(" ] <- "("
	syn[ syn == "\\)" ] <- ")"
	# substitute variables by numbers
	VV <- length(items)
	useditems <- NULL
	for (vv in 1:VV){
		ind <- which( syn == items[vv] )
		if (length(ind) > 0 ){
			syn[ ind ] <- vv 
			useditems <- c( useditems , items[vv] )   
						}               
					}
	syn <- paste0( syn , collapse="")
	mirtmodel <- mirt::mirt.model(syn)
	mirtmodel <- list("model"=mirtmodel , "syntax" = syn )
	return(mirtmodel)
		}
#########################################################################		
# split_syn_string vectorized input
split_syn_string_vec <- function( syn, vecstr ){
	for (vv in vecstr){
    syn <- split_syn_string( syn , vv  )    
            }
	return(syn)
		}
########################################################################


########################################################################
# cleans syntax in a vector from strings vv
split_syn_string <- function( syn , vv ){
	syn <- as.list(syn )
	syn.vv <- grep( vv , syn )
	LL <- length(syn.vv)
	if (LL>0){
		for (ii in 1:LL){
			ll <- syn.vv[ii]
 			syn.ll <- syn[[ll]]		
			syn[[ll]] <- split_conc( syn.ll , vv )
							}
					}
	syn <- unlist(syn)
	return(syn)
			}
########################################################################

########################################################################
# splits a string syn.ll and concatanates it with string vv
split_conc <- function( syn.ll , vv ){
	g1 <- strsplit( syn.ll , vv , perl=FALSE )[[1]] 
	Lg1 <- length(g1)
	vec <- NULL
	if (Lg1 == 1 ){ vec <- c( g1 , vv ) }
	if (Lg1 > 1 ){
		vec <- rep("" , Lg1 + (Lg1-1) )
		vec[ seq( 1 , 2*Lg1 , 2 ) ] <- g1
		vec[ seq( 2 , 2*Lg1 - 1 , 2 ) ] <- vv	
		Ls1 <- nchar(syn.ll)
		if ( substring( syn.ll , Ls1 , Ls1 ) == vv ){
		       vec <- c( vec , vv )
							}
			}	
	return(vec)
			}
########################################################################			