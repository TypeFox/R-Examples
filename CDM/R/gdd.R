#################################################################
# generalized distance discriminating method
gdd <- function( data , q.matrix , theta , b , a  , skillclasses=NULL){
	data <- as.matrix(data)
	dataresp <- as.matrix( 1 - is.na( data ) )
	data[ is.na(data) ] <- 0
	q.matrix <- as.matrix(q.matrix)
	skillspace <- skillclasses
	# compute ideal response pattern
	res <- ideal.response.pattern( q.matrix , skillspace )
	idealresp <- res$idealresp
	skillspace <- res$skillspace
	# apply generalized distance discriminating method
	res <- generalized_distance_method__Cpp( data , dataresp , 
	            idealresp, theta , a , b )
	# extract results
	distmatrix <- res$dist
	skillclass.est <- skillspace[ res$est_skill , ]
	res <- list( "skillclass.est" = skillspace , "distmatrix" = distmatrix ,
				  "skillspace" = skillspace , "theta" = theta )   
	return(res)
		}
###############################################################################
# auxiliary function for generalized distance discriminating method
#  SEXP generalized_distance_method__C( SEXP data_, SEXP dataresp_, SEXP idealresp_, 
#   SEXP theta_, SEXP a_, SEXP b_) ;
generalized_distance_method__Cpp <- 
function( data_ , dataresp_ , idealresp_ , theta_ , a_ , b_ ){
	res <- .Call("generalized_distance_method__C", 
					data_ , dataresp_ , idealresp_ , theta_ , a_ , b_ , 
					PACKAGE = "CDM")
    return(res)			
			}
###############################################################################