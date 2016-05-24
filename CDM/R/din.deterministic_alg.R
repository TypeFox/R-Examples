
##############################################################
# link to Rcpp functions
din.deterministic.devcrit <- function( dat , datresp , latresp , guess , slip ){
	res <- .Call("din_deterministic_devcrit_C", 
					dat , datresp , latresp , guess , slip , 
					PACKAGE = "CDM")
    return(res)			
			}
#**********			
# JML estimation function			
din.jml.devcrit <- function( dat , datresp , latresp , guess , slip ){
	res <- .Call("din_jml_devcrit_C", 
					dat , datresp , latresp , guess , slip , 
					PACKAGE = "CDM")
    return(res)			
			}			
#################################################################
# define different attribute pattern for dichotomous attributes
define.attribute.space <- function( q.matrix ){
    nodes <- c(0,1)
    K <- ncol(q.matrix)
    attr.patt <- as.matrix( expand.grid( as.data.frame( matrix( rep(nodes, K) , ncol = K ) ) ) )    
	if ( ! is.null( colnames(q.matrix) ) ){
		colnames(attr.patt) <- colnames(q.matrix)
						}
    return(attr.patt)
        }
#################################################################		
# compute latent responses
compute.latent.response <- function( attr.patt , q.matrix ,rule=NULL){
	AP <- nrow(attr.patt)
	I <- nrow(q.matrix)
	if ( is.null(rule) ){ rule <- rep("DINA" , I ) }
	latresp <- matrix( NA , nrow=AP , ncol=I)
	colnames(latresp) <- rownames(q.matrix)
	for (ii in 1:I){
		# ii <- 1
		comp.ii <- 1
		if ( rule[ii] == "DINA"){ comp.ii <- sum(q.matrix[ii,] ) }
		latresp[,ii] <- 1 * ( attr.patt %*% q.matrix[ii,] >= comp.ii )
			}
	return( as.matrix(latresp) )
		}
#####################################################################
# calculate deviation criterion
# loop over attributes
calc.devcrit <- function( dat , dat.resp , latresp ,slip , guess , N , I , AP){
    dev.crit <- matrix( NA , nrow=N , ncol=AP )
    for (aa in 1:AP){
        # aa <- 4
        lat.aa <- matrix( latresp[aa,] , nrow=N , ncol=I , byrow=TRUE)
        dev.crit[,aa] <- rowSums( slip * ( lat.aa - dat ) * ( lat.aa == 1 ) * dat.resp + 
                    guess * ( dat - lat.aa ) * ( lat.aa == 0 ) * dat.resp )
                    }
    return(dev.crit)
            }
#######################################################################