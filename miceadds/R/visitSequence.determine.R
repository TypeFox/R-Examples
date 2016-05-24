
########################################################
visitSequence.determine <- function( impMethod , vis ,
	data , maxit = 10){
	dat <- data
	ind <- grep( "~ I" , impMethod )
	I1 <- length(ind)
	# inits
	vis1 <- vis0 <- vis   
	# select passive variables
	pass_vars <- names(impMethod)[ind]
	iter <- 0
	while (iter < maxit ){
		vis <- vis1
		for (var.ii in pass_vars){
			vis1 <- handle.variable( var.ii , impMethod , vis1 , dat )
                    }
		if ( equal_vecs( vis , vis1 )  ){
			iter <- maxit + 1 
							}
		iter <- iter + 1
						}
	return(vis1)							
				}
########################################################				

#*********************************************+
# handle one passive variable
handle.variable <- function( var.ii , impMethod , vis1 , dat ){
    ness.ii <- colnames( stats::get_all_vars( 
				stats::as.formula( impMethod[ var.ii ] ) , data=dat[1,] ) )
    IN <- length(ness.ii)
    for (nn in 1:IN){
    #    nn <- 1
        i2 <- match( ness.ii[nn] , names(vis1) )
        for (oo in 1:length(i2)){
        vis1 <- add.entry( ness.ii.nn=ness.ii[nn] , var.ii , vis1 , occ=oo)
                    }
                    }
    vis1 <- remove.entry( vis1 , ness.ii , var.ii )
    return(vis1)
            }
#*************************************************				
# add missing entry	in visit sequence			
add.entry <- function( ness.ii.nn , var.ii , vis1 , occ=1 ){    
	index.ness.ii <- which( names(vis1) == ness.ii.nn )[occ]	
	index.var.ii <- which( names(vis1) %in% var.ii )	
	if ( ! ( (index.ness.ii+1) %in% index.var.ii ) ){ 	
        vis1 <- append( vis1 , vis1[ var.ii ][1] , index.ness.ii )
                        }												
    return(vis1)
        }				
#***************************************************		
# remove entry
remove.entry <- function( vis1 , ness.ii , var.ii ){
    index.ness.ii <- unlist( sapply( ness.ii , FUN = function(nn){ 
				match( nn , names(vis1) ) } ) )
    index.var.ii <- which( names(vis1) %in% var.ii )    
    IN <- length( index.var.ii)
    l1 <- NULL
    for (vv in 1:IN){
        if ( ! ( ( index.var.ii[vv] - 1 ) %in% index.ness.ii  ) ){
            l1 <- c( l1 , index.var.ii[vv] )                              
                        }
                    }
	if ( length(l1) > 0 ){
		vis1 <- vis1[ - l1 ]
					}
    return(vis1)
            }
#*****************************************************			
# equality of two vectors
equal_vecs <- function(v1,v2){
    val <- FALSE
    if ( length(v1) == length(v2) ){
        if ( sum( v1 != v2 ) == 0){
            val <- TRUE
                        }
                }
    return(val)
        }
#***************************************************