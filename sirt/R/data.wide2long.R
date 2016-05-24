


#####################################################################
# converts a data frame in wide format into long format
data.wide2long <- function( dat , id=NULL , X=NULL , Q=NULL){
    items <- base::setdiff( colnames(dat) , id )
    I <- length(items)
    N <- nrow(dat)
#    dat1 <- matrix( as.matrix(dat[,items]) , nrow= N*I , ncol=1 , byrow=TRUE )
    dat1 <- matrix( t(dat[,items]) , nrow= N*I , ncol=1 , byrow=FALSE )
    dat2 <- data.frame( dat1 )
    colnames(dat2) <- "resp"
    dat1 <- data.frame( "id_index" = rep( 1:N , each=I ) )
		
    if ( ! is.null(id) ){
        dat1 <- cbind( dat1 , rep( dat[, id] , each=I ) ) 
        colnames(dat1)[2] <- id    
				}  
    dat1 <- cbind( dat1 ,  "item" = rep( items , N ) ,
                "item_index" = rep(1:I , N ) , 
                "resp"= dat2$resp )
    if ( ! is.null(X) ){ 
        dat1 <- cbind( dat1 , X[ rep(1:N , each=I ) , ] )
                }
    rownames(dat1) <- NULL
	if ( ! is.null(Q) ){
	  if ( is.null(colnames(Q) ) ){
	     colnames(Q) <- paste0("q",1:ncol(Q) )
						}
	   if ( sum( colnames(Q) %in% "item" ) == 0 ){
					Q <- as.data.frame(Q)
					Q$item <- colnames(dat)
							}
		dat1 <- base::merge( x = dat1 , y = Q , by ="item" , all.x=TRUE )
					}
	dat1 <- dat1[ order( 10000*dat1$id_index + dat1$item_index ) , ]		
	dat1 <- data.frame( "rowindex" = 1:(N*I) , dat1 )
    return(dat1)
        }
#####################################################################
