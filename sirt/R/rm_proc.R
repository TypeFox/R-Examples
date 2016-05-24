
##########################################
# Data preprocessing rater models
.prep.data.rm <- function( dat , pid , rater ){
    rater <- paste( rater )
    # create table of rater indizes
    rater.index <- data.frame( "rater" = sort( unique( rater )) )
    rater.index$rater.id <- seq( 1 , nrow(rater.index) )
    RR <- nrow(rater.index)
    # create table of person indizes
    person.index <- data.frame( "pid" = sort( unique( pid )) )
    person.index$person.id <- seq( 1 , nrow(person.index) )
    PP <- nrow( person.index )
    # number of variables
    VV <- ncol(dat)
    vars <- colnames(dat)    
    # create data frame with crossed items and raters
    dat2 <- data.frame( matrix( NA , nrow=PP , ncol=RR*VV ) )
    colnames(dat2) <- paste0( rep(vars , RR ) , "-" , rep( rater.index$rater , each=VV) )
    rownames(dat2) <- person.index$pid
    for (rr in 1:RR){
        # rr <- 1
        ind.rr <- which( rater == rater.index$rater[rr] )
        dat.rr <- dat[ ind.rr  , ]
        pid.rr <- pid[ ind.rr ]
        i1 <- match(  pid.rr , person.index$pid )
        colnames(dat.rr) <- NULL
        dat2[ i1 , VV*(rr-1) + 1:VV ] <- dat.rr
            }
	# variable list
    dataproc.vars <- list( "item.index" = rep( 1:VV , RR )	 ,
			"rater.index" = rep(1:RR , each=VV ) )
	# arrange response data
    dat2.resp <- 1 - is.na(dat2)
	dat20 <- dat2
	dat2[ dat2.resp == 0 ] <- 0
    res <- list( "dat2" = dat2 , "dat2.resp" = dat2.resp , 
				"dat2.NA" = dat20 , 
				"dat" = dat , "person.index" = person.index , 
                "rater.index" = rater.index , "VV"=VV , "N" = PP , "RR" = RR ,
				"dataproc.vars" = dataproc.vars )
    return(res)
        }
#################################################