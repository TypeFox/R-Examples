
#################################################################################
# calculate model implied probabilities in GDINA models
gdina.probitem <- function( Mj, Aj , delta , rule , linkfct , delta.summary ){
    I <- length(delta)
    pjj <- as.list( 1:I )
    ljjj <- rep(0,I)    
    for (ii in 1:I){
        # ii <- 1
        pjjt <- ( Mj[[ii]][[1]] %*% delta[[ii]] )[,1]
        names(pjjt) <- paste0("A",apply( Aj[[ii]]   , 1 , FUN = function(ll){ paste(ll , collapse="") } ) )
        if (linkfct == "logit"){ pjjt <- stats::plogis( pjjt ) }
        if (linkfct == "log"){ pjjt <- exp( pjjt ) }    
        pjj[[ii]] <- pjjt
        ljjj[ii] <- length(pjjt)
            }
    
    pjj <- unlist( pjj )
    res <- data.frame( "itemno" = rep(1:I , ljjj) , "skillcomb"= names(pjj) , "prob" = pjj )
    dres <- NULL
    for (ii in 1:I){
        dii <- delta.summary[ delta.summary$itemno == ii , ]
        dii <- dii[ nrow(dii) , c("item" , "rule" , "partype.attr" ) ]
        dres <- rbind( dres , dii )
                } 
    res <- cbind( dres[ res$itemno , ] , res )
    rownames(res) <- NULL
	return(res)
        }
#################################################################################