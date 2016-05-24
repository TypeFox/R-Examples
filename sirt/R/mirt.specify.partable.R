

####################################################################### 
##   function parameter table specification
mirt.specify.partable <- function(mirt.partable, parlist, verbose=TRUE) {
    mirt.partable0 <- mirt.partable
    mirt.partable0$prior.type <- paste(mirt.partable0$prior.type)
    if(length(unique(mirt.partable$group)) == 1L) ####
        parlist <- list('all' = parlist) ####
    parlistfull <- parlist ####
    gnames <- names(parlistfull) ####
    if(!all(gnames %in% unique(mirt.partable$group))) ####
        stop('Defined group variable names in parlist do not match mirt.partable') ####
    #@@@
    items <- paste0( mirt.partable$item )        
    for(gg in gnames){ ####
        parlist <- parlistfull[[gg]] ####
        LP <- length(parlist) ####  moved this down from first line
        for (pp in 1:LP) {
            # pp <- 1 extract parameter
            par.pp <- names(parlist)[[pp]]
            # extract data frame with parameter specifications
            dfr.pp <- parlist[[pp]]
            #@@ check item names
            g1 <- sort(setdiff( rownames(dfr.pp)  , items ) )
            dfr.pp <- dfr.pp[ intersect( rownames(dfr.pp) , items ) , , drop=FALSE]      
            #@@@
            if (verbose){
                cat("*** Process group " , gg , " - parameter " , par.pp , "\n")                                                
                if ( length(g1) > 0 ){
                    cat("  - Following items do not exist in parameter table: " , paste0( g1 , collapse=" " ) , "\n")                
                                }
                           }            
                           
            # loop through items and parameters
            NI.pp <- nrow(dfr.pp)
            NP.pp <- ncol(dfr.pp)
            if (NI.pp>=1){
            for (ii in 1:NI.pp) {
                for (cc in 1:NP.pp) {
                    # ii <- 1 ; cc <- 2
                    ind <- which((paste(mirt.partable0$item) == rownames(dfr.pp)[ii]) &
                                     (paste(mirt.partable0$name) == colnames(dfr.pp)[cc]) &
                                     mirt.partable0$group == gg) ####
                    if (!is.na(dfr.pp[ii, cc])) {					
						if(!length(ind))
								warning(sprintf("Parameter \'%s\' is not relevant for item %s",
														 colnames(dfr.pp)[cc], rownames(dfr.pp)[ii]))					
                      if (par.pp != "prior.type") {
                        mirt.partable0[ind, par.pp] <- dfr.pp[ii, cc]
                      } else {
                        mirt.partable0[ind, par.pp] <- paste(dfr.pp[ii, cc])
                      }
                    }
                }
            }
            }
        }
    }
    if(!all(mirt.partable0$prior.type %in% c('none', 'norm', 'beta', 'lnorm'))) ####
        stop('Improper prior.type declared. Please only use the following:
             \'none\', \'norm\', \'beta\', \'lnorm\' ') ####
    return(mirt.partable0)
}
############################################################################
