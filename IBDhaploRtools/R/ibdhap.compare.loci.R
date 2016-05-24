ibdhap.compare.loci <- function( calls, true, data.type ){
### Compare called states against true states, per locus


    ## Check input data type 
    if(length(data.type)>1){ stop("data.type improperly indicated")}
    else if( is.element("h", data.type)){ no.ibd.ind = 15}
    else if( is.element("g", data.type)){no.ibd.ind = 9}
    else if( is.element("r", data.type)){no.ibd.ind = 4}
    else{ stop("data.type improperly indicated")}
    
    ## Classify all loci
    ibdsites <- which( true < no.ibd.ind, arr.ind = T )
    noibdsites <- which( true == no.ibd.ind, arr.ind = T )

    ## ALL SITES SUMMARIES
    all <- matrix(NA, nrow = 5, ncol = 1 )
    rownames( all ) <- c("Number of sites",
                         "Called Correctly",
                         "Called IBD Incorrectly",
                         "Called noIBD Incorrectly",
                         "No-call" )
    colnames( all ) <- " "
    all[1,] <- length( calls )
    all[2,] <- mean( true == calls )
    all[3,] <- mean( true != calls & calls > 0 & calls < no.ibd.ind )
    all[4,] <- mean( true != calls & calls > 0 & calls == no.ibd.ind )
    all[5,] <- mean( true != calls & calls == 0 ) 
    
    ## IBD SITE SUMMARIES
    ibd <- matrix(NA, nrow = 5, ncol = 1 )
    rownames( ibd ) <- c("Number of sites",
                         "Called Correctly",
                         "Called as wrong IBD",
                         "Called as no IBD",
                         "No-call")
    colnames( ibd ) = " " 
    ibd[1,] <- dim( ibdsites )[1]
    ibd[2,] <- mean( true[ ibdsites ] == calls[ ibdsites ] )*100
    ibd[3,] <- mean( true[ ibdsites ] != calls[ ibdsites ] &
                    calls[ ibdsites ] < no.ibd.ind  &  calls[ ibdsites ] > 0)*100
    ibd[4,] <- mean( true[ ibdsites ] != calls[ ibdsites ] & calls[ ibdsites ] == no.ibd.ind  )*100
    ibd[5,] <- mean( true[ ibdsites ] != calls[ ibdsites ] & calls[ ibdsites ] == 0  )*100

    ## NON-IBD SITE SUMMARIES
    nonibd <- matrix(NA, nrow = 4, ncol = 1)
    rownames( nonibd ) <- c("Number of sites",
                            "Called Correctly",
                            "Called as IBD",
                            "No-call")
    colnames( nonibd ) = " "
    nonibd[1,] <- dim( noibdsites )[1]
    nonibd[2,] <- mean( true[ noibdsites ] == calls[ noibdsites ] )*100
    nonibd[3,] <- mean( true[ noibdsites ] != calls[ noibdsites ] & calls[ noibdsites ] < no.ibd.ind & calls[ noibdsites ] >0 )*100
    nonibd[4,] <- mean( true[ noibdsites ] != calls[ noibdsites ] & calls[ noibdsites ] == 0 )*100

    ## CATEGORIES OF ALL SITES
    infcat <- hist(as.matrix(calls ), breaks = seq(-0.5,no.ibd.ind+0.5,1), plot = FALSE)$density
    trucat <-  hist(as.matrix(true ), breaks = seq(-0.5,no.ibd.ind+0.5,1), plot = FALSE)$density
    categories <- data.frame( "Category" = 0:no.ibd.ind, "Inferred" = infcat*100,  "True" =trucat*100 ) 
    rownames( categories ) = 0:no.ibd.ind
    
     return(list( all = all, ibd=ibd, nonibd=nonibd, categories=categories ) )
}
