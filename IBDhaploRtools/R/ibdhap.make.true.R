ibdhap.make.true <- function( true.filename){


    ## ## OLD FORMAT FROM fgl2ibd_4gam --> matrix rows=markers, cols=sets
    ## read in the whole data, remove names and transpose
    ##rawdat <- read.table( true.filename, sep = " ", header = FALSE )
    ##outdat <- t( rawdat[,-(1:2)] )
    ## remove the row of NA at the bottom
    ##dm <- dim( outdat )[1]
    ##outdat <- outdat[ -dm , ]

    ## read in columns [1] set [2] marker [3] cM position [4] an ordering [5] ibdlabel [6...] ibd vector
    
    rawdat <- read.table( true.filename, sep = " ", header = FALSE )
    sets <- unique( rawdat[,1] )
    num.markers <- sum( rawdat[,1]== sets[1] )

    outdat <- matrix( rawdat[,5], nrow = num.markers, ncol=length(sets), byrow=FALSE)
    colnames( outdat ) <- sets

    ## convert label to jacquard

    outdat <- apply( outdat, c(1,2), function(x){
        c(1,3,NA,4,2,5,6,9,11,10,7,13,12,14,8,15)[ x+1 ]
    })
    ## outdat is matrix with column per set, row per marker
    return( outdat )
}
