.packageName <- "cgh"

# Smith-Waterman algorithm method
# of locating copy number changes
# in microarray CGH data
#
# R 1.9.1 script
#
# (c) T.S.Price
# 2004-2008


###############################################################################################


# Smith-Waterman algorithm
sw <-
function( x, max.nIslands = NULL, trace = FALSE )
{
    .Call( "sw", x, max.nIslands, trace, NAOK = FALSE )
}

# permutation test for island scores
sw.perm.test <-
function( x, max.nIslands = 1, nIter = 1000, seed = NULL, trace = FALSE )
{
    .Call( "sw_permTest", x, max.nIslands, nIter, seed, trace, parent.frame(),
        PACKAGE = "cgh", NAOK = FALSE )
}

# robustness calculation
sw.rob <-
function( x, lo.func = function( x ) median( x ),
    hi.func = function( x ) median( x ) + .4 * mad( x ), prec = 100 )
{
    len <- length( x )
    z.lo <- lo.func( x )
    z.hi <- hi.func( x )
    z.seq <- seq( z.lo, z.hi, length = prec )
    r <- double( len )
    for ( z in z.seq ) {
        swz <- sw( x - z )
        top.seg <- swz$start[ 1 ]:( swz$start[ 1 ] + swz$length[ 1 ] - 1 )
        r[ top.seg ] <- r[ top.seg ] + 1 / prec
    }
    return( r )
}

# threshold calculation
sw.threshold <-
function( logratio, threshold.func = function( x ) median( x ) + .2 * mad( x ), sign = +1 )
{
    logratio <- sign * logratio
    threshold <- threshold.func( logratio )
    return( logratio - threshold )
}

# plot Smith-Waterman results
sw.plot <-
function( logratio, location = seq( length( logratio ) ),
    threshold.func = function( x ) median( x ) + .2 * mad( x ),
    sign = -1, highest = TRUE, expected = NULL, rob = NULL, legend = TRUE,
    xlab = "Chromosomal location", ylab = "Intensity log ratio", ... )
{
    my.line <-
    function( x, y, ... )
    {
        len <- length( x )
        run <- rle( y )[[ 1 ]]
        run.len <- length( run )
        j <- 1
        m <- 2 * x[ 1 ] - x[ 2 ]
        if ( run.len == 1 )
    	      lines(
                x = c( 3/2 * x[ 1 ] - 1/2 * x[ 2 ], 3/2 * x[ len ] - 1/2 * x[ len - 1 ] ),
                y = c( y[ 1 ], y[ len ] ), ... )
        else {
    		for ( i in 1:( run.len - 1 ) ) {
    		    k <- run[ i ]
    		    lines(
                    x = c( ( x[ j ] + m ) / 2, ( x[ j + k - 1 ] + x[ j + k ] ) / 2 ),
                    y = c( y[ j ], y[ j ] ), ... )
    		    lines(
                    x = rep( ( x[ j + k - 1 ] + x[ j + k ] ) / 2, 2 ),
                    y = c( y[ j ], y[ j + k ] ), ... )
                m <- x[ j + k - 1 ]
                j <- j + k
    		    }
    		lines(
    		    x= c( ( m + x[ j ] ) /  2, 3/2 * x[ len ] - 1/2 * x[ len - 1 ] ),
    		    y= c( y[ j ], y[ j ] ), ... )
        }
    }

    island.line <-
    function( x, y, start = 1, len = length( x ), edge = 0, ... )
    {
        if ( is.null( start ) || is.null( len ) ||
            length( start ) == 0 || length( len ) == 0 ||
            len <= 0 )
            return
        lenx <- length( x )
        x1 <- c( 2 * x[ 1 ] - x[ 2 ], x, 2 * x[ lenx ] - x[ lenx - 1 ] )
        x2 <- x1[ 1:( lenx + 1 ) ] + diff( x1 ) / 2
        lines(
            x = c( rep( x2[ start ], 2 ), rep( x2[ start + len ], 2 ) ),
            y = c( y - edge, y, y, y - edge ),
            ... )
    }

    log2 <- log( 2 )
    len <- length( logratio )

    par <- par()
    par( mar = c( 5, 4, 4, 5 ) + .1, yaxs = "r" )

    # calculate threshold
    threshold <- threshold.func( sign * logratio )

    plot(
        y = logratio,
        x = location,
        type = "n",
        xlab = "",
        ylab = "",
        ... )
    rug( location, ticksize= .01 )
    maxlr <- max( logratio )
    minlr <- min( logratio )
    axis(
        side = 4, at = seq( minlr, maxlr, length = 5 ),
        labels = c( "0", ".25"," .50", ".75", "1" ) )
    mtext( xlab, side = 1, line = 3, cex = .8 )
    mtext( ylab, side = 2, line = 3, cex = .8 )
    abline( h = maxlr, lty = 2 )
    abline( h = minlr, lty = 2 )

    if ( !is.null( rob ) ) {
        mtext( "Robustness", side = 4, line = 3, cex = .8 )
        my.line(
            y = ( maxlr - minlr ) * rob + minlr,
            x = location,
            col = "#99ffff" )
    }

    # highest scoring island
    if ( highest ) {
        x <- sign * logratio - threshold
        swx <- sw( x )
        if ( length( swx$score ) ) {
            island.line(
                x = location,
                y = maxlr + ( maxlr - minlr ) * .02,
                start = swx$start[ 1 ],
                len = swx$length[ 1 ],
                edge = ( maxlr - minlr ) * .01,
                col = "#0000ff" )
        }
    }

    # threshold
    my.line( y = rep( sign * threshold, len ), x = location, col = "#00ff00" )

    # expected values
    if ( !is.null( expected ) ) {
        my.line( y = log2( expected ) - 1, x = location, col = "#ff0000" )
    }

    # plot points
    points( y = logratio, x = location, pch = 4, cex = .5 )

    # legend
    if ( legend ) {
        sel <- c( TRUE, highest, !is.null( expected ), !is.null( rob ) )
        legend.str <- c( "threshold value", "highest-scoring island", "expected values", "robustness" )
        legend.col <- c( "#00ff00", "#0000ff", "#ff0000", "#99ffff" )
        legend( location[ 1 ], minlr + ( maxlr - minlr ) * .25, legend.str[ sel ],
            lty = rep( 1, sum( sel ) ), col = legend.col[ sel ], cex = 0.8 )
    }

    par( mar = par$mar, yaxs = par$yaxs )
}



###############################################################################################



.noGenerics <- TRUE

##.onUnload <- function( libpath ) library.dynam.unload( "cgh", libpath )
