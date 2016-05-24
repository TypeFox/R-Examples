# print method for preCKrigePolygons class object

"print.preCKrigePolygons" <- function(object)
{
    if (!inherits(object, "preCKrigePolygons"))
    {
	stop("Not a preCKrigePolygons object")
    }
    n <- length( object@covmat )
    model <- object@model
    class(model) <- "covmodel"	
    
    pix.nrows <- object@pixconfig[[1]]$nrows
    pix.ncols <- object@pixconfig[[1]]$ncols 
    pix.cw <- object@pixconfig[[1]]$colwidth
    pix.rw <- object@pixconfig[[1]]$rowwidth
    m.poly.a <- mean(unlist(lapply(object@polygons, function(x){x@area})))
    sd.m.poly.a <- sd(unlist(lapply(object@polygons, function(x){x@area})))/sqrt(n)
    m.pix.poly <- mean(unlist( lapply(object@pixconfig, function(x){mean(x$no.pix.in.poly)})))
    sd.m.pix.poly <- sd(unlist( lapply(object@pixconfig, function(x){mean(x$no.pix.in.poly)})))/sqrt(n)
    poly.as.point <- sum(unlist( lapply(object@pixconfig, function(x){sum(x$sa.polygons)})))
    m.neig.n <-   mean( unlist( lapply(object@pixconfig, function(x){ length(x$posindex) - 1 } ) ) )
    min.neig.n <-   min( unlist( lapply(object@pixconfig, function(x){ length(x$posindex) - 1 } ) ) )
    max.neig.n <-   max( unlist( lapply(object@pixconfig, function(x){ length(x$posindex) - 1 } ) ) )
    
    sd.m.neig.n <- sd( unlist( lapply(object@pixconfig, function(x){ length(x$posindex) - 1 } ) ) )/sqrt(n)
    dim.data <- dim(object@data)
    
    cat("\n", "preCKrigePolygons object of",n,"polygons")
    cat("\n", "Mean polygon area:", paste(round(m.poly.a,2), " (sd of mean: ", round(sd.m.poly.a,3), ")", sep=""),"\n")
    cat("\n", "Average number of neighbours per polygon:", paste(round(m.neig.n,2), " (sd of mean: ", round(sd.m.neig.n,3), ")", sep = "")) 
    if(m.neig.n > 0)
    {
	cat("\n", "Minimum number of neighbours:", min.neig.n,"maximum number of neighbours:", max.neig.n,"\n")
    }

    cat("\n", "The block-block variance-covariance matrices\n", "are approximated by a gird of pixels:\n")
    cat(" Pixel gird dimension:", pix.nrows, "rows,", pix.ncols,"columns \n")
    cat(" Pixel size: ", paste("width: ", pix.cw,", height: ", pix.rw, sep = ""), "\n") 
    cat(" Average pixel number per polygon:", paste(round( m.pix.poly,2), " (sd of mean: ", round(sd.m.pix.poly,3), ")", sep=""),"\n")
    if(dim.data[1] == 0 & dim.data[2] ==0)
    {
	cat("\n", "The polygons have non attributes.\n")
    }else{
	cat("\n","The object contains a", paste("(",dim.data[1], ", ", 
		dim.data[2] ,")-data", sep =""), "frame with polygon 
		attributes.\n")
    }

    invisible(object )
}
#
setMethod("show", c("preCKrigePolygons"), print.preCKrigePolygons)

# summary method for preCKrigePolygons class object
summary.preCKrigePolygons <- function(object,...)
{
    if (!inherits(object, "preCKrigePolygons"))
    {
	stop("Not a preCKrigePoints object")
    }
    n <- length( object@covmat )
    model <- object@model
    class( model) <- "covmodel"	
    
    pix.nrows <- object@pixconfig[[1]]$nrows
    pix.ncols <- object@pixconfig[[1]]$ncols 
    pix.cw <- object@pixconfig[[1]]$colwidth
    pix.rw <- object@pixconfig[[1]]$rowwidth
    m.poly.a <- mean(unlist(lapply(object@polygons, function(x){x@area})))
    sd.m.poly.a <- sd(unlist(lapply(object@polygons, function(x){x@area})))/sqrt(n)
    m.pix.poly <- mean(unlist( lapply(object@pixconfig, function(x){mean(x$no.pix.in.poly)})))
    sd.m.pix.poly <- sd(unlist( lapply(object@pixconfig, function(x){mean(x$no.pix.in.poly)})))/sqrt(n)
    poly.as.point <- sum(unlist( lapply(object@pixconfig, function(x){sum(x$sa.polygons)})))
    m.neig.n <-   mean( unlist( lapply(object@pixconfig, function(x){ length(x$posindex) - 1 } ) ) )
    min.neig.n <-   min( unlist( lapply(object@pixconfig, function(x){ length(x$posindex) - 1 } ) ) )
    max.neig.n <-   max( unlist( lapply(object@pixconfig, function(x){ length(x$posindex) - 1 } ) ) )
    
    sd.m.neig.n <- sd( unlist( lapply(object@pixconfig, function(x){ length(x$posindex) - 1 } ) ) )/sqrt(n)
    dim.data <- dim(object@data)
    
    cat("\n", "preCKrigePolygons object of",n,"polygons")
    cat("\n", "Mean polygon area:", paste(round(m.poly.a,2), " (sd of mean: ", round(sd.m.poly.a,3), ")", sep=""),"\n")
    cat("\n", "Average number of neighbours per polygon:", paste(round(m.neig.n,2), " (sd of mean: ", round(sd.m.neig.n,3), ")", sep = "")) 
    if(m.neig.n > 0)
    {
	cat("\n", "Minimum number of neighbours:", min.neig.n,"maximum number of neighbours:", max.neig.n,"\n")
    }

    cat("\n", "The block-block variance-covariance matrices\n", "are approximated by a gird of pixels:\n")
    cat(" Pixel gird dimension:", pix.nrows, "rows,", pix.ncols,"columns \n")
    cat(" Pixel size: ", paste("width: ", pix.cw,", height: ", pix.rw, sep = ""), "\n") 
    cat(" Average pixel number per polygon:", paste(round( m.pix.poly,2), " (sd of mean: ", round(sd.m.pix.poly,3), ")", sep=""),"\n")
	cat(" Polygons treated as point:", poly.as.point , "\n")
    cat("\n", "Parameter of the used Spatial-covariance model:\n\n")
    print( model, right = F)
    if(dim.data[1] == 0 & dim.data[2] ==0)
    {
	cat("\n", "The polygons have non attributes.\n")
    }else{
	cat("\n","The object contains a", paste("(",dim.data[1], ", ", dim.data[2] ,")-data", sep =""), "frame with the following polygon attributes:\n")
    }
print(str(object@data))
    invisible(object)
}

setMethod("summary", c("preCKrigePolygons"), summary.preCKrigePolygons)

# print method for preCKrigePoints class object
"print.preCKrigePoints" <- function(object)
{
    x <- object
    if (!inherits(x, "preCKrigePoints"))
    {
	stop("Not a preCKrigePolygons object")
    }
    n <- length( object@covmat )
    model <- object@model
    class(model) <- "covmodel"	
    
   
    m.neig.n <-   mean( unlist( lapply(object@posindex, function(x){ length(x) - 1 } ) ) )
    min.neig.n <-   min( unlist( lapply(object@posindex, function(x){ length(x) - 1 } ) ) )
    max.neig.n <-   max( unlist( lapply(object@posindex, function(x){ length(x) - 1 } ) ) )
    sd.m.neig.n <- sd( unlist( lapply(object@posindex, function(x){ length(x) - 1 } ) ) ) / sqrt(n)
    dim.data <- dim(object@data)
    
    cat("\n", "preCKrigePoints object of",n,"points")
    cat("\n", "Average number of neighbours per point:", paste(round(m.neig.n,2), " (sd of mean: ", round(sd.m.neig.n,3), ")", sep = "")) 
    if(m.neig.n > 0)
    {
	cat("\n", "Minimum number of neighbours:", min.neig.n,"maximum number of neighbours:", max.neig.n,"\n")
    }

 
    if(dim.data[1] == 0 & dim.data[2] ==0)
    {
	cat("\n", "The points have non attributes.\n")
    }else{
	cat("\n","The object contains a", paste("(",dim.data[1], ", ", 
		dim.data[2] ,")-data", sep =""), "frame with point attributes.\n")
    }

    invisible(object)
}
#
setMethod("show", c("preCKrigePoints"), print.preCKrigePoints)



# summary method for preCKrigePoints class object
summary.preCKrigePoints <- function(object,...)
{
    if (!inherits(object, "preCKrigePoints"))
    {
	stop("Not a preCKrigePoints object")
    }
    n <- length( object@covmat )
    model <- object@model
    class(model) <- "covmodel"	
    
   
    m.neig.n <-   mean( unlist( lapply(object@posindex, function(x){ length(x) - 1 } ) ) )
    min.neig.n <-   min( unlist( lapply(object@posindex, function(x){ length(x) - 1 } ) ) )
    max.neig.n <-   max( unlist( lapply(object@posindex, function(x){ length(x) - 1 } ) ) )
    sd.m.neig.n <- sd( unlist( lapply(object@posindex, function(x){ length(x) - 1 } ) ) ) / sqrt(n)
    dim.data <- dim(object@data)
    
    cat("\n", "preCKrigePoints object of",n,"points")
    cat("\n", "Average number of neighbours per point:", paste(round(m.neig.n,2), " (sd of mean: ", round(sd.m.neig.n,3), ")", sep = "")) 
    if(m.neig.n > 0)
    {
	cat("\n", "Minimum number of neighbours:", min.neig.n,"maximum number of neighbours:", max.neig.n,"\n")
    }


    cat("\n", "Parameter of the used Spatial-covariance model:\n\n")
    print( model, right = F)
    if(dim.data[1] == 0 & dim.data[2] ==0)
    {
	cat("\n", "The polygons have non attributes.\n")
    }else{
	cat("\n","The object contains a", paste("(",dim.data[1], ", ", dim.data[2] ,")-data", sep =""), "frame with the following point attributes:\n")
    }
print(str(object@data))
    invisible(object)
}
setMethod("summary", c("preCKrigePoints"), summary.preCKrigePoints)

# print method for CKrige.exout.polygons
"print.CKrige.exout.polygons" <- function(x,...)
{
    if (!inherits(x, "CKrige.exout.polygons"))
    {
	stop("Not a CKrige.exout.polygons object")
    }

    
    cat("\n", "List object, extended output of the CKrige function.\n\n")

    cat(" Components:\n") 
    cat(" $object: SpatialPolygonsDataFrame with kriging predictions (standard output object of CKrige)\n")
    cat(" $krig.method: used kriging method\n")
    if( length(x)  == 7)
    {
	cat(" $CMCK.par: list of CMCK parameter matrices, P1, Q1 and K\n")
    }
    cat(" $parameter: list of gls coefficients\n")
    cat(" $sk.weights: matrix of the simple kriging weights\n")
    cat(" $inv.Sigma: inverse covariance matrix of the data\n")
    cat(" $residuals: gls residuals\n")
    cat("\n use summary to get the summary of the kriging predictions\n")
   
    invisible(x)
}

## print method for CKrige.exout.polygons
"summary.CKrige.exout.polygons" <- function(object,...)
{
    x <- object
    if (!inherits(x, "CKrige.exout.polygons"))
    {
	stop("Not a CKrige.exout.polygons object")
    }
    
   if(x$krig.method == 1)
    {
	cat("\nSummary of the data frame with the universal kringing results.\n\n")
    }
    
    if(x$krig.method == 2)
    {
	cat("\nSummary of the data frame with the constrained kringing results.\n\n")
    }
    
    if(x$krig.method == 3)
    {
	cat("\nSummary of the data frame with the covariance-matching constrained kringing results.\n\n")
    }
    
    summary( x$object@data )
}
#

# print method for CKrige.exout.points
"print.CKrige.exout.points" <- function(x,...)
{
    if (!inherits(x, "CKrige.exout.points"))
    {
	stop("Not a CKrige.exout.points object")
    }

    
    cat("\n", "List object, extended output of the CKrige function.\n\n")

    cat(" Components:\n") 
    cat(" $object: SpatialPointsDataFrame with kriging predictions (standard output object of CKrige)\n")
    cat(" $krig.method: used kriging method\n")
    if( length(x)  == 7)
    {
	cat(" $CMCK.par: list of CMCK parameter matrices, P1, Q1 and K\n")
    }
    cat(" $parameter: list of gls coefficients\n")
    cat(" $sk.weights: matrix of the simple kriging weights\n")
    cat(" $inv.Sigma: inverse covariance matrix of the data\n")
    cat(" $residuals: gls residuals\n")
    cat("\n use summary to get the summary of the kriging predictions\n")
   
    invisible(x)
}
#
# summary method for CKrige.exout.points
"summary.CKrige.exout.points" <- function(object,...)
{
    x <- object
    
    if (!inherits(x, "CKrige.exout.points"))
    {
	stop("Not a CKrige.exout.points object")
    }

    if(x$krig.method == 1)
    {
	cat("\nSummary of the data frame with the universal kringing results.\n\n")
    }
    
    if(x$krig.method == 2)
    {
	cat("\nSummary of the data frame with the constrained kringing results.\n\n")
    }
    
    if(x$krig.method == 3)
    {
	cat("\nSummary of the data frame with the covariance-matching constrained kringing results.\n\n")
    }
    
}
#









