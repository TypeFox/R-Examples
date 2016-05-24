f.points.CKrige<- function(
    		formula,
    		data,
		locations,
		object,
		method = 2,
		ex.out = F
	    )
### purpose: calculate the constrained (C), covariance-matching constrained (CM)
###          and universal (U) kriging predictor for one block configuration
### ...........................................................................
### arguments:
###  pixconfig = list,
###                       1. element = coords of the lower left pixel 
###                       2. element = vector with first number = number
###                          of pixel second number = number of polygons
###                       3. element = logical vector, TURE means the
###                          appropriate polyon area is smaller than the
###                          pixel area
###                       4. element = vector, numbers are the pixel centers
###                          in the appropiate polygon  e.g 0, 2, 1 means 0
###                          pixels in first poly  2 pixel in second poly
###                          and 1 pixel in
###                       5. element = matrix, coords of the centrois of the
###                          polygons
###                       6. element = list, each element is a number, the number
###                          describe the polygon in which the appropriate
###                          pixel center is, e.g 13th element hast number
###                          3 that means that the 13th pixel center is
###                            in the 3th polygon 
###                       7. element = vector, numbers represent which polyon are
###                          considered e.g. c(6,3,10) means its a polygon
###                          configuration with the 6th, 3th and 10th polygon
### t.bb.cov = list, elements are block-block covariances matrices
### t.pred.designmat = n x number of predicton location matrix, model matrix of
###                 the prediction location
### model  = list, every element describ the autocorrelation
### t.cov.spline = Spline object of the  covariance between the sampling points
###                t.sample and the cells (rectangles) 
### t.inv.Sigma = n x n inverse of the covariance matrix of the sampling points
###              t.sample
### t.list.beta.coef = list of two elements
###                        1. element = 1 x pi vector with the parameter of the
###                           general least square regression
###                        2. element = pi x pi matrix of the covariance between
###                           the general least square parameter
### t.Y =  n x 1 vector target value at the sampling points
### t.X = design matrix of the samples
### t.sample = n X 2 matrix with the coordinates of the smapling points,
###            first row are the x values and second row are the y values
### t.eps =  is the tuning parameter. It will replace those eigen values
###          when the positive definiteness of the Q or P matrix  is not
###          satisfied.
###
### ..............................................................................
{
#
t.support.designmat <- model.matrix(object = formula, data = data)
#
t.response  <- model.response( model.frame( formula, data ) )
#
locations <- terms( x= locations)
attr(locations, "intercept") = 0
locations <- model.matrix(object = locations, data = data)
#
if( dim( object@data )[1] == 0 && dim( object@data )[2] == 0 ) # ordinary case
{
    t.pred.designmat <- matrix( rep( 1, dim( object@coords )[1] ) , ncol = 1 )
}
else
{
    # t.terms.wr = terms without the response variable
    t.terms.wr <- delete.response( termobj = terms( x = formula ) )
    t.pred.designmat <- model.matrix( object = t.terms.wr, data = object@data)
    rm( t.terms.wr )
}
#
### covmodel of measurement free error process
model.me.free <- object@model[unlist(lapply(1:length(object@model), function(i,m){m[[i]]$model != "mev"}, m = object@model))]
## Kovarianzmatrix der ST?tzpunkte
t.support.covmat <- f.covmat.support( object@model, locations )   

# Inversere Choleskymatrix L^-1
t.ichol <- forwardsolve(
    t.chol  <- t( chol(t.support.covmat ) ),
    diag( nrow( t.support.covmat ) )
    )
    
# Inverse Kovarianzmatrix
t.isupport.covmat <- crossprod( t.ichol, t.ichol )  
#orthogonalisierte Designmatrix und Messwerte der Datenpunkte
t.orth.designmat <- t.ichol %*% t.support.designmat
t.orth.data <- t.ichol %*% t.response
    
# QR  <- Zerlegung der orthogonalisierten Designmatrix
t.qr <- qr( as.data.frame( t.orth.designmat  )  )

t.R <- qr.R( t.qr )
t.iR <- backsolve( t.R, diag( nrow( t.R ) ) )
# General least square Koeffizienten (beta GLS)
t.beta.coef <- matrix( qr.coef( t.qr, as.data.frame( t.orth.data )), ncol = 1)

# (Co)-Varianzen der beta Koeffizienten
t.cov.beta.coef  <- crossprod( t( t.iR ), t( t.iR ) ) 
# GLS-Residuen
t.orth.resid <- qr.resid( t.qr, t.orth.data )
t.residuals  <-  t.chol   %*% t.orth.resid 

# lapply Laufvariable
t.index  <- as.list( 1:dim( object@coords )[1] )

t.krige.res <- lapply(t.index,
    function(i, coords, posindex, t.bb.covmat.list, t.pred.designmat,
	locations, model, t.ichol, t.isupport.covmat, t.beta.coef,
	t.cov.beta.coef , t.orth.resid, t.orth.designmat, t.iR, method)
    {

	t.indices <- posindex[[ i ]]

	coords <- matrix( coords[ t.indices, ], ncol= 2 )

	## m X m covariance matrix of the block 
	t.bb.covmat <- as.matrix( t.bb.covmat.list[[ t.indices[1] ]] ) 

	## block design matrix
	t.pred.designmat <- matrix(t.pred.designmat[ t.indices, ],
	    ncol = length( t.beta.coef ), nrow = length( t.indices ) )

	#cat("Parzelle: ", i, "\n")
	# Kovarianz zwischen den St?tzpunkten und den Zielbl?cken (n X m Matrix)
	t.dist <- f.row.dist( locations, coords )

        t.pred.covmat <- matrix(f.pp.cov( as.vector(t.dist), model), ncol = ncol(t.dist) )

	rm( t.dist )

	t.pred.covmat.ichol.trans <-  crossprod( t.pred.covmat, t( t.ichol) ) 
	
	# t(c) %*% Sigma^-1 = t(c) %*% t(L^-1) %*% L^1
	t.sk.weights  <- t.pred.covmat.ichol.trans %*% t.ichol
	

	## least square prediction

	t.lin.trend.est <-  crossprod( t( t.pred.designmat), t.beta.coef )


	## residuals weithed by the simple kriging weights

	t.weighted.resid <- crossprod( t( t.pred.covmat.ichol.trans ), t.orth.resid )

	## UK Vorhersage
	t.uk.pred <- t.lin.trend.est  + t.weighted.resid

	# UK MSPE
	t.aux <- crossprod( t( t.pred.designmat - crossprod(t(t.pred.covmat.ichol.trans), t.orth.designmat) ) , t.iR)
	t.uk.mspe <- t.bb.covmat - tcrossprod( t.pred.covmat.ichol.trans, t.pred.covmat.ichol.trans ) + t.aux %*% t( t.aux )
	
	#print(t.bb.covmat - tcrossprod( t.pred.covmat.ichol.trans, t.pred.covmat.ichol.trans ))

	# Kriging Vorhersage und MSPE
	t.result <- f.kriging.method( t.lin.trend.est, t.cov.beta.coef, t.weighted.resid, 
    					t.uk.mspe, t.bb.covmat, t.pred.designmat,
					t.orth.designmat, t.pred.covmat.ichol.trans, method)

	if(method == 1)
	{
	    return( list( prediction = t.result[[1]], sqrt.mspe = t.result[[2]], sk.weights = t.sk.weights ) )
	}
	if( method == 2 )
	{
	    return( list( prediction = t.result[[1]], sqrt.mspe = t.result[[2]],
		    P1 = t.result[[3]], Q1 =t.result[[4]], K = t.result[[5]], sk.weights = t.sk.weights) )
	}
	if( method == 3)
	{
	    return( list( prediction = t.result[[1]], mspe = t.result[[2]],
		    P1 = t.result[[3]], Q1 =t.result[[4]], K = t.result[[5]], sk.weights = t.sk.weights ) )
	}

    },
    coords = object@coords,
    posindex = object@posindex,
    t.bb.covmat.list = object@covmat,
    t.pred.designmat = t.pred.designmat,
    locations = locations,
    model = model.me.free,
    t.ichol = t.ichol,
    t.isupport.covmat = t.isupport.covmat,
    t.beta.coef = t.beta.coef,
    t.cov.beta.coef = t.cov.beta.coef,
    t.orth.resid = t.orth.resid,
    t.orth.designmat = t.orth.designmat,
    t.iR = t.iR,
    method = method
) ## end of lapply
##
## output
##
# # # 
if( method == 1 ){
   krige.result <- data.frame( matrix( NaN, nrow = nrow( t.pred.designmat ), ncol = 2 ) )
    krige.result[,1] <- unlist( lapply( t.krige.res, function( poly ){ return( poly$prediction ) } ) )
    krige.result[,2] <- unlist( lapply( t.krige.res, function( poly ){ return( poly$sqrt.mspe ) } ) )
    colnames( krige.result) = c("prediction", "prediction.se")
}
#
if( method == 2 )
{
    krige.result <- data.frame( matrix( NaN, nrow = nrow(t.pred.designmat), ncol = 5 ) )
    krige.result[,1] <- unlist( lapply( t.krige.res, function( poly ){ return( poly$prediction ) } ) )
    krige.result[,2] <- unlist( lapply( t.krige.res, function( poly ){ return(  poly$sqrt.mspe  ) } ) )
    krige.result[,3] <- unlist( lapply( t.krige.res, function( poly ){ return( poly$P1 ) } ) )
    krige.result[,4] <- unlist( lapply( t.krige.res, function( poly ){ return( poly$Q1 ) } ) )
    krige.result[,5] <- unlist( lapply( t.krige.res, function( poly ){ return( poly$K ) } ) )
    colnames( krige.result) = c("prediction", "prediction.se", "P1", "Q1", "K")
}
#
if( method == 3)
{
    krige.result <- data.frame( matrix( NaN, nrow = nrow(t.pred.designmat), ncol = 5 ) )
    krige.result[,1] <- unlist( lapply( t.krige.res, function( poly ){ return( poly$prediction[1,1] ) } ) )
    krige.result[,2] <- unlist( lapply( t.krige.res, function( poly ){ return( sqrt( poly$mspe[1,1] ) ) } ) )
    krige.result[,3] <- unlist( lapply( t.krige.res, function( poly ){ return( poly$P1[1,1] ) } ) )
    krige.result[,4] <- unlist( lapply( t.krige.res, function( poly ){ return( poly$Q1[1,1] ) } ) )
    krige.result[,5] <- unlist( lapply( t.krige.res, function( poly ){ return( poly$K[1,1] ) } ) )
    colnames( krige.result) = c("prediction", "prediction.se","P1.11", "Q1.11", "K.11") 
}
#
if(  ex.out == T)
{
    if( method == 1 | method == 2 )
    {
	#sk.weights.matrix
	if( is.null( dim(t.krige.res[[1]]$sk.weights) ) )
	{
	    sk.weights.matrix <- t(matrix( unlist( lapply( t.krige.res, function( points ){ return( points$sk.weights ) } ) ), ncol = nrow( data ), byrow = T ))
	}else{
	    sk.weights.matrix <- t(matrix( unlist( lapply( t.krige.res, function( points ){ return( points$sk.weights[1,] ) } ) ), ncol = nrow( data ), byrow = T))
	    }
	## extended output is a list
	object <- SpatialPointsDataFrame( SpatialPoints( object@coords ), data = krige.result, match.ID = F)
	res <- list(
	    object = object,
	    krig.method = method,
	    parameter = list( beta.coef = t.beta.coef, cov.beta.coef = t.cov.beta.coef),
	    sk.weights = sk.weights.matrix,
	    inv.Sigma = t.isupport.covmat,
	    residuals = t.residuals
	)
	class( res )  <- "CKrige.exout.points"
	return( res )
    }
    else
    {
	object <- SpatialPointsDataFrame( SpatialPoints(object@coords ), data = krige.result, match.ID = F)
	P1.list = lapply( t.krige.res, function( points ){ return(points$P1)})
	Q1.list = lapply( t.krige.res, function( points ){ return(points$Q1)})
	K.list = lapply( t.krige.res, function( points ){ return(points$K)})
	sk.weights.list <-lapply( t.krige.res, function( points ){ return(t(points$sk.weights))}) 
      res <- list(
	    object = object,
	    krig.method = method,
	    CMCK.par = list(P1 = P1.list, Q1 =Q1.list,K = K.list),
	    parameter = list( beta.coef = t.beta.coef, cov.beta.coef = t.cov.beta.coef),
	    sk.weights = sk.weights.list,
	    inv.Sigma = t.isupport.covmat,
	    residuals = t.residuals
	)
	class( res )  <- "CKrige.exout.points"
	return( res ) 
    }
}
else
{
    
    row.names( krige.result ) <- row.names( object@coords )
    object <- SpatialPointsDataFrame( SpatialPoints( object@coords ), krige.result )
    return( object )  
}
#    
} # end of function

