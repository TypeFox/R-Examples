f.polygons.CKrige <- function(
    		formula,
    		data,
		locations,
		object,
		method = 2,
		ex.out = F
	    )
### purpose: calculate the constrained (C), covariance-matching constrained (CM)
###          and universal (U) kriging prediction
### fromula = linear trend as formula
### data = data.frame 
### locations = x,y coordinates of the sampling sites
### object = object of class preCKrigePoints or	preCKrigePolygons
### method = numeric 1 = UK, 2= CK, 3 = CMCK
### ex.out = extended output
### ..............................................................................
{
t.support.designmat <- model.matrix( formula, data)
#
t.response  <- model.response( model.frame( formula, data ) )
# 
locations <- terms( x= locations)
attr(locations, "intercept") = 0
locations <- model.matrix(object = locations, data = data)
#
if( dim( object@data )[1] == 0 && dim( object@data )[2] == 0 ) # ordinary case
{
    t.pred.designmat <- matrix( rep( 1, length( object@polygons ) ), ncol = 1 )
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
#
# Inversere Choleskymatrix L^-1
t.ichol <- forwardsolve(
    t.chol  <- t( chol( t.support.covmat ) ),
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
#
# (Co)-Varianzen der beta Koeffizienten
t.cov.beta.coef  <- crossprod( t( t.iR ), t( t.iR ) ) 
# GLS-Residuen
t.orth.resid <- qr.resid( t.qr, t.orth.data )
t.residuals  <-  t.chol   %*% t.orth.resid 

# lapply Laufvariable
t.index  <- as.list( 1:length( object@pixconfig ) )

t.krige.res <- lapply(t.index,
    function(i, pixconfig, t.bb.covmat.list, t.pred.designmat,
	locations, model, t.ichol, t.isupport.covmat, t.beta.coef,
	t.cov.beta.coef , t.orth.resid, t.orth.designmat, t.iR, method)
    {
	pixconfig <- pixconfig[[ i ]]

	t.indices <- pixconfig$posindex

	## m X m covariance matrix of the block 
	t.bb.covmat <- as.matrix( t.bb.covmat.list[[ t.indices[1] ]] ) 

	## block design matrix
	t.pred.designmat <- matrix(t.pred.designmat[ t.indices, ], ncol = nrow( t.beta.coef ))

	#cat("Parzelle: ", i, "\n")
	# Kovarianz zwischen den St?tzpunkten und den Zielbl?cken (n X m Matrix)
	t.pred.covmat <- f.point.block.cov( pixconfig = pixconfig, locations = locations, model = model)

       #  C'(L^-1)'

	#t.pred.covmat.ichol.trans <- t( crossprod(t.pred.covmat, t( t.ichol) ) )
	t.pred.covmat.ichol.trans <-  crossprod( t.pred.covmat, t( t.ichol) ) 

	# t(c) %*% Sigma^-1 = t(c) %*% t(L^-1) %*% L^1
	
	t.sk.weights  <- t.pred.covmat.ichol.trans %*% t.ichol

	### print( dim( t.sk.weights ))
	### print( c(t(t.pred.covmat) %*% solve(t.support.covmat ))[1:2] )
	
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
	    return( list( prediction = t.result[[1]], sqrt.mspe = t.result[[2]], sk.weights = t.sk.weights[1,]  ) )
	}
	if( method == 2 )
	{
	    return( list( prediction = t.result[[1]], sqrt.mspe = t.result[[2]],
		    P1 = t.result[[3]], Q1 =t.result[[4]], K = t.result[[5]], sk.weights = t.sk.weights[1,]  ) )
	}
	if( method == 3)
	{
	    return( list( prediction = t.result[[1]], mspe = t.result[[2]],
		    P1 = t.result[[3]], Q1 =t.result[[4]], K = t.result[[5]], sk.weights = t.sk.weights  ) )
	}

    },
    pixconfig = object@pixconfig,
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
if( method == 1 ){
    krige.result <- data.frame( matrix( NaN, nrow = nrow(t.pred.designmat), ncol = 2 ) )
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
    colnames( krige.result) = c("prediction", "prediction.se", "sqrt.P", "sqrt.Q", "K")
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
	    sk.weights.matrix <- t(matrix( unlist( lapply( t.krige.res, function( poly ){ return( poly$sk.weights ) } ) ), ncol = nrow( data ), byrow = T) )
	}else{
	    sk.weights.matrix <- t(matrix( unlist( lapply( t.krige.res, function( poly ){ return( poly$sk.weights[1,] ) } ) ), ncol = nrow( data ), byrow = T ) )
	    }
	## extended output is a list
	object <- SpatialPolygonsDataFrame( SpatialPolygons( object@polygons ), data = krige.result, match.ID = F)
	
	res <- list(
	    object = object,
	    krig.method = method,
	    parameter = list( beta.coef = t.beta.coef, cov.beta.coef = t.cov.beta.coef),
	    sk.weights = sk.weights.matrix,
	    inv.Sigma = t.isupport.covmat,
	    residuals = t.residuals
	)
	class( res )  <- "CKrige.exout.polygons"
	return( res )
    }
    else
    {
	object <- SpatialPolygonsDataFrame( SpatialPolygons( object@polygons ), data = krige.result, match.ID = F)
	P1.list = lapply( t.krige.res, function( polygon ){ return(polygon$P1)})
	Q1.list = lapply( t.krige.res, function( polygon ){ return(polygon$Q1)})
	K.list = lapply( t.krige.res, function( polygon ){ return(polygon$K)})
	sk.weights.list <-lapply( t.krige.res, function( polygon ){ return(t(polygon$sk.weights))}) 
       	res <- list(
	    object = object,
	    krig.method = method,
	    CMCK.par = list(P1 = P1.list, Q1 =Q1.list,K = K.list),
	    parameter = list( beta.coef = t.beta.coef, cov.beta.coef = t.cov.beta.coef),
	    sk.weights = sk.weights.list,
	    inv.Sigma = t.isupport.covmat,
	    residuals = t.residuals
	)
	class( res )  <- "CKrige.exout.polygons"
	return( res ) 
    }
}
else
{
    
    row.names( krige.result ) <- unlist( lapply( object@polygons, function(x){ return (x@ID ) } ) )
    res <- SpatialPolygonsDataFrame( SpatialPolygons( object@polygons ), krige.result )
    return( res )  
}
#    
} # end of function

