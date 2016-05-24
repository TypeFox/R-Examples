## computeV call the function K if the class attribut in computeV is
## set as "special" for more ?computeV of the  R-package spatialCovariance
K <- function( dist, model )
{
    t.cov <- CovarianceFct( x = dist, model = model$model,
	   		    param = c(NA, variance = 1,
			    nugget = 0, scale = model$scale, model$parameter )
		    )
    return( t.cov )
}
# # # 
# # # 
# # # 
f.pp.cov <- function(t.dist, model)
### purpose: calculate the  spatial covariance for a 
###          given distance vector with the cov model in model
###
### arguments:
###            t.dist = numeric vector with distances
###            model = list with covariance parameter
###                        one element of the list = 
###                        $model = cov model
###                        $scale = range parameter
###                        $parameter = vector with the cov model parameter
###                        $variance = sill matrix of the cov model
###
### author: ch.hofer       
### date: 20.2.2006
{

### function CovarianceFct is package RandomFields

t.cov <- 0

t.part.cov.list <- lapply( model, 
    function(model, t.dist)
    {
	if( model$model != "special" )
	{
	    t.part.cor <- CovarianceFct(
		x = t.dist,
		model = model$model,
		param = c(
		    mean = NA,
		    variance = 1,
		    nugget = 0,
		    scale.par = model$scale,
		    cov.par = model$parameter
		)
	    )
	}
	else
	{
	    t.part.cor <- K(
	    dist = t.dist, 
	    model = model)
    }
    return( t.part.cov <- t.part.cor * model$variance )
}, t.dist
)

t.part.cov.mat <- matrix( unlist( t.part.cov.list ), ncol = length( model ) )
rm( t.part.cov.list ) 
return( rowSums( t.part.cov.mat ) )

} ## end function f.pp.cov
