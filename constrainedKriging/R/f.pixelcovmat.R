f.pixelcovmat<- function(
                    pixgrid,    ### grid information like mesh size positioin and so on..
                    model)### parameter of the covariance
## purpose: calculate the cov matrix for all rectangles within
##          the boundingbox
## arguments:
##
##
## author: Ch. Hofer
## date: 25. Jan 2007
{
t.dim  <- pixgrid$pixgridInfo$nrows * pixgrid$pixgridInfo$ncols
#pixcm = PixelKovarianzmatrix
pixcm <- matrix(0, nrow = t.dim , ncol = t.dim)
#
t.pixel.area <- pixgrid$pixgridInfo$rowwidth * pixgrid$pixgridInfo$colwidth
#
## build covariance matrix of the pixels defined by pixgrid
t.q <- 1
for(i in 1:length( model ) )
{
    if(model[[ i ]]$model == "nugget")
    {
	diag( pixcm ) <- rep( 1, length( diag( pixcm ) ) ) * (model[[ i ]]$variance / t.pixel.area) + diag( pixcm)
    }
    else
    {
	pixcm <-  computeV( pixgrid$pixgridInfo,
	    			  class = "special",
				  params = model[[ i ]]
			      ) * model[[ i ]]$variance + pixcm
        t.q <- t.q + 1
    } # end if
} # end for 
return( pixcm )
} ### end of function
