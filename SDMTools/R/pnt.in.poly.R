#' Point in Polygon
#' 
#' \code{pnt.in.poly} works out if 2D points lie within the boundaries of a
#' defined polygon. \cr \cr \bold{Note:} Points that lie on the boundaries of
#' the polygon or vertices are assumed to be within the polygon.
#' 
#' The algorithm implements a sum of the angles made between the test point and
#' each pair of points making up the polygon. The point is interior if the sum
#' is 2pi, otherwise, the point is exterior if the sum is 0. This works for
#' simple and complex polygons (with holes) given that the hole is defined with
#' a path made up of edges into and out of the hole. \cr \cr This sum of angles
#' is not able to consistently assign points that fall on vertices or on the
#' boundary of the polygon. The algorithm defined here assumes that points
#' falling on a boundary or polygon vertex are part of the polygon.
#' 
#' @param pnts a 2-column matrix or dataframe defining locations of the points
#' of interest
#' @param poly.pnts a 2-column matrix or dataframe defining the locations of
#' vertices of the polygon of interest
#' @return A 3-column dataframe where the first 2 columns are the original
#' locations of the points. The third column (names pip) is a vector of binary
#' values where 0 represents points not with the polygon and 1 within the
#' polygon.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
#' 
#' 
#' #define the points and polygon
#' pnts = expand.grid(x=seq(1,6,0.1),y=seq(1,6,0.1))
#' polypnts = cbind(x=c(2,3,3.5,3.5,3,4,5,4,5,5,4,3,3,3,2,2,1,1,1,1,2),
#'     y=c(1,2,2.5,2,2,1,2,3,4,5,4,5,4,3,3,4,5,4,3,2,2))
#' 
#' #plot the polygon and all points to be checked
#' plot(rbind(polypnts, pnts))
#' polygon(polypnts,col='#99999990')
#' 
#' #create check which points fall within the polygon
#' out = pnt.in.poly(pnts,polypnts)
#' head(out)
#' 
#' #identify points not in the polygon with an X
#' points(out[which(out$pip==0),1:2],pch='X')
#' 
#' 
#' @export pnt.in.poly
#' @useDynLib SDMTools pip
pnt.in.poly <- function(pnts,poly.pnts)	
{
	#check if pnts & poly is 2 column matrix or dataframe
	pnts = as.matrix(pnts); poly.pnts = as.matrix(poly.pnts)
	if (!(is.matrix(pnts) & is.matrix(poly.pnts))) stop('pnts & poly.pnts must be a 2 column dataframe or matrix')
	if (!(dim(pnts)[2] == 2 & dim(poly.pnts)[2] == 2)) stop('pnts & poly.pnts must be a 2 column dataframe or matrix')
	
	#ensure first and last polygon points are NOT the same
	if (poly.pnts[1,1] == poly.pnts[nrow(poly.pnts),1] & poly.pnts[1,2] == poly.pnts[nrow(poly.pnts),2]) poly.pnts = poly.pnts[-1,]
	
	#run the point in polygon code
	out = .Call('pip',pnts[,1],pnts[,2],nrow(pnts),poly.pnts[,1],poly.pnts[,2],nrow(poly.pnts),PACKAGE='SDMTools')
	
	#return the value
	return(data.frame(pnts,pip=out))
}
