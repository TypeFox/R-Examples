#' Landscape Patch Statistics
#' 
#' \code{PatchStat} calculates the patch statistics for individual patches
#' identified in a matrix of data. The matrix can be a raster of class 'asc'
#' (adehabitat package), 'RasterLayer' (raster package) or
#' 'SpatialGridDataFrame' (sp package).
#' 
#' The patch statistics are based on statistics calculated by fragstats
#' \url{http://www.umass.edu/landeco/research/fragstats/fragstats.html}.
#' 
#' @param mat a matrix of data with individual patches identified as with
#' \code{ConnCompLabel}; The matrix can be a raster of class 'asc' (this &
#' adehabitat package), 'RasterLayer' (raster package) or
#' 'SpatialGridDataFrame' (sp package)
#' @param cellsize cell size (in meters) is a single value representing the
#' width/height of cell edges (assuming square cells)
#' @param latlon boolean value representing if the data is geographic. If
#' latlon == TRUE, matrix must be of class 'asc', 'RasterLayer' or
#' 'SpatialGridDataFrame'
#' @return a data.frame listing \item{patchID}{the unique ID for each patch.}
#' \item{n.cell}{the number of cells for each patch, specified in square
#' meters.} \item{n.core.cell}{the number of cells in the core area, without
#' the edge area.} \item{n.edges.perimeter}{the number of outer perimeter cell
#' edges of the patch.} \item{n.edges.internal}{the number of internal cell
#' edges of the patch.} \item{area}{the area of each patch comprising a
#' landscape mosaic.} \item{core.area}{represents the interior area of the
#' patch, greater than the specified depth-of-edge distance from the
#' perimeter.} \item{perimeter}{the perimeter of the patch, including any
#' internal holes in the patch, specified in meters.}
#' \item{perim.area.ratio}{the ratio of the patch perimeter (m) to area (m2).}
#' \item{shape.index}{the shape complexity, sum of each patches perimeter
#' divided by the square root of patch area.} \item{frac.dim.index}{fractal
#' dimension index reflects shape complexity across a range of spatial scales;
#' approaches 2 times the logarithm of patch perimeter (m) divided by the
#' logarithm of patch area (m2).} \item{core.area.index}{quantifies core area
#' as a percentage of patch area.}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @seealso \code{\link{ClassStat}}, \code{\link{ConnCompLabel}}
#' @references McGarigal, K., S. A. Cushman, M. C. Neel, and E. Ene. 2002.
#' FRAGSTATS: Spatial Pattern Analysis Program for Categorical Maps. Computer
#' software program produced by the authors at the University of Massachusetts,
#' Amherst. Available at the following web site:
#' \url{www.umass.edu/landeco/research/fragstats/fragstats.html}
#' @examples
#' 
#' 
#' #define a simple binary matrix
#' tmat = { matrix(c( 0,0,0,1,0,0,1,1,0,1,
#'                    0,0,1,0,1,0,0,0,0,0,
#'                    0,1,NA,1,0,1,0,0,0,1,
#'                    1,0,1,1,1,0,1,0,0,1,
#'                    0,1,0,1,0,1,0,0,0,1,
#'                    0,0,1,0,1,0,0,1,1,0,
#'                    1,0,0,1,0,0,1,0,0,1,
#'                    0,1,0,0,0,1,0,0,0,1,
#'                    0,0,1,1,1,0,0,0,0,1,
#'                    1,1,1,0,0,0,0,0,0,1),nr=10,byrow=TRUE) }
#' 					
#' #do the connected component labelling
#' ccl.mat = ConnCompLabel(tmat)
#' ccl.mat
#' image(t(ccl.mat[10:1,]),col=c('grey',rainbow(length(unique(ccl.mat))-1)))
#' 
#' #calculate the patch statistics
#' ps.data = PatchStat(ccl.mat)
#' ps.data
#' 
#' 
#' @export PatchStat
#' @useDynLib SDMTools projectedPS geographicPS
PatchStat <- 
function(mat,cellsize=1,latlon=FALSE)	{
	##method to calculate shape index or aggregation indexes
	#a = area of the patch in number of cells
	#p is the perimeter in number of edges
	shape.index <- 
	function(a,p) {
		n = trunc(sqrt(a))
		m = a - n^2
		minp=rep(0,length(m))
		for (ii in 1:length(m)){
			if (m[ii]==0) minp[ii] = 4*n[ii]
			if (n[ii]^2<a[ii] & a[ii]<=n[ii]*(1+n[ii])) minp[ii] = 4 * n[ii] + 2
			if (a[ii] > n[ii]*(1+n[ii])) minp[ii] = 4 * n[ii] + 4
		}
		return(p/minp)
	}

	#check if raster from sp or raster package and convert if necessary
	if (any(class(mat) %in% 'RasterLayer')) mat = asc.from.raster(mat)
	if (any(class(mat) == 'SpatialGridDataFrame')) mat = asc.from.sp(mat)
	#if latlon data
	if (latlon){
		if (!any(class(mat) == 'asc')) stop('matrix must be of class asc, RasterLayer or SpatialGridDataFrame... see helpfile')
		#get the cell size info
		cellinfo = grid.info(getXYcoords(mat)$y,attr(mat,'cellsize'))
		#check to ensure matrix
		mat = try(as.matrix(mat))
		#get the unique patch ID's
		ID.vals = as.numeric(na.omit(unique(as.vector(mat))));ID.vals = ID.vals[order(ID.vals)]
		#extract the base patch info
		out = as.data.frame(.Call('geographicPS',mat,ID.vals,cellinfo$area,cellinfo$top,cellinfo$bottom,cellinfo$side,PACKAGE='SDMTools'))
		names(out) = c('patchID','n.cell','n.core.cell','n.edges.perimeter','n.edges.internal','area','core.area','perimeter')
	} else {
		#check to ensure matrix
		mat = try(as.matrix(mat))
		if (!is.matrix(mat)) stop('objects must be a matrix')
		#get the unique patch ID's
		ID.vals = as.numeric(na.omit(unique(as.vector(mat))));ID.vals = ID.vals[order(ID.vals)]
		#extract the base patch info
		out = as.data.frame(.Call('projectedPS',mat,ID.vals,PACKAGE='SDMTools'))
		names(out) = c('patchID','n.cell','n.core.cell','n.edges.perimeter','n.edges.internal')
		#calculate other stats
		out$area = out$n.cell * cellsize^2
		out$core.area = out$n.core.cell * cellsize^2
		out$perimeter = out$n.edges.perimeter * cellsize
	}
	out$perim.area.ratio = out$perimeter / out$area 
	out$shape.index = shape.index(out$n.cell,out$n.edges.perimeter)
	out$frac.dim.index = (2 * log(0.25 * out$perimeter)) / log(out$area)
	out$core.area.index = out$core.area / out$area
	return(out)
}
