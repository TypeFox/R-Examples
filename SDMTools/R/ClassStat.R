#' Landscape Class Statistics
#' 
#' \code{ClassStat} calculates the class statistics for patch types identified
#' in a matrix of data or in a raster of class 'asc' (SDMTools & adehabitat
#' packages), 'RasterLayer' (raster package) or 'SpatialGridDataFrame' (sp
#' package).
#' 
#' The class statistics are based on statistics calculated by fragstats
#' \url{http://www.umass.edu/landeco/research/fragstats/fragstats.html}.
#' 
#' @param mat a matrix of data with patches identified as classes (unique
#' integer values) as e.g., a binary lanscape of a species distribution or a
#' vegetation map. Matrix can be a raster of class 'asc' (adehabitat package),
#' 'RasterLayer' (raster package) or 'SpatialGridDataFrame' (sp package)
#' @param cellsize cell size (in meters) is a single value representing the
#' width/height of cell edges (assuming square cells)
#' @param bkgd the background value for which statistics will not be calculated
#' @param latlon boolean value representing if the data is geographic. If
#' latlon == TRUE, matrix must be of class 'asc', 'RasterLayer' or
#' 'SpatialGridDataFrame'
#' @return a data.frame listing \item{class}{a particular patch type from the
#' original input matrix (\code{mat}).} \item{n.patches}{the number of patches
#' of a particular patch type or in a class.} \item{total.area}{the sum of the
#' areas (m2) of all patches of the corresponding patch type.}
#' \item{prop.landscape}{the proportion of the total lanscape represented by
#' this class} \item{patch.density}{the numbers of patches of the corresponding
#' patch type divided by total landscape area (m2).} \item{total.edge}{the
#' total edge length of a particular patch type.} \item{edge.density}{edge
#' length on a per unit area basis that facilitates comparison among landscapes
#' of varying size.} \item{landscape.shape.index}{a standardized measure of
#' total edge or edge density that adjusts for the size of the landscape.}
#' \item{largest.patch.index}{largest patch index quantifies the percentage of
#' total landscape area comprised by the largest patch.}
#' \item{mean.patch.area}{average area of patches.}
#' \item{sd.patch.area}{standard deviation of patch areas.}
#' \item{min.patch.area}{the minimum patch area of the total patch areas. }
#' \item{max.patch.area}{the maximum patch area of the total patch areas.}
#' \item{perimeter.area.frac.dim}{perimeter-area fractal dimension equals 2
#' divided by the slope of regression line obtained by regressing the logarithm
#' of patch area (m2) against the logarithm of patch perimeter (m).}
#' \item{mean.perim.area.ratio}{the mean of the ratio patch perimeter. The
#' perimeter-area ratio is equal to the ratio of the patch perimeter (m) to
#' area (m2).} \item{sd.perim.area.ratio}{standard deviation of the ratio patch
#' perimeter.} \item{min.perim.area.ratio}{minimum perimeter area ratio}
#' \item{max.perim.area.ratio}{maximum perimeter area ratio.}
#' \item{mean.shape.index}{mean of shape index} \item{sd.shape.index}{standard
#' deviation of shape index.} \item{min.shape.index}{the minimum shape index.}
#' \item{max.shape.index}{the maximum shape index.}
#' \item{mean.frac.dim.index}{mean of fractal dimension index.}
#' \item{sd.frac.dim.index}{standard deviation of fractal dimension index.}
#' \item{min.frac.dim.index}{the minimum fractal dimension index.}
#' \item{max.frac.dim.index}{the maximum fractal dimension index.}
#' \item{total.core.area}{the sum of the core areas of the patches (m2).}
#' \item{prop.landscape.core}{proportional landscape core}
#' \item{mean.patch.core.area}{mean patch core area.}
#' \item{sd.patch.core.area}{standard deviation of patch core area.}
#' \item{min.patch.core.area}{the minimum patch core area.}
#' \item{max.patch.core.area}{the maximum patch core area.}
#' \item{prop.like.adjacencies}{calculated from the adjacency matrix, which
#' shows the frequency with which different pairs of patch types (including
#' like adjacencies between the same patch type) appear side-by-side on the map
#' (measures the degree of aggregation of patch types).}
#' \item{aggregation.index}{computed simply as an area-weighted mean class
#' aggregation index, where each class is weighted by its proportional area in
#' the landscape.} \item{lanscape.division.index}{based on the cumulative patch
#' area distribution and is interpreted as the probability that two randomly
#' chosen pixels in the landscape are not situated in the same patch}
#' \item{splitting.index}{based on the cumulative patch area distribution and
#' is interpreted as the effective mesh number, or number of patches with a
#' constant patch size when the landscape is subdivided into S patches, where S
#' is the value of the splitting index.} \item{effective.mesh.size}{equals 1
#' divided by the total landscape area (m2) multiplied by the sum of patch area
#' (m2) squared, summed across all patches in the landscape.}
#' \item{patch.cohesion.index}{measures the physical connectedness of the
#' corresponding patch type.}
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @seealso \code{\link{PatchStat}}, \code{\link{ConnCompLabel}}
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
#' #calculate the class statistics
#' cl.data = ClassStat(tmat)
#' cl.data
#' 
#' #identify background data is 0
#' cl.data = ClassStat(tmat,bkgd=0)
#' cl.data
#' 
#' 
#' 
#' @export 
ClassStat  = function(mat,cellsize=1,bkgd=NA,latlon=FALSE) {
	##method to calculate shape index or aggregation indexes
	#a = area of the patch in number of cells
	#p is the perimeter in number of edges
	#g is the number of 'internal' edges (single count)
	aggregation.index = function(a,g) {
		n = trunc(sqrt(a))
		m = a - n^2
		if (m==0) maxg = 2*n*(n-1)
		if (m<=n) maxg = 2*n*(n-1)+2*m-1
		if (m>n) maxg = 2*n*(n-1)+2*m-2
		minp=rep(0,length(m))
		for (ii in 1:length(m)){
			if (m[ii]==0) minp[ii] = 4*n[ii]
			if (n[ii]^2<a[ii] & a[ii]<=n[ii]*(1+n[ii])) minp[ii] = 4 * n[ii] + 2
			if (a[ii] > n[ii]*(1+n[ii])) minp[ii] = 4 * n[ii] + 4
		}
		return((g/maxg)*100)
	}
	shape.index = function(a,p) {
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

	######
	#check if raster from sp or raster package and convert if necessary
	if (any(class(mat) %in% 'RasterLayer')) mat = asc.from.raster(mat)
	if (any(class(mat) == 'SpatialGridDataFrame')) mat = asc.from.sp(mat)
	#check to ensure matrix
	mat = try(as.matrix(mat))
	if (!is.matrix(mat)) stop('objects must be a matrix')
	#get the uniqu classes of data
	classes = as.numeric(na.omit(unique(as.vector(mat))));classes = classes[order(classes)]
	#omit the background value
	if (!is.na(bkgd)) classes = classes[-which(classes==bkgd)]
	#out is the final object to be returned
	out = NULL
	#cycle through each of the classes
	for (cl in classes){
		#create a reclassed matrix
		mat2 = mat; mat2 = mat * 0; mat2[which(mat==cl)] = 1
		#get the patch info for the class
		out.patch = PatchStat(ConnCompLabel(mat2),cellsize=cellsize,latlon=latlon);rm(mat2)
		#define a couple constants
		L.cell = sum(out.patch$n.cell) #n cells in landscape
		L.area = sum(out.patch$area) #full area of landscape
		#remove the background patch (id = 0)
		if (0 %in% out.patch$patchID) out.patch = out.patch[-which(out.patch$patchID==0),]		
		#create a temporary variable to store output & calculate patch stats
		tout = list(class=cl)
		tout$n.patches = nrow(out.patch)
		tout$total.area = sum(out.patch$area)
		tout$prop.landscape = tout$total.area / L.area
		tout$patch.density = tout$n.patches / L.area
		tout$total.edge = sum(out.patch$perimeter)
		tout$edge.density = tout$total.edge / L.area
		tout$landscape.shape.index = shape.index(sum(out.patch$n.cell),sum(out.patch$n.edges.perimeter))
		tout$largest.patch.index = max(out.patch$area) / L.area
		tout$mean.patch.area = mean(out.patch$area)
		tout$sd.patch.area = sd(out.patch$area)
		tout$min.patch.area = min(out.patch$area)
		tout$max.patch.area = max(out.patch$area)
		tout$perimeter.area.frac.dim = 2 / (((tout$n.patches*sum(log(out.patch$perimeter)+log(out.patch$area)))-(tout$total.edge*tout$total.area))/(tout$n.patches*sum(log(out.patch$perimeter^2))-tout$total.edge^2))
		tout$mean.perim.area.ratio = mean(out.patch$perim.area.ratio)
		tout$sd.perim.area.ratio = sd(out.patch$perim.area.ratio)
		tout$min.perim.area.ratio = min(out.patch$perim.area.ratio)
		tout$max.perim.area.ratio = max(out.patch$perim.area.ratio)
		tout$mean.shape.index = mean(out.patch$shape.index,na.rm=T)
		tout$sd.shape.index = sd(out.patch$shape.index,na.rm=T)
		tout$min.shape.index = min(out.patch$shape.index,na.rm=T)
		tout$max.shape.index = max(out.patch$shape.index,na.rm=T)
		tout$mean.frac.dim.index = mean(out.patch$frac.dim.index,na.rm=T)
		tout$sd.frac.dim.index = sd(out.patch$frac.dim.index,na.rm=T)
		tout$min.frac.dim.index = min(out.patch$frac.dim.index,na.rm=T)
		tout$max.frac.dim.index = max(out.patch$frac.dim.index,na.rm=T)	
		tout$total.core.area = sum(out.patch$core.area)
		tout$prop.landscape.core = tout$total.core.area / L.area
		tout$mean.patch.core.area = mean(out.patch$core.area)
		tout$sd.patch.core.area = sd(out.patch$core.area)
		tout$min.patch.core.area = min(out.patch$core.area)
		tout$max.patch.core.area = max(out.patch$core.area)
		tout$prop.like.adjacencies = sum(out.patch$n.edges.internal) / sum(out.patch$n.edges.internal+out.patch$n.edges.perimeter*2)
		tout$aggregation.index = aggregation.index(sum(out.patch$n.cell),sum(out.patch$n.edges.internal)/2)
		tout$lanscape.division.index = 1-sum((out.patch$n.cell / L.cell)^2)
		tout$splitting.index = L.area^2 / sum(out.patch$area^2)
		tout$effective.mesh.size = sum(out.patch$area^2) / L.area 
		tout$patch.cohesion.index = ((1-(sum(out.patch$n.edges.internal)/sum(out.patch$n.edges.internal*sqrt(out.patch$n.cell))) )*((1-1/sqrt(L.cell))/10))*100
		
		#store in out 
		out = rbind(out,as.data.frame(tout))
	}
	return(out)
}
