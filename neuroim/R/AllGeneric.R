roxygen <- function() NULL

#' Generic function to print an object 
#' @param x the object to print
#' @param ... additional arguments
#' @export 
#' @rdname print-methods
setGeneric(name="print", def=function(x, ...) standardGeneric("print"))

#' Generic function to extract data values of object 
#' @param x the object to get values from
#' @param ... additional arguments
#' @export 
#' @rdname values-methods
setGeneric(name="values", def=function(x, ...) standardGeneric("values"))


#' Generic function to load data from a data source
#' 
#' @param x a data source
#' @param ... additional arguments
#' @export 
#' @rdname loadData-methods
setGeneric(name="loadData", def=function(x, ...) standardGeneric("loadData"))

#' Generic function to apply a function to an object
#' @param x the object that is mapped
#' @param m the mapping object
#' @param ... additional arguments
#' @export 
#' @rdname map-methods
setGeneric(name="map", def=function(x, m, ...) standardGeneric("map"))


#' Generic function to extract the number of dimensions of an object
#' 
#' @param x n-dimensional object
#' @param ... additional arguments
#' @export 
#' @examples 
#' x = BrainSpace(c(10,10,10), c(1,1,1))
#' ndim(x) == 3
#' x = BrainSpace(c(10,10,10,3), c(1,1,1,1))
#' ndim(x) == 4
#' 
#' @rdname ndim-methods
setGeneric(name="ndim", def=function(x, ...) standardGeneric("ndim"))

#' Generic function to add a dimension to an object
#' @param x a dimensioned object
#' @param n the size of the dimension to add
#' @export 
#' @rdname addDim-methods
#' @examples 
#' x = BrainSpace(c(10,10,10), c(1,1,1))
#' x1 <- addDim(x, 10)
#' ndim(x1) == 4
#' dim(x1)[4] == 10
setGeneric(name="addDim", def=function(x, n) standardGeneric("addDim"))

#' Generic function to drop a dimension from an object
#' @param x a dimensioned object
#' @param dimnum the index of the dimension to drop
#' @export 
#' @rdname dropDim-methods
#' @examples
#' x = BrainSpace(c(10,10,10), c(1,1,1))
#' x1 <- dropDim(x)
#' ndim(x1) == 2
#' dim(x1)[2] == 10
setGeneric(name="dropDim", def=function(x, dimnum) standardGeneric("dropDim"))

#' Generic function to extract geometric properties of an image.
#' @param x the object to query, e.g. an instance of \code{BrainVolume} or \code{BrainVector}
#' @param ... additional arguments
#' @return an object representing the geometric space of the image of type \code{\linkS4class{BrainSpace}}
#' @export space
#' @examples
#' x = BrainSpace(c(10,10,10), c(1,1,1))
#' vol <- BrainVolume(rnorm(10*10*10), x)
#' identical(x,space(vol))
#' 
#' @rdname space-methods
setGeneric(name="space", def=function(x, ...) standardGeneric("space"))

#' Generic function to fill disjoint sets of values with the output of a function
#' @param x the object to split
#' @param fac the factor to split by
#' @param FUN the function to summarize the the sets
#' @return a new object where the original values have been replaced by the function output
#' @export 
#' @rdname splitFill-methods
#' @details \code{FUN} can either return a scalar for each input vector or a vector equal to the length of the input vector. 
#' If it returns a scalar then every voxel in the set will be filled with that value in the output vector.
#' @examples 
#' 
#' ## summarize with mean -- FUN returns a scalar
#' x = BrainSpace(c(10,10,10), c(1,1,1))
#' vol <- BrainVolume(rnorm(10*10*10), x)
#' fac <- factor(rep(1:10, length.out=1000))
#' ovol.mean <- splitFill(vol, fac, mean)
#' identical(dim(ovol.mean), dim(vol))
#' length(unique(as.vector(ovol.mean))) == 10
#' ## transform by reversing vector -- FUN returns a vector.
#' ovol2 <- splitFill(vol, fac, rev)
setGeneric(name="splitFill", def=function(x, fac, FUN) standardGeneric("splitFill"))

#' Generic function to map values from one set to another using a user-supplied lookup table
#' @param x the object to map values from
#' @param lookup the lookup table. The first column is the "key" the second column is the "value".
#' @return a new object where the original values have been filled in with the values in the lookup table
#' @export 
#' @examples 
#' x <- BrainSpace(c(10,10,10), c(1,1,1))
#' vol <- BrainVolume(sample(1:10, 10*10*10, replace=TRUE), x)
#' 
#' ## lookup table is list
#' lookup <- lapply(1:10, function(i) i*10)
#' ovol <- fill(vol, lookup)
#' 
#' ## lookup table is matrix. First column is key, second column is value
#' names(lookup) <- 1:length(lookup)
#' lookup.mat <- cbind(as.numeric(names(lookup)), unlist(lookup))
#' ovol2 <- fill(vol, lookup.mat)
#' all.equal(as.vector(ovol2), as.vector(ovol))
#' 
#' @rdname fill-methods
setGeneric(name="fill", def=function(x, lookup) standardGeneric("fill"))



#' Generic function to center/scale row-subsets of a matrix or matrix-like object
#' @param x a numeric matrix or matrix-like object
#' @param f the splitting object, typically a \code{factor} or set of \code{integer} indices. must be equal to number of rows of matrix.
#' @param center should values within each submatrix be centered? (mean removed from each column of submatrix)
#' @param scale should values be scaled? (divide vector by standard deviation from each column of submatrix)
#' @return a new matrix or matrix-like object where the original rows have been grouped by \code{f} and then centered and/or scaled for each grouping
#' @docType methods
#' @export 
#' @examples 
#' 
#' M <- matrix(rnorm(1000), 10, 100)
#' fac <- factor(rep(1:2, each=5))
#' Ms <- splitScale(M, fac)
#' 
#' ## correctly centered
#' all(abs(apply(Ms[fac == 1,], 2, mean)) < .000001)
#' all(abs(apply(Ms[fac == 2,], 2, mean)) < .000001)
#' 
#' # correctly scaled
#' all.equal(apply(Ms[fac == 1,], 2, sd), rep(1, ncol(Ms)))
#' all.equal(apply(Ms[fac == 2,], 2, sd), rep(1, ncol(Ms)))
#' @rdname splitScale-methods
setGeneric(name="splitScale", def=function(x, f, center, scale) standardGeneric("splitScale"))

#' Generic function to summarize subsets of an object by first splitting by row and then "reducing" by a summary \code{function}
#' @param x a numeric matrix(like) object
#' @param fac the factor to define subsets of the object
#' @param FUN the function to apply to each subset. if \code{FUN} is missing, than the mean of each sub-matrix column is computed.
#' @return a new \code{matrix} where the original values have been reduced
#' @docType methods
#' @details if \code{FUN} is supplied it must take a vector and return a single scalar value. If it returns more than one value, an error will occur.
#' 
#' if \code{x} is a \code{BrainVector} instance then voxels (dims 1:3) are treated as columns and time-series (dim 4) as rows. 
#' The summary function then is applied to groups of voxels. However, if the goal is to apply a function to groups of time-points, 
#' then this can be achieved as follows: 
#' 
#' \code{ splitReduce(t(as.matrix(bvec)), fac) }
#' 
#'
#' @export 
#' @examples 
#' mat = matrix(rnorm(100*100), 100, 100)
#' fac = sample(1:3, nrow(mat), replace=TRUE)
#' ## compute column means of each sub-matrix
#' ms <- splitReduce(mat, fac)
#' all.equal(row.names(ms), levels(fac))
#' 
#' ## compute column medians of each sub-matrix
#' ms <- splitReduce(mat, fac, median)
#' 
#' ## compute time-series means grouped over voxels. 
#' ## Here, \code{length(fac)} must equal the number of voxels: \code{prod(dim(bvec)[1:3]}
#' bvec <- BrainVector(array(rnorm(24*24*24*24), c(24,24,24,24)), BrainSpace(c(24,24,24,24), c(1,1,1)))
#' fac <- factor(sample(1:3, prod(dim(bvec)[1:3]), replace=TRUE))
#' ms <- splitReduce(bvec, fac)
#' ms2 <- splitReduce(bvec, fac, mean)
#' all.equal(row.names(ms), levels(fac))
#' all.equal(ms,ms2)
#' 
#' @rdname splitReduce-methods
setGeneric(name="splitReduce", def=function(x, fac, FUN) standardGeneric("splitReduce"))


#' Generic function to extract the voxel dimensions of an image
#' @param x the object
#' @return a numeric vector
#' @export 
#' @examples 
#' bspace <- BrainSpace(c(10,10,10), c(2,2,2))
#' all.equal(spacing(bspace), c(2,2,2))
#' @rdname spacing-methods
setGeneric(name="spacing", def=function(x) standardGeneric("spacing"))

#' Generic function to extract the spatial bounds (origin + dim * spacing) of an image
#' param x the object
#' @export
#' @param x the object with \code{bounds} property
#' @return a \code{matrix} where each row contains the min (column 1) and max (column 2) bounds of the image dimension from 1 to \code{ndim(image)}.
#' @examples 
#' bspace <- BrainSpace(c(10,10,10), c(2,2,2))
#' b <- bounds(bspace)
#' nrow(b) == ndim(bspace)
#' ncol(b) == 2
#' @rdname bounds-methods
setGeneric(name="bounds",     def=function(x) standardGeneric("bounds"))


#' Generic getter function to extract image axes
#' @param x an object with a set of axes
#' @export 
#' @rdname axes-methods
setGeneric(name="axes",  def=function(x) standardGeneric("axes"))

#' Generic getter to extract image origin
#' @param x an object with an origin
#' @export 
#' @examples 
#' bspace <- BrainSpace(c(10,10,10), c(2,2,2))
#' origin(bspace)
#' 
#' @rdname origin-methods
setGeneric(name="origin", def=function(x) standardGeneric("origin"))

#' Generic getter to extract image coordinate transformation
#' @param x an object with a transformation
#' @export 
#' @details 
#' This function returns a transformation that can be used to go from "grid coordinates" to "real world coordinates" in millimeters.
#' @details This function returns a transformation that can be used to go from "grid coordinates" to "real world coordinates" in millimeters.
#' see \code{\linkS4class{BrainSpace}}
#' @examples 
#' bspace <- BrainSpace(c(10,10,10), c(2,2,2))
#' trans(bspace)
#' all.equal(dim(trans(bspace)), c(4,4))
#' @rdname trans-methods
setGeneric(name="trans",  def=function(x) standardGeneric("trans"))

#' Generic getter to extract inverse image coordinate transformation
#' @param x an object
#' @export 
#' @examples 
#' bspace <- BrainSpace(c(10,10,10), c(2,2,2))
#' itrans <- inverseTrans(bspace)
#' identical(trans(bspace) %*% inverseTrans(bspace), diag(4))
#' @rdname inverseTrans-methods
setGeneric(name="inverseTrans", def=function(x) standardGeneric("inverseTrans"))

#' Generic function to read a sequence of elements from an input source
#' @param x the input channel
#' @param numElements the number of elements to read
#' @return the elements as a vector
#' @export 
#' @rdname readElements-methods
setGeneric(name="readElements", def=function(x, numElements) standardGeneric("readElements"))

#' Generic function to read a set of column vector from an input source (e.g. \code{ColumnReader})
#' @param x the input channel
#' @param columnIndices the column indices
#' @return a \code{matrix} consisting of the requested column vectors 
#' @export 
#' @rdname readColumns-methods
setGeneric(name="readColumns", def=function(x, columnIndices) standardGeneric("readColumns"))


#' Generic function to write a sequence of elements from an input source
#' @param x the output channel
#' @param els the elements to write
#' @export 
#' @rdname writeElements-methods
setGeneric(name="writeElements", def=function(x, els) standardGeneric("writeElements"))


#' Generic function to write a 3D image volume to disk
#' @param x an image object, typically a \code{BrainVolume} instance.
#' @param fileName output file name
#' @param format file format string. Since "NIFTI" is the only currently supported format, this parameter can be safely ignored and omitted.
#' @param dataType output data type, If specified should be a \code{character} vector of: "BINARY", "UBYTE", "SHORT", "INT", "FLOAT", "DOUBLE". 
#' Otherwise output format will be inferred from R the datatype of the image.
#' @export 
#' @details 
#'  
#'  The output format will be inferred from file extension.
#' @details The output format will be inferred from file extension.
#'  \code{writeVolume(x, "out.nii")} outputs a NIFTI file.
#'  \code{writeVolume(x, "out.nii.gz")} outputs a gzipped NIFTI file.
#'  
#' No other file output formats are currently supported.
#' 
#' 
#' @examples 
#' 
#' bvol <- BrainVolume(array(0, c(10,10,10)), BrainSpace(c(10,10,10), c(1,1,1)))
#' \dontrun{
#' writeVolume(bvol, "out.nii")
#' writeVolume(bvol, "out.nii.gz")
#' }
#' @rdname writeVolume-methods
setGeneric(name="writeVolume",  def=function(x, fileName, format, dataType) standardGeneric("writeVolume"))


#' Generic function to write a 4D image vector to disk
#' @param x an image object, typically a \code{BrainVector} instance.
#' @param fileName output file name.
#' @param format file format string. Since "NIFTI" is the only currently supported format, this parameter can be safely ignored and omitted.
#' @param dataType the numeric data type. If specified should be a \code{character} vector of: "BINARY", "UBYTE", "SHORT", "INT", "FLOAT", "DOUBLE". 
#' Otherwise output format will be inferred from R the datatype of the image.
#' @export 
#' @examples 
#' 
#' bvec <- BrainVector(array(0, c(10,10,10,10)), BrainSpace(c(10,10,10,10), c(1,1,1)))
#' \dontrun{
#' writeVector(bvol, "out.nii")
#' writeVector(bvol, "out.nii.gz")
#' writeVector(bvec, "out.nii")
#' writeVector(bvec, "out.nii.gz")
#' }
#' @rdname writeVector-methods
setGeneric(name="writeVector",  def=function(x, fileName, format, dataType) standardGeneric("writeVector"))


#' Generic function to convert 1D indices to N-dimensional grid coordinates
#' @param x the object
#' @param idx the 1D \code{vector} of indices
#' @return a matrix of grid coordinates
#' @export 
#' @examples 
#'  
#'  bvol <- BrainVolume(array(0, c(10,10,10)), BrainSpace(c(10,10,10), c(1,1,1)))
#'  idx <- 1:10
#'  g <- indexToGrid(bvol, idx)
#'  bvol[g]
#' 
#' @rdname indexToGrid-methods
setGeneric(name="indexToGrid",   def=function(x, idx) standardGeneric("indexToGrid"))

#' Generic function to convert 1D indices to N-dimensional real world coordinates
#' @param x the object
#' @param idx the 1D indices
#' @return a matrix of real coordinates
#' @export 
#' @examples 
#' bvol <- BrainVolume(array(0, c(10,10,10)), BrainSpace(c(10,10,10), c(1,1,1)))
#' idx <- 1:10
#' g <- indexToCoord(bvol, idx)
#' idx2 <- coordToIndex(bvol, g)
#' all.equal(idx, idx2)
#' @rdname indexToCoord-methods
setGeneric(name="indexToCoord",   def=function(x, idx) standardGeneric("indexToCoord"))

#' Generic function to convert N-dimensional real world coordinates to 1D indices
#' @param x the object
#' @param coords a matrix of real world coordinates
#' @return a vector of indices
#' @export 
#' @rdname coordToIndex-methods
setGeneric(name="coordToIndex",   def=function(x, coords) standardGeneric("coordToIndex"))

#' Generic function to convert N-dimensional real world coordinates to grid coordinates
#' @param x the object
#' @param coords a matrix of real world coordinates
#' @return a matrix of grid coordinates
#' @export 
#' @rdname coordToGrid-methods
setGeneric(name="coordToGrid",   def=function(x, coords) standardGeneric("coordToGrid"))

#' Generic function to convert N-dimensional grid coordinate coordinates to real world coordinates
#' Generic function to convert N-dimensional grid coordinates to real world coordinates
#' @param x the object
#' @param coords a matrix of grid coordinates
#' @return a matrix of real coordinates
#' @export 
#' @rdname gridToCoord-methods
setGeneric(name="gridToCoord",   def=function(x, coords) standardGeneric("gridToCoord"))



#' Generic function to convert 1-dimensional real axis coordinates along a single axis dimension to an 1D index along the same axis
#' @param x the object
#' @param real the axis coordinates
#' @param dimNum the dimension number of the axis (e.g.  1, 2, 3)
#' @return a vector of axis indices
#' @export 
#' @rdname axisToIndex-methods
setGeneric(name="axisToIndex",   def=function(x, real, dimNum) standardGeneric("axisToIndex"))

#' Generic function to convert N-dimensional grid coordinate to 1D indices
#' @param x the object, typically a \code{BrainVolume} or \code{BrainSpace} instance.
#' @param coords a matrix where each row is a coordinate or a vector of length equal to \code{ndim(x)}
#' @return a vector of indices
#' @export 
#' @rdname gridToIndex-methods
setGeneric(name="gridToIndex",   def=function(x, coords) standardGeneric("gridToIndex"))


#' Generic function to apply a function to each volume of a four-dimensional image
#' @param x four-dimensional image, e.g. of class \code{BrainVector}
#' @param FUN a \code{function} taking one or two arguments (depending on the value of \code{withIndex})
#' @param withIndex whether the index of the volume supplied as the second argument to the function
#' @param mask an image mask indicating subset of volume elements to apply function over
#' @param ... additional arguments
#' @return a \code{list} of results of apply \code{FUN} to each volume.
#' @export 
#' @examples 
#' bvec <- BrainVector(array(rnorm(24*24*24*24), c(24,24,24,24)), BrainSpace(c(24,24,24,24), c(1,1,1)))
#' res <- eachVolume(bvec, mean)
#' 
#' res <- eachVolume(bvec, function(x,i) median(x), withIndex=TRUE)
#' @rdname eachVolume-methods
setGeneric(name="eachVolume", def=function(x, FUN, withIndex, mask, ...) standardGeneric("eachVolume"))

#' Generic function to extract a one or more individual volumes from a four-dimensional image
#' @param x four-dimensional image
#' @param i the indices of the volume(s) to extract
#' @param ... additional arguments
#' @return a list of \code{BrainVolume} elements
#' @export 
#' @examples 
#' bvec <- BrainVector(array(rnorm(24*24*24*24), c(24,24,24,24)), BrainSpace(c(24,24,24,24), c(1,1,1)))
#' vol <- takeVolume(bvec,1)
#' all.equal(dim(vol), c(24,24,24))
#' 
#' vol <- takeVolume(bvec,1:3)
#' length(vol) == 3
#' class(vol) == "list"
#' @rdname takeVolume-methods
setGeneric(name="takeVolume", def=function(x, i, ...) standardGeneric("takeVolume"))

#' Generic function to extract a sub-vector from a \code{BrainVector} object.
#' @param x four-dimensional image
#' @param i the indices of the volume(s) to extract
#' @param ... additional arguments
#' @return a  \code{BrainVector} object that is a sub-vector of the supplied object.
#' @export 
#' @examples 
#' bvec <- BrainVector(array(rnorm(24*24*24*24), c(24,24,24,24)), BrainSpace(c(24,24,24,24), c(1,1,1)))
#' vec <- subVector(bvec,1:2)
#' all.equal(2, dim(vec)[4])
#' 
#' vec <- subVector(bvec, c(1,3,5,7))
#' all.equal(4, dim(vec)[4])
#' 
#' mask <- LogicalBrainVolume(rep(TRUE, 24*24*24), BrainSpace(c(24,24,24), c(1,1,1)))
#' svec <- SparseBrainVector(array(rnorm(24*24*24*24), c(24,24,24,24)), 
#' BrainSpace(c(24,24,24,24), c(1,1,1)), mask)
#' vec <- subVector(svec, c(1,3,5))
#' all.equal(3, dim(vec)[4])
#' @rdname subVector-methods
setGeneric(name="subVector", def=function(x, i, ...) standardGeneric("subVector"))



#' Generic functions to apply a function to each (2D) slice of an image
#' @param x the object
#' @param FUN a \code{function} taking one or two arguments (depending on the value of \code{withIndex}
#' @param withIndex whether the index of the slice is supplied as the second argument to the function
#' @param ... additional arguments
#' @export 
#' @rdname eachSlice-methods
setGeneric(name="eachSlice", def=function(x, FUN, withIndex, ...) standardGeneric("eachSlice"))

#' Generic functions to apply a function to each series of a 4D image
#' That is, if the 4th dimension is 'time' each series is a 1D time series.
#' @param x a four dimensional image
#' @param FUN a \code{function} taking one or two arguments (depending on the value of \code{withIndex}
#' @param withIndex whether the index of the series is supplied as the second argument to the function
#' @param ... additional arguments
#' @export 
#' @examples 
#' bvec <- BrainVector(array(rnorm(24*24*24*24), c(24,24,24,24)), BrainSpace(c(24,24,24,24), c(1,1,1)))
#' res <- eachSeries(bvec, mean)
#' length(res) == 24*24*24
#' @rdname eachSeries-methods
setGeneric(name="eachSeries", def=function(x, FUN, withIndex, ...) standardGeneric("eachSeries"))

#' Generic functions to scale (center and/or normalize by standard deviation) each series of a 4D image
#' That is, if the 4th dimension is 'time' each series is a 1D time series.
#' @param x a four dimensional image
#' @param center a \code{logical} value indicating whether series should be centered
#' @param scale a \code{logical} value indicating whether series should be divided by standard deviation
#' @export 
#' @examples 
#' bvec <- BrainVector(array(rnorm(24*24*24*24), c(24,24,24,24)), BrainSpace(c(24,24,24,24), c(1,1,1)))
#' res <- scaleSeries(bvec, TRUE, TRUE)
#' @rdname scaleSeries-methods
setGeneric(name="scaleSeries", def=function(x, center, scale) standardGeneric("scaleSeries"))



#' Generic function to extract a set of series from a 4D image
#' @param x a four dimensional image
#' @param indices the indices of the series' to extract
#' @param ... additional arguments
#' @export 
#' @rdname takeSeries-methods
setGeneric(name="takeSeries", def=function(x, indices, ...) standardGeneric("takeSeries"))



#' Convert to from dense to sparse representation
#' 
#' @param x the object to make sparse, e.g. \code{DenseBrainVolume} or \code{DenseBrainVector}
#' @param mask the elements to retain
#' @param ... additional arguments
#' 
#' @details \code{mask} can be an integer vector of 1D indices or a mask volume of class \code{LogicalBrainVolume}
#' @export
#' @examples 
#' bvol <- BrainVolume(array(runif(24*24*24), c(24,24,24)), BrainSpace(c(24,24,24), c(1,1,1)))
#' indmask <- sort(sample(1:(24*24*24), 100))
#' svol <- as.sparse(bvol, indmask)
#' 
#' 
#' mask <- LogicalBrainVolume(runif(length(indmask)), space=space(bvol), indices=indmask)
#' sum(mask) == 100
#' @rdname as.sparse-methods
setGeneric(name="as.sparse", def=function(x, mask, ...) standardGeneric("as.sparse"))

#' Convert to a LogicalBrainVolume
#' @param x the object to binarize
#' @param indices the indices to set to TRUE
#' @export 
#' @rdname as.mask-methods
setGeneric(name="as.mask", def=function(x, indices) standardGeneric("as.mask"))


#' tesselate
#' @param x the object to tesselate
#' @param K the number of partitions
#' @param ... extra arguments
#' @export 
#' @rdname tesselate-methods
setGeneric(name="tesselate", def=function(x, K, ...) standardGeneric("tesselate"))

#' partition
#' @param x the object to partition
#' @param K the number of partitions
#' @param features the features used to define the partition
#' @param ... additional arguments
#' @export 
#' @rdname partition-methods
setGeneric(name="partition", def=function(x, K, features, ...) standardGeneric("partition"))

#' mergePartitions
#' @param x the object to merge
#' @param K the number of merged partitions
#' @param features the features used to define the partition
#' @param ... additional arguments
#' @export 
#' @rdname mergePartitions-methods
setGeneric(name="mergePartitions", def=function(x, K, features, ...) standardGeneric("mergePartitions"))

#' numClusters
#' @param x the object to extract number of clusters 
#' @export 
#' @rdname numClusters-methods
setGeneric(name="numClusters", def=function(x) standardGeneric("numClusters"))

#' clusterCenters
#' @param x the object to extract cluster centers from
#' @param features additional features
#' @param FUN a user-supplied function
#' @export 
#' @rdname clusterCenters-methods
setGeneric(name="clusterCenters", def=function(x, features, FUN) standardGeneric("clusterCenters"))


#' pick
#' @param x the object to pick from
#' @param mask a mask object
#' @param ... addiitonal arguments
#' @export
#' @rdname pick-methods
setGeneric(name="pick", def=function(x, mask, ...) standardGeneric("pick"))


#' Extract coordinates
#' @param x the object to extract coordinates from
#' @param ... additional arguments
#' @export 
#' @rdname coords-methods
setGeneric(name="coords", def=function(x, ...) standardGeneric("coords"))

#' Extract indices
#' @param x the object to extract indices
#' @export 
#' @rdname indices-methods
setGeneric(name="indices", def=function(x) standardGeneric("indices"))

#' Index Lookup operation
#' @param x the object to query
#' @param i the index to lookup
#' @param ... additional arguments
#' @export 
#' @rdname lookup-methods
setGeneric(name="lookup", def=function(x, i, ...) standardGeneric("lookup"))

#' Extract vector series from object
#' @param x the object
#' @param i the series index
#' @param ... additional arguments
#' @export 
#' @rdname series-methods
setGeneric(name="series", def=function(x, i, ...) standardGeneric("series"))   

#' Extract a 2D slice from an image volume
#' @param x the object
#' @param zlevel coordinate (in voxel units) along the sliced axis
#' @param along the axis along which to slice
#' @param orientation the target orientation of the 2D slice
#' @param ... additional arguments
#' @export 
#' @rdname slice-methods
setGeneric(name="slice", def=function(x, zlevel, along, orientation, ...) standardGeneric("slice"))   

#' Render an image to create a drawable image.
#' @param x the object, e.g. an instance of type \code{BrainSlice}
#' @param width width of the rendered image
#' @param height height of the rendered image
#' @param colmap the colors used to map from values to RGBA colors.
#' @param ... additional arguments
#' @rdname render-methods
setGeneric(name="render", def=function(x, width, height, colmap,...) standardGeneric("render"))


#' Render a slice at z coordinate
#' @param x the object, e.g. an instance of type \code{Layer} or \code{Overlay}
#' @param zpos the z coordinate to slice through.
#' @param width width of the rendered image
#' @param height height of the rendered image
#' @param colmap the colors used to map from values to RGBA colors.
#' @param ... additional arguments
#' @rdname renderSlice-methods
setGeneric(name="renderSlice", def=function(x, zpos, width, height, colmap,...) standardGeneric("renderSlice"))




#' Extract permutation matrix
#' @param x the object
#' @param ... additional arguments
#' @export 
#' @rdname permMat-methods
setGeneric(name="permMat", def=function(x, ...) standardGeneric("permMat"))   

#' Concatenate two objects
#' @param x the first object, typically \code{BrainVolume} or \code{BrainVector}
#' @param y the second object, typically \code{BrainVolume} or \code{BrainVector}
#' @details The \code{x} and \code{y} images must have compatible dimensions. a \code{BrainVolume} can be concatenated to \code{BrainVector}, and vice versa. See examples.
#' @param ... additional objects
#' 
#' @examples 
#' bv1 <- BrainVolume(rep(1,1000), BrainSpace(c(10,10,10), c(1,1,1)))
#' bv2 <- BrainVolume(rep(2,1000), BrainSpace(c(10,10,10), c(1,1,1)))
#' bv3 <- concat(bv1,bv2)
#' inherits(bv3, "BrainVector")
#' 
#' bv4 <- concat(bv3, bv1)
#' dim(bv4)[4] == 3
#' bv5 <- concat(bv1, bv3)
#' dim(bv4)[4] == 3
#' 
#' bv6 <- concat(bv4,bv5)
#' dim(bv6)[4] == 6
#' 
#' @export 
#' @rdname concat-methods
setGeneric(name="concat", def=function(x,y, ...) standardGeneric("concat"))

#' Find connected components
#' @name connComp
#' @param x the image object
#' @param ... additonal arguments
#' @export
#' @rdname connComp-methods
setGeneric(name="connComp", def=function(x, ...) standardGeneric("connComp"))

#' seriesIter
#' 
#' Construct a series iterator
#' @param x the object to be iterated over. This is typically an instance of class \code{\linkS4class{BrainVector}}
#' @return an \code{iter} object from the \code{iterators} package.
#' @export
#' @examples 
#' 
#' ## create a BrainVector with 10X10X10X10, where the last dimension is 
#' ## by convention the fourth dimension.
#' bvec <- BrainVector(array(rnorm(10*10*10*10), rep(10,4)), BrainSpace(rep(10,4), c(1,1,1)))
#' iter <- seriesIter(bvec)
#' 
#' ## compute mean of each series
#' library(foreach)
#' library(iterators)
#' foreach(i=iter, .combine=c) %do% { mean(i) }
#' iter <- seriesIter(bvec)
#' 
#' ## combine all series into a matrix
#' foreach(i=iter, .combine=rbind) %do% { i }
#' 
#' ## scale all series, add as columns in matrix.
#' foreach(i=seriesIter(bvec), .combine=cbind) %do% { scale(i) }
setGeneric(name="seriesIter", def=function(x) standardGeneric("seriesIter"))


#' extract voxel coordinates
#' @param x the object to extract voxels from
#' @param ... additional arguments to function
#' @export 
#' @rdname voxels-methods
setGeneric(name="voxels", def=function(x, ...) standardGeneric("voxels"))


if (!isGeneric("image"))
  setGeneric("image", function(x, ...) standardGeneric("image"))

if (!isGeneric("as.raster"))
  setGeneric("as.raster", function(x, ...) standardGeneric("as.raster"))

#' overlay two objects
#' @param x the underlay object
#' @param y the overlay object
#' @param ... additional arguments for class-specific implementations
#' @export 
#' @rdname overlay-methods
setGeneric("overlay", function(x, y, ...) standardGeneric("overlay"))


#' Generic function to test whether a file name conforms to the given \code{\linkS4class{BrainFileDescriptor}} instance.
#' Will test for match to either header file or data file
#' @param x object for which the file name is to matched to
#' @param fileName file name to be matched
#' @return TRUE for match, FALSE otherwise
#' @export fileMatches
#' @rdname fileMatches-methods
setGeneric(name="fileMatches", def=function(x, fileName) standardGeneric("fileMatches"))


#' Generic function to test whether a file name conforms to the given \code{\linkS4class{BrainFileDescriptor}} instance.
#' Will test for match to header file only
#' @param x object for which the file name is to matched to
#' @param fileName file name to be matched
#' @return TRUE for match, FALSE otherwise
#' @export headerFileMatches
#' @rdname headerFileMatches-methods
setGeneric(name="headerFileMatches", def=function(x, fileName) standardGeneric("headerFileMatches"))

#' Generic function to test whether a file name conforms to the given a \code{\linkS4class{BrainFileDescriptor}} instance.
#' Will test for match to data file only
#' @param x object for which the file name is to matched to
#' @param fileName file name to be matched
#' @return TRUE for match, FALSE otherwise
#' @export dataFileMatches
#' @rdname dataFileMatches-methods
setGeneric(name="dataFileMatches", def=function(x, fileName) standardGeneric("dataFileMatches"))

#' Generic function to get the name of the header file, given a file name and a \code{\linkS4class{BrainFileDescriptor}} instance.
#' @param x descriptor instance
#' @param fileName file name to be stripped of its extension
#' @return the correct header name
#' @export headerFile
#' @rdname headerFile-methods
setGeneric(name="headerFile", def=function(x, fileName) standardGeneric("headerFile"))

#' Generic function to get the name of the data file, given a file name and a \code{\linkS4class{BrainFileDescriptor}} instance.
#' @param x descriptor instance
#' @param fileName file name to be stripped of its extension
#' @return the correct header name
#' @export dataFile
#' @rdname dataFile-methods
setGeneric(name="dataFile", def=function(x, fileName) standardGeneric("dataFile"))

#' Generic function to strip extension from file name, given a \code{\linkS4class{BrainFileDescriptor}} instance.
#' @param x descriptor instance
#' @param fileName file name to be stripped of its extension
#' @return fileName without extension
#' @export stripExtension
#' @rdname stripExtension-methods
setGeneric(name="stripExtension", def=function(x, fileName) standardGeneric("stripExtension"))

#' Generic function to read image meta info given a file and a \code{\linkS4class{BrainFileDescriptor}} instance.
#' @param x descriptor instance
#' @param fileName file name contianing meta information
#' @export readMetaInfo
#' @rdname readMetaInfo-methods
setGeneric(name="readMetaInfo", def=function(x, fileName) standardGeneric("readMetaInfo"))


