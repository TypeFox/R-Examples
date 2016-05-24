#' Subsamples points on a mesh surface
#'
#' Subsamples surface of a triangular mesh and returns a set of points located on that mesh
#' @param mesh triangular mesh of class 'mesh3d'
#' @param SampleNum integer: number of sampled points (see \code{details} below)
#' @param type character: seclect sampling type ("mc"=MonteCarlo Sampling, "pd"=PoissonDisk Sampling,"km"=kmean clustering)
#' @param MCsamp integer: MonteCarlo sample iterations used in PoissonDisk sampling.
#' @param geodes logical: maximise geodesic distance between sample points (only for Poisson Disk sampling)
#' @param strict logical: if \code{type="pd"} and the amount of coordinates exceeds \code{SampleNum},  the resulting coordinates will be subsampled again by kmean clustering to reach the requested number.
#' @details Poisson disk subsampling will not generate the exact amount of coordinates specified in \code{SampleNum}, depending on \code{MCsamp} the result will contain more or less coordinates.
#' @return sampled points
#' @examples
#' 
#' data(humface)
#' ss <- vcgSample(humface,SampleNum = 500, type="pd")
#' \dontrun{
#' require(rgl)
#' points3d(ss)
#' }
#' @export vcgSample
vcgSample <- function(mesh, SampleNum=100,type=c("km","pd","mc"),MCsamp=20,geodes=TRUE,strict=FALSE)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        type <- match.arg(type[1],c("km","pd","mc"))
        if (type == "mc")
            type <- 1
        else if (type == "pd")
            type <- 2
        else if (type == "km")
            type <- 3
        noit <- FALSE
        if (!is.matrix(mesh$it)) {
            warning("mesh has no faces, kmean clustering on vertices is used")
            noit <- TRUE
            type <- 3
        }
        if (type %in% 1:2) {
            mesh <- meshintegrity(mesh)
            vb <- mesh$vb[1:3,,drop=FALSE]
            it <- mesh$it - 1
            storage.mode(it) <- "integer"
            dimit <- dim(it)[2]
            dimvb <- dim(vb)[2]
            type <- as.integer(type)
            SampleNum <- as.integer(SampleNum)
            MCsamp <- as.integer(MCsamp)
            if (!is.logical(geodes) || (FALSE %in% is.integer(c(it,type, MCsamp, SampleNum))) || (FALSE %in% is.numeric(vb)))
                stop("Please provide sensible arguments!")
            tmp <- .Call("Rsample", vb, it, SampleNum, type, MCsamp, geodes)
            tmp <- t(tmp)
            if (strict && nrow(tmp) > SampleNum)
                tmp <- kmeans(tmp,centers=SampleNum, iter.max=100)$centers
        } else {
            tmp <- kmeans(t(mesh$vb[1:3,]),centers=SampleNum, iter.max=100)$centers
            if (!noit)
                tmp <- t(vcgClost(tmp, mesh)$vb[1:3,])
        }
        return(tmp)
    }
