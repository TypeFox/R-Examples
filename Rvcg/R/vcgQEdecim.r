#' Performs Quadric Edge Decimation on triangular meshes.
#' 
#' Decimates a mesh by adapting the faces of a mesh either to a target face
#' number, a percentage or an approximate mesh resolution (a.k.a. mean edge
#' length
#' 
#' This is basically an adaption of the cli tridecimator from vcglib
#' 
#' @param mesh Triangular mesh of class "mesh3d"
#' @param tarface Integer: set number of target faces.
#' @param percent Numeric: between 0 and 1. Set amount of reduction relative to
#' existing face number. Overrides tarface argument. 
#' @param edgeLength Numeric: tries to decimate according to a target mean edge
#' length. Under the assumption of regular triangles, the edges are half as
#' long by dividing the triangle into 4 regular smaller triangles.
#' @param topo logical: if TRUE, mesh topology is preserved.
#' @param quality logical: if TRUE, vertex quality is considered.
#' @param bound logical: if TRUE, mesh boundary is preserved.
#' @param optiplace logical: if TRUE, mesh boundary is preserved.
#' @param scaleindi logical: if TRUE, decimatiion is scale independent.
#' @param normcheck logical: if TRUE, normal directions are considered.
#' @param safeheap logical: if TRUE, safeheap update option enabled.
#' @param qthresh numeric: Quality threshold for decimation process.
#' @param boundweight numeric: Weight assigned to mesh boundaries.
#' @param normalthr numeric: threshold for normal check in radians.
#' @param silent logical, if TRUE no console output is issued.
#' @return Returns a reduced mesh of class mesh3d.
#' @author Stefan Schlager
#' @seealso \code{\link{vcgSmooth}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#'
#' 
#' data(humface)
#' ##reduce faces to 50% 
#' decimface <- vcgQEdecim(humface, percent=0.5)
#' ## view
#' \dontrun{
#' require(rgl)
#' shade3d(decimface, col=3)
#' 
#' ## some light smoothing
#' decimface <- vcgSmooth(decimface,iteration = 1)
#' } 
#' @export vcgQEdecim
vcgQEdecim <- function(mesh,tarface=NULL,percent=NULL,edgeLength=NULL, topo=FALSE,quality=TRUE,bound=FALSE, optiplace = TRUE, scaleindi = TRUE, normcheck = FALSE, safeheap =FALSE, qthresh=0.3, boundweight = 1, normalthr = pi/2,silent=FALSE)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        doit <- TRUE
        mesh <- meshintegrity(mesh,facecheck=TRUE)
        vb <- mesh$vb[1:3,,drop=FALSE]
        it <- mesh$it-1
        dimit <- ncol(it)
        outmesh <- list()
        class(outmesh) <- "mesh3d"
        
        if (is.null(tarface) && is.null(percent)&& is.null(edgeLength))
            stop("please enter decimation option")
        if (!is.null(percent))
            {
                if (percent <= 0 || percent > 1)
                    stop ("percent must be between 0 and 1")
                tarface <- floor(percent*dimit)
            }
        if (!is.null(edgeLength))
            {
                currEL <- vcgMeshres(mesh)$res
                if (currEL >= edgeLength)
                    warning("edges already shorter than required - nothing to do")
                coef <- (currEL/edgeLength)^2
                tarface <- floor(coef*dimit)
                
            }
        ##concatenate parameters
        boolparams <- c(topo, quality, bound, optiplace, scaleindi, normcheck, safeheap)
        storage.mode(boolparams) <- "logical"
        doubleparams <- c(qthresh, boundweight, normalthr)
        storage.mode(doubleparams) <- "double"
###tmp <- .C("RQEdecim",vb,ncol(vb),it,ncol(it),tarface,vb,as.integer(topo),as.integer(quality),as.integer(bound))
        tmp <- .Call("RQEdecim", vb, it, tarface, boolparams, doubleparams,silent)
        outmesh$vb <- rbind(tmp$vb,1)
        
        outmesh$it <- tmp$it
        outmesh$normals <- rbind(tmp$normals, 1)
                                        #outmesh <- adnormals(outmesh)
        if (!is.null(edgeLength) && !silent)
            cat(paste("Mean Edge length is",vcgMeshres(outmesh)$res,"\n"))
        return(outmesh)
    }
