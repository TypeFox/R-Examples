##' A simple mesh generator for non-convex regions in n-D space
##' 
##' An unstructured simplex requires a choice of meshpoints (vertex nodes) and
##' a triangulation.  This is a simple and short algorithm that improves the
##' quality of a mesh by relocating the meshpoints according to a relaxation
##' scheme of forces in a truss structure. The topology of the truss is reset
##' using Delaunay triangulation. A (sufficiently smooth) user supplied signed
##' distance function (\code{fd}) indicates if a given node is inside or
##' outside the region. Points outside the region are projected back to the
##' boundary.
##' 
##' This is an implementation of original Matlab software of Per-Olof Persson.
##' 
##' Excerpt (modified) from the reference below:
##' 
##' \sQuote{The algorithm is based on a mechanical analogy between a triangular
##' mesh and a n-D truss structure. In the physical model, the edges of the
##' Delaunay triangles of a set of points correspond to bars of a truss. Each
##' bar has a force-displacement relationship \eqn{f(\ell, \ell_{0})}{F(L,L0)}
##' depending on its current length \eqn{\ell}{L} and its unextended length
##' \eqn{\ell_{0}}{L0}.}
##' 
##' \sQuote{External forces on the structure come at the boundaries, on which
##' external forces have normal orientations. These external forces are just
##' large enough to prevent nodes from moving outside the boundary. The
##' position of the nodes are the unknowns, and are found by solving for a
##' static force equilibrium. The hope is that (when \code{fh = function(p)
##' return(rep(1,nrow(p)))}), the lengths of all the bars at equilibrium will
##' be nearly equal, giving a well-shaped triangular mesh.}
##' 
##' See the references below for all details. Also, see the comments in the
##' source file of \code{distmesh2d}.
##' 
##' @param fdist Vectorized signed distance function, for example
##' \code{\link{mesh.dsphere}}, accepting an \code{m}-by-\code{n}
##' matrix, where \code{m} is arbitrary, as the first argument.
##' @param fh Vectorized function, for example \code{\link{mesh.hunif}},
##' that returns desired edge length as a function of position.
##' Accepts an \code{m}-by-\code{n} matrix, where \code{n} is
##' arbitrary, as its first argument.
##' @param h Initial distance between mesh nodes.
##' @param box \code{2}-by-\code{n} matrix that specifies the bounding box.
##' (See \link{distmesh2d} for an example.)
##' @param pfix \code{nfix}-by-2 matrix with fixed node positions.
##' @param \dots parameters that are passed to \code{fdist} and \code{fh}
##' @param ptol Algorithm stops when all node movements are smaller than
##' \code{dptol}
##' @param ttol Controls how far the points can move (relatively) before a
##' retriangulation with \code{\link{delaunayn}}.
##' @param deltat Size of the time step in Eulers method.
##' @param geps Tolerance in the geometry evaluations.
##' @param deps Stepsize \eqn{\Delta x} in numerical derivative computation for
##' distance function.
##' @return \code{m}-by-\code{n} matrix with node positions.
##' @section Wishlist : \itemize{ \item*Implement in C/Fortran \item*Translate
##' other functions of the matlab package }
##' @author Raoul Grasman; translated from original Matlab sources of Per-Olof
##' Persson.
##' @seealso \code{\link{distmesh2d}}, \code{\link[tripack]{tri.mesh}},
##' \code{\link{delaunayn}}, \code{\link{mesh.dsphere}},
##' \code{\link{mesh.hunif}},\cr \code{\link{mesh.diff}},
##' \code{\link{mesh.union}}, \code{\link{mesh.intersect}}
##' @references \url{http://persson.berkeley.edu/distmesh/}
##' 
##' \cite{P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB. SIAM
##' Review, Volume 46 (2), pp. 329-345, June 2004}
##' @keywords math optimize dplot graphs
##' @examples
##' 
##' \dontrun{
##' # examples distmeshnd
##' require(rgl)
##' 
##' fd = function(p, ...) sqrt((p^2)%*%c(1,1,1)) - 1
##'      # also predefined as `mesh.dsphere'
##' fh = function(p,...)  rep(1,nrow(p))
##'      # also predefined as `mesh.hunif'
##' bbox = matrix(c(-1,1),2,3)
##' p = distmeshnd(fd,fh,0.2,bbox, maxiter=100)
##'     # this may take a while:
##'     # press Esc to get result of current iteration
##' }
##'
##' @export
"distmeshnd"  <-
function (fdist, fh, h, box, pfix = array(dim = c(0, ncol(box))),
    ..., ptol = 0.001, ttol = 0.1, deltat = 0.1, geps = 0.1 *
        h, deps = sqrt(.Machine$double.eps) * h)
{
# %DISTMESHND N-D Mesh Generator using Distance Functions.
    dim = ncol(as.matrix(box))
    L0mult = 1 + 0.4/2^(dim - 1)
    rownorm2 = function(x) drop(sqrt((x^2) %*% rep(1, ncol(x))))

    # %1. Create initial distribution in bounding box
    if (dim == 1) {
        p = seq(box[1], box[2], by = h)
    }
    else {
        cbox = lapply(1:dim, function(ii) seq(box[1, ii], box[2,
            ii], by = h))
        p = do.call("expand.grid", cbox)
        p = as.matrix(p)
    }

    # %2. Remove points outside the region, apply the rejection method
    p = p[fdist(p, ...) < geps, ]
    r0 = fh(p, ...)
    p = rbind(pfix, p[stats::runif(nrow(p)) < min(r0)^dim/r0^dim, ])
    N = nrow(p)
    if (N <= dim + 1)
        stop("Not enough starting points inside boundary (is h0 too large?).")
    on.exit(return(invisible(p)))

    cat("Press esc if the mesh seems fine but the algorithm hasn't converged.\n")
    utils::flush.console()
    count = 0

    p0 = 1/.Machine$double.eps

    # mimick Matlab call ``localpairs=nchoosek(1:dim+1,2)'':
    localpairs = as.matrix(expand.grid(1:(dim + 1), 1:(dim + 1)))
    localpairs = localpairs[lower.tri(matrix(TRUE, dim + 1, dim + 1)), 2:1]

    while (TRUE) {
        if (max(rownorm2(p - p0)) > ttol * h) {
            # %3. Retriangulation by Delaunay:

            p0 = p
            t = delaunayn(p)
            pmid = matrix(0, nrow(t), dim)
            for (ii in 1:(dim + 1)) pmid = pmid + p[t[, ii],
                ]/(dim + 1)
            t = t[fdist(pmid, ...) < (-geps), ]
            pair = array(dim = c(0, 2))
            for (ii in 1:nrow(localpairs)) {
                pair = rbind(pair, t[, localpairs[ii, ]])
            }

            # %4. Describe each edge by a unique pair of nodes
            pair = Unique(pair, TRUE); # base-function `unique' is way too slow
            if (dim == 2) {
                trimesh(t, p[, 1:2])
            }
            else if (dim == 3) {
                if (count%%5 == 0) {
                  tetramesh(t, p)
                }
            }
            else {
                cat("Retriangulation #", 15, "\n")
                utils::flush.console()
            }
            count = count + 1
        }
        bars = p[pair[, 1], ] - p[pair[, 2], ]
        L = rownorm2(bars)
        L0 = fh((p[pair[, 1], ] + p[pair[, 2], ])/2, ...)
        L0 = L0 * L0mult * (sum(L^dim)/sum(L0^dim))^(1/dim)
        F = L0 - L
        F[F < 0] = 0
        Fbar = cbind(bars, -bars) * matrix(F/L, nrow = nrow(bars),
            ncol = 2 * dim)
        ii = pair[, t(matrix(1:2, 2, dim))]
        jj = rep(1, nrow(pair)) %o% c(1:dim, 1:dim)
        s = c(Fbar)
        ns = length(s)
        dp = matrix(0, N, dim)
        dp[1:(dim * N)] = rowsum(s, ii[1:ns] + ns * (jj[1:ns] -
            1))
        if (nrow(pfix) > 0)
            dp[1:nrow(pfix), ] = 0
        p = p + deltat * dp
        d = fdist(p, ...)
        ix = d > 0
        gradd = matrix(0, sum(ix), dim)
        for (ii in 1:dim) {
            a = rep(0, dim)
            a[ii] = deps
            d1x = fdist(p[ix, ] + rep(1, sum(ix)) %o% a, ...)
            gradd[, ii] = (d1x - d[ix])/deps
        }
        p[ix, ] = p[ix, ] - (d[ix] %o% rep(1, dim)) * gradd
        maxdp = max(deltat * rownorm2(dp[d < (-geps), ]))
        if (maxdp < ptol * h)
            break
    }
}
