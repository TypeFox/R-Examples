#' Finite Volume Community Ocean Model grid
#'
#' \code{fvcom.grid} provides a represenation of the unstructured triangular
#' sigma grid used by the Finite Volume Community Ocean Model (FVCOM). As a
#' disclaimer, please note that the author of this package is a user of FVCOM,
#' but is not affiliated with its development.
#'
#' @section Slots:
#'   \describe{
#'     \item{nodes.n}{The number of nodes in the grid.}
#'     \item{nodes.x}{x-coordinates of the nodes (m).}
#'     \item{nodes.y}{y-coordinates of the nodes (m).}
#'     \item{nodes.h}{z-coordinates (depth) of the nodes (m).}
#'     \item{nodes.lat}{Latitude of the nodes.}
#'     \item{nodes.lon}{Longitude of the nodes.}
#'     \item{elems.n}{Number of elements in the grid}
#'     \item{elems.v1}{1st set of node indices}
#'     \item{elems.v2}{2nd set of node indices}
#'     \item{elems.v3}{3rd set of node indices}
#'     \item{elems.size}{The area of each element (m^2)}.
#'     \item{proj}{A string that could be passed to proj4::project to convert
#'                 x and y values to latitude and longitude.}
#'   }
#' 
#' @references http://fvcom.smast.umassd.edu/FVCOM/
#' @name fvcom.grid-class
#' @rdname fvcom.grid-class
#' @exportClass fvcom.grid
setClass("fvcom.grid",
         representation(
             nodes.n="integer",
             nodes.x="numeric",
             nodes.y="numeric",
             nodes.h="numeric",
             nodes.lat="numeric",
             nodes.lon="numeric",
             elems.n="integer",
             elems.v1="integer",
             elems.v2="integer",
             elems.v3="integer",
             elems.size="numeric",
             proj="character"
             )
         )

#' An example instance of the \code{fvcom.grid} class.
#'
#' A dataset containing an example of the \code{fvcom.grid} class. This
#' example was generated for the purposes of demonstrating the plotting
#' functions of this package and is not an actual grid. The grid uses
#' bathymetry from NOAA's ETOPO1 dataset, which is available at
#' (http://www.ngdc.noaa.gov/mgg/global/global.html).
#'
#' @docType data
#' @keywords datasets
#' @format A \code{fvcom.grid} instance.
#' @name ocean.demo.grid
NULL
 
#' Load an FVCOM grid from a NetCDF output file.
#'
#' Loads enough of an FVCOM grid to use the other methods associated with
#' the \code{fvcom.grid} class. The (x,y,h) locations of the nodes and their
#' connections to form an unstructured triangular mesh are loaded.
#' @param filename The name of an output NetCDF file from FVCOM 2.7.
#' @param proj A projection string which could be passed to
#'   \code{proj4::project} and converts between the x,y and lat,lon coordinate
#'   systems for this grid.
#' @return An instance of the \code{fvcom.grid} class.
loadFVCOMGrid27 <- function(filename, proj) {
    ncid <- nc_open(filename)
    x <- as.vector(ncvar_get(ncid, 'x'))
    y <- as.vector(ncvar_get(ncid, 'y'))
    h <- as.vector(ncvar_get(ncid, 'h'))
    if(proj != '') {
        ll <- project(data.frame(x, y), proj, inverse=TRUE)
    } else {
        ll <- data.frame(x=NULL, y=NULL)
    }
    ev <- ncvar_get(ncid, 'nv') ## Element vertices
    elem.size <- sapply(seq(nrow(ev)), function(i)
        det(matrix(c(x[ev[i, 1]], x[ev[i, 2]], x[ev[i, 3]],
                     y[ev[i, 1]], y[ev[i, 2]], y[ev[i, 3]],
                     1, 1, 1),
                   3, 3, byrow=T))) / 2
    nc_close(ncid)
    return(new("fvcom.grid",
               nodes.n=length(x), nodes.x=x, nodes.y=y, nodes.h=h,
               nodes.lat=ll$y, nodes.lon=ll$x, proj=proj,
               elems.n=nrow(ev),
               elems.v1=ev[,1], elems.v2=ev[,2], elems.v3=ev[,3]))
}

#' Checks if the points (xy$x, xy$y) are in the \code{fvcom.grid} \code{grid}.
#'
#'  Given a grid and a point or set of points, \code{find.elem} returns
#'  the element number of each point in the grid. NA is returned for any
#'  point that does not lie within the grid. \code{x} and \code{y} must be
#'  the same length.
#' 
#' @param grid A \code{fvcom.grid} instance.
#' @param xy A \code{data.frame} with components \code{x} and \code{y} with
#'           the x and y locations of the points.
#' @param units Either 'll' for latitude and longitude or 'm' for meters.
#' @return A vector of integer values of length \code{nrow(xy)}. The values
#'         are indexes into the elements of grid.
#'
#' @examples {
#' # Create a regular matrix grid
#' nodes = get.nodes(ocean.demo.grid)
#' lattice.grid = expand.grid(x=seq(min(nodes$x), max(nodes$x), len=50),
#'                            y=seq(min(nodes$y), max(nodes$y), len=50))
#' elems = find.elem(ocean.demo.grid, lattice.grid, units="m")
#' # Plot the result
#' plot(lattice.grid$x, lattice.grid$y, pch=15,
#'      col=jet.colors(max(elems, na.rm=TRUE) * 2)[elems+2])
#' }
#'
#' @name find.elem
#' @aliases find.elem,fvcom.grid-method
#' @docType methods
#' @rdname find.elem-methods
setGeneric("find.elem", function(grid, xy, ...) {})
setMethod("find.elem", "fvcom.grid",
function(grid, xy, units='ll') {
    if(units == 'll') {
        grid.x <- grid@nodes.lon        
        grid.y <- grid@nodes.lat
    } else if(units == 'm') {
        grid.x <- grid@nodes.x
        grid.y <- grid@nodes.y
    } else {
        stop("Invalid units, must be 'll' or 'm'.")
    }
    elems <- integer(nrow(xy))
    ## TODO Look into ways to avoid copying the node/element data.
    elems <- .C('R_find_element', PACKAGE='ocean',
                x_pts=as.double(xy$x), y_pts=as.double(-xy$y), n_pts=nrow(xy),
                x=as.double(grid.x), y=as.double(-grid.y),
                n_grid_pts=as.integer(grid@elems.n),
                v1=as.integer(grid@elems.v1 - 1), ## TODO Subtraction in
                v2=as.integer(grid@elems.v2 - 1), ## C code
                v3=as.integer(grid@elems.v3 - 1),
                elems=as.integer(elems))$elems + 1
    ## TODO Need to modify C code to return elem + 1 except if elem = -1
    ## TODO C code should use R NA
    #print(sum(elems == 0))
    elems[elems == 0] <- NA
    return(elems)
}
)
#' Get the depth at each node in the grid.
#'
#' @param grid A \code{fvcom.grid} instance.
#' @return A vector of length \code{get.nnodes(grid)} with the depth at each node.
#' 
#' @name get.depth
#' @aliases get.depth,fvcom.grid-method
#' @docType methods
#' @rdname get.depth-methods
setGeneric("get.depth", function(grid, ...) {})
setMethod("get.depth", "fvcom.grid",
function(grid)
    return(grid@nodes.h)
)

#' Get the indices of the element vertices in the grid.
#'
#' @param grid A \code{fvcom.grid} instance.
#' @return A \code{data.frame} with \code{v1}, \code{v2}, and \code{v3}
#' elements that correspond to the indices of the nodes at vertices 1, 2, and
#' 3.
#'
#' @name get.elems
#' @aliases get.elems,fvcom.grid-method
#' @docType methods
#' @rdname get.elems-methods
setGeneric("get.elems", function(grid, ...) {})
setMethod("get.elems", "fvcom.grid",
function(grid)
    return(data.frame(v1=grid@elems.v1, v2=grid@elems.v2, v3=grid@elems.v3))
)

#' Get the number of elements in the grid.
#'
#' @param grid A \code{fvcom.grid} instance.
#' @return The number of elements in \code{grid}.
#'
#' @name get.nelems
#' @aliases get.nelems,fvcom.grid-method
#' @docType methods
#' @rdname get.nelems-methods
setGeneric("get.nelems", function(grid, ...) {})
setMethod("get.nelems", "fvcom.grid",
function(grid)
    return(grid@elems.n)
)

#' Get the number of nodes in the grid.
#'
#' @param grid A \code{fvcom.grid} instance.
#' @return The number of nodes in \code{grid}.
#'
#' @name get.nnodes
#' @aliases get.nnodes,fvcom.grid-method
#' @docType methods
#' @rdname get.nnodes-methods
setGeneric("get.nnodes", function(grid, ...) {})
setMethod("get.nnodes", "fvcom.grid",
function(grid)
    return(grid@nodes.n)
)

#' Get the values of the nodes in the grid.
#'
#' @param grid A \code{fvcom.grid} instance.
#' @return A \code{data.frame} with \code{x}, \code{y}, \code{h}, \code{lon},
#' and \code{lat} elements that correspond to the x, y, depth, longitude, and
#' latitude of each node.
#'
#' @name get.nodes
#' @aliases get.nodes,fvcom.grid-method
#' @docType methods
#' @rdname get.nodes-methods
setGeneric("get.nodes", function(grid, ...) {})
setMethod("get.nodes", "fvcom.grid",
function(grid)
    return(data.frame(x=grid@nodes.x, y=grid@nodes.y, h=grid@nodes.h,
                      lon=grid@nodes.lon, lat=grid@nodes.lat))
)

#' Get the value of the projection string.
#'
#' @param grid A \code{fvcom.grid} instance.
#' @return The value of the projection string specified when the grid was
#'         created.
#'
#' @name get.proj
#' @aliases get.proj,fvcom.grid-method
#' @docType methods
#' @rdname get.proj-methods
setGeneric("get.proj", function(grid, ...) {})
setMethod("get.proj", "fvcom.grid",
function(grid)
    return(grid@proj)
)

#' Convert a single scalar or node based quantity to element based. 
#' 
#' The length of \code{x} determines how it will be treated. If \code{x} has
#' length 1, it is returned as a single color. If the length of \code{x} is
#' the number of nodes in the grid, its value for each element is calculated
#' as the average of the value at the adjoining nodes. If the length of
#' \code{x} is the number of elements in the grid, it is returned as is. Any
#' other values throw an error.
#' 
#' @param grid An fvcom.grid instance.
#' @param x A vector of length 1 or \code{length get.nnodes(grid)}
#' @return A vector of length \code{get.nelems(grid)}
#'
#' @name interp
#' @aliases interp,fvcom.grid-method
#' @docType methods
#' @rdname interp-methods
setGeneric("interp", function(grid, ...) {})
setMethod("interp", "fvcom.grid",
function(grid, x) {
    if(length(x) == 1) {
        ## Repeat the scalar value nelems times.
        x.ret <- rep(x, get.nelems(grid))
    } else if(length(x) == get.nnodes(grid)) {
        ## Average the values at each node to get the value in each element.
        x.ret <- sapply(seq(get.nelems(grid)), function(i)
                        return(mean(c(x[grid@elems.v1[i]],
                                      x[grid@elems.v2[i]],
                                      x[grid@elems.v3[i]]))))
    } else if(length(x) == get.nelems(grid)) {
        x.ret <- x
    } else {
        stop(paste('length(x) must be equal to 1, the number of nodes, or ',
                   'the number of elements\n', sep=''))
    }
    return(x.ret)
}
)

#' Plot a \code{fvcom.grid} instance as a heatmap.
#'
#' Given a vector of data, this function plots the data as a heatmap on
#' an unstructured grid. The length of the data vector must be as long as
#' either the number nodes in the grid or the number of elements in the
#' grid. The grid is currently stored as a data.frame, but will be
#' converted to an S4 object in the future.
#' 
#' @param x A \code{fvcom.grid} instance.
#' @param z A vector to plot as a heatmap.
#' @param units Either 'll' for latitude and longitude or 'm' for meters.
#' @param col A list of colors, such as that returned by bathy.colors.
#' @param add Should the plot be added to the current plot?
#' #' @param xlim x-limits for the plot.
#' @param ylim y-limits for the plot.
#' @param zlim z-limits for the plot. 
#' @param border.col Color of the element borders. If not provided the borders
#'                   will be colored to match the adjacent polygons.
#' @param border.lwd Line width of the element borders.
#' @param bg.col Color of the background. The background is only plotted if
#'               add=F, otherwise bg.col is ignored.
#'
#' @examples {
#'   op = par(ask=TRUE)
#'   # Plot the grid in a single color
#'   image(ocean.demo.grid, col='white')
#'   # Plot the grid in bathy colors
#'   image(ocean.demo.grid)
#'   par(op)
#' }
#'
#' @name image
#' @aliases image,fvcom.grid-method
#' @docType methods
#' @rdname image-methods
setMethod("image", "fvcom.grid",
function(x, z=get.depth(x), units='ll',
         col=bathy.colors(100), add=FALSE,
         xlim=NA, ylim=NA, zlim=NA,
         border.col=NA, bg.col='gray', border.lwd=1) {
    ## TODO Set aspect ratio automatically
    grid = x
    ## Check the parameters for validity.
    ## Convert the passed in values to an element based quantity.
    z <- interp(grid, z)
    ## Calculate z limits
    if(!is.na(zlim)[1]) {
        ## If zlimits were provided, any values that lie outside of them are
        ## converted to NA.
        z[z < zlim[1]] <- NA
        z[z > zlim[2]] <- NA
    } else {
        ## If no z limits were provided, calculate them as the minimum and
        ## maximum of z ignoring missing values.
        zlim <- c(min(z, na.rm=TRUE), max(z, na.rm=TRUE))
    }
    ## Get the lat/lon or x/y locations of the nodes surrounding the elements.
    ## c(rbind()) interleaves the vectors
    if(units == "ll") { # Lat/lon
        x <- grid@nodes.lon[c(rbind(grid@elems.v1, grid@elems.v2,
                                    grid@elems.v3))]
        y <- grid@nodes.lat[c(rbind(grid@elems.v1, grid@elems.v2,
                                    grid@elems.v3))]
    } else if(units == "m") {
        x <- grid@nodes.x[c(rbind(grid@elems.v1, grid@elems.v2,
                                  grid@elems.v3))]
        y <- grid@nodes.y[c(rbind(grid@elems.v1, grid@elems.v2,
                                  grid@elems.v3))]
    } else {
        stop("Invalid units, options are 'll' or 'm'.")
    }
    ## Calculate the colors for the polygons and the borders
    n.cols <- length(col)
    cols.idx <- floor((z - zlim[1]) / (zlim[2] - zlim[1]) * n.cols)
    cols.idx[cols.idx <= 0] <- 1
    cols.idx[cols.idx > n.cols] <- n.cols
    col <- col[cols.idx]
    if(is.na(border.col)[1])
        border.col <- col
    ## cut.poly adds NAs after each polygon, ensuring that the polygons are
    ## not connected to one another by spurious lines.
    cut.poly <- function(x, n.sides=3) {
        n <- length(x) %/% n.sides
        polys <- rep(NA, n*n.sides)
        polys[rep(c(rep(T, n.sides), F), n)] <- x
        return(polys)
    }
    ## Set the x and y limits.
    if(is.na(xlim[1]))
        xlim <- c(min(x), max(x))
    if(is.na(ylim[1]))
        ylim <- c(min(y), max(y))
    ## Do the actual plotting
    if(!add)
        image(matrix(1,1,1), col=bg.col,
              xlab = "Longitude",
              ylab = "Latitude",
              xlim=xlim, ylim=ylim)
    polygon(cut.poly(x), cut.poly(y), col=col, border=border.col, lwd=border.lwd)
    #if(legend) (TODO Add legend)
    #    legend(min(x), max(y),
    #           legend=c(round(zlim[2], 2), rep("", n.cols-2),
    #           round(zlim[1], 2)),
    #           bt="n",
    #           col=rev(cols),
    #           y.intersp=10/n.cols,
    #           pch=15,
    #           cex=2)
}
)

#' Checks if the points (xy$x, xy$y) are in the \code{fvcom.grid} \code{grid}.
#'
#' @param grid A \code{fvcom.grid} instance.
#' @param xy A \code{data.frame} with components \code{x} and \code{y} with
#'           the x and y locations of the points.
#' @param units Either 'll' for latitude and longitude or 'm' for meters.
#' @return A vector of logical values of length \code{nrow(xy)}. The ith
#'         element is \code{TRUE} if (\code{xy$x[i]}, \code{xy$y[i]}) is in
#'         \code{grid} and \code{FALSE} otherwise.
#'
#' @examples {
#' # Create a regular grid of test points
#' nodes = get.nodes(ocean.demo.grid)
#' lattice.grid = expand.grid(x=seq(min(nodes$x), max(nodes$x), len=50),
#'                            y=seq(min(nodes$y), max(nodes$y), len=50))
#' # Check which points are in the grid
#' in.grid = is.in.grid(ocean.demo.grid, lattice.grid, units="m")
#' # Plot the points that are in the grid
#' with(lattice.grid[in.grid,], plot(x, y))
#' }
#'
#' @name is.in.grid
#' @aliases is.in.grid,fvcom.grid-method
#' @docType methods
#' @rdname is.in.grid-methods
setGeneric("is.in.grid", function(grid, xy, ...) {})
setMethod("is.in.grid", "fvcom.grid",
isInFVCOMGrid <- function(grid, xy, units='ll')
    return(!is.na(find.elem(grid, xy, units)))
)

#' Plot the density of x and y on grid.
#' 
#' Plots the distribution of x and y on grid. This function follows the
#' method described in Simons et al 2013: first vertically integrating the
#' data, then dividing by the number of particles spawned, and finally
#' applying an Gaussian blur filter. The coordinates should be passed in as
#' x,y in meters, then are inverse projected into lat/long using the PROJ.4
#' library.  
#' 
#' @param grid A \code{fvcom.grid} instance.
#' @param xy A \code{data.frame} with x and y location of points. NAs are not
#' supported.
#' @param npoints The number of points to scale the density plot by. This
#'                defaults to the number of points passed in, but it may be
#'                useful to set it to a different value if only a subset of
#'                the points are being plotted (e.g. some points are outside
#'                of the domain).
#' @param res The resolution of the plot in the same dimensions as \code{xy}
#'            is given. Square boxes will be plotted with each side of length
#'            \code{res}.
#' @param sigma The standard deviation of the Gaussian smoothing filter to be
#'              applied. If no filter is required, set sigma=0. The units
#'              should be the same as for \code{res}.
#' @param log Should the density be log10 transformed before plotting?
#' @param bg.col The background color.
#' @param col A list of colors, such as that returned by heat.colors.
#' @param add Should the plot be added to the current plot?
#' @param xlim x-limits for the plot.
#' @param ylim y-limits for the plot.
#' @param lim.units Units for xlim and ylim. One of 'm' (meters) or 'll'
#'                  (latitude and longitude).
#' @param zlim z-limits for the plot.
#'
#' @examples {
#' # Generate artificial data from a Gaussian distribution
#' nodes = get.nodes(ocean.demo.grid)
#' set.seed(1)
#' x = rnorm(50000, mean(nodes$x), sd(nodes$x))
#' y = rnorm(50000, mean(nodes$y), sd(nodes$y))
#' # Plot the density of the mixture
#' res = (max(nodes$x) - min(nodes$x)) / 50
#' grd = pdd(ocean.demo.grid, data.frame(x=x, y=y), res=res, sigma=5)
#' }
#'
#' @references {
#' Simons, R.D. and Siegel, D.A. and Brown K.S. 2013 Model sensitivity
#' and robustness in the estimation of larval transport: A study of
#' particle tracking parameters \emph{J. Marine Systems} 119--120:
#' 19--29.
#' }
#'
#' @name pdd
#' @aliases pdd,fvcom.grid-method
#' @docType methods
#' @rdname pdd-methods
setGeneric("pdd", function(grid, xy, ...) {})
setMethod("pdd", "fvcom.grid",
function(grid, xy, npoints=nrow(xy), res=1000, sigma=0,
         log=F, bg.col='gray', col=heat.colors(100), add=F,
         xlim=c(min(get.nodes(grid)$x), max(get.nodes(grid)$x)),
         ylim=c(min(get.nodes(grid)$y), max(get.nodes(grid)$y)),
         lim.units = 'm',
         zlim=NA) {
    ## TODO Use a matrix with rownames and colnames attrs
    ## Create a lattice grid for calculating the density.
    if(lim.units == 'll') {
        lim.proj = project(data.frame(x=xlim, y=ylim), proj=get.proj(grid))
        xlim = lim.proj$x
        ylim = lim.proj$y
        rm(lim.proj)
    }
    grd <- list(x=seq(xlim[1], xlim[2], by=res),
                y=seq(ylim[1], ylim[2], by=res))
    grd$data <- matrix(0, nrow=length(grd$x) - 1, ncol=length(grd$y) - 1)
    ## Bin the data into the grid
    bin.data <- function(x, y, grid) {
        out <- matrix(.C('R_bin_data', PACKAGE='ocean',
                         as.double(x), as.double(y), as.integer(length(x)),
                         as.double(grid$x), as.double(grid$y),
                         as.integer(nrow(grid$data)),
                         as.integer(ncol(grid$data)),
                         data=as.integer(as.vector(grid$data)))$data,
                      nrow(grid$data))
        return(out)
    }
    grd$data <- bin.data(xy$x, xy$y, grd)
    if(sum(grd$data > 0) == 0)
        stop('No points were within the grid.')
    ## Rescale by the number of particles released
    grd$data <- grd$data / npoints
    ## Apply a 2D Gaussian filter with std dev = sigma
    ## TODO Add xy.units as an option
    if(sigma > 0)
        grd$data <- filter2d(grd$data, sigma / res)
    ## Log transform the data if necessary
    if(log) {
        grd$data[grd$data == 0] <- NA ## Because log10(0) = -Inf
        grd$data <- matrix(log10(grd$data), nrow(grd$data))
    }
    ## Project the grid x,y into lat/lon
    x.proj <- c(grd$x, rep(grd$x[1], length(grd$y)))
    y.proj <- c(rep(grd$y[1], length(grd$x)), grd$y)
    p <- project(data.frame(x=x.proj, y=y.proj),
                 proj=get.proj(grid), inverse=TRUE)
    grd$x <- p$x[seq(length(grd$x))]
    grd$y <- p$y[length(p$y) - rev(seq(length(grd$y))) + 1]
    ## Project xlim and ylim if necessary
    lim.proj <- project(data.frame(x=xlim, y=ylim),
                        proj=get.proj(grid), inverse=TRUE)
    xlim <- lim.proj$x
    ylim <- lim.proj$y
    ## Set up a land mask
    ## Calculate the center of each grid cell. A cell is considered to be part
    ## of the grid if its center lies on the grid.
    ## TODO Adjust boundary cells to match grid exactly.
    grd$x.cent <- sapply(seq(length(grd$x) - 1), function(i)
                         mean(c(grd$x[i], grd$x[i + 1])))
    grd$y.cent <- sapply(seq(length(grd$y) - 1), function(i)
                         mean(c(grd$y[i], grd$y[i + 1])))
    mask <- expand.grid(x=grd$x.cent, y=grd$y.cent)
    mask$on.grid <- is.in.grid(grid, data.frame(x=mask$x, y=mask$y))
    ## TODO Why recast this? Convert row to column major? If so, just reverse
    ## the expand.grid arguments.
    grd$mask <- matrix(mask$on.grid, nrow(grd$data))
    grd$data[!grd$mask] <- NA
    rm(mask)
    
    ## Calculate the plotting limits
    if(is.na(zlim[1]))
        zlim <- c(min(grd$data, na.rm=TRUE), max(grd$data, na.rm=TRUE))
    xlim = c(max(min(grd$x), xlim[1]), min(max(grd$x), xlim[2]))
    ylim = c(max(min(grd$y), ylim[1]), min(max(grd$y), ylim[2]))
    ## Do the actual plotting
    image(matrix(1, 1, 1), xlim=xlim, ylim=ylim,
          col=bg.col, xlab='Longitude', ylab='Latitude')
    image(grd$data, x=grd$x, y=grd$y, col=col, add=TRUE, zlim=zlim)
    return(grd)
}
)

#' Plot an instance of the \code{fvcom.grid} class and overlay the 
#' trajectories given by \code{xy}.
#' 
#' @param x A \code{fvcom.grid} instance
#' @param xy A \code{list} with matrices \code{x} and \code{y} components that
#'           contain the trajectories to plot. The columns of \code{xy$x} are
#'           plotted against the columns of \code{xy$y}, so each particle
#'           trajectory should be in a column and each time index in a row.
#' @param plot.units The units for plotting. Either 'm' for meters or 'll' for
#'                    latitude and longitude.
#' @param xy.units The units of \code{xy}. Either 'm' for meters or 'll' for
#'                 latitude and longitude.
#' @param ... Additional arguments to be passed to \code{matlines}.
#'
#' @examples {
#' # Create a set of random trajectories.
#' nodes = get.nodes(ocean.demo.grid)
#' set.seed(1)
#' x = apply(matrix(rnorm(1000, 0, 0.01), 250), 2, cumsum) - 69.5
#' y = apply(matrix(rnorm(1000, 0, 0.01), 250), 2, cumsum) + 42.5
#' lines(ocean.demo.grid, list(x=x, y=y), lty=1)
#' }
#'
#' @name lines
#' @aliases lines,fvcom.grid-method
#' @docType methods
#' @rdname lines-methods
#setGeneric("lines", function(x, xy, ...) {})
setMethod("lines", "fvcom.grid",
function(x, xy, plot.units='ll', xy.units='ll', ...) {
    ## Project xy if necessary
    if((xy.units == 'm') & (plot.units == 'll')) {
        xy.proj = project(data.frame(x=as.vector(xy$x), y=as.vector(xy$y)),
                          proj=get.proj(x), inverse=TRUE)
        xy$x = matrix(xy.proj$x, nrow=nrow(xy$x))
        xy$y = matrix(xy.proj$y, nrow=nrow(xy$y))
    } else if((xy.units == 'll') & (plot.units == 'm')) {
        xy.proj = project(data.frame(x=as.vector(xy$x), y=as.vector(xy$y)),
                          proj=get.proj(x), inverse=FALSE)
        xy$x = matrix(xy.proj$x, nrow=nrow(xy$x))
        xy$y = matrix(xy.proj$y, nrow=nrow(xy$y))
    }
    ## Plot the background, then plot the trajectories
    image(x, col='white', units=plot.units)
    matlines(xy$x, xy$y, ...)
}
)

#' Check if a \code{fvcom.grid} instance is valid.
#'
#' @param object A \code{fvcom.grid} instance to check for validity.
#' @return \code{TRUE} if object is well formed, otherwise \code{FALSE.}
validFVCOMGrid <- function(object) {
    if((length(object@nodes.x) != object@nodes.n) ||
       (length(object@nodes.y) != object@nodes.n) ||
       (length(object@nodes.h) != object@nodes.n))
       return("Number of x, y, and z coordinates of nodes must be equal.")
    else if((length(object@elems.v1) != object@elems.n) ||
            (length(object@elems.v2) != object@elems.n) ||
            (length(object@elems.v3) != object@elems.n))
        return("All elements must have 3 vertices.")
    else if((sum((c(object@nodes.x, object@nodes.y, object@nodes.h) == Inf))
             != 0 ||
             sum(is.na(c(object@nodes.x, object@nodes.y, object@nodes.h)))
             != 0))
        return("Nodes may not have Inf or NA coordinates")
    else if((sum((c(object@nodes.x, object@nodes.y, object@nodes.h) == Inf))
             != 0 ||
             sum(is.na(c(object@nodes.x, object@nodes.y, object@nodes.h)))
             != 0))
        return("Elems may not have Inf or NA indices")
    else return(TRUE)
}
setValidity("fvcom.grid", validFVCOMGrid)
            
#' Create a new FVCOM grid instance from a FVCOM NetCDF output file.
#' 
#' @param filename The name of an output NetCDF file from FVCOM.
#' @param proj A string to passed to proj4::project to convert x,y locations
#' to latitude and longitude.
#' @return An instance of the \code{fvcom.grid} class.
fvcom.grid <- function(filename, proj)
    return(loadFVCOMGrid27(filename, proj))

