#' Animate the output of a RAMAS Metapop simulation.
#' 
#' Animate temporal habitat and abundance dynamics on a gridded landscape.
#'  
#' @param dat A \code{SpatialPointsDataFrame} object returned by 
#'   \code{\link{mp2xy}}.
#' @param habitat A RasterStack or RasterBrick object. The number of Raster 
#'   layers should equal the number of simulation time steps. Further, the 
#'   layers should be ordered temporally, and should correspond to the levels of
#'   \code{dat$time}. It is assumed that the time step between layers is 
#'   consistent (i.e., the interval between animation frames is constant).
#' @param outfile A character string giving the desired output path and 
#'   filename.
#' @param zlim A numeric vector of length 2 giving the lower and upper limits of
#'   the color scale indicating habitat quality. If this is not provided, a 
#'   pretty range will be calculated (though this will impact efficiency).
#' @param axes Logical. Should axes be drawn?
#' @param col.regions A \code{\link{colorRampPalette}} function that will be 
#'   used to generate the colour ramp for grids. If \code{NULL}, a default 
#'   colour ramp based on \code{\link{terrain.colors}} is used.
#' @param pt.col A \code{\link{colorRampPalette}} function that will be used
#'   to generate the colour ramp for points. These colours will be interpolated 
#'   into 100 colours, which indicate relative mean population size, ranging 
#'   from 1 (first element of the colour ramp) to the maximum mean population 
#'   size that exists in the simulation output. If \code{NULL}, a default colour
#'   ramp ranging from white to black is used.
#' @param pt.cex Point size magnifier, relative to the default.
#' @param height Numeric. The height of the animation, in pixels. Default is 
#'   800.
#' @param width Numeric. The width of the animation, in pixels. Default is 800.
#' @param interval The time interval of the animation, in seconds. Default is 
#'   0.05, i.e. 20 frames per second.
#' @param label Should a time step counter be plotted?
#' @param label.pos A vector of two numbers giving the normalised parent
#'   coordinates at which the time step counter will be plotted. The first
#'   number gives the x-coordinate (0 = left edge, 1 = right edge) and the
#'   second gives the y-coordinate (0 = bottom edge, 1 = top edge). The default
#'   value of \code{c(0.98, 0.05)} plots the label at the bottom right corner. 
#'   Ignored if \code{label} is \code{FALSE}.
#' @param label.just Justification of the time step counter label, relative to
#'   \code{label.pos}. See the description of \code{just} at
#'   \code{\link{grid.text}}. Ignored if \code{label} is \code{FALSE}.
#' @param label.cex Size of the time step counter text. Ignored if \code{label}
#'   is \code{FALSE}.
#' @param label.font Font face of the time step counter text. See
#'   \code{fontface} at \code{\link{gpar}} for available options. Ignored if
#'   \code{label} is \code{FALSE}.
#' @param overwrite Should \code{outfile} be overwritten if it already exists?
#' @return \code{NULL}. The animation is saved as an animated .gif file at the
#'   specified path (\code{outfile}).
#' @details \code{mp_animate} requires that either 
#'   \href{ImageMagick}{http://www.imagemagick.org} or 
#'   \href{GraphicsMagick}{http://www.graphicsmagick.org} are installed on the 
#'   system. See the documentation for 
#'   \code{animation::\link[animation]{im.convert}} for further details.
#'   
#'   An animated gif is created, with points indicating the location of
#'   populations with mean abundance greater than zero, overlaid upon a raster
#'   grid indicating habitat suitability. Relative population size is
#'   represented by point colour, with white corresponding to a population with 
#'   mean abundance between 0 and 1% of the maximum mean abundance across all
#'   populations and time steps, and black corresponding to the maximum mean
#'   abundance. Colours for intermediate values are scaled linearly. The colour
#'   key indicates carrying capacity, and corresponds to the colour of grid
#'   cells.
#'  
#'   An example of this function's use is provided in the vignette "Introduction 
#'   to mptools" (\code{vignette('intro', 'mptools')}).
#' @keywords spatial
#' @seealso \code{\link{mp2sp}}
#' @importFrom raster stack extent nlayers ymin ymax cellStats
#' @importFrom animation ani.options im.convert
#' @importFrom rasterVis levelplot
#' @importFrom sp sp.points
#' @importFrom grid grid.segments unit gpar grid.text grid.points
#' @importFrom grDevices colorRampPalette terrain.colors png dev.off
#' @importFrom latticeExtra layer
#' @importFrom graphics par
#' @importFrom methods is
#' @importFrom viridis viridis
#' @export
mp_animate <- function(dat, habitat, outfile, zlim, axes=FALSE, 
                       col.regions=NULL, pt.col=NULL, pt.cex=1, height=800, 
                       width=800, interval=0.05, label=TRUE, 
                       label.pos=c(0.98, 0.05), label.just='right', 
                       label.cex=1.5, label.font=2, overwrite=FALSE) {
  nl <- raster::nlayers(habitat)
  if(nl != length(unique(dat$time))) 
    stop('The number of unique levels of time in dat (',  
         length(unique(dat$time)), ') should be the same as the number of ',
         'layers of habitat (', nl, ').')
  if(label) {
    if(any(findInterval(label.pos, c(0, 1), rightmost.closed=TRUE) != 1))
      stop('label.pos should be a vector of two numeric values giving the ',
           'normalised parent coordinates (range: 0 to 1) at which time step ',
           'label text will be plotted.', call.=FALSE)  
    label.pos <- grid::unit(label.pos, 'npc')
  }
  if(!overwrite & file.exists(outfile)) 
    stop('File ', outfile, ' already exists.', call.=FALSE)
  if(!dir.exists(dirname(outfile))) 
    stop('Directory ', dirname(outfile), ' does not exist.', call.=FALSE)
  if (!methods::is(habitat, "Raster")) 
    stop('habitat must be a Raster* object.', call.=FALSE)
  if (is.null(col.regions)) col.regions <- viridis::viridis
  if (is.null(pt.col))
    pt.col <- grDevices::colorRampPalette(c("white", "black"))
  dat$N.scl <- ceiling(100 * dat$N/max(dat$N, na.rm=TRUE))
  
  message('Creating gif animation.')
  ylims <- c(raster::ymin(raster::extent(habitat) * 1.2), raster::ymax(habitat))
  if(is.null(zlim)) 
    zlim <- range(pretty(c(0, max(raster::cellStats(habitat, max)))))
  prefix <- paste(sample(letters, 6, replace=TRUE), collapse='')
  grDevices::png(sprintf('%s/%s_%%0%sd.png', tempdir(), prefix, nchar(nl)), 
                 type='cairo', width=width, height=height)
  plot_t <- function(i, z) {
    d <- dat[dat$time==z & dat$N > 0, ]
    p <- rasterVis::levelplot(
      habitat[[i]], margin=FALSE, col.regions=col.regions, ylim=ylims,
      at=seq(zlim[[1]], zlim[[2]], length.out=101), scales=list(draw=axes),
      useRaster=TRUE, colorkey=list(height=0.6)) +
      latticeExtra::layer(
        sp::sp.points(d, pch=21, col=1, fill=pt.col(101)[d$N.scl + 1], 
                      cex=pt.cex), 
        data=list(d=d, pt.cex=pt.cex, pt.col=pt.col)) +
      latticeExtra::layer(
        grid::grid.text(z, label.pos[1], label.pos[2], just=label.just, 
                        gp=grid::gpar(fontface=label.font, cex=label.cex),
                        draw=label), 
        data=list(z=z, label=label, label.pos=label.pos, 
                  label.just=label.just, label.font=label.font, 
                  label.cex=label.cex))
    print(p)
  }
  mapply(plot_t, seq_len(nl), as.character(unique(dat$time))) 
  grDevices::dev.off()
  oopt <- animation::ani.options(
    ani.width=width, ani.height=height, interval=0.05, ani.dev='png', 
    ani.type='png', nmax=nl, autobrowse=FALSE)
  on.exit(animation::ani.options(oopt), add=TRUE)
  ff <- list.files(tempdir(), sprintf('^%s.*\\.png$', prefix), full.names=TRUE)
  animation::im.convert(ff, output=outfile, extra.opts='-dispose Background', 
                        clean=TRUE)
  invisible(NULL)
}
