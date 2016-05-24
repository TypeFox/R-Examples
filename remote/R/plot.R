if ( !isGeneric('plot') ) {
  setGeneric('plot', function(x, y, ...)
    standardGeneric('plot'))
}


#' Plot an Eot* object
#' 
#' @description
#' This is the standard plotting routine for the results of \code{\link{eot}}.
#' Three panels will be drawn i) the predictor domain, ii) the response 
#' domain, iii) the time series at the identified base point
#' 
#' @param x either an object of EotMode or EotStack as returned by \code{\link{eot}}
#' @param y integer or character of the mode to be plotted (e.g. 2 or "mode_2")
#' @param pred.prm the parameter of the predictor to be plotted.\cr
#' Can be any of "r", "rsq", "rsq.sums", "p", "int" or "slp"
#' @param resp.prm the parameter of the response to be plotted.\cr
#' Can be any of "r", "rsq", "rsq.sums", "p", "int" or "slp"
#' @param show.bp logical. If \code{TRUE} a grey circle will be drawn 
#' in the predictor image to indicate the location of the base point
#' @param anomalies logical. If \code{TRUE} a reference line will be drawn
#' a 0 in the EOT time series
#' @param add.map logical. If \code{TRUE} country outlines will be added 
#' to the predictor and response images
#' @param ts.vec an (optional) time series vector of the considered 
#' EOT calculation to be shown as the x-axis in the time series plot
#' @param arrange whether the final plot should be arranged in "wide" or
#' "long" format
#' @param clr an (optional) color palette for displaying of the 
#' predictor and response fields
#' @param locations logical. If x is an EotStack, set this to TRUE to 
#' produce a map showing the locations of all modes. Ignored if x is an
#' EotMode
#' @param ... further arguments to be passed to \code{\link[raster]{spplot}}
#' 
#' @examples
#' data(vdendool)
#' 
#' ## claculate 2 leading modes
#' nh_modes <- eot(x = vdendool, y = NULL, n = 2, 
#'                 reduce.both = FALSE, standardised = FALSE, 
#'                 verbose = TRUE)
#'
#' ## default settings 
#' plot(nh_modes, y = 1) # is equivalent to
#'
#' \dontrun{
#' plot(nh_modes[[1]]) 
#' 
#' plot(nh_modes, y = 2) # shows variance explained by mode 2 only
#' plot(nh_modes[[2]]) # shows cumulative variance explained by modes 1 & 2
#' 
#' ## showing the loction of the mode
#' plot(nh_modes, y = 1, show.bp = TRUE)
#' 
#' ## changing parameters
#' plot(nh_modes, y = 1, show.bp = TRUE,
#'      pred.prm = "r", resp.prm = "p")
#'         
#' ## change plot arrangement
#' plot(nh_modes, y = 1, show.bp = TRUE, arrange = "long") 
#' 
#' ## plot locations of all base points
#' plot(nh_modes, locations = TRUE)
#' }
#' 
#' @export
#' @name plot
#' @rdname plot
#' @aliases plot,EotMode,ANY-method

# set methods -------------------------------------------------------------

setMethod('plot', signature(x = 'EotMode',
                            y = 'ANY'), 
          function(x,
                   y,
                   pred.prm = "rsq",
                   resp.prm = "r",
                   show.bp = FALSE,
                   anomalies = TRUE,
                   add.map = TRUE,
                   ts.vec = NULL,
                   arrange = c("wide", "long"),
                   clr = NULL,
                   locations = FALSE,
                   ...) {
            
            pkgs <- c("lattice", "latticeExtra", "grid", "gridExtra",
                      "RColorBrewer", "maps")
            tst <- sapply(pkgs, "requireNamespace", 
                          quietly = TRUE, USE.NAMES = FALSE)
            
            if (all(tst == TRUE)) {
              
              try(attachNamespace("gridExtra"), silent = TRUE)
              
              if (is.null(clr)) {
                clr <- colorRampPalette(
                  rev(RColorBrewer::brewer.pal(9, "Spectral")))(1000)
              }
              
              if (missing(y)) y <- 1
              
              p.prm <- paste(pred.prm, "predictor", sep = "_")
              r.prm <- paste(resp.prm, "response", sep = "_")
              
              ps <- slot(x, p.prm)
              rs <- slot(x, r.prm)
              
              if (is.null(ts.vec)) 
                ts.vec <- seq(raster::nlayers(x@resid_response))
              
              xy <- raster::xyFromCell(x@rsq_predictor, 
                                       cell = x@cell_bp)
              
              mode.location.p <- lattice::xyplot(xy[1, 2] ~ xy[1, 1], 
                                                 cex = 2,
                                                 pch = 21, fill = "grey80", 
                                                 col = "black")
              
              if (isTRUE(add.map)) {
                
                try(attachNamespace("maps"), silent = TRUE)
                
                mm180 <- maps::map("world", plot = FALSE, fill = TRUE, 
                                   col = "grey70")
                mm360 <- data.frame(maps::map(plot = FALSE, 
                                              fill = TRUE)[c("x","y")])
                mm360 <- within(mm360, {
                  x <- ifelse(x < 0, x + 360, x)
                  x <- ifelse((x < 1) | (x > 359), NA, x)
                })
                
                if (max(extent(ps)@xmax) > 180) {
                  mm.pred <- mm360
                } else {
                  mm.pred <- mm180
                }
                
                if (max(extent(rs)@xmax) > 180) {
                  mm.resp <- mm360
                } else {
                  mm.resp <- mm180
                }
              } else {
                add.map <- FALSE
                mm.pred <- NULL
                mm.resp <- NULL
              }
              
              
              px.pred <- raster::ncell(ps)
              px.resp <- raster::ncell(rs)
              
              pred.p <- sp::spplot(ps, 
                                   mm = mm.pred, maxpixels = px.pred,
                                   colorkey = list(space = "top",
                                                   width = 0.7, 
                                                   height = 0.8), 
                                   main = paste(p.prm, "mode", x@mode, 
                                                sep = " "),
                                   col.regions = clr, 
                                   panel = function(..., mm) {
                                     lattice::panel.levelplot(...)
                                     if (isTRUE(add.map)) {
                                       lattice::panel.polygon(
                                         mm$x, mm$y, 
                                         lwd = 0.5, 
                                         border = "grey20")
                                     }
                                   }, ...) 
              
              if (show.bp) pred.p <- pred.p + 
                latticeExtra::as.layer(mode.location.p)
              
              resp.p <- sp::spplot(rs, 
                                   mm = mm.resp, maxpixels = px.resp,
                                   colorkey = list(space = "top",
                                                   width = 0.7, 
                                                   height = 0.8), 
                                   main = paste(r.prm, "mode", x@mode, 
                                                sep = " "), 
                                   col.regions = clr, 
                                   panel = function(..., mm) {
                                     lattice::panel.levelplot(...)
                                     if (isTRUE(add.map)) {
                                       lattice::panel.polygon(
                                         mm$x, mm$y, 
                                         lwd = 0.5, 
                                         border = "grey20")
                                     }
                                   }, ...) 
              
              if (show.bp) resp.p <- resp.p + 
                latticeExtra::as.layer(mode.location.p)
              
              md <- x@mode
              
              ts.main <- paste("time series eot", x@mode, "\n",
                               "cumulative explained response domain variance:", 
                               round(x@cum_exp_var * 100 , 2), "%", sep = " ")
              
              eot.ts <- lattice::xyplot(x@eot ~ ts.vec,
                                        type = "b", pch = 20, col = "black", 
                                        ylab = "", xlab = "",
                                        scales = list(tck = c(0.5, 0), 
                                                      x = list(axs = "i")), 
                                        main = ts.main,
                                        panel = lattice::panel.xyplot)
              
              if (anomalies) {
                eot.ts <- eot.ts + 
                  latticeExtra::layer(lattice::panel.abline(h = 0, 
                                                            col = "grey40", 
                                                            lty = 3), 
                                      under = TRUE)
              }
              
              ### set layout to wide or long
              arrange <- arrange[1]
              if (arrange == "wide") ncls <- 2 else ncls <- 1
              
              ### amalgamate pred.p and resp.p according to layout
              c.pred.resp <- gridExtra::arrangeGrob(pred.p, resp.p, 
                                                    ncol = ncls)
              
              ### clear plot area
              grid::grid.newpage()
              
              ### combine c.pred.resp and eot time series and plot
              tst <- gridExtra::arrangeGrob(c.pred.resp, eot.ts, 
                                            heights = c(1, 0.5), 
                                            ncol = 1)
              grid::grid.draw(tst)
              
            } else {
              stop("need packages 'latticeExtra' & 'gridExtra' for plotting EOT results")
            }
            
          }
)

#' @describeIn plot

setMethod('plot', signature(x = 'EotStack',
                            y = 'ANY'), 
          function(x,
                   y,
                   pred.prm = "rsq",
                   resp.prm = "r",
                   show.bp = FALSE,
                   anomalies = TRUE,
                   add.map = TRUE,
                   ts.vec = NULL,
                   arrange = c("wide", "long"),
                   clr = NULL,
                   locations = FALSE,
                   ...) {
            
            pkgs <- c("lattice", "latticeExtra", "grid", "gridExtra",
                      "RColorBrewer", "maps")
            tst <- sapply(pkgs, "requireNamespace", 
                          quietly = TRUE, USE.NAMES = FALSE)
            
            if (all(tst == TRUE)) {
              
              try(attachNamespace("gridExtra"), silent = TRUE)              
              
              if (is.null(clr)) {
                clr <- colorRampPalette(
                  rev(RColorBrewer::brewer.pal(9, "Spectral")))(1000)
              }
              
              if (!locations) {
                
                if (missing(y)) y <- 1 else
                  if (is.character(y)) y <- which(names(x) == y)
                
                p.prm <- paste(pred.prm, "predictor", sep = "_")
                r.prm <- paste(resp.prm, "response", sep = "_")
                
                ps <- slot(x[[y]], p.prm)
                rs <- slot(x[[y]], r.prm)
                
                if (is.null(ts.vec)) 
                  ts.vec <- seq(raster::nlayers(x[[y]]@resid_response))
                
                xy <- raster::xyFromCell(x[[y]]@rsq_predictor, 
                                         cell = x[[y]]@cell_bp)
                
                mode.location.p <- lattice::xyplot(xy[1, 2] ~ xy[1, 1], 
                                                   cex = 2, pch = 21, 
                                                   fill = "grey80", 
                                                   col = "black")
                
                if (isTRUE(add.map)) {
                  
                  try(attachNamespace("maps"), silent = TRUE)
                  
                  mm180 <- maps::map("world", plot = FALSE, fill = TRUE, 
                                     col = "grey70")
                  mm360 <- data.frame(maps::map(plot = FALSE, 
                                                fill = TRUE)[c("x","y")])
                  mm360 <- within(mm360, {
                    x <- ifelse(x < 0, x + 360, x)
                    x <- ifelse((x < 1) | (x > 359), NA, x)
                  })
                  
                  if (max(extent(ps)@xmax) > 180) {
                    mm.pred <- mm360
                  } else {
                    mm.pred <- mm180
                  }
                  
                  if (max(extent(rs)@xmax) > 180) {
                    mm.resp <- mm360
                  } else {
                    mm.resp <- mm180
                  }
                } else {
                  add.map <- FALSE
                  mm.pred <- NULL
                  mm.resp <- NULL
                }
                
                
                px.pred <- raster::ncell(ps)
                px.resp <- raster::ncell(rs)
                
                pred.p <- sp::spplot(ps, 
                                     mm = mm.pred, maxpixels = px.pred,
                                     colorkey = list(space = "top",
                                                     width = 0.7, 
                                                     height = 0.8), 
                                     main = paste(p.prm, "mode", y, 
                                                  sep = " "), 
                                     col.regions = clr, 
                                     panel = function(..., mm) {
                                       lattice::panel.levelplot(...)
                                       if (isTRUE(add.map)) {
                                         lattice::panel.polygon(
                                           mm$x, mm$y, 
                                           lwd = 0.5, 
                                           border = "grey20")
                                       }
                                     }, ...) 
                
                if (show.bp) pred.p <- pred.p + 
                  latticeExtra::as.layer(mode.location.p)
                
                resp.p <- sp::spplot(rs, 
                                     mm = mm.resp, maxpixels = px.resp,
                                     colorkey = list(space = "top",
                                                     width = 0.7, 
                                                     height = 0.8), 
                                     main = paste(r.prm, "mode", y, 
                                                  sep = " "), 
                                     col.regions = clr, 
                                     panel = function(..., mm) {
                                       lattice::panel.levelplot(...)
                                       if (isTRUE(add.map)) {
                                         lattice::panel.polygon(
                                           mm$x, mm$y, 
                                           lwd = 0.5, 
                                           border = "grey20")
                                       }
                                     }, ...) 
                
                if (show.bp) resp.p <- resp.p + 
                  latticeExtra::as.layer(mode.location.p)
                
                md <- x[[y]]@mode
                
                ts.main <- paste("time series eot", x[[y]]@mode, "\n",
                                 "explained response domain variance:", 
                                 round(if (y > 1) {
                                   x[[y]]@cum_exp_var * 100 -
                                     x[[y - 1]]@cum_exp_var * 100
                                 } else {
                                   x[[y]]@cum_exp_var * 100
                                 }, 2), "%", sep = " ")
                
                eot.ts <- lattice::xyplot(x[[y]]@eot ~ ts.vec,
                                          type = "b", pch = 20, 
                                          col = "black", 
                                          ylab = "", xlab = "",
                                          scales = list(tck = c(0.5, 0), 
                                                        x = list(axs = "i")), 
                                          main = ts.main,
                                          panel = lattice::panel.xyplot)
                
                if (anomalies) {
                  eot.ts <- eot.ts + 
                    latticeExtra::layer(lattice::panel.abline(h = 0, 
                                                              col = "grey40", 
                                                              lty = 3), 
                                        under = TRUE)
                }
                
                ### set layout to wide or long
                arrange <- arrange[1]
                if (arrange == "wide") ncls <- 2 else ncls <- 1
                
                ### amalgamate pred.p and resp.p according to layout
                c.pred.resp <- gridExtra::arrangeGrob(pred.p, resp.p, 
                                                      ncol = ncls)
                
                ### clear plot area
                grid::grid.newpage()
                
                ### combine c.pred.resp and eot time series and plot
                tst <- gridExtra::arrangeGrob(c.pred.resp, eot.ts, 
                                              heights = c(1, 0.5), 
                                              ncol = 1)
                grid::grid.draw(tst)
                
              } else {
                plotLocations(x, ...)
              }
              
            } else {
              stop("need packages 'latticeExtra' & 'gridExtra' for plotting EOT results")
            }
            
          }
)


# definde function --------------------------------------------------------
plotLocations <- function(x, ...) {
  
  pkgs <- c("lattice", "latticeExtra", "grid", "gridExtra",
            "RColorBrewer", "maps")
  tst <- sapply(pkgs, "requireNamespace", 
                quietly = TRUE, USE.NAMES = FALSE)
  
  if (all(tst == TRUE)) {
    
    try(attachNamespace("gridExtra"), silent = TRUE)
    try(attachNamespace("maps"), silent = TRUE)
    
    ### plot function
    loc.df <- as.data.frame(do.call("rbind", 
                                    lapply(seq(nmodes(x)), function(i) {
                                      raster::xyFromCell(
                                        x[[i]]@rsq_predictor, 
                                        cell = x[[i]]@cell_bp)
                                    })))
    
    loc.df$eot <- paste("EOT", sprintf("%02.f", seq(x)), 
                        sep = "_")
    
    mm <- maps::map("world", plot = FALSE, fill = TRUE)
    px.pred <- raster::ncell(x[[1]]@r_predictor)
    
    pred.p <- sp::spplot(x[[1]]@rsq_predictor, 
                         mm = mm, maxpixels = px.pred,
                         colorkey = FALSE, 
                         col.regions = "grey50", 
                         panel = function(..., mm) {
                           lattice::panel.levelplot(...)
                           lattice::panel.polygon(
                             mm$x, mm$y, lwd = 0.5, 
                             border = "grey20", 
                             col = "grey70")
                         }, ...) 
    
    clrs.hcl <- function(n) {
      hcl(h = seq(270, 0, length.out = n), 
          c = 70, l = 50, fixup = TRUE)
    }
    
    n <- nmodes(x)
    clrs <- clrs.hcl(n)
    
    points.p <- lattice::xyplot(y ~ x, data = loc.df, col = "black", 
                                fill = clrs, pch = 21,
                                cex = 2)
    
    out <- pred.p + latticeExtra::as.layer(points.p)
    
    grid::grid.newpage()
    
    map.vp <- grid::viewport(x = 0, y = 0, 
                             height = 1, width = 0.85,
                             just = c("left", "bottom"))
    
    grid::pushViewport(map.vp)
    
    print(out, newpage = FALSE)
    
    grid::downViewport(lattice::trellis.vpname(name = "figure"))
    
    leg.vp <- grid::viewport(x = 1, y = 0.5, 
                             height = n / 10, width = 0.15,
                             just = c("left", "centre"))
    
    grid::pushViewport(leg.vp)  
    
    if(n == 1) ypos <- 0.5 else ypos <- seq(0.95, 0.05, length.out = n + 2)
    if(n == 1) ypos <- ypos else ypos <- ypos[-c(1, length(ypos))]
    xpos.pts <- grid::unit(0.15, "npc")
    size.pts <- 1 / n
    
    for (i in 1:n) {
      
      vp <- grid::viewport(x = xpos.pts, y = ypos[i], 
                           height = size.pts, width = 0.1,
                           just = c("left", "centre"))
      
      grid::pushViewport(vp)
      
      grid::grid.circle(gp = grid::gpar(fill = clrs[i], 
                                        col = "black"))
      
      grid::upViewport()
      
    }
    
    xpos.txt <- grid::unit(0.25, "npc")
    width.txt <- 0.7
    
    for (i in 1:n) {
      
      vp <- grid::viewport(x = xpos.txt, y = ypos[i], 
                           height = size.pts, width = width.txt,
                           just = c("left", "centre"))
      
      grid::pushViewport(vp)
      
      txt <- grid::textGrob(x = 0.2, sort(names(x))[i],
                            just = "left")
      
      grid::grid.draw(txt)
      
      grid::popViewport()
      
    }
    
    grid::upViewport(0)
    
  } else {    
    stop("need packages 'gridExtra', 'latticeExtra' & 'maps' to plot locations")
  }
  
} 