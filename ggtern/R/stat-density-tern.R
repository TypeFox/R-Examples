#'Density Estimate
#' 
#' Perform a 2D kernel density estimatation using kde2d and display the results with contours. This can be 
#' useful for dealing with overplotting. Additional weight aesthetic (see aesthetic section below) permits better weighting if desired
#' 
#'@section Aesthetics:
#' \Sexpr[results=rd,stage=build]{ggtern:::rd_aesthetics("stat", "density_tern")}
#' @inheritParams ggplot2::stat_density2d
#' @name stat_density_tern
#' @rdname stat_density_tern
#' @param base the base transformation of the data, options include 'identity' (ie direct on the cartesian space), or 'ilr'
#' which means to use the isometric log ratio transformation.
#' @param h Bandwidth (vector of length two) as a multiple of the best estimate, estimated using \code{\link[MASS]{bandwidth.nrd}}. 
#' @param weight weighting for weighted kde2d esimate, default's to 1, which is non-weighted and equivalent to the usual kde2d calculation
#' @param expand Calculate on a mesh which extends beyond the grid of the plot region by this amount
#' If \code{NULL}, estimated using \code{\link[MASS]{bandwidth.nrd}}.
#' @examples 
#' #Plot Density Estimate, on isometric log ratio transformation of original data
#' data(Feldspar)
#' ggtern(Feldspar,aes(Ab,An,Or)) + 
#' stat_density_tern(aes(fill=..level..),geom='polygon',base='ilr')
#' @export
stat_density_tern <- function(mapping = NULL, data = NULL, geom = "density_tern",position = "identity",
                              ...,
                              contour = TRUE,
                              n = 100, h = NULL, na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE,weight=1,base='ilr',expand = c(.5,.5)) {
  layer(
    data        = data,
    mapping     = mapping,
    stat        = StatDensityTern,
    geom        = geom,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(
      na.rm     = na.rm,
      contour   = contour,
      n         = n,
      h         = h,
      base      = base,
      expand    = expand,
      weight    = weight,
      ...
    )
  )
}

#' @rdname stat_density_tern
#' @format NULL
#' @usage NULL
#' @export
StatDensityTern  = ggproto("StatDensityTern", 
  Stat,
  retransform    = FALSE,
  required_aes   = c("x",'y',"z"),
  setup_data = function(data, params) {
    data$weight = data$weight %||% params$weight %||% numeric(nrow(data)) + 1
    data
  },
  compute_group  = function(self,data, scales, na.rm = FALSE, 
                            h = NULL, contour = TRUE, n = 100, base='ilr', expand = c(.5,.5), weight=NULL) {
    
    if(!base %in% c('identity','ilr')) stop('base must be either identity or ilr',call.=FALSE)
    data   = remove_missing(data,vars=self$required_aes,na.rm=TRUE,name='StatConfidence',finite=TRUE)
    if(empty(data)) return(data.frame())
    
    coord = coord_tern()
    raes  = self$required_aes
    
    f    = get(base,mode='function')
    fInv = get(sprintf("%sInv",base),mode='function')

    if(base == 'identity') data = tlr2xy(data,coord,inverse=FALSE,scale=TRUE) 
    data[,c('x','y')] = f( as.matrix(data[,which(colnames(data) %in% raes)]) )

    h.est = c(MASS::bandwidth.nrd(data$x), MASS::bandwidth.nrd(data$y)) 
    h = ifthenelse(is.null(h),h.est, h.est*if(length(h)!=2){ rep(h[1],2) } else{ h })
    
    .getLims = function(){
      if(base == 'identity'){
        rngx = coord$limits$x; rngy = coord$limits$y
      }else{
        rngx = range(data$x);  rngy = range(data$y)
      }
      expand = ifthenelse(length(expand) != 2,rep(expand[1],2),expand)
      c(expand_range(rngx,expand[1]),expand_range(rngy,expand[2]))
    }
    
    dens <- kde2d.weighted(data$x, data$y, h=h, n=n, lims=.getLims(), w=data$weight)
    data <- data.frame(expand.grid(x = dens$x, y = dens$y), z = as.vector(dens$z),group=data$group[1])
    if (contour) {
      data = StatContour$compute_panel(data, scales)
    } else {
      names(data) <- c("x", "y", "density", "group")
      data$level <- 1
      data$piece <- 1
    }
    if(base == 'identity') data = tlr2xy(data,coord,inverse=TRUE,scale=TRUE)
    data[,self$required_aes] = fInv( as.matrix(data[,which(colnames(data) %in% raes)]) )
    data
  }
)
