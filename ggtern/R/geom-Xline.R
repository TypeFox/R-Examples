#' Fixed Value Lines
#' 
#' Plot fixed value lines, for the top, left and right axis, analagous to the \code{\link{geom_hline}}
#' and \code{\link{geom_vline}} geometries in \code{\link[=ggplot]{ggplot2}}
#' 
#' @name geom_Xline
#' @aliases geom_Tline geom_Lline geom_Rline Tline tline Lline lline Rline rline
#' @author Nicholas Hamilton
#' @inheritParams ggplot2:::geom_hline
#' @param Tintercept,Lintercept,Rintercept the intercepts for the T, L and R axis respectively
#' @examples 
#' ggtern() + 
#' geom_Tline(Tintercept=.5,arrow=arrow(), colour='red') + 
#' geom_Lline(Lintercept=.2, colour='green') + 
#' geom_Rline(Rintercept=.1, colour='blue')
#' @rdname geom_Xline
NULL

#' @export
#' @name geom_Xline
#' @rdname geom_Xline
geom_Tline <- function(mapping = NULL, data = NULL,
                       ...,
                       Tintercept,
                       na.rm = FALSE, show.legend = NA) {
  
  # Act like an annotation
  if (!missing(Tintercept)) {
    data <- data.frame(Tintercept = Tintercept)
    mapping <- aes(Tintercept = Tintercept)
    show.legend=FALSE
  }
  
  layer(
    data        = data,
    mapping     = mapping,
    stat        = StatIdentity,
    geom        = GeomTline,
    position    = PositionIdentity,
    show.legend = show.legend,
    inherit.aes = FALSE,
    params      = list(
      na.rm     = na.rm,
      ...
    )
  )
}

#' @rdname  geom_Xline
#' @format NULL
#' @usage NULL
#' @export
GeomTline <- ggproto("GeomTline",Geom,
                     draw_panel = function(self,data, panel_scales,coord){
                       .drawTRLLinesX(self,data,panel_scales,coord,'T')
                     },
                     default_aes = aes(colour = "black", size = 0.5, linetype = 1, alpha = NA,arrow=NULL),
                     required_aes = c("Tintercept"),
                     draw_key = draw_key_Tline
)

#' @export
#' @rdname geom_Xline
Tline <- geom_Tline

#' @export
#' @name geom_Xline
tline <- Tline

#' @export
#' @rdname geom_Xline
geom_Lline <- function(mapping = NULL, data = NULL,
                       ...,
                       Lintercept,
                       na.rm = FALSE, show.legend = NA) {
  
  # Act like an annotation
  if (!missing(Lintercept)) {
    data <- data.frame(Lintercept = Lintercept)
    mapping <- aes(Lintercept = Lintercept)
    show.legend=FALSE
  }
  
  layer(
    data        = data,
    mapping     = mapping,
    stat        = StatIdentity,
    geom        = GeomLline,
    position    = PositionIdentity,
    show.legend = show.legend,
    inherit.aes = FALSE,
    params      = list(
      na.rm     = na.rm,
      ...
    )
  )
}

#' @export
#' @rdname geom_Xline
#' @format NULL
#' @usage NULL
GeomLline <- ggproto("GeomLline",Geom,
                     draw_panel = function(self,data, panel_scales,coord){
                       .drawTRLLinesX(self,data,panel_scales,coord,'L')
                     },
                     default_aes = aes(colour = "black", size = 0.5, linetype = 1, alpha = NA,arrow=NULL),
                     required_aes = c("Lintercept"),
                     draw_key = draw_key_Lline
)

#' @export
#' @rdname geom_Xline
Lline <- geom_Lline

#' @export
#' @name geom_Xline
lline <- Lline

#' @export
#' @rdname geom_Xline
geom_Rline <- function(mapping = NULL, data = NULL,
                       ...,
                       Rintercept,
                       na.rm = FALSE, show.legend = NA) {
  
  # Act like an annotation
  if (!missing(Rintercept)) {
    data <- data.frame(Rintercept = Rintercept)
    mapping <- aes(Rintercept = Rintercept)
    show.legend=FALSE
  }
  
  layer(
    data        = data,
    mapping     = mapping,
    stat        = StatIdentity,
    geom        = GeomRline,
    position    = PositionIdentity,
    show.legend = show.legend,
    inherit.aes = FALSE,
    params      = list(
      na.rm     = na.rm,
      ...
    )
  )
}

#' @export
#' @rdname geom_Xline
#' @format NULL
#' @usage NULL
GeomRline <- ggproto("GeomRline",Geom,
                     draw_panel = function(self,data, panel_scales,coord){
                       .drawTRLLinesX(self,data,panel_scales,coord,'R')
                     },
                     default_aes = aes(colour = "black", size = 0.5, linetype = 1, alpha = NA,arrow=NULL),
                     required_aes = c("Rintercept"),
                     draw_key = draw_key_Rline
)

#' @export
#' @rdname geom_Xline
Rline <- geom_Rline

#' @export
#' @name geom_Xline
rline <- Rline




#internal function
.drawTRLLinesX <- function(self,data,panel_scales, coord, feat){
  if(!'CoordTern' %in% class(coord)) return(zeroGrob())
  axisNames = names(coord$mapping)
  if(!feat %in% axisNames) stop(sprintf("Invalid 'feat' variable ('%s'), please use %s",
                                        feat,
                                        joinCharacterSeries(axisNames,'or')),call.=FALSE)
  data      = remove_missing(data,vars=paste(axisNames,'intercept',sep=""),na.rm=TRUE,name=class(self)[1],finite=TRUE)
  if(empty(data)) return(zeroGrob())
  
  ranges    = coord$range(panel_scales);
  mapping   = coord$mapping
  
  #Get the correct sequence of other axes, relative to the featured axis
  getOthers = function(mapping,feat){
    others  = rep(names(mapping),2)
    ix.feat = which(others == feat)
    mapping[ (others[-ix.feat])[ ix.feat[1]+c(0,1) ] ]
  }; others = getOthers(mapping,feat)
  
  featIntercept = sprintf("%sintercept",feat)
  for(x in c(0:1) ){
    s = if(x == 0) "" else "end"
    data[,sprintf("%s%s",mapping[[feat]],s) ] = data[,featIntercept]
    data[,sprintf("%s%s",  others[[1+x]],s) ] = 1 - data[, mapping[[feat]] ] - min(coord$scales[[ names(others)[2-x] ]]$limits)
    data[,sprintf("%s%s",  others[[2-x]],s) ] = min(coord$scales[[ names(others)[2-x] ]]$limits)
  }
  scale_details = list( x.range = ranges[['x']], y.range = ranges[['y']] )
  data          = coord$transform(data,scale_details)
  grob          = zeroGrob()
  tryCatch({
    cw   = calc_element('tern.axis.clockwise',coord$theme) ##Clockwise
    grob = segmentsGrob(if(!cw) data$x else data$xend, 
                        if(!cw) data$y else data$yend, 
                        if( cw) data$x else data$xend, 
                        if( cw) data$y else data$yend,
                        default.units     = "npc",
                        gp = gpar(col     = alpha(data$colour, data$alpha),
                                  fill    = alpha(data$colour, data$alpha),
                                  lwd     = data$size*find_global_tern(".pt"),
                                  lty     = data$linetype,
                                  lineend = 'butt'
                        ),
                        arrow = data$arrow)
  },error=function(e){message(as.character(e))})
  grob
}
