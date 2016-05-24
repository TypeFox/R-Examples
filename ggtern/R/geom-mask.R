#' Apply Manual Clipping Mask
#' 
#' This function creates a manual clipping mask, which in turn suppresses the standard clipping mask that would otherwise
#' be rendered in the foregound rendering procedure, giving the user control over the exact placement with respect to 
#' other layers. For example, the user may wish to have the clipping mask placed after 
#' the \code{geom_point(...)} layer, but before the \code{geom_label(...)} layer, this situation has been
#' demonstrated in the example below. In the event that the user wishes to suppress the mask altogether, then a convenience
#' function has been provided, \code{theme_nomask()}.
#' @author Nicholas Hamilton
#' @examples 
#' data(Feldspar)
#' x = ggtern(Feldspar,aes(Ab,An,Or,label=Experiment)) + geom_point()
#' 
#' #Default Behaviour
#' x + geom_label()
#' 
#' #Insert manual mask before the labels, to prevent them being truncated
#' x + geom_point(size=6) + geom_mask() + geom_label()
#' @export
#' @rdname geom_mask
geom_mask <- function() {
  layer(
    data        = data.frame(x=1,y=1,z=1),
    mapping     = NULL,
    stat        = "identity",
    geom        = GeomMask,
    position    = "identity",
    show.legend = FALSE,
    inherit.aes = FALSE,
    params      = list(
      na.rm     = TRUE
    )
  )
}


#' @rdname geom_mask
#' @format NULL
#' @usage NULL
GeomMask <- ggproto("GeomMask", Geom,
  default_aes = aes("x","y","z"),
  draw_panel  = function(self, data, panel_scales, coord){
    items = gList()
    if(!inherits(coord,'CoordTern'))return(items)
    tryCatch({
      themeElements = c('tern.plot.background','tern.panel.background')
      for(ixEl in c(1:2)){
        e  = calc_element(themeElements[ixEl],coord$theme %||% theme_get(),verbose=F)
        if(!identical(e,element_blank())){
          a  = c(0,1); b = 0.5; if(ixEl == 2){ a = expand_range(a,0.01) }; #EXPAND THE TOP MASK SLIGHTLY
          if(ixEl == 1){
            ex  = data.frame(diag(1,3,3)); colnames(ex) = as.character(coord$mapping)
          }else{
            ex  = .get.tern.extremes(coord,panel_scales,transform=FALSE)
          }
          ex  = coord$transform(ex,scale_details = panel_scales)
          ex = rbind(ex,ex[1,,drop=F])
          for(ix in c(1:2)){
            xvals = c(a[1],a[1],b[1],if(ix==1){ ex$x }else{NULL},b[1],a[2],a[2],a[1])
            yvals = c(a[1],a[2],a[2],if(ix==1){ ex$y }else{NULL},a[2],a[2],a[1],a[1])
            grob     <- polygonGrob(  x = xvals,
                                      y = yvals,
                                      default.units = "npc",
                                      id   = rep(1,length(xvals)),
                                      gp   = gpar(  col  = if(ix==1       | is.null(e$colour)){NA}else{ e$colour },
                                                    fill = alpha(if(ix==2 | is.null(e$fill)){NA}else{e$fill},is.numericor(e$alpha,1)),
                                                    lwd  = if(ix==1){0}else{ifthenelse(!is.numeric(e$size),0,e$size)*find_global_tern(".pt")},
                                                    lty  = e$linetype)
            )
            
            vp <- viewport(x = 0.5, y = 0.5, width = 1, height = 1, just = c("center","center"),clip=if(ixEl == 1) 'inherit' else 'off')
            grob = editGrob(grob, vp = vp, name = sprintf("mask-%i-%i",ixEl,ix))
            items[[length(items) + 1]] = grob
          }
        }
      }
    },error=function(e){writeLines(as.character(e))})
    items
  },
  draw_key = FALSE
)