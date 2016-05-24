## #old code for wrapping

## sectionview <- function(model, ...) {
##   if (class(model)[1]=="km") {
##     sectionview.km(model=model,...)
##   } else if (class(model)[1]=="list") {
##     sectionview.dm(model=model,...)
##   }
## }


## sectionview3d <- function(model, ...) {
##   if (class(model)[1]=="km") {
##     sectionview3d.km(model=model,...)
##   } else if (class(model)[1]=="list") {
##     sectionview3d.dm(model=model,...)
##   }
## }

## Wrapper for view
if(!isGeneric("view")) {
    setGeneric(name = "view",
               def = function(type="auto", model, ...) standardGeneric("view")
               )
}

setMethod("view", c("character","km"), 
          function(type = "auto",
                   model, model.type = "UK",
                   center = NULL, axis = NULL,
                   npoints = NULL,
                   col_points = "red",
                   col_surf = "blue",
                   col_needles = NA,
                   conf_lev = NULL,
                   conf_blend = NULL,
                   bg_blend = NULL,
                   mfrow = NULL,
                   nlevels = 10,
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = NULL,
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              if(is.null(type)) type="auto"
              if(type=="auto"){
                  if(model@d==1) type="xy"
                  if(model@d==2) type="3d"
                  if(model@d>2) {
                      type="section"
                      if(is.null(center)) center=colMeans(model@X)
                  }
              }
              
              if (type == "section" || type == "xy")
                  sectionview.km(model = model, type = model.type,
                                 center = center, axis = axis,
                                 npoints = ifelse(is.null(npoints),100,npoints),
                                 col_points = col_points,
                                 col_surf = col_surf,
                                 conf_lev = ifelse(is.null(conf_lev),c(0.5,0.8,0.9,0.95,0.99),conf_lev),
                                 conf_blend = conf_blend,
                                 bg_blend = ifelse(is.null(bg_blend),5,bg_blend),
                                 mfrow = mfrow,
                                 Xname = Xname,
                                 yname=yname,
                                 Xscale = Xscale,
                                 yscale = yscale,
                                 xlim = xlim,
                                 ylim = ylim,
                                 title = title,
                                 ...)
              
              if (type == "section3d" || type == "3d")
                  sectionview3d.km(model = model, type = model.type,
                                   center = center, axis = axis,
                                   npoints = ifelse(is.null(npoints),20,npoints),
                                   col_points = col_points,
                                   col_surf = col_surf,
                                   col_needles = col_needles,
                                   conf_lev = ifelse(is.null(conf_lev),c(0.95),conf_lev),
                                   conf_blend = conf_blend,
                                   bg_blend = ifelse(is.null(bg_blend),5,bg_blend),
                                   Xname = Xname,
                                   yname = yname,
                                   Xscale = Xscale,
                                   yscale = yscale,
                                   xlim = xlim,
                                   ylim = ylim,
                                   title = title,
                                   ...)
              
              if (type == "contour")
                  contourview.km(model = model, type = model.type,
                                 center = center, axis = axis,
                                 npoints = ifelse(is.null(npoints),20,npoints),
                                 col_points = col_points,
                                 col_surf = col_surf,
                                 bg_blend = ifelse(is.null(bg_blend),1,bg_blend),
                                 mfrow = mfrow,
                                 nlevels = nlevels,
                                 Xname = Xname,
                                 yname = yname,
                                 Xscale = Xscale,
                                 yscale = yscale,
                                 xlim = xlim,
                                 ylim = ylim,
                                 title = title,
                                 ...)
          }
          )

setMethod("view", c("character","list"), 
          function(type = "auto",
                   model,
                   center = NULL, axis = NULL,
                   npoints = NULL,
                   col_points = "red",
                   col_surf = "blue",
                   col_needles = NA,
                   bg_blend = NULL,
                   mfrow = NULL,
                   nlevels = 10,
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = NULL,
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              if(is.null(type)) type="auto"
              if(type=="auto"){
                  if(ncol(model$data$X)==1) type="xy"
                  if(ncol(model$data$X)==2) type="3d"
                  if(ncol(model$data$X)>2) {
                      type="section"
                      if(is.null(center)) center=colMeans(model$data$X)
                  }
              }
              
              if (type == "section" || type == "xy")
                  sectionview.list(model = model,
                                   center = center, axis = axis,
                                   npoints = ifelse(is.null(npoints),100,npoints),
                                   col_points = col_points,
                                   col_surf = col_surf,
                                   bg_blend = ifelse(is.null(bg_blend),5,bg_blend),
                                   mfrow = mfrow,
                                   Xname = Xname,
                                   yname=yname,
                                   Xscale = Xscale,
                                   yscale = yscale,
                                   xlim = xlim,
                                   ylim = ylim,
                                   title = title,
                                   ...)
              
              if (type == "section3d" || type == "3d")
                  sectionview3d.list(model = model,
                                     center = center, axis = axis,
                                     npoints = ifelse(is.null(npoints),20,npoints),
                                     col_points = col_points,
                                     col_surf = col_surf,
                                     col_needles = col_needles,
                                     bg_blend = ifelse(is.null(bg_blend),5,bg_blend),
                                     Xname = Xname,
                                     yname = yname,
                                     Xscale = Xscale,
                                     yscale = yscale,
                                     xlim = xlim,
                                     ylim = ylim,
                                     title = title,
                                     ...)
              
              if (type == "contour")
                  contourview.list(model = model,
                                   center = center, axis = axis,
                                   npoints = ifelse(is.null(npoints),20,npoints),
                                   col_points = col_points,
                                   col_surf = col_surf,
                                   bg_blend = ifelse(is.null(bg_blend),1,bg_blend),
                                   mfrow = mfrow,
                                   nlevels = nlevels,
                                   Xname = Xname,
                                   yname = yname,
                                   Xscale = Xscale,
                                   yscale = yscale,
                                   xlim = xlim,
                                   ylim = ylim,
                                   title = title,
                                   ...)
          }
          )

setMethod("view", c("character","function"), 
          function(type = "auto",
                   model , dim,
                   center = NULL, axis = NULL,
                   npoints = NULL,
                   col = "blue",
                   mfrow = NULL,
                   nlevels = 10,
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = c(0,1),
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              if(is.null(type)) type="auto"
              if(type=="auto"){
                  if(dim==1) type="xy"
                  if(dim==2) type="3d"
                  if(dim>2) {
                      type="section"
                      if(is.null(center)) if(!is.null(xlim)) center=colMeans(xlim) else center=rep(0.5,dim)
                  }
              }
              
              if (type == "section" || type == "xy")
                  sectionview.fun(fun = model, dim = dim,
                                  center = center, axis = axis,
                                  npoints = ifelse(is.null(npoints),100,npoints),
                                  col = col,
                                  mfrow = mfrow,
                                  Xname = Xname,
                                  yname=yname,
                                  Xscale = Xscale,
                                  yscale = yscale,
                                  xlim = xlim,
                                  ylim = ylim,
                                  title = title,
                                  ...)
              
              if (type == "section3d" || type == "3d")
                  sectionview3d.fun(fun = model, dim = dim,
                                    center = center, axis = axis,
                                    npoints = ifelse(is.null(npoints),20,npoints),
                                    col = col,
                                    Xname = Xname,
                                    yname = yname,
                                    Xscale = Xscale,
                                    yscale = yscale,
                                    xlim = xlim,
                                    ylim = ylim,
                                    title = title,
                                    ...)
              
              if (type == "contour")
                  contourview.fun(fun = model, dim = dim,
                                  center = center, axis = axis,
                                  npoints = ifelse(is.null(npoints),20,npoints),
                                  col = col,
                                  mfrow = mfrow,
                                  nlevels = nlevels,
                                  Xname = Xname,
                                  yname = yname,
                                  Xscale = Xscale,
                                  yscale = yscale,
                                  xlim = xlim,
                                  ylim = ylim,
                                  title = title,
                                  ...)
          }
          )


## Wrapper for sectionview
if(!isGeneric("sectionview")) {
    setGeneric(name = "sectionview",
               def = function(model, ...) standardGeneric("sectionview")
               )
}

setMethod("sectionview", "km", 
          function(model,
                   type = "UK",
                   center = NULL,
                   npoints = 100,
                   col_points = "red",
                   col_surf = "blue",
                   conf_lev = c(0.5,0.8,0.9,0.95,0.99),
                   conf_blend = NULL,
                   bg_blend = 5,
                   mfrow = NULL,
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = NULL,
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              sectionview.km(model = model,
                             type = type,
                             center = center,
                             npoints = npoints,
                             col_points = col_points,
                             col_surf = col_surf,
                             conf_lev = conf_lev,
                             conf_blend = conf_blend,
                             bg_blend = bg_blend,
                             mfrow = mfrow,
                             Xname = Xname,
                             yname=yname,
                             Xscale = Xscale,
                             yscale = yscale,
                             xlim = xlim,
                             ylim = ylim,
                             title = title,
                             ...)		
          }
          )

setMethod("sectionview", "list", 
          function(model,
                   center = NULL,
                   npoints = 100,
                   col_points = "red",
                   col_surf = "blue",
                   bg_blend = 5,
                   mfrow = NULL,
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = NULL,
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              sectionview.list(model = model,
                               center = center,
                               npoints = npoints,
                               col_points = col_points,
                               col_surf = col_surf,
                               bg_blend = bg_blend,
                               mfrow = mfrow,
                               Xname = Xname,
                               yname = yname,
                               Xscale = Xscale,
                               yscale = yscale,
                               xlim = xlim,
                               ylim = ylim,
                               title = title,
                               ...)		
          }
          )

setMethod("sectionview", "function", 
          function(model,dim,
                   center = NULL, axis = NULL,
                   npoints = 100,
                   col = "blue",
                   mfrow = NULL,
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = c(0,1),
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              sectionview.fun(fun = model,
                              dim = dim,
                              center = center, axis = axis,
                              npoints = npoints,
                              col = col,
                              mfrow = mfrow,
                              Xname = Xname,
                              yname = yname,
                              Xscale = Xscale,
                              yscale = yscale,
                              xlim = xlim,
                              ylim = ylim,
                              title = title,
                              ...)        
          }
          )

#Wrapper for sectionview3d
if(!isGeneric("sectionview3d")) {
    setGeneric(name = "sectionview3d",
               def = function(model, ...) standardGeneric("sectionview3d")
               )
}

setMethod("sectionview3d", "km", 
          function(model, type = "UK",
                   center = NULL, axis = NULL,
                   npoints = 20,
                   col_points = "red",
                   col_surf = "blue",
                   col_needles = NA,
                   conf_lev = c(0.95),
                   conf_blend = NULL,
                   bg_blend = 5,
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = NULL,
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              sectionview3d.km(model = model, type = type,
                               center = center, axis = axis,
                               npoints = npoints,
                               col_points = col_points,
                               col_surf = col_surf,
                               col_needles = col_needles ,
                               conf_lev = conf_lev,
                               conf_blend = conf_blend,
                               bg_blend = bg_blend,
                               Xname = Xname,
                               yname = yname,
                               Xscale = Xscale,
                               yscale = yscale,
                               xlim = xlim,
                               ylim = ylim,
                               title = title,
                               ...)		
          }
          )

setMethod("sectionview3d", "list", 
          function(model,
                   center = NULL, axis = NULL,
                   npoints = 20,
                   col_points = "red",
                   col_surf = "blue",
                   bg_blend = 5,
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = NULL,
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              sectionview3d.list(model = model,
                                 center = center, axis = axis,
                                 npoints = npoints,
                                 col_points = col_points,
                                 col_surf = col_surf,
                                 bg_blend = bg_blend,
                                 Xname = Xname,
                                 yname = yname,
                                 Xscale = Xscale,
                                 yscale = yscale,
                                 xlim = xlim,
                                 ylim = ylim,
                                 title = title,
                                 ...)		
          }
          )

setMethod("sectionview3d", "function", 
          function(model,dim,
                   center = NULL, axis = NULL,
                   npoints = 20,
                   col = "blue",
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = c(0,1),
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              sectionview3d.fun(fun = model, dim = dim,
                                center = center, axis = axis,
                                npoints = npoints,
                                col = col,
                                Xname = Xname,
                                yname = yname,
                                Xscale = Xscale,
                                yscale = yscale,
                                xlim = xlim,
                                ylim = ylim,
                                title = title,
                                ...)        
          }
          )

#Wrapper for contourview
if(!isGeneric("contourview")) {
    setGeneric(name = "contourview",
               def = function(model, ...) standardGeneric("contourview")
               )
}

setMethod("contourview", "km", 
          function(model, type = "UK",
                   center = NULL, axis = NULL,
                   npoints = 20,
                   col_points = "red",
                   col_surf = "blue",
                   bg_blend = 1,
                   nlevels = 10,
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = NULL,
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              contourview.km(model = model, type = type,
                             center = center, axis = axis,
                             npoints = npoints,
                             col_points = col_points,
                             col_surf = col_surf,
                             bg_blend = bg_blend,
                             nlevels = nlevels,
                             Xname = Xname,
                             yname = yname,
                             Xscale = Xscale,
                             yscale = yscale,
                             xlim = xlim,
                             ylim = ylim,
                             title = title,
                             ...)        
          }
          )

setMethod("contourview", "list", 
          function(model,
                   center = NULL, axis = NULL,
                   npoints = 20,
                   col_points = "red",
                   col_surf = "blue",
                   bg_blend = 1,
                   nlevels = 10,
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = NULL,
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              contourview.list(model = model,
                               center = center, axis = axis,
                               npoints = npoints,
                               col_points = col_points,
                               col_surf = col_surf,
                               bg_blend = bg_blend,
                               nlevels = nlevels,
                               Xname = Xname,
                               yname = yname,
                               Xscale = Xscale,
                               yscale = yscale,
                               xlim = xlim,
                               ylim = ylim,
                               title = title,
                               ...)        
          }
          )

setMethod("contourview", "function", 
          function(model,dim,
                   center = NULL, axis = NULL,
                   npoints = 20,
                   col = "blue",
                   nlevels = 10,
                   Xname = NULL,
                   yname = NULL,
                   Xscale = 1,
                   yscale = 1,
                   xlim = c(0,1),
                   ylim = NULL,
                   title = NULL,
                   ...){
              
              contourview.fun(fun = model, dim = dim,
                              center = center, axis = axis,
                              npoints = npoints,
                              col = col,
                              nlevels = nlevels,
                              Xname = Xname,
                              yname = yname,
                              Xscale = Xscale,
                              yscale = yscale,
                              xlim = xlim,
                              ylim = ylim,
                              title = title,
                              ...)        
          }
          )
