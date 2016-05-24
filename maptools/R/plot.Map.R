# Copyright 2000-2001 (c) Nicholas Lewin-Koh 
# modifications 2001-2004 Roger Bivand

plot.Map <- function(x, recs, auxvar=NULL, add=FALSE, fg ='gray', 
                   ol='black', prbg=NULL, glyph=16, color='red',
                   type='q', nclass=5, ...) 
{
  .Deprecated("plot.Spatial", package="maptools")
  theMap <- x
  if(!inherits(theMap, "Map"))
  stop("Map.obj must be of class Map")

  if(missing(recs)) recs <- 1:attr(theMap$Shapes,'nshps')

  if (length(fg) != length(recs)) fg <- rep(fg[1], length(recs))

  xylims <- Map2maplim(theMap)

  if(!add){
     plot(xylims$x, xylims$y, asp=1, type='n',...)
  }
  if(!is.null(prbg)) {
    plim <- par()$usr
    rect(plim[1], plim[2], plim[3], plim[4], col=prbg) #,border=par()$bg)
   }

  ret <- NULL
  if(attr(theMap$Shapes,'shp.type') == 'point' ) {
    for(i in 1:length(recs)) {
      points(theMap$Shapes[[recs[i]]]$verts, pch=glyph, col=fg[i])
    }
  }
  if(attr(theMap$Shapes,'shp.type') == 'arc'){
    for(i in 1:length(recs)) {
      if(attr(theMap$Shapes[[recs[i]]], 'nParts') == 1) {
        lines(theMap$Shapes[[recs[i]]]$verts, col=ol)
      }
      if(attr(theMap$Shapes[[recs[i]]], 'nParts') > 1){
        for(j in 1:attr(theMap$Shapes[[recs[i]]], 'nParts')) {
	  if(j < attr(theMap$Shapes[[recs[i]]], 'nParts'))
             lines(theMap$Shapes[[recs[i]]]$verts[j:(j+1)-1], col= ol)
          else
             lines(theMap$Shapes[[recs[i]]]$verts[j:attr(theMap$Shapes[[recs[i]]], 'nVerts')], col= ol)  
        }
      }
    }
  }
  if(attr(theMap$Shapes,'shp.type') == 'poly'){
    if(!is.null(auxvar) && nclass > 1) {
      if (length(auxvar) != attr(theMap$Shapes,'nshps'))
        stop("lengths conflict")
      col.rmp <- color.ramp(nclass, color=color, nvec=auxvar[recs], type=type)
      for(i in 1:length(recs)) {
        ii <- recs[i]
        if(attr(theMap$Shapes[[ii]], 'nParts') == 1) {
          polygon(theMap$Shapes[[ii]]$verts,
                col=col.rmp$ramp[col.rmp$col.class[i]],
                border= ol, ...)
        }
        if(attr(theMap$Shapes[[ii]], 'nParts') > 1) {
          for(j in 1:attr(theMap$Shapes[[ii]], 'nParts')) {
	    if(j < attr(theMap$Shapes[[ii]], 'nParts')) {
              polygon(theMap$Shapes[[ii]]$verts[(theMap$Shapes[[ii]]$Pstart[j]+1):
                    theMap$Shapes[[ii]]$Pstart[j+1],],
		    col=col.rmp$ramp[col.rmp$col.class[i]], border=ol, ...)
            } else {
              polygon(theMap$Shapes[[ii]]$verts[(theMap$Shapes[[ii]]$Pstart[j]+1):
                    attr(theMap$Shapes[[ii]],'nVerts'),],
                    col=col.rmp$ramp[col.rmp$col.class[i]],border= ol, ...)
            }
          }
        }
      }
      ret <- col.rmp
    } else {
      for(i in 1:length(recs)) {
        ii <- recs[i]
        if(attr(theMap$Shapes[[ii]],'nParts') == 1) {
          polygon(theMap$Shapes[[ii]]$verts, col=fg[i], border= ol, ...)
        }
        if(attr(theMap$Shapes[[ii]],'nParts') > 1) {
          for(j in 1:attr(theMap$Shapes[[ii]], 'nParts')) {
	    if(j<attr(theMap$Shapes[[ii]], 'nParts')) {
              polygon(theMap$Shapes[[ii]]$verts[(theMap$Shapes[[ii]]$Pstart[j]+1):
                    theMap$Shapes[[ii]]$Pstart[j+1],], col=fg[i] ,border= ol, ...)
            } else {
              polygon(theMap$Shapes[[ii]]$verts[(theMap$Shapes[[ii]]$Pstart[j]+1):
                    attr(theMap$Shapes[[ii]],'nVerts'),],
                    col=fg[i], border=ol, ...)
            }
          }
        }
      }
    }
  }
  if(attr(theMap$Shapes,'shp.type')=='multipoint'){
    stop("Multipoint shape type not yet plotted")
  }
  invisible(ret)
}


