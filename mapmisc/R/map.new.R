map.new = function(x, legendRight=FALSE, buffer=0, mar=c(0,0,0,0), ...) {
	
	if(!is.null(attributes(x)$ellipse)) {
		ellipse = x = attributes(x)$ellipse
	} else {
		ellipse = NULL
	}
	
  xpoints = as(extend(extent(x), buffer), 'SpatialPoints')

	thecrs = try(proj4string(x), silent=TRUE) 
	if(class(thecrs)!="try-error" & !is.na(thecrs)) {
			proj4string(xpoints) = CRS(thecrs)
	}
  
  oldpar = par(c('mar','xaxs','yaxs'))

  if(is.logical(legendRight))
    legendRight = c(1,0.8)[1+ legendRight ]
  
  if(!is.na(legendRight))
    mar[4] = (1-legendRight) * par('fin')[1]/par('csi')

  par(mar=mar, xaxs='i',yaxs='i',xpd=FALSE)
  
	plot(xpoints,pch=NA, axes=TRUE, ...)

	if(!is.null(ellipse)) {
		plot(ellipse, add=TRUE, ...)
	}
  

  return(invisible(oldpar))
		
}
