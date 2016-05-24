
scaleBar = function(crs, 
    pos="bottomright",
    cex=par('cex'),
    pt.cex = 1.1*cex,
    seg.len=5*cex,
    title.cex=cex,
    outer=TRUE, ...) {

  forLegend = list(cex=cex,pt.cex=pt.cex,seg.len=seg.len,...)
  

  if(length(forLegend$scale.cex)){
    warning('use seg.len instead of scale.cex')
    forLegend$seg.len = forLegend$scale.cex
    forLegend$scale.cex=NULL
  }

  
  
  # if text.width=0, no north arrow
  if(pt.cex==0) {
    forLegend$pt.cex = 0.1*forLegend$cex
    noArrow=TRUE
  } else {
    noArrow=FALSE
  }
  if(seg.len==0){
    forLegend$seg.len = 0.1
    noScale = TRUE
  } else {
    noScale = FALSE
  }

  if(!length(forLegend$text.width)){
    forLegend$text.width = 2*strwidth('m', cex=forLegend$pt.cex)
  }
  

	if(is.character(crs))
		crs = CRS(crs)
	if(class(crs) != "CRS")
		crs = CRS(proj4string(crs))
	

#	dash = "\u2517\u2501\u2501\u2501\u2537\u2501\u2501\u2501\u251B"

  # we'll target a label this wide
  widthTargetUsr = strwidth('m', cex=forLegend$cex)*forLegend$seg.len
  
  
	xpoints = t(bbox(extent(par("usr"))))
	xcentre = apply(xpoints, 2, mean)
	xpoints = rbind(centre=xcentre, 
			dashright = xcentre + c(widthTargetUsr,0)
	)
	
	xpoints = SpatialPoints(xpoints, proj4string=crs)

	if(requireNamespace('rgdal', quietly=TRUE)) {	
		xll = spTransform(xpoints, crsLL)
	} else {
		xll= xpoints
		if(!isLonLat(crs))
			warning('rgdal not intalled, assuming the plot is long-lat')
	}

	up = matrix(coordinates(xll)["centre",]+c(0,0.1),ncol=2,
			dimnames=list("up",NULL))
	
	
	xll=rbind(xll, SpatialPoints(up, 
					proj4string=CRS(proj4string(xll))))
	
# how long (in m) is our dashtemplate
	dashdist = spDists(xll[c("centre","dashright"),])[1,2]*1000

	theb = log10(dashdist)
	candidates = 10^c(floor(theb), ceiling(theb))
	candidates = c(candidates[1]*c(1,2,5), candidates[2])
	segdist=candidates[order(abs(candidates - dashdist))[1]]
	
	forLegend$seg.len = forLegend$seg.len *segdist / dashdist
	
	
	if(segdist >1100) {
		lunits="km"
		segdistPrint = segdist / 1000
	} else {
		lunits="m"
    segdistPrint = segdist
	}

	eps = 0.175
	
	dimIn = par("pin")
	dimUser = par("usr")
	dimUser = c(dimUser[2]-dimUser[1], dimUser[4]-dimUser[3])
	InPerUnit = dimIn/dimUser
	

	theN = c(0, 0+1i, eps+1i, (1-eps)+1.5*eps*1i,
			(1-eps)+1i, 1+1i, 1+0i,
			(1-eps)+0i, eps+(1-1.5*eps)*1i,
			eps) - 0.5 - 0.5*1i
	theHat = c(-0.25+1i, 0.5+1.6i, 1.25+1i)
  
  Nwidth = strwidth('N', cex=forLegend$pt.cex)
	theHat = c(theHat, rev(theHat) + 1.5*eps*1i)- 0.5-0.5*1i
	theN =  Nwidth*Re(theN) + 
			1i*Nwidth*InPerUnit[1]/InPerUnit[2]*
			Im(theN)
	theHat =  Nwidth* Re(theHat)+ 
			1i*Nwidth*InPerUnit[1]/InPerUnit[2]*
			Im(theHat)
	
	
	
if(!noScale) {	
	thelabel = paste(segdistPrint, lunits,sep="")
} else {
  thelabel=''
}

	
	defaults = list(col='black', 
			xjust=0.7, 
			x=pos)
	
	if(outer){
		par(xpd=TRUE)
		fromEdge = matrix(par("plt"), 2, 2, 
				dimnames=list(c("min","max"), c("x","y")))
		propIn = apply(fromEdge, 2, diff)
		if(is.character(pos)) {
			inset = c(0,0)
			if(length(grep("^bottom", pos))){
				inset[2] = -fromEdge["min","y"]					
			} else if(length(grep("^top", pos))){
				inset[2] = fromEdge["max","y"]-1					
			}
			if(length(grep("left$", pos))){
				inset[1] = -fromEdge["min","x"]					
			} else if(length(grep("right$", pos))){
				inset[1] = fromEdge["max","x"]-1					
			}
			inset = inset/propIn
			defaults$inset = inset+ 0.01
		}
	} else {
		defaults$inset = 0.01
	}

	for(D in names(defaults)) {
		if(is.null(forLegend[[D]]))
			forLegend[[D]] = defaults[[D]]			
	}
  if(!length(forLegend$text.col))
    forLegend$text.col = forLegend$col

	onem = par("cxy")[1]
	
	forLegend$lty = as.integer(!noScale)
	forLegend$title=' '
	forLegend$legend = ' '
	forLegend$lwd=3

	
	
	thelegend = do.call(graphics::legend, forLegend)
  
  thelegend$title = thelabel
  
  thelegend$textxy =  
       mean(c(thelegend$rect$left, thelegend$text$x)) + 
#       max(c(dashdist, sum(strwidth(c('n',thelabel), cex=cex))))/2 + 
          1i*thelegend$text$y
			
      
      
	if(forLegend$lty)
		text(Re(thelegend$textxy),
        Im(thelegend$textxy),
				label=thelegend$title,  
        pos=3, cex=title.cex,
        offset=0.5,
        col=forLegend$text.col)
		
  # if there's no scale bar or box and pos is numeric,
# put the N at the point spcified
  if(noScale & !nchar(forLegend$title) & all(forLegend$bty=='n') & is.numeric(forLegend$x)) {	
    thecentre =  c(forLegend$x,forLegend$y)[1:2]
  } else {
    thecentre =  c(thelegend$text$x + forLegend$text.width/2, thelegend$text$y)
  }
  
  xpoints = SpatialPoints(t(thecentre),
      proj4string=crs)
  
  if(requireNamespace('rgdal', quietly=TRUE)) {	
    xll = spTransform(xpoints, crsLL)
  } else {
    xll= xpoints
    if(!length(grep("longlat", projection(xpoints))))
      warning('rgdal not intalled, assuming the plot is long-lat')
  }
  xll = rbind(xll,
      SpatialPoints(
          xll@coords+c(0,1),
      proj4string=crs(xll))
  )
  
  if(requireNamespace('rgdal', quietly=TRUE)) {	
    xpoints2 = spTransform(xll, crs)
  } else{
    xpoints2 = xll
  }
  thediff=apply(coordinates(xpoints2), 2,diff)
  north=atan(thediff[1]/thediff[2])+pi*(thediff[2]<0)
  
  theN = theN * exp(-1i*north)
  theHat = theHat * exp(-1i*north)
  
  thecentre = thecentre[1] + 1i*thecentre[2]
  
  if(!noArrow){
	 polygon(forLegend$pt.cex*theN +thecentre, 
			 col=forLegend$text.col,border=NA)
	 polygon(forLegend$pt.cex*theHat + thecentre, 
			 col=forLegend$text.col,border=NA)
  }
 
 
	return(invisible(list(
              out=thelegend, 
              call=c(forLegend, list(label=thelabel)), 
              centre=c(Re(thecentre),Im(thecentre))
  )))
}