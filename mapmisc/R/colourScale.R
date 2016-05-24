roundForBreaks = function(breaks, dec){
	if(is.null(dec)) return(breaks)
dec = 10^dec
breaks = breaks * dec
breaks = c(floor(breaks[1]), 
		round(breaks[seq(2,by=1,len=length(breaks)-2)]),
		ceiling(breaks[length(breaks)]))
breaks / dec
}

radarCol <- c("#FFFFFF", "#99CCFF", "#0099FF", "#00FF66", "#00CC00", "#009900", 
				"#006600", "#FFFF33", "#FFCC00", "#FF9900", "#FF6600", "#FF0000", 
				"#FF0299", "#9933CC", "#660099")

colourScale = function(x=NULL, breaks=5, 
	style=c("quantile","equal","unique", "fixed"),
	col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
	transform=NULL, revCol=FALSE, exclude=NULL, 
	labels=NULL,...) {


	UseMethod("colourScale")	

}

colorScale = function(...) {
	colourScale(...)
}

colourScale.character =   function(x=NULL, breaks=5, 
		style=c("quantile","equal","unique", "fixed"),
		col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
		transform=NULL, revCol=FALSE, exclude=NULL, labels=NULL,...) {

	colourScale(factor(x), breaks=breaks, 
	style=style,
	col=col, opacity=opacity, dec=dec, firstBreak=firstBreak, 
	transform=transform, revCol=revCol, exclude=exclude, labels=labels,...)
}

colourScale.factor = function(x=NULL, breaks=5, 
		style=c("quantile","equal","unique", "fixed"),
		col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
		transform=NULL, revCol=FALSE, exclude=NULL, labels=NULL,...) {

	res=colourScale(as.integer(x), 
			breaks=breaks,
			style="unique",
			col=col, opacity=opacity, dec=dec, firstBreak=firstBreak, 
			transform=transform, revCol=revCol, exclude=exclude, labels=as.character(x),
			...)
	
	res
}

colourScale.Raster = function(x=NULL, breaks=5, 
style=c("quantile","equal","unique", "fixed"),
col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
transform=NULL, revCol=FALSE, exclude=NULL, labels=NULL,...) {

style = style[1]
weights=NULL

NforSample = 5e+05

if(is.character(col)){
if(col[1]=='radar'){
	col = radarCol
}
}


	if(style == "equal") {
		if(length(exclude)) {
			x = unique(x)
		} else {
			x = try(range(c(minValue(x), maxValue(x))), silent=TRUE)
      if(class(x)=='try-error')
        x = range(quantile(sampleRegular(x, size=NforSample)))
		}
	} else if(style=='fixed') {
		x = NULL
} else if(style=='unique') {

    if(is.null(labels)){
      levelsx = levels(x)[[1]]
    } else {
      # if labels is a data frame, use it
      if(is.data.frame(labels)) {
        levelsx = labels
      } else if(
          length(labels) == length(breaks)
          ){
        levelsx = data.frame(
            ID=breaks,
            label=as.character(labels)
            )
      } else { # different number of labels and breaks
        levelsx = data.frame(
            ID=sort(unique(x))
            )
        if(ncol(levelsx)==length(labels)) {
          levelsx$label = as.character(labels)
        } else {
          warning('labels must be same length as either breaks or unique(x)')
          levelsx$label = as.character(levelsx$ID)
        }
      } # end different numbers of levels and breaks
      
  	} # labels is not null
    
    if(ncell(x)<1e+06) {
      x = freq(x)
      weights = x[,2]
      x=x[,1]
    } else {
      x = table(
              sampleRegular(x, NforSample)
      )
      weights = x
      x = as.numeric(names(x))
    }
    
    if(is.null(levelsx)){
      levelsx = data.frame(
          ID=sort(x))
    }
    notInLevels = which(! x %in% levelsx$ID)
    if(length(notInLevels)){
      # add more values to ID
      toAdd = matrix(NA,
          length(notInLevels),ncol(levelsx), 
          dimnames=list(NULL, colnames(levelsx)))
      toAdd[,1] = x[notInLevels]
    levelsx = rbind(levelsx, toAdd)
    levelsx = levelsx[order(levelsx$ID),]
    } # end add more ID in levels
    if(is.vector(labels)){
      if(length(labels)==nrow(levelsx))
        levelsx$label = labels
    } # end labels are vector
    if(is.null(levelsx$label))
      levelsx$label = as.character(levelsx$ID)
    
    levelsx$freq = weights[match(
            levelsx$ID,
            x
            )]
  # if more breaks than ID's have been requested
  # set breaks to all ID's
  if(length(breaks)==1 & all(breaks > nrow(levelsx)))
    breaks = levelsx$ID
  
  if((
        length(breaks)==nrow(levelsx)
        ) & (
        'col' %in% names(levelsx)
        )
  ){
   # colours have been provided 
    col = levelsx$col
  }
  
  weights = levelsx$freq
  x=levelsx$ID
  labels = levelsx$label
		
} else { # not unique or equal or fixed, take a sample
    Nsample = 20000
    xVec= c()
    Diter = 0
    while(length(xVec)<Nsample & Diter < 5){
      xVec = c(xVec,
          na.omit(sampleRegular(x, min(c(ncell(x),Nsample))))
      )
      Diter = Diter + 1
    }
    x = c(xVec, maxValue(x), minValue(x))
} # end if style== stuff

	res=colourScale(x, breaks, 
			style,
			col=col, opacity=opacity, dec=dec, firstBreak=firstBreak, 
			transform=transform, revCol=revCol, exclude=exclude, labels=labels,
			weights=weights,...)
	res[!names(res)%in% 'plot']
}



colourScale.numeric = function(x=NULL, breaks=5, 
		style=c("quantile","equal","unique", "fixed"),
		col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
		transform=NULL, revCol=FALSE, exclude=NULL, labels=NULL, 
		weights=NULL,...) {

	xOrig = x

	if(is.character(col)){
		if(col[1]=='radar'){
			col = radarCol
		}
	}
	
	
	style = style[1]
	if(!is.function(col)){		
		if(is.matrix(col)| is.data.frame(col)) {
			redCol = grep("^[[:space:]]*r(ed)?[[:space:]]*$", colnames(col), value=TRUE,ignore.case=TRUE)
			blueCol = grep("^[[:space:]]*b(lue)?[[:space:]]*$", colnames(col), value=TRUE,ignore.case=TRUE)
			greenCol = grep("^[[:space:]]*g(reen)?[[:space:]]*$", colnames(col), value=TRUE,ignore.case=TRUE)
			if(!all(c(length(redCol),length(greenCol), length(blueCol))==1)) {
				warning("col is a matrix but it's not clear which columns are red, green and blue")
			}		
			if(any(col[,c(redCol,greenCol,blueCol)] > 2 )) {
				theMax = 255
			} else{
				theMax = 1
			}
			col = rgb(col[,redCol], col[,greenCol], col[,blueCol], maxColorValue=theMax)
		}
		
		colString = col
		if(length(colString)==1){
			if(requireNamespace('RColorBrewer',quietly=TRUE)) {
				col = function(n) RColorBrewer::brewer.pal(n, colString)[1:n]
			} else {
				col = function(n) heat.colors(n)
			}
		} else {
			col = function(n) colString[1:n]
		}
	}

	eps=0.01
	
	
	if(style=='unique'){
		# unique breaks	

		if(length(weights)!=length(x)) {
			thetable = as.data.frame(table(ID=x))
			thetable$ID = as.numeric(as.character(thetable$ID)) 
		} else {
			thetable = tapply(weights, x,sum)
			thetable = data.frame(ID=as.numeric(names(thetable)),
					Freq=thetable)
		}
    
    if(length(breaks)==1){
      # breaks is the maximum number of breaks
      thetable = thetable[order(thetable$Freq, decreasing=TRUE),]
      
      shouldExclude = which(thetable$ID %in% exclude)
      if(length(shouldExclude))
        thetable = rbind(
          thetable[-shouldExclude,],
          thetable[shouldExclude,]
          )
      breaks = thetable[
          seq(1,min(breaks,nrow(thetable))),'ID'
          ]
      
    }
    
		notInX = which(! breaks %in% x)
		if(length(notInX))
			thetable = rbind(thetable,
				data.frame(ID=breaks[notInX],
						Freq=rep(0,length(notInX)))
				)
				
		if(length(labels) == length(breaks)) {
			thetable$label = labels[match(thetable$ID,breaks)]
		} else {
			if(length(labels) == length(xOrig)) {
				thetable$label = labels[match(thetable$ID,xOrig)]
			} else if(length(labels)==nrow(thetable)){
				thetable$label = labels[order(thetable$ID)]
			} else {
				thetable$label = as.character(thetable$ID)
			}
		}	

		if(length(breaks)==1) { 
			# assume breaks is the number of unique values to assign colours to
			thetableExc = thetable[!( thetable$ID %in% exclude |
										thetable$label %in% exclude),]
			ncol = min(c(breaks, nrow(thetableExc)))
			breaks = thetableExc[order(thetableExc$Freq,decreasing=TRUE)[1:ncol] , 'ID']
		}  else {
			breaks = breaks[! ( breaks %in% exclude | 
								breaks %in% thetable[thetable$label %in% exclude,'ID']
								) ]
		}
		
		thetable[match( breaks, thetable$ID),'col'] = col(length(breaks))
		thetable = thetable[order(thetable$ID),]

		
		colVec = thetable$col
		names(colVec) =as.character(thetable$ID)
		breaks = thetable$ID
		breaks = c(breaks[1]-1/2, breaks+c(diff(breaks)/2,1/2))
	} else { # not unique breaks
		if(length(exclude) & length(x)) {
			toexclude = which(x %in% exclude)
			x[toexclude]	= NA
		}
		
		
		thetable=NULL
		if(!is.null(transform)) {
			if(is.numeric(transform)) {
				# assume it's a box=cox parameter
				transform = transform[1]
				if(abs(transform)<eps){
					transform="log"
				} else if(abs(transform-1)<eps) {
					transform=NULL
				} else if(abs(transform-0.5)<eps) {
					transform = "sqrt"
				} else if(abs(transform-2)<eps) {
					transform = "square"
				} else {
					if(any(x<0,na.rm=TRUE) ) {
						warning("negative x's given with Box-Cox parameter")
					}
					boxcox = transform		
					transform = list(
							function(x) {
								(-1)^(boxcox<0)*(x^boxcox - 1)/boxcox
							},
							function(x) {
								((-1)^(boxcox<0)*x*boxcox+1)^(1/boxcox)									
							}
					)						
				}
			} # end transform numeric 
			if(is.character(transform)){
				if(transform=="exp") { 
					tranform = list(exp, log)
				} else if(transform=="log") {
					if(any(x<=0,na.rm=TRUE) ) {
						warning("negative or zero x's given with log transform")
					}
					transform = list(log,exp)
				} else if(transform=="sqrt") {
					if(any(x<0,na.rm=TRUE) ) {
						warning("negative x's given with root transform")
					}
					transform = list(sqrt, function(x) x^2)						
				} else if(transform=="square") {
					if(any(x<0,na.rm=TRUE) ) {
						warning("negative x's given with square transform")
					}
					transform = list(function(x) x^2, sqrt)						
				}
			} # end transform character 	
		} else {# end transform not null
			transform = list(function(x) x)[c(1,1)]		
		}
		xOrig = x
		x = transform[[1]](xOrig)
		
		if(style=="quantile"){
			breaks = quantile(x, prob=seq(0,1,len=breaks[1]), na.rm=TRUE)
		} else if(style=="equal"){
			startHere = min(x, na.rm=TRUE)
			if(is.null(transform) & !is.null(firstBreak) )	{
				startHere = firstBreak
				firstBreak = NULL
			}					
			breaks = seq(startHere, max(x, na.rm=TRUE),len=breaks[1])
		} else if(style!='fixed') {
		 # style is passed to classInt
			if (requireNamespace("classInt", quietly = TRUE)) { 
				if(length(x)>2100)
					x = sample(length(x), 2000)
				breaks = classInt::classIntervals(x, n=breaks, 
					style=style, ...)$brks
			} else {
				warning("Install the classInt package to use style=", style)
			}
		} # end classint
		
		breaks = transform[[2]](breaks)
		
		# round
		if(!is.null(dec)) {
			breaks = roundForBreaks(breaks, dec)
		} # end rounding
		if(!is.null(firstBreak))
			breaks[1] = firstBreak
		breaks = sort(unique(breaks))			
		
		colVec = col(length(breaks)-1)		
		if(revCol) colVec = rev(colVec)
	} # end style not unique

		
	if(any(opacity < 1-eps)) {
			if(length(opacity)==1){
				opacity = rep(opacity, length(colVec))
			} else if(length(opacity)==2){
				opacity = seq(opacity[1], opacity[2], len=length(colVec))
			} else if(length(opacity)==3){
				opacity = c(opacity[1], 
						seq(opacity[2], opacity[3],
								len=length(colVec)-1))
			} else {
				opacity = opacity[round(
								seq(1,length(opacity), 
										len=length(colVec)))]
			}
			opacity = toupper(as.hexmode(round(opacity*255)))
			hasOpacity = grep("^#[[:xdigit:]]{8}$", colVec)
			colVec[hasOpacity] = substr(colVec, 1, 7)
			
			isRGB = grep("^#[[:xdigit:]]{6}$", colVec)
			colForPlot = colVec
			colForPlot[isRGB] = paste(colForPlot[isRGB], opacity[isRGB],sep="")
			
	} else {
			colForPlot = colVec
	}
	
	result = list(col=colVec, breaks=breaks, colOpacity=colForPlot)
	if(style=="unique") {
		thetable$colOpacity = colForPlot
		
		result$colourtable = rep(NA, max(thetable$ID)+1)
		result$colourtable[1+thetable$ID] = thetable$colOpacity
		result$colortable = result$colourtable

		if(length(xOrig))
			result$plot = thetable[match(xOrig, thetable$ID),'colOpacity']
		result$levels = thetable
		if(length(thetable$label))
			result$legend = as.character(thetable$label)
	} else if (length(xOrig)){		
		result$plot = as.character(cut(
						xOrig, 
						breaks=breaks,
						labels=colForPlot,
						include.lowest=TRUE
			))
		if(length(labels)==length(result$col))
			result$legend = labels
	}
	
	result
		
}

colourScale.NULL = colourScale.numeric




colourScale.logical = function(x=NULL, breaks=5, 
		style=c("quantile","equal","unique", "fixed"),
		col="YlOrRd", opacity=1, dec=NULL, firstBreak=NULL, 
		transform=NULL, revCol=FALSE, exclude=NULL, labels=NULL, ...) {

	colourScale(as.numeric(x), breaks=breaks, 
			style=style,
			col=col, opacity=opacity, dec=dec, firstBreak=firstBreak, 
			transform=transform, revCol=revCol, exclude=exclude, labels=labels,...)
	
}
		