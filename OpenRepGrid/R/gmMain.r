#####################################################################
#  	GRAPHICAL MODULES (gm) script collection						
#  	graphical moduls is a collection of simple graphical templates 
#  	that can be used to construct complex custom graphics.
#   All functions from this collection start with the letters gm.
#	Similar to the grid package when a grob is given as output
#	the same functions name ends with Grob and has a corresponding
#	cuntion that does not have grob at the end. e.g. gmFoo and gmFooGrob
# 	
#  	by Mark Heckmann 2009
#####################################################################
# TODO: naming of function and corresponding grob function

# package dependencies:
#require(grid)
#require(colorspace)

#####################################################################
#####################################################################

###### FUNCTION DEFINITION ########
# extract luminance value from hex color value. 
# done by conversion to rgb to lch space
# luminance is returned
# works vectorwise
# default.na is value returned if hex contains NAs

gmGetHexLuminanceValues <- function(hex, default.na =NA)
{
	sapply(hex, function(x){						# check if strings are hex, that is if the start with a "#"
		if(substr(x, 1, 1)!="#" & !is.na(x)) 
			stop("hex is not hexadecimal. getHexLuminanceValues() needs a hex color value!")	
	} )
	NAs <- is.na(hex)								# hex2RGB needs hex values, NAs not accepted 
	hex[NAs] <- "#FFFFFF"							# replace NAs by dummy hex value (white)
	
	lch <- as(hex2RGB(hex), "polarLUV")				# convert hex to rgb to lch	 , thanks to A. Zeileis r-help 20100129
	lum <- as.vector(lch@coords[,"L"])				# return the luminance values only
	lum[NAs] <- default.na
	lum
}

## NOT RUN
# getHexLuminanceValues(c("#3B5715", "#910322"))


###### FUNCTION DEFINITION ########
# select a color from supplied vector corresponding to the 
# luminance value of given hex colors and given breaks
# if hex contains NAs a default hex value can be passed and 
# works vectorwise
# default.na is value returned if hex contains NAs

gmSelectTextColorByLuminance <- function(hex, breaks=c(-1,50,101), breakColors=c("white", "black"), default.na=NA)
{
	luminanceVec <- gmGetHexLuminanceValues(hex)				# get luminance values from hex color				
	indices <- as.integer(cut(luminanceVec, breaks=breaks))		# cut by breaks and get indices
	breakColors[indices]										# return color by index
}

## NOT RUN
# selectTextColorByLuminance(c("#3B5715", "#910322"))

### SHOW EXAMPLE ###
## plot with random background and corresponding textcolor
#library(RColorBrewer)
#bgColors <- c(brewer.pal(8,"Purples"), brewer.pal(8,"YlOrRd"))
#textColors <- gmSelectTextColorByLuminance(bgColors)

#pushViewport(viewport(layout=grid.layout(4, 4, respect=TRUE)))
#	for(i in 1:4){
#		for(j in 1:4){
#			grid.rect(gp=gpar(col="white", fill=bgColors[4*(i-1) + j]), vp=viewport(layout.pos.col=j, layout.pos.row=i))
#			grid.text(paste("Zelle (", i, ",", j, ")", sep=""), gp=gpar(col=textColors[4*(i-1) + j]), vp=viewport(layout.pos.col=j, layout.pos.row=i))
#		}
#	}	
#popViewport()

#####################################################################



#####################################################################
#####################################################################
# like a gmTextBox
# like Murrels example that is rezisable but also viewport 
# rotation enabled allowed



#####################################################################
#####################################################################
# gmTextBox		fill a viewport with color, and text or two texts
#
# evtl mit gpar Objekten??
# wenn zwei Texte übergeben werden, so werden diese, je nachdem, ob horiz=TRUE oder FALSE ist neben 
# oder untereinander dargestellt.
# TODO: ggf. das plotten von borderlines integrieren?
#		Probleme mit den Rändern.

gmTextBox <- function(text=c("text 1", "text 2"), textCol = c("black", "black"), 
					  bgCol = c(grey(.9), grey(.9)), vp=viewport(), textsize= c(.8, .8), 
					  fontface=c("bold", "plain"), horiz=FALSE, vAdjust = c(.4, .65))
{
	#vp=viewport()
	#text=c("text 1", "text 2")
	#textCol = c("black", "black")
	#bgCol = "grey"
	#horiz=FALSE
	#borderCol <- c("black", "black", "black","black")
	#twoTexts <- TRUE
	#####

	if(length(text) == 2 ) {twoTexts <- TRUE} else {twoTexts <- FALSE}				# two text bodies?	
	gpText_1 <- gpar(col=textCol[1], cex=textsize[1], fontface=fontface[1]) 		# make gpar objects
	gpText_2 <- gpar(col=textCol[2], cex=textsize[2], fontface=fontface[2]) 		# 
	gpFill_1 <- gpar(fill=bgCol[1], col=bgCol[1])#, col=NA)											# no border
	gpFill_2 <- gpar(fill=bgCol[2], col= bgCol[2])#, col=NA)											# no border

	# three options: 1 text, 2 texts vertical, 2 texts horizontal
	
	if(!twoTexts){																	# just one text body
		pushViewport(vp)
			grid.rect(gp=gpFill_1)
			grid.text(text[1], gp=gpText_1)
		popViewport()
	}
	
	# TODO: dieser Ansatz ist noch nicht perfekt. Speziell, da ich verschiedene fontfaces für
	# die obere und untere Zelle haben möchte müssen zwei textGrobs gebaut werden. Hier muss 
	# noch ein wenig Arbeit geleistet werden, um deren Größe zu messen unds ie nebeneinander 
	# sauber zu platzieren. Vlt. in der nächsten Version.
	if(twoTexts){																			# two text bodies
		pushViewport(vp)																	# outer viewport
			if(horiz){ nRow <- 1; nCol <- 2; yOffset=.5 } else {nRow <- 2; nCol <- 1; yOffset=vAdjust[1]}					# define layout with respect to orientation (horiz T/F)
			pushViewport(viewport(layout=grid.layout(nRow,nCol)))							# split viewport horizontally or vertically (horiz T/F)
				posRow <- 1; posCol <- 1 													# define row and column position of first viewport
				pushViewport(viewport(layout.pos.row=posRow, layout.pos.col=posCol))		# push upper viewport
					grid.rect(gp=gpFill_1)
					grid.text(y=yOffset, text[1], gp=gpText_1, just=c("center", "center"))
				popViewport()
				if(horiz){ posRow <- 1; posCol <- 2; yOffset=.5} else { posRow <- 2; posCol <- 1; yOffset=vAdjust[2]}		# define row and column position of second viewport
				pushViewport(viewport(layout.pos.row=posRow, layout.pos.col=posCol))		# push lower viewport
					grid.rect(gp=gpFill_2)
					grid.text(y=yOffset, text[2], gp=gpText_2, just=c("center", "center"))
				popViewport()		
			popViewport()
			grid.rect(gp=gpar(col="white", lwd=2))											# border around whole vp
		popViewport()
	}											
}


### NOT RUN ###
#gmTextBox()
#gmTextBox(c(12,"(14 %)"))
#gmTextBox(c(12))

## make a grid of gmTextBoxes
#pushViewport(viewport(layout=grid.layout(4, 4, respect=FALSE)))
#	for(i in 1:4){
#		for(j in 1:4){
#			gmTextBox(vp=viewport(layout.pos.row=i, layout.pos.col=j))
#		}
#	}	
#popViewport()

## use data from data frame
#script <- matrix(sample(1:100, 16), ncol=4)
#subscript <- matrix(paste("(", sample(1:100, 16), "%)", sep=""), ncol=4)
#nCol <- ncol(script); nRow <- nrow(script) 
#pushViewport(viewport(layout=grid.layout(nRow, nCol, respect=FALSE)))
#	for(i in 1:nRow){
#		for(j in 1:nCol){
#			gmTextBox(c(script[i,j], subscript[i,j]), horiz=FALSE, vp=viewport(layout.pos.row=i, layout.pos.col=j))
#		}
#	}	
#popViewport()


## plot with random background colors and corresponding textcolor overlay
#library(RColorBrewer)

#bgColors <- brewer.pal(8,"YlOrRd")
#script <- matrix(sample(1:100, 16), ncol=4)
#subscript <- matrix(paste("(", sample(1:100, 16), "%)", sep=""), ncol=4)
#bgColorsScript <- matrix(sample(bgColors, 16, rep=T), ncol=4)
#textColorsScript <- matrix(gmSelectTextColorByLuminance(bgColorsScript), ncol=4)  # benutzt gmSelectTextColorByLuminance

#nCol <- ncol(script); nRow <- nrow(script) 
#pushViewport(viewport(layout=grid.layout(nRow, nCol, respect=TRUE)))
#	for(i in 1:nRow){
#		for(j in 1:nCol){
#			gmTextBox(text=c(script[i,j], subscript[i,j]), textCol=c(textColorsScript[i,j], textColorsScript[i,j]), 
#					 bgCol=c(bgColorsScript[i,j], bgColorsScript[i,j]), vp=viewport(layout.pos.row=i, layout.pos.col=j))
#		}
#	}	
#popViewport()

## uses gmRandomColor
#nRow <- 5; nCol <- 5
#script <- matrix(sample(1:100, nRow*nCol), ncol=nCol)
#subscript <- matrix(paste("(", sample(1:100, nRow*nCol), "%)", sep=""), ncol=nCol)
#pushViewport(viewport(layout=grid.layout(nRow, nCol, respect=TRUE)))
#	for(i in 1:nRow){
#		for(j in 1:nCol){
#			randColor <- gmRandomColor()
#			gmTextBox(text=c(script[i,j], subscript[i,j]), textCol=rep(gmSelectTextColorByLuminance(randColor),2), 
#					 bgCol=rep(randColor, 2), vp=viewport(layout.pos.row=i, layout.pos.col=j))
#		}
#	}	
#popViewport()


#####################################################################
#####################################################################
# gmSplitTextGrob
# text grob that automatically does line breaks in text, allows resizing
# and vertical orientation of
# TODO: - ggf. muss man sich überlegen, ob die momentane Form bei horiz=F
# 		geeignet ist. Denn man muss da umdenken. x muss nun .5 sein anstatt 0
# 		wenn der text nicht auch z.B an den linken Rand soll.
# 		- no vectorized form available yet

# adopted from Murrell(2008) R Graphics, p...
gmSplitString  <- function(text, horiz=TRUE, splitWidth=unit(.98, "npc"))			# function to split grobtext
{				
	#require(grid)
	if(is.expression(text)){			# Expressions können nicht weiter verarbeitet werden
		return(text)
		break
	}
	if(is.null(text)) text <- ""
	if(length(text) ==1 & is.na(text)) text <- ""
	if(is.character(text) & length(text)==0) text <- ""
	if(text==""){
		return(paste(text))
		break
	}
	strings <- strsplit(as.character(text), " ")[[1]]
	if(length(strings)==1){
		return(paste(strings))
		break
	}
	newstring <- strings[1]
	linewidth <- stringWidth(newstring)
	gapwidth <-  stringWidth(" ") 

	if(!horiz){
		availwidth <- convertHeight(splitWidth, "inches", valueOnly=TRUE)
	}else{
		availwidth <- convertWidth(splitWidth, "inches", valueOnly=TRUE)
	}
	#print(availwidth);
	for (i in 2:length(strings)){
		width <- stringWidth(strings[i])
		if (convertWidth(linewidth + gapwidth + width,
			"inches", valueOnly=TRUE) < availwidth){
			sep <- " "
			linewidth <- linewidth + gapwidth + width
		} else {
			sep <- "\n"
			linewidth <- width
		}
		newstring <- paste(newstring, strings[i], sep=sep)
	}
	newstring
}


# make text grob
gmSplitTextGrob <- function(text, x=unit(0.5, "npc"), y=unit(0.5, "npc"), just=c("center", "center"),  gp=gpar(), horiz=TRUE, splitWidth=unit(.98, "npc"), ...)
{						
    if (!is.unit(splitWidth)) splitWidth <- unit(splitWidth, "npc")

	if(!horiz) rot <- 90 else rot <- 0
	#print(horiz); print(rot);
	grob(text=text, cl="gmSplitTextGrob", x=x, y=y, just=just, rot = rot, horiz=horiz, gp=gp, splitWidth=splitWidth, ...)
}


# variation to explore
drawDetails.gmSplitTextGrob <- function(x, recording)    	# drawdetails method is called when resizing window
{				
    #str(x);
	if(!x$horiz) {
		grid.text(label=gmSplitString(x$text, horiz=x$horiz, splitWidth=x$splitWidth), 
				  rot=x$rot, just=x$just, x=x$x, y=x$y, gp=x$gp)
	} else {
		grid.text(label=gmSplitString(x$text, horiz=x$horiz, splitWidth=x$splitWidth), 
				  rot=x$rot, just=x$just, x=x$x, y=x$y, gp=x$gp,)	
	}
}


# printing wrapper for gmSplitTextGrob 
gmSplitTextBox <- function(text, x=unit(0.5, "npc"), y=unit(0.5, "npc"), just=c("center", "center"),  gp=gpar(), horiz=TRUE, splitWidth=unit(.98, "npc"), ...)
{
	tg <- gmSplitTextGrob(text, x=x, y=y, just=just, gp=gp, horiz=horiz, splitWidth=splitWidth, ...)  			# gmSplitTextGrob
	grid.draw(tg)																								# print gmSplitTextGrob
}	



### NOT RUN
#text <- "some random longer text that might be the label of an item"
#grid.draw(gmSplitTextGrob(text, horiz=T, just=c("center", "center")))
#grid.draw(gmSplitTextGrob(text, horiz=T, splitWidth=.9, x=0.05, y=.5, just=c("left", "center"), gp=gpar(col="darkgrey")))
#gmSplitTextBox(text, horiz=T, splitWidth=.9, x=0.05, y=.5, just=c("left", "center"), gp=gpar(col="darkgrey", lineheight=.8))

#gmSplitTextBox(text, h=F)

#splitText <- gmSplitTextGrob(text, horiz=F, class="gmSplitTextGrob", gp=gpar(fontsize=12, lineheight=.9))
#grid.draw(splitText)

## matrix of text with random orientation
#grid.newpage()
#text <- "some random longer text that might be the label of an item"
#textOrientation <- matrix(sample(c(T,F), 16, rep=TRUE), ncol=4)
#pushViewport(viewport(layout=grid.layout(4, 4, respect=FALSE)))
#	for(i in 1:4){
#		for(j in 1:4){
#			grid.draw(gmSplitTextGrob(text, horiz=textOrientation[i,j],
#								vp=viewport(layout.pos.row=i, layout.pos.col=j),
#								gp=gpar(fontsize=12, lineheight=.9)))
#			grid.rect(vp=viewport(layout.pos.row=i, layout.pos.col=j))											
#		}
#	}	
#popViewport()



## matrix of text with random orientation and random fore- and background color
#grid.newpage()
#text <- "some random longer text that might be the label of an item"
#textOrientation <- matrix(sample(c(T,F), 16, rep=TRUE), ncol=4)
#pushViewport(viewport(layout=grid.layout(4, 4, respect=FALSE)))
#	for(i in 1:4){
#		for(j in 1:4){
#			randColor <- gmRandomColor()
#			grid.rect(vp=viewport(layout.pos.row=i, layout.pos.col=j), 
#					  gp=gpar(fill=randColor, col="lightgrey")))
#			gmSplitTextBox(text, splitWidth=.9, horiz=textOrientation[i,j], vp=viewport(layout.pos.row=i, layout.pos.col=j),
#						   gp=gpar(fontsize=12, lineheight=.9, col=gmSelectTextColorByLuminance(randColor)))							
#		}
#	}	
#popViewport()


#####################################################################
#####################################################################
# gmMakeVpBorders
# uses the current vp and adds border lines at specified places
# is useful for the construction of tables etc. and an alternative
# to do this afterwards by whole lines.
# be careful with the visually adequate order of the sides as they are printed in the order given by side
# TODO: - recycle vector in lwd etc.? Evtl. kann ein NA stattdessen lieber 
# 		dazu genutzt werden, dass die Linie nicht gezeichnet wird.
# 		- evtl. noch nicht ganz perfekt in bezug auf das clipping, da dies nicht
# 		als Funtkionsargument implementiert ist

gmMakeVpBorders <- function(side, col, lwd, ...)
{
	#col <- gmRandomColor(4)
	#side <- 1:4
	#lwd <- 50:54
	for(i in side){
		if(i==1 | i=="bottom") grid.lines(x=c(0,1), y=c(0,0), gp=gpar(lwd=lwd[side %in% i], col=col[side %in% i], lineend="square"), ...)
		if(i==2 | i=="left") grid.lines(x=c(0,0), y=c(0,1), gp=gpar(lwd=lwd[side %in% i], col=col[side %in% i], lineend="square"), ...)
		if(i==3 | i=="top") grid.lines(x=c(0,1), y=c(1,1), gp=gpar(lwd=lwd[side %in% i], col=col[side %in% i], lineend="square"), ...)
		if(i==4 | i=="right") grid.lines(x=c(1,1), y=c(0,1), gp=gpar(lwd=lwd[side %in% i], col=col[side %in% i], lineend="square"), ...)
	}	
}

#gmSplitTextBox("Some long text is written here", splitWidth=.9)
#gmMakeVpBorders(1:4, gmRandomColor(4), lwd=rep(30,4))      

## matrix of text with random orientation and random fore- and background color
#grid.newpage()
#text <- "some random longer text that might be the label of an item"
#textOrientation <- matrix(sample(c(T,F), 16, rep=TRUE), ncol=4)
#pushViewport(viewport(layout=grid.layout(4, 4, respect=FALSE)))
#	for(i in 2:3){
#		for(j in 2:3){
#			tmpVp <- viewport(layout.pos.row=i, layout.pos.col=j)
#			grid.rect(gp=gpar(fill=gmRandomColor()), vp=tmpVp)
#			gmSplitTextBox(text, splitWidth=.9, horiz=textOrientation[i,j], vp=tmpVp)
#			gmMakeVpBorders(1:4, rep("grey", 4), lwd=rep(10,4), vp=tmpVp)      						
#		}
#	}	
#popViewport()

## booktab like look
#grid.newpage()
#text <- "some random longer text that might be the label of an item"
#textOrientation <- matrix(sample(c(T,F), 16, rep=TRUE), ncol=4)
#pushViewport(viewport(layout=grid.layout(4, 4, respect=FALSE)))
#	for(i in 1:4){
#		for(j in 2:3){
#			tmpVp <- viewport(layout.pos.row=i, layout.pos.col=j, clip=T)			# im Moment noch clipping bei vp definition. Lieber direkt in Funktion
#			grid.rect(gp=gpar(fill=gmRandomColor(v=.8), col=NA), vp=tmpVp)
#			gmSplitTextBox(text, splitWidth=.9, horiz=textOrientation[i,j], vp=tmpVp)
#			gmMakeVpBorders(1:4, rep("black", 4), lwd=c(4,NA, 4, NA), vp=tmpVp)      						
#		}
#	}	
#popViewport()



#####################################################################
#####################################################################
# gmBulletPointsBox					
# A function that prints a list of text elements as bullet points
# Bullets can be chosen any pch, numbers, letters or any other vector.



#####################################################################
#####################################################################
# gmProfileLines
# ask Hadley first if he already implicitly has it...



#####################################################################
#####################################################################
# gmRandomColor 
# small convenience wrapper that returns a vector of random colors as hex
# unsing HSV scheme (hue, saturation, value)
# (requires RColorBrewer package).
# hue (1-360), saturation 0-1 and value 0-1 can be fixed or restricted to a range
# shuffle = shuffle the outoput vector, so patterns are destroyed

gmRandomColor <- function(n=1, h=runif(n)*360, s=runif(n), v=runif(n), shuffle=TRUE, plot=FALSE )
{
	#require(colorspace)
	#hexColorVec <- hex(HSV(runif(n), runif(n), runif(n)))
	hexColorVec <- hex(HSV(h, s, v))
	if(shuffle) hexColorVec <- hexColorVec[sample(seq_along(hexColorVec), length(hexColorVec))]
	if(plot){
		pal <- function(col, border = "light gray", ...) 
		{
			n <- length(col) 
			plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
				axes = FALSE, xlab = "", ylab = "", ...) 
			rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
		}
		pal(hexColorVec)
	}
	return(hexColorVec)
}

### NOT RUN
# gmRandomColor() 
# gmRandomColor(20, plot=T)
# gmRandomColor(30, h=100:200, v=3:10/10, p=T, shuffle=F)
# gmRandomColor(30, h=100:200, v=3:10/10, p=T)



#####################################################################
#####################################################################
# gmArrowIndicator
# an arrow of given size, angle, filling and background color which can be 
# used to visuallize changes or rates.
# angle
# initAngle
# size in mm
# fill
# border
# background

## deprecated due to grob version below
#gmArrowIndicator <- function(angle=0, col="black", size=5, circle=FALSE, ...)
#{
#	#gp <- modifyList(gpar(fill="black", col=NA, lwd=1, lineend ="square",   		# set default gpar and overwrite if provided
#	#				linejoin ="mitre", linemitre=1), gp)
#
#	if(hasArg(vp)) vp <- list(...)$vp else vp <- viewport()							# if vp is passed use it else create empty viewport
#		
#	pushViewport(vp)
#		pushViewport(viewport(angle=angle))
#		if(circle) grid.circle(x=0.5, y=0.5, r=unit(size+2,"mm"), gp=gpar(fill="lightgrey", col=NA))
#		
#		arrow <- arrow(angle = 35, length = unit(size, "mm"), type = "closed")			# make arrow object to be passed to grid.lines
#		grid.lines( x = unit(c(.5, .5), "npc") + unit(c(0,size), "mm"), 
#					y = unit(c(.5, .5), "npc"), 
#					arrow = arrow, 
#					gp=gpar(fill=col, lwd=1, col=NA, lineend ="square",
#							linejoin ="mitre", linemitre=1))              
#		grid.lines( x = unit(c(.5, .5), "npc") + unit(c(-size/3,0), "mm"), 
#					y = unit(c(.5, .5), "npc"), 
#					gp= gpar(col=col, lwd=2*size, lineend ="square", linejoin ="mitre"))
#		popViewport()
#	popViewport()
#}



gmArrowIndicatorGrob <- function(angle=0, col="black", size=5, circle=FALSE, initangle=0,  ...)
{
	#gp <- modifyList(gpar(fill="black", col=NA, lwd=1, lineend ="square",   		# set default gpar and overwrite if provided
	#				linejoin ="mitre", linemitre=1), gp)
	vp <- viewport(angle=angle + initangle)											# created rotated viewport for arrow direction
	if(hasArg(vp)) vp <- vpStack(list(...)$vp, vp)									# if vp is passed use it and stack the two viewports
	if(circle) 
		circleBackgroundGrob <- circleGrob(x=0.5, y=0.5, r=unit(size+2,"mm"), gp=gpar(fill="lightgrey", col=NA))
	else
		circleBackgroundGrob <- nullGrob()
	arrow <- arrow(angle = 35, length = unit(size, "mm"), type = "closed")			# make arrow description object to be passed to grid.lines
	gTree(children=gList(
		circleBackgroundGrob,
		linesGrob(  x = unit(c(.5, .5), "npc") + unit(c(0,size), "mm"), 
				    y = unit(c(.5, .5), "npc"), 
				    arrow = arrow, 
				    gp=gpar(fill=col, lwd=1, col=NA, lineend ="square",
				    	   	linejoin ="mitre", linemitre=1),
				    vp=vp),              
		linesGrob(  x = unit(c(.5, .5), "npc") + unit(c(-size/3,0), "mm"), 
				    y = unit(c(.5, .5), "npc"), 
					gp= gpar(col=col, lwd=2*size, lineend ="square", linejoin ="mitre"),
					vp=vp)
					)
	)	
}

gmArrowIndicator <- function(angle=0, col="black", size=5, circle=FALSE, initangle=0, ...){
	aiGrob <- gmArrowIndicatorGrob( angle=angle, col=col, size=size, 
									circle=circle, initangle=initangle, ...)
	grid.draw(aiGrob)
}

## NOT RUN:
# gmArrowIndicator2(angle=10, vp=viewport(x=.9))
## array of arrows
# grid.newpage()
# angleMatrix <- matrix(sample(1:360, 16, rep=T), ncol=4)
# pushViewport(viewport(layout=grid.layout(4, 4, respect=FALSE)))
# 	for(i in 1:4){
# 		for(j in 1:4){
# 			vpTmp <- viewport(layout.pos.row=i, layout.pos.col=j)
# 			grid.rect(vp=vpTmp)
# 			gmArrowIndicator(angleMatrix[i,j], vp=vpTmp)						
# 		}
# 	}	
# popViewport()
# 
# ## array of arrows colored by angle
# grid.newpage()
# angleMatrix <- matrix(sample(1:360, 100, rep=T), ncol=10)
# pushViewport(viewport(layout=grid.layout(10, 10, respect=FALSE)))
# 	for(i in 1:10){
# 		for(j in 1:10){
# 			vpTmp <- viewport(layout.pos.row=i, layout.pos.col=j)
# 			grid.rect(vp=vpTmp)
# 			gmArrowIndicator(angleMatrix[i,j], col=gmRandomColor(), size=7, vp=vpTmp)						
# 		}
# 	}	
# popViewport()
# 
# ## array of arrows colored by angle
# grid.newpage()
# angleMatrix <- matrix(1:100*3.6, ncol=10, byrow=T)
# pushViewport(viewport(layout=grid.layout(10, 10, respect=FALSE)))
# 	for(i in 1:10){
# 		for(j in 1:10){
# 			vpTmp <- viewport(layout.pos.row=i, layout.pos.col=j)
# 			grid.rect(vp=vpTmp)
# 			col <- gmSelectColorByValue(angleMatrix[i,j], seq(0, 360, by=10))
# 			gmArrowIndicator(angleMatrix[i,j], col=col, size=7, vp=vpTmp)						
# 		}
# 	}	
# popViewport()
# 
# ## a frameGrob example
# rows <- 10; cols <- 10
# fg <- frameGrob(layout=grid.layout(rows,cols, widths=unit(rep(1.5,cols), "cm"), heights=unit(rep(1.5,rows), "cm")))
# for(i in 1:rows) for (j in 1:cols)	fg <- placeGrob(fg, gmArrowIndicatorGrob(), i, j)
# grid.draw(fg)
# 
# ## a frameGrob example
# rows <- 10; cols <- 10
# angleMatrix <- matrix(1:(rows*cols)*3.6, ncol=cols, byrow=T)
# gmSelectColorByValue(angleMatrix[i,j], seq(0, 360, by=10))
# fg <- frameGrob(layout=grid.layout(rows,cols, widths=unit(rep(1.5,cols), "cm"), heights=unit(rep(1.5,rows), "cm")))
# for(i in 1:rows) {
# 	for (j in 1:cols){
# 		col <- gmSelectColorByValue(angleMatrix[i,j], seq(0, 360, by=10))
# 		fg <- placeGrob(fg, gmArrowIndicatorGrob(angle=angleMatrix[i,j], col=col), i, j)
# 	}	
# }
# grid.draw(fg)


#####################################################################
#####################################################################
# gmShowPalette 
# convenient wrapper to look at a palett, taken from colorspace vignette

#gmShowPalette <- function(col, border = "light gray", ...) 
#{
#	n <- length(col) 
#	plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),axes = FALSE, xlab = "", ylab = "", ...) 
#	rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
#}

# gmShowPalette(diverge_hcl(30, h = c(120, 20), c = 70, l = c(55, 98)))

# grid based version for better placement
gmShowPalette <- function(col, border = "light gray", ...) 
{
	layout <- grid.layout(ncol=length(col))
	fg <- frameGrob(layout=layout, ...)
	for(i in seq_along(col))
		fg <- placeGrob(fg, rectGrob(gp=gpar(fill=col[i], col=border)), col=i)
	grid.draw(fg)		
}
# gmShowPalette(diverge_hcl(30, h = c(120, 20), c = 70, l = c(55, 98)))


#####################################################################
#####################################################################
# gmSelectColorByValue
# TODO: vectorize, work on matrix, df etc.
# if value does not lie within any interval defined by breaks NA is returned
# if x contains NAs default.na is returned. The default is NA
# but a color can be specified in case needed.

gmSelectColorByValue <- function(x, breaks= seq(0, 100, by=10), 
								 colors=diverge_hcl(length(breaks)-1, h = c(120, 20), c = 70, l = c(55, 98)),
								 default.na=NA)
{
	if(length(breaks)!=(length(colors)+1)) 
			stop("breaks and colors have to be the same length!") 	# check if vectors have same length
	is.na(x) <- is.na(x)											# replacec NaNs by NA

	x <- as.vector(x)
	col <- cut(x, breaks=breaks, labels =colors)
	col <- colors[as.integer(col)]
	col[is.na(col)] <- default.na									# replace NAs by default.na
	col
}
# evtl. mal cut2 aus Hmisc anschauen

#gmSelectColorByValue(1:100)
#gmSelectColorByValue(c(NA, NA, 1:9))
#gmSelectColorByValue(c(NA, NA, 1:9), default.na="#EDEBEB")

#tmp <- gmSelectColorByValue(1:100, breaks= seq(0, 100, by=5))
#print(tmp)
#gmShowPalette(tmp)

#tmp <- gmSelectColorByValue(1:100, c(0,50,100), c("black", "white"))
#print(tmp)
#gmShowPalette(tmp)


#####################################################################
#####################################################################
# gmLegends
# there is a function for grid legends in vcd package, but it does
# not allow for multiple characters. The code is slightly modified
# allowing multiple characters in the first row

# WORKING PARTLY BUT STILL UNDER CONSTRUCTION!!!

#library(vcd)
#grid_legend(0.8, 0.9, c("aa","bb"), c("blue", "blue"), c("Port", "Starboard"), title = "SIDE")

#grid_legend(0.8, 0.9, c(1, 19), c("red", "blue"),
#  c("Port", "Starboard"), title = "SIDE")

#x=.5
#y=.5
#pch=c(1,2)
#col="black"
#labels=c("Text 1", "Text2")
#frame = TRUE
#hgap = unit(0.5, "lines")
#vgap = unit(0.3, "lines")
#default_units = "lines"
#gp = gpar()
#draw = TRUE
#title = "Legend:"

# TODO: automatic deterination of wFirstRow by max stringwidth

gmLegend <- function (x, y, pch, symbol=FALSE, col, labels, hgap = unit(0.5, 
    "lines"), wFirstCol=unit(2,"lines"), vgap = unit(0.3, "lines"), default_units = "lines", 
    gpRect = gpar(), gpText=gpar(), draw = TRUE, title = "Legend:") 
{
    labels <- as.character(labels)
    if (is.logical(title) && !title) 
        title <- NULL
    if (!is.null(title)) {
        labels <- c(title, labels)
        pch <- c(NA, pch)
        col <- c(NA, col)
    }
    nkeys <- length(labels)
    if (length(pch) != nkeys) 				
        stop("pch and labels not the same length")
    if (!is.unit(hgap)) 
        hgap <- unit(hgap, default_units)
    if (length(hgap) != 1) 
        stop("hgap must be single unit")
    if (!is.unit(vgap)) 
        vgap <- unit(vgap, default_units)
    if (length(vgap) != 1) 
        stop("vgap must be single unit")
	
    legend.layout <- grid.layout(nkeys, 3, 
			widths = unit.c(wFirstCol, max(unit(rep(1, nkeys), 
							"strwidth", as.list(labels))), hgap), 
			heights = unit.pmax(unit(1, "lines"), vgap + unit(rep(1, nkeys), 
								"strheight", as.list(labels))))
    fg <- frameGrob(layout = legend.layout, gp = gpText)
	# background col
	fg <- placeGrob(fg, rectGrob(gp = gpRect)) 
	
    for (i in 1:nkeys) {
        tit <- !is.null(title) && i == 1
        if (!tit)
			if(symbol) { # print text if symbol is FALSE (default) 
				fg <- placeGrob(fg, pointsGrob(0.5, 0.5, pch = pch[i], 
				                gp = gpar(col = col[i])), col = 1, row = i)
			} else {
	        	fg <- placeGrob(fg, textGrob(label= pch[i], x=0.1, y=0.5, 
	                gp = gpar(col = col[i]), just = c("left", "center")), 
					col = 1, row = i)		
			}
        fg <- placeGrob(fg, textGrob(labels[i], x = 0 + 0.3 * 
            tit, y = 0.5, just = c("left", "center")), col = 2 - 
            tit, row = i)
    }
    pushViewport(viewport(x, y, height = unit(nkeys, "lines"), 
        width = grobWidth(fg)))
# 	if (frame)
#		fg <- placeGrob(fg, rectGrob(gp = gpar(fill = "transparent"))) 
	if (draw) 
        grid.draw(fg)
    popViewport(1)
    invisible(fg)
}

#labels=1:20
# legend with symbols
#gmLegend(x=0.25, y=0.5, symbol=TRUE, pch=seq_along(labels), col=rep("blue", length(labels)),
#  		 labels=labels, gpText=gpar(cex=.7), title = "Test 1")
# without frame
#gmLegend(x=0.5, y=0.5, symbol=TRUE, pch=seq_along(labels), col=rep("brown", length(labels)),
#  		 labels=labels, gpText=gpar(cex=.7, col=grey(.7)), gpRect=gpar(col=NA), title = "Test 2")
# legend with multi-character index column
#gmLegend(x=0.75, y=0.5, wFirstCol=unit(3, "lines"), hgap = unit(1, "lines"), 
#		 pch=paste(LETTERS[1:20], labels, sep=""), labels=LETTERS[1:20], col=rep(rainbow(20), length(labels)), 
#		 gpRect=gpar(col=1, fill=grey(.95), lty=3), gpText=gpar(col=grey(.5), cex=.7), title = NULL)



#####################################################################
#####################################################################
# gmLegends_2
# there is a function for grid legends in vcd package, but it does
# not allow for multiple characters. The code is slightly modified
# allowing multiple characters in the first row

# TODO: placeGrob nutzen sowie Spalten und Zeilenhöhe berechnen!
# 		background wird im moment einfach noch gezeichnet ohne args!

gmLegend2 <- function(colors, labels, ncol=NA, nrow=NA, byrow=TRUE, 
					  symbolSize=unit(3, "mm"), symbolMargin=unit(2, "mm"),
					  bg=1, na.bg =TRUE, force.height=FALSE, dynamic=TRUE)
{
#args:
#byrow=T
#ncol=NA
#nrow=1
#colors <- rainbow(10)
##labels <- LETTERS[1:10]
#labels <- sapply(1:10, getRandString)
##args:
#symbolSize <- unit(3, "mm")
#symbolMargin <- unit(2, "mm")
#bg <- 1
#na.bg =TRUE				# draw background in NA (unused) cells?

	# input check
	if(length(labels) != length(colors))					# do colors match labels?
		stop("Same length of colors and labels required!")
	if(length(colors)==1 & length(labels) > 1)				# one color many labels -> recycle color
		colors <- rep(colors, length(labels))
	if(sum(is.na(c(ncol, nrow))) != 1 | 
		!is.numeric(c(ncol, nrow))) stop("Please specify ncol OR nrow as positive integer.")
	noCells <- length(labels)

	# calcs: determine matrix size
	if(is.na(nrow)) nrow <- noCells %/% ncol + (noCells %% ncol != 0) 	# needed no of rows
	if(is.na(ncol)) ncol <- noCells %/% nrow + (noCells %% nrow != 0) 	# needed no of cols 

	# missing cells are given NAs
	if(noCells != ncol * nrow){
		labels <- c(labels, rep(NA, ncol * nrow - noCells))
		colors <- c(colors, rep(NA, ncol * nrow - noCells))
	}
	labelsMat <- matrix(labels, ncol=ncol, nrow=nrow, byrow=byrow)
	colorsMat <- matrix(colors, ncol=ncol, nrow=nrow, byrow=byrow)
	#labelsMat; colorsMat

	# filling the layout
	labelCell <- function(label){
		gTree(children=gList(
			  gmSplitTextGrob(label,
			    				x=unit(2, "mm"), 
								y=unit(.5, "npc"), 
								just=c("left", "center"), 
								gp=gpar(lineheight=.7, cex=.8))
			  ))
	}

	symbolCell <- function(fill, col="black"){
		gTree(children=gList(
			  rectGrob(width=symbolSize, height=symbolSize, 
					   gp=gpar(fill=fill, col=col))
		))
	}	

	backgroundCell <- function(gp=gpar()){
		gTree(children=gList(
			  rectGrob(width=1, height=1, gp=gp)
		))
	}	
	
	# make layout and frame
	layout <- grid.layout(nrow=nrow, ncol = ncol*2,
	 					  widths=unit(rep(c(7,1), ncol),  rep(c("mm", "null"), ncol)),
						  heights=unit(rep(1, nrow), "lines"))
	fg <- frameGrob(layout=layout, name="topFrame")

	# make and add background object
	bgCell <- rectGrob(gp=gpar(fill=grey(0.6), col="white", lwd=5))
	fg <- packGrob(fg, bgCell, dynamic=dynamic, force.height=force.height)	

	# fill frame
	for(i in 1:nrow){
		for(j in 1:(ncol)){
			bgGrob <- backgroundCell(gpar(fill=grey(.95), col="white"))
			draw.bg <- !(is.na(labelsMat[i,j]) & !na.bg)
			if(bg==1 & draw.bg)
				fg <- packGrob(fg, bgGrob, col=(2*j-1):(2*j), row=i, dynamic=dynamic, force.height=force.height)	
			if(bg==2 & draw.bg){
				fg <- packGrob(fg, bgGrob, col=(2*j-1), row=i, dynamic=dynamic, force.height=force.height)	
				fg <- packGrob(fg, bgGrob, col=(2*j), row=i, dynamic=dynamic, force.height=force.height)	
			}
			
			if(!is.na(colorsMat[i,j])){
				symbolGrob <- symbolCell(colorsMat[i,j], col=NA) 
				fg <- packGrob(fg, symbolGrob, col=2*j-1, row=i, dynamic=dynamic, force.height=force.height)	
			}
			if(!is.na(labelsMat[i,j])){
				cellGrob <- labelCell(labelsMat[i,j])	
				fg <- packGrob(fg, cellGrob, col=(2*j), row=i, dynamic=dynamic, force.height=force.height)
			}
		}
	}
	return(fg)
}


#fg <- gmLegend2(rainbow(7), letters[1:7], ncol=3, bg=0)
#pushViewport(viewport(y=.3,height=unit(10, "mm"), width=.4))
#	grid.draw(fg)
#popViewport()

#getRandString <- function(len=12) return(paste(sample(c(rep(0:9, each=5), LETTERS,letters, rep(c(" "), 10)),len,replace=TRUE),collapse=''))
#fg <- gmLegend2(rainbow(10), sapply(1:10, getRandString), ncol=3, byrow=F)
#pushViewport(viewport(y=.3,height=unit(10, "mm"), width=.4))#
#	grid.draw(fg)
#popViewport()



















