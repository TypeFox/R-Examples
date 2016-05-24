

### daltonize ###


#' @export
#' @name daltonize
#' @aliases LSD.daltonize
#' @title Dichromat vision simulation for colorpalettes
#' @description Dichromat vision simulation and enhancement according to \url{http://www.daltonize.org}.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param show logical: if \code{TRUE} (by default), the resulting colorpalettes are depicted in an R plot.
#' @return \code{daltonize} returns a list, where each entry is a vector containing R built-in colors in hexadecimal representation: \item{simulated}{vector of simulated colors} \item{enhanced}{vector of enhanced colors}
#' @author Bjoern Schwalb
#' @seealso \code{\link{disco}}, \code{\link{colorpalette}}, \code{\link{demotour}}
#' @references \url{http://www.daltonize.org}
#' @examples daltonize("heat",cvd = "d")
#' daltonize("colorblind",cvd = "p")
#' @keywords daltonize, colorblind


daltonize = function(colpal,cvd = "p",show = TRUE)
{
	
	# Color Vision Deficiency Transformation Matrices #
	
	cvdtransform = switch(cvd,
			
			# protanope: reds are greatly reduced (1% men, 0,02% women) #
			
			p = matrix(data = c(0, 2.02344, -2.52581, 0, 1, 0, 0, 0, 1),nrow = 3,ncol = 3,byrow = TRUE),
			
			# deuteranope: greens are greatly reduced (1% men, 0,01% women) #
			
			d = matrix(data = c(1, 0, 0, 0.494207, 0, 1.24827, 0, 0, 1),nrow = 3,ncol = 3,byrow = TRUE),
			
			# tritanope: blues are greatly reduced (0,002% men, 0,001% women) #
			
			t = matrix(data = c(1, 0, 0, 0, 1, 0, -0.395913, 0.801109, 0),nrow = 3,ncol = 3,byrow = TRUE)	
	)
	
	# color management ###
	
	rgbvec = colorpalette(colpal)
	rgbmat = col2rgb(rgbvec)
	
	# conversion of RGB coordinates into LMS, a color space suitable for calculating color blindness as it's represented by the three types of cones of the human eye, named after their sensitivity at wavelengths; Long (564-580nm), Medium (534-545nm) and Short (420-440nm) #
	
	rgb2lms = matrix(data = c(17.8824, 43.5161, 4.11935, 3.45565, 27.1554, 3.86714, 0.0299566, 0.184309, 1.46709),nrow = 3,ncol = 3,byrow = TRUE)
	lmsmat = rgb2lms%*%rgbmat
	
	# simulation of color blindness by reducing the colors along a dichromatic confusion line, the line parallel to the axis of the missing photoreceptor, to a single color #

	simmat = cvdtransform%*%lmsmat
	
	# conversion of LMS coordinates into RGB #
	
	lms2rgb = solve(rgb2lms)
	simrgbmat = lms2rgb%*%simmat

	# adjust for numerical instability #
	
	simrgbmat = apply(simrgbmat,c(1,2),function(x){if (x>255){return(255)} else if (x<0){return(0)} else return(x)})
	
	# conversion to hexadecimal color code #
			
	simpal = apply(round(simrgbmat),2,function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)})
	
	# compensation for color blindness by shifting wavelengths away from the portion of the spectrum invisible to the dichromat, towards the visible portion #

		# isolate invisible colors to color vision deficiency (calculate error matrix) #
		
		errormat = rgbmat - simrgbmat
		
		# shift colors towards visible spectrum (apply error modifications) #
		
		adjustmat = matrix(data = c(0, 0, 0, 0.7, 1, 0, 0.7, 0, 0),nrow = 3,ncol = 3,byrow = TRUE)
		compensationmat = adjustmat%*%errormat
		
		# add compensation to original values #
		
		enhancedrgbmat = compensationmat + rgbmat
	
	# adjust for numerical instability #
	
	enhancedrgbmat = apply(enhancedrgbmat,c(1,2),function(x){if (x>255){return(255)} else if (x<0){return(0)} else return(x)})
	
	# conversion to hexadecimal color code #
	
	enhancedpal = apply(enhancedrgbmat,2,function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)})
	
	# simulation of color blindness for the compensation values #
	
		# conversion of RGB coordinates into LMS #
		
		testlmsmat = rgb2lms%*%enhancedrgbmat
		
		# simulation of color blindness #
		
		testsimmat = cvdtransform%*%testlmsmat
		
		# conversion of LMS coordinates into RGB #
		
		testsimrgbmat = lms2rgb%*%testsimmat
		
		# adjust for numerical instability #
		
		testsimrgbmat = apply(testsimrgbmat,c(1,2),function(x){if (x>255){return(255)} else if (x<0){return(0)} else return(x)})
		
		# conversion to hexadecimal color code #
		
		testsimpal = apply(round(testsimrgbmat),2,function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)})
	
	# plot resulting colpals #
	
	if (show){
		nr = length(rgbvec)
		
		# The nonlinear relations for L, a, and b are intended to mimic the nonlinear response of the eye. Furthermore, uniform changes of components in the Lab color space aim to correspond to uniform changes in perceived color, so the relative perceptual differences between any two colors in Lab can be approximated by treating each color as a point in a three dimensional space and taking the Euclidean distance between them #
	
			# linearize RGB #
	
			linsimmat = apply(simrgbmat/255,c(1,2),function(x){if (x <= 0.04045)return(x/12.92) else return(((x + 0.055)/(1 + 0.055))^(2.4))})
					
			# conversion of RGB coordinates into XYZ (sRGB D65) ###
		
			rgb2xyz = matrix(data = c(0.4124564, 0.3575761, 0.1804375, 0.2126729, 0.7151522, 0.0721750, 0.0193339, 0.1191920, 0.9503041),nrow = 3,ncol = 3,byrow = TRUE)
			simxyzmat = rgb2xyz%*%linsimmat
			
			# tristimulus values of the reference white point: standard light source (D65) chosen according to 6500K color temperature (average daylight) #
		
			x = 0.312713
			y = 0.329016
			z = 1-x-y
			
			xn = x/y
			yn = 1
			zn = z/y
			
			# conversion of XYZ coordinates into Lab ###
			
			sqfkt = function(x){if (x > (6/29)^3)return(x^(1/3)) else return((1/3)*(29/6)^2*x + (4/29))}
			simLabmat = apply(simxyzmat,2,function(x){L = 116*sqfkt(x[2]/yn) - 16;a = 500*(sqfkt(x[1]/xn) - sqfkt(x[2]/yn));b = 200*(sqfkt(x[2]/yn) - sqfkt(x[3]/zn));return(c(L,a,b))})
				
			euclids = c()
			for (i in 1:nr){
				diffmat = simLabmat - simLabmat[,i]
				euclids = c(euclids,apply(diffmat[,setdiff(1:nr,i),drop = FALSE],2,function(x){return(sqrt(sum(x^2)))}))
			}
			euclids = unique(euclids)
		
			# Euclidean perception difference score #
	
			Epds = min(euclids)/sqrt(100^2 + (2*128)^2 + (2*128)^2)
			
		npal = 4
		plot(1,1,xlim=c(0,nr),ylim=c(0,npal*2),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
		rect(xleft = 0:(nr - 1),ybottom = 1 - 1,xright = 1:nr,ytop = 1 - 0.2,col = testsimpal,border = NA)
		rect(xleft = 0,ybottom = 1 - 1,xright = nr,ytop = 1 - 0.2,col = "transparent",border = "darkgrey")
		rect(xleft = 0:(nr - 1),ybottom = 3 - 1,xright = 1:nr,ytop = 3 - 0.2,col = enhancedpal,border = NA)
		rect(xleft = 0,ybottom = 3 - 1,xright = nr,ytop = 3 - 0.2,col = "transparent",border = "darkgrey")
		rect(xleft = 0:(nr - 1),ybottom = 5 - 1,xright = 1:nr,ytop = 5 - 0.2,col = simpal,border = NA)
		rect(xleft = 0,ybottom = 5 - 1,xright = nr,ytop = 5 - 0.2,col = "transparent",border = "darkgrey")
		rect(xleft = 0:(nr - 1),ybottom = 7 - 1,xright = 1:nr,ytop = 7 - 0.2,col = rgbvec,border = NA)
		rect(xleft = 0,ybottom = 7 - 1,xright = nr,ytop = 7 - 0.2,col = "transparent",border = "darkgrey")
		cvdtype = switch(cvd,p = "protanope",d = "deuteranope",t = "tritanope")
		text(rep(nr/2,npal),c(8,6,4,2)-0.6,labels = c(paste("Original colors"),paste("Simulated colors as seen by",cvdtype,"individuals"),paste("Enhanced colors for",cvdtype,"individuals"),paste("Enhanced simulated colors as seen by",cvdtype,"individuals")))
		mtext(paste("Euclidean perception difference score (Epds):",signif(Epds,1)),3,2,cex=1.25)
		mtext(paste("( should be >= 0.01 for sufficient perception differences )"),3,0,col="darkgrey")
	}
	
	# return colpals #

	return(list("simulated" = simpal,"enhanced" = enhancedpal))
}


# alias #


LSD.daltonize = daltonize


### colorpalette ###


#' @export
#' @name colorpalette
#' @aliases LSD.colorpalette
#' @title Provides colorpalettes containing R built-in colors
#' @description Provides pre-designed colorpalettes (character vectors containing R built-in colors) of this and several other R packages (grDevices, RColorBrewer, colorRamps) as well as custom-made ones.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}).
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to length of \code{colpal}, if not specified).
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is returned to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is returned to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param rev logical: if \code{TRUE} (\code{FALSE} by default), a reversed colorpalette is returned.
#' @return \code{colorpalette} returns a vector containing R built-in colors in hexadecimal representation.
#' @author Achim Tresch, Bjoern Schwalb
#' @seealso \code{\link{disco}}, \code{\link{demotour}}
#' @examples colorpalette("heat")
#' colorpalette(c("darkred","grey","darkblue"),10)
#' @keywords color, alpha


colorpalette = function(colpal,nrcol = NULL,simulate = FALSE,daltonize = FALSE,cvd = "p",alpha = NULL,rev = FALSE)
{
	if (length(colpal) > 1){palette = colpal}
	else{palette = switch(colpal,			
				
				# custom-made palettes #
				
				heat = c("grey","darkblue","red","orange","gold"),
				crazyred = c( "#940000","#A50000","#FF5C5C","#FFB9B9"),
				crazygreen = c("dark green","#009700","green","#C0F5D0"),
				crazyblue = c("dark blue","blue","#7390EE","light blue"),
				mountain = c("light green","dark green","black","dark grey","#F0F0F0"),
				girly = c("violet","violetred","violetred1","purple","purple3"),
				jamaica = c("red","yellow","green"),
				standard = c("brown","gold","yellow","lightyellow","white"),
				colorblind = c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"),
				
				# palettes from the RColorBrewer package #
				
				ylorrd = rev(c("#FFFFCC","#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#800026")),
				ylorbr = rev(c("#FFFFE5","#FFF7BC","#FEE391","#FEC44F","#FE9929","#EC7014","#CC4C02","#993404","#662506")),
				ylgnbu = rev(c("#FFFFD9","#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8","#253494","#081D58")),
				ylgn = rev(c("#FFFFE5","#F7FCB9","#D9F0A3","#ADDD8E","#78C679","#41AB5D","#238443","#006837","#004529")),
				reds = rev(c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")),
				rdpu = rev(c("#FFF7F3","#FDE0DD","#FCC5C0","#FA9FB5","#F768A1","#DD3497","#AE017E","#7A0177","#49006A")),
				purples = rev(c("#FCFBFD","#EFEDF5","#DADAEB","#BCBDDC","#9E9AC8","#807DBA","#6A51A3","#54278F","#3F007D")),
				purd = rev(c("#F7F4F9","#E7E1EF","#D4B9DA","#C994C7","#DF65B0","#E7298A","#CE1256","#980043","#67001F")),
				pubugn = rev(c("#FFF7FB","#ECE2F0","#D0D1E6","#A6BDDB","#67A9CF","#3690C0","#02818A","#016C59","#014636")),
				pubu = rev(c("#FFF7FB","#ECE7F2","#D0D1E6","#A6BDDB","#74A9CF","#3690C0","#0570B0","#045A8D","#023858")),
				orrd = rev(c("#FFF7EC","#FEE8C8","#FDD49E","#FDBB84","#FC8D59","#EF6548","#D7301F","#B30000","#7F0000")),
				oranges = rev(c("#FFF5EB","#FEE6CE","#FDD0A2","#FDAE6B","#FD8D3C","#F16913","#D94801","#A63603","#7F2704")),
				greys = rev(c("#FFFFFF","#F0F0F0","#D9D9D9","#BDBDBD","#969696","#737373","#525252","#252525","#000000")),
				greens = rev(c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")),
				gnbu = rev(c("#F7FCF0","#E0F3DB","#CCEBC5","#A8DDB5","#7BCCC4","#4EB3D3","#2B8CBE","#0868AC","#084081")),
				bupu = rev(c("#F7FCFD","#E0ECF4","#BFD3E6","#9EBCDA","#8C96C6","#8C6BB1","#88419D","#810F7C","#4D004B")),
				bugn = rev(c("#F7FCFD","#E5F5F9","#CCECE6","#99D8C9","#66C2A4","#41AE76","#238B45","#006D2C","#00441B")),
				blues = rev(c("#F7FBFF","#DEEBF7","#C6DBEF","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C","#08306B")),
				spectral = c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"),
				rdylgn = c("#A50026","#D73027","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#D9EF8B","#A6D96A","#66BD63","#1A9850","#006837"),
				rdylbu = c("#A50026","#D73027","#F46D43","#FDAE61","#FEE090","#FFFFBF","#E0F3F8","#ABD9E9","#74ADD1","#4575B4","#313695"),
				rdgy = c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#FFFFFF","#E0E0E0","#BABABA","#878787","#4D4D4D","#1A1A1A"),
				rdbu = c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061"),
				puor = c("#7F3B08","#B35806","#E08214","#FDB863","#FEE0B6","#F7F7F7","#D8DAEB","#B2ABD2","#8073AC","#542788","#2D004B"),
				prgn = c("#40004B","#762A83","#9970AB","#C2A5CF","#E7D4E8","#F7F7F7","#D9F0D3","#A6DBA0","#5AAE61","#1B7837","#00441B"),
				piyg = c("#8E0152","#C51B7D","#DE77AE","#F1B6DA","#FDE0EF","#F7F7F7","#E6F5D0","#B8E186","#7FBC41","#4D9221","#276419"),
				brbg = c("#543005","#8C510A","#BF812D","#DFC27D","#F6E8C3","#F5F5F5","#C7EAE5","#80CDC1","#35978F","#01665E","#003C30"),
				
				# palettes from the grDevices package #
				
				standardterrain = terrain.colors(9),
				standardtopo = topo.colors(9),
				standardheat = heat.colors(9),
				standardrainbow = rainbow(9,start=0.7,end=0.1),
				standardcm = cm.colors(9),
				
				# palettes from the colorRamps package #

				bl2gr = c("#0000FF","#001AE6","#0033CC","#004DB3","#006699","#008080","#009966","#00B34C","#00CC33","#00E619","#00FF00"),
				bl2gr2rd = c("#0000BF","#0000FF","#0080FF","#00FFFF","#40FFBF","#80FF80","#BFFF40","#FFFF00","#FF8000","#FF0000","#BF0000"),
				bl2rd = c("#0000FF","#0033FF","#0066FF","#0099FF","#00CCFF","#00FFFF","#FFCC00","#FF9900","#FF6600","#FF3300","#FF0000"),
				bl2yl = c("#0000FF","#1A1AE6","#3333CC","#4D4DB3","#666699","#808080","#999966","#B3B34C","#CCCC33","#E6E619","#FFFF00"),
				cy2yl = c("#00FFFF","#1AFFE6","#33FFCC","#4DFFB3","#66FF99","#80FF80","#99FF66","#B3FF4C","#CCFF33","#E6FF19","#FFFF00"),
				gr2rd = c("#00FF00","#1AE600","#33CC00","#4DB300","#669900","#808000","#996600","#B34C00","#CC3300","#E61900","#FF0000"),
				ma2gr = c("#FF00FF","#E61AE6","#CC33CC","#B34DB3","#996699","#808080","#669966","#4CB34C","#33CC33","#19E619","#00FF00"),
				matlablike = c("#0000AA","#0040FF","#0080FF","#40BFFF","#80FFFF","#BFFFBF","#FFFF80","#FFBF40","#FF8000","#FF4000","#AA0000"),
				matlablike2 = c("#0000BF","#0000FF","#0080FF","#00FFFF","#40FFBF","#80FF80","#BFFF40","#FFFF00","#FF8000","#FF0000","#BF0000") 
		)}
	if (is.null(palette)){stop("'colpal' should be a valid colorpalette name (see 'disco()') or a character vector of at least two R built-in color names for interpolation!")}
	if (is.null(nrcol)){nrcol = length(palette)}
	palette = try(colorRampPalette(palette)(nrcol),silent = TRUE)
	if (class(palette) == "try-error"){stop("'colpal' should be a valid colorpalette name (see 'disco()') or a character vector of at least two R built-in color names for interpolation!")}
	if (rev){palette = rev(palette)}
	if (simulate){palette = daltonize(palette,cvd = cvd,show = FALSE)$simulated}
	if (daltonize){palette = daltonize(palette,cvd = cvd,show = FALSE)$enhanced}
	if (!is.null(alpha)){palette = convertcolor(palette,alpha)}
	return(palette)
}


# alias #


LSD.colorpalette = colorpalette


### disco ###


#' @export
#' @name disco
#' @aliases LSD.disco
#' @aliases display.colorpalette
#' @aliases LSD.display.colorpalette
#' @title Disco (DISplays COlorpalettes)
#' @description Displays pre-designed colorpalettes as well as custom-made ones (see \code{colorpalette}).
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (displays all colorpalettes, if not specified).
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to length of \code{colpal}, if not specified).
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @author Bjoern Schwalb
#' @seealso \code{\link{colorpalette}}, \code{\link{demotour}}
#' @examples disco()
#' disco("rdbu",10)
#' @keywords color, alpha


disco = function(colpal = NULL,nrcol = NULL,alpha = NULL)
{
	if (is.null(colpal)){
		if (is.null(nrcol)){nrcol = 11}
		
		# plot custom-made palettes #
		
		ownpals = c("colorblind","standard","crazyred","crazygreen","crazyblue","girly","jamaica","mountain","heat")
		pals = ownpals
		npal = length(ownpals)
		plot(1,1,xlim=c(0,nrcol),ylim=c(0,npal),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
		for (i in 1:npal){rect(xleft = 0:(nrcol - 1),ybottom = i - 1,xright = 1:nrcol,ytop = i - 0.2,col = colorpalette(pals[i],nrcol),border = NA)}
		for (i in 1:npal){rect(xleft = 0,ybottom = i - 1,xright = nrcol,ytop = i - 0.2,col = "transparent",border = "darkgrey")}
		text(rep(-0.1,npal),(1:npal)-0.6,labels = pals,xpd = TRUE,adj = 1,cex=0.8)
		mtext(paste("palettes from the LSD package"),3,2,cex=1.25)
		mtext(paste("( character strings can be passed to LSD functions as 'colpal' )"),3,0,col="darkgrey")
		devAskNewPage(ask = TRUE)
		
		# plot palettes from the RColorBrewer package (part 1) #
		
		brewerpals = c("spectral","rdylgn","rdylbu","rdgy","rdbu","puor","prgn","piyg","brbg")
		pals = brewerpals
		npal = length(brewerpals)
		plot(1,1,xlim=c(0,nrcol),ylim=c(0,npal),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
		for (i in 1:npal){rect(xleft = 0:(nrcol - 1),ybottom = i - 1,xright = 1:nrcol,ytop = i - 0.2,col = colorpalette(pals[i],nrcol),border = NA)}
		for (i in 1:npal){rect(xleft = 0,ybottom = i - 1,xright = nrcol,ytop = i - 0.2,col = "transparent",border = "darkgrey")}
		text(rep(-0.1,npal),(1:npal)-0.6,labels = pals,xpd = TRUE,adj = 1,cex=0.8)
		mtext(paste("palettes from the RColorBrewer package"),3,2,cex=1.25)
		mtext(paste("( character strings can be passed to LSD functions as 'colpal' )"),3,0,col="darkgrey")
		devAskNewPage(ask = TRUE)
		
		# plot palettes from the RColorBrewer package (part 2) #
		
		brewerpals = c("ylorrd","ylorbr","ylgnbu","ylgn","reds","rdpu","purples","purd","pubugn","pubu","orrd","oranges","greys","greens","gnbu","bupu","bugn","blues")
		pals = brewerpals
		npal = length(brewerpals)
		plot(1,1,xlim=c(0,nrcol),ylim=c(0,npal),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
		for (i in 1:npal){rect(xleft = 0:(nrcol - 1),ybottom = i - 1,xright = 1:nrcol,ytop = i - 0.2,col = colorpalette(pals[i],nrcol),border = NA)}
		for (i in 1:npal){rect(xleft = 0,ybottom = i - 1,xright = nrcol,ytop = i - 0.2,col = "transparent",border = "darkgrey")}
		text(rep(-0.1,npal),(1:npal)-0.6,labels = pals,xpd = TRUE,adj = 1,cex=0.8)
		mtext(paste("palettes from the RColorBrewer package"),3,2,cex=1.25)
		mtext(paste("( character strings can be passed to LSD functions as 'colpal' )"),3,0,col="darkgrey")
		devAskNewPage(ask = TRUE)
		
		# plot palettes from the colorRamps package #

		rampspals = c("bl2gr","bl2gr2rd","bl2rd","bl2yl","cy2yl","gr2rd","ma2gr","matlablike","matlablike2")
		pals = rampspals
		npal = length(rampspals)
		plot(1,1,xlim=c(0,nrcol),ylim=c(0,npal),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
		for (i in 1:npal){rect(xleft = 0:(nrcol - 1),ybottom = i - 1,xright = 1:nrcol,ytop = i - 0.2,col = colorpalette(pals[i],nrcol),border = NA)}
		for (i in 1:npal){rect(xleft = 0,ybottom = i - 1,xright = nrcol,ytop = i - 0.2,col = "transparent",border = "darkgrey")}
		text(rep(-0.1,npal),(1:npal)-0.6,labels = pals,xpd = TRUE,adj = 1,cex=0.8)
		mtext(paste("palettes from the colorRamps package"),3,2,cex=1.25)
		mtext(paste("( character strings can be passed to LSD functions as 'colpal' )"),3,0,col="darkgrey")
		devAskNewPage(ask = TRUE)
		
		# plot palettes from the grDevices package #

		grpals = c("standardterrain","standardtopo","standardheat","standardrainbow","standardcm")
		pals = grpals
		npal = length(grpals)
		plot(1,1,xlim=c(0,nrcol),ylim=c(0,npal),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
		for (i in 1:npal){rect(xleft = 0:(nrcol - 1),ybottom = i - 1,xright = 1:nrcol,ytop = i - 0.2,col = colorpalette(pals[i],nrcol),border = NA)}
		for (i in 1:npal){rect(xleft = 0,ybottom = i - 1,xright = nrcol,ytop = i - 0.2,col = "transparent",border = "darkgrey")}
		text(rep(-0.1,npal),(1:npal)-0.6,labels = pals,xpd = TRUE,adj = 1,cex=0.8)
		mtext(paste("palettes from the grDevices package"),3,2,cex=1.25)
		mtext(paste("( character strings can be passed to LSD functions as 'colpal' )"),3,0,col="darkgrey")
	} else {
		
		# plot user-specified colorpalette if 'colpal != NULL' #
		
		if (is.null(nrcol)){nrcol = length(colorpalette(colpal,nrcol,alpha = alpha))}
		plot(1,1,xlim=c(0,nrcol),ylim=c(0,8),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
		rect(xleft = 0:(nrcol - 1),ybottom = 3,xright = 1:nrcol,ytop = 5,col = colorpalette(colpal,nrcol,alpha = alpha),border = NA)
		rect(xleft = 0,ybottom = 3,xright = nrcol,ytop = 5,col = "transparent",border = "darkgrey")
		if (length(colpal) > 1){mtext(paste("colorpalette"),3,2,cex=1.25)} else {mtext(paste(colpal,"colorpalette"),3,2,cex=1.25)}
		mtext(paste("( ",paste(colorpalette(colpal,nrcol,alpha = alpha),collapse = "  ")," )"),3,0,col="darkgrey",cex=1 + ((1-0.25)*(nchar(paste("( ",paste(colorpalette(colpal,nrcol,alpha = alpha),collapse = "  ")," )"))-50))/(50-185))
	}
}


# aliases #


LSD.disco = disco
display.colorpalette = disco
LSD.display.colorpalette = disco


### distinctcolors ###


#' @export
#' @name distinctcolors
#' @aliases LSD.distinctcolors
#' @title Find preferably distinct R built-in colors
#' @description Find a vector of distinct R built-in colors for a pre-defined length ('nrcol').
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to the length of \code{10}, if not specified).
#' @param method character string implying the method for color selection to be used ("RGB" uses a grid in the RGB space (default), "Lab" uses a grid in the Lab space, "goldenratio" uses the golden ratio as spacing between colors in the HSV color space).
#' @param bw logical: if \code{TRUE} (\code{FALSE} by default), the colors "black" and "white" are removed from the resulting colorpalette.
#' @param show logical: if \code{TRUE} (by default), the resulting colorpalettes are depicted in an R plot.
#' @param simulate logical: if \code{TRUE} (by default), a converted colorpalettes are additionally depicted to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @return \code{distinctcolors} returns a vector containing R built-in colors in hexadecimal representation.
#' @author Bjoern Schwalb
#' @seealso \code{\link{disco}}, \code{\link{colorpalette}}, \code{\link{demotour}}
#' @examples distinctcolors()
#' @keywords color


distinctcolors = function(nrcol = 10,method = "RGB",bw = FALSE,show = TRUE,simulate = TRUE)
{	
	# get.permutations function #
	
	get.permutations = function(n,r,v = 1:n){
		v = unique(sort(v))
		v0 = vector(mode(v),0)
		sub = function(n,r,v){
			if (r == 1){
				matrix(v,n,1)
			} else if (n == 1) {
				matrix(v,1,r)
			} else {
				inner = Recall(n,r - 1,v)
				cbind(rep(v,rep(nrow(inner),n)),matrix(t(inner),ncol = ncol(inner),nrow = nrow(inner)*n,byrow = TRUE))
			}
		}
		sub(n, r, v[1:n])
	}
	
	if (method == "RGB"){
		
		# 'RGB' uses a grid in the RGB space #
		
		# lay grid into the RGB space #
		
		lower = 1
		while ((nrcol+2) > lower^3) lower = lower + 1
		valseq = seq(0,255,length.out = lower)
        spl = get.permutations(length(valseq),3,valseq)
		
		# remove black and white #
		
		spl = spl[-c(1,nrow(spl)),]
		
		# sample nrcol of colors needed #
		
		spl = spl[sample(nrow(spl),nrcol),,drop=FALSE]
		
		# conversion to hexadecimal color code #
		
		colpal = apply(spl,1,function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)})
		
		# black and white? #
		
		if (!bw){colpal[1] = "#000000";colpal[length(colpal)] = "#FFFFFF"}
		
	} else if (method == "Lab"){
		
		# 'Lab' uses a grid in the Lab space #
		
		# lay grid into the Lab space #
		
		lower = 1
		while ((nrcol+2) > lower^3) lower = lower + 1
		valseq = seq(0,1,length.out = lower)
		spl = get.permutations(length(valseq),3,valseq)
		spl[,1] = spl[,1]*100
		spl[,2:3] = (spl[,2:3]*2 - 1)*150
		
		# The nonlinear relations for L, a, and b are intended to mimic the nonlinear response of the eye. Furthermore, uniform changes of components in the Lab color space aim to correspond to uniform changes in perceived color, so the relative perceptual differences between any two colors in Lab can be approximated by treating each color as a point in a three dimensional space and taking the Euclidean distance between them #
		
		# tristimulus values of the reference white point: standard light source (D65) chosen according to 6500K color temperature (average daylight) #
		
		x = 0.312713
		y = 0.329016
		z = 1-x-y
		
		xn = x/y
		yn = 1
		zn = z/y
		
		# conversion of Lab coordinates into XYZ ###
		
		sqfkt = function(x){if (x > 6/29)return(x^3) else return(3*(6/29)^2*(x - (4/29)))}
		xyzspl = apply(t(spl),2,function(x){X = xn*sqfkt((x[1] + 16)/116 + x[2]/500);Y = yn*sqfkt((x[1] + 16)/116);Z = zn*sqfkt((x[1] + 16)/116 - x[3]/200);return(c(X,Y,Z))})
		
		# conversion of XYZ coordinates into RGB (sRGB D65) ###
		
		rgb2xyz = matrix(data = c(0.4124564, 0.3575761, 0.1804375, 0.2126729, 0.7151522, 0.0721750, 0.0193339, 0.1191920, 0.9503041),nrow = 3,ncol = 3,byrow = TRUE)
		linrgbspl = solve(rgb2xyz)%*%xyzspl
		
		# adjust for numerical instability #
		
		linrgbspl = apply(linrgbspl,c(1,2),function(x){if (x>1){return(1)} else if (x<0){return(0)} else return(x)})
		
		# reverse linearization in RGB #
		
		revrgbspl = t(apply(linrgbspl,c(1,2),function(x){if (x <= 0.0031308)return(x*12.92) else return(x^(1/2.4)*(1 + 0.055) - 0.055)})*255)
		
		# remove black and white #
		
		revrgbspl = revrgbspl[-which(apply(round(revrgbspl),1,sum) %in% c(0,3*255)),]
		
		# sample nrcol of colors needed #
		
		revrgbspl = revrgbspl[sample(nrow(revrgbspl),nrcol),,drop=FALSE]
		
		# conversion to hexadecimal color code #
		
		colpal = apply(revrgbspl,1,function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)})
		
		# black and white? #
		
		if (!bw){colpal[1] = "#000000";colpal[length(colpal)] = "#FFFFFF"}
		
	} else if (method == "goldenratio"){
		
		# 'goldenratio' uses the golden ratio as spacing between colors in the HSV color space #
	
		# use random start value #
		
		h = runif(1)
		
		# subsequently add the golden ratio conjugate (1/golden ratio) modulo 1 #
		
		for (i in 1:(nrcol-1)){h = c(h,(h[length(h)] + 0.618033988749895)%%1)}
		
		# conversion to hexadecimal color code #
		
		colpal = hsv(h = h,s = 1,v = 1)
		
		# black and white? #
		
		if (!bw){colpal[1] = "#000000";colpal[length(colpal)] = "#FFFFFF"}
		
		# method switch trap #
		
	} else {stop("method need to be 'RGB', 'Lab' or 'goldenratio' !")}	
	
	# plot results? #
	
	if (show){
		if (simulate){
			nrcolpal = 4
			plot(1,1,xlim=c(0,nrcol),ylim=c(0,nrcolpal*2),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
			rect(xleft = 0:(nrcol - 1),ybottom = 1 - 1,xright = 1:nrcol,ytop = 1 - 0.2,col = daltonize(colpal,cvd = "t",show = FALSE)$simulated,border = NA)
			rect(xleft = 0,ybottom = 1 - 1,xright = nrcol,ytop = 1 - 0.2,col = "transparent",border = "darkgrey")
			rect(xleft = 0:(nrcol - 1),ybottom = 3 - 1,xright = 1:nrcol,ytop = 3 - 0.2,col = daltonize(colpal,cvd = "d",show = FALSE)$simulated,border = NA)
			rect(xleft = 0,ybottom = 3 - 1,xright = nrcol,ytop = 3 - 0.2,col = "transparent",border = "darkgrey")
			rect(xleft = 0:(nrcol - 1),ybottom = 5 - 1,xright = 1:nrcol,ytop = 5 - 0.2,col = daltonize(colpal,cvd = "p",show = FALSE)$simulated,border = NA)
			rect(xleft = 0,ybottom = 5 - 1,xright = nrcol,ytop = 5 - 0.2,col = "transparent",border = "darkgrey")
			rect(xleft = 0:(nrcol - 1),ybottom = 7 - 1,xright = 1:nrcol,ytop = 7 - 0.2,col = colpal,border = NA)
			rect(xleft = 0,ybottom = 7 - 1,xright = nrcol,ytop = 7 - 0.2,col = "transparent",border = "darkgrey")
			text(rep(nrcol/2,nrcolpal),c(8,6,4,2)-0.6,labels = c(paste("Original colors"),paste("Simulated colors as seen by protanope individuals"),paste("Simulated colors as seen by deuteranope individuals"),paste("Simulated colors as seen by tritanope individuals")))
			mtext(paste("distinctcolors"),3,2,cex=1.25)
			mtext(paste("( ",paste(colpal,collapse = "  ")," )"),3,0,col="darkgrey",cex=1 + ((1-0.25)*(nchar(paste("( ",paste(colpal,collapse = "  ")," )"))-50))/(50-185))
		} else {
			plot(1,1,xlim=c(0,nrcol),ylim=c(0,8),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
			rect(xleft = 0:(nrcol - 1),ybottom = 3,xright = 1:nrcol,ytop = 5,col = colpal,border = NA)
			rect(xleft = 0,ybottom = 3,xright = nrcol,ytop = 5,col = "transparent",border = "darkgrey")
			mtext(paste("distinctcolors"),3,2,cex=1.25)
			mtext(paste("( ",paste(colpal,collapse = "  ")," )"),3,0,col="darkgrey",cex=1 + ((1-0.25)*(nchar(paste("( ",paste(colpal,collapse = "  ")," )"))-50))/(50-185))
		}
	}
	
	# return colpal #
	
	return(colpal)
}


# alias #


LSD.distinctcolors = distinctcolors


### convertcolor ###


#' @export
#' @name convertcolor
#' @aliases LSD.convertcolor
#' @title Map R colors to hexadecimal representation
#' @description Convert R built-in colors to hexadecimal representation.
#' @param cols a character vector containing R built-in colors.
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @return \code{convertcolor} returns a vector containing R built-in colors in hexadecimal representation.
#' @author Bjoern Schwalb
#' @seealso \code{\link{disco}}, \code{\link{colorpalette}}, \code{\link{demotour}}
#' @examples convertcolor(c("red","green","blue"))
#' @keywords color


convertcolor = function(cols,alpha = NULL)
{
	colmat = col2rgb(cols)
	cols = apply(colmat,2,function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)})
	if (!is.null(alpha)){cols = paste(cols,alpha,sep="")}
	return(cols)
}


# alias #


LSD.convertcolor = convertcolor


### complementarycolor ###


#' @export
#' @name complementarycolor
#' @aliases LSD.complementarycolor
#' @title Complement R colors
#' @description Convert R built-in colors to their color complement
#' @param cols a character vector containing R built-in colors.
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @return \code{complementarycolor} returns a vector containing R built-in colors in hexadecimal representation.
#' @author Bjoern Schwalb
#' @seealso \code{\link{disco}}, \code{\link{colorpalette}}, \code{\link{demotour}}
#' @examples complementarycolor(c("red","green","blue"))
#' @keywords color


complementarycolor = function(cols,alpha = NULL)
{
	colmat = c(255,255,255) - col2rgb(cols)
	cols = apply(colmat,2,function(x){rgb(x[1],x[2],x[3],maxColorValue = 255)})
	if (!is.null(alpha)){cols = paste(cols,alpha,sep="")}
	return(cols)
}


# alias #


LSD.complementarycolor = complementarycolor


### convertgrey ###


#' @export
#' @name convertgrey
#' @aliases LSD.convertgrey
#' @title Convert R colors to greyscale
#' @description Greyscale R built-in colors.
#' @param cols a character vector containing R built-in colors.
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @return \code{convertgrey} returns a vector containing R built-in colors in hexadecimal representation.
#' @author Bjoern Schwalb
#' @seealso \code{\link{disco}}, \code{\link{colorpalette}}, \code{\link{demotour}}
#' @examples convertgrey(c("red","green","blue"))
#' @keywords greyscale, color


convertgrey = function(cols,alpha = NULL)
{
	colmat = col2rgb(cols)
	cols = apply(colmat,2,function(x){rgb(0.3*x[1]+0.59*x[2]+0.11*x[3],0.3*x[1]+0.59*x[2]+0.11*x[3],0.3*x[1]+0.59*x[2]+0.11*x[3],maxColorValue = 255)})
	if (!is.null(alpha)){cols = paste(cols,alpha,sep="")}
	return(cols)
}


# alias #


LSD.convertgrey = convertgrey



