#
# Release Version 130602 build - V1.0 (micromapST)
#


panelFill <-
function(col="#D0D0D0",border=NA,...)
# fill a panel specified by "usr" with a fill color of "col=" and a outline color of "border="
{
    xy <- par("usr")               # get usr data (x&y scaling) - points on the panel box
    #  draw polygon (box) with color "border" and fill "col"
    polygon(xy[c(1, 2, 2, 1)], xy[c(3, 3, 4, 4)],col=col,border=border,xpd=TRUE,...)
}


panelGrid <-
function(x = NULL, y = NULL, col = 2, lwd = 1, lty = 1)
#
# defaults = col = red,  lwd = 1 pt,  lty = solid
# place grids in panel.  If x present = vertical grids
#                        if y present = horizontal grids  
#                        if x and y present = both grids
{
     if(!is.null(x))
		abline(v = x, lwd = lwd, lty = lty, col = col)
     if(!is.null(y))
		abline(h = y, lwd = lwd, lty = lty, col = col)
}


panelInbounds <-
function(bnds)
#  bnds = min and max of panel boundaries.
#  times potential pretty to panel limits.
{
    grid = pretty(bnds)
    return(grid[bnds[1] < grid & grid < bnds[2]])
}


panelLengthen <-
function(x, n=1)
#   x = original vector to lengthen 
#   n = number of entries in new vector. 
#  expand number of columns or rows to "n", but replicated data may not be right..  WATCH!.
#    data in original vector is to be replicated to new elements. 
#  If no original data, zeroes are provided.
{
   if(n<1)stop("panelLengthen: invalid required length")       
   if(length(x)==0) return(rep(0,n))    # original vector length = 0,  return vector with "n" zeros.
   newx = rep(x,ceiling(n/length(x)))   # repeat x into new space. (multiple length of original)
   length(newx) = n                     # set results to length of n  (trim to "n" length)
   return(newx)
}


panelOutline <-
function(col = "black", lwd = 1, lty = 1)
# Outline panel in "col". Current panel = "usr", col="black", lwd = 1, lty = solid
{
	xy <- par("usr")            # get window size (save to reuse)
	polygon(xy[c(1, 2, 2, 1)], xy[c(3, 3, 4, 4)], density=0, col=col, xpd=TRUE)
}

panelScale <-
function(rx = c(0, 1), ry = c(0, 1),firstp=FALSE, inches = FALSE)
#  Set scale of panel.
#  If inches - set rx and ry to current inches values and return.
#  Do "new=TRUE" plot to set the scale.  (could use plot.window?)
{
	if(inches) {
	   pin = par("pin")
	   rx = c(0, pin[1])
	   ry = c(0, pin[2])
	}
	warn = unlist(options('warn'))
        options(warn=-1)
        
        par(new=TRUE)
        options(warn=warn)
	plot(rx, ry, type = "n", axes = FALSE, xaxs = "i", yaxs = "i", 
		xlab = "", ylab = "", main = "")
	return(list(rx = rx, ry = ry))
}


panelSelect <-
function(layout, i = 1, j = 1, margin = NULL)
#
# Panel Select
#    Layout = panel structure
#       dim = dimensions of panel = c(i,j).  If i or j > respective dimension - fatal.
#      If no margin specified---
#       datfig = par(fig <- data)
#       pad = par(mai = layout$pad[c(4,1,3,2)] # reorder
#      if margin specified---
#       labfig = par(fig = layout$labfig[ind,]) # based on margin 
#       brd = par(mai = layout$brd[c(4,1,3,2)] # reorder 
#         par(fig = c(a,b,c,d)...  NDC coordinates of figure region in the display.  
#                                  new plot, unless new=TRUE?
#         par(mai = c(b,l,t,r)...  numerical vector margin size (bot, left, top, right) in inches.
#    i =     column index
#    j =     row index
#    margin = left, right, top, bottom, bot,...
#

{
	dim = layout$dim
	if(i > dim[1] || j > dim[2])
	    stop("PS-01 Dimension error. Program error - index 'i' or 'j' is out of bounds.")
	
	if(is.null(margin)) {               # "margin" is missing.
            k = dim[2] * (i - 1) + j          # #col * (rowInx-1) + colInx 
                     # datfig layout as C1R1, C2R1, C3R1, C4R1, C2R1, ...
	    par(fig = layout$datfig[k,  ], 
	        mai = layout$pad[c(4, 1, 3, 2)] )  # pad is c(L(1), R(2), T(3), B(4))
	               # 
	           # fig is c(x1, x2, y1, y2)     # no units.
	           # mai is c(bot, left, top, right) <- pad[c(4,1,3,2)] in inches
	           
	}
	else {
	    vec = 1:4
	    nam = c("left", "right", "top", "bottom", "bot")
	    ind = match(margin, nam)
	    if(is.na(ind))
		stop("Bad label region name")
	    if(ind == 5)   ind = 4      # "bot" -> "bottom"
	    par(fig = layout$labfig[ind,  ], 
	        mai = layout$brd[c(4, 1, 3, 2)] )
	}
#	"done"
}


panelLayout <-
function(nrow = 1, 
      ncol = 1, 
      leftMargin = 0,                         # inches
      rightMargin = 0,                        # inches
      topMargin = 1,                          # inches, leave room for page titles. 
      bottomMargin = 0,                       # inches
      borders = rep(0.5, 4),                  # inches 
      # The figure borders are left(1), right, top, bottom(4)
      colSize = 1, 
      rowSize = 1, 
      colSep = 0, 
      rowSep = 0, 
      pad = NULL)
{
  # Note fig matrices rounded to 6 places in an effort of avoid a R problem with fig when
  #  values appear in exponential notation.

	oldpar = par()                       # save original par values.
	din = oldpar$din                     # get device dimensions (inches)
	din.x = din[1]                       #  x = width
	din.y = din[2]                       #  y = height

	plotX = din.x - borders[1] - borders[2] - leftMargin - rightMargin  # usable width inches
	plotY = din.y - borders[3] - borders[4] - bottomMargin - topMargin  # usable height inches
	
	# bounds (x1, x2, y1, y2)
	
	#   bounds (edge left, margin left, margin right, edge right) shifted right by "borders[1]"
	xbnds = c(0, leftMargin, leftMargin + plotX, leftMargin + plotX + 
		rightMargin) + borders[1]  #shift all by the left border
	
	#   bounds (edge bottom, margin bottom, margin top, edge top) shifted up by "borders[4]"
	ybnds = c(0, bottomMargin, bottomMargin + plotY, bottomMargin + 
		plotY   + topMargin) + borders[4] # shift all by bottom border
	
	# the right and top borders are handled in the first calculation.
	
	#  fig.scale = inch coordinates of device space.
	fig.scale = c(din.x, din.x, din.y, din.y)	

        # left figure is in the left margin space of the plot area from top to bottom
	leftfig = c(xbnds[1] - borders[1], xbnds[2] + borders[2], ybnds[1] - 
		borders[4], ybnds[4] + borders[3])
	
	# right figure is in the right margin space of the plot area from top to bottom
	rightfig = c(xbnds[3] - borders[1], xbnds[4] + borders[2], ybnds[1] - 
		borders[4], ybnds[4] + borders[3])
	
	# top figure is in the top margin space from from left to right
	topfig = c(xbnds[1] - borders[1], xbnds[4] + borders[2], ybnds[3] - 
		borders[4], ybnds[4] + borders[3])
	
	# bottom figure is in the bottom margin space from left to right
	botfig = c(xbnds[1] - borders[1], xbnds[4] + borders[2], ybnds[1] - 
		borders[4], ybnds[2] + borders[3])
	
	# these figure areas are from the devices left to right and top to bottom limits.
	# 
	
	colSep = panelLengthen(colSep, ncol + 1)      # initially 1 element of zero - now "ncol" elements
	rowSep = panelLengthen(rowSep, nrow + 1)      # nrow elements.
	
	if(is.null(pad)) 
	  {  # no pad, initialize
	     pad = c(borders[1] + colSep[1] + leftMargin, 
		     borders[2] + colSep[ncol + 1] + rightMargin, 
		     borders[3] + rowSep[1] + topMargin, 
		     borders[4] + rowSep[nrow + 1] + bottomMargin)	
	  }
	# The borders should align around the edge.
	
	colSep = cumsum(colSep)             # convert individual spaces to cumulative sums.           
	rowSep = cumsum(rowSep)
	
	plotX = plotX - colSep[ncol + 1]    # subtract room for separators (value of last colSep or rowSep
	plotY = plotY - rowSep[nrow + 1]
	
	# the colSize and rowSize values are relative to a projected sum of units.
	
	# sum(colSize) = ncol (number of columns)
	# sum(rowSize) = 71.64 as coded for 10 full groups and the median group
	
	relx = panelLengthen(colSize, ncol)  # size required == 3 * 1    relx is 3 elements of 1
	rely = panelLengthen(rowSize, nrow)  # size required == 10 * 7 + 1.64  
	relx = relx/sum(relx)                #   Factional each element / sum
	rely = rely/sum(rely)
	
	xinc = plotX * cumsum(c(0, relx))
	yinc = plotY * cumsum(c(0, rely))
	
	fig = matrix(0, nrow = nrow * ncol, ncol = 4)
	k = 0
	for(i in 1:nrow) {
           for(j in 1:ncol) {
              k = k + 1
	      fig[k, 1] = xbnds[2] + xinc[j] + colSep[j] - pad[1]             # x.pos k<=1:- start
	      fig[k, 2] = xbnds[2] + xinc[j + 1] + colSep[j] + pad[2]         # x pos      - end
	      fig[k, 4] = ybnds[3] - yinc[i] - rowSep[i] + pad[3]             # y pos - start
	      fig[k, 3] = ybnds[3] - yinc[i + 1] - rowSep[i] - pad[4]         # y pos - end 
	   }
	}
	
	# fig now has the x1, x2, y1, y2 physical position of each panel on the page.
	
	labfig = rbind(leftfig, rightfig, topfig, botfig)
	# lab figure has four rows - one for each margin/border space around the plot area. 
	
	#  Scale the "inch" coordinates to coordinates 0 to 1 (relative)
	fig    = abs(t(t(fig)/fig.scale))
	labfig = t(t(labfig)/fig.scale)	
	
	# coltabs are in inches and start inside the left border
	coltabs = cbind(c(0, colSep + xinc + leftMargin), leftMargin + c(0,colSep) + c(xinc, 
	                  xinc[length(xinc)] + rightMargin))	
	
	# rowtabs are in inches and start below the upper border
        
        # yinc is padded with a leading 0.  
	rowtabs = cbind(c(ybnds[3], ybnds[3] - rowSep - c(yinc[-1], yinc[nrow + 1] + bottomMargin)), 
	                c(ybnds[4], ybnds[3] - rowSep - yinc)) - borders[4]
	
	# The tabs provide the physical points for each panel.
	
	list(dim = c(nrow, ncol), datfig = round(fig,6), labfig = round(labfig,6), brd = borders,
		pad = pad, coltabs = coltabs, rowtabs = rowtabs, 
		figsize = c(din.x, din.y))
}

