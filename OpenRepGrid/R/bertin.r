
###############################################################################
###  							        BERTIN DISPLAYS 				                          ###
###############################################################################



constructCellGrob <- function(text, gp=gpar(), horiz=TRUE){
	gp <-  modifyList(gpar(fill=grey(.95)), gp)
	col <- gmSelectTextColorByLuminance(gp$fill)
	gTree(children=gList( rectGrob(width=1, height=1, 
							   		gp=gpar(fill=gp$fill, col="white")),
					  	  gmSplitTextGrob(text=text, horiz=horiz, gp=modifyList(gp, gpar(col=col)))
						 ))
}



bertin1 <- function(x, draw=TRUE){
	if(!inherits(x, "repgrid")) 
		stop("Object must be of class 'repgrid'")
	
  # determine color range (shades of grey)
  nrow <- nrow(x@ratings)
  ncol <- ncol(x@ratings)
  	
  # settings
	height.element.label <- 5
	height.cell <- unit(3, "mm")
	height.fg.top <- unit(ncol * height.element.label, "mm")
	

	bertinCell <- function(label, fill){
		textColor <- gmSelectTextColorByLuminance(fill)
		gTree(children=gList(
					  rectGrob(width=1, height=1, 
							   gp=gpar(fill=fill, col="white")),
					  textGrob(label=label, gp=gpar(lineheight=.7, cex=.6, col=textColor))
		))		
	}
	
	# rating framegrob
	dp.fg <- frameGrob(grid.layout(nrow=nrow, ncol=ncol, respect=F))
	scale.range <- x@scale$max - x@scale$min
	for (row in seq_len(nrow)){
		for (col in seq_len(ncol)){
			score <- x@ratings[row, col, 1]
			rg <- bertinCell(label=score, fill=grey((score - x@scale$min) / scale.range))
			dp.fg <- placeGrob(dp.fg , rg, row = row, col = col)
		}
	}
	
	# left framegrob (initial pole)
	left.c.fg <- frameGrob(grid.layout(nrow=nrow, ncol=1))
	for(row in seq_len(nrow)){
	  label <- x@constructs[[row]]$leftpole$name
		tg <- textGrob(label=label, gp=gpar(cex=.6))
		left.c.fg <- placeGrob(left.c.fg, tg, row=row)
	}
	
	# top framegrob (elements)
	top.e.fg <- frameGrob(grid.layout(ncol=ncol, nrow=ncol + 1, respect=F))
	rg <- rectGrob( gp=gpar(fill="black", col="white"), 
					        vp=viewport(width=unit(1, "points")))
	for(row in seq_len(ncol)){
	  label <- x@elements[[row]]$name
		tg <- textGrob(label=label, x=.4, just="left", gp=gpar(cex=.6))
		top.e.fg <- placeGrob(top.e.fg, tg, row=row, col=row)					
		top.e.fg <- placeGrob(top.e.fg, rg, row=row:ncol + 1, col=row)			
	}
	
	# combine framegrobs
	main.fg <- frameGrob(grid.layout(nrow=4, ncol=3, heights=c(.1, 2 ,2,.2), widths=c(1,2,1)))
	main.fg <- placeGrob(main.fg, top.e.fg, row= 2, col = 2)	
	main.fg <- placeGrob(main.fg, left.c.fg, row= 3, col = 1)
	main.fg <- placeGrob(main.fg, dp.fg, row= 3, col = 2)
	main.fg <- placeGrob(main.fg, left.c.fg, row=3, col = 3)
	if(draw) grid.draw(main.fg) else main.fg
}





bertin2 <- function(x, ratings=TRUE, top=unit(40, "mm"), sides=unit(40, "mm") , 
					left=sides, right=sides, 
					cell=unit(6, "mm"), cell.height=cell, cell.width=cell, 
					gp.cells=gpar(), gp.constructs=gpar(), gp.elements=gpar(), 
					bg.col=grey(.95), colors=c("white", "black"), draw=TRUE){
	if(!inherits(x, "repgrid")) 
		stop("Object must be of class 'repgrid'")
	
	gp.cells <- modifyList(gpar(lineheight=.7, cex=.6, fill=bg.col), gp.cells) 
	gp.constructs <- modifyList(gpar(lineheight=.7, cex=.8, fill=bg.col), gp.constructs) 
	gp.elements <- modifyList(gpar(lineheight=.7, cex=.8, fill=bg.col), gp.elements) 
	
	# determine color range (shades of grey)
	nrow <- nrow(x@ratings)
	ncol <- ncol(x@ratings)
	
	height.top <- top
	width.left <- left
	width.right <- right
	height.cell <- cell.height
	width.cell <- cell.width
	height.body <- 	nrow * height.cell
	width.body <- ncol * width.cell
	
	bertinCell <- function(label, fill, gp=gpar(), ratings=TRUE){
		textColor <- gmSelectTextColorByLuminance(fill)
		gp <- modifyList(gp, gpar(col=textColor))
		if(ratings) tg <- textGrob(label=label, gp=gp) else tg <- nullGrob()
		gTree(children=gList(
					  rectGrob(width=1, height=1, 
							   gp=gpar(fill=fill, col="white")),
					  tg
		))		
	}
	
	# rating framegrob
	colorFun <- makeStandardRangeColorRamp(colors)
	dp.fg <- frameGrob(grid.layout(nrow=nrow, ncol=ncol, respect=F))
	scale.range <- x@scale$max - x@scale$min
	scale.min <- x@scale$min
	for (row in seq_len(nrow)){
		for (col in seq_len(ncol)){
			score <- x@ratings[row, col, 1]
			rg <- bertinCell(label=score, fill=colorFun((score-scale.min)/scale.range), gp=gp.cells, ratings=ratings)
			dp.fg <- placeGrob(dp.fg , rg, row = row, col = col)
		}
	}
	
	# left framegrob (initial pole)
	left.c.fg <- frameGrob(grid.layout(nrow=nrow, ncol=1))
	for(row in seq_len(nrow)){
	  text <- x@constructs[[row]]$leftpole$name
		tg <- constructCellGrob(text=text, gp=gp.constructs)
		left.c.fg <- placeGrob(left.c.fg, tg, row=row)	
	}
	
	# right framegrob (contrast pole)
	right.c.fg <- frameGrob(grid.layout(nrow=nrow, ncol=1))
	for(row in seq_len(nrow)){
	  text <- x@constructs[[row]]$rightpole$name
		tg <- constructCellGrob(text=text, gp=gp.constructs)
		right.c.fg <- placeGrob(right.c.fg, tg, row=row)	
	}
	
	# top framegrob (elements)
	top.e.fg <- frameGrob(grid.layout(ncol=ncol, nrow=1))
	for(col in seq_len(ncol)){
	  text <- x@elements[[col]]$name
		tg <- constructCellGrob(text=text, horiz=FALSE, gp=gp.elements)
		top.e.fg <- placeGrob(top.e.fg, tg, row=NULL, col=col)					
	}
	
	# combine framegrobs
	main.fg <- frameGrob(grid.layout(nrow=2, ncol=3, heights=unit.c(height.top, height.body), widths=unit.c(width.left, width.body, width.right)))
	main.fg <- placeGrob(main.fg, top.e.fg, row= 1, col = 2)	
	main.fg <- placeGrob(main.fg, left.c.fg, row= 2, col = 1)
	main.fg <- placeGrob(main.fg, dp.fg, row= 2, col = 2)
	main.fg <- placeGrob(main.fg, right.c.fg, row=2, col = 3)
	if(draw) grid.draw(main.fg) else main.fg
}



bertin2PlusLegend <- function(x, ratings=TRUE, top=unit(40, "mm"), 
								sides=unit(40, "mm"), left=sides, right=sides, 
								cell=unit(6, "mm"), cell.height=cell, cell.width=cell, 
								gp.cells=gpar(), gp.constructs=gpar(), gp.elements=gpar(), 
								bg.col=grey(.95), colors=c("white", "black"), draw=TRUE,
								vspace=unit(2,"mm"), legend.just="left", legend.height=unit(10, "mm"),
								legend.width=unit(40, "mm"))
{
		fg.bertin <- bertin2(	x=x, ratings=ratings, top=top, 
								sides=sides, left=left, right=right, 
								cell=cell, cell.height=cell.height, cell.width=cell.width, 
								gp.cells=gp.cells, gp.constructs=gp.constructs, gp.elements=gp.elements, 
								bg.col=bg.col, colors=colors, draw=FALSE)
		
		widths <- fg.bertin$framevp$layout$widths
		heights <- fg.bertin$framevp$layout$heights
		nrow <- fg.bertin$framevp$layout$nrow
		ncol <- fg.bertin$framevp$layout$ncol

		colorFun <- makeStandardRangeColorRamp(colors)	
		lg <- gmLegend2(colorFun(c(0,1)), c("left pole", "right pole"), ncol=2, byrow=F)
		fg.legend <- frameGrob(grid.layout(widths=legend.width, just=legend.just))
		fg.legend <- placeGrob(fg.legend, lg)
		fg.main <- frameGrob(grid.layout(nrow=nrow + 2, heights=unit.c(heights, vspace, legend.height),
										 ncol=ncol, widths=widths))
		fg.main <- placeGrob(fg.main, fg.bertin, row=1:nrow)
		fg.main <- placeGrob(fg.main, fg.legend , row=nrow + 2)
		
		if(draw) grid.draw(fg.main)	else fg.main
}

# bertin2PlusLegend(rg2, colors=c("darkred", "white"))
# bertin2PlusLegend(rg2, colors=c("darkred", "white"), top=unit(4, "cm"), sides=unit(4, "cm"))




# TODO: -may work with closures here to store old row and column when marking 
#        rows and columns?
#       -splitString has a bug, breaks too late
#       -trimming of elements and constructs
#
#' Workhorse for the biplot printing. 
#'
#' Prints a bertin to the output 
#' device. It uses the R base graphics system and 
#' this is very fast. This is useful for working with grids. Not so much for
#' producing high-quality output.
#'
#' @param x         \code{repgrid} object. 
#' @param ratings   Vector. rating scores are printed in the cells
#' @param margins   Vector of length three (default \code{margins=c(0,1,1)}). 
#'                  1st element denotes the left, 2nd the upper and 3rd the 
#'                  right margin in npc coordinates (i.e. 0 to zero).
#' @param trim      Vector (default \code{trim=c(F,F)}).If a number the string
#'                  is trimmed to the given number of characters. If set 
#'                  to TRUE the labels are trimmed to the available space
#' @param add       Logical. Wether to add bertin to existent plot (default is 
#'                  \code{FALSE}). If \code{TRUE, plot.new()} will not be called
#'                  \code{par(new=TRUE)}.
#' @return \code{NULL} just for printing.
#'
#' @export
#' @keywords internal
#' @author Mark Heckmann
#' 
bertinBase <- function(nrow, ncol, labels="", labels.elements="", 
                       labels.left="", labels.right="", 
                       col.text=NA, cex.text=.6, cex.elements=.7, 
                       cex.constructs=.7, col.fill=grey(.8), border="white", 
                       xlim=c(0,1), ylim=c(0,1), margins=c(0,1,1), lheight=.75,
                       text.margin=0.005, elements.offset=c(0.002, 0.002), 
                       id=c(T,T), cc=0, cr=0, cc.old=0, cr.old=0, 
                       col.mark.fill="#FCF5A4", print=TRUE, byrow=FALSE, add=FALSE)
{
  if (byrow)
    labels <- as.vector(matrix(labels, nrow=nrow, ncol=ncol, byrow=TRUE))
  col.fill <- recycle(col.fill, nrow*ncol)    # recycle col.fill if too short e.g. one color
  if (identical(col.text, NA))                # if not explicitly defined replace col.text according to bg color
    col.text <- gmSelectTextColorByLuminance(col.fill)
  else recycle(col.text, nrow*ncol)
  #if (length(trim) == 1)    # if only one parameter given, extend to the other
  #   trim <- recycle(trim, 2)
  if (length(id) == 1)
    id <- recycle(id, 2)    

  makeMain <- function(){
    rect(x1, y1, x2, y2, col = col.fill, border = border)
    text(x1 + cell.width/2, y1 + cell.height/2, labels=labels, col=col.text, cex=cex.text)
  }
  
  makeElements <- function(){       #### elements
    index <- cascade(ncol, type=2)
    if (id[2]){
      labels.elements[index$left] <- paste(labels.elements[index$left], 
                                           "-", index$left)
      labels.elements[index$right] <- paste(index$right, "-", 
                                            labels.elements[index$right])
    }

    height.strokes <- (margins[2] - ylim[2]) / (max(cascade(ncol) + 1))
    x.lines <- xlim[1] + x1.o * diff(xlim) + cell.width / 2
    y1.lines <- ylim[2]
    y2.lines <- y1.lines + cascade(ncol) * height.strokes   # upper end of bertin main plus offset
    segments(x.lines, y1.lines, x.lines, y2.lines)
    text(x.lines[index$left] + elements.offset[1], 
        y2.lines[index$left] + elements.offset[2], 
        labels=labels.elements[index$left], adj=c(1,0), cex=cex.elements, xpd=T)
    text(x.lines[index$right] - elements.offset[1], 
         y2.lines[index$right] + elements.offset[2], 
         labels=labels.elements[index$right], adj=c(0,0), cex=cex.elements, xpd=T)
  }
  
  makeConstructs <- function(){     ### constructs
    if (id[1]){
      labels.left <- paste(labels.left, " (", 1:nrow, ")", sep="")
      labels.right <- paste("(", 1:nrow, ") ", labels.right, sep="")
    }  
    labels.left <- baseSplitString(labels.left, availwidth= (xlim[1] - margins[1])* .95, cex=cex.text)
    labels.right <- baseSplitString(labels.right, availwidth=(margins[3] - xlim[2]) * .95, cex=cex.text)
    par(lheight=lheight)    # set lineheight
    text(xlim[1] - text.margin, y1[1:nrow] + cell.height/2, labels=labels.left, 
         cex=cex.constructs, adj=1, xpd=T)
    text(xlim[2] + text.margin, y1[1:nrow] + cell.height/2, labels=labels.right, 
         cex=cex.constructs, adj=0, xpd=T)
  }
  
  colorRow <- function(cr){
    par(new=TRUE)   # next plot will overplot not earse the old one, necessary for setting the same regions
    plot.new()
    #plot.window(xlim=0:1, ylim=0:1) #, xaxs="i", yaxs="i")#, asp =nrow/ncol)
    if (cr >= 1 & cr <= nrow){      # color current row cr
      labels.rows <- labels[(1:ncol-1)*nrow + cr]
      col.mark.text=gmSelectTextColorByLuminance(col.mark.fill)
      rect(x1.rc, y1.rc[cr], x2.rc, y2.rc[cr], 
           col = col.mark.fill, border = border)
      text(x1.rc + cell.width/2, y1.rc[cr] + cell.height/2, 
           labels=labels.rows, col=col.mark.text, cex=cex.text)
    }
  }
  
  colorColumn <- function(cc){
    par(new=TRUE)   # next plot will overplot not earse the old one, necessary for setting the same regions
    plot.new()
    #plot.window(xlim=0:1, ylim=0:1) #, xaxs="i", yaxs="i")#, asp =nrow/ncol)
    if (cc >= 1 & cc <= ncol){      # color current column cc
      labels.cols <- labels[1:nrow + (cc-1)*nrow]
      #col.fill <- col.fill[1:nrow + (cc-1)*nrow]
      #col.text=gmSelectTextColorByLuminance(col.fill)
      col.mark.text=gmSelectTextColorByLuminance(col.mark.fill)
      rect(x1.rc[cc], y1.rc, x2.rc[cc], y2.rc, 
           col = col.mark.fill, border = border)
      text(x1.rc[cc] + cell.width/2, y1.rc + cell.height/2, 
           labels=labels.cols, col=col.mark.text, cex=cex.text)  
      # color vertical stroke
      height.strokes <- (1 - ylim[2]) / (max(cascade(ncol) + 1))
      x.lines <- xlim[1] + x1.o * diff(xlim) + cell.width / 2
      y1.lines <- ylim[2]
      y2.lines <- y1.lines + cascade(ncol) * height.strokes  
      segments(x.lines[cc], y1.lines, x.lines[cc], y2.lines[cc], lwd=3, col="white")  # overplot old stroke in white
      segments(x.lines[cc], y1.lines, x.lines[cc], y2.lines[cc], col=col.mark.fill) 
    }
  }
  
  renewColumn <- function(cc){
    if (cc >= 1 & cc <= ncol){
      # vertical stroke
      height.strokes <- (1 - ylim[2]) / (max(cascade(ncol) + 1))
      x.lines <- xlim[1] + x1.o * diff(xlim) + cell.width / 2
      y1.lines <- ylim[2]
      y2.lines <- y1.lines + cascade(ncol) * height.strokes  
      segments(x.lines[cc], y1.lines, x.lines[cc], y2.lines[cc], lwd=3, col="white")  # overplot old stroke in white
      segments(x.lines[cc], y1.lines, x.lines[cc], y2.lines[cc], col="black") 
      
      # plot rects and text
      labels.cols <- labels[1:nrow + (cc-1)*nrow]
      col.fill <- col.fill[1:nrow + (cc-1)*nrow]
      col.text=gmSelectTextColorByLuminance(col.fill)
      rect(x1.rc[cc], y1.rc, x2.rc[cc], y2.rc, 
           col = col.fill, border = border)
      text(x1.rc[cc] + cell.width/2, y1.rc + cell.height/2, 
           labels=labels.cols, col=col.text, cex=cex.text)
    }
  }
  
  renewRow <- function(cr){
    if (cr >= 1 & cr <= nrow){
      # plot rects and text
      labels.rows <- labels[(1:ncol-1)*nrow + cr]
      col.fill <- col.fill[(1:ncol-1)*nrow + cr]
      col.text=gmSelectTextColorByLuminance(col.fill)
      rect(x1.rc, y1.rc[cr], x2.rc, y2.rc[cr], 
           col = col.fill, border = border)
      text(x1.rc + cell.width/2, y1.rc[cr] + cell.height/2, 
           labels=labels.rows, col=col.text, cex=cex.text)
    }
  }
  
  # make basic calculations
  x1.o <- 0:(ncol - 1)/ncol
  x2.o <- 1:ncol/ncol
  y1.o <- rev(0:(nrow - 1)/nrow)
  y2.o <- rev(1:nrow/nrow)
  
  x1 <- rep(x1.o, each=nrow)
  x2 <- rep(x2.o, each=nrow)
  y1 <- rep(y1.o, ncol)
  y2 <- rep(y2.o, ncol)

  x1 <- xlim[1] + x1 * diff(xlim)      # rescale coordinates according to given limits
  x2 <- xlim[1] + x2 * diff(xlim) 
  y1 <- ylim[1] + y1 * diff(ylim)
  y2 <- ylim[1] + y2 * diff(ylim)
  
  cell.width <- diff(xlim) / ncol
  cell.height <- diff(ylim) / nrow
  
  x1.rc <- x1[(1:ncol)*nrow]      # calc coords for row and col starts and ends
  x2.rc <- x2[(1:ncol)*nrow]
  y1.rc <- y1[1:nrow]
  y2.rc <- y2[(1:nrow)]
  
  # set plotting parameters
  #old.par <- par(no.readonly = TRUE)    # save parameters
  #on.exit(par(old.par))                 # reset old par when done
  op <- par(oma=rep(0,4), mar=rep(0,4), xaxs="i", yaxs="i")
  if (print)                  # in case no new printing should occur
    par(new=FALSE)
  else 
    par(new=TRUE)
  if (add)                    # will bertin be added to existent plot?
    par(new=TRUE)
  
  plot.new()
  #plot.window(xlim=0:1, ylim=0:1) #, xaxs="i", yaxs="i")#, asp =nrow/ncol)
  
  # plotting
  if (print) {
    makeMain()
    makeElements()
    makeConstructs()
    colorRow(cr)        # color current row or column
    colorColumn(cc)
  } else {
    renewColumn(cc.old)
    renewRow(cr.old)
    colorRow(cr)
    colorColumn(cc)
  }     
  #par(op)
  invisible(NULL)
}


#bertinBase(20, 70, xlim=c(.2,.8), ylim=c(0,.4))
#bertinBase(10,20)
#bertinBase(10,20, xlim=c(0.1, .9), ylim=c(.2, .8), cex.text=.8)
#bertinBase(20, 30, grey(runif(13)), cex.text=.6)
#labels <- randomSentences(20, 6)
#bertinBase(20, 70, xlim=c(.25,.75), ylim=c(.1,.4), margins=c(.03,.9,.97), id=F, 
#           labels.l=labels, labels.ri=labels, labels.el=rep(labels, 4))
           
#x <- randomGrid(20, 40)
#nc <- length(x@constructs)
#ne <- length(x@elements)
#color <- c("darkred", "white", "darkgreen")
#colorFun <- makeStandardRangeColorRamp(color)
#scale.min <- x@scale$min
#scale.max <- x@scale$max
#scores <- as.vector(x@ratings[,,1])
#col.fill <- colorFun((scores-scale.min)/(scale.max-scale.min))
#bertinBase(nc, ne, col.fill, scores , xlim=c(.2, .8), ylim=c(0,.6), cex.text=.6, border="white")
#bertinBase(nc, ne, col.fill, scores , xlim=c(.2, .8), ylim=c(0,.6), cex.text=.6, border="white", cc=10, cr=10, pri=F)
#bertinBase(nc, ne, col.fill, scores , xlim=c(.2, .8), ylim=c(0,.6), cex.text=.6, border="white", cc.old=10, pri=F)
#bertinBase(nc, ne, col.fill, scores , xlim=c(.2, .8), ylim=c(0,.6), cex.text=.6, border="white", cr.old=10, pri=F)

#bertinBase(nc, ne, col.fill, scores , xlim=c(.2, .8), ylim=c(0,.6), cex.text=.6, border="white")
#for (row in 1:10){
#  for (col in 1:15) {
#    bertinBase(nc, ne, col.fill, scores , xlim=c(.2, .8), ylim=c(0,.6), cex.text=.6, 
#               border="white", cc=col, cr=row, cc.old=col -1, cr.old=row-1, pri=F)
#    Sys.sleep(.2)
#  }
#}


#' Make Bertin display of grid data.
#'
#' One of the most popular ways of displaying grid data has been adopted 
#' from Bertin's (1974) graphical proposals, which have had an immense 
#' influence onto data visualization. One of the most appealing 
#' ideas presented by Bertin is the concept of the reordable matrix. 
#' It is comprised of graphical displays for each cell, allowing to 
#' identify structures by eye-balling reordered versions of the data matrix 
#' (see Bertin, 1974). In the context of repertory grids, 
#' the display is made up of a simple colored rectangle where the color 
#' denotes the corresponding score. Bright values correspond to low, dark 
#' to high scores. For an example of how to analyze a Bertin display see 
#' e.g. Dick (2000) and Raeithel (1998).
#'
#' @param x               \code{repgrid} object.
#' @param colors          Vector. Two or more colors defininig the color ramp for 
#'                        the bertin (default \code{c("white", "black")}).
#' @param showvalues      Logical. Wether scores are shown in bertin
#' @param xlim            Vector. Left and right limits inner bertin (default 
#'                        \code{c(.2, .8)}).
#' @param ylim            Vector. Lower and upper limits of inner bertin
#'                        default(\code{c(.0, .6)}).
#' @param margins         Vector of length three (default \code{margins=c(0,1,1)}). 
#'                        1st element denotes the left, 2nd the upper and 3rd the 
#'                        right margin in npc coordinates (i.e. 0 to zero).
#' @param cex.elements    Numeric. Text size of element labels (default \code{.7}).
#' @param cex.constructs  Numeric. Text size of construct labels (default \code{.7}).
#' @param cex.text        Numeric. Text size of scores in bertin cells (default \code{.7}).
#' @param col.text        Color of scores in bertin (default \code{NA}). By default
#'                        the color of the text is chosen according to the 
#'                        background color. If the background ist bright the text 
#'                        will be black and vice versa. When a color is specified
#'                        the color is set independetn of background.
#' @param border          Border color of the bertin cells (default \code{white}).
#' @param lheight         Line height for constructs.
#' @param id              Logical. Wheteher to print id number for constructs and elements
#'                        respectively (default \code{c(T,T)}).
#' @param cc              Numeric. Current column to mark.
#' @param cr              Numeric. Current row to mark.
#' @param cc.old          Numeric. Column to unmark.
#' @param cr.old          Numeric. Row to unmark.
#' @param col.mark.fill   Color of marked row or column (default \code{"#FCF5A4"}).
#' @param print           Print whole bertin. If \code{FALSE} only current and old
#'                        row and column are printed.
#' @param ...             Optional arguments to be passed on to \code{bertinBase}.
#'
#' @return \code{NULL} just for the side effects, i.e. printing.
#'
#' @export
#' @references    Bertin, J. (1974). \emph{Graphische Semiologie: Diagramme, Netze, 
#'                  Karten}. Berlin, New York: de Gruyter.
#'
#'                Dick, M. (2000). The Use of Narrative Grid Interviews in 
#'                  Psychological Mobility Research. \emph{Forum Qualitative 
#'                  Sozialforschung / Forum: Qualitative Social Research, 1}(2). 
#'                
#'                Raeithel, A. (1998). Kooperative Modellproduktion von 
#'                  Professionellen und Klienten - erlauetert am Beispiel des 
#'                  Repertory Grid. \emph{Selbstorganisation, Kooperation, Zeichenprozess: 
#'                  Arbeiten zu einer kulturwissenschaftlichen, anwendungsbezogenen 
#'                  Psychologie} (pp. 209-254). Opladen: Westdeutscher Verlag.
#'
#' @examples \dontrun{
#' 
#'    bertin(feixas2004)
#'    bertin(feixas2004, c("white", "darkblue"))
#'    bertin(feixas2004, showvalues=F)
#'    bertin(feixas2004, border="grey")
#'    bertin(feixas2004, cex.text=.9)
#'    bertin(feixas2004, id=c(F, F))
#'    
#'    bertin(feixas2004, cc=3, cr=4)
#'    bertin(feixas2004, cc=3, cr=4, col.mark.fill="#e6e6e6")
#' }
#'                       
bertin <- function(x, colors=c("white", "black"), showvalues=TRUE, 
                   xlim=c(.2, .8), ylim=c(0,.6), margins=c(0,1,1),
                   cex.elements=.7, cex.constructs=.7, cex.text=.6, col.text=NA, 
                   border="white", lheight=.75, id=c(T,T),
                   cc=0, cr=0, cc.old=0, cr.old=0, col.mark.fill="#FCF5A4", print=TRUE, 
                   ...){
  if (!inherits(x, "repgrid")) 							      # check if x is repgrid object
  	stop("Object x must be of class 'repgrid'")
  	
  nc <- length(x@constructs)
  ne <- length(x@elements)
  colorFun <- makeStandardRangeColorRamp(colors)
  scale.min <- x@scale$min
  scale.max <- x@scale$max
  scores <- as.vector(x@ratings[,,1])
  scores.standardized <- (scores-scale.min)/(scale.max-scale.min) 
  col.fill <- colorFun(scores.standardized)
  if (!showvalues)
    scores <- recycle("", nc * ne)
  bertinBase(nrow=nc, ncol=ne, labels=scores, labels.elements=getElementNames(x),
             labels.left=getConstructNames(x)$leftpole, 
             labels.right=getConstructNames(x)$rightpole,
             col.fill=col.fill,
             xlim=xlim, ylim=ylim, margins=margins,
             cex.elements=cex.elements, cex.constructs=cex.elements,
             cex.text=cex.text, col.text=col.text, 
             border=border, lheight=lheight, id=id, cc=cc, cr=cr, cc.old=cc.old, cr.old=cr.old,
             col.mark.fill=col.mark.fill, print=print, ...)
  invisible(NULL)
}

#x <- randomGrid(10,20)
#x





#' Bertin display with corresponding cluster anaylsis. 
#'
#' Element columns and 
#' constructs rows are ordered according to cluster criterion. Various 
#' distance measures as well as cluster methods are supported.
#'
#' @param x           \code{repgrid} object.
#' @param  dmethod    The distance measure to be used. This must be one of 
#'                    \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, 
#'                    \code{"canberra"}, \code{"binary"}, or \code{"minkowski"}.
#'                    Default is \code{"euclidean"}.
#'                    Any unambiguous substring can be given (e.g. \code{"euc"} 
#'                    for \code{"euclidean"}). 
#'                    A vector of length two can be passed if a different distance measure for
#'                    constructs and elements is wanted (e.g.\code{c("euclidean", "manhattan")}).
#'                    This will apply euclidean distance to the constructs and
#'                    manhattan distance to the elements.
#'                    For additional information on the different types see
#'                    \code{?dist}. 
#' @param  cmethod    The agglomeration method to be used. This should be (an
#'                    unambiguous abbreviation of) one of \code{"ward"}, 
#'                    \code{"single"}, \code{"complete"}, \code{"average"}, 
#'                    \code{"mcquitty"}, \code{"median"} or \code{"centroid"}.
#'                    Default is \code{"ward"}.
#'                    A vector of length two can be passed if a different cluster method for
#'                    constructs and elements is wanted (e.g.\code{c("ward", "euclidean")}).
#'                    This will apply ward clustering to the constructs and
#'                    single linkage clustering to the elements. If only one of either
#'                    constructs or elements is to be clustered the value \code{NA}
#'                    can be supplied. E.g. to cluster elements only use \code{c(NA, "ward")}.
#' @param  p          The power of the Minkowski distance, in case \code{"minkowski"}
#'                    is used as argument for \code{dmethod}. \code{p} can be a vector
#'                    of length two if different powers are wanted for constructs and
#'                    elements respectively (e.g. \code{c(2,1)}).
#' @param align       Whether the constructs should be aligned before clustering
#'                    (default is \code{TRUE}). If not, the grid matrix is clustered 
#'                    as is. See Details section in function \code{\link{cluster}} for more information.
#' @param trim        The number of characters a construct is trimmed to (default is
#'                    \code{10}). If \code{NA} no trimming is done. Trimming
#'                    simply saves space when displaying the output.
#' @param type        Type of dendrogram. Either or \code{"triangle"} (default) 
#'                    or \code{"rectangle"} form.
#' @param xsegs       Numeric vector of normal device coordinates (ndc i.e. 0 to 1) to mark
#'                    the widths of the regions for the left labels, for the
#'                    bertin display, for the right labels and for the 
#'                    vertical dendrogram (i.e. for the constructs).
#' @param ysegs       Numeric vector of normal device coordinates (ndc i.e. 0 to 1) to mark
#'                    the heights of the regions for the horizontal dendrogram 
#'                    (i.e. for the elements), for the bertin display and for 
#'                    the element names.
#' @param x.off       Horizontal offset between construct labels and construct dendrogram and 
#                     between the outer right margin and the dendrogram 
#'                    (default is \code{0.01} in normal device coordinates).
#' @param y.off       Vertical offset between bertin display and element dendrogram and 
#                     between the lower margin and the dendrogram
#'                    (default is \code{0.01} in normal device coordinates).
#' @param cex.axis    \code{cex} for axis labels, default is \code{.6}.
#' @param col.axis    Color for axis and axis labels, default is \code{grey(.4)}.
#' @param draw.axis   Whether to draw axis showing the distance metric for the dendrograms 
#'                    (default is \code{TRUE}).
#' @param ...         additional parameters to be passed to function \code{\link{bertin}}.
#'
#' @return            A list of two \code{\link{hclust}} object, for elements and constructs
#'                    respectively.
#'
#' @author        Mark Heckmann
#' @export
#' @seealso  \code{\link{cluster}}
#'
#' @examples \dontrun{
#'
#'    # default is euclidean distance and ward clustering 
#'    bertinCluster(bell2010)                                     
#'
#'    ### applying different distance measures and cluster methods
#'
#'    # euclidean distance and single linkage clustering 
#'    bertinCluster(bell2010, cmethod="single")
#'    # manhattan distance and single linkage clustering             
#'    bertinCluster(bell2010, dmethod="manhattan", cm="single") 
#'    # minkowksi distance with power of 2 = euclidean distance  
#'    bertinCluster(bell2010, dm="mink", p=2)                     
#' 
#'    ### using different methods for constructs and elements
#'
#'    # ward clustering for constructs, single linkage for elements
#'    bertinCluster(bell2010, cmethod=c("ward", "single"))        
#'    # euclidean distance measure for constructs, manhatten 
#'    # distance for elements
#'    bertinCluster(bell2010, dmethod=c("euclidean", "man"))
#'    # minkowski metric with different powers for constructs and elements    
#'    bertinCluster(bell2010, dmethod="mink", p=c(2,1)))          
#'
#'    ### clustering either constructs or elements only
#'    # euclidean distance and ward clustering for constructs no 
#'    # clustering for elements
#'    bertinCluster(bell2010, cmethod=c("ward", NA))  
#'    # euclidean distance and single linkage clustering for elements 
#'    # no clustering for constructs            
#'    bertinCluster(bell2010, cm=c(NA, "single"))                 
#'
#'    ### changing the appearance
#'    # different dendrogram type        
#'    bertinCluster(bell2010, type="rectangle")  
#'    # no axis drawn for dendrogram                 
#'    bertinCluster(bell2010, draw.axis=F)                        
#'
#'    ### passing on arguments to bertin function via ...
#'     # grey cell borders in bertin display
#'    bertinCluster(bell2010, border="grey")  
#'    # omit printing of grid scores, i.e. colors only                  
#'    bertinCluster(bell2010, showvalues=FALSE)                   
#'
#'    ### changing the layout
#'    # making the vertical dendrogram bigger
#'    bertinCluster(bell2010, xsegs=c(0, .2, .5, .7, 1))
#'    # making the horizontal dendrogram bigger          
#'    bertinCluster(bell2010, ysegs=c(0, .3, .8, 1))              
#' }
#'
bertinCluster <- function(x, dmethod=c("euclidean", "euclidean"), 
                          cmethod=c("ward", "ward"), p=c(2,2), align=TRUE, 
                          trim=NA, type=c("triangle"), 
                          xsegs = c(0, .2, .7, .9, 1), ysegs = c(0, .1, .7, 1),
                          x.off=0.01, y.off=0.01,
                          cex.axis =.6, col.axis =  grey(.4), draw.axis=TRUE, ...)
{
  if (length(dmethod) == 1)       # if only one value is passed
    dmethod <- rep(dmethod, 2)
  if (length(cmethod) == 1)       # if only one value is passed
    cmethod <- rep(cmethod, 2)
  if (length(p) == 1)             # if only one value is passed
    p <- rep(p, 2)

  cex.dend <- 0.001       # size text dendrogram, only needed for sanity 
                          # check purposes, otherwise 0.001 so no dend labels are drawn
                                                
  inr.x <- xsegs[4]       # inner figure region (bertin) ndc x coordinate range
                          # range goes from left side to y dendrogram region                                            
  inr.y <- 1 - ysegs[2]   # bertin fig region range as ndc coords
                          # range goes from end of x dendrogram region to end of device (i.e. 1)              
  
  # transform xsegs and ysegs coordinates (ndc) into 
  # ndc coordinates for inner figure region used by bertin plot
  xlim.bertin <- xsegs[2:3] / inr.x
  ylim.bertin <- c(0, (ysegs[3] - ysegs[2]) / inr.y)
  
  if (align)               # align grid if promoted, uses dmethod etc. for constructs, i.e. [1]
    x <- align(x, along = 0, dmethod = dmethod[1], 
               cmethod = cmethod[1], p = p[1])  
    
  r <- getRatingLayer(x, trim=trim)    # get ratings

  # dendrogram for constructs
  if (is.na(cmethod[1])){
    con.ord <- seq_len(getNoOfConstructs(x))          # no change in order
    fit.constructs <- NULL
  } else {
    dc <- dist(r, method = dmethod[1], p=p[1])        # make distance matrix for constructs
    fit.constructs <- hclust(dc, method=cmethod[1])   # hclust object for constructs
    dend.con <- as.dendrogram(fit.constructs)
    con.ord <- order.dendrogram(rev(dend.con))
  }
  
  # dendrogram for elements
  if (is.na(cmethod[2])){
    el.ord <- seq_len(getNoOfConstructs(x))          # no change in order
    fit.elements <- NULL
  } else {
    de <- dist(t(r), method = dmethod[2], p=p[2])     # make distance matrix for elements
    fit.elements <- hclust(de, method=cmethod[2])     # hclust object for elements
    dend.el <- as.dendrogram(fit.elements)
    el.ord <- order.dendrogram(dend.el)
  }
  
  x <- x[con.ord, el.ord]   # reorder repgrid object
  
  plot.new()
  par(fig = c(xsegs[c(1,4)], ysegs[c(2,4)]) , new=T)
  #par(fig = c(0, .8, .2, 1), new=T)
  
  bertin(x, xlim=xlim.bertin, ylim=ylim.bertin, add=T, ...)           # print reordered bertin
  
  # x dendrogram (horizontal) elements
  if (!is.na(cmethod[2])){
    dend.x.fig <- c(xsegs[2:3], ysegs[1:2]) + c(0,0, y.off, -y.off)     # adjust for offsets
    par(fig = dend.x.fig, new=T, mar=c(0,0,0,0))
    ymax.el <- attr(dend.el, "height")
    plot(dend.el, horiz=F, xlab="", xaxs="i", yaxs="i", yaxt="n",
          nodePar=list(cex=0, lab.cex=cex.dend), ylim=c(ymax.el,0), type=type)
    if (draw.axis)              # wether to draw axis
      axis(2, las=1, cex.axis=cex.axis, col=col.axis, col.axis=col.axis)
  }
  
  # y dendrogram (vertical) constructs
  if (!is.na(cmethod[1])){
    dend.y.fig <- c(xsegs[4:5], ysegs[2:3]) + c(x.off, -x.off, 0, 0)    # adjust for offsets           
    par(fig = dend.y.fig, new=T, mar=c(0,0,0,0))
    xmax.con <- attr(dend.con, "height")
    plot(dend.con, horiz=T, xlab="", xaxs="i", yaxs="i", yaxt="n",
          nodePar=list(cex=0, lab.cex=cex.dend), xlim=c(0,xmax.con), type=type)
    if (draw.axis)              # wether to draw axis
      axis(1, las=1, cex.axis=cex.axis, col=col.axis, col.axis= col.axis)
  }
  # return hclust objects for elements and constructs
  invisible(list(constructs=fit.constructs, elements=fit.elements))
}

# TODO: use of layout does not work with bertinCluster
# a future version could use layout
# layout (matrix(1:4), 2)
# bertinCluster(bell2010)

# bertinCluster(bell2010, type="t", bor=grey(.5))
# dev.new()
# bertinCluster(bell2010, type="t", dm="manhattan", cm="single")
# dev.new()
# bertinCluster(bell2010, type="t", dm="manhattan", cm="centroid")



# base graphics for quick clustering
# # make dendrogram
#x <- bell2010

# compute new layout of bertin










