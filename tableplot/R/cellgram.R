# mod MF: use TRUE/FALSE, not T/F
# mod MF: changed t.col to label.col for consistency

cellgram <- function(

	#-- Arguments that may be vectorized:

	cell,  		       #-- cell value(s) 
	
	shape=0,		   		#-- 1. shape of cell value(s); 0=circle, 1=diamond, 2=square, 3=cross
	shape.col="black", 	#-- 2. color of shape(s), outline only
	shape.lty=1,	   		#-- 3. line type used for shape(s)
	
	shape.neg=0,		   	#-- 4. See 1.
	shape.col.neg="red",	#-- 5. 
	shape.lty.neg=1,		#-- 6. 

	#-- Arguments that can not be vectorized:
	
	cell.fill="white", 	#-- 7. inside color of |smallest| shape in a cell
	back.fill="white", 	#-- 8. background color of cell 
	
	label=0,           	#-- 9. how many cell values will be printed; max is 4
	label.size=0.7, 		#-- 10.		
#	t.col="black",     	#-- 11. color of cell label(s)
	label.col="black",  #-- 11. color of cell label(s)
	
	ref.lines=FALSE,   	#-- 12. to draw ref lines or not
	ref.col="grey80",  	#-- 13. color for ref lines
	
	scale.max=1,			  #-- 14.
	
	shape.lwd=1,	   		#-- line width of shape(s)
	frame.col="black", 	#-- color of frame around cell
	frame.lwd="0.5"    	#-- line width of frame around cell
	)
	
	{

	## A function that sorts a vector according to the absolute value of its elements:
	absort <- function(y,decreasing=TRUE){
		y[sort(abs(y),decreasing=decreasing,index.return=TRUE)$i]
		}
		
	grid.rect(gp=gpar(fill=back.fill, col=frame.col, lwd=frame.lwd))

	if (length(cell)>length(shape))     shape=rep(shape, length(cell))
	if (length(cell)>length(shape.col)) shape.col=rep(shape.col, length(cell))
	if (length(cell)>length(shape.lty)) shape.lty=rep(shape.lty, length(cell))

	if (length(cell)>length(shape.neg))     shape.neg=rep(shape.neg, length(cell))
	if (length(cell)>length(shape.col.neg)) shape.col.neg=rep(shape.col.neg, length(cell))
	if (length(cell)>length(shape.lty.neg)) shape.lty.neg=rep(shape.lty.neg, length(cell))
	
	#-- Draw reference lines:

	if (ref.lines==TRUE) {
		grid.segments(x0=0,y0=.5,x1=1,y1=.5, gp=gpar(col=ref.col, lwd=0.5))
		grid.segments(x0=.5,y0=0,x1=.5,y1=1, gp=gpar(col=ref.col, lwd=0.5))
		}

	#-- Rescale cell values:

	s.cell = cell / scale.max

	#-- Draw cell values:

	k.to.fill <- which.min(abs(cell)) #-- Which shape to fill?
	
	for (k in 1:length(cell)){
		
		if (!is.na(cell[k])) { ## This is to allow missing values; but if all missing, then error ensues!

			if (cell[k] < 0) this.col   <- shape.col.neg[k] else this.col <- shape.col[k]
			if (cell[k] < 0) this.lty   <- shape.lty.neg[k] else this.lty <- shape.lty[k]
			if (cell[k] < 0) this.shape <- shape.neg[k] else this.shape <- shape[k]
			
			if (k==k.to.fill) fill <- cell.fill else fill <- "transparent"
			
			if (this.shape==0 || this.shape=="circle") 
				grid.circle( 
					r=abs(s.cell[k]/2), 
					gp=gpar(col=this.col, lty=this.lty, lwd=shape.lwd, fill=fill))

			if (this.shape==1 || this.shape=="diamond") {
				r1 = 0.5 - 0.5*abs(s.cell[k])
				r2 = 0.5*abs(s.cell[k]) + 0.5
				grid.polygon( 
					x=c(r1, .5, r2, .5), y=c(.5, r2, .5, r1), 
					gp=gpar(col=this.col, lty=this.lty, lwd=shape.lwd, fill=fill)) }

			if (this.shape==2 || this.shape=="square") 
				grid.rect(
					height=abs(s.cell[k]), width=abs(s.cell[k]), 
					gp=gpar(col=this.col, lty=this.lty, lwd=shape.lwd, fill=fill))

			if (this.shape==3 || this.shape=="cross"){
				s3 = (1 - s.cell[k])/2 # can't think of any other name
				grid.segments(x0=s3, y0=.5, x1=1-s3, y1=.5, gp=gpar(col=this.col, lwd=shape.lwd))   
				grid.segments(x0=.5, y0=s3, x1=.5, y1=1-s3, gp=gpar(col=this.col, lwd=shape.lwd)) }  		
		}}

	# For black and white plots, do the following:

	#if (min(cell[!is.na(cell)]) < 0) cell.fill = "white"

	#-- Labels
	
	if (label > 0){
		cell = absort(cell,decreasing=TRUE) # omit step if making best fit residual plots; see above for "absort"
		d = matrix(c(1,1,0,0,0,1,1,0),4,2)
		for (k in 1:min(label,4,length(cell))){
			grid.text(cell[k], gp=gpar(cex=label.size, col=label.col),
			    x=unit(c(.97,.97,.03,.03)[k],"npc"), 
			    y=unit(c(.03,.97,.97,.03)[k],"npc"), just=d[k,])
			}
		}
	}