# mod MF: use TRUE/FALSE, not T/F
# mod MF: changed t.col to label.col for consistency
# modified 1-14-2010 MF: set default for assign.sets
# modified 1-16-2010 MF: replaced ugly call to cellgram with do.call; cleaned up code
# modified 1-19-2010 MF: title default=NULL, better defaults for top.space, left.space
# modified 1-25-2010 MF: added grid.rect() for text.m display

# DONE: Calculate better default values for top.space, left.space
# TODO: Perhaps should use mar = c(left, top, right, bottom) instead to allow space on all margins
# TODO: There should be some default for cell.specs:  formals(cellgram)[-1], but scale.max more sensible
# TODO: Allow better positioning of title (which should probably be called 'main')
# TODO: Rename args for easier understanding: label.size -> label.cex, empty.text.size -> text.cex

tableplot <-
	function(values,  ...) UseMethod("tableplot")

tableplot.default <- function(

	values, 			# Array of values to plot; can also be a matrix or data frame
	assign.sets,			# Matrix of specification assignments 
	cell.specs, 			# List of lists; each list is one specification
	
	v.parts	=0, 	# List of column clusters
	h.parts	=0, 	# List of row clusters
	gap		=2, 	# Width of partitions (in millimeters)
	
	text.m 	=0,		# Matrix of text for insertion into empty cell(s)
	empty.text.size 	= 0.8,
	empty.text.col 	= "grey30", 
	
	title		=NULL, 

	table.label=TRUE,	
	label.size =0.8,	
	side.rot	=0,		# Degree of rotation (positive for counter-clockwise)
	
	left.space	=10,		# Millimeters between left of tableplot and left edge of drawing region
	top.space	=10+10*(!is.null(title)),	# Analogous to above
	...                 # Required for consistency with generic
	){

	require(grid)
#	require(lattice)    # does not require lattice

	rows <- dim(values)[1]
	cols <- dim(values)[2]
	# default for assign.sets
	if (missing(assign.sets)) assign.sets <- matrix(1, rows, cols)
	n.sets <- max(sets <- unique(as.vector(assign.sets)))
	
	## A function to construct a gap list:
	gap.list <- function(partitions=0,x){
		if (length(partitions)==1) rep(0,x) else {
		rep(1:length(partitions), partitions)-1}
		}
		
	grid.newpage()

	#---Constructing vectors of gaps if partitions provided.

	v.gaps <- gap.list(partitions=v.parts, x=cols)
	h.gaps <- gap.list(partitions=h.parts, x=rows)

	#---Constructing labels, if no row or column names in values.

	if (length(rownames(values)) == rows) side.label <- rownames(values) 
	else side.label <- (1:rows)
	
	if (length(colnames(values)) == cols) top.label <- colnames(values) 
	else top.label <- (1:cols)
	
	#---Add on extra dimension to values if values only has two dimensions.
	#---This step cannot occur before previous step.
	
	if (length(dim(values))==2) {dim(values) <- c(rows, cols, 1)}
	
	#---Create a text.m matrix if none provided.
	#---This would be filled with NAs; useful if empty cells requested but do not insert text
	
	if (!is.matrix(text.m)) text.m <- matrix("", rows, cols)
	
	#---Create Layout 1 and write main title.

	L1 <- grid.layout(2,1,heights=unit(c(top.space,1),c("mm","null"))) # Layout has 2 rows and 1 column
	
	pushViewport(viewport(layout=L1, width=1, height=1, x=unit(0,"npc"), y=unit(1,"npc"), just=c("left","top")))
	
	#---Push row 1 of Layout 1.
	
	pushViewport(viewport(layout.pos.row=1)) 	
	if (!is.null(title)) grid.text(title, x=0.02, just=c("left", "bottom"))
	upViewport()

	#---Create Layout 2.

	L2 <- grid.layout(1,2,widths=unit(c(left.space,1),c("mm","null")))
	
	#---Push row 2 of Layout 1.
	
	pushViewport(viewport(layout.pos.row=2, layout=L2)) 

	#---Create Layout 3.

	L3 <- grid.layout(rows,cols, respect=TRUE, just=c("left","top"))
	
	#---Push col 2 of Layout 2.
	
	pushViewport(viewport(layout.pos.col=2))		
	#---Push Layout 3, but with adjustments to accomodate possible partitions.

	pushViewport(viewport(layout=L3,  
				    x=0, y=1, just=c(0,1), 
				    width =unit(1,"npc")-unit(gap,"mm")*(length(v.parts)-1),
				    height=unit(1,"npc")-unit(gap,"mm")*(length(h.parts)-1)))
	
	#---Draw cellgrams.

	for (i in 1:rows){
		for (j in 1:cols){

			pushViewport(viewport(layout.pos.row=i, layout.pos.col=j))
			pushViewport(viewport(just=c(0,1), height=1, width=1, 
							x=unit(gap,"mm")*v.gaps[j], 
							y=unit(1,"npc")-unit(gap,"mm")*h.gaps[i]))				

			if (assign.sets[i,j]==0) {
				grid.text(text.m[i,j],gp=gpar(cex=empty.text.size, col=empty.text.col))
				grid.rect(gp=gpar(col="black", lwd=1))
			}
			else 
			{
			spec <- cell.specs[[assign.sets[i,j]]]
#browser()
            do.call(cellgram, c(list(cell=values[i,j,]), spec))
#			cellgram(cell  		= values[i,j,], 
#				   shape 	 		= spec[[1]],
#				   shape.col 		= spec[[2]],
#				   shape.lty 		= spec[[3]],
#				   shape.neg		= spec[[4]],
#				   shape.col.neg = spec[[5]],
#				   shape.lty.neg = spec[[6]],
#				   cell.fill 		= spec[[7]],
#				   back.fill 		= spec[[8]],
#				   label	 		= spec[[9]],
#				   label.size		= spec[[10]],
#				   label.col			= spec[[11]],
#				   ref.lines		= spec[[12]],
#				   ref.col	 		= spec[[13]],
#				   scale.max 		= spec[[14]]
#				   ) 
			}

			##grid.rect()
			if ((j==1) && (table.label==TRUE)) {
					if (side.rot==0) {grid.text(side.label[i], x=-0.14, just=1, gp=gpar(cex=label.size))}
					else 			{grid.text(side.label[i], x=-0.35, just=c("center"), rot=side.rot, gp=gpar(cex=label.size))}
					}
			if ((i==1) && (table.label==TRUE)) {grid.text(top.label[j],  y=1.15, vjust=0, gp=gpar(cex=label.size))}
			upViewport()
			upViewport()
			
			}
		}
	popViewport(0)
	}
