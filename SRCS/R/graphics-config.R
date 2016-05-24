## Returns a tagged list with the arguments to be passed to graphics::layout(...) to set the general layout
## of the figure
## Arguments:
## out.X: tagged list with the parameters of the outer X (horizontal, on top of the figure) variable and labels.
## out.Y: tagged list with the parameters of the outer Y variable and labels
## inner.X.par: tagged list with the parameters of the outer X variable and labels
## inner.Y.par: tagged list with the parameters of the outer Y variable and labels
## outerColumns: number of different levels of the outer X variable
## outerRows: number of different levels of the outer Y variable
## heatmaps.per.row: number of levels that should be displayed in each row of heatmaps. Should be smaller than outerColumns.
## show.colorbar: logical, whether a legend colorbar should be displayed on the right of every row of the figure
## longest.inner.Y.lev.width: numeric, number of inches that the longest string from the inner X levels would occupy. Used here to create an additional dummy blank
##           column so that the inner X levels of the left-most column of heatmaps, which are drawn outside the cell margin, have enough space 
## highest.inner.X.level: the same for the height (in inches)
## annotate: boolean, whether annotations should be included in the top left corner of the plot or not
# ----------------------------------------------------------------------------
.compose.layout <- function(out.X, out.Y, inner.X.par, inner.Y.par, outerColumns, outerRows, heatmaps.per.row, show.colorbar, 
                            longest.inner.Y.lev.width, highest.inner.X.level, annotate){
	
	left.cols = 0;		# additional columns on the left for outer Y label and Y levels labels
	right.cols = 0;		# additional column on the right for the legend colorbar
	top.rows = 0;			# additional rows on top of the figure for the outer X label and X levels labels
	bottom.rows = 0;	# additional bottom row for inner X label
	
	outY.levlabs = FALSE;
	outX.levlabs = FALSE;
	outX.lab = FALSE;
	outY.lab = FALSE;
	outY.blank.column = FALSE;    # provide space for the levels of the inner Y variable displayed in the Y axis of every cell in the left-most column
	inner.xlab = FALSE;
	inner.ylab = FALSE;	

	if(out.Y[["lab"]] != FALSE)					      { 	left.cols = left.cols + 1;	outY.lab = TRUE; 				    }
	if(out.Y[["levels.lab"]] != FALSE)	      {  	left.cols = left.cols + 1;	outY.levlabs = TRUE; 		    }
	if(inner.Y.par[["lab"]] != FALSE)         {   left.cols = left.cols + 1;  inner.ylab = TRUE; 			    }
	if(inner.Y.par[["levels.loc"]] == "left") {   left.cols = left.cols + 1;  outY.blank.column = TRUE;  }
	
	if(out.X[["lab"]] != FALSE)					{		top.rows = top.rows + 1;	outX.lab = TRUE; 					}
	if(out.X[["levels.lab"]] != FALSE)	{  	top.rows = top.rows + 1;	outX.levlabs = TRUE; 			}
		
	if(show.colorbar)										          {		right.cols = right.cols + 1;            						      }
	if(inner.X.par[["lab"]] != FALSE)		          {		bottom.rows = bottom.rows + 1;	inner.xlab = TRUE; 	      }		
	if(inner.X.par[["levels.loc"]] == "bottom")		{		bottom.rows = bottom.rows + 1;	inX.blank.row = TRUE; 	}		

	counter = 1;	
	if(outY.lab){ counter = 2; } 				# first row starts with counter=2, because 1
																							# is used repeatedly for the outer Y variable label (left-most columns)
	
	# Compose the layout matrix by rows
		
	row1 = NA; 
	row2 = NA;
	left.most.col = NULL;
	dummy.blank.col = NULL;
	
	layout.matrix = NA;
	bottom.row = NA;
	if(!is.null(heatmaps.per.row)){ 
		if(is.na(heatmaps.per.row)){ 
			heatmaps.per.row = NULL; 
		}
		else if(heatmaps.per.row == outerColumns){
			heatmaps.per.row = NULL;
		} 
	}
	
	width.Ylab              = character(0); # inches of the left most column containing the name of the outer Y variable
	width.Ylevels           = character(0); # inches of the second left-most column with the levels of the outer Y variable
	width.innerYlab		      = character(0); # inches of the third left-most column with the (repeated) name of the inner Y variable
	width.blank.margin.col  = character(0); # inches of the (dummy) fourth left-most column that gives space for inner Y levels when inner.Y.par[["levels.loc"]] == "left"
	width.colorbar 		      = character(0); # inches of the right-most column that contains the legend colorbar
	height.Xlab             = character(0); # inches of the top row containing the name of the outer X variable
	height.Xlevels          = character(0); # inches of the second top row containing the levels of the outer X variable
	height.innerXlab	      = character(0); # inches of the bottom row containing the name of the inner X variable
	height.blank.margin.row = character(0); # inches of the (dummy) second bottom row that gives space for inner X levels when inner.X.par[["levels.loc"]] == "bottom"

	actual.columns = min(outerColumns, heatmaps.per.row);     # this also works OK when heatmaps.per.row is NULL
	
	# First row of the layout: outer X variable label plus additional columns
	if(outX.lab){	
		row1 = c( rep(0,left.cols), rep(counter, actual.columns), rep(0, right.cols) ); 
		counter = counter + 1;
		height.Xlab = out.X[["lab.width"]];
	}
	
	# Second row: outer X level names (if the user has asked so)
	if(outX.levlabs){ 
		row2 = c( rep(0, left.cols), seq(from = counter, to = counter + actual.columns - 1), rep(0, right.cols) );
		counter = (counter + actual.columns - 1) + 1;		# starting value for correlative layout integers		
		height.Xlevels = out.X[["levels.lab.width"]];
	}
	
	# Left-most column: vertical label of the outer Y axis (all 1 in the layout)
	if(outY.lab){
		if(is.null(heatmaps.per.row)) { left.most.col = rep(1, outerRows); }
		else													{ left.most.col = rep(1, outerRows * ceiling(outerColumns / heatmaps.per.row)); }
		width.Ylab = out.Y[["lab.width"]];
	}
	# Second left-most column: levels of the outer Y axis
	if(outY.levlabs){
    width.Ylevels = out.Y[["levels.lab.width"]];
	}
	# Third left-most column: name of the inner Y variable
	if(inner.ylab){
		width.innerYlab = inner.Y.par[["lab.width"]];
	}
	# Fourth left-most column: blank space for the levels of the inner Y variable displayed in the inner Y axis of every heatmap in the left part
	if(outY.blank.column){
    width.blank.margin.col = lcm(longest.inner.Y.lev.width * 2.54); # Convert longest.inner.Y.lev.width from inches to centimeters
	}

	# Second bottom row: blank space for the levels of the inner X variable displayed below the inner X axis of every heatmap at the bottom 
	if(inX.blank.row){
    height.blank.margin.row = lcm(highest.inner.X.level * 2.54); # Convert highest.inner.X.lev.width from inches to centimeters
	}
	
	# Bottom row
	if(inner.xlab){
		height.innerXlab = inner.X.par[["lab.width"]];
	}
	
	if(show.colorbar){
		width.colorbar = lcm(2);
	}
			
	if(is.null(heatmaps.per.row)){	# General case: the layout matrix is easy to generate 	
		# This is the number of columns of a submatrix of the total matrix that does not contain the 2 upper rows
		# nor the left-most column and the additional bottom row
		
		ncols = actual.columns;	                        
		if(outY.levlabs)			{ ncols = ncols + 1;  }	
		if(show.colorbar)			{ ncols = ncols + 1;  } 			# not counting the left-most Y variable label column	
		if(inner.ylab)		    {	ncols = ncols + 1;	}
		layout.matrix = matrix(data = counter : (counter + ncols*outerRows - 1), nrow = outerRows, ncol = ncols, byrow = TRUE);
	}
	else{					# Special case: some cells must be set to 0 because no plot will be displayed on them
	
		rows.per.ylevel = ceiling(outerColumns / heatmaps.per.row);
		
		## Create a block for a complete ylevel (with rows.per.ylevel rows)
		inner.block = matrix(data = counter : (counter + heatmaps.per.row * rows.per.ylevel - 1), nrow = rows.per.ylevel, ncol = heatmaps.per.row, byrow = TRUE);
		extra.col.Y.levlabs = NA;
		if(outY.levlabs){ 																	    # Append left column if necessary
			extra.col.Y.levlabs = rep(counter, rows.per.ylevel);	# extra column on the left with the appropriate counter value
			inner.block = 1 + inner.block; 										    # add an offset of 1 
		}
				
    if(inner.ylab){ 																	      # Redefine the inner block since inner Y label must be repeated in every heatmap row
                                                            ## so it goes consecutive with the cell numbering 
      counter = min(inner.block);
      inner.block = matrix(data = counter : (counter + (heatmaps.per.row + 1)*rows.per.ylevel - 1), nrow = rows.per.ylevel, ncol = heatmaps.per.row + 1, byrow = TRUE);
		}
		
		## Set to 0 those cells of the last row of the block that should not display any heatmap
		last.filled.col = outerColumns %% heatmaps.per.row;
		to.zero = (1:heatmaps.per.row) > last.filled.col;		# General case, when there is no inner.ylab column
		if(inner.ylab){ to.zero = c(FALSE, to.zero); }      # Particular case when there is one more column on the left for inner Y label
		      
		inner.block[rows.per.ylevel, to.zero] = 0;

		## Append right column
		colorbar.counter = max(inner.block) + 1;		
		if(show.colorbar){																	# Append right column if necessary
			inner.block = cbind(inner.block, rep(colorbar.counter, rows.per.ylevel));
		}	
		## Append left column calculated at the beginning
		if(outY.levlabs){
			inner.block = cbind(extra.col.Y.levlabs, inner.block);			# append the column on the left
			counter = counter - 1;
		}
		
		## Duplicate this block and add the offset corresponding to each of the outer Y levels
		block.offset = max(inner.block) - counter + 1;
		all.blocks = lapply(X = 0:(outerRows-1), FUN = function(x, mat, offset){ 
				mask = (mat != 0); 
				return( (mat + (x*offset)) * mask ); 
			}, 
			mat = inner.block, offset = block.offset
		);
		layout.matrix = do.call(rbind, all.blocks);
	}

  actual.rows = nrow(layout.matrix);

  # Dummy blank column providing space for inner Y levels of the inner Y axis
  if(outY.blank.column){    
    layout.matrix = cbind(rep(NA, nrow(layout.matrix)), layout.matrix);
    current.col = 1;
    if(outY.levlabs) { layout.matrix[,current.col] = layout.matrix[,current.col + 1];  current.col = current.col + 1; }
    if(inner.ylab)   { layout.matrix[,current.col] = layout.matrix[,current.col + 1];  current.col = current.col + 1; }
    layout.matrix[,current.col] = 0;    # dummy column
  }
  
	if(outY.lab)					{ layout.matrix = cbind(left.most.col, layout.matrix);	}
	if(outX.levlabs)			{ layout.matrix = rbind(row2, layout.matrix);						}
	if(outX.lab)					{	layout.matrix = rbind(row1, layout.matrix);						}
	
	# Dummy blank row providing space for inner X levels of the inner X axis
	if(inX.blank.row){	
    dummy.row = rep(0, ncol(layout.matrix));
    layout.matrix = rbind(layout.matrix, dummy.row);
	}
	
	# Bottom row: inner X variable name
	if(inner.xlab){
		counter = max(layout.matrix) + 1;
		bottom.row = c( rep(0, left.cols), seq(from = counter, to = counter + actual.columns - 1), rep(0, right.cols) );
		layout.matrix = rbind(layout.matrix, bottom.row);
	}
	
	# Annotations?
	if(annotate){	# The cells on the top left corner are fussioned to allow for an annotation label
		layout.matrix[layout.matrix > 0] = layout.matrix[layout.matrix > 0] + 1;
		layout.matrix[1:top.rows, 1:left.cols] = 1;
	}
	
	colnames(layout.matrix) = NULL;
	rownames(layout.matrix) = NULL;

  # Set column widths and row heights. Variables that remain as character(0) are ignored in practice when doing the concatenation.
  layout.col.widths = c(width.Ylab,	width.Ylevels, width.innerYlab, width.blank.margin.col, rep(1, actual.columns), width.colorbar);
  layout.row.heights = c(height.Xlab, height.Xlevels, rep(1, actual.rows), height.blank.margin.row, height.innerXlab); 

  arg.list = list(mat = layout.matrix, widths = layout.col.widths, heights = layout.row.heights);	
	return(arg.list);
}

# _____________________________________________________________________________________________________________

.compute.default.margins <- function(longest.innerX, longest.innerY, heat.cell.par, heat.axes.par, heatmap.titles, inner.X.levels.loc, 
                                    inner.Y.levels.loc, inner.X.las){

		# ---------------------------------------------------------------------------------------------------------
		# Compute the length (in inches) of the longest label of the inner Y levels (preliminary computations)
		# ---------------------------------------------------------------------------------------------------------

		actual.cex = par("cex");
		if(!is.null(heat.cell.par$cex))						{ 			heat.cell.par$cex; }
		
		myaxis.font = par("font.axis");	
		myaxis.cex = par("cex.axis");	
		if(!is.null(heat.axes.par$font))					{			myaxis.font = heat.axes.par$font;				}
		else{	
			if(!is.null(heat.cell.par$font.axis))		{			myaxis.font = heat.cell.par$font.axis;	}		
		}
		if(!is.null(heat.cell.par$cex.axis))			{			myaxis.cex = heat.cell.par$cex.axis;		}
			
		# ---------------------------------------------------------------------------------------------------------
		# Compute the height (in inches) of the individual heatmap titles
		# ---------------------------------------------------------------------------------------------------------
		mytitles.font = par("font.main");
		mytitles.cex = par("cex.main");
		if(!is.null(heatmap.titles)){
			if(!is.null(heat.cell.par$font.main))		{			mytitles.font = heat.cell.par$font.main;		}		
			if(!is.null(heat.cell.par$cex.main))		{			mytitles.cex	= heat.cell.par$cex.main;			}
		}

		old.cex = par("cex");	
		par(cex = 0.66);									# Just temporarily so the next calls to strwidth() and strheight() behave exactly 
																			# the same way as when heatmaps are displayed later on
    inches.innerY.maxlen = strwidth(longest.innerY, cex = myaxis.cex, font = myaxis.font, units = "inches");
    
    # vertical labels in the inner X axis
    inches.innerX.maxheight = strwidth(longest.innerX, cex = myaxis.cex, font = myaxis.font, units = "inches");
    if( inner.X.las == 0 || inner.X.las == 1){
      # horizontal labels in the inner X axis
			inches.innerX.maxheight = strheight("M", cex = myaxis.cex, font = myaxis.font, units = "inches");
		}
		
		inches.left.mar = 0;
		inches.heatmap.titles.height = 0;
		inches.heatmap.inner.X.levels.height = 0;
		
		if(inner.Y.levels.loc == "all"){ inches.left.mar = inches.innerY.maxlen;				}
		if(!is.null(heatmap.titles)){ inches.heatmap.titles.height = strheight("M", cex = mytitles.cex, font = mytitles.font, units = "inches"); }		
		if(inner.X.levels.loc == "all"){ inches.heatmap.inner.X.levels.height = strheight("M", cex = myaxis.cex, font = myaxis.font, units = "inches"); }
		
		left.mar = inches.left.mar + 0.1;
		top.mar = inches.heatmap.titles.height + 0.1;
		bottom.mar = inches.heatmap.inner.X.levels.height + 0.1;		
		
		par(cex = old.cex);												# Restore original value
		
		return(c(bottom.mar, left.mar, top.mar, inches.innerY.maxlen, inches.innerX.maxheight)); # Note the fourth value is NOT the right margin !!
}