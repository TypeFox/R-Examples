#' Heatmap plot of the ranking achieved by a target variable levels after all statistical pairwise comparisons in multi-parameter problem instances.
#'
#' @description \code{plot.SRCS}: Function to display a grid of heatmaps representing the statistical ranking of one level of the target factor (usually, the algorithm)
#'	vs the rest of levels of the target factor, over several problem configurations characterized by (at most) 3 parameters in addition to the target factor.
#' @param x An SRCS object containing columns for the names of the problem parameters (including the algorithm) and the rank obtained
#' 	by that algorithm when compared with the rest over the same problem configuration. Typically this is the object returned by a call to \code{\link{SRCSranks}}.
#' @param yOuter,xOuter Names of the variables in \code{x} that will be placed vertically (in the left-most part) and horizontally (on the top), respectively. 
#'	Each level of \code{yOuter} (resp. \code{xOuter}) corresponds to a complete row (complete column) of heatmaps in the grid.
#' @param yInner,xInner Names of the variables in \code{x} that will be placed on the Y axis and on the X axis of every heatmap, respectively. 
#'	Each level of \code{yInner} (resp. \code{xInner}) corresponds to a row (column) inside a heatmap.
#' @param zInner Name of the variable in \code{x} that will be represented with a color code inside every heatmap. Usually corresponds to the ranking column of \code{x},
#'	which will most often contain integer values (negatives are allowed). When the SRCS object being plotted has been returned by a call to \code{SRCSranks}, 
#'	this column is called "rank".
#'	For \code{\link{animatedplot}}, it should be a vector of strings containing the names of the ranks
#'  columns that will be depicted, each at a time, sorted by time instant (from the earliest to the most recent).
#' @param out.Y.par,out.X.par A tagged list with parameters to configure how the labels of the outer Y and X variables and their levels are displayed. 
#'	Valid parameters and their default values are as follows:
#' \itemize{
#'	\item{\code{lab = TRUE}}{ Label with the name of the variable. Will be displayed vertically
#'	for the outer Y variable and horizontally for the outer X variable. Valid values: FALSE for no label; 
#'	NULL or TRUE for the name of the outer variable; and any string for a specific label. Defaults to TRUE.}
#'	\item{\code{lab.width = lcm(1)}}{ Width of the left-most column (for the outer Y variable) or top row (for the outer X variable) containing the name of the variable.}
#'	\item{\code{lab.textpar = list(cex = 1.6)}}{ Rest of parameters that will be passed to \code{\link{text}} to display this label.
#'		Parameter \code{cex} (text magnification factor) will be set 1.6 by default when not provided. }
#'	\item{\code{levels.lab = TRUE}}{ Whether a label should be displayed (or not) for every level of the variable. }
#'	\item{\code{levels.lab.width = lcm(1)}}{	Width of the row or column containing the levels of the variable.}
#'	\item{\code{levels.lab.textpar = list(cex = 1.4)}}{ Tagged list with more parameters that will be passed directly to \code{\link{text}} to display this label. 
#'		Parameter \code{cex} (text magnification factor) will be set 1.4 by default when not provided.
#'		NOTE: if present, the value of parameter \code{str} will always be overwritten by 0 (horizontal text) for the outer X, and 90 (vertical text) for the outer Y variable.}
#'	\item{\code{lab.bg = NULL}}{ Background color of the rectangle where the label is placed. Default is transparent. 
#'		No additional checks will be done concerning the validity of this parameter.}
#'	\item{\code{levels.bg = NULL}}{ Background color of the rectangle where the levels of the label are placed. Default is transparent. 
#'		No additional checks will be done concerning the validity of this parameter.}
#'	\item{\code{lab.border = NULL}}{ Border color of the rectangle where the label is placed. Defaults to NULL (no line).
#'		No additional checks will be done concerning the validity of this parameter.}
#'	\item{\code{levels.border	= NULL}}{ Line color of the rectangle border where the levels of this label are placed. Defaults to NULL (no line). 
#'		No additional checks will be done concerning the validity of this parameter.}
#' }
#' @param inner.X.par,inner.Y.par A tagged list with parameters to configure how the labels of the innter Y and X variables and their levels are displayed.
#'	Valid parameters and their default values are the following:
#'	\itemize{
#'		\item{\code{lab = TRUE}}{	Inner label to be shown. Valid values are FALSE for no label, NULL or TRUE for the name of variable passed as argument to \code{\link{plot.SRCS}},
#'															or any string for a specific label. Defaults to TRUE.}
#'		\item{\code{lab.width = lcm(1)}}{ Width of the optional space for the label of the inner Y variable. 
#'			The label will be repeated along the rows of the left-most column of heatmaps.}
#'		\item{\code{lab.textpar = list(cex = 1)}}{	Rest of parameters passed to \code{\link{text}} to display this label. }
#' 		\item{\code{levels.loc = c("bottom", "left", "all", "none")}}{ Location of the inner level labels: only in heatmaps of the left-most column or the bottom row, 
#'			or in every heatmap of the plot, or none. Defaults to "bottom" for the inner X variable and "left" for the inner Y variable.
#'			When levels.loc is set to "none", the value of params[["levels.at"]] is ignored.}
#' 		\item{\code{levels.at = NULL}}{ Levels of the inner variable where the label will be shown. Defaults to all the levels. They can be provided in any order, 
#'			since the order in which they will be displayed only depends on the order defined by the levels argument when that factor column of the data was created.}
#'		\item{\code{levels.las = 1}}{ Orientation of the level labels of this variable, defined as in \code{\link{axis}}. 1 for horizontal, 2  for vertical.}
#'	}
#' @param color.function A custom function that receives one argument (number of colors to be generated, (maxrank - minrank + 1) in our case) and returns a vector of that length
#'	with the hexadecimal codes of the colors to be used in the heatmaps, see \code{heat.colors} or \code{terrain.colors} for instance. Defaults to the \code{heat.colors} function.
#' @param colorbar.par Tagged list to configure the aspect of the colorbar legend displayed on the right part of the figure:
#'	\itemize{
#'		\item{\code{levels.at = NULL}}{ String vector: Levels at which the Y axis ticks of the colorbar will be shown. By default, three levels will be labeled: 
#'			0, the \code{min(x[[zInner]])} and \code{max(x[[zInner]])}.}
#'		\item{\code{hlines = TRUE}}{ Logical: whether black horizontal lines should be displayed in the colorbar to separate the colors. Defaults to TRUE. }
#'	}
#' @param heatmaps.per.row Maximum number of heatmaps displayed in a row of the grid. Useful when variable \code{xOuter} has too many levels so they can be
#'	splitted in two or more sub-rows of heatmaps, with all the sub-rows corresponding to a single level of the \code{yOuter} variable.
#' @param heatmaps.titles A vector of the same length as the total number of heatmaps, i.e. unique(x[[yOuter]]) * unique(x[[xOuter]]), containing the
#'	titles to be displayed on top of each heatmap. The elements of the vector will be associated to the heatmaps of the grid from left to right and from top to bottom.
#' @param show.colorbar Logical: whether a colorbar legend will be shown on the right of the figure (one for each row of heatmaps) or not. Defaults to TRUE
#' @param annotation.lab String with the annotation title that will be displayer on the top left corner. Defaults to NULL, indicating no annotation will be shown.
#' @param heat.cell.par Tagged list that will be passed to \code{\link{par}} just before displaying each heatmap. This way expert users can configure exactly
#'	the appearance of the heatmaps. No additional checks will be done concerning the validity of this list.
#' @param heat.axes.par Tagged list that will be passed to \code{\link{axis}} when creating the heatmap axes. No additional validity checks are done.
#'	The values of the arguments \code{side, at, labels} will always be replaced by suitable ones according to \code{inner.X.par[["levels.at"]]} or \code{inner.Y.par[["levels.at"]]}.
#' @param colorbar.cell.par Tagged list that will be passed to \code{\link{par}} just before showing each colorbar. No additional validity checks are done.
#' @param colorbar.axes.par Tagged list that will be passed to \code{\link{axis}} to draw the axes of the colorbar. No additional validity checks are done.
#' @param annotation.text.par Tagged list that will be passed to \code{\link{text}} to show an additional title on the top left corner. No additional validity checks are done.
#' @details \code{plot.SRCS} plots a grid with the results over all problem configurations, and should be applied to the object returned by \code{\link{SRCSranks}} with
#'	only one \code{performance} column. 
#'	
#'	\code{singleplot} is used for plotting only one heatmap for a subset of problem configurations 
#'	in which the outer X and Y parameters take a fixed value, and should be applied to the object returned by \code{\link{SRCScomparison}}. 
#'
#'	\code{animatedplot} creates a video from a sequence of plots, intended to show the temporal evolution of the ranking over time. 
#'	It should be applied only to the object returned by \code{\link{SRCSranks}} when the \code{performance} argument passed to it was a vector of strings, 
#'	each of them being the performance column of the data at a given time instant.
#' @note The function uses the base graphics system.
#' @examples
#'	# Example from a Machine Learning problem with noisy data
#'	ranks = SRCSranks(ML1, params = c("Dataset", "Noise type", "Noise ratio"), 
#' 	  target = "Algorithm", performance="Performance", maximize = TRUE, ncores = 2, 
#'	  paired = TRUE, pairing.col = "Fold");
#'	singleplot(ranks, yInner = "Noise type",
#'    xInner = "Noise ratio", Algorithm = "C4.5", Dataset = "glass")
#' 	plot(x = ranks, yOuter = "Dataset", xOuter = "Algorithm", yInner = "Noise type", 
#'	  xInner = "Noise ratio", out.X.par = list(levels.lab.textpar = 
#'	  list(col = "white"), levels.bg = "black", levels.border = "white"), 
#'	  out.Y.par = list(levels.bg = "gray"), colorbar.axes.par = list(cex.axis = 0.8), 
#'	  show.colorbar = TRUE)
#'	SRCScomparison(ranks, "Algorithm", Dataset = "automobile", `Noise type` = "ATT_GAUS", 
#'	  `Noise ratio`= 10, pvalues = FALSE)
#' @seealso \code{\link{text}, \link{par}, \link{axis}, \link{SRCSranks}, \link{animatedplot}, \link{singleplot}},
#' \link[RColorBrewer]{brewer.pal}, \link[colorspace]{RGB}
plot.SRCS <- function(x, yOuter, xOuter, yInner, xInner, 
											zInner = "rank",
											out.Y.par = list(), 	out.X.par = list(),
											inner.X.par = list(),	inner.Y.par = list(), 
											colorbar.par = list(),			
											color.function = heat.colors,		# function that will be called with (maxValue - minValue + 1) to generate heatmap colors
											heatmaps.per.row = NULL,		# in case there are many values for the outer X variable
											heatmaps.titles = NULL, 		# in case the user needs distinct titles for each heatmap.  
											show.colorbar = TRUE,				# TRUE to show a legend colorbar on the right of the plot
											annotation.lab = NULL,			# whether an annotation label will be displayer on the top left corner
											heat.cell.par = list(), 		# extra arguments for the par() function which only affect heatmaps cells. NO CHECKS
											heat.axes.par = list(),			# extra arguments for the axis() function which controls heatmaps axes. NO CHECKS
											colorbar.cell.par = list(),	# extra arguments for the par() function which only affect colorbars cells. NO CHECKS
											colorbar.axes.par = list(), # extra arguments for the axis() function which controls colorbars heatmaps. NO CHECKS
											annotation.text.par = list(), #extra arguments for text() in the annotation label. NO CHECKS
											...
	)
{
	data = x;
	out.Y.par = .check.outer.graphical.params(out.Y.par, yOuter);
	out.X.par = .check.outer.graphical.params(out.X.par, xOuter);
	inner.Y.par = .check.inner.graphical.params(inner.Y.par, "y");
	inner.X.par = .check.inner.graphical.params(inner.X.par, "x");
	colorbar.par = .check.colorbar.params(colorbar.par);

	colNames = names(data);
	for(param in c(yOuter, xOuter, yInner, xInner, zInner)){
		if(!is.null(param)){
			if(sum(param %in% colNames) < length(param)){
				stop(paste("ERROR:",param[which.min(param %in% colNames)],"not found in the frame data"));
			}
		}
  }
	
  if(is.null(xInner) || is.null(yInner) || is.null(zInner)){
		stop("inner X, Y and Z variables cannot be null");
  }
  if(is.null(yOuter) && is.null(xOuter)){
		stop("outer x and y variables cannot be both null at the same time. Use singleplot() instead");
  }
  
  if(length(show.colorbar) > 1 || typeof(show.colorbar) != "logical"){
		stop("show.colorbar argument must be a 1-element vector either TRUE or FALSE");
  }
  
  if(is.null(color.function)){
		color.function = heat.colors;
  }

	# --------------------------------------------------------------------------
	# Check validity of heatmaps.per.row
	# --------------------------------------------------------------------------	
  if(!is.null(heatmaps.per.row)){
		if(typeof(heatmaps.per.row) != "numeric" || length(heatmaps.per.row) > 1){
			stop("heatmaps.per.row must be a single positive integer value");
		}
		else if(heatmaps.per.row <= 0){
			stop("heatmaps.per.row must be a single positive integer value");
		}
  }

	# --------------------------------------------------------------------------  
  
  #TODO: check that the target and the rank variables are not the xInner and yInner
  #TODO: check that the rank variable matches zInner, and the target matches either yOuter or xOuter


	maxValue <- max(data[[zInner]], na.rm = TRUE);
	minValue <- min(data[[zInner]], na.rm = TRUE);
	mycolors = do.call(color.function, list(maxValue - minValue + 1)); # heat.colors palette by default

	## Count the number of values in each dimension
		
	outerRows = 1; outerColumns = 1; 
	uniqueOuterRows = NA;
	uniqueOuterColumns = NA;
	
	if(!is.null(xOuter)){ 
		uniqueOuterColumns = sort(unique(data[[xOuter]]));
		outerColumns = length(uniqueOuterColumns);	
	} else{ 
		uniqueOuterColumns = c("");	
		outerColumns = 1;
	}	
	
	if(!is.null(yOuter)){ 
		uniqueOuterRows = sort(unique(data[[yOuter]]), decreasing = TRUE);
		outerRows = length(uniqueOuterRows);	
	} else{	
		uniqueOuterRows = c("");
		outerRows = 1;
	}

	# --------------------------------------------------------------------------
	# Check validity of heatmaps.titles
	# --------------------------------------------------------------------------	
	if(!is.null(heatmaps.titles)){
		if(typeof(heatmaps.titles) != "character" || length(heatmaps.titles) != outerRows*outerColumns){
			stop(paste0("argument heatmaps.titles must be a character-string vector of the same length as the number of heatmaps to be displayed (",outerRows*outerColumns,")"));
		}
	}

	# --------------------------------------------------------------------------

	par(oma = c(0,0,0,0));
	par(mar = c(0,0,0,0));
	
	# ---------------------------------------------------------------------------------------------------------
	# Compute the width of inner Y levels and height of the individual heatmap titles (in inches) and store
	#	the results in variables left.mar and top.mar that will later be used for calling par(mai = ...) ) in heatmap cells
	# ---------------------------------------------------------------------------------------------------------	
	top.mar = 0;
	left.mar = 0;
	bottom.mar = 0;
	inches.longest.innerY = 0;
	inches.highest.innerX = 0;
	if(is.null(heat.cell.par[["mar"]]) && is.null(heat.cell.par[["mai"]])){		# Otherwise, use the heat.cell.par[["mar"]] (or "mai") vector provided by the user
			string.innerY = as.character(data[[yInner]]);
			string.innerX = as.character(data[[xInner]]);
			longest.innerY = "";									# Longest of the inner Y levels (either in the whole data or within the user-specified subset of levels to be labeled)
			longest.innerX = "";									# Longest of the inner X levels (either in the whole data or within the user-specified subset of levels to be labeled)
			
			if(!is.null(inner.Y.par[["levels.at"]]))		{		longest.innerY = inner.Y.par[["levels.at"]][which.max(nchar(inner.Y.par[["levels.at"]]))]; }
			else																				{		longest.innerY = string.innerY[which.max(nchar(string.innerY))];	}
			
			if(!is.null(inner.X.par[["levels.at"]]))		{		longest.innerX = inner.X.par[["levels.at"]][which.max(nchar(inner.X.par[["levels.at"]]))]; }
			else																				{		longest.innerX = string.innerX[which.max(nchar(string.innerX))];	}
			
			default.mar = .compute.default.margins(longest.innerX, longest.innerY, 
																						heat.cell.par, heat.axes.par, heatmaps.titles, inner.X.par[["levels.loc"]], inner.Y.par[["levels.loc"]],
																						inner.X.par[["levels.las"]]);
			
			bottom.mar = default.mar[1];
			left.mar = default.mar[2]; 
			top.mar = default.mar[3];			
      inches.longest.innerY = default.mar[4];
      inches.highest.innerX = default.mar[5];
	}
	
	# ---------------------------------------------------------------------------------------------------------	
	
	if(!is.null(heat.cell.par[["bg"]])){
		par(bg = heat.cell.par[["bg"]]);
	}
	
	annotate = !is.null(annotation.lab);	
	
	# Layout of the complete plot, including variable names and levels
	layout.args = .compose.layout(out.X.par, out.Y.par, inner.X.par, inner.Y.par, outerColumns, outerRows, heatmaps.per.row, show.colorbar, 
                               inches.longest.innerY, inches.highest.innerX, annotate);
	do.call(layout, layout.args);	
	#layout.show(n = max(layout.args[["mat"]]));
	#stop("END");
	
  # -----------------------------------------------------------------------
  # 													(OPTIONAL) ANNOTATION LABEL
  # -----------------------------------------------------------------------
	if(annotate){
		par(mar = c(0,0,0,0));
		plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');

		textargs = annotation.text.par;
		textargs[["labels"]] = annotation.lab;
		if(is.null(textargs[["x"]])){ textargs[["x"]] = 0.5; }
		if(is.null(textargs[["y"]])){ textargs[["y"]] = 0.5; }
		if(is.null(textargs[["cex"]])){ textargs[["cex"]] = 1.2; }
		do.call(text, textargs);
	}

  # -----------------------------------------------------------------------
  # 													OUTER Y LABEL
  # -----------------------------------------------------------------------
  if( ifelse(!is.null(out.Y.par[["lab"]]), out.Y.par[["lab"]] != FALSE, FALSE) ){
		label = "";
		if(out.Y.par[["lab"]] == TRUE){ 
			if(!is.null(yOuter))				{		label = yOuter; 					}
		}
		else 													{ label = out.Y.par[["lab"]]; }
		par(mar = c(0,0,0,0));
		plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
		
		# Plot background rectangle
		if(!is.null(out.Y.par[["lab.bg"]]) || !is.null(out.Y.par[["lab.border"]])){
			usr = par("usr");
			bordercolor = NA;
			bg = NA;
			if(!is.null(out.Y.par[["lab.border"]])){
				bordercolor = out.Y.par[["lab.border"]];
			}
			if(!is.null(out.Y.par[["lab.bg"]])){
				bg = out.Y.par[["lab.bg"]];
			}
			rect(usr[1], usr[3], usr[2], usr[4], border = bordercolor, col = bg);
		}
		
		textargs = out.Y.par[["lab.textpar"]];
		textargs[["labels"]] = label; 	
		textargs[["srt"]] = 90;		
		do.call(text, textargs);	# by default, font color is black
  }
  # -----------------------------------------------------------------------
  # 													OUTER X LABEL
  # -----------------------------------------------------------------------
  if( ifelse(!is.null(out.X.par[["lab"]]), out.X.par[["lab"]] != FALSE, FALSE) ){
		label = "";
		if(out.X.par[["lab"]] == TRUE){ 
			if(!is.null(xOuter))				{ label = xOuter; 						} 
		}
		else 													{ label = out.X.par[["lab"]]; }
		par(mar = c(0,0,0,0));
		plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');

		# Plot background rectangle
		if(!is.null(out.X.par[["lab.bg"]]) || !is.null(out.X.par[["lab.border"]])){
			usr = par("usr");
			bordercolor = NA;
			bg = NA;
			if(!is.null(out.X.par[["lab.border"]])){
				bordercolor = out.X.par[["lab.border"]];
			}
			if(!is.null(out.X.par[["lab.bg"]])){
				bg = out.X.par[["lab.bg"]];
			}
			rect(usr[1], usr[3], usr[2], usr[4], border = bordercolor, col = bg);
		}

		textargs = out.X.par[["lab.textpar"]];
		textargs[["labels"]] = label; 		
		textargs[["srt"]] = 0;		
		do.call(text, textargs);	# by default, font color is black
  }
  # -----------------------------------------------------------------------
  # 													OUTER X LEVELS
  # -----------------------------------------------------------------------
  if( ifelse(!is.null(out.X.par[["levels.lab"]]), out.X.par[["levels.lab"]] == TRUE, FALSE) ){		
		par(mar = c(0,0,0,0));
		names.list = uniqueOuterColumns[1:min(heatmaps.per.row, length(uniqueOuterColumns))];
		lapply(X = names.list, FUN = function(label, args){
			plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
			
			# Plot background rectangles
			if(!is.null(out.X.par[["levels.bg"]]) || !is.null(out.X.par[["levels.border"]])){
				usr = par("usr");
				bordercolor = NA;
				bg = NA;
				if(!is.null(out.X.par[["levels.border"]])){
					bordercolor = out.X.par[["levels.border"]];
				}
				if(!is.null(out.X.par[["levels.bg"]])){
					bg = out.X.par[["levels.bg"]];
				}
				rect(usr[1], usr[3], usr[2], usr[4], border = bordercolor, col = bg);
			}			
			
			textargs = args;
			textargs[["labels"]] = label; 		
			textargs[["srt"]] = 0;				
			do.call(text, textargs);
		}, args = out.X.par[["levels.lab.textpar"]]);
  }
	# -----------------------------------------------------------------------
  	
	for (i in 1:outerRows)
	{
		outY = uniqueOuterRows[i];						# label of the outer Y level
		boolVector = rep(TRUE, nrow(data));
		if(outY != ""){
			boolVector = boolVector & (data[[yOuter]] == outY);	# retain those with the i-th value in the outer row variable
		}
		
		original.boolVector = boolVector;
		
		# -----------------------------------------------------------------------
		# 													OUTER Y LEVEL "outY"
		# -----------------------------------------------------------------------
		if( ifelse(!is.null(out.Y.par[["levels.lab"]]), out.Y.par[["levels.lab"]] == TRUE, FALSE) ){
			par(mar = c(0,0,0,0));
			plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');

			# Plot background rectangle
			if(!is.null(out.Y.par[["levels.bg"]]) || !is.null(out.Y.par[["levels.border"]])){
				usr = par("usr");
				bordercolor = NA;
				bg = NA;
				if(!is.null(out.Y.par[["levels.border"]])){
					bordercolor = out.Y.par[["levels.border"]];
				}
				if(!is.null(out.Y.par[["levels.bg"]])){
					bg = out.Y.par[["levels.bg"]];
				}
				rect(usr[1], usr[3], usr[2], usr[4], border = bordercolor, col = bg);
			}
	
			textargs = out.Y.par[["levels.lab.textpar"]];
			textargs[["labels"]] = outY; 	
			textargs[["srt"]] = 90;
			do.call(text, textargs);
		}
		# -----------------------------------------------------------------------		
		
		for (j in 1:outerColumns)
		{	
			in.left.most.column = ifelse(!is.null(heatmaps.per.row), j %% heatmaps.per.row == 1, j == 1);
			if( in.left.most.column )  {	# if we are starting a row 
					# -----------------------------------------------------------------------
					# 													INNER Y LABEL (repeated across all rows)
					# -----------------------------------------------------------------------
					label = "";
					if(inner.Y.par[["lab"]] == TRUE)	{ 
						if(!is.null(yInner))						{	label = yInner; 							}
					}
					else 															{ label = inner.Y.par[["lab"]]; }
					par(mar = c(0,0,0,0));
					plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
					
					textargs = inner.Y.par[["lab.textpar"]];
					textargs[["labels"]] = label; 
					textargs[["srt"]] = 90;
					do.call(text, textargs);
					# -----------------------------------------------------------------------		
			}
		
			outX = uniqueOuterColumns[j];
			if(outX != ""){
				boolVector = original.boolVector & (data[[xOuter]] == outX); # retain those with the j-th value in the outer column variable
			}			
			dataSubset = data[boolVector, c(xInner, yInner, zInner)];
		
			ordered.indices = order(dataSubset[[xInner]], dataSubset[[yInner]]);
			dataSubset = dataSubset[ordered.indices, ];

			uniqueInnerY = unique(dataSubset[[yInner]]);
			uniqueInnerX = unique(dataSubset[[xInner]]);

			heatMatrix = matrix(dataSubset[[zInner]], nrow = length(uniqueInnerY), ncol = length(uniqueInnerX), byrow=FALSE);
			rownames(heatMatrix) = uniqueInnerY;
			colnames(heatMatrix) = uniqueInnerX;

			# Show inner X and Y levels only in those levels specified by the user: delete the rest of names
			if(!is.null(inner.X.par[["levels.at"]])){				
				notfound = !(uniqueInnerX %in% inner.X.par[["levels.at"]]);		# inner X labels not found among those provided by the user
				colnames(heatMatrix)[notfound] = "";
			}
			if(!is.null(inner.Y.par[["levels.at"]])){				
				notfound = !(uniqueInnerY %in% inner.Y.par[["levels.at"]]);		# inner Y labels not found among those provided by the user
				rownames(heatMatrix)[notfound] = "";
			}			
			
			# -----------------------------------------------------------------------
			# 								SPECIFIC SETTINGS FOR HEATMAP CELLS
			# -----------------------------------------------------------------------								
			oldpar = par(names(heat.cell.par));				# Retrieve the original par values
			names(oldpar) = names(heat.cell.par);	
					
			if(length(heat.cell.par) > 0){			par(heat.cell.par);			}
			if(is.null(heat.cell.par[["mar"]]) && is.null(heat.cell.par[["mai"]])){			
				# Set default margins for every heatmap depending on the length of the longest inner Y level and whether there are titles on the heatmaps or not 				
				par(mai = c(bottom.mar, left.mar, top.mar, 0.1));
			}
			# -----------------------------------------------------------------------						
			
			mainTitle = "";
			if( !is.null(heatmaps.titles) ){		mainTitle = heatmaps.titles[(i-1)*outerColumns + j]; }
			
			# Separation from axes label to the axis, from ticks label to the tick marker, and from tick marker to axis line. Default is (3,1,0)
			if(is.null(heat.cell.par[["mgp"]])){	par(mgp = c(4,1,0));	}

			image(x = seq.int(1:length(colnames(heatMatrix))), 
						y = seq.int(1:length(rownames(heatMatrix))), 
						z = t(heatMatrix), zlim = c(minValue, maxValue),
						col=mycolors, xlab = "", ylab = "", xaxt="n", yaxt="n", main = mainTitle, xaxs = "i", yaxs="i");	

			par(oldpar);		# Restore the original par values
			
			# -----------------------------------------------------------------------
			# 								HEATMAP AXES IF NECESSARY
			# -----------------------------------------------------------------------
			# 								Y axis
			# -----------------------------------------------------------------------			
			# If we are in the left part of the grid (or if all heatmaps should have Y axis), add a vertical axis with inner Y labels to every plot
			# Otherwise, add a vertical axis with no labels	
			mylabels = NULL;
			if( in.left.most.column || inner.Y.par[["levels.loc"]] == "all"){			mylabels = rownames(heatMatrix);		}
			else																														{ 		mylabels = rep("",length(rownames(heatMatrix)));	} 
			axisparams = heat.axes.par;
			axisparams[["side"]] = 2;
			axisparams[["at"]] = seq.int(1:length(rownames(heatMatrix)));
			axisparams[["labels"]] = mylabels;
			axisparams[["las"]] = inner.Y.par[["levels.las"]];
			do.call(axis, axisparams); 
			# -----------------------------------------------------------------------
			# 								X axis
			# -----------------------------------------------------------------------						
			mylabels = NULL;
			axisparams = NULL;
			# If we are in the last outer row of the grid (or if all heatmaps should have X axis), add the horizontal axis to every plot
			if(i==outerRows || inner.X.par[["levels.loc"]] == "all"){					mylabels = labels = colnames(heatMatrix);					}
			else 																										{					mylabels = rep("", length(colnames(heatMatrix)))	 }
			axisparams = heat.axes.par;
			axisparams[["side"]] = 1;
			axisparams[["at"]] = seq.int(1:length(colnames(heatMatrix)));
			axisparams[["labels"]] = mylabels;
			axisparams[["las"]] = inner.X.par[["levels.las"]];
			do.call(axis, axisparams);
		}
		
		if(show.colorbar){
			# -----------------------------------------------------------------------
			# 		SPECIFIC SETTINGS FOR COLORBAR LEGEND (repeated across the rows)
			# -----------------------------------------------------------------------		
			oldpar = par(names(colorbar.cell.par));		# Retrieve the original par values
			names(oldpar) = names(colorbar.cell.par);
			#oldpar = par(no.readonly=T);
			if(length(colorbar.cell.par) > 0){			par(colorbar.cell.par);			}
			if(is.null(colorbar.cell.par[["mar"]]) && is.null(colorbar.cell.par[["mai"]])){
				# Compute default margins for every heatmap depending on the length of the longest inner Y level and whether there are titles on the heatmaps or not 				
				par(mai = c(bottom.mar,0.1,top.mar,0.5));
			}
			
			axisargs = colorbar.axes.par;
			colorbar.data = matrix(data=seq(from=minValue,to=maxValue), nrow = 1,ncol = (maxValue - minValue + 1), byrow = TRUE);
			colnames(colorbar.data) = seq(from=minValue, to=maxValue);

			if(!is.null(colorbar.par[["levels.at"]])){ 	# delete those names not found among user provided labels
				names = colnames(colorbar.data);
				notfound = !(names %in% colorbar.par[["levels.at"]]);		# inner X labels not found among those provided by the user
				names[notfound] = "";
				axisargs[["at"]] = names[names != ""]; #seq.int(minValue, maxValue);
				axisargs[["labels"]] = axisargs[["at"]]; #colnames(colorbar.data);
			}
			
			# Now plot the legend on the right (ramp colorbar)		
			image(y = seq(from=minValue,to=maxValue),
						z = colorbar.data,					
						col = mycolors, xaxt="n", yaxt="n", ylab="");

			if(colorbar.par[["hlines"]]){	# Horizontal lines to separate the colors
				sapply(X = (minValue-0.5):(maxValue+0.5), FUN = function(value){ abline(h=value); });
			}
			
			axisargs[["side"]] = 4;
			axisargs[["las"]] = 1;				
			do.call(axis, axisargs);
			par(oldpar);		# Restore the original par values
			# -----------------------------------------------------------------------		
	  }

	}
	
  if( ifelse(!is.null(inner.X.par[["lab"]]), inner.X.par[["lab"]] != FALSE, FALSE) ){
    for (j in 1:min(heatmaps.per.row, outerColumns)){	  # this also works if heatmaps.per.row is NULL
      # -----------------------------------------------------------------------
      # 													INNER X LABEL
      # -----------------------------------------------------------------------
      label = "";
      if(inner.X.par[["lab"]] == TRUE)	{ 
        if(!is.null(xInner))						{	label = xInner; 							}
      }
      else 															{ label = inner.X.par[["lab"]]; }
      par(mar = c(0,0,0,0));
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
      textargs = inner.X.par[["lab.textpar"]];
      textargs[["labels"]] = label; 
      textargs[["srt"]] = 0;
      do.call(text, textargs);
      #text(x = 0.56, y = 0.45, xInner, cex = 1, col = "black");	    
      # -----------------------------------------------------------------------			
    }
	}	
}

# ___________________________________________________________________________
#' Video sequence composed of plots of successive data, generated with a call to external software ImageMagick
#'
#' @rdname plot.SRCS
#' @description \code{animatedplot}: Function to generate an animated video consisting of a temporal sequence of grid plots like those generated by \code{\link{plot.SRCS}}.
#' 	The function requires software ImageMagick has been installed.
#' @param filename Name of the output video file, including the extension. It is strongly recommended that the name ends in ".gif" to preserve most of image quality.
#' @param path.to.converter String with the full path to the converter program delivered with ImageMagick, e.g. "C:/Program Files/ImageMagick-<version>/convert.exe"
#' @param width,height Width and height, in pixels, of the result video. Both default to 800
#' @param res Nominal resolution (in ppi) of the output video. Used to set text size and line widths. Defaults to 100. See \code{\link{png}, \link{jpeg}, \link{bmp}, \link{tiff}}.
#' @param pointsize Point size argument to be passed to the functions that print to image. Defaults to 16.
#' @param delay Time delay (in 1/100th of a second) spent in each of the images that compose the video. Defaults to 30, i.e. 0.3 seconds.
#' @param type The type of image file generated for each frame. The image files will be then joined together into a video. Should be one of \code{"png", "jpeg", "bmp", "tiff"}.
#' @param quality The quality of the images, in a scale from 1 to 100. The less the quality, the more the compression and the smaller the file size.
#' @param compression (For TIFF format only) Used to indicate the kind of compression. Must be one of \code{"none", "rle", "lzw", "jpeg"}. Ignored if \code{type} is not \code{"tiff"}.
#' @param annotations Vector of strings with the annotation label of every image of the video. Should have the same length as \code{zInner}. Defaults to NULL (no annotations).
#' @param ... (In \code{animatedplot}): Rest of optional parameters that will be passed to \code{\link{plot.SRCS}} to plot every frame.
#' (In \code{singleplot}): A number of named arguments of the form \code{variable = value}, where \code{variable} is a column in \code{x}, for subsetting
#' \code{x} in a way that there exists exactly one occurrence of all the levels of \code{zInner} for each combination of yInner,xInner.
#' @examples
#' # ---------------------------------------------------
#' \dontrun{
#' mat = matrix(NA, nrow = nrow(MPBall), ncol = ncol(MPBall))
#' # First, take the average of the previous performance columns up to each change point
#' for(j in 6:ncol(MPBall)){
#'   mat[,j] = rowSums(MPBall[,5:j])/(j-5+1)
#' }
#' MPBall[,6:ncol(MPBall)] = mat[,6:ncol(MPBall)]
#'
#' ranksall = SRCSranks(MPBall, params = c("Dim", "CF", "Severity"), target="Algorithm", 
#'    test = "tukeyHSD", performance=paste("OffError", seq(from=1, to = 100, by = 24), 
#'    sep = "_"), maximize = FALSE, ncores = 2)
#' 
#' # Adjust argument path.to.converter to point to ImageMagick convert utility
#' animatedplot(x = ranksall, filename = "MPBconv_reduced.gif", 
#'	             path.to.converter = "C:/Program Files/ImageMagick-6.8.8-Q8/convert.exe",
#'	             yOuter = "Algorithm", xOuter = "Dim", yInner = "CF", xInner = "Severity",
#'	             zInner = paste0("rank",1:5), delay = 30,
#'	             annotations = paste0("At change ",seq.int(from = 1, to = 100, by = 24)),
#'	             inner.Y.par = list(levels.at = c("40", "200", "400", "600", "800", "1000"), 
#'               lab = "Change\nfrequency", levels.loc = "left"), 
#'	             heat.cell.par = list(pty = "s"),
#'	             inner.X.par = list(levels.at = c("2", "8", "14")), 
#'	             out.Y.par = list(levels.lab.textpar = list(cex = 1, col = "white"), 
#'               levels.bg = "black", levels.border = "white"), 
#'	             out.X.par = list(lab = "Dimension", levels.bg = "gray"),
#'	             colorbar.par = list(levels.at = c("-2", "0", "2")),
#'	             colorbar.axes.par = list(cex.axis = 0.8), 
#'	             show.colorbar = TRUE, height = 500
#'             )
#' # The full dataset (20 MB) can be downloaded from 
#' # http://decsai.ugr.es/~pjvi/SRCSfiles/MPBall.RData
#' # (the average must still be computed before plotting, just as in the example above)
#' # Check the script in http://decsai.ugr.es/~pjvi/SRCSfiles/DOPvideoScript.R
#' }
animatedplot <- function(x, filename, path.to.converter,
											yOuter, xOuter, yInner, xInner, zInner,
											width = 800, height = 800, res = 100, pointsize = 16,
											delay = 30, type = c("png", "jpeg", "bmp", "tiff"), quality = 75, compression = c("none", "rle", "lzw", "jpeg", "zip"),
											annotations = NULL,
											...
	)
{
	data = x;
	
	colNames = names(data);
	for(param in c(yOuter, xOuter, yInner, xInner, zInner)){
		if(!is.null(param)){
			if(!(param %in% colNames)){
				stop(paste("ERROR:",param,"not found in the frame data"));
			}
		}
  }
	
	type = match.arg(type);
	if(type == "tiff"){ 
		compression = match.arg(compression); 
	}
	if(type == "jpeg" && (quality <= 0 || quality > 100)){
		stop("quality must be a positive number greater than 0 and not greater than 100");
	}
	nplots = length(zInner);
	
	if(!is.null(annotations)){
		if(length(annotations) != nplots){
			msg = paste0("the annotation vector must have the same length as zInner, the number of columns to be plotted (", nplots, ")");
			stop(msg);
		}
	}	
	
	currdev = dev.cur();
	
	if(type == "png"){					png(filename = paste0(filename, "%02d.png"), width = width, height = height, pointsize = pointsize); }
	else if(type == "jpeg"){   jpeg(filename = paste0(filename, "%02d.jpeg"), width = width, height = height, pointsize = pointsize, quality = quality); }
	else if(type == "bmp"){		  bmp(filename = paste0(filename, "%02d.bmp"), width = width, height = height, pointsize = pointsize); }
	else{											 tiff(filename = paste0(filename, "%02d.tiff"), width = width, height = height, pointsize = pointsize, compression = compression); }
	
	mypar = par(no.readonly = TRUE);
	
	plotargs = list(...);
	plotargs[["x"]] = data;
	plotargs[["yOuter"]] = yOuter;
	plotargs[["xOuter"]] = xOuter;
	plotargs[["yInner"]] = yInner;
	plotargs[["xInner"]] = xInner;	
	
	for(i in 1:nplots){		
		plotargs[["zInner"]] = zInner[i];
		plotargs[["annotation.lab"]] = annotations[i];
		do.call(plot.SRCS, plotargs);
		par(mypar);
	}	
	
	while(dev.cur() != currdev){	
    dev.off();    # close device		 
  }
	
	# Delete the file to make sure it does not exist previously
	unlink(filename);
	# Compose the shell command 	

	# Now call ImageMagick
	#shell(shQuote(cmdstring));
	if(.Platform$OS.type == "windows"){
		wincmdstring = paste0(shQuote(path.to.converter), " -delay ",delay," -compress None -quality 100% ",filename,"*.png ", filename);	
		shell(wincmdstring);	
	}
	else{
		unixcmdstring = paste0(path.to.converter, " -delay ",delay," -compress None -quality 100% ",filename,"*.png ", filename);
		system(unixcmdstring);
	}
}

# ___________________________________________________________________________
#' Heatmap plot of the ranking achieved by target variable levels after statistical pairwise comparisons for a given combination of the remaining parameters
#' in multi-parameter problem instances.
#'
#' @rdname plot.SRCS
#' @description \code{singleplot}: Function to display either a single heatmap representing the statistical ranking of one level of the target factor (usually, the algorithm)
#'	vs the rest of levels of the target factor, over one single problem configurations defined by a combination of values for the problem configuration parameters.
#' @param labels.par Tagged list to configure how the labels will be displayed:
#' \itemize{
#'	\item{\code{xlab	= TRUE}}{ Label with the name of the variable for the X axis. Will be displayed horizontally. Valid values:
#'		FALSE for no label; NULL or TRUE for the name of the outer variable; and any string for a specific label. Defaults to TRUE.}
#'	\item{\code{ylab	= TRUE}}{ Analogous for  the Y axis.}
#' 	\item{\code{xlevels.at = NULL}}{ Levels of the X axis variable where the label will be shown. Defaults to all the levels.
#'		The levels can be provided in any order, since the order in which they will be depicted only depends on the original order defined 
#'	  when the corresponding factor column of the data was created.}
#'	\item{\code{ylevels.at = NULL}}{ Analogous for Y axis variable.}
#' }
#' @param haxis,vaxis Whether the X and the Y axes should be displayed or not. Defaults to TRUE for both.
#' @param title Title of the plot.
singleplot <- function(x, yInner, xInner, zInner = "rank",
														color.function = heat.colors,
                            labels.par = list(),                             
                            colorbar.par = list(), 
                            heat.axes.par = list(),		# no check
                            colorbar.axes.par = list(),	# no check
                            haxis = TRUE,
                            vaxis = TRUE,
                            title = "",
                            show.colorbar = TRUE,
                            ...
	)
{
	par(c(5,4,4,2)+0.1);
	data = x;
	labels.par = .check.labels.params(labels.par);
	colorbar.par = .check.colorbar.params(colorbar.par);
	
	params = list(...);
	nms = names(params);
	if(is.null(yInner) || is.null(xInner)){
		stop("both yInner and xInner must be non-null");
	}

	if(!(zInner %in% names(data))){
		stop("zInner column not found on this object");
	}
	if(sum(!(nms %in% names(data))) > 0){
		bad = params[[which.min(nms %in% names(data))]];
		stop(paste0("column ", bad, " was not found on this object"));
	} 
	vbool = rep(TRUE, nrow(data));
	
	# Create boolean vector of rows that fulfill all conditions simultaneously
	for(i in 1:length(params)){
		temp = (data[[ nms[i] ]] == params[[i]]);
		vbool = vbool & temp;
	};
	
	if(sum(vbool) == 0){
		stop("the data subset meeting all the specified column values is empty");
	}
	
	maxValue <- max(data[[zInner]]);
	minValue <- min(data[[zInner]]);
	mycolors = do.call(color.function, args = list(maxValue - minValue + 1)); # heat.colors palette by default
	
	dataSubset = data[vbool, c(xInner, yInner, zInner)];
	ordered.indices = order(dataSubset[[xInner]], dataSubset[[yInner]]);
	dataSubset = dataSubset[ordered.indices, ];

	uniqueInnerY = unique(dataSubset[[yInner]]);
	uniqueInnerX = unique(dataSubset[[xInner]]);
	
	if(nrow(dataSubset) > length(uniqueInnerY)*length(uniqueInnerX)){
		stop("too few values provided: the <column = value> pairs specified do not restrict the zInner sufficiently.\nTry passing more <column = value> pairs as arguments");
	}

	string.innerY = as.character(uniqueInnerY);
	longest.innerY = NULL;
	if(!is.null(labels.par[["ylevels.at"]]))		{		longest.innerY = labels.par[["ylevels.at"]][which.max(nchar(labels.par[["ylevels.at"]]))]; }
	else																				{		longest.innerY = string.innerY[which.max(nchar(string.innerY))];	}
	
	heatMatrix = matrix(dataSubset[[zInner]], nrow = length(uniqueInnerY), ncol = length(uniqueInnerX), byrow=FALSE);
	rownames(heatMatrix) = uniqueInnerY;
	colnames(heatMatrix) = uniqueInnerX;	
	
	# Show inner X and Y levels only in those levels specified by the user: delete the rest of names
	if(!is.null(labels.par[["xlevels.at"]])){				
		notfound = !(uniqueInnerX %in% labels.par[["xlevels.at"]]);		# inner X labels not found among those provided by the user
		colnames(heatMatrix)[notfound] = "";
	}
	if(!is.null(labels.par[["ylevels.at"]])){				
		notfound = !(uniqueInnerY %in% labels.par[["ylevels.at"]]);		# inner Y labels not found among those provided by the user
		rownames(heatMatrix)[notfound] = "";
	}			

	# minimum width needed for the left margin
	mywidth = strwidth(longest.innerY, units = "inches");

	if(show.colorbar){	layout(matrix(c(1,2), nrow = 1, ncol = 2), widths = c(1, 0.3)); }

	if(labels.par[["xlab"]] == TRUE){ 	myxlab = xInner;		}
	else{		
		if(labels.par[["xlab"]] == FALSE){ myxlab = ""; }
		else{ 														myxlab = labels.par[["xlab"]];	}
	}
	
	if(labels.par[["ylab"]] == TRUE){ 	myylab = yInner;		}
	else{		
		if(labels.par[["ylab"]] == FALSE){ myylab = ""; }
		else{ 														 myylab = labels.par[["ylab"]];	}
	}

	# set left margin 
	mymai = par("mai");	
	if(2*mywidth > mymai[2]){
		mymai[2] = 2*mywidth;
		par(mai = mymai);
		mymgp = par("mgp");
		par(mgp = c(as.integer(mywidth/par("csi"))+2, mymgp[2:3]));
	}	
		
	image(x = seq.int(1:length(colnames(heatMatrix))), 
				y = seq.int(1:length(rownames(heatMatrix))), 
				t(heatMatrix),col=mycolors, xlab = myxlab, ylab = myylab, xaxt="n", yaxt="n", main = title, xaxs = "i", yaxs="i");
	
	# -----------------------------------------------------------------------
	# 								HEATMAP AXES IF NECESSARY
	# -----------------------------------------------------------------------
	# 								Y axis
	# -----------------------------------------------------------------------			
	if(vaxis){			
		mylabels = rownames(heatMatrix);									
		heat.axes.par[["side"]] = 2; 
		heat.axes.par[["at"]] = seq.int(1:length(rownames(heatMatrix)));
		heat.axes.par[["labels"]] = mylabels;
		heat.axes.par[["las"]] = 1;
		do.call(axis, heat.axes.par);		
	}
	# -----------------------------------------------------------------------
	# 								X axis
	# -----------------------------------------------------------------------						
	if(haxis){
		mylabels = labels = colnames(heatMatrix);
		heat.axes.par[["side"]] = 1;
		heat.axes.par[["at"]] = seq.int(1:length(colnames(heatMatrix)));
		heat.axes.par[["labels"]] = mylabels;
		heat.axes.par[["las"]] = 1;
		do.call(axis, heat.axes.par);
	}
	
	if(show.colorbar){
			# -----------------------------------------------------------------------
			# 							SPECIFIC SETTINGS FOR COLORBAR LEGEND 
			# -----------------------------------------------------------------------		
			mymar = par("mar");			
			par(mar = c(mymar[1], 2.1, mymar[3], mymar[4])); 	
			axisargs = colorbar.axes.par;
			colorbar.data = matrix(data=seq(from=minValue,to=maxValue), nrow = 1,ncol = (maxValue - minValue + 1), byrow = TRUE);
			colnames(colorbar.data) = seq(from=minValue, to=maxValue);

			if(!is.null(colorbar.par[["levels.at"]])){ 	# delete those names not found among user provided labels
				names = colnames(colorbar.data);
				notfound = !(names %in% colorbar.par[["levels.at"]]);		# inner X labels not found among those provided by the user
				names[notfound] = "";
				axisargs[["at"]] = names[names != ""]; #seq.int(minValue, maxValue);
				axisargs[["labels"]] = axisargs[["at"]]; #colnames(colorbar.data);
			}
			
			# Now plot the legend on the right (ramp colorbar)		
			image(y = seq(from=minValue,to=maxValue),
						z = colorbar.data,					
						col = mycolors, xaxt="n", yaxt="n", ylab="");

			if(colorbar.par[["hlines"]]){	# Horizontal lines to separate the colors
				sapply(X = (minValue-0.5):(maxValue+0.5), FUN = function(value){ abline(h=value); });
			}
			
			axisargs[["side"]] = 4;
			axisargs[["las"]] = 1;				
			do.call(axis, axisargs);			
			# -----------------------------------------------------------------------		
			# Restore old margins
			par(mar = mymar);
	  }	
}