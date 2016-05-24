##  =================================================================
## 	Possible graphical options for the outer variable names and text:
##  =================================================================
##
##  lab	= TRUE							# - Label with the name of the variable. Will be displayed vertically 
##													# 	for the outer Y variable and horizontally for the outer X variable. Valid values:
##													# 	FALSE for no label; NULL or TRUE for the name of the outer variable; and any string for a specific label. Defaults to TRUE
##	lab.width = lcm(1)			# - Width of the left-most column or bottom row containing the name of the yOuter variable 
##	lab.textpar = list(cex = 1.6, x = 0.5, y = 0.5)	 # - Rest of parameters that will be passed to the graphics::text() function call to display this label.
##													#   Parameter cex defaults to 1.6 for outer X and Y labels if not specified.
##													#   NOTE: the value of parameter str will be ignored if specified by the user, and overwritten by 
##													#   str = 0 (for outer X) or str = 90 (for outer Y)
##	levels.lab = TRUE				# - Whether a label should be displayed (or not) for every level of the variable
##													# 	FALSE means no label while NULL or TRUE means level labels will be shown. Defaults to TRUE.
##	levels.lab.width = lcm(1)	# - Width of the row or column containing the levels of the variable 
##	levels.lab.textpar = list(cex = 1.4, x = 0.5, y = 0.5) # - Rest of parameters that will be passed to the graphics::text() function call to display this label.
##													#   Parameter cex defaults to 1.4 for outer X and Y labels if not specified.
##													#   NOTE: the value of parameter str will be ignored if specified by the user, and overwritten by 
##													#   str = 0 (for outer X) or str = 90 (for outer Y)
##	lab.bg = NULL						#   Background color of the rectangle where the label is placed. Default is transparent. 	NO CHECK
##	levels.bg = NULL				#		Background color of the rectangle where the levels of the label are placed. Default is transparent. NO CHECK
##	lab.border = NULL				#		Line color of the rectangle border where the label is placed. Defaults to NULL (no line). NO CHECK
##	levels.border	= NULL		#		Line color of the rectangle border where the levels of this label are placed. Defaults to NULL (no line). NO CHECK
.check.outer.graphical.params <- function(params, var.name){
	
	# Check for NAs	
	if(length(params) > 0){
		if(sum(is.na(params)) > 0){
			nms = names(params);
			stop(paste("NAs not allowed for any parameter. Invalid value (NA) for parameter", nms[which.max(is.na(params))]));
		}
	}

	result.list = list();
	
	# ---------------------------------------------------------------- params[["lab"]]
	
	if(!is.null(params[["lab"]])){
		if( typeof(params[["lab"]]) != "character" && typeof(params[["lab"]]) != "logical" ){
			stop("invalid type for graphic parameter lab: must be either a character string or TRUE/FALSE");
		}
	}
	else params[["lab"]] = TRUE;				# defaults to TRUE (display the name of the variable passed in the xOuter or yOuter arguments)	
	result.list[["lab"]] = params[["lab"]];
	params[["lab"]] = NULL;

  if(is.null(var.name)){    result.list[["lab"]] = FALSE;    }

	# ---------------------------------------------------------------- params[["lab.width"]]

	if(!is.null(params[["lab.width"]])){
		if(!is(params[["lab.width"]], "numeric") || length(params[["lab.width"]]) > 1){
			stop("invalid type for graphic parameter lab.width: must be numeric and positive");
		}		
		if(params[["lab.width"]] <= 0){
			stop("invalid value for graphic parameter lab.width: must be numeric and positive");
		}	
		else{		result.list[["lab.width"]] = lcm(params[["lab.width"]]);		}
	}
	else{ result.list[["lab.width"]] = lcm(1); }	# lab.width defaults to lcm(1), i.e. 1 cm
	params[["lab.width"]] = NULL;	

	# ---------------------------------------------------------------- params[["lab.textpar"]]

	if(!is.null(params[["lab.textpar"]])){
		if(typeof(params[["lab.textpar"]]) != "list"){
			stop("invalid type for graphic parameter lab.textpar: must be a tagged list");
		}
	}
	else{ params[["lab.textpar"]] = list(); }

	if(is.null(params[["lab.textpar"]][["cex"]])){		params[["lab.textpar"]][["cex"]] = 1.6;			}
	if(is.null(params[["lab.textpar"]][["x"]]))	 {    params[["lab.textpar"]][["x"]] = 0.5; 			}
	if(is.null(params[["lab.textpar"]][["y"]]))	 { 		params[["lab.textpar"]][["y"]] = 0.5;				}
	
	result.list[["lab.textpar"]] = params[["lab.textpar"]];
	params[["lab.textpar"]] = NULL;
	
	# ---------------------------------------------------------------- params[["levels.lab"]]
	
	if(!is.null(params[["levels.lab"]])){
		if(typeof(params[["levels.lab"]]) != "logical"){
			stop("invalid type for graphic parameter levels.lab: must be logical (either TRUE or FALSE)");
		}
		else{		
			result.list[["levels.lab"]] = params[["levels.lab"]];
		}
	}
	else{ result.list[["levels.lab"]] = TRUE; }	# lev.labs defaults to TRUE
	params[["levels.lab"]] = NULL;
	
	if(is.null(var.name)){    result.list[["levels.lab"]] = FALSE;    }
	
	# ---------------------------------------------------------------- params[["levels.lab.width"]]
	
	if(!is.null(params[["levels.lab.width"]])){
		if(!is(params[["levels.lab.width"]], "numeric") || length(params[["levels.lab.width"]]) > 1){
			stop("invalid type for graphic parameter levels.lab.width: must be numeric and positive");
		}	
		if(params[["levels.lab.width"]] <= 0){
			stop("invalid value for graphic parameter levels.lab.width: must be numeric and positive");
		}
		else{		result.list[["levels.lab.width"]] = lcm(params[["levels.lab.width"]]);		}
	}
	else{ result.list[["levels.lab.width"]] = lcm(1); }	# levels.lab.width defaults to lcm(1), i.e. 1 cm
	params[["levels.lab.width"]] = NULL;	
	
	# ---------------------------------------------------------------- params[["levels.lab.textpar"]]

	if(!is.null(params[["levels.lab.textpar"]])){
		if(typeof(params[["levels.lab.textpar"]]) != "list"){
			stop("invalid type for graphic parameter levels.lab.textpar: must be a tagged list");
		}
	}
	else{ params[["levels.lab.textpar"]] = list(); }
	
	if(is.null(params[["levels.lab.textpar"]][["cex"]])) {		params[["levels.lab.textpar"]][["cex"]] = 1.4;		}
	if(is.null(params[["levels.lab.textpar"]][["x"]]))	 {    params[["levels.lab.textpar"]][["x"]] = 0.5; 			}
	if(is.null(params[["levels.lab.textpar"]][["y"]]))	 { 		params[["levels.lab.textpar"]][["y"]] = 0.5;			}
	
	result.list[["levels.lab.textpar"]] = params[["levels.lab.textpar"]];
	params[["levels.lab.textpar"]] = NULL;
	
	result.list[["lab.bg"]] 				= params[["lab.bg"]];
	result.list[["levels.bg"]] 			= params[["levels.bg"]];				
	result.list[["lab.border"]] 		= params[["lab.border"]];
	result.list[["levels.border"]] 	= params[["levels.border"]];
	params[["lab.bg"]] = NULL;
	params[["levels.bg"]] = NULL;
	params[["lab.border"]] = NULL;
	params[["levels.border"]] = NULL;

#		Line color of the rectangle border where the label is placed. Defaults to NULL (no line). NO CHECK
##	levels.border	= NULL		#		Line color of the rectangle border where the levels of this label are placed. Defaults to NULL (no line). NO CHECK
	
	# ----------------------------------------------------------------	
	## Check for unrecognized tags in the list
	
	if(sum(!sapply(params, is.null)) > 0){
		nms = names(params);
		stop(paste("unrecognized graphical parameter:",nms[which.max(!sapply(params, is.null))]));
	}	
		
	return(result.list);
}

# ___________________________________________________________________________

##  ===========================================================================
## 	Possible graphical options for the heatmaps, inner variable names and text:
##  ===========================================================================
##
##	lab = TRUE,										# Inner label to be shown. Valid values are FALSE for no label, NULL or TRUE for the name of variable passed as argument to SRCS.plot,
##																# or any string for a specific label. Defaults to TRUE.
##	lab.width = lcm(1),						# width of the optional space for the label of the inner Y variable. 
##																# The label will be repeated along the rows of the left-most column of heatmaps
##	lab.textpar = list(cex = 1)	# Rest of parameters passed to function graphics::text() to display this label
##	levels.loc = c("bottom", "left", "all", "none"), 	# Location of the inner level labels: only in heatmaps of the left-most column or the bottom row, 
##																# or in every heatmap of the plot, or none. Defaults to "bottom" for the inner X variable and "left" for the inner Y variable
##																# When levels.loc is set to "none", the value of params[["levels.at"]] is ignored
##	levels.at = NULL,							# levels of the inner variable where the label will be shown. Defaults to all the levels.
##																# They can be provided in any order, since the order in which they will be depicted only depends on the
##																# original order defined by the levels argument when the corresponding factor column of the data was created
##	levels.las = 1,								# Orientation of the level labels of the variable. 1 for horizontal, 2 for vertical.
## inner.variable is either "y" or "x"
.check.inner.graphical.params <- function(params, inner.variable){

	result.list = list();

	# ---------------------------------------------------------------- params[["lab"]]
	
	if(!is.null(params[["lab"]])){
		if( typeof(params[["lab"]]) != "character" && typeof(params[["lab"]]) != "logical" ){
			stop("invalid type for graphic parameter lab: must be either a character string or TRUE/FALSE");
		}
	}
	else params[["lab"]] = TRUE;				# defaults to TRUE (display the name of the variable passed in the yInner argument)
	result.list[["lab"]] = params[["lab"]];
	params[["lab"]] = NULL;

	# ---------------------------------------------------------------- params[["lab.width"]]

	if(!is.null(params[["lab.width"]])){
		if(!is(params[["lab.width"]], "numeric") || length(params[["lab.width"]]) > 1){
			stop("invalid type for graphic parameter lab.width: must be numeric and positive");
		}
		if(params[["lab.width"]] <= 0){
			stop("invalid value for graphic parameter lab.width: must be numeric and positive");
		}
		else{		result.list[["lab.width"]] = lcm(params[["lab.width"]]);		}
	}
	else{ result.list[["lab.width"]] = lcm(1); }	# lab.width defaults to lcm(1), i.e. 1 cm
	params[["lab.width"]] = NULL;	

	# ---------------------------------------------------------------- params[["lab.textpar"]]

	if(!is.null(params[["lab.textpar"]])){
		if(typeof(params[["lab.textpar"]]) != "list"){
			stop("invalid type for graphic parameter lab.textpar: must be a tagged list");
		}
	}
	else{ params[["lab.textpar"]] = list(); }
	if(is.null(params[["lab.textpar"]][["cex"]])){		params[["lab.textpar"]][["cex"]] = 1;		   	}
	if(is.null(params[["lab.textpar"]][["x"]]))	 {    params[["lab.textpar"]][["x"]] = 0.5; 			}
	if(is.null(params[["lab.textpar"]][["y"]]))	 { 		params[["lab.textpar"]][["y"]] = 0.5;				}

	result.list[["lab.textpar"]] = params[["lab.textpar"]];
	params[["lab.textpar"]] = NULL;

	# ---------------------------------------------------------------- params[["levels.loc"]]
	
	if(!is.null(params[["levels.loc"]])){
		if( ifelse(inner.variable == "x", params[["levels.loc"]] != "bottom" && params[["levels.loc"]] != "all" && params[["levels.loc"]] != "none", 
                                      params[["levels.loc"]] != "left" && params[["levels.loc"]] != "all" && params[["levels.loc"]] != "none") ){
       stop("invalid value for graphic parameter levels.loc: must be either \"left\" (for inner Y var), \"bottom\" (for inner X var) , \"all\" or \"none\"");
		}
	}
	else{ 
		if(inner.variable == "x"){		params[["levels.loc"]] = "bottom"; }		# Inner "x" defaults to "bottom": show the inner x levels only in the heatmaps of the bottom row
		else{													params[["levels.loc"]] = "left"; 	}		# Inner "y" defaults to "left": show the inner y levels only in the heatmaps of the left-most column
	}
	result.list[["levels.loc"]] = params[["levels.loc"]];
	params[["levels.loc"]] = NULL;
	
	# ---------------------------------------------------------------- params[["levels.at"]]
	
	if(!is.null(params[["levels.at"]])){
		if(!is(params[["levels.at"]], "character")){
			stop("invalid value for graphic parameter levels.at: must be a vector of text strings");
		}
	}
	# Ignore the user-specified value for this parameter if params[["levels.loc"]] == "none"
	if(result.list[["levels.loc"]] == "none"){  params[["levels.at"]] = c("");  }	
	result.list[["levels.at"]] = params[["levels.at"]];
	params[["levels.at"]] = NULL;

	# ---------------------------------------------------------------- params[["levels.las"]]
	
	if(!is.null(params[["levels.las"]])){
		if(!(params[["levels.las"]] %in% c(1,2,3,4))){
			stop("invalid value for graphic parameter levels.las: must be an integer of 1, 2, 3 or 4");
		}
		params[["levels.las"]] = strtoi(params[["levels.las"]]);
		if(is.na(params[["levels.las"]])){ 
      stop("invalid value for graphic parameter levels.las: must be an integer of 1, 2, 3 or 4");
		}
	}
	else{	params[["levels.las"]] = 1; 	}	
	result.list[["levels.las"]] = params[["levels.las"]];
	params[["levels.las"]] = NULL;
		
	# ----------------------------------------------------------------	
	## Check for unrecognized tags in the list
	
	if(sum(!sapply(params, is.null)) > 0){
		nms = names(params);
		stop(paste("unrecognized graphical parameter:",nms[which.max(!sapply(params, is.null))]));
	}

	return(result.list);		
}

# ___________________________________________________________________________

##  ===========================================================================
## 	Possible graphical options for the colorbar
##  ===========================================================================
##
##	levels.at = NULL,												# Levels at which the Y axis ticks of the colorbar will be shown. Defaults to 0, the minimum value and the maximum value.
##	hlines = TRUE,													# Whether black horizontal lines should be displayed in the colorbar to separate the colors. Defaults to TRUE

.check.colorbar.params <- function(params){
	
	# Check for NAs	
	if(length(params) > 0){		
		if(sum(is.na(params)) > 0){
			nms = names(params);
			stop(paste("NAs not allowed for any parameter. Invalid value (NA) for parameter", nms[which.max(is.na(params))]));
		}
	}

	result.list = list();
	
	# ---------------------------------------------------------------- params[["levels.at"]]	
	
	if(!is.null(params[["levels.at"]])){
		if(!is(params[["levels.at"]], "character")){
			stop("invalid value for graphic parameter levels.at: must be a vector of text strings");
		}
	}
	result.list[["levels.at"]] = params[["levels.at"]];
	params[["levels.at"]] = NULL;
	
	# ---------------------------------------------------------------- params[["levels.at"]]		
	
	if(!is.null(params[["hlines"]])){ 
		if(typeof(params[["hlines"]]) != "logical"){
			stop("invalid value for graphic parameter hlines: must be either TRUE or FALSE");
		}
	}
	else{ params[["hlines"]] = TRUE; }
	result.list[["hlines"]] = params[["hlines"]];
	params[["hlines"]] = NULL;

	# ----------------------------------------------------------------	
	## Check for unrecognized tags in the list
	
	if(sum(!sapply(params, is.null)) > 0){
		nms = names(params);
		stop(paste("unrecognized graphical parameter:",nms[which.max(!sapply(params, is.null))]));
	}
	
	return(result.list);
}	

# ___________________________________________________________________________

##  ===========================================================================
## 	Possible graphical options for the singleplot labels
##  ===========================================================================
##
##  xlab	= TRUE							# - Label with the name of the variable. Will be displayed horizontally. Valid values:
##														# 	FALSE for no label; NULL or TRUE for the name of the outer variable; and any string for a specific label. Defaults to TRUE
##  ylab	= TRUE							# - Label with the name of the variable. Will be displayed vertically. Valid values:
##														# 	FALSE for no label; NULL or TRUE for the name of the outer variable; and any string for a specific label. Defaults to TRUE
##	xlevels.at = NULL,				# levels of the inner variable where the label will be shown. Defaults to all the levels.
##														# They can be provided in any order, since the order in which they will be depicted only depends on the
##														# original order defined by the levels argument when the corresponding factor column of the data was created
##	ylevels.at = NULL,				# levels of the inner variable where the label will be shown. Defaults to all the levels.
##														# They can be provided in any order, since the order in which they will be depicted only depends on the
##														# original order defined by the levels argument when the corresponding factor column of the data was created
.check.labels.params <- function(params){
	
	# Check for NAs	
	if(length(params) > 0){		
		if(sum(is.na(params)) > 0){
			nms = names(params);
			stop(paste("NAs not allowed for any parameter. Invalid value (NA) for parameter", nms[which.max(is.na(params))]));
		}
	}

	result.list = list();

	# ---------------------------------------------------------------- params[["xlab"]]
	if(!is.null(params[["xlab"]])){
		if( typeof(params[["xlab"]]) != "character" && typeof(params[["xlab"]]) != "logical" ){
			stop("invalid type for graphic parameter lab: must be either a character string or TRUE/FALSE");
		}
	}
	else params[["xlab"]] = TRUE;				# defaults to TRUE (display the name of the variable passed in the xOuter or yOuter arguments)	
	result.list[["xlab"]] = params[["xlab"]];
	params[["xlab"]] = NULL;

	# ---------------------------------------------------------------- params[["ylab"]]
	
	if(!is.null(params[["ylab"]])){
		if( typeof(params[["ylab"]]) != "character" && typeof(params[["ylab"]]) != "logical" ){
			stop("invalid type for graphic parameter lab: must be either a character string or TRUE/FALSE");
		}
	}
	else params[["ylab"]] = TRUE;				# defaults to TRUE (display the name of the variable passed in the xOuter or yOuter arguments)	
	result.list[["ylab"]] = params[["ylab"]];
	params[["ylab"]] = NULL;

	# ---------------------------------------------------------------- params[["xlevels.at"]]
	
	if(!is.null(params[["xlevels.at"]])){
		if(!is(params[["xlevels.at"]], "character")){
			stop("invalid value for graphic parameter xlevels.at: must be a vector of text strings");
		}
	}
	result.list[["xlevels.at"]] = params[["xlevels.at"]];
	params[["xlevels.at"]] = NULL;

	# ---------------------------------------------------------------- params[["ylevels.at"]]
	
	if(!is.null(params[["ylevels.at"]])){
		if(!is(params[["ylevels.at"]], "character")){
			stop("invalid value for graphic parameter ylevels.at: must be a vector of text strings");
		}
	}
	result.list[["ylevels.at"]] = params[["ylevels.at"]];
	params[["ylevels.at"]] = NULL;
	
	return(result.list);
}