plot.ice = function(x, plot_margin = 0.05, frac_to_plot = 1, plot_points_indices = NULL,
					plot_orig_pts_preds = TRUE, pts_preds_size = 1.5,
					colorvec, color_by = NULL, x_quantile = FALSE, plot_pdp = TRUE, centered = FALSE, 
					prop_range_y = TRUE, rug_quantile = seq(from = 0, to = 1, by = 0.1), 
					centered_percentile = 0.01, point_labels = NULL, point_labels_size = NULL,
					prop_type = "sd", ...){
	
	DEFAULT_COLORVEC = c("firebrick3", "dodgerblue3", "gold1", "darkorchid4", "orange4", "forestgreen", "grey", "black")
	#think of x as x. needs to be 'x' to match R's generic.

	#list of passed arguments, including the ...
	#arg_list = as.list(match.call(expand.dots = TRUE))

	#some argument checking
	if (class(x) != "ice"){ 
		stop("object is not of class \"ice\"")
	}
	if (frac_to_plot <= 0 || frac_to_plot > 1 ){
		stop("frac_to_plot must be in (0,1]")
	}
	if(!(prop_type %in% c("sd","range"))){
		stop("prop_type must be either 'sd' or 'range'")
	}

	#extract the grid and lines to plot	
	grid = x$gridpts 
	n_grid = length(grid)
	ecdf_fcn = NULL
	if (x_quantile){
		ecdf_fcn = ecdf(grid)
		grid = ecdf_fcn(grid)
	}
	ice_curves = x$ice_curves
	N = nrow(ice_curves)
	
	if (!is.null(point_labels)){
		if (length(point_labels) != N){
			stop("point_labels must be same length as number of ICE curves: ", N)
		}
	}

	#### figure out the colorvec.
	legend_text = NULL #default is no legend.
	#case 1: random
	if (missing(colorvec) && missing(color_by)){
		#we're going to choose dark grey and randomly alpha the lines
		colorvec = sort(rgb(rep(0.4, N), rep(0.4, N), rep(0.4, N), runif(N, 0.4, 0.8)))
	} 
	#case 2: both colorvec and color_by specified, so print a warning but use colorvec.
	if (!missing(colorvec) && !missing(color_by)){
		if (!missing(colorvec) && length(colorvec) < N){
			stop("color vector has length ", length(colorvec), " but there are ", N, " lines to plot")
		}
#		warning("Both colorvec and color_by_predictor are specified...using colorvec.")
	}	
	#case 3: colorvec missing but color_by is specified.
	if (!missing(color_by) && missing(colorvec)){
		#argument checking first:
		arg_type = class(color_by)
		if(!(arg_type %in% c("character", "numeric", "factor"))){
			stop("color_by must be a column name in X or a column index")
		}
		if(class(color_by) == "character"){
			if(!(color_by %in% names(x$Xice))){
				stop("The predictor name given by color_by was not found in the X matrix")
			}
			x_color_by = x$Xice[, color_by]
		} else if (length(color_by) > N){ #tell the user the thing they passed in doesn't line up
			stop("The color_by_data vector you passed in has ", length(color_by), " entries but the ICEbox object only has ", N, " curves.")
		} else if (length(color_by) == N){ #it's an actual data vector
			x_color_by = color_by
		}		
		else{  #check numeric
			if( color_by < 1 || color_by > ncol(x$Xice) || (color_by%%1 !=0)){
				stop("color_by must be a column name in X or a column index")
			}
			x_color_by = x$Xice[, color_by]
		}
		x_unique = unique(x_color_by)
		num_x_color_by = length(x_unique)		
		
		
		#if there are 10 or fewer unique values of this x value, we use the
		#same color in DEFAULT_COLORVEC for each. Otherwise, we use a rainbow.
		if (num_x_color_by <= 10){
			
			which_category = match(x_color_by, x_unique)
			colorvec = DEFAULT_COLORVEC[which_category]
			
			#now make the legend.
			legend_text = as.data.frame(cbind(x_unique, DEFAULT_COLORVEC[1 : num_x_color_by]))
			x_column_name = ifelse(length(color_by) == N, "data vector level", ifelse(is.character(color_by), color_by, paste("x_", color_by, sep = "")))
			names(legend_text) = c(x_column_name,"color")
			cat("ICE Plot Color Legend\n")
			print(legend_text, row.names = FALSE)			
		} else {
			if (is.factor(x_color_by)){
				warning("color_by is a factor with greater than 10 levels: coercing to numeric.")
				x_color_by = as.numeric(x_color_by)
			}			

			alpha_blend_colors = matrix(0, nrow = N, ncol = 3)
			
			alpha_blend_colors[, 1] = seq(from = 1, to = 0, length.out = N)
			alpha_blend_colors[, 2] = seq(from = 0, to = 1, length.out = N)
			alpha_blend_colors[, 3] = 0
			
			
			rgbs = array(NA, N)
			for (i in 1 : N){
				rgbs[i] = rgb(alpha_blend_colors[i, 1], alpha_blend_colors[i, 2], alpha_blend_colors[i, 3])
			}
			
			colorvec = rgbs[sort(x_color_by, index.return = T)$ix]
			#now make the legend.
			cat("ICE Plot Color Legend: red = low values of the color_by variable and green = high values\n")
		}
	}

	
	#pull out a fraction of the lines to plot
	if (is.null(plot_points_indices)){
		plot_points_indices = which(as.logical(sample(1 : N, round(frac_to_plot * N))))
	} else {
		if (frac_to_plot < 1){
			stop("frac_to_plot has to be 1 when plot_points_indices is passed to the plot function.")
		}
	}
	
	ice_curves = ice_curves[plot_points_indices, ]
	if (nrow(ice_curves) == 0){
		stop("no rows selected: frac_to_plot too small.")
	}
	if (centered){
		centering_vector = ice_curves[, ceiling(ncol(ice_curves) * centered_percentile + 0.00001)]
		ice_curves = ice_curves - centering_vector
	}
	colorvec = colorvec[plot_points_indices]
	
	##### now start plotting
	min_ice_curves = min(ice_curves)
	max_ice_curves = max(ice_curves)
	range_ice_curves = max_ice_curves - min_ice_curves
	min_ice_curves = min_ice_curves - plot_margin * range_ice_curves
	max_ice_curves = max_ice_curves + plot_margin * range_ice_curves

  
  arg_list = list(...)
  #add the x and y values
  arg_list = modifyList(arg_list, list(x = grid, y = ice_curves[1, ]))
  
	#get the xlabel if it wasn't already passed explicitly.
	if( is.null(arg_list$xlab)){
		xlab = x$xlab
    	arg_list = modifyList(arg_list, list(xlab = xlab))
	}
	if (x_quantile){
		xlab = paste("quantile(", xlab, ")", sep = "")
		arg_list = modifyList(arg_list, list(xlab = xlab))
	}
	if (!missing(color_by)){
		xlab = paste(xlab, "colored by", ifelse(length(color_by) == N, "a provided data vector", color_by))
		arg_list = modifyList(arg_list, list(xlab = xlab))
	}
	
	#same for y label
	if( is.null(arg_list$ylab)){	
		if (x$logodds){
			ylab = "partial log-odds"
			arg_list = modifyList(arg_list, list(ylab = ylab))
		} else if(x$probit){
			ylab = "partial probit"
			arg_list = modifyList(arg_list, list(ylab = ylab))		
		}else {
			ylab = paste("partial yhat", ifelse(centered, "(centered)", ""))
			arg_list = modifyList(arg_list, list(ylab = ylab))
		}
	}

	#set xact if not passed explicitly 
	if( is.null(arg_list$xaxt) ){
		xaxt = ifelse(x$nominal_axis, "n", "s")
		arg_list = modifyList(arg_list, list(xaxt = xaxt))
	}

	#set ylim if not passed explicitly
	if (is.null(arg_list$ylim)){
		ylim = c(min_ice_curves, max_ice_curves) 
		arg_list = modifyList(arg_list, list(ylim = ylim))
	}
	#set type if not passed explicitly
	if (is.null(arg_list$type)){
		type = "n"
		arg_list = modifyList(arg_list, list(type = type))
	}

  
	
	#plot all the prediction lines
	do.call("plot", arg_list)
  
  
	if (x$nominal_axis){
		axis(1, at = sort(x$xj), labels = sort(x$xj), cex.axis = arg_list$cex.axis)
	}	
	if (centered && prop_range_y && !x$logodds && !x$probit){ #don't draw this axis for logodds since it makes no sense
		at = seq(min(ice_curves), max(ice_curves), length.out = 5)
		#we need to organize it so it's at zero
		at = at - min(abs(at))
		
		#check prop type.
		if(prop_type == "range"){
			labels = round(at / x$range_y, 2)  #as a fraction of range of y
		}else{
			labels = round(at / x$sd_y, 2)     #as a fraction of sd(y)
		}
		axis(4, at = at, labels = labels, cex.axis = arg_list$cex.axis)
	}
	
	for (i in 1 : nrow(ice_curves)){
		points(grid, ice_curves[i, ], col = colorvec[i], type = "l")
	}

	if (plot_orig_pts_preds){ #indicate the fitted values associated with observed xj values
		yhat_actual = x$actual_prediction[plot_points_indices]
		if (centered){
			yhat_actual = yhat_actual - centering_vector
		}
				
		if (x_quantile){
			xj = ecdf_fcn(x$xj)[plot_points_indices]
		} else {
			xj = x$xj[plot_points_indices]
		}
		for (i in 1 : length(xj)){
			points(xj[i], yhat_actual[i], col = rgb(0.1, 0.1, 0.1), pch = 16, cex = pts_preds_size)
			points(xj[i], yhat_actual[i], col = colorvec[i], pch = 16, cex = round(pts_preds_size * 0.7))
		}
	}
	
	if (!is.null(point_labels)){
		text(xj, yhat_actual, 
				pos = 4,
				labels = point_labels[plot_points_indices], 
				cex = ifelse(is.null(point_labels_size), pts_preds_size, point_labels_size))
	}
	
	if (!is.null(rug_quantile) && !x_quantile){
		axis(side = 1, line = -0.1, at = quantile(x$xj, rug_quantile), lwd = 0, tick = T, tcl = 0.4, lwd.ticks = 2, col.ticks = "blue4", labels = FALSE, cex.axis = arg_list$cex.axis)
	}
	
	#if plot_pdp is true, plot actual pdp (in the sense of Friedman '01)
	#Ensure this is done after all other plotting so nothing obfuscates the PDP
	if (plot_pdp){
		pdp = apply(ice_curves, 2, mean) # pdp = average over the columns (we don't use the one from the ICE object since plot_points_indices may have been passed)
		#cat("pdp has", nrow(ice_curves), "rows\n")
		if (centered){
#			ice_curves[, ceiling(ncol(ice_curves) * centered_percentile + 0.00001)]
			pdp = pdp - pdp[ceiling(length(pdp) * centered_percentile + 0.00001)]
		}		

		#calculate the line thickness based on how many lines there are
		num_lines = length(plot_points_indices)
		points(grid, pdp, col = "yellow", type = "l", lwd = min(5.5 + (num_lines / 100) * 0.75, 8)) #every 100 lines we get 0.5 more highlight up to 8
		points(grid, pdp, col = "BLACK", type = "l", lwd = 4)
	}
		
	invisible(list(plot_points_indices = plot_points_indices, legend_text = legend_text, pdp = pdp))
}

