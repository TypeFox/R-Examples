plot.dice = function(x, plot_margin = 0.05, frac_to_plot = 1, plot_sd = TRUE, plot_orig_pts_deriv = TRUE,
 						pts_preds_size = 1.5, colorvec, color_by = NULL, x_quantile = FALSE, plot_dpdp = TRUE, 
						rug_quantile = seq(from = 0, to = 1, by = 0.1), ...){
	#think of x as 'dice_obj'
	
	DEFAULT_COLORVEC = c("firebrick3", "dodgerblue3", "gold1", "darkorchid4", "orange4", "forestgreen", "grey", "black")
	
	arg_list = list(...)

	#list of passed arguments, including the ...
#	arg_list = as.list(match.call(expand.dots = TRUE))

	#some argument checking
	if (class(x) != "dice"){ 
		stop("object is not of class 'dice'")
	}
	if (frac_to_plot <= 0 || frac_to_plot > 1 ){
		stop("frac_to_plot must be in (0,1]")
	}
	if(!is.null(arg_list$ylim) && plot_sd == TRUE){
		stop("Cannot specify both ylim and plot_sd=TRUE.")
	}
	
	#extract the grid and lines to plot	
	grid = x$gridpts 
	n_grid = length(grid)
	ecdf_fcn = NULL
	if (x_quantile){
		ecdf_fcn = ecdf(grid)
		grid = ecdf_fcn(grid)
	}
	d_ice_curves = x$d_ice_curves
	N = nrow(d_ice_curves)

	#### figure out the colorvec.
	legend_text = NULL #default is no legend.
	#case 1: random
	if (missing(colorvec) && missing(color_by)){
		colorvec = sort(rgb(rep(0.4, N), rep(0.4, N), rep(0.4, N), runif(N, 0.4, 0.8)))
	} 
	#case 2: both colorvec and color_by specified, so print a warning but use colorvec.
	if(!missing(colorvec) && !missing(color_by)){
		if (!missing(colorvec) && length(colorvec) < N){
			stop("color vector has length ", length(colorvec), " but there are ", N, " lines to plot")
		}
#		warning("Both colorvec and color_by_predictor are specified...using colorvec.")
	}	
	#case 3: colorvec missing but color_by is specified.
	if(!missing(color_by) && missing(colorvec)){
		#argument checking first:
		arg_type = class(color_by)
		if(!(arg_type %in% c("character", "numeric"))){
			stop("color_by must be a column name in X or a column index")
		}
		if(class(color_by) == "character"){
			if(!(color_by %in% names(x$Xice))){
				stop("The predictor name given by color_by was not found in the X matrix")
			}
		} else{  #check numeric
			if( color_by < 1 || color_by > ncol(x$Xice) || (color_by%%1 !=0)){
				stop("color_by must be a column name in X or a column index")
			}
		}
		x_color_by = x$Xice[, color_by]
		x_unique = unique(x_color_by)
		num_x_color_by = length(x_unique)		
		
		
		#if there are 10 or fewer unique values of this x value, we use the
		#same color in DEFAULT_COLORVEC for each. Otherwise, we use a rainbow.
		if (num_x_color_by <= 10){
			
			which_category = match(x_color_by, x_unique)
			colorvec = DEFAULT_COLORVEC[which_category]
			
			#now make the legend.
			legend_text = as.data.frame(cbind(x_unique, DEFAULT_COLORVEC[1 : num_x_color_by]))
			x_column_name = ifelse(is.character(color_by), color_by, paste("x_", color_by, sep = ""))
			names(legend_text) = c(x_column_name,"color")
			cat("dICE Color Legend\n")
			print(legend_text)			
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
			
#			plot(1:200,1:200, type = "n")
#			for (i in 1 : N){
#				abline(a = i, b = 0, col = rgbs[i])
#			}
			colorvec = rgbs[sort(x_color_by, index.return = T)$ix]
		}
	}

	
	#pull out a fraction of the lines to plot
	plot_points_indices = which(as.logical(rbinom(N, 1, frac_to_plot)))
	d_ice_curves = d_ice_curves[plot_points_indices, ]
	if (nrow(d_ice_curves) == 0){
		stop("no rows selected: frac_to_plot too small.")
	}
	colorvec = colorvec[plot_points_indices]
	
	##### now start plotting
	min_dice = min(d_ice_curves)
	max_dice = max(d_ice_curves)
	range_dice = max_dice - min_dice
	min_dice = min_dice - plot_margin * range_dice
	max_dice = max_dice + plot_margin * range_dice

   #add the x and y values
   arg_list = modifyList(arg_list, list(x = grid, y = d_ice_curves[1, ]))
  
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
		xlab = paste(xlab, "colored by", color_by)
		arg_list = modifyList(arg_list, list(xlab = xlab))
	}
	
	#same for y label
	if( is.null(arg_list$ylab)){	
		if (x$logodds){
			ylab = "partial log-odds"
			arg_list = modifyList(arg_list, list(ylab = ylab))
		} else {
			ylab = "derivative f-hat"
			arg_list = modifyList(arg_list, list(ylab = ylab))
		}
	}

	#set xact if not passed explicitly 
	if( is.null(arg_list$xaxt) ){
		xaxt = ifelse(x$nominal_axis, "n", "s")
		arg_list = modifyList(arg_list, list(xaxt = xaxt))
	}

	#set ylim if not passed explicitly
	if( is.null(arg_list$ylim) ){
		if(plot_sd){
			offset = 1.5 * max(x$sd_deriv)
			ylim = c(min_dice - offset, max_dice)	
		}else{
			ylim = c(min_dice, max_dice) 
		}
		arg_list = modifyList(arg_list, list(ylim = ylim))
	}
	#set type if not passed explicitly
	if( is.null(arg_list$type) ){
		type = "n"
		arg_list = modifyList(arg_list, list(type = type))
	}


	## if plot_sd = TRUE, set up the layout to have
    ## the dpdp above and the sd plot below.
	do.call("plot", arg_list)
  
  
	if (x$nominal_axis){
		axis(1, at = sort(x$xj), labels = sort(x$xj), cex.axis = arg_list$cex.axis)
	}	
	
	for (i in 1 : nrow(d_ice_curves)){
		points(grid, d_ice_curves[i, ], col = colorvec[i], type = "l")
	}

	if (plot_orig_pts_deriv){ #indicate the fitted values associated with observed xj values
		deriv_actual = x$actual_deriv[plot_points_indices]
				
		if (x_quantile){
			xj = ecdf_fcn(x$xj)[plot_points_indices]
		} else {
			xj = x$xj[plot_points_indices]
		}
		for (i in 1 : length(xj)){
			points(xj[i], deriv_actual[i], col = "black", pch = 16, cex = pts_preds_size)
			points(xj[i], deriv_actual[i], col = colorvec[i], pch = 16)
		}
	}


	#if plot_dpdp is true, plot actual dpdp (in the sense of Friedman '01)
	if (plot_dpdp){
		friedman_dpdp = x$dpdp
		
		#calculate the line thickness based on how many lines there are
		num_lines = length(plot_points_indices)
		#every 100 lines we get 0.5 more highlight up to 8
		points(grid, friedman_dpdp, col = "yellow", type = "l", lwd = min(5.5 + (num_lines / 100) * 0.75, 8)) 
		points(grid, friedman_dpdp, col = "BLACK", type = "l", lwd = 4)
	}
	
	if (!is.null(rug_quantile) && !x_quantile){
		axis(side = 1, line = -0.1, at = quantile(x$xj, rug_quantile), lwd = 0, tick = T, tcl = 0.4, lwd.ticks = 2, col.ticks = "blue4", labels = FALSE,
		     cex.axis = arg_list$cex.axis)
	}

	#do the sd plot if required.
	if (plot_sd){
		abline(h = ylim[1] + offset, col = rgb(0.8,0.8,0.8))
		at = seq(ylim[1], ylim[1] + max(x$sd_deriv), length.out = 2)	

		labels = round(seq(0, max(x$sd_deriv), length.out = 2), 1)
		axis(4, at = at, labels = labels, cex.axis = arg_list$cex.axis)
    print(arg_list$cex)
		mtext("sd(deriv)", side = 4,line = 0.5, cex = arg_list$cex.lab)

		points(x= grid, y = (x$sd_deriv+ylim[1]),type='l')
	}
		
	if (is.null(legend_text)){
		invisible(list(plot_points_indices = plot_points_indices, legend_text = legend_text))
	} else {
		invisible(list(plot_points_indices = plot_points_indices, legend_text = legend_text))
	}
}
