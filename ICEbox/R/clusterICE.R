clusterICE = function(ice_obj, nClusters, plot = TRUE, plot_margin = 0.05, colorvec, plot_pdp = FALSE,
			x_quantile = FALSE, avg_lwd = 3, centered = FALSE, plot_legend = FALSE, ...){

	DEFAULT_COLORVEC = c("green", "red", "blue", "black", "green", "yellow", "pink", "orange", "forestgreen", "grey")
	
	if(class(ice_obj) != "ice" && class(ice_obj) != "dice"){
		stop("'ice_obj' must be of class 'ice' or 'dice'.")
	}

	objclass = class(ice_obj)
	
	if (missing(nClusters) || !(nClusters %% 1 == 0 ) || (nClusters <= 0)){
		stop("nClusters must be a positive integer.")
	}

	#make all curves have avg value = 0.
	if(objclass == "ice"){
		curves = t(scale(t(ice_obj$ice_curves), center = T, scale = F))  #sum(ice_curves[i,]) = 0 for all i
	}
	else{
		curves = t(scale(t(ice_obj$d_ice_curves), center = T, scale = F))
	}

	#cluster
	cl = kmeans(curves, iter.max = 20, centers = nClusters, ...)
	if (missing(colorvec)){
		colorvec = DEFAULT_COLORVEC
		if(length(colorvec) < nClusters){
			colorvec = c(colorvec, rgb(runif(nClusters - 10, 0, 0.7), runif(nClusters - 10, 0, 0.7), runif(nClusters - 10, 0, 0.7)))
		}
	}
	
	cluster_centers = cl$centers
	
	if (centered){
		for (k in 1 : nrow(cluster_centers)){
			cluster_centers[k, ] = cluster_centers[k, ] - cluster_centers[k, 1] #ifelse(cluster_centers[o, 1] < 0, -cluster_centers[o, 1], cluster_centers[o, 1]
		}		
	}
		
	if (plot){
		#y limits
		rg = range(cluster_centers) 
		dist = rg[2] - rg[1]
		rg_min = rg[1] - plot_margin * dist
		rg_max = rg[2] + plot_margin * dist

		#x grid and xlab
		xlab = ice_obj$xlab
		grid = ice_obj$gridpts 
		if (x_quantile){
			xlab = paste("quantile(", xlab, ")", sep = "")
			ecdf_fcn = ecdf(grid)
			grid = ecdf_fcn(grid)
		}
		plot(grid, as.numeric(cluster_centers[1, ]), 
				type = 'n', 
				ylim = c(rg_min, rg_max), 
				ylab = paste("cluster yhat"), 
				xlab = xlab, 
				xaxt = ifelse(ice_obj$nominal_axis, "n", "s"))

		cluster_order = order(cl$centers[, 1]) #use original, non-centered object only
		cluster_size = cl$size / sum(cl$size)
		total_line_width = avg_lwd * nClusters
		
		centers_to_plot = array(NA, nrow(cluster_centers))
		for(i in 1 : nrow(cluster_centers)){		
			#we re-order it so that when the code is rerun, randomness in kmeans
			#doesn't switch which cluster goes with which color.
			center_to_plot = cluster_order[i]
			centers_to_plot[i] = center_to_plot
			points(grid, cluster_centers[center_to_plot, ], col = colorvec[i], type = "l", 
					lwd = cluster_size[center_to_plot] * total_line_width)
		}
				
		if(plot_pdp){  
			# this is done after all other plotting so nothing obfuscates the PDP
			if(objclass == "ice"){
				pdp = ice_obj$pdp
			}
			else{
				pdp = ice_obj$dpdp #really want the dpdp
			}
			#now mean center it:
			pdp = pdp - mean(pdp)

			#if centered=TRUE, start everything at 0.
			if(centered){
				pdp = pdp - pdp[1]
			}		

			#calculate the line thickness based on how many lines there are
			#every add'l cluster we get 0.5 more highlight up to 8
			points(grid, pdp, col = "yellow", type = "l", lwd = min(4 +  nClusters*0.5, 5),8)  #yellow
			points(grid, pdp, col = "BLACK", type = "l", lwd = 3, lty=3)  #black, dotted.
		}#end of pdp

		#if (prop_range_y){
		#	at = seq(min(curves), max(curves), length.out = 5)
		#	#we need to organize it so it's at zero
		#	at = at - min(abs(at))
		#	
		#	labels = round(at / ice_obj$range_y, 2)
		#	axis(4, at = at, labels = labels)
		#}
		
		if (ice_obj$nominal_axis){
			axis(1, at = sort(ice_obj$xj), labels = sort(ice_obj$xj))
		}

	
		#legend
		if (plot_legend){
			#fix ordering
			prop_data_in_clusters = round(cl$size[cluster_order] / sum(cl$size), 2)
			
			legend("topleft", title = "Prop.", inset = 0.01, 
				legend = as.character(prop_data_in_clusters), fill = colorvec, cex = 0.8)
		}
	}

	invisible(cl)
}
