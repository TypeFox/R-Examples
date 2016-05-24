autoKrige.cv = function(formula, input_data, data_variogram = input_data,
                          model = c("Sph", "Exp", "Gau", "Ste"), kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), 
						  fix.values = c(NA,NA,NA), verbose = c(FALSE, TRUE), GLS.model = NA,
                          start_vals = c(NA,NA,NA), miscFitOptions = list(),...)
# Automatically fits a variogram model to the data using autofitVariogram and performs
# crossvalidation by calling krige.cv
{
    if(length(verbose) != 2) stop("Argument 'verbose' should be of length 2, e.g. c(TRUE,FALSE).")
	kr.cv = krige.cv(formula, input_data, 
                      model = autofitVariogram(formula, 
                                  data_variogram, 
                                  model = model, 
                                  kappa = kappa,
                                  fix.values = fix.values,
                                  verbose = verbose[1],
                                  GLS.model = GLS.model,
                                  start_vals = start_vals,
                                  miscFitOptions = miscFitOptions)$var_model,
                         verbose = verbose[2],
						 ...)
	kr.cv = list(krige.cv_output = kr.cv)
	class(kr.cv) = "autoKrige.cv"
	return(kr.cv)
}

summary.autoKrige.cv = function(object, ..., digits = 4) {
# The summary function for the autoKrige.cv object. Returns some
# statistics about the cross-validation
	obj = object$krige.cv_output
	out = list()
	# mean error, ideally 0:
	out$mean_error = mean(obj$residual)
    # mean error divided by the mean of the observed values, measure for how large the mean_error is in contrast to the mean of the dataset
    out$me_mean = out$mean_error / mean(obj$observed)
	# mean absolute error, ideally 0, less vulnerable to outliers
	out$MAE = mean(abs(obj$residual))
	# MSE, ideally small
	out$MSE = mean(obj$residual^2)
	# Mean square normalized error, ideally close to 1
    out$MSNE = mean(obj$zscore^2)
	# correlation observed and predicted, ideally 1
    out$cor_obspred = cor(obj$observed, obj$observed - obj$residual)
	# correlation predicted and residual, ideally 0
    out$cor_predres = cor(obj$observed - obj$residual, obj$residual)
	# RMSE, ideally small
	out$RMSE = sqrt(sum(obj$residual^2) / length(obj$residual))
    # RMSE / sd(observed), measure for how much the residuals vary to the total variation in the dataset
    out$RMSE_sd = out$RMSE / sd(obj$observed)
	# URMSE, ideally zero
	out$URMSE = sqrt((sum(obj$residual^2) / length(obj$residual)) - mean(obj$residual)^2)
	# Inter quartile range, ideally small
	out$iqr = IQR(obj$residual)

    out = lapply(out, signif, digits = digits)
    out = t(t(out))
	return(out)
}

compare.cv = function(..., col.names, bubbleplots = FALSE, zcol = "residual", 
						   layout, key.entries, reference = 1, plot.diff = FALSE,
                           digits = 4, ggplot = FALSE, addPoly = NULL) 
# A function to compare cross-validations to each other in both statistics (using summary.autoKrige.cv) or
# in bubble plots (using cv.compare.bubble). '...' can be both output from krige.cv or autoKrige.cv.
{
	dots = list(...)
    
    # If dots consists of only one entry, which is a list, that list is assigned to dots
    if(length(dots) == 1 & inherits(dots, "list")) {
      dots = dots[[1]]
      if(missing(col.names)) col.names = names(dots)
    }

    # If a user passes output directly from krige.cv, change it to an autoKrige.cv object
    dots = checkIfautokrige.cv(dots)

	out = do.call("cbind", lapply(dots, summary, digits = digits))
	if(missing(col.names)) {
        col.names = as.character(as.list(match.call(expand.dots=TRUE))[-1])[1:length(dots)]
	}
	dimnames(out)[[2]] = col.names

	if(bubbleplots) {
		if(missing(layout)) {
			layout1 = unlist(sapply(1:10, function(x) return(rep(x, 2*x - 1))))
			layout2 = unlist(sapply(1:10, function(x) return(rep(x, 2*x))))
			layout = c(layout1[length(dots)], layout2[length(dots)])
		}
        if(ggplot) {
            print(cv.compare.ggplot(dots, zcol = zcol, layout = layout, addPoly = addPoly, col.names = col.names))
        } else {
            if(missing(key.entries)) key.entries = quantile(dots[[reference]]$krige.cv_output[[zcol]])
            cv.compare.bubble(dots, zcol = zcol, col.names = col.names, 
                              layout = layout, key.entries = key.entries, 
                              reference = reference, plot.diff = plot.diff)
        }
	}	

	return(data.frame(out))
}

# Checks wheter the list contains autoKrige.cv objects,
# if they are spatial, convert to autoKrige.cv object.
checkIfautokrige.cv = function(l) {
  l = lapply(l, function(x) { 
          if(inherits(x, "autoKrige.cv")) {                   # autoKrige.cv output
              return(x)
          } else if(inherits(x, "SpatialPointsDataFrame")){   # krige.cv output
              x = list(krige.cv_output = x)
              class(x) = "autoKrige.cv"
              return(x)
          } else {                                            # unknown output
              stop(sprintf("One of the objects in \'...\' has class \'%s\'. Class of objects in \'...\' should be one of \n  \'autoKrige.cv\' or \'SpatialPointsDataFrame\'", class(x)[1]))
          }
      })
  return(l)
}


cv.compare.bubble = function(objs, zcol, key.entries, layout, col.names, plot.diff, reference,...) 
# A function to create bubble plots of multiple cross validation in order to compare them.
{
  ref_data = objs[[reference]]$krige.cv_output[[zcol]]
  plot.column = zcol
  plot_list = list()
  it = 1
  for(obj in objs) {
	main = col.names[it]
    if(plot.diff) {
		obj$krige.cv_output$diff = ref_data^2 - obj$krige.cv_output[[zcol]]^2
		plot.column = "diff"
		if(it == reference) main = "Reference" else main = as.expression(parse(text = paste("Ref^2 -", col.names[it],"^2", sep = "")))
	}
	plot_list[[it]] = bubble(obj$krige.cv_output, 
    					zcol = plot.column, 
						key.entries = key.entries,
						main = main,
						sub = zcol,
						...)
    it = it + 1     	  
    }

  pos_list = cbind(rep(0:(layout[1]-1) / layout[1], layout[2]),
		   rep((layout[2]-1):0/layout[2], each = layout[1]),
		   rep(1:layout[1] / layout[1], layout[2]),
		   rep(layout[2]:1/layout[2],each = layout[1]))
  more = TRUE
  it = 1
  while(more) {
    if(it == length(plot_list)) more = FALSE
    print(plot_list[[it]], position = pos_list[it,], more = more)
    it = it + 1   
  }
}

# Same as cv.compare.bubble, only now the plot is made using ggplot
# instead of lattice. I find this version preferable.
cv.compare.ggplot = function(objs, zcol = "residual", layout, col.names, addPoly = NULL) {
  require(ggplot2, quietly = TRUE)
  objs = checkIfautokrige.cv(objs)
  if(!missing(col.names)) names(objs) = col.names

  # Build data for ggplot
  coornames = names(data.frame(coordinates(objs[[1]]$krige.cv_output)))
  dat = melt(lapply(objs, function(x) {
      dum = data.frame(x$krige.cv_output)[c(coornames, zcol)]
      return(dum)
    }), measure.vars = NULL)
  zcol_abs = sprintf("%s_abs", zcol)
  dat[[zcol_abs]] = abs(dat[[zcol]])
  if(!missing(col.names)) dat$L1 = factor(dat$L1, levels = col.names)

  # Make plot
  ggobj = ggplot(aes_string(x = coornames[1], y = coornames[2], size = zcol_abs, color = zcol), data = dat) +
      facet_wrap(~ L1, ncol = layout[1], nrow = layout[2], as.table = TRUE) + 
      scale_x_continuous(name = "", breaks = NA, labels = NA) + 
      scale_y_continuous(name = "", breaks = NA, labels = NA) + coord_equal() +
      scale_color_gradient2("Under- or \noverestimation", high = "green", low = "red", mid = "grey80") + 
      scale_size_continuous("Amount")
  if(!is.null(addPoly)) {
    require(gpclib)
    require(maptools)
    if(inherits(addPoly, "SpatialPolygons")) addPoly = fortify(addPoly) else stop("addPoly should be SpatialPolygons* object")
    addPoly[[zcol]] = min(dat[[zcol]])
    addPoly[[zcol_abs]] = min(dat[[zcol_abs]])
    ggobj = ggobj + geom_path(aes_string(x = "long", y = "lat", group = "group"), color = "lightgrey", size = 0.5, data = addPoly, legend = FALSE)
  }
  return(ggobj + geom_point())
}








