
### BOXPLOT of each variable ####

"boxplot_variables" = function(dataset, variables = NULL, samples = NULL, horizontal = T, 
                               col = "lightblue", nchar.label = 10, cex.axis = 0.8, ...) {
  
  if (is.null(variables)) { # assume all
    variables = rownames(dataset$data)
  } 
  
  if (is.null(samples)) {
    samples = colnames(dataset$data)
  }
  if (is.numeric(variables))
    names.short = substr(get_x_values_as_text(dataset)[variables], 1, nchar.label)
  else 
    names.short = substr(variables, 1, nchar.label)

  if (length(variables) > 1) {
      boxplot(t(dataset$data[variables,samples]), names = names.short, 
          horizontal = horizontal, las = 2, 
          col = col, cex.axis = cex.axis, ...)
  }
  else {
      boxplot(dataset$data[variables,samples], xlab = names.short,  
            horizontal = F, las = 2, 
            col = col, cex.axis = cex.axis, ...)
  }
}

boxplot_vars_factor = function(dataset, meta.var, variables = NULL, samples = NULL, 
                               horizontal = F, nchar.label = 10, col = NULL,
                               vec.par = NULL, cex.axis = 0.8, ylabs = NULL, ...)
{
  if (is.null(variables)) { # assume all
    variables = rownames(dataset$data)
  } 
  
  if (is.null(samples)) {
    samples = colnames(dataset$data)
  }
  
  if (is.numeric(variables))
    names.short = substr(get_x_values_as_text(dataset)[variables], 1, nchar.label)
  else 
    names.short = substr(variables, 1, nchar.label)
  
  if (is.null(vec.par))
    vec.par = c(length(variables), 1)
  
  par(mfrow = vec.par)
  for (i in 1:length(variables)) {
    if (is.null(col)) coli = i+1
    else coli = col
	if (!is.null(ylabs)){
		ylab = ylabs[i]
	} else {
		ylab = NULL
	}
    boxplot(dataset$data[variables[i],samples] ~ dataset$metadata[,meta.var], 
            horizontal = horizontal, las = 2, main = names.short[i], 
            col = coli, cex.axis = cex.axis, ylab = ylab, ...)
  }
  par(mfrow = c(1,1))
}

plotvar_twofactor = function(dataset, variable, meta.var1, meta.var2, colour = "darkblue", title = "", 
                             xlabel = NULL, ylabel = NULL)
{

  df = data.frame(dataset$data[variable,], dataset$metadata[,meta.var1], 
                  dataset$metadata[,meta.var2])
  if (is.numeric(variable)) n1 = rownames(dataset$data)[variable]
  else n1 = variable
  if (is.numeric(meta.var1)) n2 = colnames(dataset$metadata)[meta.var1]
  else n2 = meta.var1
  if (is.numeric(meta.var2)) n3 = colnames(dataset$metadata)[meta.var2]
  else n3 = meta.var2
  colnames(df) = c("name1", "name2", "name3")

  g = ggplot2::ggplot(data = df)
  g = g + ggplot2::aes_string('name2', 'name1')
  g = g + ggplot2::geom_boxplot(fill= colour)
  g = g + ggplot2::facet_grid(. ~ name3)
  if (is.null(xlabel)) g = g + ggplot2::xlab(n2)
  else g = g + ggplot2::xlab(xlabel)
  if (is.null(ylabel)) g = g + ggplot2::ylab(n1)
  else g = g + ggplot2::ylab(ylabel)
  if (title != "") g = g + ggplot2::ggtitle(title)
  g
}

##############################SPECTRA PLOTS###############################

# samples - list of samples to plot
# variable.bounds - interval of x values to plot: [1] - minimum value; [2] - maximum value
# lty, lwd, col - parameters to pass to matplot
# ... - extra parameters passed to matplot function
"plot_spectra_simple" = function(dataset, samples = NULL, variable.bounds = NULL, xlab = NULL,
                               ylab = NULL, lty = 1, lwd = 1, col = 1, reverse.x = F, ...) {
  
  if (is.null(xlab)) xlab = get_x_label(dataset)
  if (is.null(ylab)) ylab = get_value_label(dataset)
  
  if (!is.null(dataset$labels$x) && !is.expression(dataset$labels$x) && dataset$labels$x == "mz/rt"){
	if (is.null(variable.bounds)){
		variables = 1:length(get_x_values_as_text(dataset))
		vars = as.numeric(gsub("/.*", '', get_x_values_as_text(dataset)))

	} else {
		x.vars = as.numeric(gsub("/.*", '', get_x_values_as_text(dataset)))
		variables = which(x.vars > variable.bounds[1] & x.vars < variable.bounds[2])
		vars = x.vars[variables]
	}
  } else {
	  if (is.null(variable.bounds)){
		variables = rownames(dataset$data)
		vars = variables
	  } 
	  else {
		x.vars = get_x_values_as_num(dataset)
		variables = rownames(dataset$data)[x.vars > variable.bounds[1] & x.vars < variable.bounds[2]] 
		vars = variables
	  }
  }
  if (reverse.x) xlim = c( max(as.numeric(vars)), min(as.numeric(vars)) )
  else xlim = range(as.numeric(vars))
  
  if (is.null(samples)){
    samples = colnames(dataset$data)
  } 

  matplot(vars, dataset$data[variables,samples,drop=F], type="l", lty=lty, col = col,
            xlab = xlab, ylab = ylab, xlim = xlim, ...)
}


# plots spectra with colors per groups, given by a metadata variable

# legend.place = "none" if no legend or other accepted string in legend function
# 
"plot_spectra" = function(dataset, column.class, func = NULL, samples = NULL, 
                        variable.bounds = NULL, xlab = NULL, ylab = NULL, lty = 1,
                        legend.place = "topright", cex = 0.8, reverse.x = F, ...) {
  
  if (is.null(xlab)) xlab = get_x_label(dataset)
  if (is.null(ylab)) ylab = get_value_label(dataset)
  
 if (!is.null(dataset$labels$x) && !is.expression(dataset$labels$x) && dataset$labels$x == "mz/rt"){
	if (is.null(variable.bounds)){
		variables = 1:length(get_x_values_as_text(dataset))
		vars = as.numeric(gsub("/.*", '', get_x_values_as_text(dataset)))
	} else {
		x.vars = as.numeric(gsub("/.*", '', get_x_values_as_text(dataset)))
		variables = which(x.vars > variable.bounds[1] & x.vars < variable.bounds[2])
		vars = x.vars[variables]
	}
  } else {
	  if (is.null(variable.bounds)){
		variables = rownames(dataset$data)
		vars = variables
	  } 
	  else {
		x.vars = get_x_values_as_num(dataset)
		variables = rownames(dataset$data)[x.vars > variable.bounds[1] & x.vars < variable.bounds[2]] 
		vars = variables
	  }
  }
  
  if (reverse.x) xlim = c( max(as.numeric(vars)), min(as.numeric(vars)) )
  else xlim = range(as.numeric(vars))
	
  if (is.null(samples)){
		samples = colnames(dataset$data)
		metadata = dataset$metadata[,column.class]
	} 
  else {
		metadata = factor(dataset$metadata[samples, column.class])
	}
	if (is.null(func)){
		matplot(vars, dataset$data[variables,samples], type="l", lty=lty, col=as.integer(metadata), 
            xlab = xlab, ylab = ylab, xlim = xlim, ...)
    if (legend.place != "none")
		  legend(legend.place, levels(metadata), cex=cex, fill = sort(as.integer(factor(levels(metadata))))) 
	} 
  else {
		aggregate.result = aggregate(t(dataset$data[variables,samples]), by = list(metadata), func)
		matplot(vars, t(aggregate.result[-1]), type = "l", lty=1, col=as.integer(aggregate.result[,1]), 
            xlab = xlab, ylab = ylab, xlim = xlim, ...)
		if (legend.place != "none")
		  legend(legend.place, levels(metadata), cex=cex, fill = sort(as.integer(factor(levels(metadata)))))
	}
}


##########################################################################################################
## Multiplot from ggplot2 - function taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

multiplot <- function(plots, plotlist=NULL, file, cols=1, layout=NULL) {


  # Make a list from the ... arguments and plotlist
  plots <- c(plots, plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
