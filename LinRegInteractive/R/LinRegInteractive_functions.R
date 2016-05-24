fxInteractive <- function(model, ...) UseMethod("fxInteractive")

fxInteractive.default <- function(...)
	{
	cat("This class is not supported yet. \nPlease refer to the documentation for which classes this holds.", fill = TRUE)
	}

fxInteractive.lm <- function(model, # object of class lm is mandatory
		
		## Control initial appearance ##
		initial.values           = as.list(NULL), # Initial values for the metric covariates in a named list, default to the means. See details.
		preselect.var            = NA, # Name of continuous variable to be displayed as character or NA for menu selection, default to NA.
		preselect.type           = "effect",
		preselect.groups 		 = NULL, # Index of groups to be preselectes. If set to NULL (default) all groups are preselected.
		
		## Parameters for plot layout ##
		dev.height          = 18,	 # Height of graphic device in cm.
		dev.width           = 18,	 # Width of plot area in graphic device in cm.
		dev.width.legend    = 8,	 # Width of legend area in graphic device in cm.
		dev.pointsize		= 10,	 # Character pointsize of graphic device.
		dev.defined         = FALSE, # Graphic device predefined, e. g. customized for printing? Default to FALSE, see details.
		ylim                = NA,    # With a numeric vector of length 2 the plot limits in y-direction can be set. If NA (the default) these are determined automatically.
		col					= NA,	 # Vector of color specifications for the groups. Passed to the line commands and to the legend. Actual palette and consecutive sequence if NA (default).
		lty					= 1,     # Vector of line type specifications for the groups. Passed to the line commands and to the legend, default to solid lines.
		lwd					= 1,     # Vector of line width specifications for the groups. Passed to the line commands and to the legend, default to 1.
		main				= NA,	 # Label for the plot title.
		main.line			= 1.5,	 # Height in lines for plot title which is passed to title(), default to 1.5.
		xlab				= NA,	 # Label for the x-axis. Name of the selected covariate, if NA.
		ylab				= NA,	 # Label for the y-axis. Name of the selected plot type (see argument label.types), if NA.
		legend.add          = TRUE,  # Should a legend be added to the plot? Default to TRUE, but too many groups can cause problems.
		legend.space        = legend.add, # Should the space for the legend be reserved? Default to the value of legend.add. Setting legend.add to FALSE and legend.space to TRUE plots white space instead of the legend. This can be useful if different termplots will be arranged in a document.  
		legend.only         = FALSE,	  # Plot the legend alone.
		legend.pos			= "center",   # Position of the legend, see legend() for details.
		legend.cex          = 1,	      # Relative size of legend text, reduce for many groups.
		legend.width.factor = 1,	   # Factor by which the width of legend box is manipulated.	
		rug.ticksize		= 0.02,    # Length of rugplot tickmarks. Set to 0, if no rugplot should be drawn.
		rug.col				= "black", # Color of rugplot tickmarks, default to black.	
		vline.actual        = TRUE,    # Add vertical line at actual postion of selceted metric covariate? Default to TRUE.
		pos.hlines       	= c(0,0),  # Positon of horizontal line for [1] effect plot and [2] for marginal effects plot, NA for no lines.
		n.effects           = 100,     # Number of equally space points over the span of the metric covariate used for plotting.
		
		## Parameters for plot snapshot
		autosave.plot       = FALSE,    # Save the initial plot?
		snapshot.plot       = FALSE,    # Save plot as PDF when snapshot button is pressed? Default to FALSE.
		graphics.filename   = "LinRegIntPlot",  # Filename as character for graphic file.
		graphics.numbering  = !autosave.plot,   # Automatically append a 3-digits-number to the filenname to avoid that existing graphic files are overwritten.
		graphics.type       = "pdf",    # Graphics file type.
		
		## Parameters for text-output ##
		factor.sep          = "|",   # Character by which the factors are separated in the groupnames.
		level.sep           = ".",   # Character by which the levels are separated in the groupnames.
		latex2console       = FALSE, # Should the textoutput triggered by the snapshot button be printed as LaTeX?
		xtable.big.mark 	= ".",   # Bigmark character for LaTeX output, argument passed to print.xtable().
		xtable.decimal.mark = ",",   # Decimal character for LaTeX output, argument passed to print.xtable().
		xtable.digits 		= NULL,  # Number of digits, argument passed to xtable().
		xtable.display 		= NULL,  # Display style, argument passed to xtable().
		xtable.booktabs 	= FALSE, # Use LaTeX package booktabs for horizontal lines, argument passed to print.xtable().
				
		## Annotations for panel 
		panel.title         = "Linear Model" ,					# Title used in the panel.
		label.button        = "Snapshot" ,						# Label for the snapshot-button.
		label.slider.act    = "Variable displayed: " ,			# Additional label for the slider of the selected metric covariate. 
		label.box.type      = "Type" ,							# Title for the radiogroup box.
		label.types         = c("effect", "marginal effect"),   # Lables for radiogroup buttons (character vector of length 2).
		label.box.groups    = "Groups" ,						# Title for the checkbox.
		
		## Parameters to control the size of the panel.
		slider.width                = 200 , # Width of each slider.    
		slider.height               = 60  , # Height of each slider.   
		button.height               = 30  , # Height of snapshot button.
		box.type.height             = 70  , # Height of radiobox for type selection.
		box.group.character.width   = 7   , # The width of the boxes is basically a this value times the number of characters.
		box.group.line.height       = 28  , # The height of the checkbox is this value times the number of groups.
		dist.obj.width              = 20  , # Vertical distance between sliders and boxes and vertical margins. 
		dist.obj.height             = 10  , # Horizontal distance between panel objects.
		
		... ) # other graphical parameters passed to par()	
{
	#################################################################	
	# Pick covariates employed, check variables for factors, assign #
	# corresponding objects and initialize graphic window			#
	#################################################################	
	
	# pick covariates which are employed in the model
	if(!is.null(model$data))
		{# extract data from model$data if possible
		X <- get_all_vars(model$terms, model$data)[,-1,drop=FALSE]
		}else
		{# try to extract data from model$model
		X <- get_all_vars(model$terms, model$model)[,-1,drop=FALSE]	
		}
	
	# identify factors from model matrix
	logicalindex.factor <- sapply(X, is.factor)
	factors.present      <- any(logicalindex.factor)
	
	if(factors.present)
	{
		num.level <- sapply(X[,logicalindex.factor,drop=FALSE], function(x) length(levels(x)))
	}
	
	# separate factors and continuous covariates 
	X.continuous   <- X[,!logicalindex.factor, drop=FALSE]
	num.continuous <- dim(X.continuous)[2]
	
	# switch for special treatment of single metric covariate
	single.covariate <- num.continuous==1 
	
	if(factors.present)
	{
		# separate factors 
		X.factor    <- X[,logicalindex.factor, drop=FALSE]
		num.groups  <- prod(num.level)
		
		# build all factor combinations
		factor.comb <- factorCombinations(X.factor, factor.sep=factor.sep, level.sep=level.sep, count=FALSE)
	}
	
	# If not specified in the function call: if there is only one continuous covariate, choose it. 
	# Otherwise select the continuous covariate to be displayed from popup list.
	if(is.na(preselect.var))
	{
		if(single.covariate)
		{
			var.selection <- colnames(X.continuous)[1]
		}else
		{	
			var.selection <- select.list(colnames(X.continuous))
		}
	}else
	{# in no variable is preselected
		if(single.covariate)
		{# if there is only one metric covariate, chose it
			var.selection <- colnames(X.continuous)[1]
		}else
		{# Check if preselected covariate exists, if not, show the menu	
			if(any(colnames(X.continuous)==preselect.var))
			{
				var.selection <- preselect.var		
			}else
			{
				var.selection <- select.list(colnames(X.continuous))	
			}	
		}
	}
	
	# If no variable is selected from popup list leave function.
	if(var.selection=="") return()
	
	# pick continuous covariate under invstigation
	x.actual           <- X.continuous[,var.selection]
	x.actual.sequence  <- seq(min(x.actual, na.rm=TRUE), max(x.actual, na.rm=TRUE), length=n.effects)
	X.pred             <- data.frame(x=x.actual.sequence)
	colnames(X.pred)   <- var.selection
	
	# pick other continuous covariates if present
	if(!single.covariate)
	{	
		X.continuous.other     <- X.continuous[,-match(var.selection, colnames(X.continuous)), drop=FALSE]
		nam.X.continuous.other <- colnames(X.continuous.other)
		num.continuous.other   <- length(nam.X.continuous.other)
		val.X.continuous.other <- colMeans(X.continuous.other)
	}
	
	# If not predefined nor switched off via the argument: set up a new graphic device 
	if(!dev.defined)
	{	
		if(legend.only)
		{# Initialize device for legend only
			dev.new(width=dev.width.legend/2.54, height=dev.height/2.54, pointsize=dev.pointsize, buffered=TRUE, noRStudioGD = TRUE)
			par(cex=1, ...)
		}else
		{# If factors are used as covariates and legend should be printed, split plot region via layout()	
			if(factors.present & legend.space)
			{
				dev.new(width=(dev.width + dev.width.legend)/2.54, height=dev.height/2.54, pointsize=dev.pointsize, buffered=TRUE, noRStudioGD = TRUE)
				fraction.plot <- round((dev.width/(dev.width + dev.width.legend))*100, digits=0)
				layoutmatrix  <- matrix(c(rep(2, fraction.plot), rep(1, (100-fraction.plot))),1,100)
				layout(layoutmatrix)
				par(cex=1, ...)
			}else
			{
				dev.new(width=dev.width/2.54, height=dev.height/2.54, pointsize=dev.pointsize, buffered=TRUE, noRStudioGD = TRUE)
				par(cex=1, ...)
			}
		}
	}
	
	#############################################	
	#  Action when panel items are manipulated 	#
	#  (except action for snapshot-button)		#	
	#############################################
	func.panel.action <- function(panel) 
	{# Main panel function, must be defined within the namespace where panel is established
		# read plot type from panel
		if(length(panel$type)==0) # in the first call of the panel the slots are sometimes NULL   
		{
			type.actual <- "effect"
		}else
		{
			type.actual <- panel$type
		}
		
		# read group selection from panel if groups are present
		if(factors.present)
		{
			if(length(panel$groups)==0) 
			{# in the first call of the panel the slots are sometimes NULL 
				if(is.null(preselect.groups))
				{
					logical.index.groups  <- rep(TRUE, times=num.groups)
				}else
				{
					logical.index.groups  <- rep(FALSE, times=num.groups)
					logical.index.groups[preselect.groups]  <- TRUE
				}
				index.groups.selected <- c(1:num.groups)[logical.index.groups]
				# if all groups are deselected choose the first
				if(sum(logical.index.groups)==0)
				{
					index.groups.selected  <- 1  
				}
			}else
			{
				logical.index.groups  <- panel$groups
				index.groups.selected <- c(1:num.groups)[logical.index.groups]
				
				# if all groups are deselected choose the first
				if(sum(logical.index.groups)==0)
				{
					index.groups.selected  <- 1  
				}
			}
		}# end if factors.present	
		
		# read value of actual variable from panel
		if(length(panel$var.selected.act)==0)
		{
			var.selected.act <- 1
		}else
		{
			var.selected.act <- panel$var.selected.act      
		}   
		
		# build designmatrix for predict-method from continuous covariates
		if(single.covariate)
		{
			X.pred.metr <- X.pred			
		}else
		{# read actual values from panel and add to designmatrix for predict-method
			for(i in 1:num.continuous.other)
			{
				if(length(panel[[nam.X.continuous.other[i]]])==0) # initialization fails
				{
					val.X.continuous.other[i] <- i
				}else
				{
					val.X.continuous.other[i] <- panel[[nam.X.continuous.other[i]]]
				}
			}
			X.pred.metr <- cbind(X.pred, t(val.X.continuous.other))
		}
		
		# If factors are used as covariates, calculate predictions for each group.
		if(factors.present)
		{
			# calculate predictions for each selected group
			list.groups <- as.list(NULL)
			for(i in seq(along=index.groups.selected))
			{
				index.group.actual <- index.groups.selected[i]  
				
				# read actual factor combination
				factor.group.actual <- factor.comb$comb[index.group.actual,,drop=FALSE]
				rownames(factor.group.actual) <- ""
				
				# combine with continuous covariates to new designmatrix
				X.pred <- cbind(X.pred.metr, factor.group.actual)
				
				# calculate effect or marginal effect
				switch(type.actual,
						effect = pred.x.actual <- predict(model, newdata=X.pred),
						marginal = pred.x.actual <- splinefun(x=x.actual.sequence, y=predict(model, newdata=X.pred))(x.actual.sequence, deriv=1))
				
				
				# store results for each and every group 
				sub.list <- list(name=factor.comb$names[index.group.actual],
						prog=pred.x.actual)
				
				list.groups[[i]] <- sub.list 
			}
		}else
		{
			# calculate effect or marginal effect
			switch(type.actual,
					effect = pred.x.actual <- predict(model, newdata=X.pred.metr),
					marginal = pred.x.actual <- splinefun(x=x.actual.sequence, y=predict(model, newdata=X.pred.metr))(x.actual.sequence, deriv=1))
			
			index.groups.selected <- 1
			
			# store results
			list.groups <- as.list(NULL)
			list.groups[[1]] <- list(name="default", prog=pred.x.actual)
		}
		
		# determine plot limits in y-direction
		if(any(is.na(ylim))) ylim <- c(min(unlist(lapply(list.groups,"[",2)), na.rm=TRUE), max(unlist(lapply(list.groups,"[",2)), na.rm=TRUE))		
		
		
		### Plotting commands ###
		func.ploteffects <- function()
		{
			# draw effects for each and every group
			# when no factors are present set number of groups to 1 for correct color handling
			if(!factors.present) num.groups <- 1
			# specify color scheme
			if(all(is.na(col))) 
			{
				col.types <- c(1:num.groups)
			}else
			{	
				col.types <- rep(col, times=ceiling(num.groups/length(col)))
			}
			
			lty.types <- rep(lty, times=ceiling(num.groups/length(lty)))
			lwd.types <- rep(lwd, times=ceiling(num.groups/length(lwd)))
			
			if(factors.present & legend.only)
			{# When legend.only is active plot legend only
				old.mar <- par("mar")
				par(mar=c(old.mar[1],0,old.mar[3],0))
				plot(1,1, type="n", axes=FALSE, xlab=NA, ylab=NA)
				legend.labels <- paste(factor.comb$names[index.groups.selected],"", sep="")
				legend(x=legend.pos, legend=legend.labels, text.width=max(strwidth(legend.labels, cex=legend.cex))*legend.width.factor, col=col.types[index.groups.selected], lty=lty.types[index.groups.selected], lwd=lwd.types[index.groups.selected], cex=legend.cex)
				par(mar=old.mar)    
			}else	
			{# when inactive incorporate other legend parameters
				# when factors are used as covariates and switch is on: add legend first to allow adding elements to the main plot	
				if(factors.present & legend.space)
				{
					old.mar <- par("mar")
					par(mar=c(old.mar[1],0,old.mar[3],0))
					plot(1,1, type="n", axes=FALSE, xlab=NA, ylab=NA)
					if(legend.add)
					{
						legend.labels <- paste(factor.comb$names[index.groups.selected],"", sep="")
						legend(x=legend.pos, legend=legend.labels, text.width=max(strwidth(legend.labels, cex=legend.cex))*legend.width.factor, col=col.types[index.groups.selected], lty=lty.types[index.groups.selected], lwd=lwd.types[index.groups.selected], cex=legend.cex)
					}
					par(mar=old.mar)    
				}		
				
				for(i in seq(along=index.groups.selected))
				{
					# initial plot in first loop	
					if(i==1)
					{
						if(is.na(xlab)) xlab <- var.selection
						if(is.na(ylab))
							switch(type.actual,
									effect   = ylab <- label.types[1],
									marginal = ylab <- label.types[2])
						
						
						plot(x.actual.sequence, rep(1, times=length(x.actual.sequence)), type="n", ylim=ylim, xlab=xlab, ylab=ylab, main=NA)
						title(main=main, line=main.line)
						if(!((rug.ticksize==0)||is.na(rug.ticksize))) rug(x.actual,  ticksize = rug.ticksize, col= rug.col) 
					}
					
					index.group.actual <- index.groups.selected[i]  
					
					# if there is only a litte bit variation in the effect draw a horizontal line
					if(sd(list.groups[[i]]$prog, na.rm = TRUE)<10^-10)
					{
						abline(h=mean(list.groups[[i]]$prog, na.rm = TRUE), col=col.types[index.group.actual], lty=lty.types[index.group.actual], lwd=lwd.types[index.group.actual])	
					}else	
					{ 
						lines(x.actual.sequence, list.groups[[i]]$prog, type="l", col=col.types[index.group.actual], lty=lty.types[index.group.actual], lwd=lwd.types[index.group.actual])
					}     
				}
				
				# vertical line at actual position
				if(vline.actual) abline(v=(var.selected.act))
				
				# draw horizontal line at specified positions
				switch(type.actual,
						effect 	 = abline(h=pos.hlines[1]),
						marginal = abline(h=pos.hlines[2]))	
			}# end if legend.only		
		}# end func.ploteffects
		
		# In conjunction with dev.flush() and buffered=TRUE animation is more fluent.
		dev.hold()
		func.ploteffects()
		dev.flush()
		
		# When autosave is activated directly save the plot
		if(autosave.plot)
		{
			# Create filename for plot
			if(graphics.numbering)
			{
				for(i in 1:1000)
				{
					graphics.filename.candidate <- paste(graphics.filename,"-", formatC(i, width=3, flag="0"), sep="")
					if(file.exists(paste(graphics.filename.candidate, ".", graphics.type, sep=""))) next else break
				}
			}else
			{
				graphics.filename.candidate <- paste(graphics.filename, ".", graphics.type, sep="") 
			}
			
			# Platform dependent save operation
			if(.Platform$OS.type != "windows") 
			{ # Mac OS, Linux
				if (any(graphics.type == c("png","jpeg","jpg","tiff","bmp")))
				{
					sptype <- graphics.type
					if (graphics.type == "jpg") {sptype <- "jpeg"}
					savePlot(filename=graphics.filename.candidate, type=sptype , ... )      
				} 
				if(graphics.type == "pdf")
				{
					dev.copy2pdf(file=graphics.filename.candidate)
				}
				if(graphics.type == "eps")
				{
					dev.copy2eps(file=graphics.filename.candidate)
				}
			}else
			{ # Windows OS
				savePlot(filename = graphics.filename.candidate, type = graphics.type)	
			}
		}
		
		panel
	} # end body action()
	
	#############################################	
	#    Action when snapshot-button is used  	#
	#############################################
	internal.snapshot <- function(panel)
	{# function must be defined within the namespace where panel is established
		# save plot if snapshot.plot is TRUE
		if(snapshot.plot)      
		{               
			# Create filename for plot
			if(graphics.numbering)
			{
				for(i in 1:1000)
				{
					graphics.filename.candidate <- paste(graphics.filename,"-", formatC(i, width=3, flag="0"), sep="")
					if(file.exists(paste(graphics.filename.candidate, ".", graphics.type, sep=""))) next else break
				}
			}else
			{
				graphics.filename.candidate <- paste(graphics.filename, ".", graphics.type, sep="") 
			}
			
			# Platform dependent save operation
			if(.Platform$OS.type != "windows") 
			{ # Mac OS, Linux
				if (any(graphics.type == c("png","jpeg","jpg","tiff","bmp")))
				{
					sptype <- graphics.type
					if (graphics.type == "jpg") {sptype <- "jpeg"}
					savePlot(filename=graphics.filename.candidate, type=sptype , ... )      
				} 
				if(graphics.type == "pdf")
				{
					dev.copy2pdf(file=graphics.filename.candidate)
				}
				if(graphics.type == "eps")
				{
					dev.copy2eps(file=graphics.filename.candidate)
				}
			}else
			{ # Windows OS
				savePlot(filename = graphics.filename.candidate, type = graphics.type)	
			}
		}
		
		# read values from panel and build design matrix for predict-method
		var.selected.act <- panel$var.selected.act
		if(single.covariate)
		{
			X.pred.metr <- data.frame(var.selected.act)
			colnames(X.pred.metr)[1] <- var.selection
			vek.x.continuous.actual <- var.selected.act
			names(vek.x.continuous.actual) <- var.selection
		}else
		{	
			for(i in 1:num.continuous.other)
			{
				val.X.continuous.other[i] <- panel[[nam.X.continuous.other[i]]]
			}
			
			# build designmatrix from continuous covariates
			X.pred.metr <- cbind(var.selected.act, t(val.X.continuous.other))
			colnames(X.pred.metr)[1] <- var.selection
			vek.x.continuous.actual <- as.vector(X.pred.metr)
			names(vek.x.continuous.actual) <- colnames(X.pred.metr)
		}
		
		if(factors.present)
		{
			X.pred <- cbind(X.pred.metr, factor.comb$comb)
		}else
		{
			X.pred <- as.data.frame(X.pred.metr)
		}
				
		output.effect.resp <- cbind(predict(model, newdata=X.pred),1)[,-2, drop=FALSE]		
		colnames(output.effect.resp)  <- c("effect")
		
		if(factors.present)
		{
			rownames(output.effect.resp)  <- factor.comb$names
		}else
		{
			rownames(output.effect.resp)  <- ""   
		}
		
		# calculate ECDF-value for each continuous covariate and  marginal effects in each group
		# Order corresponds to appearance in the original data.frame
		output.vek.x.continuous.actual <- NULL
		F.x.continuous 				   <- NULL
		output.me 		               <- NULL
		
		for(i in 1:num.continuous)
		{
			nam.x.marginal.actual       <- colnames(X.continuous)[i]
			
			x.marginal.actual           <- X.continuous[,nam.x.marginal.actual, drop=TRUE]
			x.marginal.actual.sequence  <- seq(min(x.marginal.actual, na.rm=TRUE), max(x.marginal.actual, na.rm=TRUE), length=n.effects)
			X.marginal.pred             <- data.frame(x=x.marginal.actual.sequence)
			colnames(X.marginal.pred)   <- nam.x.marginal.actual
			
			if(!single.covariate)
			{	
				nam.X.marginal.other        <- colnames(X.continuous)[-i]
				X.marginal.pred             <- cbind(X.marginal.pred, X.pred.metr[,nam.X.marginal.other, drop=FALSE])
				colnames(X.marginal.pred)   <- c(nam.x.marginal.actual, nam.X.marginal.other)   
			}
			
			x.marginal.actual.selected     <- X.pred.metr[,nam.x.marginal.actual, drop=TRUE]

			output.vek.x.continuous.actual <- c(output.vek.x.continuous.actual, x.marginal.actual.selected)
      
      		# ECDF-value	
			F.x.continuous.actual        <- sum(x.marginal.actual <= x.marginal.actual.selected)/length(x.marginal.actual)
			F.x.continuous               <- c(F.x.continuous, F.x.continuous.actual)
			
			# calculate marginals for actual continuous covariate in each group
			marginal.actual.groups <- NULL
			if(factors.present)
			{
				for(i in 1:num.groups)
				{
					# read actual factor combination
					factor.group.actual <- factor.comb$comb[i,,drop=FALSE]
					rownames(factor.group.actual) <- ""
					
					# combine with continuous covariates to new designmatrix
					X.marginal.pred.actual <- cbind(X.marginal.pred, factor.group.actual)
					
					marginal.actual.groups.actual   <- splinefun(x=x.marginal.actual.sequence, y=predict(model, newdata=X.marginal.pred.actual))(x.marginal.actual.selected, deriv=1)
					marginal.actual.groups          <- c(marginal.actual.groups, marginal.actual.groups.actual)
				}   
			}else
			{
				marginal.actual.groups  <-	splinefun(x=x.marginal.actual.sequence, y=predict(model, newdata=X.marginal.pred))(x.marginal.actual.selected, deriv=1) 
			}
			
			output.me <- cbind(output.me, marginal.actual.groups)
		}#end for-loop continuous covariates
		
		output.x.continuous                   <- rbind(output.vek.x.continuous.actual, F.x.continuous)
		colnames(output.x.continuous)         <- colnames(X.continuous)
		rownames(output.x.continuous)[c(1:2)] <- c("value", "ECDF(value)")
		
		if(factors.present)
		{
			rownames(output.me) <- factor.comb$names
		}else
		{
			rownames(output.me) <- "marginal effect"   
		}
		
		# Add actual values to table of marginal effects
		colnames(output.me) <- colnames(X.continuous)
		output.me           <- rbind(output.x.continuous, output.me)
		
		if(!latex2console)
		{
			print(summary(model))
			cat(fill=TRUE)
			cat("Selected values of metric covariates", fill=TRUE)
			print(output.x.continuous, fill=TRUE)
			
			cat(fill=TRUE)
			cat("Effects in different groups for selected values of metric covariates", fill=TRUE)
			print(output.effect.resp)
			
			cat(fill=TRUE)
			cat("Marginal effects in different groups for selected values of metric covariates", fill=TRUE)
			print(output.me)
			
			}else
			{
			vec.align <- c("l", rep("r", times=dim(output.x.continuous)[2]))
			output.x.continuous.xtable <- xtable(output.x.continuous, digits = xtable.digits, display = xtable.display,
					align = vec.align,
					caption = "Selected values of metric covariates",
					label = "tab-values")
			print(output.x.continuous.xtable, 
					format.args=list(big.mark = xtable.big.mark, decimal.mark = xtable.decimal.mark),
					booktabs = xtable.booktabs)
			
			vec.align <- c("l", rep("r", times=dim(output.effect.resp)[2]))
			output.effect.resp.xtable <- xtable(output.effect.resp, digits = xtable.digits, display = xtable.display,
					align = vec.align,
					caption = "Effects in different groups for selected values of metric covariates",
					label = "tab-effects")
			print(output.effect.resp.xtable, 
					format.args=list(big.mark = xtable.big.mark, decimal.mark = xtable.decimal.mark),
					booktabs = xtable.booktabs,
					include.rownames = factors.present)
			
			vec.align <- c("l", rep("r", times=dim(output.me)[2]))
			output.me.xtable <- xtable(output.me, digits = xtable.digits, display = xtable.display,
					align = vec.align,
					caption = "Marginal effects in different groups for selected values of metric covariates",
					label = "tab-marginaleffects")
			print(output.me.xtable, hline.after=c(-1,0,2,nrow(output.me.xtable)), 
					format.args=list(big.mark = xtable.big.mark, decimal.mark = xtable.decimal.mark),
					booktabs = xtable.booktabs)
			
		}# end if latex2console
		panel
	}# end body internal.snapshot
	
	#############################################	
	#  Determine panel dimensions and layout 	#
	#############################################	
	
	# Calculate slider positions and sizes
	y.pos.slider.other <- 2.5*dist.obj.height + slider.height
	
	if(factors.present)
	{
		# Determine maximum length of group name
		length.factorgroup.name <- NULL 
		for(i in 1:num.groups)
		{
			length.factorgroup.name <- c(length.factorgroup.name, length(strsplit(factor.comb$names[i], split=NULL)[[1]]))
		}
		max.length.factorgroup.name <- max(length.factorgroup.name)
	}else
	{#if no groups are present set to fix value
		max.length.factorgroup.name <- 20   
	}
	
	# Calculate button position and size
	button.width        <- max(box.group.character.width*max.length.factorgroup.name, box.group.character.width*18)
	button.pos.x        <- 2*dist.obj.width + slider.width
	button.pos.y        <- dist.obj.height
	
	# Calculate box positions and sizes
	boxes.width         <- button.width
	box.pos.x           <- button.pos.x
	box.type.pos.y      <- button.pos.y + button.height + dist.obj.height
	box.groups.pos.y    <- box.type.pos.y  + box.type.height + dist.obj.height
	
	if(factors.present)
	{
		box.group.height <- box.group.line.height*length(factor.comb$names)
	}else
	{
		box.group.height <- 0   
	}
	
	# Calculate overall panel size
	overall.width   <- box.pos.x + dist.obj.width + boxes.width
	overall.height  <- max((4*dist.obj.height + box.type.height + box.group.height + button.height), ((num.continuous-1)*slider.height + y.pos.slider.other + dist.obj.height))
	
	#########################################################	
	#  Define main panel and populate it with GUI controls 	#
	#########################################################
	mainpanel <- rp.control(title=panel.title, size = c(overall.width,overall.height))
	
	# Slider for displayed continuous variable
	initial.act <- initial.values[[var.selection]]
	if(is.null(initial.act)) initial.act <- mean(x.actual, na.rm = TRUE) 
	
	eval(wrapRpSlider(panel="mainpanel", variable="var.selected.act", title=paste(label.slider.act, var.selection, sep=""), 
					from = min(x.actual, na.rm=TRUE), to=max(x.actual, na.rm=TRUE), initval=initial.act,
					pos= c(dist.obj.width, dist.obj.height, slider.width, slider.height) , showvalue =TRUE, action = "func.panel.action"))
	
	# Sliders for other continuous variables (if present)
	if(!single.covariate)
	{
		for(i in 1:num.continuous.other)
		{
			y.pos.slider <- y.pos.slider.other + (i-1)*slider.height
			
			initial.act <- initial.values[[nam.X.continuous.other[i]]]
			if(is.null(initial.act)) initial.act <- mean(X.continuous.other[,i], na.rm = TRUE)
			
			slidercall <- wrapRpSlider(panel="mainpanel", variable=nam.X.continuous.other[i],title=nam.X.continuous.other[i], 
					from = min(X.continuous.other[,i], na.rm=TRUE), to=max(X.continuous.other[,i], na.rm=TRUE), initval=initial.act,
					pos= c(dist.obj.width, y.pos.slider, slider.width, slider.height) , showvalue =TRUE, action = "func.panel.action")
			eval(slidercall)
		}
	}
	
	# Button for snapshot
	rp.button(panel=mainpanel, action = internal.snapshot, title = label.button , pos = c(button.pos.x, button.pos.y, button.width, button.height))
	
	# Radiogroup for type
	# if type is preselected check for allowed value
	if(any(preselect.type==c("effect", "marginal")))
	{# if allowed use it
		type.initval <- preselect.type[1]
	}else
	{# if not allowed set as effect
		type.initval <- "effect"	
	}
	
	# Cheat R CMD check	
	radiogroup.call <- "rp.radiogroup(panel=mainpanel, variable=type, vals=c('effect', 'marginal'), initval=type.initval, labels=label.types, 
			title = label.box.type, action = func.panel.action, pos = c(box.pos.x, box.type.pos.y, boxes.width, box.type.height))"
	eval(parse(text=radiogroup.call)) 
	
	# Checkbox for groups, if factors are employed in the model
	if(factors.present)
	{
		if(is.null(preselect.groups))
		{
			init.logical.index.groups  <- rep(TRUE, times=num.groups)
		}else
		{
			init.logical.index.groups  <- rep(FALSE, times=num.groups)
			init.logical.index.groups[preselect.groups]  <- TRUE
		}
		
		# Cheat R CMD check	
		checkbox.call <- "rp.checkbox(panel=mainpanel, variable=groups, action = func.panel.action, initval=init.logical.index.groups,
				labels = factor.comb$names, title = label.box.groups, pos = c(box.pos.x, box.groups.pos.y, boxes.width, box.group.height))"
		eval(parse(text=checkbox.call)) 
	}
	
	# Initalize panel
	rp.do(panel=mainpanel, func.panel.action)
	
	# when autosave is activated close panel after initialization
	if(autosave.plot) rp.control.dispose(panel=mainpanel)
}


fxInteractive.glm <- function(model, # object of class glm is mandatory
		
		## Control initial appearance ##
		initial.values           = as.list(NULL), # Initial values for the metric covariates in a named list, default to the means. See details.
		preselect.var            = NA, # Name of continuous variable to be displayed as character or NA for menu selection, default to NA.
		preselect.type           = "link",
		preselect.groups 		 = NULL, # Index of groups to be preselectes. If set to NULL (default) all groups are preselected.
		
		## Parameters for plot layout ##
		dev.height          = 18,	 # Height of graphic device in cm.
		dev.width           = 18,	 # Width of plot area in graphic device in cm.
		dev.width.legend    = 8,	 # Width of legend area in graphic device in cm.
		dev.pointsize		= 10,	 # Character pointsize of graphic device.
		dev.defined         = FALSE, # Graphic device predefined, e. g. customized for printing? Default to FALSE, see details.
		ylim                = NA,    # With a numeric vector of length 2 the plot limits in y-direction can be set. If NA (the default) these are determined automatically.
		col					= NA,	 # Vector of color specifications for the groups. Passed to the line commands and to the legend. Actual palette and consecutive sequence if NA (default).
		lty					= 1,     # Vector of line type specifications for the groups. Passed to the line commands and to the legend, default to solid lines.
		lwd					= 1,     # Vector of line width specifications for the groups. Passed to the line commands and to the legend, default to 1.
		main				= NA,	 # Label for the plot title.
		main.line			= 1.5,	 # Height in lines for plot title which is passed to title(), default to 1.5.
		xlab				= NA,	 # Label for the x-axis. Name of the selected covariate, if NA.
		ylab				= NA,	 # Label for the y-axis. Name of the selected plot type (see argument label.types), if NA.
		legend.add          = TRUE,  # Should a legend be added to the plot? Default to TRUE, but too many groups can cause problems.
		legend.space        = legend.add, # Should the space for the legend be reserved? Default to the value of legend.add. Setting legend.add to FALSE and legend.space to TRUE plots white space instead of the legend. This can be useful if different termplots will be arranged in a document.  
		legend.only         = FALSE,	  # Plot the legend alone.
		legend.pos			= "center",   # Position of the legend, see legend for details.
		legend.cex          = 1,	      # Relative size of legend text, reduce for many groups.
		legend.width.factor = 1, 		  # Factor by which the width of legend box is manipulated.		
		rug.ticksize		= 0.02,       # Length of rugplot tickmarks. Set to 0, if no rugplot should be drawn.
		rug.col				= "black",	  # Color of rugplot tickmarks, default to black.
		vline.actual        = TRUE,       # Add vertical line at actual postion of selceted metric covariate? Default to TRUE.
		pos.hlines       	= c(0,0.5,0), # Positon of horizontal line for [1] link plot, [2] response plot and [3] for marginal effects plot, NA for no lines.
		n.effects           = 100,        # Number of equally spaced points over the span of the selected metric covariate used for plotting.
		
		## Parameters for plot snapshot
		autosave.plot       = FALSE,    # Save the initial plot?
		snapshot.plot       = FALSE,    # Save plot when snapshot button is pressed? Default to FALSE.
		graphics.filename   = "LinRegIntPlot",  # Filename as character for graphic file.
		graphics.numbering = !autosave.plot,    # Automatically append a 3-digits-number to the filenname to avoid that existing graphic files are overwritten.
		graphics.type       = "pdf",	# Graphics file type.
		
		## Parameters for text-output ##
		factor.sep          = "|",   # Character by which the factors are separated in the groupnames.
		level.sep           = ".",   # Character by which the levels are separated in the groupnames.
		latex2console       = FALSE, # Should the textoutput triggered by the snapshot button be printed as LaTeX?
		xtable.big.mark 	= ".", # Bigmark character for LaTeX output, argument passed to print.xtable().
		xtable.decimal.mark = ",", # Decimal character for LaTeX output, argument passed to print.xtable().
		xtable.digits 		= NULL, # Number of digits, argument passed to xtable().
		xtable.display 		= NULL, # Display style, argument passed to xtable().
		xtable.booktabs 	= FALSE, # Use LaTeX package booktabs for horizontal lines, argument passed to print.xtable().
		
		## Annotations for panel 
		panel.title         = "Generalized Linear Model" ,                           # Title used in the panel.
		label.button        = "Snapshot" ,                                            # Label for the snapshot-button.
		label.slider.act    = "Variable displayed: " ,                                # Additional label for the slider of the selected metric covariate. 
		label.box.type      = "Type" ,                                                # Title for the radiogroup box.
		label.types         = c("linear predictor", "response", "marginal effect"),   # Lables for radiogroup buttons (character vector of length 3).
		label.box.groups    = "Groups" ,                                              # Title for the checkbox.
		
		## Parameters to control the size of the panel.
		slider.width                = 200 , # Width of each slider.    
		slider.height               = 60  , # Height of each slider.   
		button.height               = 30  , # Height of snapshot button.
		box.type.height             = 90  , # Height of radiobox for type selection.
		box.group.character.width   = 7   , # The width of the boxes is basically a this value times the number of characters.
		box.group.line.height       = 28  , # The height of the checkbox is this value times the number of groups.
		dist.obj.width              = 20  , # Vertical distance between sliders and boxes and vertical margins. 
		dist.obj.height             = 10  , # Horizontal distance between panel objects.
		
		... ) # other graphical parameters passed to par()	
{
	#################################################################	
	# Pick covariates employed, check variables for factors, assign #
	# corresponding objects and initialize graphic window			#
	#################################################################
		
	# pick covariates which are employed in the model
	if(!is.null(model$data))
		{# extract data from model$data
		X <- get_all_vars(model$terms, model$data)[,-1,drop=FALSE]
		}else
		{# try to extract data from model$model
		X <- get_all_vars(model$terms, model$model)[,-1,drop=FALSE]	
		}
	
	# identify factors from model matrix
	logicalindex.factor <- sapply(X, is.factor)
	factors.present      <- any(logicalindex.factor)
	
	if(factors.present)
	{
		num.level <- sapply(X[,logicalindex.factor,drop=FALSE], function(x) length(levels(x)))
	}
	
	# separate factors and continuous covariates 
	X.continuous   <- X[,!logicalindex.factor, drop=FALSE]
	num.continuous <- dim(X.continuous)[2]
	
	# switch for special treatment of single metric covariates
	single.covariate <- num.continuous==1 
	
	if(factors.present)
	{
		# separate factors 
		X.factor    <- X[,logicalindex.factor, drop=FALSE]
		num.groups  <- prod(num.level)
		
		# build all factor combinations
		factor.comb <- factorCombinations(X.factor, factor.sep=factor.sep, level.sep=level.sep, count=FALSE)
	}
	
	# If not specified in the function call: if there is only one continuous covariate, choose it. 
	# Otherwise select the continuous covariate to be displayed from popup list.
	if(is.na(preselect.var))
	{
		if(single.covariate)
		{
			var.selection <- colnames(X.continuous)[1]
		}else
		{	
			var.selection <- select.list(colnames(X.continuous))
		}
	}else
	{# in no variable is preselected
		if(single.covariate)
		{# if there is only one metric covariate, chose it
			var.selection <- colnames(X.continuous)[1]
		}else
		{# Check if preselected covariate exists, if not, show the menu	
			if(any(colnames(X.continuous)==preselect.var))
			{
				var.selection <- preselect.var		
			}else
			{
				var.selection <- select.list(colnames(X.continuous))	
			}	
		}
	}
	
	# If no variable is selected from popup list leave function.
	if(var.selection=="") return()
	
	# pick continuous covariate under invstigation
	x.actual           <- X.continuous[,var.selection]
	x.actual.sequence  <- seq(min(x.actual, na.rm=TRUE), max(x.actual, na.rm=TRUE), length=n.effects)
	X.pred             <- data.frame(x=x.actual.sequence)
	colnames(X.pred)   <- var.selection
	
	# pick other continuous covariates if present
	if(!single.covariate)
	{	
		X.continuous.other     <- X.continuous[,-match(var.selection, colnames(X.continuous)), drop=FALSE]
		nam.X.continuous.other <- colnames(X.continuous.other)
		num.continuous.other   <- length(nam.X.continuous.other)
		val.X.continuous.other <- colMeans(X.continuous.other)
	}
	
	# If not predefined nor switched off via the argument: set up a new graphic device 
	if(!dev.defined)
	{	
		if(legend.only)
		{# Initialize device for legend only
			dev.new(width=dev.width.legend/2.54, height=dev.height/2.54, pointsize=dev.pointsize, buffered=TRUE, noRStudioGD = TRUE)
			par(cex=1, ...)
		}else
		{# If factors are used as covariates and legend should be printed, split plot region via layout()	
			if(factors.present & legend.space)
			{
				dev.new(width=(dev.width + dev.width.legend)/2.54, height=dev.height/2.54, pointsize=dev.pointsize, buffered=TRUE, noRStudioGD = TRUE)
				fraction.plot <- round((dev.width/(dev.width + dev.width.legend))*100, digits=0)
				layoutmatrix  <- matrix(c(rep(2, fraction.plot), rep(1, (100-fraction.plot))),1,100)
				layout(layoutmatrix)
				par(cex=1, ...)
			}else
			{
				dev.new(width=dev.width/2.54, height=dev.height/2.54, pointsize=dev.pointsize, buffered=TRUE, noRStudioGD = TRUE)
				par(cex=1, ...)
			}
		}
	}
	
	#############################################	
	#  Action when panel items are manipulated 	#
	#  (except action for snapshot-button)		#	
	#############################################
	func.panel.action <- function(panel) 
	{# Main panel function, must be defined within the namespace where panel is established
		# read plot type from panel
		if(length(panel$type)==0) # in the first call of the panel the slots are sometimes NULL   
		{
			type.actual <- "link"
		}else
		{
			type.actual <- panel$type
		}
		
		# read group selection from panel if groups are present
		if(factors.present)
		{
			if(length(panel$groups)==0) 
			{# in the first call of the panel the slots are sometimes NULL 
				if(is.null(preselect.groups))
				{
					logical.index.groups  <- rep(TRUE, times=num.groups)
				}else
				{
					logical.index.groups  <- rep(FALSE, times=num.groups)
					logical.index.groups[preselect.groups]  <- TRUE
				}
				index.groups.selected <- c(1:num.groups)[logical.index.groups]
				# if all groups are deselected choose the first
				if(sum(logical.index.groups)==0)
				{
					index.groups.selected  <- 1  
				}
			}else
			{
				logical.index.groups  <- panel$groups
				index.groups.selected <- c(1:num.groups)[logical.index.groups]
				
				# if all groups are deselected choose the first
				if(sum(logical.index.groups)==0)
				{
					index.groups.selected  <- 1  
				}
			}
		}# end if factors.present
		
		# read value of actual variable from panel
		if(length(panel$var.selected.act)==0)
		{
			var.selected.act <- 1
		}else
		{
			var.selected.act <- panel$var.selected.act      
		}   
		
		# build designmatrix for predict-method from continuous covariates
		if(single.covariate)
		{
			X.pred.metr <- X.pred			
		}else
		{# read actual values from panel and add to designmatrix for predict-method
			for(i in 1:num.continuous.other)
			{
				if(length(panel[[nam.X.continuous.other[i]]])==0) # initialization fails
				{
					val.X.continuous.other[i] <- i
				}else
				{
					val.X.continuous.other[i] <- panel[[nam.X.continuous.other[i]]]
				}
			}
			X.pred.metr <- cbind(X.pred, t(val.X.continuous.other))
		}
		
		# If factors are used as covariates, calculate predictions for each group.
		if(factors.present)
		{
			# calculate predictions for each selected group
			list.groups <- as.list(NULL)
			for(i in seq(along=index.groups.selected))
			{
				index.group.actual <- index.groups.selected[i]  
				
				# read actual factor combination
				factor.group.actual <- factor.comb$comb[index.group.actual,,drop=FALSE]
				rownames(factor.group.actual) <- ""
				
				# combine with continuous covariates to new designmatrix
				X.pred <- cbind(X.pred.metr, factor.group.actual)
				
				# calculate linear predictor, probability or marginal effect
				switch(type.actual,
						link = pred.x.actual <- predict(model, newdata=X.pred, type="link"),
						response = pred.x.actual <- predict(model, newdata=X.pred, type="response"),
						marginal = pred.x.actual <- splinefun(x=x.actual.sequence, y=predict(model, newdata=X.pred, type="response"))(x.actual.sequence, deriv=1))
				
				# store results for each and every group 
				sub.list <- list(name=factor.comb$names[index.group.actual],
						prog=pred.x.actual)
				
				list.groups[[i]] <- sub.list 
			}
		}else
		{
			# calculate probability, linear predictor or marginal effect
			switch(type.actual,
					link     = pred.x.actual <- predict(model, newdata=X.pred.metr, type="link"),
					response = pred.x.actual <- predict(model, newdata=X.pred.metr, type="response"),
					marginal = pred.x.actual <- splinefun(x=x.actual.sequence, y=predict(model, newdata=X.pred.metr, type="response"))(x.actual.sequence, deriv=1))
			
			index.groups.selected <- 1
			
			# store results
			list.groups <- as.list(NULL)
			list.groups[[1]] <- list(name="default", prog=pred.x.actual)
		}
		
		# determine plot limits in y-direction
		if(any(is.na(ylim))) ylim <- c(min(unlist(lapply(list.groups,"[",2)), na.rm=TRUE), max(unlist(lapply(list.groups,"[",2)), na.rm=TRUE))		
		
		### Plotting commands ###
		func.ploteffects <- function()
		{
			# draw effects for each and every group
			# when no factors are present set number of groups to 1 for correct color handling
			if(!factors.present) num.groups <- 1
			# specify color scheme
			if(all(is.na(col))) 
			{
				col.types <- c(1:num.groups)
			}else
			{	
				col.types <- rep(col, times=ceiling(num.groups/length(col)))
			}
			
			lty.types <- rep(lty, times=ceiling(num.groups/length(lty)))
			lwd.types <- rep(lwd, times=ceiling(num.groups/length(lwd)))
			
			if(factors.present & legend.only)
			{# When legend.only is active plot legend only
				old.mar <- par("mar")
				par(mar=c(old.mar[1],0,old.mar[3],0))
				plot(1,1, type="n", axes=FALSE, xlab=NA, ylab=NA)
				legend.labels <- paste(factor.comb$names[index.groups.selected],"", sep="")
				legend(x=legend.pos, legend=legend.labels, text.width=max(strwidth(legend.labels, cex=legend.cex))*legend.width.factor, col=col.types[index.groups.selected], lty=lty.types[index.groups.selected], lwd=lwd.types[index.groups.selected], cex=legend.cex)
				par(mar=old.mar)    
			}else	
			{# when inactive incorporate other legend parameters
				# when factors are used as covariates and switch is on: add legend first to allow adding elements to the main plot	
				if(factors.present & legend.space)
				{
					old.mar <- par("mar")
					par(mar=c(old.mar[1],0,old.mar[3],0))
					plot(1,1, type="n", axes=FALSE, xlab=NA, ylab=NA)
					if(legend.add)
					{
						legend.labels <- paste(factor.comb$names[index.groups.selected],"", sep="")
						legend(x=legend.pos, legend=legend.labels, text.width=max(strwidth(legend.labels, cex=legend.cex))*legend.width.factor, col=col.types[index.groups.selected], lty=lty.types[index.groups.selected], lwd=lwd.types[index.groups.selected], cex=legend.cex)
					}
					par(mar=old.mar)    
				}
				
				for(i in seq(along=index.groups.selected))
				{
					# initial plot in first loop	
					if(i==1)
					{
						if(is.na(xlab)) xlab <- var.selection
						if(is.na(ylab)) switch(type.actual,
									link     = ylab <- label.types[1],
									response = ylab <- label.types[2],
									marginal = ylab <- label.types[3])
						
						plot(x.actual.sequence, rep(1, times=length(x.actual.sequence)), type="n", ylim=ylim, xlab=xlab, ylab=ylab, main=NA)
						title(main=main, line=main.line)
						if(!((rug.ticksize==0)||is.na(rug.ticksize))) rug(x.actual,  ticksize = rug.ticksize, col= rug.col) 
					}
					
					index.group.actual <- index.groups.selected[i]  
					
					# if there is only a litte bit variation in the effect draw a horizontal line
					if(sd(list.groups[[i]]$prog, na.rm = TRUE)<10^-10)
					{
						abline(h=mean(list.groups[[i]]$prog, na.rm = TRUE), col=col.types[index.group.actual], lty=lty.types[index.group.actual], lwd=lwd.types[index.group.actual])	
					}else	
					{ 
						lines(x.actual.sequence, list.groups[[i]]$prog, type="l", col=col.types[index.group.actual], lty=lty.types[index.group.actual], lwd=lwd.types[index.group.actual])
					}    
				}	
				
				# vertical line at actual position
				if(vline.actual) abline(v=(var.selected.act))
				
				# draw horizontal line at specified positions
				switch(type.actual,
						link 	 = abline(h=pos.hlines[1]),
						response = abline(h=pos.hlines[2]),
						marginal = abline(h=pos.hlines[3]))
				
			}# end if legend.only
			
		}# end func.ploteffects
		
		# in conjunction with dev.flush() and buffered=TRUE animation is more fluent
		dev.hold()
		func.ploteffects()
		dev.flush()
		
		# When autosave is activated directly save the plot
		if(autosave.plot)
		{
			# Create filename for plot
			if(graphics.numbering)
			{
				for(i in 1:1000)
				{
					graphics.filename.candidate <- paste(graphics.filename,"-", formatC(i, width=3, flag="0"), sep="")
					if(file.exists(paste(graphics.filename.candidate, ".", graphics.type, sep=""))) next else break
				}
			}else
			{
				graphics.filename.candidate <- paste(graphics.filename, ".", graphics.type, sep="") 
			}
			
			# Platform dependent save operation
			if(.Platform$OS.type != "windows") 
			{ # Mac OS, Linux
				if (any(graphics.type == c("png","jpeg","jpg","tiff","bmp")))
				{
					sptype <- graphics.type
					if (graphics.type == "jpg") {sptype <- "jpeg"}
					savePlot(filename=graphics.filename.candidate, type=sptype , ... )      
				} 
				if(graphics.type == "pdf")
				{
					dev.copy2pdf(file=graphics.filename.candidate)
				}
				if(graphics.type == "eps")
				{
					dev.copy2eps(file=graphics.filename.candidate)
				}
			}else
			{ # Windows OS
				savePlot(filename = graphics.filename.candidate, type = graphics.type)	
			}
			
		}
		
		panel
	} # end body action()
	
	#############################################	
	#    Action when snapshot-button is used  	#
	#############################################
	internal.snapshot <- function(panel)
	{# function must be defined within the namespace where panel is established
		# save plot if snapshot.plot is TRUE
		if(snapshot.plot)      
		{
			# Create filename for plot
			if(graphics.numbering)
			{
				for(i in 1:1000)
				{
					graphics.filename.candidate <- paste(graphics.filename,"-", formatC(i, width=3, flag="0"), sep="")
					if(file.exists(paste(graphics.filename.candidate, ".", graphics.type, sep=""))) next else break
				}
			}else
			{
				graphics.filename.candidate <- paste(graphics.filename, ".", graphics.type, sep="") 
			}
			
			# Platform dependent save operation
			if(.Platform$OS.type != "windows") 
			{ # Mac OS, Linux
				if (any(graphics.type == c("png","jpeg","jpg","tiff","bmp")))
				{
					sptype <- graphics.type
					if (graphics.type == "jpg") {sptype <- "jpeg"}
					savePlot(filename=graphics.filename.candidate, type=sptype , ... )      
				} 
				if(graphics.type == "pdf")
				{
					dev.copy2pdf(file=graphics.filename.candidate)
				}
				if(graphics.type == "eps")
				{
					dev.copy2eps(file=graphics.filename.candidate)
				}
			}else
			{ # Windows OS
				savePlot(filename = graphics.filename.candidate, type = graphics.type)	
			}
		}       
		
		# read values from panel and build design matrix for predict-method
		var.selected.act <- panel$var.selected.act
		if(single.covariate)
		{
			X.pred.metr <- data.frame(var.selected.act)
			colnames(X.pred.metr)[1] <- var.selection
			vek.x.continuous.actual <- var.selected.act
			names(vek.x.continuous.actual) <- var.selection
		}else
		{	
			for(i in 1:num.continuous.other)
			{
				val.X.continuous.other[i] <- panel[[nam.X.continuous.other[i]]]
			}
			
			# build designmatrix from continuous covariates
			X.pred.metr <- cbind(var.selected.act, t(val.X.continuous.other))
			colnames(X.pred.metr)[1] <- var.selection
			vek.x.continuous.actual <- as.vector(X.pred.metr)
			names(vek.x.continuous.actual) <- colnames(X.pred.metr)
		}
		
		if(factors.present)
		{
			X.pred <- cbind(X.pred.metr, factor.comb$comb)
		}else
		{
			X.pred <- as.data.frame(X.pred.metr)
		}
		
		pred.link.x.actual      <- predict(model, newdata=X.pred, type="link")
		pred.response.x.actual  <- predict(model, newdata=X.pred, type="response")
		
		output.effect.resp            <- cbind(pred.link.x.actual, pred.response.x.actual)
		colnames(output.effect.resp)  <- c("link", "response")
		
		if(factors.present)
		{
			rownames(output.effect.resp)  <- factor.comb$names
		}else
		{
			rownames(output.effect.resp)  <- ""  
		}
		
		# calculate ECDF-value for each continuous covariate and  marginal effects in each group
		output.vek.x.continuous.actual <- NULL
		F.x.continuous 				   <- NULL
		output.me 					   <- NULL
		
		for(i in 1:num.continuous)
		{
			nam.x.marginal.actual       <- colnames(X.continuous)[i]
			
			x.marginal.actual           <- X.continuous[,nam.x.marginal.actual, drop=TRUE]
			x.marginal.actual.sequence  <- seq(min(x.marginal.actual, na.rm=TRUE), max(x.marginal.actual, na.rm=TRUE), length=n.effects)
			X.marginal.pred             <- data.frame(x=x.marginal.actual.sequence)
			colnames(X.marginal.pred)   <- nam.x.marginal.actual
			
			if(!single.covariate)
			{	
				nam.X.marginal.other        <- colnames(X.continuous)[-i]
				X.marginal.pred             <- cbind(X.marginal.pred, X.pred.metr[,nam.X.marginal.other, drop=FALSE])
				colnames(X.marginal.pred)   <- c(nam.x.marginal.actual, nam.X.marginal.other)   
			}
			
			x.marginal.actual.selected  <- X.pred.metr[,nam.x.marginal.actual, drop=TRUE]
			
			output.vek.x.continuous.actual <- c(output.vek.x.continuous.actual, x.marginal.actual.selected)
			
			# ECDF-value	
			F.x.continuous.actual        <- sum(x.marginal.actual <= x.marginal.actual.selected)/length(x.marginal.actual)
			F.x.continuous               <- c(F.x.continuous, F.x.continuous.actual)
			
			# calculate marginals for actual continuous covariate in each group
			marginal.actual.groups <- NULL
			if(factors.present)
			{
				for(i in 1:num.groups)
				{
					# read actual factor combination
					factor.group.actual <- factor.comb$comb[i,,drop=FALSE]
					rownames(factor.group.actual) <- ""
					
					# combine with continuous covariates to new designmatrix
					X.marginal.pred.actual <- cbind(X.marginal.pred, factor.group.actual)
					
					marginal.actual.groups.actual   <- splinefun(x=x.marginal.actual.sequence, y=predict(model, newdata=X.marginal.pred.actual, type="response"))(x.marginal.actual.selected, deriv=1)
					marginal.actual.groups          <- c(marginal.actual.groups, marginal.actual.groups.actual)
				}   
			}else
			{								
				marginal.actual.groups  <- splinefun(x=x.marginal.actual.sequence, y=predict(model, newdata=X.marginal.pred, type="response"))(x.marginal.actual.selected, deriv=1)
			}
			
			output.me <- cbind(output.me, marginal.actual.groups)
		}#end for-loop continuous covariates
		
		output.x.continuous                   <- rbind(output.vek.x.continuous.actual, F.x.continuous)
		colnames(output.x.continuous)         <- colnames(X.continuous)
		rownames(output.x.continuous)[c(1:2)] <- c("value", "ECDF(value)")
		
		if(factors.present)
		{
			rownames(output.me) <- factor.comb$names
		}else
		{
			rownames(output.me) <- "marginal effect"   
		}
		
		# Add actual values to table of marginal effects
		colnames(output.me) <- colnames(X.continuous)
		output.me           <- rbind(output.x.continuous, output.me)
		
		if(!latex2console)
		{
			print(summary(model))
			cat(fill=TRUE)
			cat("Selected values of metric covariates", fill=TRUE)
			print(output.x.continuous, fill=TRUE)
			
			cat(fill=TRUE)
			cat("Effects in different groups for selected values of metric covariates", fill=TRUE)
			print(output.effect.resp)
			
			cat(fill=TRUE)
			cat("Marginal effects in different groups for selected values of metric covariates", fill=TRUE)
			print(output.me)
			
		}else
		{
			vec.align <- c("l", rep("r", times=dim(output.x.continuous)[2]))
			output.x.continuous.xtable <- xtable(output.x.continuous, digits = xtable.digits, display = xtable.display,
					align = vec.align,
					caption = "Selected values of metric covariates",
					label = "tab-values")
			print(output.x.continuous.xtable, 
					format.args=list(big.mark = xtable.big.mark, decimal.mark = xtable.decimal.mark),
					booktabs = xtable.booktabs)
			
			vec.align <- c("l", rep("r", times=dim(output.effect.resp)[2]))
			output.effect.resp.xtable <- xtable(output.effect.resp, digits = xtable.digits, display = xtable.display,
					align = vec.align,
					caption = "Effects in different groups for selected values of metric covariates",
					label = "tab-effects")
			print(output.effect.resp.xtable, 
					format.args=list(big.mark = xtable.big.mark, decimal.mark = xtable.decimal.mark),
					booktabs = xtable.booktabs,
					include.rownames = factors.present)
			
			vec.align <- c("l", rep("r", times=dim(output.me)[2]))
			output.me.xtable <- xtable(output.me, digits = xtable.digits, display = xtable.display,
					align = vec.align,
					caption = "Marginal effects in different groups for selected values of metric covariates",
					label = "tab-marginaleffects")
			print(output.me.xtable, hline.after=c(-1,0,2,nrow(output.me.xtable)), 
					format.args=list(big.mark = xtable.big.mark, decimal.mark = xtable.decimal.mark),
					booktabs = xtable.booktabs)
			
		}# end if latex2console
		
		panel	
	} # end body internal.snapshot
	
	#############################################	
	#  Determine panel dimensions and layout 	#
	#############################################	
	
	# Calculate slider positions and sizes
	y.pos.slider.other <- 2.5*dist.obj.height + slider.height
	
	if(factors.present)
	{
		# Determine maximum length of group name
		length.factorgroup.name <- NULL 
		for(i in 1:num.groups)
		{
			length.factorgroup.name <- c(length.factorgroup.name, length(strsplit(factor.comb$names[i], split=NULL)[[1]]))
		}
		max.length.factorgroup.name <- max(length.factorgroup.name)
	}else
	{#if no groups are present set to fix value
		max.length.factorgroup.name <- 20   
	}
	
	# Calculate button position and size
	button.width        <- max(box.group.character.width*max.length.factorgroup.name, box.group.character.width*18)
	button.pos.x        <- 2*dist.obj.width + slider.width
	button.pos.y        <- dist.obj.height
	
	# Calculate box positions and sizes
	boxes.width         <- button.width
	box.pos.x           <- button.pos.x
	box.type.pos.y      <- button.pos.y + button.height + dist.obj.height
	box.groups.pos.y    <- box.type.pos.y  + box.type.height + dist.obj.height
	
	if(factors.present)
	{
		box.group.height <- box.group.line.height*length(factor.comb$names) + 10
	}else
	{
		box.group.height <- 0   
	}
	
	# Calculate overall panel size
	overall.width   <- box.pos.x + dist.obj.width + boxes.width
	overall.height  <- max((4*dist.obj.height + box.type.height + box.group.height + button.height), ((num.continuous-1)*slider.height + y.pos.slider.other + dist.obj.height))
	
	#########################################################	
	#  Define main panel and populate it with GUI controls 	#
	#########################################################
	mainpanel <- rp.control(title=panel.title, size = c(overall.width,overall.height))
	
	# Slider for displayed continuous variable
	initial.act <- initial.values[[var.selection]]
	if(is.null(initial.act)) initial.act <- mean(x.actual, na.rm = TRUE) 
	
	eval(wrapRpSlider(panel="mainpanel", variable="var.selected.act", title=paste(label.slider.act, var.selection, sep=""), 
					from = min(x.actual, na.rm=TRUE), to=max(x.actual, na.rm=TRUE), initval=initial.act,
					pos= c(dist.obj.width, dist.obj.height, slider.width, slider.height) , showvalue =TRUE, action = "func.panel.action"))
	
	# Sliders for other continuous variables (if present)
	if(!single.covariate)
	{
		for(i in 1:num.continuous.other)
		{
			y.pos.slider <- y.pos.slider.other + (i-1)*slider.height
			
			initial.act <- initial.values[[nam.X.continuous.other[i]]]
			if(is.null(initial.act)) initial.act <- mean(X.continuous.other[,i], na.rm = TRUE)
			
			slidercall <- wrapRpSlider(panel="mainpanel", variable=nam.X.continuous.other[i],title=nam.X.continuous.other[i], 
					from = min(X.continuous.other[,i], na.rm=TRUE), to=max(X.continuous.other[,i], na.rm=TRUE), initval=initial.act,
					pos= c(dist.obj.width, y.pos.slider, slider.width, slider.height) , showvalue =TRUE, action = "func.panel.action") 
			eval(slidercall)
		}
	}
	
	# Button for snapshot
	rp.button(panel=mainpanel, action = internal.snapshot, title = label.button , pos = c(button.pos.x, button.pos.y, button.width, button.height))
	
	# Radiogroup for type
	# if type is preselected check for allowed value
	if(any(preselect.type==c("link","response","marginal")))
	{# if allowed use it
		type.initval <- preselect.type[1]
	}else
	{# if not allowed set as link
		type.initval <- "link"	
	}
	
	# Cheat R CMD check	
	radiogroup.call <- "rp.radiogroup(panel=mainpanel, variable=type, vals=c('link', 'response', 'marginal'), initval=type.initval, labels=label.types,
			title = label.box.type, action = func.panel.action, pos = c(box.pos.x, box.type.pos.y, boxes.width, box.type.height))"
	eval(parse(text=radiogroup.call))
	
	# Checkbox for groups, if factors are employed in the model
	if(factors.present)
	{
		if(is.null(preselect.groups))
		{
			init.logical.index.groups  <- rep(TRUE, times=num.groups)
		}else
		{
			init.logical.index.groups  <- rep(FALSE, times=num.groups)
			init.logical.index.groups[preselect.groups]  <- TRUE
		}	
		
		# Cheat R CMD check	
		checkbox.call <- "rp.checkbox(panel=mainpanel, variable=groups, action = func.panel.action, initval=init.logical.index.groups,
				labels = factor.comb$names, title = label.box.groups, pos = c(box.pos.x, box.groups.pos.y, boxes.width, box.group.height))"
		eval(parse(text=checkbox.call)) 
	}		
	
	# Initalize panel
	rp.do(panel=mainpanel, func.panel.action)
	
	# when autosave is activated close panel after initialization
	if(autosave.plot) rp.control.dispose(panel=mainpanel)
} 


wrapRpSlider <- function(panel, variable, from, to, action, title=NA, 
		log = FALSE, showvalue = FALSE, resolution = 0, initval = NULL, 
		pos = NULL, horizontal = TRUE)
{
	if(is.na(title))
	{
		title.call <- paste("'",variable,"'", sep="")
	}else
	{
		title.call <- paste("'",title,"'", sep="")
	}
	
	slider.call <- paste("rp.slider(panel=", panel,", variable=",variable, ", from=",from,", to=",to,", action=", action,", title=",title.call,", pos=",deparse(pos),
			", log =",deparse(log),", showvalue = ",deparse(showvalue),", resolution = ",resolution,", initval = ",deparse(initval),
			", horizontal = ",deparse(horizontal),")", sep="")
	
	return(parse(text=slider.call)) 
}



factorCombinations <- function(X, factor.sep="|", level.sep=".", count=TRUE)
{
	logicalindex.factors <- sapply(X, is.factor)
	if(!any(logicalindex.factors)) stop("no factors in data.frame")
	
	X.factor <- X[,logicalindex.factors, drop=FALSE]
	if(dim(X.factor)[2]==1)
	{# just a single factor in the covariates, groups are the factor levels 
		factorgroup <- data.frame(levels(X.factor[,1]))
		colnames(factorgroup) <- colnames(X.factor)
		
		groups.names <- paste(colnames(factorgroup)[1], factorgroup[,1], sep=level.sep)
		
		# count occurences of levels
		groups.counts <- as.vector(table(X.factor))
		
		output <- list(comb=factorgroup, names=groups.names, counts=groups.counts)
		return(output)
	}
	
	factor.levels 	  <- lapply(X.factor, levels)
	num.factor.levels <- sapply(factor.levels, length)
	num.factors       <- dim(X.factor)[2]
	
	# first factor
	fac.act <- rep(factor.levels[[1]], each=prod(tail(num.factor.levels,-1)))
	factorgroups <- data.frame(fac.act)
	
	for(j in c(2:num.factors))
	{
		if(j==num.factors)
		{# for the last factor
			num.rep <- 1	
		}else
		{	
			num.rep      <- prod(tail(num.factor.levels,-j))
		}
		
		fac.act      <- rep(factor.levels[[j]], each=num.rep)
		factorgroups <- cbind(factorgroups, fac.act)
	}
	
	colnames(factorgroups) <- colnames(X.factor)
	
	# check for orderd factors and if so assign order to corresponding column
	for(j in 1:dim(X.factor)[2])
	{
		if(is.ordered(X.factor[,j]))
		{
			factorgroups[,j] <- factor(factorgroups[,j], levels=levels(X.factor[,j]), ordered = TRUE)
		}
	}
	
	# build group names from factor combinations
	groups.names <- as.vector(NULL)
	for(j in 1:dim(factorgroups)[2])
	{
		groups.names.actual <- paste(colnames(factorgroups)[j], factorgroups[,j], sep=level.sep)
		groups.names <- paste(groups.names, groups.names.actual, sep=factor.sep)
	}
	groups.names <- substring(groups.names,2)
	
	# count occurences of groups
	if(count)
	{	
		X.factor.complete <- X.factor[complete.cases(X.factor),,drop=FALSE]
		
		groups.data <- as.vector(NULL)
		for(j in 1:dim(X.factor.complete)[2])
		{
			groups.data.actual <- paste(colnames(X.factor.complete)[j], X.factor.complete[,j], sep=level.sep)
			groups.data <- paste(groups.data, groups.data.actual, sep=factor.sep)
		}
		groups.data <- substring(groups.data,2)
		
		groups.counts <- sapply(groups.names, function(x) sum(x==groups.data), USE.NAMES = FALSE)
	}else
	{
		groups.counts <- NULL
	}
	
	output <- list(combinations=factorgroups, names=groups.names, counts=groups.counts)
	return(output)
}


fxInteractive.lme <- function(model, # object of class lm is mandatory
		
		predict.lme.level = 0, # Levels of grouping to be displayed, passed to predict.lme. Only one level can be visualized, default to 0. 
				
		## Control initial appearance ##
		initial.values           = as.list(NULL), # Initial values for the metric covariates in a named list, default to the means. See details.
		preselect.var            = NA, # Name of continuous variable to be displayed as character or NA for menu selection, default to NA.
		preselect.type           = "effect",
		preselect.groups 		 = NULL, # Index of groups to be preselectes. If set to NULL (default) all groups are preselected.
		
		## Parameters for plot layout ##
		dev.height          = 18,	 # Height of graphic device in cm.
		dev.width           = 18,	 # Width of plot area in graphic device in cm.
		dev.width.legend    = 8,	 # Width of legend area in graphic device in cm.
		dev.pointsize		= 10,	 # Character pointsize of graphic device.
		dev.defined         = FALSE, # Graphic device predefined, e. g. customized for printing? Default to FALSE, see details.
		ylim                = NA,    # With a numeric vector of length 2 the plot limits in y-direction can be set. If NA (the default) these are determined automatically.
		col					= NA,	 # Vector of color specifications for the groups. Passed to the line commands and to the legend. Actual palette and consecutive sequence if NA (default).
		lty					= 1,     # Vector of line type specifications for the groups. Passed to the line commands and to the legend, default to solid lines.
		lwd					= 1,     # Vector of line width specifications for the groups. Passed to the line commands and to the legend, default to 1.
		main				= NA,	 # Label for the plot title.
		main.line			= 1.5,	 # Height in lines for plot title which is passed to title(), default to 1.5.
		xlab				= NA,	 # Label for the x-axis. Name of the selected covariate, if NA.
		ylab				= NA,	 # Label for the y-axis. Name of the selected plot type (see argument label.types), if NA.
		legend.add          = TRUE,  # Should a legend be added to the plot? Default to TRUE, but too many groups can cause problems.
		legend.space        = legend.add, # Should the space for the legend be reserved? Default to the value of legend.add. Setting legend.add to FALSE and legend.space to TRUE plots white space instead of the legend. This can be useful if different termplots will be arranged in a document.  
		legend.only         = FALSE,	  # Plot the legend alone.
		legend.pos			= "center",   # Position of the legend, see legend() for details.
		legend.cex          = 1,	      # Relative size of legend text, reduce for many groups.
		legend.width.factor = 1,	   # Factor by which the width of legend box is manipulated.	
		rug.ticksize		= 0.02,    # Length of rugplot tickmarks. Set to 0, if no rugplot should be drawn.
		rug.col				= "black", # Color of rugplot tickmarks, default to black.	
		vline.actual        = TRUE,    # Add vertical line at actual postion of selceted metric covariate? Default to TRUE.
		pos.hlines       	= c(0,0),  # Positon of horizontal line for [1] effect plot and [2] for marginal effects plot, NA for no lines.
		n.effects           = 100,     # Number of equally space points over the span of the metric covariate used for plotting.
		
		## Parameters for plot snapshot
		autosave.plot       = FALSE,    # Save the initial plot?
		snapshot.plot       = FALSE,    # Save plot as PDF when snapshot button is pressed? Default to FALSE.
		graphics.filename   = "LinRegIntPlot",  # Filename as character for graphic file.
		graphics.numbering  = !autosave.plot,   # Automatically append a 3-digits-number to the filenname to avoid that existing graphic files are overwritten.
		graphics.type       = "pdf",    # Graphics file type.
		
		## Parameters for text-output ##
		factor.sep          = "|",   # Character by which the factors are separated in the groupnames.
		level.sep           = ".",   # Character by which the levels are separated in the groupnames.
		latex2console       = FALSE, # Should the textoutput triggered by the snapshot button be printed as LaTeX?
		xtable.big.mark 	= ".",   # Bigmark character for LaTeX output, argument passed to print.xtable().
		xtable.decimal.mark = ",",   # Decimal character for LaTeX output, argument passed to print.xtable().
		xtable.digits 		= NULL,  # Number of digits, argument passed to xtable().
		xtable.display 		= NULL,  # Display style, argument passed to xtable().
		xtable.booktabs 	= FALSE, # Use LaTeX package booktabs for horizontal lines, argument passed to print.xtable().
		
		## Annotations for panel 
		panel.title         = "Linear Mixed-Effects Model",		# Title used in the panel.
		label.button        = "Snapshot" ,						# Label for the snapshot-button.
		label.slider.act    = "Variable displayed: " ,			# Additional label for the slider of the selected metric covariate. 
		label.box.type      = "Type" ,							# Title for the radiogroup box.
		label.types         = c("effect", "marginal effect"),   # Lables for radiogroup buttons (character vector of length 2).
		label.box.groups    = "Groups" ,						# Title for the checkbox.
		
		## Parameters to control the size of the panel.
		slider.width                = 200 , # Width of each slider.    
		slider.height               = 60  , # Height of each slider.   
		button.height               = 30  , # Height of snapshot button.
		box.type.height             = 70  , # Height of radiobox for type selection.
		box.group.character.width   = 7   , # The width of the boxes is basically a this value times the number of characters.
		box.group.line.height       = 28  , # The height of the checkbox is this value times the number of groups.
		dist.obj.width              = 20  , # Vertical distance between sliders and boxes and vertical margins. 
		dist.obj.height             = 10  , # Horizontal distance between panel objects.
		
		... ) # other graphical parameters passed to par()	
{
	#################################################################	
	# Pick covariates employed, check variables for factors, assign #
	# corresponding objects and initialize graphic window			#
	#################################################################	
	
	# pick covariates which are employed in the model
	if(!is.null(model$data))
	{# extract data from model$data if possible
		X <- get_all_vars(model$terms, model$data)[,-1,drop=FALSE]
	}else
	{# try to extract data from model$model
		X <- get_all_vars(model$terms, model$model)[,-1,drop=FALSE]	
	}
	
	X <- cbind(X, model$data[,colnames(model$groups), drop=FALSE])
	
	# identify factors from model matrix
	logicalindex.factor <- sapply(X, is.factor)
	factors.present      <- any(logicalindex.factor)
	
	if(factors.present)
	{
		num.level <- sapply(X[,logicalindex.factor,drop=FALSE], function(x) length(levels(x)))
	}
	
	# separate factors and continuous covariates 
	X.continuous   <- X[,!logicalindex.factor, drop=FALSE]
	num.continuous <- dim(X.continuous)[2]
	
	# switch for special treatment of single metric covariate
	single.covariate <- num.continuous==1 
	
	if(factors.present)
	{
		# separate factors 
		X.factor    <- X[,logicalindex.factor, drop=FALSE]
		num.groups  <- prod(num.level)
		
		# build all factor combinations
		factor.comb <- factorCombinations(X.factor, factor.sep=factor.sep, level.sep=level.sep, count=FALSE)
	}
	
	# If not specified in the function call: if there is only one continuous covariate, choose it. 
	# Otherwise select the continuous covariate to be displayed from popup list.
	if(is.na(preselect.var))
	{
		if(single.covariate)
		{
			var.selection <- colnames(X.continuous)[1]
		}else
		{	
			var.selection <- select.list(colnames(X.continuous))
		}
	}else
	{# in no variable is preselected
		if(single.covariate)
		{# if there is only one metric covariate, chose it
			var.selection <- colnames(X.continuous)[1]
		}else
		{# Check if preselected covariate exists, if not, show the menu	
			if(any(colnames(X.continuous)==preselect.var))
			{
				var.selection <- preselect.var		
			}else
			{
				var.selection <- select.list(colnames(X.continuous))	
			}	
		}
	}
	
	# If no variable is selected from popup list leave function.
	if(var.selection=="") return()
	
	# pick continuous covariate under invstigation
	x.actual           <- X.continuous[,var.selection]
	x.actual.sequence  <- seq(min(x.actual, na.rm=TRUE), max(x.actual, na.rm=TRUE), length=n.effects)
	X.pred             <- data.frame(x=x.actual.sequence)
	colnames(X.pred)   <- var.selection
	
	# pick other continuous covariates if present
	if(!single.covariate)
	{	
		X.continuous.other     <- X.continuous[,-match(var.selection, colnames(X.continuous)), drop=FALSE]
		nam.X.continuous.other <- colnames(X.continuous.other)
		num.continuous.other   <- length(nam.X.continuous.other)
		val.X.continuous.other <- colMeans(X.continuous.other)
	}
	
	# If not predefined nor switched off via the argument: set up a new graphic device 
	if(!dev.defined)
	{	
		if(legend.only)
		{# Initialize device for legend only
			dev.new(width=dev.width.legend/2.54, height=dev.height/2.54, pointsize=dev.pointsize, buffered=TRUE, noRStudioGD = TRUE)
			par(cex=1, ...)
		}else
		{# If factors are used as covariates and legend should be printed, split plot region via layout()	
			if(factors.present & legend.space)
			{
				dev.new(width=(dev.width + dev.width.legend)/2.54, height=dev.height/2.54, pointsize=dev.pointsize, buffered=TRUE, noRStudioGD = TRUE)
				fraction.plot <- round((dev.width/(dev.width + dev.width.legend))*100, digits=0)
				layoutmatrix  <- matrix(c(rep(2, fraction.plot), rep(1, (100-fraction.plot))),1,100)
				layout(layoutmatrix)
				par(cex=1, ...)
			}else
			{
				dev.new(width=dev.width/2.54, height=dev.height/2.54, pointsize=dev.pointsize, buffered=TRUE, noRStudioGD = TRUE)
				par(cex=1, ...)
			}
		}
	}
	
	#############################################	
	#  Action when panel items are manipulated 	#
	#  (except action for snapshot-button)		#	
	#############################################
	func.panel.action <- function(panel) 
	{# Main panel function, must be defined within the namespace where panel is established
		# read plot type from panel
		if(length(panel$type)==0) # in the first call of the panel the slots are sometimes NULL   
		{
			type.actual <- "effect"
		}else
		{
			type.actual <- panel$type
		}
		
		# read group selection from panel if groups are present
		if(factors.present)
		{
			if(length(panel$groups)==0) 
			{# in the first call of the panel the slots are sometimes NULL 
				if(is.null(preselect.groups))
				{
					logical.index.groups  <- rep(TRUE, times=num.groups)
				}else
				{
					logical.index.groups  <- rep(FALSE, times=num.groups)
					logical.index.groups[preselect.groups]  <- TRUE
				}
				index.groups.selected <- c(1:num.groups)[logical.index.groups]
				# if all groups are deselected choose the first
				if(sum(logical.index.groups)==0)
				{
					index.groups.selected  <- 1  
				}
			}else
			{
				logical.index.groups  <- panel$groups
				index.groups.selected <- c(1:num.groups)[logical.index.groups]
				
				# if all groups are deselected choose the first
				if(sum(logical.index.groups)==0)
				{
					index.groups.selected  <- 1  
				}
			}
		}# end if factors.present	
		
		# read value of actual variable from panel
		if(length(panel$var.selected.act)==0)
		{
			var.selected.act <- 1
		}else
		{
			var.selected.act <- panel$var.selected.act      
		}   
		
		# build designmatrix for predict-method from continuous covariates
		if(single.covariate)
		{
			X.pred.metr <- X.pred			
		}else
		{# read actual values from panel and add to designmatrix for predict-method
			for(i in 1:num.continuous.other)
			{
				if(length(panel[[nam.X.continuous.other[i]]])==0) # initialization fails
				{
					val.X.continuous.other[i] <- i
				}else
				{
					val.X.continuous.other[i] <- panel[[nam.X.continuous.other[i]]]
				}
			}
			X.pred.metr <- cbind(X.pred, t(val.X.continuous.other))
		}
		
		# If factors are used as covariates, calculate predictions for each group.
		if(factors.present)
		{
			# calculate predictions for each selected group
			list.groups <- as.list(NULL)
			for(i in seq(along=index.groups.selected))
			{
				index.group.actual <- index.groups.selected[i]  
				
				# read actual factor combination
				factor.group.actual <- factor.comb$comb[index.group.actual,,drop=FALSE]
				rownames(factor.group.actual) <- ""
				
				# combine with continuous covariates to new designmatrix
				X.pred <- cbind(X.pred.metr, factor.group.actual)
				
				# calculate effect or marginal effect
				switch(type.actual,
						effect = pred.x.actual <- predict(model, newdata=X.pred, level=predict.lme.level[1]),
						marginal = pred.x.actual <- splinefun(x=x.actual.sequence, y=predict(model, newdata=X.pred, level=predict.lme.level[1]))(x.actual.sequence, deriv=1))
				
				
				# store results for each and every group 
				sub.list <- list(name=factor.comb$names[index.group.actual],
						prog=pred.x.actual)
				
				list.groups[[i]] <- sub.list 
			}
		}else
		{
			# calculate effect or marginal effect
			switch(type.actual,
					effect = pred.x.actual <- predict(model, newdata=X.pred.metr, level=predict.lme.level[1]),
					marginal = pred.x.actual <- splinefun(x=x.actual.sequence, y=predict(model, newdata=X.pred.metr, level=predict.lme.level[1]))(x.actual.sequence, deriv=1))
			
			index.groups.selected <- 1
			
			# store results
			list.groups <- as.list(NULL)
			list.groups[[1]] <- list(name="default", prog=pred.x.actual)
		}
		
		# determine plot limits in y-direction
		if(any(is.na(ylim))) ylim <- c(min(unlist(lapply(list.groups,"[",2)), na.rm=TRUE), max(unlist(lapply(list.groups,"[",2)), na.rm=TRUE))		
		
		
		### Plotting commands ###
		func.ploteffects <- function()
		{
			# draw effects for each and every group
			# when no factors are present set number of groups to 1 for correct color handling
			if(!factors.present) num.groups <- 1
			# specify color scheme
			if(all(is.na(col))) 
			{
				col.types <- c(1:num.groups)
			}else
			{	
				col.types <- rep(col, times=ceiling(num.groups/length(col)))
			}
			
			lty.types <- rep(lty, times=ceiling(num.groups/length(lty)))
			lwd.types <- rep(lwd, times=ceiling(num.groups/length(lwd)))
			
			if(factors.present & legend.only)
			{# When legend.only is active plot legend only
				old.mar <- par("mar")
				par(mar=c(old.mar[1],0,old.mar[3],0))
				plot(1,1, type="n", axes=FALSE, xlab=NA, ylab=NA)
				legend.labels <- paste(factor.comb$names[index.groups.selected],"", sep="")
				legend(x=legend.pos, legend=legend.labels, text.width=max(strwidth(legend.labels, cex=legend.cex))*legend.width.factor, col=col.types[index.groups.selected], lty=lty.types[index.groups.selected], lwd=lwd.types[index.groups.selected], cex=legend.cex)
				par(mar=old.mar)    
			}else	
			{# when inactive incorporate other legend parameters
				# when factors are used as covariates and switch is on: add legend first to allow adding elements to the main plot	
				if(factors.present & legend.space)
				{
					old.mar <- par("mar")
					par(mar=c(old.mar[1],0,old.mar[3],0))
					plot(1,1, type="n", axes=FALSE, xlab=NA, ylab=NA)
					if(legend.add)
					{
						legend.labels <- paste(factor.comb$names[index.groups.selected],"", sep="")
						legend(x=legend.pos, legend=legend.labels, text.width=max(strwidth(legend.labels, cex=legend.cex))*legend.width.factor, col=col.types[index.groups.selected], lty=lty.types[index.groups.selected], lwd=lwd.types[index.groups.selected], cex=legend.cex)
					}
					par(mar=old.mar)    
				}		
				
				for(i in seq(along=index.groups.selected))
				{
					# initial plot in first loop	
					if(i==1)
					{
						if(is.na(xlab)) xlab <- var.selection
						if(is.na(ylab))
							switch(type.actual,
									effect   = ylab <- label.types[1],
									marginal = ylab <- label.types[2])
						
						
						plot(x.actual.sequence, rep(1, times=length(x.actual.sequence)), type="n", ylim=ylim, xlab=xlab, ylab=ylab, main=NA)
						title(main=main, line=main.line)
						if(!((rug.ticksize==0)||is.na(rug.ticksize))) rug(x.actual,  ticksize = rug.ticksize, col= rug.col) 
					}
					
					index.group.actual <- index.groups.selected[i]  
					
					# if there is only a litte bit variation in the effect draw a horizontal line
					if(sd(list.groups[[i]]$prog, na.rm = TRUE)<10^-10)
					{
						abline(h=mean(list.groups[[i]]$prog, na.rm = TRUE), col=col.types[index.group.actual], lty=lty.types[index.group.actual], lwd=lwd.types[index.group.actual])	
					}else	
					{ 
						lines(x.actual.sequence, list.groups[[i]]$prog, type="l", col=col.types[index.group.actual], lty=lty.types[index.group.actual], lwd=lwd.types[index.group.actual])
					}     
				}
				
				# vertical line at actual position
				if(vline.actual) abline(v=(var.selected.act))
				
				# draw horizontal line at specified positions
				switch(type.actual,
						effect 	 = abline(h=pos.hlines[1]),
						marginal = abline(h=pos.hlines[2]))	
			}# end if legend.only		
		}# end func.ploteffects
		
		# In conjunction with dev.flush() and buffered=TRUE animation is more fluent.
		dev.hold()
		func.ploteffects()
		dev.flush()
		
		# When autosave is activated directly save the plot
		if(autosave.plot)
		{
			# Create filename for plot
			if(graphics.numbering)
			{
				for(i in 1:1000)
				{
					graphics.filename.candidate <- paste(graphics.filename,"-", formatC(i, width=3, flag="0"), sep="")
					if(file.exists(paste(graphics.filename.candidate, ".", graphics.type, sep=""))) next else break
				}
			}else
			{
				graphics.filename.candidate <- paste(graphics.filename, ".", graphics.type, sep="") 
			}
			
			# Platform dependent save operation
			if(.Platform$OS.type != "windows") 
			{ # Mac OS, Linux
				if (any(graphics.type == c("png","jpeg","jpg","tiff","bmp")))
				{
					sptype <- graphics.type
					if (graphics.type == "jpg") {sptype <- "jpeg"}
					savePlot(filename=graphics.filename.candidate, type=sptype , ... )      
				} 
				if(graphics.type == "pdf")
				{
					dev.copy2pdf(file=graphics.filename.candidate)
				}
				if(graphics.type == "eps")
				{
					dev.copy2eps(file=graphics.filename.candidate)
				}
			}else
			{ # Windows OS
				savePlot(filename = graphics.filename.candidate, type = graphics.type)	
			}
		}
		
		panel
	} # end body action()
	
	#############################################	
	#    Action when snapshot-button is used  	#
	#############################################
	internal.snapshot <- function(panel)
	{# function must be defined within the namespace where panel is established
		# save plot if snapshot.plot is TRUE
		if(snapshot.plot)      
		{               
			# Create filename for plot
			if(graphics.numbering)
			{
				for(i in 1:1000)
				{
					graphics.filename.candidate <- paste(graphics.filename,"-", formatC(i, width=3, flag="0"), sep="")
					if(file.exists(paste(graphics.filename.candidate, ".", graphics.type, sep=""))) next else break
				}
			}else
			{
				graphics.filename.candidate <- paste(graphics.filename, ".", graphics.type, sep="") 
			}
			
			# Platform dependent save operation
			if(.Platform$OS.type != "windows") 
			{ # Mac OS, Linux
				if (any(graphics.type == c("png","jpeg","jpg","tiff","bmp")))
				{
					sptype <- graphics.type
					if (graphics.type == "jpg") {sptype <- "jpeg"}
					savePlot(filename=graphics.filename.candidate, type=sptype , ... )      
				} 
				if(graphics.type == "pdf")
				{
					dev.copy2pdf(file=graphics.filename.candidate)
				}
				if(graphics.type == "eps")
				{
					dev.copy2eps(file=graphics.filename.candidate)
				}
			}else
			{ # Windows OS
				savePlot(filename = graphics.filename.candidate, type = graphics.type)	
			}
		}
		
		# read values from panel and build design matrix for predict-method
		var.selected.act <- panel$var.selected.act
		if(single.covariate)
		{
			X.pred.metr <- data.frame(var.selected.act)
			colnames(X.pred.metr)[1] <- var.selection
			vek.x.continuous.actual <- var.selected.act
			names(vek.x.continuous.actual) <- var.selection
		}else
		{	
			for(i in 1:num.continuous.other)
			{
				val.X.continuous.other[i] <- panel[[nam.X.continuous.other[i]]]
			}
			
			# build designmatrix from continuous covariates
			X.pred.metr <- cbind(var.selected.act, t(val.X.continuous.other))
			colnames(X.pred.metr)[1] <- var.selection
			vek.x.continuous.actual <- as.vector(X.pred.metr)
			names(vek.x.continuous.actual) <- colnames(X.pred.metr)
		}
		
		if(factors.present)
		{
			X.pred <- cbind(X.pred.metr, factor.comb$comb)
		}else
		{
			X.pred <- as.data.frame(X.pred.metr)
		}
		
		output.effect.resp <- cbind(predict(model, newdata=X.pred, level=predict.lme.level[1]),1)[,-2, drop=FALSE]		
		colnames(output.effect.resp)  <- c("effect")
		
		if(factors.present)
		{
			rownames(output.effect.resp)  <- factor.comb$names
		}else
		{
			rownames(output.effect.resp)  <- ""   
		}
		
		# calculate ECDF-value for each continuous covariate and  marginal effects in each group
		# Order corresponds to appearance in the original data.frame
		output.vek.x.continuous.actual <- NULL
		F.x.continuous 				   <- NULL
		output.me 		               <- NULL
		
		for(i in 1:num.continuous)
		{
			nam.x.marginal.actual       <- colnames(X.continuous)[i]
			
			x.marginal.actual           <- X.continuous[,nam.x.marginal.actual, drop=TRUE]
			x.marginal.actual.sequence  <- seq(min(x.marginal.actual, na.rm=TRUE), max(x.marginal.actual, na.rm=TRUE), length=n.effects)
			X.marginal.pred             <- data.frame(x=x.marginal.actual.sequence)
			colnames(X.marginal.pred)   <- nam.x.marginal.actual
			
			if(!single.covariate)
			{	
				nam.X.marginal.other        <- colnames(X.continuous)[-i]
				X.marginal.pred             <- cbind(X.marginal.pred, X.pred.metr[,nam.X.marginal.other, drop=FALSE])
				colnames(X.marginal.pred)   <- c(nam.x.marginal.actual, nam.X.marginal.other)   
			}
			
			x.marginal.actual.selected     <- X.pred.metr[,nam.x.marginal.actual, drop=TRUE]
			
			output.vek.x.continuous.actual <- c(output.vek.x.continuous.actual, x.marginal.actual.selected)
			
			# ECDF-value	
			F.x.continuous.actual        <- sum(x.marginal.actual <= x.marginal.actual.selected)/length(x.marginal.actual)
			F.x.continuous               <- c(F.x.continuous, F.x.continuous.actual)
			
			# calculate marginals for actual continuous covariate in each group
			marginal.actual.groups <- NULL
			if(factors.present)
			{
				for(i in 1:num.groups)
				{
					# read actual factor combination
					factor.group.actual <- factor.comb$comb[i,,drop=FALSE]
					rownames(factor.group.actual) <- ""
					
					# combine with continuous covariates to new designmatrix
					X.marginal.pred.actual <- cbind(X.marginal.pred, factor.group.actual)
					
					marginal.actual.groups.actual   <- splinefun(x=x.marginal.actual.sequence, y=predict(model, newdata=X.marginal.pred.actual, level=predict.lme.level[1]))(x.marginal.actual.selected, deriv=1)
					marginal.actual.groups          <- c(marginal.actual.groups, marginal.actual.groups.actual)
				}   
			}else
			{
				marginal.actual.groups  <-	splinefun(x=x.marginal.actual.sequence, y=predict(model, newdata=X.marginal.pred, level=predict.lme.level[1]))(x.marginal.actual.selected, deriv=1) 
			}
			
			output.me <- cbind(output.me, marginal.actual.groups)
		}#end for-loop continuous covariates
		
		output.x.continuous                   <- rbind(output.vek.x.continuous.actual, F.x.continuous)
		colnames(output.x.continuous)         <- colnames(X.continuous)
		rownames(output.x.continuous)[c(1:2)] <- c("value", "ECDF(value)")
		
		if(factors.present)
		{
			rownames(output.me) <- factor.comb$names
		}else
		{
			rownames(output.me) <- "marginal effect"   
		}
		
		# Add actual values to table of marginal effects
		colnames(output.me) <- colnames(X.continuous)
		output.me           <- rbind(output.x.continuous, output.me)
		
		if(!latex2console)
		{
			print(summary(model))
			cat(fill=TRUE)
			cat("Selected values of metric covariates", fill=TRUE)
			print(output.x.continuous, fill=TRUE)
			
			cat(fill=TRUE)
			cat("Effects in different groups for selected values of metric covariates", fill=TRUE)
			print(output.effect.resp)
			
			cat(fill=TRUE)
			cat("Marginal effects in different groups for selected values of metric covariates", fill=TRUE)
			print(output.me)
			
		}else
		{
			vec.align <- c("l", rep("r", times=dim(output.x.continuous)[2]))
			output.x.continuous.xtable <- xtable(output.x.continuous, digits = xtable.digits, display = xtable.display,
					align = vec.align,
					caption = "Selected values of metric covariates",
					label = "tab-values")
			print(output.x.continuous.xtable, 
					format.args=list(big.mark = xtable.big.mark, decimal.mark = xtable.decimal.mark),
					booktabs = xtable.booktabs)
			
			vec.align <- c("l", rep("r", times=dim(output.effect.resp)[2]))
			output.effect.resp.xtable <- xtable(output.effect.resp, digits = xtable.digits, display = xtable.display,
					align = vec.align,
					caption = "Effects in different groups for selected values of metric covariates",
					label = "tab-effects")
			print(output.effect.resp.xtable, 
					format.args=list(big.mark = xtable.big.mark, decimal.mark = xtable.decimal.mark),
					booktabs = xtable.booktabs,
					include.rownames = factors.present)
			
			vec.align <- c("l", rep("r", times=dim(output.me)[2]))
			output.me.xtable <- xtable(output.me, digits = xtable.digits, display = xtable.display,
					align = vec.align,
					caption = "Marginal effects in different groups for selected values of metric covariates",
					label = "tab-marginaleffects")
			print(output.me.xtable, hline.after=c(-1,0,2,nrow(output.me.xtable)), 
					format.args=list(big.mark = xtable.big.mark, decimal.mark = xtable.decimal.mark),
					booktabs = xtable.booktabs)
			
		}# end if latex2console
		panel
	}# end body internal.snapshot
	
	#############################################	
	#  Determine panel dimensions and layout 	#
	#############################################	
	
	# Calculate slider positions and sizes
	y.pos.slider.other <- 2.5*dist.obj.height + slider.height
	
	if(factors.present)
	{
		# Determine maximum length of group name
		length.factorgroup.name <- NULL 
		for(i in 1:num.groups)
		{
			length.factorgroup.name <- c(length.factorgroup.name, length(strsplit(factor.comb$names[i], split=NULL)[[1]]))
		}
		max.length.factorgroup.name <- max(length.factorgroup.name)
	}else
	{#if no groups are present set to fix value
		max.length.factorgroup.name <- 20   
	}
	
	# Calculate button position and size
	button.width        <- max(box.group.character.width*max.length.factorgroup.name, box.group.character.width*18)
	button.pos.x        <- 2*dist.obj.width + slider.width
	button.pos.y        <- dist.obj.height
	
	# Calculate box positions and sizes
	boxes.width         <- button.width
	box.pos.x           <- button.pos.x
	box.type.pos.y      <- button.pos.y + button.height + dist.obj.height
	box.groups.pos.y    <- box.type.pos.y  + box.type.height + dist.obj.height
	
	if(factors.present)
	{
		box.group.height <- box.group.line.height*length(factor.comb$names)
	}else
	{
		box.group.height <- 0   
	}
	
	# Calculate overall panel size
	overall.width   <- box.pos.x + dist.obj.width + boxes.width
	overall.height  <- max((4*dist.obj.height + box.type.height + box.group.height + button.height), ((num.continuous-1)*slider.height + y.pos.slider.other + dist.obj.height))
	
	#########################################################	
	#  Define main panel and populate it with GUI controls 	#
	#########################################################
	mainpanel <- rp.control(title=panel.title, size = c(overall.width,overall.height))
	
	# Slider for displayed continuous variable
	initial.act <- initial.values[[var.selection]]
	if(is.null(initial.act)) initial.act <- mean(x.actual, na.rm = TRUE) 
	
	eval(wrapRpSlider(panel="mainpanel", variable="var.selected.act", title=paste(label.slider.act, var.selection, sep=""), 
					from = min(x.actual, na.rm=TRUE), to=max(x.actual, na.rm=TRUE), initval=initial.act,
					pos= c(dist.obj.width, dist.obj.height, slider.width, slider.height) , showvalue =TRUE, action = "func.panel.action"))
	
	# Sliders for other continuous variables (if present)
	if(!single.covariate)
	{
		for(i in 1:num.continuous.other)
		{
			y.pos.slider <- y.pos.slider.other + (i-1)*slider.height
			
			initial.act <- initial.values[[nam.X.continuous.other[i]]]
			if(is.null(initial.act)) initial.act <- mean(X.continuous.other[,i], na.rm = TRUE)
			
			slidercall <- wrapRpSlider(panel="mainpanel", variable=nam.X.continuous.other[i],title=nam.X.continuous.other[i], 
					from = min(X.continuous.other[,i], na.rm=TRUE), to=max(X.continuous.other[,i], na.rm=TRUE), initval=initial.act,
					pos= c(dist.obj.width, y.pos.slider, slider.width, slider.height) , showvalue =TRUE, action = "func.panel.action")
			eval(slidercall)
		}
	}
	
	# Button for snapshot
	rp.button(panel=mainpanel, action = internal.snapshot, title = label.button , pos = c(button.pos.x, button.pos.y, button.width, button.height))
	
	# Radiogroup for type
	# if type is preselected check for allowed value
	if(any(preselect.type==c("effect", "marginal")))
	{# if allowed use it
		type.initval <- preselect.type[1]
	}else
	{# if not allowed set as effect
		type.initval <- "effect"	
	}
	
	# Cheat R CMD check	
	radiogroup.call <- "rp.radiogroup(panel=mainpanel, variable=type, vals=c('effect', 'marginal'), initval=type.initval, labels=label.types, 
			title = label.box.type, action = func.panel.action, pos = c(box.pos.x, box.type.pos.y, boxes.width, box.type.height))"
	eval(parse(text=radiogroup.call)) 
	
	# Checkbox for groups, if factors are employed in the model
	if(factors.present)
	{
		if(is.null(preselect.groups))
		{
			init.logical.index.groups  <- rep(TRUE, times=num.groups)
		}else
		{
			init.logical.index.groups  <- rep(FALSE, times=num.groups)
			init.logical.index.groups[preselect.groups]  <- TRUE
		}
		
		# Cheat R CMD check	
		checkbox.call <- "rp.checkbox(panel=mainpanel, variable=groups, action = func.panel.action, initval=init.logical.index.groups,
				labels = factor.comb$names, title = label.box.groups, pos = c(box.pos.x, box.groups.pos.y, boxes.width, box.group.height))"
		eval(parse(text=checkbox.call)) 
	}
	
	# Initalize panel
	rp.do(panel=mainpanel, func.panel.action)
	
	# when autosave is activated close panel after initialization
	if(autosave.plot) rp.control.dispose(panel=mainpanel)
}
