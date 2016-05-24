# some helpful threads
# https://stat.ethz.ch/pipermail/r-help/2008-September/172641.html
# http://tolstoy.newcastle.edu.au/R/e4/help/08/02/4875.html
# http://tolstoy.newcastle.edu.au/R/e2/help/07/01/8598.html

# http://www.r-statistics.com/wp-content/uploads/2011/01/boxplot-add-label-for-outliers.r.txt

bp.with.outlier.label <- function(y, label_name, ..., spread_text = T, data, plot = T, range = 1.5, label.col = "blue",
										jitter_if_duplicate = T, jitter_only_positive_duplicates = F)
{
	# change log:
	# 19.04.2011 - added support to "names" and "at" parameters.


	# jitter_if_duplicate - will jitter (Actually just add a bit of numbers) so to be able to decide on which location to plot the label when having identical variables...
	#require(plyr) #for is.formula and ddply

	# a function to jitter data in case of ties in Y's
	jitter.duplicate <- function(x, only_positive = F)
	{
		if(only_positive) {
			ss <- x > 0
		} else {
			ss <- T
		}
		ss_dup <- duplicated(x[ss])
		# ss <- ss & ss_dup
		temp_length <- length(x[ss][ss_dup])
		x[ss][ss_dup] <- x[ss][ss_dup] + seq(from = 0.00001, to = 0.00002, length.out = temp_length)
		x
	}
	# jitter.duplicate(c(1:5))
	# jitter.duplicate(c(1:5,5,2))
	# duplicated(jitter.duplicate(c(1:5,5,2)))
	# jitter.duplicate(c(0,0,1:5,5,2))
	# duplicated(jitter.duplicate(c(0,0,1:5,5,2)))



	# handle cases where
	if(jitter_if_duplicate) {
		# warning("duplicate jutter of values in y is ON")
		if(!missing(data)) {	#e.g: we DO have data
			# if(exists("y") && is.formula(y)) {		# F && NULL # F & NULL
			y_name <- as.character(substitute(y))	# I could have also used as.list(match.call())
												# credit to Uwe Ligges and Marc Schwartz for the help
												# https://mail.google.com/mail/?shva=1#inbox/12dd7ca2f9bfbc39
			if(length(y_name) > 1) {	# then it is a formula (for example: "~", "y", "x"
				model_frame_y <- model.frame(y, data = data)
				temp_y <- model_frame_y[,1]
				temp_y  <- jitter.duplicate(temp_y, jitter_only_positive_duplicates)	# notice that the default of the function is to work only with positive values...
				# the_txt <- paste(names(model_frame_y)[1], "temp_y", sep = "<<-") # wrong...
				the_txt <- paste("data['",names(model_frame_y)[1],"'] <- temp_y", sep = "")
				eval(parse(text = the_txt))	# jutter out y var so to be able to handle identical values.
			} else {	# this isn't a formula
				data[,y_name] <- jitter.duplicate(data[,y_name], jitter_only_positive_duplicates)
				y <- data[,y_name]	# this will make it possible for boxplot(y, data) to work later (since it is not supposed to work with data when it's not a formula, but now it does :))
			}
		} else {	# there is no "data"
			if(is.formula(y)) { # if(exists("y") && is.formula(y)) {		# F && NULL # F & NULL
				temp_y <- model.frame(y)[,1]
				temp_y  <- jitter.duplicate(temp_y, jitter_only_positive_duplicates)	# notice that the default of the function is to work only with positive values...
				temp_y_name <- names(model.frame(y))[1]	# we must extract the "names" before introducing a new enbironment (or there will be an error)
				environment(y) <- new.env()
				assign(temp_y_name, temp_y, environment(y))
					# Credit and thanks for doing this goes to Niels Richard Hansen (2 Jan 30, 2011)
					# http://r.789695.n4.nabble.com/environment-question-changing-variables-from-a-formula-through-model-frame-td3246608.html
				# warning("Your original variable (in the global environemnt) was just jittered.")	# maybe I should add a user input before doing this....
				# the_txt <- paste(names(model_frame_y)[1], "temp_y", sep = "<<-")
				# eval(parse(text = the_txt))	# jutter out y var so to be able to handle identical values.
			} else {
				y <- jitter.duplicate(y, jitter_only_positive_duplicates)
			}
		}
	}
	# the_txt <- paste("print(",names(model_frame_y)[1], ")")
	# eval(parse(text = the_txt))	# jutter out y var so to be able to handle identical values.
	# print(ls())


	# y should be a formula of the type: y~x, y~a*b
	# or it could be simply y
	if(missing(data)) {
			boxdata <- boxplot(y, plot = plot,range = range ,...)
		} else {
			boxdata <- boxplot(y, plot = plot,data = data, range = range ,...)
		}
	if(length(boxdata$names) == 1 && boxdata$names =="") boxdata$names <- 1	# this is for cases of type: boxplot(y) (when there is no dependent group)
	if(length(boxdata$out) == 0 ) stop("No outliers detected for this boxplot")

	#### if(!missing(data)) attach(data)	# this might lead to problams I should check out for alternatives for using attach here...


	# creating a data.frame with information from the boxplot output about the outliers (location and group)
	boxdata_group_name <- factor(boxdata$group)
	levels(boxdata_group_name) <- boxdata$names[as.numeric(levels(boxdata_group_name))]	# the subseting is for cases where we have some sub groups with no outliers
	if(!is.null(list(...)$at))	{	# if the user chose to use the "at" parameter, then we would like the function to still function (added on 19.04.2011)
		boxdata$group <- list(...)$at[boxdata$group]
		}
	boxdata_outlier_df <- data.frame(group = boxdata_group_name, y = boxdata$out, x = boxdata$group)


	# Let's extract the x,y variables from the formula:
	if(is.formula(y))
	{
		model_frame_y <- model.frame(y)
			# old solution: (which caused problems if we used the names parameter when using a 2 way formula... (since the order of the names is different then the levels order we get from using factor)
			# y <- model_frame_y[,1]
			# x <- model_frame_y[,-1]
		splited_model_frame_y <- split(model_frame_y[,1], model_frame_y[-1])
		y <- unlist(splited_model_frame_y)
		x <- rep(names(splited_model_frame_y), times = lapply(splited_model_frame_y, length)	)
		if(!is.null(dim(x))) {	# then x is a matrix/data.frame of the type x1*x2*..and so on - and we should merge all the variations...
			x <- apply(x,1, paste, collapse = ".")
		}
	} else {
		# if(missing(x)) x <- rep(1, length(y))
		x <- rep(1, length(y))	# we do this in case y comes as a vector and without x
	}

	# and put all the variables (x, y, and outlier label name) into one data.frame
	DATA <- data.frame(label_name, x ,y)

	if(!is.null(list(...)$names))	{	# if the user chose to use the names parameter, then we would like the function to still function (added on 19.04.2011)
		DATA$x <- factor(DATA$x, levels = unique(DATA$x))
		levels(DATA$x) = list(...)$names	# enable us to handle when the user adds the "names" parameter # fixed on 19.04.11	# notice that DATA$x must be of the "correct" order (that's why I used split above
		# warning("Careful, the use of the 'names' parameter is experimental.  If you notice any errors please e-mail me at: tal.galili@gmail.com")
		}

	#### if(!missing(data)) detach(data)	# we don't need to have "data" attached anymore.

	# let's only keep the rows with our outliers
	boxplot.outlier.data <- function(xx, y_name = "y")
	{
		y <- xx[,y_name]
		boxplot_range <- range(boxplot.stats(y, coef = range )$stats)
		ss <- (y < boxplot_range[1]) | (y > boxplot_range[2])
		return(xx[ss,])
	}
	outlier_df <-ddply(DATA, .(x), boxplot.outlier.data)


	# create propor x/y locations to handle over-laping dots...
	if(spread_text) {
		# credit: Greg Snow
	# require(TeachingDemos)
		temp_x <- boxdata_outlier_df[,"x"]
		temp_y1 <- boxdata_outlier_df[,"y"]
		temp_y2 <- temp_y1
		for(i in unique(temp_x))
		{
			tmp <- temp_x == i
			temp_y2[ tmp ] <- spread.labs( temp_y2[ tmp ], 1.3*strheight('A'), maxiter=6000, stepsize = 0.05) #, min=0 )
		}

	}



	# max(strwidth(c("asa", "a"))
	# move_text_right <- max(strwidth(outlier_df[,"label_name"]))

	# plotting the outlier labels :)  (I wish there was a non-loop wise way for doing this)
	for(i in seq_len(dim(boxdata_outlier_df)[1]))
	{
		# ss <- (outlier_df[,"x"]  %in% boxdata_outlier_df[i,]$group) & (outlier_df[,"y"] %in% boxdata_outlier_df[i,]$y)

		# if(jitter_if_duplicate) {
			# ss <- (outlier_df[,"x"]  %in% boxdata_outlier_df[i,]$group) & closest.number(outlier_df[,"y"]  boxdata_outlier_df[i,]$y)
		# } else {
		ss <- (outlier_df[,"x"]  %in% boxdata_outlier_df[i,]$group) & (outlier_df[,"y"] %in% boxdata_outlier_df[i,]$y)
		# }

		current_label <- outlier_df[ss,"label_name"]
		temp_x <- boxdata_outlier_df[i,"x"]
		temp_y <- boxdata_outlier_df[i,"y"]
		# cbind(boxdata_outlier_df,		temp_y2)
		# outlier_df



		if(spread_text) {
			temp_y_new <- temp_y2[i] # not ss
			move_text_right <- strwidth(current_label)
			text( temp_x+move_text_right, temp_y_new, offset = 0.1, current_label, cex=0.7, col = label.col)
			# strwidth
			segments( temp_x+(move_text_right/6), temp_y, temp_x+(move_text_right*.47), temp_y_new )
		} else {
			text(temp_x, temp_y, current_label, pos = 4, offset = 0.1, col = label.col, cex=0.7)
		}
	}

	# outputing some of the information we collected
	list(boxdata = boxdata, boxdata_outlier_df = boxdata_outlier_df, outlier_df=outlier_df)
}
