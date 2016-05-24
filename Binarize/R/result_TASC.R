#############################################Class TrinarizationResult############################################
#This is the base class of all results of the trinarization functions. It provides the basic methods show, print and 
#a method called plotBinarization. It also checks all created object for validity.
setClass(
	Class = "TrinarizationResult", 
	representation = representation(
		originalMeasurements = "numeric",
		trinarizedMeasurements = "integer", 
		threshold1 = "numeric",
		threshold2 = "numeric",
		method = "character",
		p.value = "numeric"
	),
	validity = function(object){
		#extract object slots
		omeasure <- object@originalMeasurements
		bmeasure <- object@trinarizedMeasurements
		thresh1 <- object@threshold1
		thresh2 <- object@threshold2
		meth <- object@method
		p.value <- object@p.value
		
		
		#initialize the basic strings
		valid_methods <- c(
			"TASC A",
			"TASC A (min)",
			"TASC B",
			"TASC B (min)"
		)
		for(i in seq(1, length(valid_methods))){
			valid_methods_string <- ifelse(i==1, sprintf("\"%s\"", valid_methods[i]), sprintf("%s, \"%s\"", valid_methods_string, valid_methods[i]))
		}
		
		#initialize the critical error messages
		critical_invalid_strings <- c(
			"'originalMeasurements' isn't set!",
			"'trinarizedMeasurements' isn't set!",
			"'threshold1' isn't set!",
			"'threshold2' isn't set!",
			"'method' isn't set!",
			"'p.value' isn't set!"
		)
		#check object for critical errors
		critical_invalid <- c(
			!length(omeasure),
			!length(bmeasure),
			!length(thresh1),
			!length(thresh2),
			!length(meth),
			!length(p.value)
		)
		#if critical error occured return the corresponding error messages
		if (sum(as.integer(critical_invalid))){
			return(critical_invalid_strings[which(critical_invalid)])
		}
		
		#initialize the weak error messages
		weak_invalid_strings <- c(
			"Only zeros, ones and twos are valid values for 'trinarizedMeasurements'.",
			sprintf("'method' must be element of {%s}, but it is \"%s\".", valid_methods_string, as.character(meth)),
			"Length of original and trinarized Measurements must be the same.",
			sprintf("'threshold1' and 'threshold2' must be within the borders of the original values, which is the interval [%f, %f], but they are %f and %f.", min(omeasure), max(omeasure), thresh1, thresh2),
			"'p.value' must be in range [0,1]."
		)
		#check object for weak errors
		weak_invalid <- c(
			length(which(bmeasure > 2)) || length(which(bmeasure < 0)),
			length(which(valid_methods == meth)) < 1,
			length(bmeasure) != length(omeasure),
			thresh1 < min(omeasure) || thresh1 > max(omeasure) || thresh2 < min(omeasure) || thresh2 > max(omeasure),
			(!is.na(p.value) && (p.value < 0 || p.value > 1))
		)
		#if weak error occured return the corresponding error messages
		if (sum(as.integer(weak_invalid))){
			return(weak_invalid_strings[which(weak_invalid)])
		}
		
		#object is valid
		return(TRUE)
	}
)

#This method prints the last three slots out to console (trinarizedMeasurements is limited to 10 values). It is called
#when creating an object without an assignment or by only typing the name of a TrinarizationResult-object at console.
setMethod(
	f = "show",
	signature = "TrinarizationResult",
	definition = function(object){
		cat("Method: ", object@method, "\n",sep="")
		if(length(object@trinarizedMeasurements) <= 10){
			cat("\nTrinarized vector: [ ", paste(object@trinarizedMeasurements, collapse=" "),
				" ]\n",sep="")
		}else{
			cat("\nTrinarized vector: [ ",paste(object@trinarizedMeasurements[1:10], collapse=" "),
				" ...]\n",sep="")
		}
		cat("\nThreshold1: ", object@threshold1, "\n", sep="")
		cat("\nThreshold2: ", object@threshold2, "\n", sep="")
		if(!is.na(object@p.value)){
			cat("\np value: ", object@p.value, "\n", sep="")
		}
	}
)

#This method prints the last three slots out to console
setMethod(
	f = "print",
	signature = "TrinarizationResult",
	definition = function(x){
		cat("Method: ", x@method, "\n", sep="")
		cat("\nThreshold1: ", x@threshold1, "\n", sep="")
		cat("\nThreshold2: ", x@threshold2, "\n", sep="")
		cat("\nTrinarized vector: [ ", paste(x@trinarizedMeasurements, collapse=" "),
			" ]\n", sep="")
		if (!is.na(x@p.value)){
			cat("\np value: ",x@p.value,"\n", sep="")
		}
	}
)

setGeneric("plot", useAsDefault = plot)

#This Method plots the computed binarization in a one- or two-dimensional way.
setMethod(
	f = "plot",
	signature = c("TrinarizationResult"),
	definition = function(x, twoDimensional=FALSE, showLegend=TRUE, showThreshold=TRUE, ...)
	{
		if (twoDimensional){
			plot(1:length(x@trinarizedMeasurements), x, showLegend=showLegend, showThreshold=showThreshold, ...)
		}else{
			#extract the base values of x
			vect_length <- length(x@originalMeasurements)
			min_val <- min(x@originalMeasurements) #floor(min(c(x@originalMeasurements,0)))
			max_val <- max(x@originalMeasurements) #ceiling(max(c(x@originalMeasurements,0)))

			#get the ... argument into a list
			args <- list(...)

			#check for several standard graphic parameters and if they aren't set, set them to default values
			if(is.null(args$ylab)){
				args$ylab <- ""
			}
			if(is.null(args$xlab)){
				args$xlab <- ""
			}
			if(is.null(args$lty)){
				args$lty <- 2
			}
			if(is.null(args$pch)){
				args$pch <- x@trinarizedMeasurements
			}else if(length(args$pch) == 3){
				args$pch <- args$pch[x@trinarizedMeasurements+1]
			}

			col <- args$col

			if(is.null(col)){
				col <- c("red","green","blue","black")
			}
			if(length(col) < 3){
				col <- rep(col, 3)[1:3]
			}
			if(length(col) == 3){
				col <- c(col, "black")
			}
			if(length(col) == 4){
				args$col <- col[x@trinarizedMeasurements+1]
			}

			if(is.null(args$type)){
				args$type <- "p"
			}

			if(is.null(args$yaxt)){
				args$yaxt="n"
			}

			#plotting the axes shouldn't be controlled by standard plot function
			#this method does it later
			#args$axes <- FALSE

			#check for the limit standard graphic parameters and if they aren't set, set them to default values
			if(is.null(args$xlim)){
				args$xlim <- c(min_val,max_val)
			}
			if(is.null(args$ylim)){
				args$ylim <- c(-0.1,0.1)
			}
			#set the point coordinates
			args$x <- x@originalMeasurements
			args$y <- rep(0,vect_length)

			#plot them
			do.call("plot", args)

			#plot the threshold as line
			if(as.logical(showThreshold)){
				par(new=TRUE)
				largs <- list(...)

				if(is.null(largs$lty)){
					largs$lty <- 2
				}

				largs$col <- col[4]

				do.call("abline", c(largs,v=x@threshold1))
				do.call("abline", c(largs,v=x@threshold2))
			}

			#if axes isn't set or TRUE plot the x-axis
			#if (is.null(list(...)$axes) || as.logical(list(...)$axes) || list(...)$yaxt != "n"){
			#    if (is.null(args$lwd)){
			#        lwd <- 1
			#    }else{
			#        lwd <- args$lwd
			#    }
			#    at <- round(seq(min_val,max_val,by=(max_val-min_val)/5),1)
			#    axis(1, at=at, lwd=lwd, pos=-0.01)
			#    #axis(1, at=at, lwd=lwd, pos=-0.05)#c(min_val,-10))
			#    #axis(1, at=at, lwd=lwd, pos=-0.1)#c(min_val,-10))
			#}

			if(as.logical(showLegend)){
				if(is.null(args$lwd)){
					lwd <- 1
				}else{
					lwd <- args$lwd
				}
				if(as.logical(showThreshold)){
					if(is.null(args$pch)){
						pch <- c(0,1,2,NA)
					}else if(length(args$pch) > 3){
						pch <- c(15, 16, 17, NA)
					}else{
						pch <- c(unique(args$pch), NA)
					}
					names <- c("zero", "one", "two", "threshold")
					lty <- c(NA, NA, NA, args$lty[1])
				}else{
					if(is.null(args$pch)){
						pch <- c(0,1,2)
					}else if(length(args$pch) > 3){
						pch <- c(15, 16, 17)
					}else{
						pch <- unique(args$pch)
					}
					names <- c("zero", "one", "two")
					lty <- c(NA, NA, NA)
					#if(is.null(args$col)){
					#    col <- "black"
					#}else if (length(args$col) < 3){
					#    col <- args$col
					#}else{
					#    col <- args$col[1:2]
					#}
				}

				legend("topleft", names, pch=pch,
					lty=lty, inset=c(0.05, 0.05), bty="n", cex=0.8, lwd=lwd, col=col)
			}
		}
	}
)

setMethod(
	f = "plot",
	signature = c("numeric","TrinarizationResult"),
	definition = function(x, y, showLegend=TRUE, showThreshold=TRUE, ...)
	{
		#extract the base values of y
		vect_length <- length(y@originalMeasurements)
		min_val <- min(y@originalMeasurements) 
		max_val <- max(y@originalMeasurements) 
		
		#get the ... argument into a list
		args <- list(...)
		
		#check for several standard graphic parameters and if they aren't set, set them to default values
		if(is.null(args$ylab)){
			args$ylab <- ""
		}
		if(is.null(args$xlab)){
			args$xlab <- ""
		}
		if(is.null(args$lty))
		   args$lty <- 2
		if(is.null(args$cex.axis))
		   args$cex.axis <- par("cex.axis")
		if(is.null(args$cex.lab))
		   args$cex.lab <- par("cex.lab")   
			
		if(is.null(args$pch))
		{
			args$pch <- y@trinarizedMeasurements
		}
		else 
		if(length(args$pch) == 3)
		{
			args$pch <- args$pch[y@trinarizedMeasurements+1]
		}
		
		col <- args$col
		
		if(is.null(col)){
			col <- c("red","green","blue","black")
		}
		if(length(col) < 3){
			col <- rep(col, 3)[1:3]
		}
		if(length(col) == 3){
			col <- c(col, "black")
		}
		if(length(col) == 4){
			args$col <- col[y@trinarizedMeasurements+1]
		}
		
		if(is.null(args$type))
			args$type <- "p"

		#plotting the axes shouldn't be controlled by standard plot function
		#this method does it later
		#args$axes <- FALSE
		args$xaxt="n"
		
		#maxx is the minimal value >= vect_length and dividable by vect_length DIV 5
		#(for example: if vect_length = 11 => maxx is 12 and if vect_length = 19 => maxx is 21)
		#maxx <- ifelse(vect_length%%5==0, vect_length, vect_length%/%5*6)
		#while(maxx < vect_length)
		#    maxx <- maxx + vect_length%/%5
		
		#check for the limit standard graphic arguments. if not set set them to default values
		#if (is.null(args$xlim))
		#    args$xlim <- c(0, maxx)
		#if (is.null(args$ylim))
		#    args$ylim <- c(min_val, max_val)
		
		#plot the binarization
		args$x <- x
		#seq(along = x@originalMeasurements)
		args$y <- y@originalMeasurements
		print(args)
		do.call("plot", args)
		
		
		#plot the threshold as line
		if(as.logical(showThreshold)){
			par(new=TRUE)
			largs <- list(...)
			
			if (is.null(largs$lty))
			  largs$lty <- 2
			
			largs$col <- col[4]
				
			do.call("abline", c(largs,h=y@threshold1))
			do.call("abline", c(largs,h=y@threshold2))
		}
		
		#if axes isn't set or TRUE plot the x and y axis according to maxx, min_val, max_val
		if(is.null(list(...)$axes) || as.logical(list(...)$axes) || list(...)$xaxt != "n")
		{
		  if(is.null(args$lwd))
		  {
				lwd <- 1
		  }
		  else
		  {
			  lwd <- args$lwd
		  }
		  axis(1, at=x, lwd=lwd, cex.axis=args$cex.axis, cex.lab=args$cex.lab)
		}

	  if(as.logical(showLegend))
	  {
		  if(is.null(args$lwd))
		  {
			  lwd <- 1
		  }
		  else{
			  lwd <- args$lwd
		  }
		  if(as.logical(showThreshold))
		  {
			  if(is.null(args$pch)){
				  pch <- c(0,1,2,NA)
			  }
			  else
			  if(length(args$pch) > 3)
			  {
				pch <- c(15, 16, 17, NA)
			  }
			  else{
				  pch <- c(unique(args$pch), NA)
			  }
			  names <- c("zero", "one", "two", "threshold")
			  lty <- c(NA, NA, NA, args$lty[1])
		  }
		  else
		  {
			  if(is.null(args$pch))
			  {
				  pch <- c(0,1,2)
			  }
			  else
			  if(length(args$pch) > 3)
			  {
				pch <- c(15, 16, 17)
			  }
			  else{
				  pch <- unique(args$pch)
			  }
			  names <- c("zero", "one", "two")
			  lty <- c(NA, NA, NA)
		  }

		  legend("topleft", names, pch=pch,
			  lty=lty, inset=c(0.05, 0.05), bty="n", cex=0.8, lwd=lwd, col=col)
		}
	}
)


#############################################Class TASCResult##############################################
#This is the result class for the TASC algorithm. It provides an additional method called plotStepFunctions and is
#derived from the BinarizationResult class.
setClass(
	Class = "TASCResult", 
	representation = representation(
		intermediateSteps = "matrix",
		intermediateHeights1 = "matrix",
		intermediateHeights2 = "matrix",
		intermediateStrongestSteps = "matrix"),
	contains = "TrinarizationResult",
	validity = function(object){
		#extract relevant object slots
		isteps <- object@intermediateSteps
		iheights1 <- object@intermediateHeights1
		iheights2 <- object@intermediateHeights2
		istrsteps <- object@intermediateStrongestSteps
		omeasure <- object@originalMeasurements
		
		#initialize the critical error messages
		critical_invalid_strings <- c(
			"'intermediateSteps' isn't set!",
			"'intermediateHeights1' isn't set!",
			"'intermediateHeights2' isn't set!",
			"'intermediateStrongestSteps' isn't set!"
		)
	
		#check object for critical errors
		critical_invalid <- c(
			!length(isteps),
			!length(iheights1),
			!length(iheights2),
			!length(istrsteps)
		)
		#if critical error occured return the corresponding error messages
		if(sum(as.integer(critical_invalid))){
			return(critical_invalid_strings[which(critical_invalid)])
		}
		
		#initialize weak error messages
		weak_invalid_strings <- c(
			"'intermediateSteps' and 'intermediateHeights' must have the same dimensionality.",
			"'intermediateStrongestSteps' must have the same number of rows as 'intermediateSteps'.",
			"The values of 'intermediateSteps' must be in range [0, #Measurements].",
			"The values of 'intermediateStrongestSteps' must be in range [1, #Measurements]."
		)
		
		#check object for weak errors
		weak_invalid <- c(
			as.logical(sum(dim(isteps) != dim(iheights1) || dim(isteps) != dim(iheights2))),
			nrow(istrsteps) != nrow(isteps),
			(sum(isteps < 0) || sum(isteps > length(omeasure))),
			(sum(istrsteps[,1] < 1) || sum(istrsteps[,1] > length(omeasure)) || sum(istrsteps[,2] < 1) || sum(istrsteps[,2] > length(omeasure)))
		)
		#if weak error occured return the corresponding error messages
		if(sum(as.integer(weak_invalid))){
			return(weak_invalid_strings[which(weak_invalid)])
		}
		
		#object is valid
		return(TRUE)
	}
)

#This method plots all the computed optimal step functions with n steps in one diagram. These step functions are formed
#by the two BASC algorithms and are used to determine the optimal jumping point and are also used to calculate the
#P-Value.
setMethod(
	f = "plotStepFunctions",
	signature = "TASCResult",
	definition = function(x, showLegend=TRUE, connected=FALSE, withOriginal=TRUE, ...){
		#check the input BASCResult-Object
		if(ncol(x@intermediateSteps) == 0 || nrow(x@intermediateSteps) == 0)
			stop("intermediateSteps has no values to plot.")
		if(ncol(x@intermediateHeights1) == 0 || nrow(x@intermediateHeights1) == 0)
			stop("intermediateHeights has no values to plot.")
		if(ncol(x@intermediateHeights2) == 0 || nrow(x@intermediateHeights2) == 0)
			stop("intermediateHeights has no values to plot.")
		if(ncol(x@intermediateStrongestSteps) == 0 || nrow(x@intermediateStrongestSteps) == 0)
			stop("intermediateStrongestSteps has no values to plot.")
			
		#get the value-count
		vect_count <- length(x@originalMeasurements)
		
		#steps is a matrix with all the jump indices computed by the C-function concatenated with
		#1:vect_count which is used for plotting the original step-function
		if(as.logical(withOriginal)){
			steps <- matrix(nrow=nrow(x@intermediateSteps)+1, ncol = vect_count, data = rep(0,(nrow(x@intermediateSteps)+1) * vect_count))
			steps[1:(nrow(steps)-1),1:ncol(x@intermediateSteps)] <- x@intermediateSteps
			steps[nrow(steps),] <- seq(along=x@originalMeasurements)
		}else{
			steps <- matrix(nrow=nrow(x@intermediateSteps), ncol = vect_count, data = rep(0,nrow(x@intermediateSteps) * vect_count))
			steps[1:nrow(steps),1:ncol(x@intermediateSteps)] <- x@intermediateSteps
		}
		
		#heights is a matrix with all the jump heights computed by the C-function concatenated with
		#the jump heights of the original step-function
		if(as.logical(withOriginal)){
			heights1 <- matrix(nrow=nrow(x@intermediateSteps)+1, ncol = vect_count, data = rep(0,(nrow(x@intermediateSteps)+1) * vect_count))
			heights1[1:(nrow(heights1)-1),1:ncol(x@intermediateSteps)] <- x@intermediateHeights1
			heights1[nrow(heights1),] <- c(diff(sort(x@originalMeasurements)), 0)

			heights2 <- matrix(nrow=nrow(x@intermediateSteps)+1, ncol = vect_count, data = rep(0,(nrow(x@intermediateSteps)+1) * vect_count))
			heights2[1:(nrow(heights2)-1),1:ncol(x@intermediateSteps)] <- x@intermediateHeights2
			heights2[nrow(heights2),] <- c(diff(sort(x@originalMeasurements)), 0)
		}else{
			heights1 <- matrix(nrow=nrow(x@intermediateSteps), ncol = vect_count, data = rep(0,nrow(x@intermediateSteps) * vect_count))
			heights1[1:nrow(heights1),1:ncol(x@intermediateSteps)] <- x@intermediateHeights1

			heights2 <- matrix(nrow=nrow(x@intermediateSteps), ncol = vect_count, data = rep(0,nrow(x@intermediateSteps) * vect_count))
			heights2[1:nrow(heights1),1:ncol(x@intermediateSteps)] <- x@intermediateHeights2
		}
		
		heights1 <- t(apply(heights1,1,function(x)x/sum(x)))
		heights2 <- t(apply(heights2,1,function(x)x/sum(x)))
		
		#the maximal y-value is calculated. y starts at 1, all the individual jump-heights are added and
		#between every single step-function there's 0.5 free space
		maxy <- nrow(steps) * 0.5 + sum(heights1) + 1
		
		# #maxx is the minimal value >= vect_length and dividable by vect_length DIV 5
		# #(for example: if vect_length = 11 => maxx is 12 and if vect_length = 19 => maxx is 21)
		#maxx <- ifelse((vect_count%%5)==0, vect_count, (vect_count%/%5)*6)
		#while(maxx < vect_count)
		#    maxx <- maxx + (vect_count%/%5)
		
		maxx <- vect_count
			
		#calculate the coordinates of the lines of the step-functions
		lines <- sapply(
			#loop over the rows of steps from last row to first row
			rev(seq(along = steps[,1])),
			function(i, st, he){
				#calculate the base y-value of the current "line"
				#it is calculated like maxy but only for the first i lines.
				cury <- ifelse(i < nrow(he), sum(he[seq(i + 1, nrow(he)),]) + (nrow(he) - i + 1) * 0.5, 0.5)
				#cury <- ifelse(i < nrow(he), sum(he[seq(i + 1, nrow(he)),]), 0.5)
				
				#get the current steps and heights row
				cur_steps <- st[i, st[i,] > 0]
				cur_heights <- he[i,]
				
				#count is the current number of single lines of the current step functions
				#except the last line which is always added directly before the return statement 
				count <- min(vect_count-1, length(cur_steps))
				if(!as.logical(connected))
					lines <- matrix(nrow=2, ncol=2+4*count)
				else
					lines <- matrix(nrow=2, ncol=2+2*count)
				#construct the coordinates of the lines first and last x,y-pair will be added 
				#after the next block
				lines[,seq(2, ncol(lines)-1)] <- matrix(
					sapply(
						seq(1, count), 
						function(j,s,h,base){
							#the NAs are neccessary because vertical lines direct at a step shouldn't be 
							#plotted
							if(!as.logical(connected)){
								result <- matrix(ncol=4, nrow=2, rep(NA,8))
								result[1,c(1,4)] <- rep(s[j], 2)
								result[2, 1] <- ifelse(j==1, base, base + sum(h[1:j-1]))
								result[2, 4] <- base + sum(h[1:j])
							}
							else{
								result <- matrix(ncol=2, nrow=2, rep(NA,4))
								result[1,c(1,2)] <- rep(s[j], 2)
								result[2, 1] <- ifelse(j==1, base, base + sum(h[1:j-1]))
								result[2, 2] <- base + sum(h[1:j])
							}
							return(result)
						}, 
						cur_steps, 
						cur_heights, 
						cury
					)
				)
				#set the first and the last coordinates pair and return all the coordinates for the
				#current step-function
				lines[, 1] <- c(0, cury)
				lines[, ncol(lines)] <- c(vect_count, lines[2, ncol(lines) - 1])               
				return(lines)
			},
			steps,
			heights1
		)
		
		#calculate the coordinates for the lines of the respective strongest steps
		if(as.logical(withOriginal)){
			ncol <- 3 * (nrow(steps) - 1)
			sequence <- rev(seq(1, nrow(steps) - 1))
		}
		else{
			ncol <- 3 * nrow(steps)
			sequence <- rev(seq(1, nrow(steps)))
		}
		strongestLines <- matrix(
			nrow = 2, 
			ncol = ncol,
			data = sapply(
				sequence,
				function(i, st, he, x, l){
					#get the current values
					cur_steps <- st[i, st[i,] > 0]
					cur_heights <- he[i,]
					cur_l <- l[[length(l) - i + 1]]
							
					#get the coordinates from the current values
					result <- matrix(nrow = 2, ncol = 3, data = rep(NA, 6))
					result[1,c(2,3)] <- rep(x@intermediateStrongestSteps[i,1], 2)
					result[2,2] <- cur_l[2, which(cur_l[1,] == x@intermediateStrongestSteps[i,1])[1]]
					result[2,3] <- cur_l[2, which(cur_l[1,] == x@intermediateStrongestSteps[i,1])[2]]
					
					return(result)
				}, 
				steps, 
				heights1, 
				x, 
				lines
			)
		)
		strongestLines2 <- matrix(
			nrow = 2, 
			ncol = ncol,
			data = sapply(
				sequence,
				function(i, st, he, x, l){
					#get the current values
					cur_steps <- st[i, st[i,] > 0]
					cur_heights <- he[i,]
					cur_l <- l[[length(l) - i + 1]]
							
					#get the coordinates from the current values
					result <- matrix(nrow = 2, ncol = 3, data = rep(NA, 6))
					result[1,c(2,3)] <- rep(x@intermediateStrongestSteps[i,2], 2)
					result[2,2] <- cur_l[2, which(cur_l[1,] == x@intermediateStrongestSteps[i,2])[1]]
					result[2,3] <- cur_l[2, which(cur_l[1,] == x@intermediateStrongestSteps[i,2])[2]]
					
					return(result)
				}, 
				steps, 
				heights2, 
				x, 
				lines
			)
		)
		
		strongestLines <- cbind(strongestLines, strongestLines2)

		#insert NA's at the strongestStep positions, because this lines are plotted seperatly
		if(as.logical(connected)){
			for(i in seq(along=lines)){
				l <- lines[[i]]
				matched <- match(strongestLines[2,],l[2,])
				matched <- matched[!is.na(matched)]
				if (length(matched) > 0){
					ind <- max(matched)
					l <- matrix(nrow=2,data=c(l[,1:(ind-1)],NA,NA,l[,-(1:(ind-1))]))
					lines[[i]] <- l
				}
			}
		}
		
		#put the additional arguments in args
		args <- list(...)
		
		#check several standard graphics parameter and set them to default values if they aren't
		#set yet
		if(is.null(args$xlim))
			args$xlim <- c(0, maxx)
		if(is.null(args$ylim))
			args$ylim <- c(0, maxy*1.01)
		if(is.null(args$pch))
			args$pch <- c(1,20)
		if(is.null(args$type))
			args$type <- "o"
		else
			args$type <- args$type[1]   #if args$type is a vector then take first element for standard-lines
										#and the second element for strongest-step lines others will be ignored
		if(is.null(args$cex))
			args$cex <- c(1,1.2)
		if(is.null(args$ylab))
			args$ylab <- ""
		if(is.null(args$xlab))
			args$xlab <- ""
		if(is.null(args$lty))
			args$lty <- 1
		else
			args$lty <- args$lty[1]     #same handling as args$type
		if(!is.null(args$col))
			args$col <- args$col[1]     #same handling as args$type
		#drawing axes will be handled later by this function and not by the standard plot function
		args$axes <- FALSE
		
		#plot the step functions
		lapply(
			lines, 
			function(l){
				args$x <- l[1,]
				args$y <- l[2,]
				
				do.call("plot", args)
				#par(new=TRUE) is neccessary beacause the old lines shouldn't be deleted
				par(new=TRUE)
			}
		)
		
		#setup args for plotting strongest steps        
		args$x <- strongestLines[1,]
		args$y <- strongestLines[2,]
		if(is.null(list(...)$type))             
			args$type <- "l"
		else if(length(list(...)$type) > 1)
			args$type <- list(...)$type[2]
		else
			args$type <- list(...)$type
			
		if(is.null(list(...)$lty))
			args$lty <- 2
		else if(length(list(...)$lty) > 1)
			args$lty <- list(...)$lty[2]
		else
			args$lty <- list(...)$lty
			
		if(!is.null(list(...)$col) & length(list(...)$col) > 1){
			args$col <- list(...)$col[2]
		}
		
		#plot the strongest steps of the step functions
		do.call("plot", args)
		
		if(is.null(list(...)$lwd))
			lwd <- 1
		else
			lwd <- list(...)$lwd
		
		#if axes isn't set or set to TRUE then plot the x-axes
		if(is.null(list(...)$axes) || as.logical(list(...)$axes)){
			axis(1, pos=0, at=seq(0,maxx,by=(vect_count%/%5)), lwd = lwd)
		}
		
		#if showLegend is TRUE plot a legend
		if(as.logical(showLegend)){
			#if lty wasn't set take the default values for the line types else take the first
			#two values (if possible) for the line type
			if(is.null(list(...)$lty))
				lty <- c(1,2) 
			else if(length(list(...)$lty) == 1)
				lty <- list(...)$lty
			else
				lty <- list(...)$lty[c(1,2)]
				
			if(is.null(list(...)$col))
				col <- "black"
			else if(length(list(...)$col) == 1)
				col <- list(...)$col
			else
				col <- list(...)$col[1:2]
			
			legend("topleft", c("steps","strongest steps"), lty=lty, col=col, inset=c(0.05,0), bty="n", cex=0.8, lwd=lwd)
		}
	}
)



