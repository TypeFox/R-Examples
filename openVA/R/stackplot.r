#' plot grouped CSMF from a "insilico" object
#' 
#' Produce bar plot of the CSMFs for a fitted \code{"insilico"} object in broader groups.
#' 
#' 
#' @param x one or a list of fitted object from \code{codeVA} function
#' @param grouping C by 2 matrix of grouping rule. If set to NULL, make it default.
#' @param type type of the plot to make
#' @param order.group list of grouped categories. If set to NULL, make it default.
#' @param order.sub Specification of the order of sub-populations to plot 
#' @param err indicator of inclusion of error bars
#' @param CI confidence interval for error bars.
#' @param sample.size.print Logical indicator for printing also the sample size for each sub-population labels.
#' @param xlab Labels for the causes.
#' @param ylab Labels for the CSMF values.
#' @param ylim Range of y-axis.
#' @param title Title of the plot.
#' @param horiz Logical indicator indicating if the bars are plotted
#' horizontally.
#' @param angle Angle of rotation for the texts on x axis when \code{horiz} is
#' set to FALSE
#' @param border The color for the border of the bars.
#' @param err_width Size of the error bars.
#' @param err_size Thickness of the error bar lines.
#' @param point_size Size of the points.
#' @param bw Logical indicator for setting the theme of the plots to be black
#' and white.
#' @param \dots Not used.
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @importFrom InSilicoVA stackplot
#' @importFrom ggplot2 ggplot aes geom_bar position_dodge geom_errorbar ggtitle coord_flip theme_bw theme element_text
#' @examples
#' \donttest{
#' data(RandomVA3)
#' test <- RandomVA3[1:200, ]
#' train <- RandomVA3[201:400, ]
#' fit1 <- codeVA(data = test, data.type = "customize", model = "InSilicoVA",
#'                     data.train = train, causes.train = "cause",
#'                     Nsim=1000, auto.length = FALSE)
#'
#' fit2 <- codeVA(data = test, data.type = "customize", model = "InterVA",
#'                data.train = train, causes.train = "cause",
#'                version = "4.02", HIV = "h", Malaria = "l")
#'
#' fit3 <- codeVA(data = test, data.type = "customize", model = "Tariff",
#'                data.train = train, causes.train = "cause", 
#'                nboot.sig = 100)
#'
#' data(SampleCategory3)
#' stackplotVA(fit1, grouping = SampleCategory3, type ="dodge", 
#'             ylim = c(0, 1), title = "InSilicoVA")
#' stackplotVA(fit2, grouping = SampleCategory3, type = "dodge", 
#'             ylim = c(0, 1), title = "InterVA4.02")
#' stackplotVA(fit3, grouping = SampleCategory3, type = "dodge", 
#'             ylim = c(0, 1), title = "Tariff")
#' }
#' @export stackplotVA
stackplotVA <- function(x, grouping = NULL,
    type = c("stack", "dodge")[1], 
    order.group = NULL, order.sub = NULL,
    err = TRUE, CI = 0.95, sample.size.print = FALSE,
	xlab = "Group", ylab = "CSMF", ylim = NULL,
	title = "CSMF by broader cause categories", 
	horiz = FALSE, angle = 60,  
	err_width = .4, err_size = .6, point_size = 2, 
	border = "black", bw = FALSE, ...){


	if(is.null(grouping)){
		data("SampleCategory", envir = environment())
		SampleCategory <- get("SampleCategory", envir  = environment())
		grouping <- SampleCategory
		order.group <- c("TB/AIDS", 
						"Communicable",
						"NCD",
						"External",
						"Maternal",
						"causes specific to infancy") 
	}
	if(is.null(order.group)){
		order.group <- unique(grouping[, 2])
	}
 

	csmf <- NULL
	
	if(type == "stack"){
		err <- FALSE
	}

	if(class(x) == "list"){
		counts <- rep(0, length(x))
		for(i in 1:length(x)){
			if(class(x[[i]]) == "insilico"){
				if(is.null(x[[i]]$subpop)){
					# for InSilicoVA, since the error bar needs to be calculated
					# based on the grouping, need the raw CSMF here
					# instead of summarized version as others
					csmf[[i]] <- x[[i]]$csmf
					counts[i] <- length(x[[i]]$id)
				}else{
					stop("Sub-population specification exists in InSilicoVA fit, please rerun with only the InSilicoVA object\n")
				}
			}else{
				csmf[[i]] <- getCSMF(x[[i]], CI = CI)	
				if(class(x[[i]]) == "interVA"){
					counts[i] <- length(x[[i]]$VA) 
				}else if(class(x[[i]]) == "tariff"){
					counts[i] <- dim(x[[i]]$causes.test)[1]
				}
			}
			if(!is.null(names(x)[i])){
				names(csmf)[i] <- names(x)[i]
			}else{
				names(csmf)[i] <- paste("Input", i)
			}
		}
		barNames <- names(x)
		if(is.null(barNames)){
			barNames <- paste("model", 1:length(x))
		}
	}else if(class(x) == "insilico"){
			return(stackplot(x, grouping = grouping,
					    type = type, 
					    order.group = order.group, order.sub = order.sub,
					    err = err, CI = CI, sample.size.print = sample.size.print,
						xlab = xlab, ylab = ylab, ylim = ylim, 
						title = title, 
						horiz = horiz, angle = angle,  
						err_width = err_width, err_size = err_size, 
						point_size = point_size, 
						border = border, bw = bw, ...))
	}else{
		csmf[[1]] <- getCSMF(x, CI = CI)
		if(class(x) == "interVA"){
			counts <- length(x$VA) 
		}else if(class(x) == "tariff"){
			counts <- dim(x$causes.test)[1]
		}
		barNames <- ""
		err <- FALSE
	}


    # to fix dodge and stack gives reverse ordering
    if(type == "dodge"){
      order.group <- rev(order.group)
    }
	
	csmf.group <- NULL
	group <- NULL
	label <- NULL
	low <- (1 - CI) / 2
	high <- 1 - low

 	for(index in 1:length(csmf)){
		
 		# now calculate grouped CSMF
		if(class(x[[index]]) == "insilico"){
			csmf.group.temp <- matrix(0, dim(csmf[[index]])[1], length(order.group))
			colnames(csmf.group.temp) <- order.group
			for(i in 1:length(order.group)){
				which.names <- grouping[which(grouping[, 2] == order.group[i]),1]
				which <- which(colnames(csmf[[index]]) %in% which.names)
				if(length(which) > 1){
					csmf.tmp <- apply(csmf[[index]][, which], 1, sum)
				}else if(length(which) == 1){
					csmf.tmp <- csmf[[index]][, which]
				}else{
					csmf.tmp <- 0
				}
				csmf.group.temp[, i] <- csmf.tmp
			}	

			csmf.mean <- apply(csmf.group.temp, 2, mean)
			csmf.lower <- apply(csmf.group.temp, 2,function(x){quantile(x, low)})
			csmf.upper <- apply(csmf.group.temp, 2,function(x){quantile(x, high)})

		}else{
			csmf.mean <- csmf.lower <- csmf.upper <- rep(0, length(order.group))
			# only grouping
			for(i in 1:length(order.group)){
				which.names <- grouping[which(grouping[, 2] == order.group[i]), 1]
				which <- which(names(csmf[[index]]) %in% which.names)
				if(length(which) > 0){
					csmf.tmp <- sum(csmf[[index]][which])
				}else{
					csmf.tmp <- 0
				}
				csmf.mean[i] <- csmf.lower[i] <- csmf.upper[i] <- csmf.tmp
			}	
		}

		csmf.group <- rbind(csmf.group, cbind(csmf.mean, csmf.lower, csmf.upper))
		group <- c(group, order.group)
		if(sample.size.print){
			label <- c(label, rep(paste(barNames[index], "\n", "n = ", counts[index], sep = ""), length(order.group)))
		}else{
			label <- c(label, rep(barNames[index], length(order.group)))
		}
	}

	rownames(csmf.group) <- NULL
	toplot <- data.frame(csmf.group)
	subpop <- Causes <- NULL
	if(type == "stack" ){
		toplot$Causes <- factor(group, levels = rev(order.group))
	}else if(type == "dodge"){
		toplot$Causes <- factor(group, levels = (order.group))
	}
	toplot$subpop <- label

	# initialize ggplot, force order of bars
	if(horiz){
		g <- ggplot(toplot, aes(x=reorder(subpop, seq(1:length(subpop))), y=csmf.mean, fill=Causes, order = as.numeric(Causes)))	
	}else{
		g <- ggplot(toplot, aes(x=reorder(subpop, seq(length(subpop):1)), y=csmf.mean, fill=Causes, order = -as.numeric(Causes)))
	}
	if(type == "stack"){
		g <- g + geom_bar(stat='identity', color=border, size = .3)
	}else if(type == "dodge"){
		g <- g + geom_bar(stat='identity', color=border, size = .3, position=position_dodge(0.9))
	}
	if(err && type == "dodge"){
		g <- g + geom_errorbar(aes(ymin = csmf.lower, ymax = csmf.upper), size = err_size, width = err_width,  position = position_dodge(.9))
	}
	g <- g + xlab(xlab) + ylab(ylab) 
	if(!is.null(ylim)){
		g <- g + ylim(ylim)
	}
	g <- g + ggtitle(title)
	if(horiz) g <- g + coord_flip()
	if(bw) g <- g + theme_bw()
	hjust <- NULL
	if(angle != 0) hjust = 1
	if(!horiz) g <- g + theme(axis.text.x = element_text(angle = angle, hjust = hjust))	
	return(g)

}