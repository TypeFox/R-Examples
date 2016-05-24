#' plot grouped CSMF from a "insilico" object
#' 
#' Produce bar plot of the CSMFs for a fitted \code{"insilico"} object in broader groups.
#' 
#' 
#' @param x fitted \code{"insilico"} object
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
#' @seealso \code{\link{insilico}}, \code{\link{summary.insilico}}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark Probabilistic cause-of-death
#' assignment using verbal autopsies, \emph{arXiv preprint arXiv:1411.3042}
#' \url{http://arxiv.org/abs/1411.3042} (2014)
#' @keywords InSilicoVA
#' @examples
#' \donttest{
#'   data(RandomVA1) 
#'   ##
#'   ## Scenario 1: without sub-population specification
#'   ##
#'   fit1<- insilico(RandomVA1, subpop = NULL,  
#'                 Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#'                 auto.length = FALSE)
#'   # stack bar plot for grouped causes
#'   # the default grouping could be seen from
#'   data(SampleCategory)
#'   stackplot(fit1, type = "dodge", xlab = "")
#'   
#'   ##
#'   ## Scenario 2: with sub-population specification
#'   ##
#'   data(RandomVA2)
#'   fit2<- insilico(RandomVA2, subpop = list("sex"),  
#'                 Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#'                 auto.length = FALSE)
#'   stackplot(fit2, type = "stack", angle = 0)
#'   stackplot(fit2, type = "dodge", angle = 0)
#'   # Change the default grouping by separating TB from HIV
#'   data(SampleCategory)
#'   SampleCategory[c(3, 9), ]
#'   SampleCategory[3, 2] <- "HIV/AIDS"
#'   SampleCategory[9, 2] <- "TB"
#'   stackplot(fit2, type = "stack", grouping = SampleCategory, 
#'             sample.size.print = TRUE, angle = 0)
#'   stackplot(fit2, type = "dodge", grouping = SampleCategory,
#'             sample.size.print = TRUE, angle = 0)
#'   
#'   # change the order of display for sub-population and cause groups
#'   groups <- c("HIV/AIDS", "TB", "Communicable", "NCD", "External",
#'               "Maternal", "causes specific to infancy") 
#'   subpops <- c("Women", "Men")
#'   stackplot(fit2, type = "stack", grouping = SampleCategory, 
#'             order.group = groups, order.sub = subpops, 
#'             sample.size.print = TRUE, angle = 0)	
#' } 
#' @export stackplot
stackplot <- function(x, grouping = NULL,
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
    # to fix dodge and stack gives reverse ordering
    if(type == "dodge"){
      order.group <- rev(order.group)
    }

	csmf <- NULL
	
	if(type == "stack"){
		err <- FALSE
	}

	if(class(x) == "list"){
		counts <- rep(0, length(x))
		for(i in 1:length(x)){
			csmf[[i]] <- x[[i]]$csmf
			if(!is.null(names(x)[i])){
				names(csmf)[i] <- names(x)[i]
			}else{
				names(csmf)[i] <- paste("Input", i)
			}
			counts[i] <- length(x[[i]]$id)
		}
	}else if(class(x) == "insilico"){
		if(is.null(x$subpop)){
				if(sample.size.print){
					title <- paste(title, "\n", "n = ", length(x$id), sep = "")
				}
				csmf[[1]] <- x$csmf
				names(csmf)[1] <- " "
				# avoid adding one more counts at xlab
				sample.size.print <- FALSE
			}else{
				if(is.null(order.sub)){
					order.sub <- names(x$csmf)
				}
				sx <- summary(x)
				counts <- rep(0, length(order.sub))
				for(i in 1:length(order.sub)){
					csmf[[i]] <- x$csmf[[which(names(x$csmf) == order.sub[i])]]
					names(csmf)[i] <- order.sub[i]
					counts[i] <- sx$subpop_counts[which(names(sx$subpop_counts) == order.sub[i])]
				}	
			}
	}

	csmf.group <- NULL
	group <- NULL
	label <- NULL
	low <- (1 - CI) / 2
	high <- 1 - low

 	for(index in 1:length(csmf)){
		csmf.group.temp <- matrix(0, dim(csmf[[index]])[1], length(order.group))
		colnames(csmf.group.temp) <- order.group

		for(i in 1:length(order.group)){
			which.names <- grouping[which(grouping[, 2] == order.group[i]), 1]
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

		csmf.group <- rbind(csmf.group, cbind(csmf.mean, csmf.lower, csmf.upper))
		group <- c(group, order.group)
		if(sample.size.print){
			label <- c(label, rep(paste(names(csmf)[index], "\n", "n = ", counts[index], sep = ""), length(order.group)))
		}else{
			label <- c(label, rep(names(csmf)[index], length(order.group)))
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