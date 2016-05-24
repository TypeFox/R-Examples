#' plot aggregated COD distribution
#' 
#' Produce a bar plot of the aggregated COD distribution as approximate CSMFs for a fitted \code{"insilico"} object.
#' 
#' 
#' @param x object from \code{get.indiv} function. 
#' @param type An indicator of the type of chart to plot. "errorbar" for line
#' plots of only the error bars on single population; "bar" for bar chart with
#' error bars on single population.
#' @param top The number of top causes to plot. If multiple sub-populations are
#' to be plotted, it will plot the union of the top causes in all
#' sub-populations.
#' @param causelist The list of causes to plot. It could be a numeric vector
#' indicating the position of the causes in the InterVA cause list (see
#' \code{\link{causetext}}), or a vector of character string of the cause
#' names. The argument supports partial matching of the cause names. e.g.,
#' "HIV/AIDS related death" could be abbreviated into "HIV"; "Other and
#' unspecified infect dis" could be abbreviated into "Other and unspecified
#' infect".
#' @param which.plot Specification of which group to plot if there are
#' multiple.
#' @param xlab Labels for the causes.
#' @param ylab Labels for the CSMF values.
#' @param title Title of the plot.
#' @param horiz Logical indicator indicating if the bars are plotted
#' horizontally.
#' @param angle Angle of rotation for the texts on x axis when \code{horiz} is
#' set to FALSE
#' @param fill The color to fill the bars when \code{type} is set to "bar".
#' @param border The color to color the borders of bars when \code{type} is set
#' to "bar".
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
#' 
#' \dontrun{
#' # Toy example with 1000 VA deaths
#' data(RandomVA1)
#' fit1<- insilico(RandomVA1, subpop = NULL,  
#'               Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#'               auto.length = FALSE)
#' summary(fit1, id = "d199")
#' 
#' # update credible interval for individual probabilities to 90%
#' indiv.new <- get.indiv(fit1, CI = 0.9)
#' fit1$indiv.prob.lower <- indiv.new$lower
#' fit1$indiv.prob.upper <- indiv.new$upper
#' fit1$indiv.CI <- 0.9
#' summary(fit1, id = "d199")
#' 
#' 
#' # get empirical aggregated COD distribution 
#' agg.csmf <- get.indiv(data = RandomVA2, fit1, CI = 0.95, 
#'                       is.aggregate = TRUE, by = NULL)
#' head(agg.csmf)
#' 

#' # aggregate individual COD distribution by sex and age
#' # note the model was fitted assuming the same CSMF for all deaths
#' # this aggregation provides an approximate CSMF for each sub-groups
#' agg.by.sex.age <- get.indiv(data = RandomVA2, fit1, CI = 0.95, 
#'                         is.aggregate = TRUE, by = list("sex", "age"))
#' head(agg.by.sex.age$mean)
#' 
#' # plot of aggregated individual COD distribution
#' # 0. plot for all data
#' indivplot(agg.csmf, top = 10)
#' # 1. plot for specific one group
#' indivplot(agg.by.sex.age, which.plot = "Men 60-", top = 10)
#' # 2. comparing multiple groups
#' indivplot(agg.by.sex.age, which.plot = list("Men 60+", "Men 60-"), 
#'                           top = 5)
#' # 3. comparing multiple groups on selected causes
#' indivplot(agg.by.sex.age, which.plot = list("Men 60-", "Women 60-"), 
#'                           top = 0, causelist = c(
#'                             "HIV/AIDS related death", 
#'                             "Pulmonary tuberculosis", 
#'                             "Other and unspecified infect dis", 
#'                             "Other and unspecified NCD"))
#' } 
#' @export plot.insilico
indivplot <- function(x, type = c("errorbar", "bar")[1], 
	top = 10, causelist = NULL, which.plot = NULL, 
	xlab = "Causes", ylab = "COD distribution", 
	title = "COD distributions for the top causes", 
	horiz = TRUE, angle = 60, fill = "lightblue", 
	err_width = .4, err_size = .6, point_size = 2, 
	border = "black", bw = FALSE, ...){
	

	
	if(class(x) == "list" && is.null(which.plot)){
		compare <- TRUE
		which.plot <- colnames(x$mean)
	}else if(class(which.plot) == "list" && length(which.plot) > 1){
		compare <- TRUE
	}else{
		compare <- FALSE
	}

	if(class(x) != "list"){
		compare <- FALSE
		which.plot <- "All"
	}
	# to please R CMD Check with ggplot
	Group <- Mean <- Lower <- Upper <- Causes <- NULL
	
	## reshape data into: 
	## [mean, median, lower, upper, group, cause]
	if(class(x) == "list"){
		dist <- cbind(as.vector(x$mean), 
					  as.vector(x$median),
					  as.vector(x$lower),
					  as.vector(x$upper))
		group <- rep(colnames(x$mean), each = dim(x$mean)[1]) 
		dist <- data.frame(dist, check.names = FALSE)
		colnames(dist) <- c("Mean", "Median", "Lower", "Upper")
		dist$Group = group
		dist$Causes <- rep(rownames(x$mean), dim(x$mean)[2]) 
		group.names <- unique(group)
		if(!is.null(which.plot)){
			group.names <- intersect(group.names, which.plot)
		}		
	}else{
		dist <- data.frame(x)
		group.names <- "All"
		which.toplot <- "All"
		dist$Group <- "All"
		dist$Causes <- rownames(x)
	}

	# get which causes to plot
	# if causelist is provided then fine, 
	# otherwise construct from top argument
	if(top > 0 && is.null(causelist)){
		for(i in 1:length(group.names)){
			sub <- which(dist$Group == group.names[i])
			topcause <- order(dist[sub, "Mean"], decreasing=TRUE)[1:round(top)]
			causelist <- c(causelist, 
					setdiff(dist[sub[topcause], "Causes"], causelist))
		}
	}
	
	csmf.toplot <- NULL
	for(i in 1:length(which.plot)){
		csmf.toplot.tmp <- subset(dist, (Group == which.plot[i]))
		csmf.toplot.tmp <- csmf.toplot.tmp[match(causelist, csmf.toplot.tmp$Causes), ]
		csmf.toplot <- rbind(csmf.toplot, csmf.toplot.tmp)
	}
	csmf.toplot <- data.frame(csmf.toplot)
	
	# making bar plot for comparison
	if(compare){
		# initialize ggplot, force order of bars
		if(horiz){
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(length(Causes), 1)),
									 y=Mean, ymax=max(Upper)*1.05,
									 fill = Group))
		}else{
			g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(1:length(Causes))),
									 y=Mean, ymax=max(Upper)*1.05, 
									 fill = Group))
		}
		g <- g + geom_point( aes(color=Group), position=position_dodge(0.5), size = point_size) 
		g <- g + geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Group), size = err_size, width = err_width, position = position_dodge(.5))
		g <- g + xlab(xlab) + ylab(ylab) 
		g <- g + ggtitle(title)
		g <- g + scale_y_continuous() 
		if(horiz) g <- g + coord_flip()
		if(bw) g <- g + theme_bw()
		if(!horiz) g <- g + theme(axis.text.x = element_text(angle = angle, hjust = 1))
		return(g)
	}else{		
		# making bar plot
		# require: number of top causes or which causes to plot;
		#		   color list
		if(type == "bar"){
			# initialize ggplot, force order of bars
			if(horiz){
				g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(length(Causes), 1)),
										 y=Mean))	
			}else{
				g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(1:length(Causes))),
										 y=Mean))
			}
			g <- g + geom_bar( stat="identity", 
					 colour=border, fill=fill,  size = .3)
			## todo		
			g <- g + geom_errorbar(aes(ymin = Lower, ymax = Upper), size = err_size, width = err_width,  position = position_dodge(.9))
			g <- g + xlab(xlab) + ylab(ylab) 
			g <- g + ggtitle(title)
			if(horiz) g <- g + coord_flip()
			if(bw) g <- g + theme_bw()
			if(!horiz) g <- g + theme(axis.text.x = element_text(angle = angle, hjust = 1))	
			return(g)
		}

		# making error bar plot
		# require: number of top causes or which causes to plot;
		#		   color list
		if(type == "errorbar"){
			# initialize ggplot, force order of bars
			if(horiz){
				g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(length(Causes), 1)),
										 y=Mean))	
			}else{
				g <- ggplot(csmf.toplot, aes(x=reorder(Causes, seq(1:length(Causes))),
										 y=Mean))
			}
			g <- g + geom_point( stat="identity", 
					 colour=border, fill=fill,  size = point_size)
			g <- g + geom_errorbar(aes(ymin = Lower, ymax = Upper), size = err_size, width = err_width, position = position_dodge(.9))
			g <- g + xlab(xlab) + ylab(ylab) 
			g <- g + ggtitle(title)
			if(horiz) g <- g + coord_flip()
			if(bw) g <- g + theme_bw()
			if(!horiz) g <- g + theme(axis.text.x = element_text(angle = angle, hjust = 1))	
			return(g)
		}
	}   

}
