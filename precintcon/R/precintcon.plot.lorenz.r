#' @export
precintcon.plot.lorenz <- function(..., 
  interval        = 1,
	grouped         = FALSE, 
	xlab            = expression(sum(n[i]), i==1),
	ylab            = expression(sum(P[i]), i==1), 
	legend.title    = "Legend",
	legend          = NULL,
	fontsize        = 10, 
	axis.text.color = "black", 
	export          = FALSE, 
	export.name     = "lorenz_plot.png", 
  width           = 8.6, 
  height          = 7.5, 
  units           = "cm",
  args            = NA
) {
	
	l <- list(...)

	##
	# Check all input types
	#
	if (checkType(l)) {

		##
		# Activating break before to show a plot
		#
		if (length(l) > 1 && !grouped)
			par(ask=T)
		
		##
		# Getting parameters name
		#
	  varl <- ifelse(is.na(args), as.character(match.call()[1:length(l)+1]), args)
		
		if (!is.null(legend) && length(varl) != length(legend))
			stop(paste("legend should has length equals to the number of input data. legend parameter length", 
							length(legend), ": number of input data", length(varl)))
		
		else if (!is.null(legend))
			varl <- as.list(legend)
		
		##
		# Should plot be grouped?
		#
		if (grouped) {
			
			l <- lapply(l, FUN=as.precintcon.fd, interval=interval)
			
			data <- data.frame()
			
			for(i in 1:length(l)) {
				
				d <- l[[i]]
				
				d <- rbind(data.frame(
					initial.class = 0, 
					final.class   = 0, 
					midpoint      = 0,
					n             = 0,
					sum.n         = 0,
					P             = 0,
					sum.P         = 0,
					p.sum.n       = 0,
					p.sum.P       = 0), d)
			
				data <- rbind(data, data.frame(dataset=paste(varl[[i]]), x=d$p.sum.n, y=d$p.sum.P, s=sum(d$p.sum.P)))
			}
			
			data <- data[order(data$s), ]

			plot <- ggplot(data) +
					geom_ribbon(aes_string(x="x", ymin="y", ymax="x", colour=factor("dataset", levels=unique("dataset")),
						fill=factor("dataset", levels=unique("dataset")), linetype=NA), alpha=.3) + 
					xlab(xlab) + ylab(ylab) +
					theme(text = element_text(size=fontsize), 
						axis.text = element_text(color = axis.text.color),
						axis.title.x = element_text(vjust = -.5),
						axis.title.y = element_text(angle=0, hjust=1)) +
					scale_x_continuous(labels = percent_format()) +
					scale_y_continuous(labels = percent_format()) + 
					scale_color_hue(legend.title) + scale_fill_hue(legend.title)
				
			if (!export)
				print(plot)
			
			else
				ggsave(export.name, plot, width=width, height=height, units=units)
			
		} else {
					
			plotl <- mapply(precintcon.plot.lorenz.ungrouped, 
				l, varl,  
					fontsize = fontsize, 
					axis.text.color = axis.text.color, 
					interval = interval,
					MoreArgs = list(xlab=xlab,
						ylab = ylab),
					SIMPLIFY = FALSE)
		
			for (i in 1:length(plotl)) {
				
				if (!export)
					print(plotl[[i]])
					
				else
					ggsave(paste(varl[[i]], export.name, sep="_"), plotl[[i]], 
							width=width, height=height, units=units)
			}
		}
		
	} else
		stop("all input data should be either of type precintcon.daily, precintcon.fd, or precintcon.classified")

	par(ask=F)
}

checkType <- function(f) {
	
	n <- length(f)
	
	if (n > 0) {
		type <- is.element("precintcon.daily", class(f[[1]])) || 
				is.element("precintcon.fd", class(f[[1]]))    || 
				is.element("precintcon.classified", class(f[[1]]))
		
		if (n > 1)
			return(type && checkType(f[-1]))
		else
			return(type)
	} else
		stop("empty input in precintcon.plot.lorenz function")
}