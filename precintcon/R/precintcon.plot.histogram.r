#' @export
precintcon.plot.histogram <- function(
   ..., 
	density         = FALSE, 
	grouped         = FALSE,
	xlab            = "Precipitation (mm)",
	ylab            = "Frequency",
   legend.title   = "Legend", 
   legend         = NULL,
	fontsize        = 10,
	axis.text.color = "black", 
	export          = FALSE, 
	export.name     = "histogram_plot.png", 
	width           = 10, 
	height          = 10, 
	units           = "cm",
	args            = NA
) {
	
	l <- list(...)
	
	if (length(l) > 0) {
	
		if (length(l) > 1 && !export && !grouped)
			par(ask = TRUE)
		
		varl <- ifelse(is.na(args), as.character(match.call()[1:length(l)+1]), args)

		if (!is.null(legend) && length(varl) != length(legend)) {
			stop(paste("legend should has length equals to the number of input data. legend parameter length", 
				length(legend), ": number of input data", length(varl)))
		} else if (!is.null(legend)) {
			varl <- as.list(legend)
		}
		
		#######
		#
		# Defining function for x data
		#
		f <- function(x, n) {
			
			ddd <- NULL
			
			if (is.element("precintcon.daily", class(x))) {
				ddd <- as.vector(as.matrix(x[ ,3:33]))
			}  else if (is.element("precintcon.monthly", class(x))) {
				ddd <- as.vector(x$precipitation)
			} else if (is.element("precintcon.spi", class(x))) {
				ddd <-as.vector(x$spi)
			} else if (is.element("precintcon.pci", class(x))) {
				ddd <-as.vector(x$pci)
			} else {
				stop("input data should be either of type precintcon.daily, precintcon.monthly, precintcon.spi, or precintcon.pci")
			}

			ddd <- ddd[!is.na(ddd)]
			
			data.frame(dataset = paste(n, sep = ""), x = ddd)
		}
		
		l <- mapply(f, l, varl, SIMPLIFY = FALSE)
		
		if (!grouped) {
			
			#######
			#
			# Defining function for graph generation
			#
			g <- function(d, n, 
					density, ylab, xlab, export.name, export, width, height, units) {	

				graph <- ggplot(d, aes_string(x = "x")) + 
					(if (!density) geom_histogram() else geom_density(fill = "grey", alpha = .3, show_guide = F)) + 
					xlab(xlab) + ylab(ylab) +
					theme(text = element_text(size = fontsize),  
						  axis.text = element_text(color = axis.text.color)) + 
					facet_grid(. ~ dataset)
				
				if (!export) {
					suppressMessages(print(graph))
				} else {
					export.name <- paste(n, export.name, sep = "_")
					ggsave(export.name, graph, width = width, height = height, units = units)
				}
			}
			
			#########
			#
			# Generating graphs for each input data
			#
			mapply(g, l, varl, density = density, export = export, 
				export.name = export.name, width = width, height = height, units = units,
				MoreArgs = list(ylab = ylab, xlab = xlab),
			   	SIMPLIFY = FALSE)
   
   		} else {
			
			######
			#
			# list of data.frame to data.frame (rbind)
			#
			l <- do.call(rbind.data.frame, l)
			
			graph <- ggplot(l, aes_string(x = "x", fill = "dataset")) + 
					 (if (!density) geom_histogram() else geom_density(alpha = .3)) + 
					 xlab(xlab) + ylab(ylab) +
					 scale_fill_discrete(legend.title) +
					 theme(text = element_text(size = fontsize),
						axis.text = element_text(color = axis.text.color))
			 
			 if (!export) {
				 suppressMessages(print(graph))
			 } else {
				 ggsave(export.name, graph, height = height, width = width, units = units)
			 }
		}
	} else {
		stop("empty input data in precitcon.plot.histogram function.")
	}
	
	par(ask = FALSE)
}
