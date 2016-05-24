#' @export
precintcon.plot.pci <- function(..., 
  xlab            = "Years", 
  ylab            = "PCI", 
  legend          = NULL,
  fontsize        = 10, 
	axis.text.color = "black",  
	export          = FALSE, 
  export.name     = "pci_plot.png", 
	width           = 10, 
  height          = 10, 
  units           = "cm",
  args            = NA
) {
	
	l <- list(...)
	
	varl <- ifelse(is.na(args), as.character(match.call()[1:length(l)+1]), args)
	
	if (length(l) > 1 && !export)
		par(ask=T)
	
	if (!is.null(legend) && length(varl) != length(legend)) {
		
		stop(paste("legend should has length equals to the number of input data. legend parameter length", 
						length(legend), ": number of input data", length(varl)))
		
	} else if (!is.null(legend))
		varl <- as.vector(legend)

	mapply(function(d, n, axis.text.color, fontsize, 
		xlab, ylab, export, export.name, width, height, units) {
	
		if (is.element("precintcon.pci", class(d)) || 
				is.element("precintcon.monthly", class(d)) || 
				is.element("precintcon.daily", class(d))) {
			
			if (is.element("precintcon.daily", class(d))) {
				d <- as.precintcon.monthly(d)
			}
			
			if (is.element("precintcon.monthly", class(d))) {
				d <- precintcon.pci.analysis(d)
			}
		
				d <- cbind(d, data.frame(dataset=paste(n)))
			
				graph <- ggplot(d, aes_string("year", "pci")) + geom_line(size=.5) +
					geom_point(size=2) + xlab(xlab) + ylab(ylab) +
					theme(text = element_text(size = fontsize), 
							axis.text = element_text(color = axis.text.color),
							axis.text.x = element_text(angle = 25),
							axis.title.x = element_text(vjust = .1)) +
					scale_x_continuous(expand = c(.02, .02), 
							breaks = seq(min(d$year), max(d$year), by = 2)) +
					facet_grid(. ~ dataset)
			
				if (!export) {
					print(graph)
				} else {
					export.name <- paste(n, export.name, sep="_") 
					ggsave(export.name, plot=graph, height=height, width=width, units=units)
				}
		
		}
	}, l, varl,
	   axis.text.color = axis.text.color, fontsize = fontsize, width = width, height = height, units = units, MoreArgs = list(xlab = xlab, ylab = ylab, 
	   export = export, export.name = export.name), SIMPLIFY = FALSE)

	par(ask=F)
}
