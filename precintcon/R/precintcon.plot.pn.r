#' @export
precintcon.plot.pn <- function(..., 
  interval        = 30, 
  scale           = "a",
	xlab            = NA, 
  ylab            = "PN", 
  fontsize        = 10, 
	axis.text.color = "black", 
  legend          = NULL, 
	export          = FALSE, 
  export.name     = "pn_plot.png", 
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
	
	mapply(function(d, n, interval, scale, axis.text.color, fontsize, 
					xlab, ylab, export, export.name, width, height, units) {
				
		if (is.element("precintcon.pn", class(d))      || 
			is.element("precintcon.monthly", class(d)) || 
			is.element("precintcon.daily", class(d))) {
			
			d <- precintcon.pn.analysis(d, interval, scale)
			
			d <- cbind(d, data.frame(dataset=paste(n)))
			
			graph <- NA
			
			if ((scale == "a") || (scale == "d")) {
				
				if (is.na(xlab))
					xlab <- "Year"
				
				graph <- ggplot(d, aes_string("year", "pn")) + geom_point(size=2) +
					scale_x_continuous(expand = c(.02, .02), breaks = seq(min(d$year), max(d$year), by = 2))
			
			} else {
				
				ddd <- NA
				
				if (scale == "s") {
					
					if (is.na(xlab))
						xlab = "Season"
					
					ddd <- 2 * (d$season - 1) + d$season 
				} else if (scale == "m"){
					
					if (is.na(xlab))
						xlab = "Month"
					
					ddd <- d$month 
				}
				
				d <- 
					cbind(
						d, data.frame(
							date = as.Date(paste("01", ddd, d$year, sep = "/"), "%d/%m/%Y")));

				graph <- 
					ggplot(d, aes_string("date", "pn")) + geom_point(size = 1.1) +
					scale_x_date(labels = date_format("%b %y"))
			}
			
			graph <- 
				graph + 
				geom_line(size=.5) + 
				xlab(xlab) +
				ylab(ylab) +
				theme(text = element_text(size = fontsize), 
					axis.text = element_text(color = axis.text.color),
					axis.text.x = element_text(angle = 25),
					axis.title.x = element_text(vjust = .1)) +
				facet_grid(. ~ dataset)
			
			if (!export) {
				print(graph)
			} else {
				export.name <- paste(n, export.name, sep="_") 
				ggsave(export.name, plot = graph, height = height, width = width, units = units)
			}
		}
	}, l, varl, interval = interval, scale = scale,
	axis.text.color = axis.text.color, MoreArgs = list(fontsize = fontsize, 
	width = width, height = height, units = units, xlab = xlab, ylab = ylab, 
	export = export, export.name = export.name), SIMPLIFY = FALSE)
	
	par(ask=F)
}
