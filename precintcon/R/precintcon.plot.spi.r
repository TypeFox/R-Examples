#' @export
precintcon.plot.spi <- function(..., 
  period          = 3, 
  distribution    = "Gamma",
	xlab            = "Months",
	ylab            = "SPI", 
	ylim            = c(-3,3),
	legend          = NULL,
	fontsize        = 10, 
	axis.text.color = "black", 
	export          = FALSE, 
	export.name     = "spi_plot.png", 
	width           = 8.6, 
  height          = 7.5, 
  units           = "cm",
  args            = NA
) {

	l <- list(...)
	
	if (length(l) > 0) {
		
		if (length(l) > 1)
			par(ask = T)
		
		l <- lapply(l, FUN = precintcon.spi.analysis, period = period, distribution = distribution)
		
		varl <- ifelse(is.na(args), as.character(match.call()[1:length(l)+1]), args)
		
		if (!is.null(legend) && length(varl) != length(legend))
			stop(paste("legend should has length equals to the number of input data. legend parameter length", 
							length(legend), ": number of input data", length(varl)))
		
		else if (!is.null(legend))
			varl <- as.list(legend)
		
		plotl <- mapply(p.plot.spi, l, varl, 
			fontsize=fontsize, 
			axis.text.color=axis.text.color, 
			MoreArgs = list(xlab=xlab,
				ylab=ylab, ylim=ylim),
			SIMPLIFY=FALSE)
		
		for (i in 1:length(plotl)) {
			
			if (!export)
				print(plotl[[i]])
			
			else
				ggsave(paste(varl[[i]], export.name, sep="_"), plotl[[i]], width=width, height=height, units=units)
		}
					
		par(ask = F)
		
	} else stop("empty input data in precintcon.plot.spi function.")
}

p.plot.spi <- function(d, n,
		xlab            = "Periods",
		ylab            = "SPI", 
		ylim            = c(-3,3),
		fontsize        = 10, 
		axis.text.color = "black"
) {
	
	if (is.element("precintcon.spi", class(d))) {
		
		d[which(nchar(d[,2]) < 2), 2] <- paste("0", d[which(nchar(d[,2]) < 2), 2], sep = "")
		
		data <- data.frame(
			x = as.Date(paste(d[ ,1], d[ ,2], "01", sep = "/"), "%Y/%m/%d"), 
			y = d[ ,3], dataset = paste(n, sep=""))
		
		p <- ggplot(data) + geom_line(aes_string(x = "x", y = "y"), size = .3) + 
				xlab(xlab) + 
		    ylab(ylab) + 
				scale_x_date(expand   = c(1/48, 1/48), 
						limits            = as.Date(c(data[1,1], tail(data$x, n = 1))), 
						date_breaks       = "20 months", labels = date_format("%b %y"), 
						date_minor_breaks = "1 month") +
				scale_y_continuous(breaks = seq(-3, 3, by = 1), limits = c(-3.5, 3.5)) +
				theme(text           = element_text(size = fontsize), 
					axis.text          = element_text(color = axis.text.color),
					axis.text.x        = element_text(angle = 25),
					axis.title.x       = element_text(vjust = .1),
					panel.grid.minor.x = element_blank()) +
			facet_grid(. ~ dataset)

	return(p)
		
	} else {
		stop("invalid input data type. It should be of type precintcon.spi")		
	}
}