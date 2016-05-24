#' @export 
precintcon.plot.deciles <- function(
  ...,
	ylab            = "Precipitation",
	legend.title    = "Legend",
	legend          = NULL,
	fontsize        = 10, 
	axis.text.color = "black", 
	export          = FALSE, 
	export.name     = "deciles_plot.png", 
	width           = 8.6, 
	height          = 7.5, 
	units           = "cm",
	ylim            = NA,
	grouped         = FALSE,
	args            = NA
) {
	
	l <- list(...)
	
	if (length(l) > 1 && !export && !grouped)
		par(ask = T)
	
	varl <- ifelse(is.na(args), as.character(match.call()[1:length(l)+1]), args)

	if (!is.null(legend) && length(varl) != length(legend))
		stop(paste("legend should has length equals to the number of input data. legend parameter length", 
						length(legend), ": number of input data", length(varl)))
	
	else if (!is.null(legend))
		varl <- as.list(legend)

	if (!grouped) {
		min <- min(unlist(lapply(l, FUN = function(d) return(as.precintcon.deciles(d)[2]))))
		max <- max(unlist(lapply(l, FUN = function(d) return(as.precintcon.deciles(d)[11]))))

      #######
      #
      # Function for plotting
      #
      f <- function(d, legend, max, min,
            ylab, fontsize, axis.text.color, 
            export, export.name, width, height, units) {
         
         d <- as.precintcon.deciles(d)
         
         d[1] <- legend
         
         graph <- ggplot(d) + 
               geom_boxplot(aes_string(x = "dataset", ymin = "D1", lower = "D2",
                           middle = "D5", upper = "D9", ymax = "D10"), stat = "identity", show.legend = FALSE) +
               ylab(ylab) +
               theme(text = element_text(size = fontsize), 
                     axis.text = element_text(color = axis.text.color),
                     axis.title.x = element_blank(),
                     axis.text.x  = element_blank(),
                     axis.ticks.x = element_blank()) +
               ylim(min, max) +
               facet_grid(. ~ dataset)
         
         if (!export) {
            print(graph)
         } else {
            export.name <- paste(export.name, sep = "_")
            ggsave(export.name, graph, width = width, height = height, units = units)
         }
      }
      
		###########
		#
		# Generating graphs for each input data
		#	
		mapply(f, l, varl, max = max, min = min, fontsize = fontsize, axis.text.color = axis.text.color, 
		   export = export, export.name = export.name, width = width, height = height, units = units,
		   MoreArgs = list(ylab = ylab), SIMPLIFY = FALSE)
	}
	
	par(ask = F)
}
