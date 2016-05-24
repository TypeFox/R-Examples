plot_precision <-
function(SE, N, labels = NULL,
                           save_name = "Graph_precision", save_dir = getwd(), ...) {
  plot_table <- data.frame(SE = SE, N = N, la = if(is.null(labels)) FALSE else labels, stringsAsFactors = FALSE)
  if(sum(!is.na(plot_table$SE) & !is.na(plot_table$N)) > 1L) {
    plot_table <- subset(plot_table, !is.na(plot_table$SE) & !is.na(plot_table$N))
    plot_table$SE <- 1/plot_table$SE
    plot_table$N <- sqrt(plot_table$N)
    
    png(paste0(save_dir, "/", save_name, ".png") )
    plot(plot_table$N, plot_table$SE, main = "Precision by Sample Size", xlab = "sqrt(sample size)", ylab = "1 / median(SE)", ...)
    if(!is.null(labels)) text(plot_table$N, plot_table$SE, labels = plot_table$la, pos = 4)
    dev.off()
  } else { print(" - - Insufficient data to create scatterplot Standard Error vs Sample Size", quote = FALSE) }
  return(invisible())
}
