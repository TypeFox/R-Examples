plot_distribution <-
function(data_table, names = 1:ncol(data_table), include = TRUE, plot_order = 1:ncol(data_table),
                              quantile_lines = FALSE, save_name = "Graph_distribution", save_dir = getwd(), ...) {
  
  table1 <- data.frame(name = names, use = include, order = plot_order, N = colSums(!is.na(data_table)))
  table2 <- data_table[ , table1$use & table1$N > 0]
  table1 <- table1[table1$use & table1$N > 0, ]
  
  plotN <- nrow(table1)
  
  if(plotN > 1L) {
    
    if(quantile_lines) {
      pointN <- nrow(table2)
      quantL <- matrix(data = 0.0, ncol = 3, nrow = plotN,
                       dimnames = list(NULL, c("quantile25", "median", "quantile75")))
      for(medI in 1:plotN) quantL[medI, ] <- quantile(table2[ , medI], names = FALSE, na.rm = table1$N[medI] < pointN)[2:4]
    }
    
    png(paste0(save_dir, "/", save_name, ".png"), width = 400 + 80 * plotN, height = 480)
    boxplot(table2[ , order(table1$order, na.last = FALSE)], names = table1$name[order(table1$order, na.last = FALSE)], ...)
    if(quantile_lines) abline(h = c(median(quantL[ ,1]), median(quantL[ ,2]), median(quantL[ ,3])), lty = 3)
    dev.off()
  } else { print(" - - Insufficient data for effect-size comparisons", quote = FALSE) }
  return(invisible())
}
