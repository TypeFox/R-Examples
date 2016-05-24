plot_skewness <-
function(skewness, kurtosis, labels = paste("Study", 1:length(skewness)), plot_labels = "outliers",
                          save_name = "Graph_skewness_kurtosis", save_dir = getwd(), ...) {
  
  if(is.logical(plot_labels)) { if(plot_labels) { plot_labels <- "all" } else { plot_labels <- "none" } }
  
  data <- data.frame(sk = skewness, ku = kurtosis, la = labels, stringsAsFactors = FALSE)
  data <- subset(data, !is.na(data$sk) & !is.na(data$ku)) 
  if(nrow(data) <= 1L) { print(" - - Insufficient data to create skewness/kurtosis plot.", quote = FALSE)
  } else {
    to_label <- logical(length = nrow(data))
    if(plot_labels == "all" ) to_label <- TRUE
    if(plot_labels == "outliers" ) to_label <- abs(data$sk) > 0.1 | data$ku > 10
    png(paste0(save_dir, "/", save_name, ".png") )
    plot(data$sk, data$ku, main = "Skewness vs. Kurtosis plot", xlab = "Skewness", ylab = "kurtosis", ...)
    abline(h = 0, v = 0, lwd = 2)
    if(any(to_label)) text(data$sk[to_label], data$ku[to_label], labels = data$la[to_label], pos = 4)
    dev.off()
  }
  return(invisible())
}
