plot_regional <-
function(dataset, chr, start_pos, end_pos,
                          plot_cutoff_p = 1, name_cutoff_p,
                          data_name = NULL, save_name = "regional_association_plot", save_dir = getwd(), header_translations,
                          main = "Regional association plot", ...) {
  plot_names <- !missing(name_cutoff_p)
  header_std <- c("PVALUE", "CHR", "POSITION", "MARKER")[c(TRUE, TRUE, TRUE, plot_names)]
  # marker-names goes last so that we know the index of p, chr and position
  
  if(missing(header_translations)) {
    dataset <- dataset[which(dataset$PVALUE <= plot_cutoff_p & dataset$CHR == chr & dataset$POSITION >= start_pos & dataset$POSITION <= end_pos), header_std[-2]]
  } else {
    header_data <- toupper(colnames(dataset))
    header_col <- integer(length = length(header_std))
    for(forI in 1:length(header_std)) {
      header_current <- identify_column(header_std[forI], header_translations, header_data)
      if(length(header_current) != 1L) { if(length(header_current) == 0L) { stop(paste("Cannot identify data column:", header_std[forI])) } else { stop(paste("Multiple data columns identified as:", header_std[forI])) } }
      header_col[forI] <- header_current
    }
    dataset <- dataset[which(dataset[ , header_col[1]] <= plot_cutoff_p & dataset[ , header_col[2]] == chr & dataset[ , header_col[3]] >= start_pos & dataset[ , header_col[3]] <= end_pos), header_col[-2]]
    colnames(dataset) <- header_std[-2]
  }
  
  if(nrow(dataset) < 10L) stop("Insufficient usuable entries")
  
  png(paste0(save_dir, "/", save_name, ".png"),
       width = 960, height = 480)
  par(mgp = c(2.5, 0.9, 0))
  plot(dataset$POSITION, -log10(dataset$PVALUE), pch = 20, xaxs = "r",
       xlab = paste("Position on chromosome", chr), ylab = "Observed -log10(p-value)",
       main = main, cex.lab = 1.2, sub = data_name, cex.sub = 1.3, ...)
  abline(h = -log10(5e-8), lty = 3, col="red")
  if(plot_names) { if(any(dataset$PVALUE <= name_cutoff_p)) {
    dataset <- dataset[dataset$PVALUE <= name_cutoff_p, ]
    text(dataset$POSITION, -log10(dataset$PVALUE), labels = dataset$MARKER, pos = 4)
  } }
  dev.off()
  return(invisible())
}
