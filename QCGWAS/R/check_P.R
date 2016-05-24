check_P <-
function(dataset, HQ_subset,
                    plot_correlation = FALSE, plot_if_threshold = FALSE, threshold_r = 0.99,
                    save_name = "dataset", save_dir = getwd(), header_translations,
                    use_log = FALSE, dataN = nrow(dataset), ...) {
  
  if(!missing(header_translations)) {
    header_test <- translate_header(header = colnames(dataset), standard = c("PVALUE", "EFFECT", "STDERR"), alternative = header_translations)
    if(any(duplicated(header_test$header_h))) stop("cannot translate header - duplicate column names")
    if(header_test$missing_N > 0L) stop(paste("Cannot identify data column(s):", paste(header_test$missing_h, collapse = ", ")))
    colnames(dataset) <- header_test$header_h
  }
  
  if(missing(HQ_subset)) {
    HQ_subset <- logical(length = nrow(dataset))
    plot_HQ <- FALSE
  } else {
    if(is.numeric(HQ_subset)) {
      temp <- logical(length = nrow(dataset))
      temp[HQ_subset] <- TRUE
      HQ_subset <- temp
    }
    plot_HQ <- TRUE
  }
  
  goodOnes <- !(is.na(dataset$PVALUE) | is.na(dataset$EFFECT) | is.na(dataset$STDERR))
  if(sum(goodOnes) > 10L) {
    p_obs <- -log10(dataset$PVALUE[goodOnes])
    p_exp <- -log10(pchisq((dataset$EFFECT[goodOnes]/dataset$STDERR[goodOnes])^2, 1, lower.tail=FALSE))
    HQ_subset <- HQ_subset[goodOnes]
    if(any(p_obs > 300)) {
      p_obs[p_obs > 300] <- 300
      if(use_log) { save_log(phaseL = 4L, checkL = "p-values", typeL = "extreme values", SNPL = sum(p_obs > 300), allSNPs = dataN, actionL = "-", noteL = "Extreme p-values in dataset - temporarly set to 1e-300 for correlation calculation and plot", fileL = paste(save_dir, save_name, sep = "/"))
      } else { print(" - - warning: Extreme p-values in dataset - temporarly set to 1e-300 for correlation calculation and plot", quote = FALSE) }
    }
    if(any(p_exp > 300)) {
      p_exp[p_exp > 300] <- 300
      if(use_log) { save_log(phaseL = 4L, checkL = "p-values", typeL = "extreme values", SNPL = sum(p_exp > 300), allSNPs = dataN, actionL = "-", noteL = "Expexted p-values capped at 1e-300", fileL = paste(save_dir, save_name, sep = "/"))		
      } else { print(" - - warning: expected p-values capped at 1e-300", quote = FALSE) }
    }
    p_cor <- cor(p_obs, p_exp, use = "everything")
    if(is.na(p_cor)){
      if(use_log) save_log(phaseL = 4L, checkL = "p-value check", typeL = "unable to calculate", SNPL = dataN, allSNPs = dataN, actionL = "-", noteL = "Unable to calculate correlation - check warnings()", fileL = paste(save_dir, save_name, sep = "/"))
      print(" - - warning: unable to calculate p-value correlation. Check warnings().", quote = FALSE)
      return(NA)      
    }
    if(p_cor < threshold_r) {
      if(use_log) save_log(phaseL = 4L, checkL = "p-value check", typeL = "poor correlation", SNPL = dataN, allSNPs = dataN, actionL = "-", noteL = paste("Reported p-values correlate poorly with recalculated p-values ( r =", round(p_cor, digits = 3),")"), fileL = paste(save_dir, save_name, sep = "/"))
      print(paste0(" - - warning: reported p-values correlate poorly to expected values(r = ", round(p_cor, digits = 3), ")"), quote = FALSE)
    }
    if(plot_correlation & (p_cor < threshold_r | !plot_if_threshold )) {
      p_max <- max(p_obs, p_exp)
      png(paste0(save_dir, "/", save_name, "_graph_p-correlation.png"),
           width = 720, height = 720, res = 144)
      plot(x = 0, y = 0, col = "white",
           xlim = c(0, p_max), ylim = c(0, p_max),
           main="P-value correlation", xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
           sub = save_name, cex.sub = 1.3, ...)
      lines(x = c(0, 1.1 * p_max), y = c(0, 1.1 * p_max))
      points(p_exp[!HQ_subset],p_obs[!HQ_subset], col = if(plot_HQ) "grey" else "black", pch = 20, cex = 0.8)      
      points(p_exp[HQ_subset], p_obs[HQ_subset], col = "black", pch = 20, cex = 0.8)
      if (plot_HQ) legend(0.1 * p_max, 0.9 * p_max, c(" Low quality", "High quality"), pch = 20, col = c("grey", "black"), cex = 0.8)
      text(0.1 * p_max, 0.95 * p_max, paste("r =", round(p_cor, digits = 3)), pos = 4, cex=1.0, col = ifelse(p_cor < threshold_r, "red", "black") )
      dev.off()
    }
  } else {
    if(use_log) { save_log(phaseL = 4L, checkL = "p-value check", typeL = "insufficient data", SNPL = dataN, allSNPs = dataN, actionL = "Test skipped", noteL = "Less than 10 entries with sufficient data to calculate corelation", fileL = paste(save_dir, save_name, sep = "/"))
    } else { print(" - - p-test aborted: less than 10 entries with the data required to calculate corelation", quote = FALSE) }
    p_cor <- NA
  }
  return(p_cor)
}
