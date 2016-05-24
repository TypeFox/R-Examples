QC_histogram <-
function(
  dataset, data_col = 1, save_name = "dataset", save_dir = getwd(), export_outliers = FALSE,
  filter_FRQ = NULL, filter_cal = NULL, filter_HWE = NULL, filter_imp = NULL,
  filter_NA = TRUE, filter_NA_FRQ = filter_NA, filter_NA_cal = filter_NA, filter_NA_HWE = filter_NA, filter_NA_imp = filter_NA,
  breaks = "Sturges", graph_name = colnames(dataset)[data_col], header_translations,
  check_impstatus = FALSE, ignore_impstatus = FALSE, T_strings = c("1", "TRUE", "yes", "YES", "y", "Y"), F_strings = c("0", "FALSE", "no", "NO", "n", "N"), NA_strings = c(NA, "NA", ".", "-"), ...
) {
  
  skip_FRQ <- if(is.null(filter_FRQ)) { TRUE } else { is.na(filter_FRQ) & !filter_NA_FRQ }
  skip_cal <- if(is.null(filter_cal)) { TRUE } else { is.na(filter_cal) & !filter_NA_cal }
  skip_HWE <- if(is.null(filter_HWE)) { TRUE } else { is.na(filter_HWE) & !filter_NA_HWE }
  skip_imp <- if(is.null(filter_imp)) { TRUE } else { is.na(filter_imp) & !filter_NA_imp }
  
  # This is to ensure that HQ_filter won't be looking at missing columns
  if(skip_FRQ) filter_FRQ <- NULL
  if(skip_cal) filter_cal <- NULL
  if(skip_HWE) filter_HWE <- NULL
  if(skip_imp) filter_imp <- NULL
  
  if(is.vector(dataset)) {
    if(check_impstatus | ignore_impstatus) stop("cannot check or ignore imp-status: vector dataset!")
    if(skip_FRQ + skip_cal + skip_HWE + skip_imp > 2L) {
      dataset <- data.frame(EFFECT = dataset)
      if(is.character(data_col)) {
        colnames(dataset) <- data_col
        data_col <- 1L
      } else {
        if(data_col != 1L) stop("Invalid column specified")
      }
    } else { stop("Insufficient data to apply filters: dataset is single column!") }
  } else {
    if(is.character(data_col)) {
      data_col <- which(colnames(dataset) == data_col)
      if(length(data_col) != 1L) stop("Invalid column specified")
    } else {
      if(is.na(colnames(dataset)[data_col])) stop("Invalid column specified")
  } }
  
  if(length(graph_name) != 1L) stop("Argument 'graph_name' has invalid length")
  # This line was added not to test graph_name, but to "fix" it before
  #	the header of dataset is checked/translated
  
  if(skip_FRQ & skip_cal & skip_HWE & skip_imp) {
    goodOnes <- !is.na(dataset[ , data_col])
    clarf <- "No filter applied"
  } else {
    header_std <- c("EFF_ALL_FREQ", "HWE_PVAL", "CALLRATE", "IMP_QUALITY", "IMPUTED")[c(!skip_FRQ, !skip_HWE, !skip_cal, !skip_imp, check_impstatus | (!ignore_impstatus & !(skip_cal & skip_HWE & skip_imp)))]
    if(missing(header_translations)) {
      if(!any(colnames(dataset) == "IMPUTED")) {
        if(!check_impstatus & (ignore_impstatus | skip_imp | (skip_HWE & skip_cal)) ) {
          if(!ignore_impstatus) {
            if(skip_imp & !skip_HWE & !skip_cal) {
              dataset$IMPUTED <- 0L
              print("Warning: no imputation-status specified - all SNPs set to genotyped") }
            if(skip_HWE & skip_cal & !skip_imp) {
              dataset$IMPUTED <- 1L
              print("Warning: no imputation-status specified - all SNPs set to imputed") }
          }
        } else { stop("Missing imputation status") }
      }
      if(!all(header_std %in% colnames(dataset))) stop("Cannot apply filter: missing or unidentified columns")
    } else {
      header_test <- translate_header(header = colnames(dataset), standard = header_std, alternative = header_translations)
      if(any(duplicated(header_test$header_h))) stop("cannot translate header - duplicate column names")
      if(header_test$missing_N > 1L) stop("cannot translate header - missing columns")
      if(header_test$missing_N == 1L) {
        if(header_test$missing_h == "IMPUTED" & !check_impstatus & (ignore_impstatus | skip_imp | (skip_HWE & skip_cal)) ) {
          if(!ignore_impstatus) {
            if(skip_imp) {
              dataset$IMPUTED <- 0L
              print("Warning: no imputation-status specified - all SNPs set to genotyped")
            } else {
              dataset$IMPUTED <- 1L
              print("Warning: no imputation-status specified - all SNPs set to imputed") }
            header_test$header_h <- c(header_test$header_h, "IMPUTED")
          }
        } else { stop(paste("cannot translate header - missing column:", paste(header_test$missing_h, collapse = ", "))) }
      }
      colnames(dataset) <- header_test$header_h
    }
    
    if(check_impstatus) {
      dataset$IMPUTED <- convert_impstatus(dataset$IMPUTED, T_strings, F_strings, NA_strings, use_log = FALSE)
      if(all(is.na(dataset$IMPUTED))) stop("imputation status missing or untranslated")
    }
    
    goodOnes <- !is.na(dataset[ , data_col]) & HQ_filter(data = dataset, ignore_impstatus = ignore_impstatus,
                                                         FRQ_val = filter_FRQ,	cal_val = filter_cal,	 HWE_val = filter_HWE,	imp_val = filter_imp,
                                                         FRQ_NA	= filter_NA_FRQ, cal_NA = filter_NA_cal, HWE_NA = filter_NA_HWE, imp_NA = filter_NA_imp)
    clarf <- "Filtered for"
    if(!skip_FRQ) {
      if(is.na(filter_FRQ)) { clarf <- paste(clarf, "missing allele frequency;")
      } else {
        if(filter_NA_FRQ) {
          clarf <- paste(clarf, "MAF <", filter_FRQ, "or missing;")
        } else {
          clarf <- paste(clarf, "MAF <", filter_FRQ, ";")
        } } }
    if(!skip_cal) {
      if(is.na(filter_cal)) { clarf <- paste(clarf, "missing call rates;")
      } else {
        if(filter_NA_cal) {
          clarf <- paste(clarf, "call rate <", filter_cal, "or missing;")
        } else {
          clarf <- paste(clarf, "call rate <", filter_cal, ";")
        } } }
    if(!skip_HWE) {
      if(is.na(filter_HWE)) { clarf <- paste(clarf, "missing HWE p-value;")
      } else {
        if(filter_NA_HWE) {
          clarf <- paste(clarf, "HWE p <", filter_HWE, "or missing;")
        } else {
          clarf <- paste(clarf, "HWE p <", filter_HWE, ";")
        } } }
    if(!skip_imp) {
      if(is.na(filter_imp)) { clarf <- paste(clarf, "missing imputation quality;")
      } else {
        if(filter_NA_imp) {
          clarf <- paste(clarf, "imp. qual. <", filter_imp, "or missing;")
        } else {
          clarf <- paste(clarf, "imp. qual. <", filter_imp, ";")
        } } }
    clarf <- substr(clarf, 1L, nchar(clarf) - 1L) # removes the final semi-colon
  }
  
  goodN <- sum(goodOnes)
  if(goodN < 4L) { print("Insufficient non-missing, non-filtered effect sizes")
  } else {
    min_dat <- min(dataset[goodOnes, data_col])
    max_dat <- max(dataset[goodOnes, data_col])
    min_N <- 0L
    max_N <- 0L
    
    png(paste0(save_dir, "/", save_name, ".png"), width = 1440, height = 480)
    par(mfrow = c(1, 2))
    (( h1<-hist(mean(dataset[goodOnes, data_col]) + (qnorm(ppoints(goodN)) * sd(dataset[goodOnes, data_col])),
                freq = FALSE, plot = TRUE, main = paste("Expected distribution of", graph_name), xlab = graph_name, breaks = breaks, sub = save_name, font.sub = 3, ...) ))
    h2_breaks <- h1$breaks
    minbreaks <- h2_breaks[1]
    maxbreaks <- h2_breaks[length(h2_breaks)]
    if (minbreaks > min_dat) {
      h2_breaks <- c(min_dat, h2_breaks)
      min_N <- sum(dataset[goodOnes, data_col] < minbreaks)
    }
    if (maxbreaks < max_dat) {
      h2_breaks <- c(h2_breaks, max_dat)
      max_N <- sum(dataset[goodOnes, data_col] > maxbreaks)
    }
    (( h2 <- hist(dataset[goodOnes, data_col], breaks = h2_breaks, xlim = c(minbreaks, maxbreaks), 
                  freq = FALSE, plot = TRUE, main = paste("Observed distribution of", graph_name), xlab = graph_name, sub = clarf, font.sub = 3, ...) ))
    if(min_N > 0L) {
      text(minbreaks, 0.6 * max(h2$density), pos = 4,
           label = paste(min_N, "values outside min. range"), cex = 1, col = "red")
    }
    if(max_N > 0L) {
      text(maxbreaks, 0.6 * max(h2$density), pos = 2,
           label = paste(max_N, "values outside max. range"), cex = 1, col = "red")
    }
    dev.off()
    
    if(export_outliers > 0L & min_N + max_N > 0L) {
      if(min_N + max_N <= export_outliers | export_outliers == 1) {
        write.table(dataset[goodOnes & (dataset[ , data_col] < minbreaks | dataset[ , data_col] > maxbreaks), ],
                    paste0(save_dir, "/", save_name, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
      } else {
        write.table(dataset[goodOnes & (dataset[ , data_col] < minbreaks | dataset[ , data_col] > maxbreaks), ][1:export_outliers, ],
                    paste0(save_dir, "/", save_name, ".txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  } } }
  return(invisible())
}
