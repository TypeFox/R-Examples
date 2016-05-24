QQ_plot <-
function(
  dataset, save_name = "dataset", save_dir = getwd(),
  filter_FRQ = NULL, filter_cal = NULL, filter_HWE = NULL, filter_imp = NULL,
  filter_NA = TRUE, filter_NA_FRQ = filter_NA, filter_NA_cal = filter_NA, filter_NA_HWE = filter_NA, filter_NA_imp = filter_NA,
  p_cutoff = 0.05, plot_QQ_bands = FALSE,
  header_translations,
  check_impstatus = FALSE, ignore_impstatus = FALSE, T_strings = c("1", "TRUE", "yes", "YES", "y", "Y"), F_strings = c("0", "FALSE", "no", "NO", "n", "N"), NA_strings = c(NA, "NA", ".", "-"), ... ) {
  
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
    if(skip_FRQ & skip_cal & skip_HWE & skip_imp) dataset <- data.frame(PVALUE = dataset) else stop("Insufficient data to apply filters: dataset is single column!")
    if(check_impstatus) stop("cannot check impstatus - dataset is a single column") }
    
  header_std <- c("PVALUE", "EFF_ALL_FREQ", "HWE_PVAL", "CALLRATE", "IMP_QUALITY", "IMPUTED")[c(TRUE, !skip_FRQ, !skip_HWE, !skip_cal, !skip_imp, check_impstatus | (!ignore_impstatus & !(skip_cal & skip_HWE & skip_imp)))]
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
    if(!all(header_std %in% colnames(dataset))) { stop("Cannot apply filter: missing or unidentified columns") }
  } else {
    header_test <- translate_header(header = colnames(dataset), standard = header_std, alternative = header_translations)
    if(any(duplicated(header_test$header_h))) { stop("cannot translate header - duplicate column names") }
    if(header_test$missing_N > 1) { stop("cannot translate header - missing columns") }
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
    if(all(is.na(dataset$IMPUTED))) { stop("imputation status missing or untranslated") } }
  
  if(skip_FRQ & skip_cal & skip_HWE & skip_imp) {
    goodOnes <- !is.na(dataset$PVALUE)
    clarf <- "No filter applied"
  } else {
    goodOnes <- !is.na(dataset$PVALUE) & HQ_filter(data = dataset, ignore_impstatus = ignore_impstatus,
                                                   FRQ_val = filter_FRQ,	cal_val = filter_cal,	 HWE_val = filter_HWE,	imp_val = filter_imp,
                                                   FRQ_NA  = filter_NA_FRQ, cal_NA = filter_NA_cal, HWE_NA = filter_NA_HWE, imp_NA = filter_NA_imp)
    clarf <- "Filtered for"
    if(!skip_FRQ) {
      if(is.na(filter_FRQ)) { clarf <- paste(clarf, "missing allele-FRQ;")
      } else {
        if(filter_NA_FRQ) {
          clarf <- paste(clarf, "FRQ <", filter_FRQ, "or missing;")
        } else {
          clarf <- paste(clarf, "FRQ <", filter_FRQ, ";")
        } } }
    if(!skip_cal) {
      if(is.na(filter_cal)) { clarf <- paste(clarf, "missing callrates;")
      } else {
        if(filter_NA_cal) {
          clarf <- paste(clarf, "callrate <", filter_cal, "or missing;")
        } else {
          clarf <- paste(clarf, "callrate <", filter_cal, ";")
        } } }
    if(!skip_HWE) {
      if(is.na(filter_HWE)) { clarf <- paste(clarf, "missing HWE-p;")
      } else {
        if(filter_NA_HWE) {
          clarf <- paste(clarf, "HWE-p <", filter_HWE, "or missing;")
        } else {
          clarf <- paste(clarf, "HWE-p <", filter_HWE, ";")
        } } }
    if(!skip_imp) {
      if(is.na(filter_imp)) { clarf <- paste(clarf, "missing imputation-Q;")
      } else {
        if(filter_NA_imp) {
          clarf <- paste(clarf, "impQ <", filter_imp, "or missing;")
        } else {
          clarf <- paste(clarf, "impQ <", filter_imp, ";")
        } } }
    clarf <- substr(clarf, 1, nchar(clarf) - 1) # removes the final semi-colon
  }
  
  P_obs_N <- sum(goodOnes)
  if(P_obs_N < 10) {
    print("Error - insuficient non-missing, non-filtered p-values")
    return(invisible(NULL))
  }
  
  P_obs <- dataset$PVALUE[goodOnes]
  P_obs <- sort(-log10(P_obs))
  incl  <- P_obs >= -log10(p_cutoff)
  
  if(sum(incl) < 10) {
    print("Error - insufficient significant p-values after filtering")
    return(invisible(NULL))
  }
  
  P_exp	<- sort(-log10(ppoints(P_obs_N)))
  P_exp_short <- P_exp[incl]
  P_obs_short <- P_obs[incl]
  P_exp_min	<- P_exp_short[1]
  P_exp_max	<- P_exp_short[length(P_exp_short)]
  P_obs_min	<- P_obs_short[1]
  P_obs_max	<- P_obs_short[length(P_obs_short)]
  
  P_min <- if(P_exp_min < P_obs_min) P_exp_min else P_obs_min
  P_max <- if(P_exp_max > P_obs_max) P_exp_max else P_obs_max
  
  if(plot_QQ_bands) {
    temp <- (1:P_obs_N)
    i1000 <- c(1, (1:1000) * floor(P_obs_N / 1000), P_obs_N)
    QQ_band_upper <- sort(-log10(qbeta( 1 - 0.05 / 2, temp, P_obs_N - temp + 1 ) ) )[i1000]
    QQ_band_lower <- sort(-log10(qbeta(     0.05 / 2, temp, P_obs_N - temp + 1 ) ) )[i1000]
    P_exp <- P_exp[i1000]
  }
  
  png(paste0(save_dir, "/", save_name, ".png"), width = 720, height = 720)
  plot(c(P_exp_min, P_exp_max), c(P_obs_min, P_obs_max), xlim = c(0, P_exp_max), ylim = c(0, P_obs_max),
       main = "QQ plot", xlab = "Expected -log10(p-value)", ylab = "Observed -log10(p-value)",
       pch = 20, sub = clarf, font.sub = 3, ...)
  if(plot_QQ_bands) {
    polygon( c(P_exp, rev(P_exp)), c( QQ_band_upper, rev(P_exp)), col="grey", border = NA )
    polygon( c(P_exp, rev(P_exp)), c( QQ_band_lower, rev(P_exp)), col="grey", border = NA )
  }
  
  lines(c(P_min, P_max), c(P_min, P_max))
  points(P_exp_short, P_obs_short, pch = 20, col = "black")
  dev.off()
  return(invisible(NULL))
}
