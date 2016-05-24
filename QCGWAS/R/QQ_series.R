QQ_series <-
function(dataset, save_name = "dataset", save_dir = getwd(),
                      filter_FRQ = NULL, filter_cal = NULL, filter_HWE = NULL, filter_imp = NULL,
                      filter_NA = TRUE, filter_NA_FRQ = filter_NA, filter_NA_cal = filter_NA, filter_NA_HWE = filter_NA, filter_NA_imp = filter_NA,
                      p_cutoff = 0.05, plot_QQ_bands = FALSE,
                      header_translations,
                      check_impstatus = FALSE, ignore_impstatus = FALSE, T_strings = c("1", "TRUE", "yes", "YES", "y", "Y"), F_strings = c("0", "FALSE", "no", "NO", "n", "N"), NA_strings = c(NA, "NA", ".", "-"), ... ) {
  
  useFRQ <- if(is.null(filter_FRQ)) FALSE else !all(is.na(filter_FRQ) & !filter_NA_FRQ)
  useCal <- if(is.null(filter_cal)) FALSE else !all(is.na(filter_cal) & !filter_NA_cal)
  useHWE <- if(is.null(filter_HWE)) FALSE else !all(is.na(filter_HWE) & !filter_NA_HWE)
  useImp <- if(is.null(filter_imp)) FALSE else !all(is.na(filter_imp) & !filter_NA_imp)
  
  if(missing(header_translations)) {
    if(!any(colnames(dataset) == "IMPUTED")) {
      if(!check_impstatus & (ignore_impstatus | (!useImp | (!useHWE & !useCal)))) {
        if(!ignore_impstatus) {
          if(!useImp & useHWE & useCal) {
            dataset$IMPUTED <- 0
            print("Warning: no imputation-status specified - all SNPs set to genotyped") }
          if(!useHWE & !useCal & useImp) {
            dataset$IMPUTED <- 1
            print("Warning: no imputation-status specified - all SNPs set to imputed") }
        }
      } else { stop("Missing imputation status") }
    }
  } else {    
    header_std <- c("PVALUE", "EFF_ALL_FREQ", "HWE_PVAL", "CALLRATE", "IMP_QUALITY", "IMPUTED")[c(TRUE, useFRQ, useHWE, useCal, useImp, check_impstatus | ( (useCal | useHWE | useImp) & !ignore_impstatus))]
    header_test <- translate_header(header = colnames(dataset), standard = header_std, alternative = header_translations)
    if(any(duplicated(header_test$header_h))) { stop("cannot translate header - duplicate column names") }
    if(header_test$missing_N > 1) { stop("cannot translate header - missing columns") }
    if(header_test$missing_N == 1) {
      if(header_test$missing_h == "IMPUTED" & !check_impstatus & (ignore_impstatus | (!useImp | (!useHWE & !useCal))) ) {
        if(!ignore_impstatus) {
          if(!useImp) {
            dataset$IMPUTED <- 0
            print("Warning: no imputation-status specified - all SNPs set to genotyped")
          } else {
            dataset$IMPUTED <- 1
            print("Warning: no imputation-status specified - all SNPs set to imputed") }
          header_test$header_h <- c(header_test$header_h, "IMPUTED")
        }
      } else { stop(paste("cannot translate header - missing column:", paste(header_test$missing_h, collapse = ", "))) }
    }
    colnames(dataset) <- header_test$header_h
  }
  
  if(check_impstatus) {
    dataset$IMPUTED <- convert_impstatus(dataset$IMPUTED, T_strings, F_strings, NA_strings, use_log = FALSE)
    if(all(is.na(dataset$IMPUTED))) stop("Imputation status missing or untranslated")
  }
  
  plot_table <- data.frame(name = "",
                           FRQ = if(useFRQ) filter_FRQ else NA, cal = if(useCal) filter_cal else NA,
                           HWE = if(useHWE) filter_HWE else NA, imp = if(useImp) filter_imp else NA,
                           NA_FRQ = filter_NA_FRQ, NA_cal = filter_NA_cal, NA_HWE = filter_NA_HWE, NA_imp = filter_NA_imp,
                           stringsAsFactors = FALSE)
  plot_table$name <- if(length(save_name) < nrow(plot_table)) paste0(save_name, 1:nrow(plot_table)) else save_name
  
  
  for(ij in 1:nrow(plot_table)) {
    QQ_plot(dataset = dataset, save_name = plot_table$name[ij], save_dir = save_dir,
            filter_FRQ = if(!useFRQ) NULL else plot_table$FRQ[ij], filter_cal = if(!useCal) NULL else plot_table$cal[ij],
            filter_HWE = if(!useHWE) NULL else plot_table$HWE[ij], filter_imp = if(!useImp) NULL else plot_table$imp[ij],
            filter_NA_FRQ = plot_table$NA_FRQ[ij], filter_NA_cal = plot_table$NA_cal[ij],
            filter_NA_HWE = plot_table$NA_HWE[ij], filter_NA_imp = plot_table$NA_imp[ij],
            p_cutoff = p_cutoff, plot_QQ_bands = plot_QQ_bands,
            check_impstatus = FALSE, ignore_impstatus = ignore_impstatus, ... ) }
  return(invisible(NULL))
}
