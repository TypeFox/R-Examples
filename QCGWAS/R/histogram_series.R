histogram_series <-
function(
  dataset, data_col = 1, save_name = paste0("dataset_F", 1:nrow(plot_table)), save_dir = getwd(), export_outliers = FALSE,
  filter_FRQ = NULL, filter_cal = NULL, filter_HWE = NULL, filter_imp = NULL,
  filter_NA = TRUE, filter_NA_FRQ = filter_NA, filter_NA_cal = filter_NA, filter_NA_HWE = filter_NA, filter_NA_imp = filter_NA,
  breaks = "Sturges", 
  header_translations, ignore_impstatus = FALSE, check_impstatus = FALSE, T_strings = c("1", "TRUE", "yes", "YES", "y", "Y"), F_strings = c("0", "FALSE", "no", "NO", "n", "N"), NA_strings = c(NA, "NA", ".", "-"), ...
) {
  
  if(length(data_col) != 1L) { stop("Multiple columns specified") }
  
  useFRQ <- if(is.null(filter_FRQ)) FALSE else !all(is.na(filter_FRQ) & !filter_NA_FRQ)
  useCal <- if(is.null(filter_cal)) FALSE else !all(is.na(filter_cal) & !filter_NA_cal)
  useHWE <- if(is.null(filter_HWE)) FALSE else !all(is.na(filter_HWE) & !filter_NA_HWE)
  useImp <- if(is.null(filter_imp)) FALSE else !all(is.na(filter_imp) & !filter_NA_imp)
  use_impstatus <- if(ignore_impstatus) check_impstatus else check_impstatus | useCal | useHWE | useImp

  if(is.vector(dataset)) {
    if(ignore_impstatus | check_impstatus) stop("Vector dataset - cannot ignore or check imp-status")
    if(is.numeric(data_col) & data_col != 1L) stop("Vector dataset - only one column")
    if(sum(useFRQ, useCal, useHWE, useImp) == 1L) {
      dataset <- data.frame(EFFECT = dataset)
      if(is.character(data_col)) {
        colnames(dataset) <- data_col
      } else {
        if(useFRQ) { colnames(dataset) <- "EFF_ALL_FREQ" }
        if(useCal) { colnames(dataset) <- "CALLRATE" }
        if(useHWE) { colnames(dataset) <- "HWE_PVAL" }
        if(useImp) { colnames(dataset) <- "IMP_QUALITY" }
      }
    } else { stop("Insufficient data to apply filters: dataset is single column!") } }
  
  if(is.character(data_col)) {
    graph_name <- data_col
    data_col <- which(colnames(dataset) == data_col)
  } else {
    graph_name <- colnames(dataset)[data_col]
  }
  if(!graph_name %in% colnames(dataset)) stop("Invalid column specified")
  
  if(missing(header_translations)) {
    if(use_impstatus & !"IMPUTED" %in% colnames(dataset)) {
      if(check_impstatus) stop("cannot check imputation status - no such column in dataset")
      if(!ignore_impstatus) {
        if(useImp & (useHWE | useCal)) stop("Missing imputation status")
        if( useHWE &  useCal & !useImp) {
          dataset$IMPUTED <- 0L
          print("Warning: no imputation-status specified - all SNPs set to genotyped") }
        if(!useHWE & !useCal &  useImp) {
          dataset$IMPUTED <- 1L
          print("Warning: no imputation-status specified - all SNPs set to imputed") }
    } }
  } else {
    header_std <- c("EFF_ALL_FREQ", "HWE_PVAL", "CALLRATE", "IMP_QUALITY", "IMPUTED")[c(useFRQ, useHWE, useCal, useImp, use_impstatus)]
    header_test <- translate_header(header = colnames(dataset), standard = header_std, alternative = header_translations)
    if(any(duplicated(header_test$header_h))) stop("cannot translate header - duplicate column names")
    if(header_test$missing_N > 1L) stop("cannot translate header - missing columns")
    if(header_test$missing_N == 1L) {
      if(header_test$missing_h == "IMPUTED" & !check_impstatus & (ignore_impstatus | (!useImp | (!useHWE & !useCal)))) {
        if(!ignore_impstatus) {
          if(!useImp) {
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
  
  plot_table <- data.frame(name = "temp",
                           FRQ = if(useFRQ) filter_FRQ else NA, cal = if(useCal) filter_cal else NA,
                           HWE = if(useHWE) filter_HWE else NA, imp = if(useImp) filter_imp else NA,
                           NA_FRQ = filter_NA_FRQ, NA_cal = filter_NA_cal, NA_HWE = filter_NA_HWE, NA_imp = filter_NA_imp,
                           stringsAsFactors = FALSE)
  if(length(save_name) < nrow(plot_table)) {
    plot_table$name <- paste0(save_name, 1:nrow(plot_table))
  } else {  plot_table$name <- save_name }
  
  
  for(ij in 1:nrow(plot_table)) {
    QC_histogram(dataset = dataset, data_col = data_col, save_name = plot_table$name[ij], save_dir = save_dir, export_outliers = export_outliers,
                 filter_FRQ = if(useFRQ) plot_table$FRQ[ij] else NULL, filter_cal = if(useCal) plot_table$cal[ij] else NULL,
                 filter_HWE = if(useHWE) plot_table$HWE[ij] else NULL, filter_imp = if(useImp) plot_table$imp[ij] else NULL,
                 filter_NA_FRQ = plot_table$NA_FRQ[ij], filter_NA_cal = plot_table$NA_cal[ij],
                 filter_NA_HWE = plot_table$NA_HWE[ij], filter_NA_imp = plot_table$NA_imp[ij],
                 breaks = breaks, graph_name = graph_name, check_impstatus = FALSE, ignore_impstatus = ignore_impstatus, ... ) }
  return(invisible())
}
