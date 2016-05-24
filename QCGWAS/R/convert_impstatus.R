convert_impstatus <-
function(impstatus, T_strings = c("1", "TRUE", "T"), F_strings = c("0", "FALSE", "F"), NA_strings = c(NA, "NA", ".", "-"),
                              use_log = FALSE, ...) {
  if(any(duplicated(c(T_strings, F_strings, NA_strings)))) stop("duplicate strings in string-arguments")
  temp_imp <- ifelse(impstatus %in% T_strings, 1L, 2L)
  temp_imp <- ifelse(impstatus %in% F_strings, 0L, temp_imp)
  if(length(NA_strings) > 0L) temp_imp <- ifelse(impstatus %in% NA_strings, NA, temp_imp)
  if(any(temp_imp == 2L, na.rm = TRUE)) {
    if(use_log) save_log(phaseL = 2L, checkL = "data integrity", typeL = "imputation status", SNPL = sum(temp_imp == 2L, na.rm = TRUE), actionL = "set to NA", noteL = "Cannot parse imputation status: unidentified strings", ...)
    print(" - - warning: untranslatable values in imputation-status column", quote = FALSE)
    temp_imp <- ifelse(temp_imp == 2L, NA, temp_imp)
  }
  return(temp_imp)
}
