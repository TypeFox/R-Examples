HQ_filter <-
function(data, ignore_impstatus = FALSE,
                      FRQ_val = NULL, HWE_val = NULL, cal_val = NULL, imp_val = NULL,
                      filter_NA = TRUE, FRQ_NA = filter_NA, HWE_NA = filter_NA, cal_NA = filter_NA, imp_NA = filter_NA) {
  goodOnes <- !logical(length = nrow(data) )
  if(!is.null(FRQ_val)) {
    if(is.na(FRQ_val)) { if(FRQ_NA) { goodOnes <- goodOnes & !is.na(data$EFF_ALL_FREQ) }
    } else {             if(FRQ_NA) { goodOnes <- ifelse( is.na(data$EFF_ALL_FREQ) |  data$EFF_ALL_FREQ < FRQ_val | data$EFF_ALL_FREQ > 1 - FRQ_val , FALSE, goodOnes)
    } else {                          goodOnes <- ifelse(!is.na(data$EFF_ALL_FREQ) & (data$EFF_ALL_FREQ < FRQ_val | data$EFF_ALL_FREQ > 1 - FRQ_val), FALSE, goodOnes) } } }
  
  if(ignore_impstatus) {
    
    if(!is.null(HWE_val)) {
      if(is.na(HWE_val)) { if(HWE_NA) { goodOnes <- goodOnes & !is.na(data$HWE_PVAL) }
      } else {             if(HWE_NA) { goodOnes <- ifelse((is.na(data$HWE_PVAL) | data$HWE_PVAL < HWE_val), FALSE, goodOnes)
      } else {                          goodOnes <- ifelse(!is.na(data$HWE_PVAL) & data$HWE_PVAL < HWE_val , FALSE, goodOnes) } } }
    if(!is.null(cal_val)) {
      if(is.na(cal_val)) { if(cal_NA) { goodOnes <- goodOnes & !is.na(data$CALLRATE) }
      } else {             if(cal_NA) { goodOnes <- ifelse((is.na(data$CALLRATE) | data$CALLRATE < cal_val), FALSE, goodOnes)
      } else {                          goodOnes <- ifelse(!is.na(data$CALLRATE) & data$CALLRATE < cal_val , FALSE, goodOnes) } } }
    if(!is.null(imp_val)) {
      if(is.na(imp_val)) { if(imp_NA) { goodOnes <- goodOnes & !is.na(data$IMP_QUALITY) }
      } else {             if(imp_NA) { goodOnes <- ifelse((is.na(data$IMP_QUALITY) | data$IMP_QUALITY < imp_val), FALSE, goodOnes)
      } else {                          goodOnes <- ifelse(!is.na(data$IMP_QUALITY) & data$IMP_QUALITY < imp_val , FALSE, goodOnes) } } }
    
  } else {
    
    if(!is.null(HWE_val)) {
      if(is.na(HWE_val)) { if(HWE_NA) { goodOnes <- ifelse( is.na(data$HWE_PVAL) & data$IMPUTED == 0 & !is.na(data$IMPUTED), FALSE, goodOnes) }
      } else {             if(HWE_NA) { goodOnes <- ifelse((is.na(data$HWE_PVAL) | data$HWE_PVAL < HWE_val) & data$IMPUTED == 0 & !is.na(data$IMPUTED), FALSE, goodOnes)
      } else {                          goodOnes <- ifelse(!is.na(data$HWE_PVAL) & data$HWE_PVAL < HWE_val  & data$IMPUTED == 0 & !is.na(data$IMPUTED), FALSE, goodOnes) } } }
    if(!is.null(cal_val)) {
      if(is.na(cal_val)) { if(cal_NA) { goodOnes <- ifelse( is.na(data$CALLRATE) & data$IMPUTED == 0 & !is.na(data$IMPUTED), FALSE, goodOnes) }
      } else {             if(cal_NA) { goodOnes <- ifelse((is.na(data$CALLRATE) | data$CALLRATE < cal_val) & data$IMPUTED == 0 & !is.na(data$IMPUTED), FALSE, goodOnes)
      } else {                          goodOnes <- ifelse(!is.na(data$CALLRATE) & data$CALLRATE < cal_val  & data$IMPUTED == 0 & !is.na(data$IMPUTED), FALSE, goodOnes) } } }
    if(!is.null(imp_val)) {
      if(is.na(imp_val)) { if(imp_NA) { goodOnes <- ifelse( is.na(data$IMP_QUALITY) & data$IMPUTED == 1 & !is.na(data$IMPUTED), FALSE, goodOnes) }
      } else {             if(imp_NA) { goodOnes <- ifelse((is.na(data$IMP_QUALITY) | data$IMP_QUALITY < imp_val) & data$IMPUTED == 1 & !is.na(data$IMPUTED), FALSE, goodOnes)
      } else {                          goodOnes <- ifelse(!is.na(data$IMP_QUALITY) & data$IMP_QUALITY < imp_val  & data$IMPUTED == 1 & !is.na(data$IMPUTED), FALSE, goodOnes) } } }
  }
  return(goodOnes)
}
