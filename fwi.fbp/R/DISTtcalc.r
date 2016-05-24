.DISTtcalc <- function(FUELTYPE,ROSeq,HR,CFB){
  alpha <- ifelse(FUELTYPE %in% c("C1","O1A","O1B","S1","S2","S3","D1"),
                  0.115,                                                                                                                         # /* page 41 */
                  0.115 - 18.8 * (CFB**2.5) * exp(-8* CFB))                                                                                      # /* 72 */
  
  DISTt  <- ROSeq * (HR + exp(-alpha * HR)/alpha - 1/alpha)               # /* 71 */  
  
  DISTt}
