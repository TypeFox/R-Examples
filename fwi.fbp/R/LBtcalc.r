  .LBtcalc  <- function(FUELTYPE,LB,HR,CFB){                                                                                            # /* this is a new subroutine to account for the discussion in 3.6.1 */
    alpha  <- ifelse(FUELTYPE %in% c("C1","O1A","O1B","S1","S2","S3","D1"),
        0.115,                                                                                                                         # /* page 41 */
        0.115 - 18.8 * (CFB**2.5) * exp(-8* CFB))                                                                                      # /* 72 */
  	LBt    <- (LB -1) * (1 - exp(-alpha*HR)) + 1						                                                                           # /* 81 - 2009 */
    LBt}
