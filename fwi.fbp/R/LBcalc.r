  .LBcalc   <- function(FUELTYPE,WSV){
     LB    <- ifelse(FUELTYPE %in% c("O1A","O1B"),
       ifelse(WSV >= 1.0,1.1 * (WSV**0.464),                                                                                           # /* corrected from "+" to "*" in the errata 80 */
           1.0),                                                                                                                       # /* 81 */
           1.0 + 8.729 * (1-exp(-0.030*WSV))**(2.155))                                                                                 # /* 79 */
     LB}
