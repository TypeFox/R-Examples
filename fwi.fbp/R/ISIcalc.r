  .ISIcalc <- function(FFMC,WSV){
     m    <- 147.2*(101.-FFMC)/(59.5+FFMC)                                                                                             # /* 46 */
     fF   <- 91.9*exp(-0.1386*m)*(1.+m**5.31/4.93e7)                                                                                   # /* 45 */
     fW   <- ifelse(WSV < 40,exp(0.05039*WSV),12 * (1-exp(-0.0818 * (WSV-28))))                                                        # /* 53 */
     ISI  <- 0.208*fW*fF                                                                                                               # /* 52 */
     ISI}
