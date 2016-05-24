  .BROScalc <- function(FUELTYPE,FFMC,BUI,WSV,FMC,SFC,PC,PDF,CC,CBH){
     m     <- 147.2*(101-FFMC)/(59.5+FFMC)                                                                                             # /* 46 */
     fF    <- 91.9*exp(-0.1386*m)*(1.+(m**5.31)/4.93e7)                                                                                # /* 45 */
     BfW   <- exp(-0.05039*WSV)                                                                                                        # /* 75 */
     BISI  <- 0.208*BfW*fF                                                                                                             # /* 76 */
    # /* Note the BUI effect is captured in ROScalc */
     BROS  <-  .ROScalc(FUELTYPE,BISI,BUI,FMC,SFC,PC,PDF,CC,CBH)                                                                        # /* 77 */
     BROS}
