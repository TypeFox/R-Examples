  .C6calc<-function(FUELTYPE,ISI,BUI,FMC,SFC,CBH,ROS,CFB,RSC,option="CFB"){                                                             #options include(&ROS, &CFB, &RSC, &RSI)
     FMEavg <- 0.778                                                                                                                   # /* page 37 */
     tt     <- 1500 - 2.75 * FMC                                                                                                       # /* 59 */
     H      <- 460 + 25.9 * FMC                                                                                                        # /* 60 */
     FME    <- ((1.5 - 0.00275 * FMC)**4.)/(460 + 25.9*FMC) * 1000                                                                     # /* 61 */
     RSI    <- 30 * (1 - exp(-0.08 * ISI))**3.0                                                                                        # /* 62 */
     RSS    <- RSI * .BEcalc(FUELTYPE, BUI)                                                                                             # /* 63 */
     RSC    <- 60 * (1 - exp(-0.0497*ISI)) * FME/FMEavg                                                                                # /* 64 */

     CFB    <- ifelse(RSC > RSS,.CFBcalc(FUELTYPE, FMC, SFC, RSS, CBH),0)
     ROS    <- ifelse(RSC > RSS,RSS + (CFB)*(RSC-RSS),RSS)                                                                             # /* 65 */

     if (option=="CFB"){CFB}else
     if (option=="ROS"){ROS}else
     if (option=="RSC"){RSC}else
     if (option=="RSI"){RSI}}
