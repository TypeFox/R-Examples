  .CFBcalc <- function(FUELTYPE,FMC,SFC,ROS,CBH,option="CFB"){
     CFB  <- 0
     CSI  <- 0.001 * (CBH**1.5) * (460 + 25.9*FMC)**1.5                                                                                #/* 56 */
     RSO  <- CSI/(300*SFC)                                                                                                             # /* 57 */
     CFB  <- ifelse(ROS > RSO, 1 - exp(-0.23*(ROS - RSO)),CFB)                                                                         # /* 58 */
     if (option=="CFB"){CFB} else
     if (option=="CSI"){CSI} else
     if (option=="RSO"){RSO}}
