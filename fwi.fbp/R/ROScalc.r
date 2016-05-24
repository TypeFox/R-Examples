.ROScalc     <- function(FUELTYPE, ISI, BUI, FMC, SFC, PC, PDF, CC, CBH){
     NoBUI    <- rep(-1,length(ISI))
     d        <- c("C1","C2","C3","C4","C5","C6","C7","D1","M1","M2","M3","M4","S1","S2","S3","O1A","O1B")
     a        <- c(90,110,110,110,30,30,45,30,0,0,120,100,75,40,55,190,250)
     b        <- c(0.0649,0.0282,0.0444,0.0293,0.0697,0.0800,0.0305,0.0232,0,0,0.0572,0.0404,0.0297,0.0438,0.0829,0.0310,0.0350)
     c0       <- c(4.5,1.5,3.0,1.5,4.0,3.0,2.0,1.6,0,0,1.4,1.48,1.3,1.7,3.2,1.4,1.7)                                                   # This is to avoid using "c" since it is a function in R.
     names(a) <-names(b)<-names(c0)<-d

     RSI      <- rep(-1,length(ISI))
     RSI      <- ifelse(FUELTYPE %in% c("C1","C2","C3","C4","C5","C7","D1","S1","S2","S3"),
                 as.numeric(a[FUELTYPE] * (1 - exp(-b[FUELTYPE] * ISI))**c0[FUELTYPE]),RSI)                                            #  /* 26 */
     RSI      <- ifelse(FUELTYPE %in% c("M1"), PC/100 * .ROScalc(rep("C2",length(ISI)),ISI,NoBUI,FMC,SFC,PC,PDF,CC,CBH)
                + (100-PC)/100 *.ROScalc(rep("D1",length(ISI)),ISI,NoBUI,FMC,SFC,PC,PDF,CC, CBH),RSI)                                   # /* 27 */
     RSI      <- ifelse(FUELTYPE %in% c("M2"), PC/100 * .ROScalc(rep("C2",length(ISI)),ISI,NoBUI,FMC,SFC,PC,PDF,CC,CBH)
                + 0.2*(100-PC)/100 *.ROScalc(rep("D1",length(ISI)),ISI,NoBUI,FMC,SFC,PC,PDF,CC, CBH),RSI)                               # /* 27 */
     RSI_m3   <- rep(-99,length(ISI))
     RSI_m3   <- ifelse(FUELTYPE %in% c("M3"), as.numeric(a[["M3"]] * ((1 - exp(-b[["M3"]] * ISI))**c0[["M3"]])),RSI_m3)               # /* 30 - 2009 */
     RSI      <- ifelse(FUELTYPE %in% c("M3"),PDF/100* RSI_m3 +
                (1-PDF/100)* .ROScalc(rep("D1",length(ISI)),ISI,NoBUI,FMC,SFC,PC,PDF,CC,CBH),RSI)                                       # /* 29 - 2009 */
     RSI_m4   <- rep(-99,length(ISI))
     RSI_m4   <- ifelse(FUELTYPE %in% c("M4"),as.numeric(a[["M4"]] * ((1 - exp(-b[["M4"]] * ISI))**c0[["M4"]])),RSI_m4)                # /* 30 - 2009 */
     RSI      <- ifelse(FUELTYPE %in% c("M4"),PDF/100* RSI_m4 +
                0.2*(1-PDF/100)* .ROScalc(rep("D1",length(ISI)),ISI,NoBUI,FMC,SFC,PC,PDF,CC,CBH),RSI)                                   # /* 31 - 2009 */
     CF       <- rep(-99,length(ISI))
     CF       <- ifelse(FUELTYPE %in% c("O1A","O1B"),ifelse(CC < 58.8,0.005*(exp(0.061*CC)-1),0.176 + 0.02*(CC-58.8)),CF)              # /* 35b - 2009 */
     RSI      <- ifelse(FUELTYPE %in% c("O1A","O1B"),a[FUELTYPE] * ((1 - exp(-b[FUELTYPE] * ISI))**c0[FUELTYPE])* CF ,RSI)             # /* 36 */
#     RSI      <- ifelse(FUELTYPE %in% c("C6"),.C6calc(FUELTYPE,ISI,BUI,FMC,SFC,CBH,option="RSI"),RSI)
     ROS      <- ifelse(FUELTYPE %in% c("C6"),.C6calc(FUELTYPE,ISI,BUI,FMC,SFC,CBH,option="ROS"),.BEcalc(FUELTYPE,BUI)*RSI )
#     ROS      <- BEcalc(FUELTYPE,BUI)*RSI
     ROS      <- ifelse(ROS <= 0,0.000001,ROS)
     ROS
}
