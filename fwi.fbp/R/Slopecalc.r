   .Slopecalc <- function(FUELTYPE,FFMC,BUI,WS,WAZ,GS,SAZ,FMC,SFC,PC,PDF,CC,CBH,ISI,output="RAZ"){                                          # output options include: RAZ and WSV
     NoBUI   <- rep(-1,length(FFMC))
     SF      <- ifelse (GS >= 70,10,exp(3.533 * (GS/100)^1.2))                                                                         # /* 39 */
     ISZ     <- .ISIcalc(FFMC, 0)
     RSZ     <- .ROScalc(FUELTYPE,ISZ,BUI=NoBUI,FMC,SFC,PC,PDF,CC,CBH)
     RSF     <- RSZ * SF                                                                                                               # /* 40 */

     d       <- c("C1","C2","C3","C4","C5","C6","C7","D1","M1","M2","M3","M4","S1","S2","S3","O1A","O1B")
     a       <- c(90,110,110,110,30,30,45,30,0,0,120,100,75,40,55,190,250)
     b       <- c(0.0649,0.0282,0.0444,0.0293,0.0697,0.0800,0.0305,0.0232,0,0,0.0572,0.0404,0.0297,0.0438,0.0829,0.0310,0.0350)
     c0      <- c(4.5,1.5,3.0,1.5,4.0,3.0,2.0,1.6,0,0,1.4,1.48,1.3,1.7,3.2,1.4,1.7)                                                    # This is to avoid using "c" since it is a function in R.
     names(a)<-names(b)<-names(c0)<-d

     RSZ     <- rep(-99,length(FFMC))
     RSF_C2  <- rep(-99,length(FFMC))
     RSF_D1  <- rep(-99,length(FFMC))
     RSF_M3  <- rep(-99,length(FFMC))
     RSF_M4  <- rep(-99,length(FFMC))
     CF      <- rep(-99,length(FFMC))
     ISF     <- rep(-99,length(FFMC))
     ISF_C2  <- rep(-99,length(FFMC))
     ISF_D1  <- rep(-99,length(FFMC))
     ISF_M3  <- rep(-99,length(FFMC))
     ISF_M4  <- rep(-99,length(FFMC))
     # FUELTYPE c("C1","C2","C3","C4","C5","C6","C7","D1","S1","S2","S3")
     options(warn=-1)
     ISF     <- ifelse(FUELTYPE %in% c("C1","C2","C3","C4","C5","C6","C7","D1","S1","S2","S3"),
            ifelse((1 - (RSF/a[FUELTYPE])**(1/c0[FUELTYPE])) >= 0.01,
                  log(1 - (RSF/a[FUELTYPE])**(1/c0[FUELTYPE]))/(-b[FUELTYPE]),                                                         # /* 41b - 2009 */
                  log(0.01)/(-b[FUELTYPE])),ISF)                                                                                       # /* 41a - 2009 */

     options(warn=1)
   # FUELTYPE c("M1","M2")
     RSZ    <- ifelse(FUELTYPE %in% c("M1","M2"),.ROScalc(rep("C2",length(ISZ)), ISZ, BUI=NoBUI, FMC, SFC, PC, PDF, CC, CBH),RSZ)
     RSF_C2 <- ifelse(FUELTYPE %in% c("M1","M2"),RSZ * SF,RSF_C2)                                                                      # /* 40 */
     RSZ    <- ifelse(FUELTYPE %in% c("M1","M2"),.ROScalc(rep("D1",length(ISZ)), ISZ, BUI=NoBUI, FMC, SFC, PC, PDF, CC, CBH),RSZ)
     RSF_D1 <- ifelse(FUELTYPE %in% c("M1","M2"),RSZ * SF,RSF_D1)                                                                      # /* 40 */
     RSF0   <- 1 - (RSF_C2/a[["C2"]])^(1/c0[["C2"]])
     ISF_C2 <- ifelse(FUELTYPE %in% c("M1","M2")&RSF0 >= 0.01,
           log(1 - (RSF_C2/a[["C2"]])**(1/c0[["C2"]]))/(-b[["C2"]]),ISF_C2)                                                            # /* 41a - 2009 */
     ISF_C2 <- ifelse(FUELTYPE %in% c("M1","M2")&RSF0 < 0.01,
           log(0.01)/(-b[["C2"]]),                                                                                                     # /* 41b - 2009 */
           ISF_C2)
     RSF0   <- 1 - (RSF_D1/a[["D1"]])^(1/c0[["D1"]])
     ISF_D1 <- ifelse(FUELTYPE %in% c("M1","M2")&RSF0 >= 0.01,
           log(1 - (RSF_D1/a[["D1"]])**(1/c0[["D1"]]))/(-b[["D1"]]),ISF_D1)                                                            # /* 41a - 2009 */
     ISF_D1 <- ifelse(FUELTYPE %in% c("M1","M2")&RSF0 < 0.01,
           log(0.01)/(-b[["D1"]]),                                                                                                     # /* 41b - 2009 */
           ISF_D1)
     ISF    <- ifelse(FUELTYPE %in% c("M1","M2"),PC/100*ISF_C2+(1-PC/100)*ISF_D1,ISF)                                                  # /* 42a - 2009 */
   # FUELTYPE c("M3")
     PDF100 <- rep(100,length(ISI))
     RSZ    <- ifelse(FUELTYPE %in% c("M3"),.ROScalc(rep("M3",length(FMC)),ISI=ISZ,BUI=NoBUI,FMC,SFC,PC,PDF100,CC,CBH),RSZ)
     RSF_M3 <- ifelse(FUELTYPE %in% c("M3"),RSZ*SF,RSF_M3)                                                                             # /* 40 */
     RSZ    <- ifelse(FUELTYPE %in% c("M3"),.ROScalc(rep("D1",length(ISZ)),ISZ,BUI=NoBUI,FMC,SFC,PC,PDF100,CC,CBH),RSZ)
     RSF_D1 <- ifelse(FUELTYPE %in% c("M3"),RSZ*SF,RSF_D1)
     RSF0   <- 1 - (RSF_M3/a[["M3"]])^(1/c0[["M3"]])                                                                                   # /* 40 */
     ISF_M3 <- ifelse(FUELTYPE %in% c("M3")&RSF0 >= 0.01,
                log(1 - (RSF_M3/a[["M3"]])**(1/c0[["M3"]]))/(-b[["M3"]]),ISF_M3)                                                       # /* 41a - 2009 */
     ISF_M3 <- ifelse(FUELTYPE %in% c("M3")&RSF0 < 0.01,
                log(0.01)/(-b[["M3"]]),                                                                                                # /* 41b - 2009 */
                ISF_M3)
     RSF0   <- 1 - (RSF_D1/a[["D1"]])^(1/c0[["D1"]])                                                                                   # /* 40 */
     ISF_D1 <- ifelse(FUELTYPE %in% c("M3")&RSF0 >= 0.01,
                log(1 - (RSF_D1/a[["D1"]])**(1/c0[["D1"]]))/(-b[["D1"]]),ISF_D1)                                                       # /* 41a - 2009 */
     ISF_D1 <- ifelse(FUELTYPE %in% c("M3")&RSF0 < 0.01,
                log(0.01)/(-b[["D1"]]),
                ISF_D1)
     ISF    <- ifelse(FUELTYPE %in% c("M3"),PDF/100*ISF_M3 + (1-PDF/100)*ISF_D1,ISF)                                                   # /* 42b - 2009 */
   # FUELTYPE c("M4")
     RSZ    <- ifelse(FUELTYPE %in% c("M4"),.ROScalc(rep("M4",length(FMC)),ISI=ISZ,BUI=NoBUI,FMC,SFC,PC,PDF100,CC,CBH),RSZ)
     RSF_M4 <- ifelse(FUELTYPE %in% c("M4"),RSZ * SF,RSF_M4)
     RSZ    <- ifelse(FUELTYPE %in% c("M4"),.ROScalc(rep("D1",length(ISZ)),ISZ,BUI=NoBUI,FMC,SFC,PC,PDF100,CC,CBH),RSZ)
     RSF_D1 <- ifelse(FUELTYPE %in% c("M4"),RSZ * SF,RSF_D1)                                                                           # /* 40 */
     RSF0   <- 1 - (RSF_M4/a[["M4"]])^(1/c0[["M4"]])                                                                                   # /* 40 */
     ISF_M4 <- ifelse(FUELTYPE %in% c("M4")&RSF0 >= 0.01,
                log(1 - (RSF_M4/a[["M4"]])**(1/c0[["M4"]]))/(-b[["M4"]]),ISF_M4)                                                                    # /* 41a - 2009 */
     ISF_M4 <- ifelse(FUELTYPE %in% c("M4")&RSF0 < 0.01,
                log(0.01)/(-b[["M4"]]),
                ISF_M4)                                                                                                                # /* 41b - 2009 */
     RSF0   <- 1 - (RSF_D1/a[["D1"]])^(1/c0[["D1"]])                                                                                   # /* 40 */
     ISF_D1 <- ifelse(FUELTYPE %in% c("M4")&RSF0 >= 0.01,
                log(1 - (RSF_D1/a[["D1"]])**(1/c0[["D1"]]))/(-b[["D1"]]),ISF_D1)                                                                    # /* 41a - 2009 */
     ISF_D1 <- ifelse(FUELTYPE %in% c("M4")&RSF0 < 0.01,
                log(0.01)/(-b[["D1"]]),
                ISF_D1)                                                                                                                # /* 41b - 2009 */
     ISF    <- ifelse(FUELTYPE %in% c("M4"),PDF/100*ISF_M4 + (1-PDF/100.)*ISF_D1,ISF)                                                  # /* 42c - 2009 */
   # FUELTYPE c("O1A","O1B")
     CF     <- ifelse(FUELTYPE %in% c("O1A","O1B"),ifelse(CC < 58.8,0.005*(exp(0.061*CC)-1),                                           # /* 35a - 2009 */
            0.176 + 0.02*(CC-58.8)),                                                                                                   # /* 35b - 2009 */
            CF)

     ISF    <- ifelse(FUELTYPE %in% c("O1A","O1B"),ifelse((1 - (RSF/(CF*a[FUELTYPE]))**(1/c0[FUELTYPE])) >= 0.01,
             log(1 - (RSF/(CF*a[FUELTYPE]))**(1/c0[FUELTYPE]))/(-b[FUELTYPE]),                                                         # /* 43a - 2009 */
             log(0.01)/(-b[FUELTYPE])),                                                                                                # /* 43b - 2009 */
             ISF)

     m      <- 147.2*(101-FFMC)/(59.5+FFMC)                                                                                            # /* 46 */
     fF     <- 91.9*exp(-.1386*m)*(1+(m**5.31)/4.93e7)                                                                                 # /* 45 */

     #//   WSE <- log(ISF/(0.208 * fF))/0.05039                                                                                        # /* 44 */
     WSE    <- 1/0.05039*log(ISF/(0.208 * fF))
     WSE    <- ifelse(WSE>40&ISF < (0.999*2.496*fF),28-(1/0.0818*log(1-ISF/(2.496*fF))),WSE)                                           # /* 44e/44b,/44c- 2009 */
     WSE    <- ifelse(WSE>40&ISF >= (0.999*2.496*fF),112.45,WSE)                                                                       # /* 44e/44b,/44c- 2009 */

     WSX    <- WS*sin(WAZ) + WSE*sin(SAZ)                                                                                              # /* 47 */

     WSY    <- WS*cos(WAZ) + WSE*cos(SAZ)                                                                                              # /* 48 */
     WSV    <- sqrt(WSX*WSX + WSY*WSY)                                                                                                 # /* 49 */


     RAZ    <- acos(WSY/(WSV))                                                           # /* in radians */                            # /* 50 */
     RAZ    <- ifelse(WSX < 0,2*pi - RAZ,RAZ)                                                                                          # /* 51 */
     if (output=="RAZ"){
     RAZ} else
     if (output=="WAZ"){
     WAZ} else
     if (output=="WSV"){
     WSV}}
