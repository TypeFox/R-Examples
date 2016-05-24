   .SFCcalc <- function(FUELTYPE,FFMC, BUI, PC, GFL){
     SFC   <- rep(-999,length(FFMC))
     SFC   <- ifelse(FUELTYPE=="C1",ifelse(FFMC>84,0.75 + 0.75*(1-exp(-0.23*(FFMC-84)))**0.5,0.75 - 0.75*(1-exp(-0.23*(84-FFMC)))**0.5),SFC)    # 9a & 9b
     SFC   <- ifelse(FUELTYPE=="C2"|FUELTYPE=="M3"|FUELTYPE=="M4",5.0 * (1 - exp(-0.0115 * BUI)),SFC)
     SFC   <- ifelse(FUELTYPE=="C3"|FUELTYPE=="C4",5.0 * (1 - exp(-0.0164*BUI))**2.24,SFC)
     SFC   <- ifelse(FUELTYPE=="C5"|FUELTYPE=="C6",5.0 * (1 - exp(-0.0149*BUI))**2.48,SFC)
     SFC   <- ifelse(FUELTYPE=="C7",ifelse(FFMC>70,2*(1-exp(-0.104*(FFMC-70))),0)+1.5*(1-exp(-0.0201*BUI)),SFC)
     SFC   <- ifelse(FUELTYPE=="D1",1.5 * (1 - exp(-0.0183 * BUI)),SFC)
     SFC   <- ifelse(FUELTYPE=="M1"|FUELTYPE=="M2",PC/100*(5.0*(1-exp(-0.0115*BUI)))+((100-PC)/100*(1.5*(1-exp(-0.0183*BUI)))),SFC)
     SFC   <- ifelse(FUELTYPE=="O1A"|FUELTYPE=="O1B",GFL,SFC)
     SFC   <- ifelse(FUELTYPE=="S1",4.0*(1-exp(-0.025*BUI))+4.0*(1-exp(-0.034*BUI)),SFC)                                               # 19,20,25
     SFC   <- ifelse(FUELTYPE=="S2",10.0*(1-exp(-0.013*BUI))+6.0*(1-exp(-0.060*BUI)),SFC)                                              # 19,20,25
     SFC   <- ifelse(FUELTYPE=="S3",12.0*(1-exp(-0.0166*BUI))+20.0*(1-exp(-0.0210*BUI)),SFC)                                           # 19,20,25
     SFC   <- ifelse(SFC<=0,0.000001,SFC)
     SFC}
