  .TFCcalc <-  function(FUELTYPE,CFL,CFB,SFC,PC,PDF,option="TFC"){
     CFC  <- CFL * CFB                                                                                                                 # /* 66  66a - 2009 */
     CFC  <- ifelse(FUELTYPE %in% c("M1","M2"),PC/100*CFC,CFC)                                                                         # /* 66b - 2009 */
     CFC  <- ifelse(FUELTYPE %in% c("M3","M4"),PDF/100*CFC,CFC)                                                                        # /* 66c - 2009 */
     TFC  <- SFC + CFC                                                                                                                 # /* 67 */
     if (option=="TFC"){TFC} else
     if (option=="CFC"){CFC}}

