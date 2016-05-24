  .BEcalc<-function(FUELTYPE,BUI){
     d    <- c("C1","C2","C3","C4","C5","C6","C7","D1","M1","M2","M3","M4","S1","S2","S3","O1A","O1B")
     BUIo <- c(72,64,62,66,56,62,106,32,50,50,50,50,38,63,31,01,01)
     Q    <- c(0.9,0.7,0.75,0.8,0.8,0.8,0.85,0.9,0.8,0.8,0.8,0.8,0.75,0.75,0.75,1.0,1.0)                                               # Use Q instead of q
     names(BUIo)<-names(Q)<-d
     BE   <- ifelse(BUI > 0 & BUIo[FUELTYPE] > 0,
     exp(50*log(Q[FUELTYPE])*(1/BUI - 1/BUIo[FUELTYPE])),                                                                              # /* 54 */
     1)
     as.numeric(BE)}
