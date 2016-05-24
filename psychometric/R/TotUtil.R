"TotUtil" <-
function(Rxy, Sy, MXg, COST, Nselected)
 {
 TU <- Nselected*Rxy*Sy*MXg - COST
 return(TU)
 }

