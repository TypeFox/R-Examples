"MargUtil" <-
function(Rxy, Sy, MXg, COST, Nselected)
 {
 MU <- Rxy*Sy*MXg - COST/Nselected
 return(MU)
 }

