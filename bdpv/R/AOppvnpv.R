AOppvnpv <-
function(se, sp){
oddsppv <- sqrt((1-se)*(1-sp)/(se*sp))
oddsnpv <- sqrt((se*sp)/((1-se)*(1-sp)))
Pppv<-oddsppv/(1+oddsppv)
Pnpv<-oddsnpv/(1+oddsnpv)
return(c(Pppv=Pppv, Pnpv=Pnpv))
}

