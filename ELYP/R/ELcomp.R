ELcomp <- function(Haz, Sur, gam) {
####input is in the order of y. No need of d?? (since Haz>0 same as d>0)
#### output is LOG Emplik (not Emplik)
if ( any(Sur >1) || any(Sur <0) ) 
  {warning("Sur out of [0,1]"); Sur <- pmax(0, pmin(1, Sur) ) }
Rti <- (1-Sur)/Sur
valu <- Haz/(gam[,2]+(gam[,1]-gam[,2])*Sur)
valu <- valu[valu > 0]
part1 <- sum( log(valu) )
part2 <- sum( (log( gam[,1]/(gam[,1]+gam[,2]*Rti) ))/gam[,2] )
return(part1+part2)
}

