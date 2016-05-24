AT.SPC.spectrum.at.depth.g.cm2 <- function(spc, depth.g.cm2, interpolate = TRUE)
{
    depth.step.of.spc      <- unique(spc$depth.step)
    depth.g.cm2.of.spc     <- unique(spc$depth.g.cm2)
    depth.step.interp      <- approx( x    = depth.g.cm2.of.spc,
                  	              y    = depth.step.of.spc,
                    	              xout = depth.g.cm2,
                    	              rule = 2)$y
    
    depth.step.int         <- floor(depth.step.interp)
    depth.step.frac        <- depth.step.interp - depth.step.int
    spc.before             <- spc[spc$depth.step == depth.step.int,] 
    if(interpolate){
	 spc.after              <- spc[spc$depth.step == (depth.step.int+1),]  

	 spc.interp             <- spc.before
	 spc.interp$depth.step  <- depth.step.interp
	 spc.interp$depth.g.cm2 <- depth.g.cm2
	 spc.interp$fluence.cm2 <- (1 - depth.step.frac) * spc.before$dN.dE.per.MeV.u.per.primary * spc.before$DE.MeV.u+ depth.step.frac * spc.after$dN.dE.per.MeV.u.per.primary * spc.after$DE.MeV.u

	 return(spc.interp)
    }else{
	 return(spc.before)
    }
}
