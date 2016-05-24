t_prog <-
function(N, m, cdiagn, creg, cdist)
 {
  d0   = 0.01
  dmax = 130
  V0 <- pi*(d0^3)/6
  Vmax <- pi*(dmax^3)/6# = 1150347
  repeat {
  vdiagn <- 1000*rlnorm(N, cdiagn[1], cdiagn[2])# we multiply by 1000 so as to get the volume instead of ml in cubic mm and calculate the diameter in mm
  vreg   <- 1000*rlnorm(N, cdist[1],  cdist[2])
  vdist  <- 1000*rlnorm(N, creg[1],   creg[2])
  if (V0<vdiagn & vdiagn<Vmax & V0<vdist & vdist<Vmax & V0<vreg & vreg<Vmax & vreg<vdist)
  break  }
  Treg   <- t_vol(vreg, m)
  Tdist  <- t_vol(vdist, m)
  Tdiagn <- t_vol(vdiagn, m)
  Ddiagn <- ((6/pi)*vdiagn)^(1/3)
  if (Tdiagn < Treg) {stage = "localized"} else if (Tdiagn < Tdist) {stage = "regional"} else stage = "distant"  
  return(list("Treg"=Treg, "Tdist"=Tdist, "Tdiagn"=Tdiagn, "Ddiagn"=Ddiagn, "stage" = stage))
 }
