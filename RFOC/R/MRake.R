`MRake` <-
function(M)
{
############  prepare a MEC structure for a focal mechanism later used in other routines
#####   print(M)
  
ang2 = GetRakeSense(M$uaz, M$ud, M$vaz, M$vd, M$paz, M$pd, M$taz, M$td)

MEC = GetRake(M$az1-90, M$d1,   M$az2-90,  M$d2, ang2)

MEC$P = list(az=M$paz, dip=M$pd)
MEC$T = list(az=M$taz, dip=M$td)

MEC$U = list(az=M$uaz, dip=M$ud)
MEC$V = list(az=M$vaz, dip=M$vd)

MEC$F = list(az=M$az1, dip=M$d1)
MEC$G = list(az=M$az2, dip=M$d2)

MEC$sense = ang2

MEC$M = M

return(MEC)

}

