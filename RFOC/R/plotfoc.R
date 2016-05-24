`plotfoc` <-
function(MEC)
{
par(mfrow=c(3,1))

radiateP(MEC )

if(length(MEC$PTS)>0)
  {
    focpoint(MEC$PTS$gaz, MEC$PTS$angP,  lab=MEC$PTS$name, UP=MEC$UP)
  }
# print(paste(" ", "P=", MEC$P$az, MEC$P$dip, "T=", MEC$T$az, MEC$T$dip))
focpoint(MEC$P$az, MEC$P$dip,  lab="P", UP=MEC$UP)
focpoint(MEC$T$az, MEC$T$dip,  lab="T", UP=MEC$UP)

# dev.set(which = dev.next())
radiateSV(MEC )
if(length(MEC$PTS)>0)
  {
focpoint(MEC$PTS$gaz, MEC$PTS$angS,  lab=MEC$PTS$name, UP=MEC$UP)
}
# dev.set(which = dev.next())
radiateSH(MEC )
if(length(MEC$PTS)>0)
  {
focpoint(MEC$PTS$gaz, MEC$PTS$angS,  lab=MEC$PTS$name, UP=MEC$UP)
}
}

