`radiateP` <-
function(MEC, SCALE=FALSE, col=col, TIT=FALSE)
{
  if(missing(SCALE)) { SCALE=FALSE }
  if(missing(col)) { col=heat.colors(20) }
  if(missing(TIT)) {  TIT=FALSE }

  updown = "Lower"
  if(MEC$UP==TRUE) { updown = " Upper" } else { updown = " Lower" }

  imageP(MEC$az1, MEC$dip1, MEC$rake1, SCALE=SCALE, UP=MEC$UP, col=col )

  net(add=TRUE)
  PlotPlanes(MEC)

  if(TIT==TRUE)
    {
   
  title(main=paste(" ", "P-wave",MEC$name, updown) ,
        sub=paste(sep="", "dip=", MEC$dip1, " strike=",
          MEC$az1, " rake=", format.default(MEC$rake1, digits=3) ))
}
}

