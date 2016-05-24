`SDRfoc` <-
function(s,d,r, u=FALSE, ALIM=c(-1,-1, +1, +1), PLOT=TRUE)
{
  #############   these are the strike dip and rake from the Harvard CMT catalogues
  if(missing(u)) u = FALSE
    if(missing(ALIM))  {ALIM=c(-1,-1, +1, +1) }
    if(missing(PLOT))  {PLOT=TRUE }

  mc = CONVERTSDR(s,d,r )
  MEC = MRake(mc$M)
  MEC$UP = u
  MEC$icol =  foc.icolor(MEC$rake1)
  MEC$ileg =  focleg(MEC$icol)
  MEC$fcol =   foc.color(MEC$icol)
  MEC$CNVRG = NA
  MEC$LIM = ALIM
  
  if(PLOT)
    {
      ##  one()

      
      Beachfoc(MEC, fcol=MEC$fcol, fcolback="white")
      net(add=TRUE)
      
      ## addPT(MEC)
      addmecpoints(MEC)
      PlotPlanes(MEC, col1="blue", col2=grey(.6) )

      
      
      tit1 = paste(sep=" ",
        paste(sep="","Strike=", formatC(MEC$az1, digits=5) ),
        paste(sep="","Dip=", formatC(MEC$dip1, digits=5)   ),
        paste(sep="","Rake=", formatC(MEC$rake1)) )

      
      tit2 = focleg(MEC$icol)
      if(MEC$UP) { tit3 = "UPPER HEMI" } else {   tit3 = "LOWER HEMI" }
      
      text(0.3420201, 0.9396926, tit3, pos=4, font=2, xpd=TRUE)
      text(-0.3420201, 0.9396926, tit1, pos=2, font=2, xpd=TRUE)
      text(-.5,0.8660254, tit2, pos=2, font=2, xpd=TRUE)
      
      text(.95,-.9, paste(sep='', "CNVRG=", MEC$CNVRG), font=2, pos=2, xpd=TRUE)
    }
  return(MEC)
  
}

