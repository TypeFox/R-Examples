markseis24<-function(pjj, pix =list(yr=2009, jd=1, hr=0, mi=0, sec=0, dur=0), col='red', LEGON=3, BARON=TRUE, ARROWS=TRUE, lwd=1 )
{
############  after plotting plotseis24 use this to add marks on the plot
  ###   duration is set in seconds
  if(missing(BARON)) { BARON = TRUE }
  if(missing(LEGON)){ LEGON= 3 }
  if(missing(col)) { col='red' }
  if(missing(ARROWS)) { ARROWS=TRUE }
  if(missing(lwd)) {lwd=1  }

  alen = 0.06
  jd = pjj$jd
  yr = pjj$yr

  ydivs = pjj$y


wy = which( is.na(ydivs) )


  ydivs[wy] = -(  (wy-1)/24 )

     ydivs = c(0, ydivs , -1)


  
  ydiff = ydivs[3] - ydivs[2]

  
 
  if(is.null(pix$dur) )
    {
      dur=rep(0, length(pix$hr))
    }
  else
    {
      dur = pix$dur
    }
  
  if(length(dur)<length(pix$hr)) {
    print("warning: (markseis24) replace dur with dur[1]")
    dur=rep(dur[1], length(pix$hr))
    
  }
  
  pix = recdatel(pix)
  pix$dur = dur
  
  w = which(pix$yr==yr & pix$jd==jd)

##   print(w)
  
  if(length(w)<1)
    {
      print("No pix times are on the day displayed")
      return(0)

    }


  
  h1 = pix$hr[w]+pix$mi[w]/60+pix$sec[w]/3600
  h2 = h1 + pix$dur[w]/(3600)

  ypos1 = floor(h1)
  ypos2 = floor(h2)

  xpos1 = (h1-ypos1)*3600
  xpos2 = (h2-ypos2)*3600

  ypos1[ypos1>23] = 23
  ypos1[ypos1<0] = 0

  ypos2[ypos1>23] = 23
  ypos2[ypos1<0] = 0

  yy1 = ydivs[ypos1+2]+ydiff*.4
  yy2 = ydivs[ypos2+2]+ydiff*.2
  
###  print(yy1)
###    print(yy2)
###  

  nonwrapper = which(ypos1==ypos2)
  
  winmark(xpos1[nonwrapper],xpos2[nonwrapper] , side=1, bar=yy1[nonwrapper], leg=yy2[nonwrapper], col=col, lwd = lwd, arrows = ARROWS, alen = 0.06, LEGON=LEGON, BARON=BARON )


wrapper = which(ypos1!=ypos2)

  ###  these are win marks that span the hour cutoff and so need to be broken up

  if(length(wrapper)>0)
    {

      WY1 = ypos1[wrapper]
      WY2 = WY1
      yh1 = h1[wrapper]
      
      yh2 = rep(3600, length(yh1))
      
      xpos1 = (yh1-WY1)*3600
      xpos2 = (yh2-WY2)*3600

      yy1 = ydivs[WY1+2]+ydiff*.4

      yy2 = ydivs[WY2+2]+ydiff*.2


      winmark(xpos1,xpos2 , side=1, bar=yy1, leg=yy2, col=col, lwd = lwd,  alen = 0.06, LEGON=LEGON, BARON=BARON, arrows = ARROWS )


      WY1 = ypos2[wrapper]
      WY2 = WY1
      yh1 = h2[wrapper]
      
      yh2 = rep(0, length(yh1))
      
      xpos1 = (yh2-WY1)*3600
      xpos2 = (yh1-WY2)*3600

      yy1 = ydivs[WY1+2]+ydiff*.4

      yy2 = ydivs[WY2+2]+ydiff*.2

      winmark(xpos1,xpos2 , side=1, bar=yy1, leg=yy2, col=col, lwd = lwd, alen = 0.06, LEGON=2, BARON=BARON, arrows = ARROWS )

    }

  
}

