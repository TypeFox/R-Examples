`PLOT.TTCURVE` <-
function(GH, STAXY=NULL, DIST=c(0,10), DY=.1, tim=1, dt=1,  sel=c(1:4), WIN=c(1,0), labs=c("CE1"), notes="CE1.V", tags="CE1.V", sfact=1,  COL='red', add=1, pts=FALSE, YAX=FALSE, TIT=NULL, SHIFT=NULL , rm.mean=TRUE, UNITS="volts", MARK=TRUE)
{
  ### plot a matrix of seismograms on a simple panel display
  ###   GH = structure of traces

  ###  add = 1,2,3  if add=1 plot and add traces
   ###                  add =2 plot, but no traces
  ###                   add = 3 no plot, but add traces
 ###   
  ###   sfact >= 2 = scale by window
###   source("/home/lees/Progs/R_PAX/RSEIS/R/PLOT.TTCURVE.R")


  if(missing(STAXY)) { print("NO STATION LIST"); return(0) }
  if(missing(DIST)) {DIST = range(STAXY$dist) }
 
 
  if(missing(sel)) { sel = 1:length(GH$JSTR) }

 if(missing(DY)) { DY =  diff(DIST) * (1/length(sel)) }
  
  
  if(missing(sfact)) { sfact=1}
  
  if(missing(dt)) { dt=rep(GH$info$dt[1],length(GH$JSTR)) }
  
  if(missing(add)) { add=1 }
  if(missing(pts)) {  pts=FALSE  }
  if(missing(YAX)) {  YAX=FALSE  }
  if(missing(TIT)) { TIT=NULL }
  if(missing(SHIFT)) { SHIFT=NULL }
  if(missing(UNITS)) { UNITS="volts" }
  if(missing(rm.mean)) { rm.mean=TRUE }
  if(missing(MARK)) { MARK=TRUE  }

  ulen = unlist(lapply(GH$JSTR, length))
  ireftrace = which.max(ulen)
  if(length(ireftrace)<1) return(-1)

  
  if(missing(tim))
    { 

      unlist(lapply(GH$JSTR, length))
      tim = GH$dt[ireftrace]*seq(from=0,to=length(GH$JSTR[[ireftrace]])-1)
    }

  if(missing(WIN))
    {
      WIN =range(tim)
    }
  if(missing(notes))
    {
      note.flag = FALSE
    }
  else
    {
      note.flag = TRUE
    }
  if(missing(tags))
    {
      tag.flag = FALSE
    }
  else
    {
      tag.flag = TRUE
      tag.col=1
    }

  if(is.list(WIN)==TRUE)
    {
      WIN = WIN$x
    }

  if(is.null(WIN)==TRUE){ WIN =range(tim) }
  
  ###  print(paste(sep=' ', "WIN", WIN[1], WIN[2]))
  ###  this following does not work as I expected it to.
  ### if(exists(deparse(substitute(WIN)))==FALSE){ WIN =range(tim) }
  
  ###  print(paste(sep=' ', "WIN", WIN[1], WIN[2]))
  mindist=DIST[1]; maxdist= DIST[2];

  
  tflag = tim>=WIN[1]&tim<=WIN[2]
  
  tr1 = 0.05
  tr2 = .9

 ###  print("plot.seisn sel")
  ### print(sel)
  
  nn = length(sel)

  
  
  if(missing(COL)) { COL=rep(1, nn)  }
if(length(COL)<nn) {  COL=c(COL, rep(1, nn-length(COL))) }
  
    if(missing(labs)) { labs=rep(NA, nn) }
  
  ttics = pretty(tim[tflag], n=10 )
  atics = ttics
 
  dy =  DY
  maxS = rep(0,nn)
  minS = rep(0,nn)
  diffS = rep(0,nn)
  meanS = rep(0,nn)
##############################
##############################
  #####   loop all selected traces to get min max and other stats for plotting:
  for(i in 1:length(sel))
    {
      ii = sel[i]

      LEN = length(GH$JSTR[[ii]])
        if(LEN<2)
          {
            next;
          }
       tim = GH$dt[ii]*seq(from=0,to=(LEN-1))
     
      tflag = tim>=WIN[1]&tim<=WIN[2]
      amp = GH$JSTR[[ii]][tflag]
      meanS[i] = mean(amp, na.rm =TRUE)
      if(rm.mean==TRUE) amp = amp -meanS[i]
     ###  print(range(amp[!is.na(amp)]))
      lamp = length(amp[!is.na(amp)])
     
      if(lamp<1)
        {
          maxS[i] = 0
          minS[i] = 0
          diffS[i] = 0
          meanS[i] = 0
        }
      else
        {
          
          maxS[i] = max(amp, na.rm=TRUE)
          minS[i] = min(amp, na.rm=TRUE)
          diffS[i] = maxS[i]-minS[i]
         
        }
    }
##############################
##############################
##############################  set up the plotting region:
  
      ##  abs weighting using only COMP
  KDIFF = which.max(diffS)
  
  if(sfact>=2)
    {
      MAXy = max(maxS, na.rm=TRUE)
      MINy = min(minS, na.rm=TRUE)
      
      maxS =rep(MAXy, nn)
      minS =rep(MINy, nn)
    }

  

  if(add==1 | add==2 )
    {
      
      plot(range(tim[tflag]), c(mindist-(maxdist-mindist)*.05, maxdist+(maxdist-mindist)*.05), type='n', axes=FALSE, xlab="", ylab="")
       
    }
   ##############################
   ##############################
   ##############################
   #########################  actual plotting  all traces ####################
   box(col=grey(0.8))

  upar = par("usr")


  MSTA = match( GH$STNS[sel], STAXY$name)
  ###   STAXY$distxyz[MSTA]
 ###ADIST =RPMG::RESCALE(STAXY$dist[MSTA], 0.05, 0.95, mindist, maxdist)
ADIST =STAXY$dist[MSTA]

###  sel =  which(GH$COMPS == "V")

  
 ###  plot(range(tim[tflag]), c(0,1), type='n', axes=FALSE, xlab="", ylab="")
       
  for(i in 1:length(sel))
    {
      ii = sel[i]
      tim = GH$dt[ii]*seq(from=0,to=length(GH$JSTR[[ii]])-1)
      tflag = tim>=WIN[1]&tim<=WIN[2]
      if(!is.null(SHIFT)==TRUE)
        {
          print(paste(sep=" ", i, ii, SHIFT[ii]))
          tim = tim-SHIFT[ii]
        }
      tflag = tim>=WIN[1]&tim<=WIN[2]      
      amp = GH$JSTR[[ii]][tflag]
      if(rm.mean==TRUE) amp = amp - meanS[i]
      tcol = COL[i]
      lamp = length(amp[!is.na(amp)])
      ###### print(paste(sep=' ',i, ii, lamp))
      if(lamp<1)
        {
          next;
        }

      ###  here I tried removing the mean value before plotting....this is wrong
      ###  amp = amp-mean(amp[!is.na(amp)])
      y3 = ADIST[i]
      
      if(sfact==1)
        {
          minamp =  min(amp[!is.na(amp)]);
          maxamp= max(amp[!is.na(amp)]);
        }
      else
        {
          minamp =  minS[i];
          maxamp= maxS[i];

        }
           ###  print( paste(sep=' ', "IN PLOT.SEISN", minamp, maxamp))

     ### if(add!=3) addtix(side=3, pos=y3+dy,   tck=0.005, at=ttics, labels=FALSE, col=gray(0.8) )
      
      z = RPMG::RESCALE(amp, y3-(dy/2), y3+(dy/2), minamp, maxamp )
      
    ###  if(add!=3)abline(h=y3, lty=2, col=grey(0.8))
      #############  draw seismograms
    if(add!=2)
      {
         lines(tim[tflag], z, col=tcol)
       }
      
      ylab = labs[i]
          
      if(note.flag==TRUE)
        {
                                        #  print( paste(sep=' ', "IN PLOT.MATN", notes[i]))
          
          if(add!=3)text(max(tim[tflag]), y3-dy*0.1, notes[i], adj=1, col=tcol)
          
        }

      if(add!=3 & tag.flag==TRUE)
        {
          text(max(tim[tflag]), y3, labels=tags[i], pos=4 , col=tcol)
        }
      
      if(add!=3 & tag.flag==FALSE)
        text(max(tim[tflag]), y3, labels=i, pos=4, col=gray(0.8))
      
    }


###  end plotting

  
  if(add!=3)
    {
      axis(side=1, tck=0.01, at=ttics, labels=FALSE)
      axis(side=1, tick=FALSE,  at=ttics, labels=atics, line=-1)
      
      
      moretics = seq(from=min(ttics), to=max(ttics), by=1)
      if(length(moretics)<500)
        {
          axis(side=3, tck=0.01, at=moretics, labels=FALSE)
        }
      title(xlab='Time (s)', ylab="Distance, km",  line=1.4, cex=1.2) 
    }
  
  u = par("usr")
  if(MARK)
    {
      
      ftime = Zdate(GH$info, sel[1], WIN[1])
      mtext( ftime, side = 3, at = u[1] , line=0.5, adj=0)
      
      ftime = Zdate(GH$info, sel[1], 0)
      mtext( ftime, side = 1, at = u[1] , line=1.5, adj=0)
      
    }
  
  
  
  if(!is.null(TIT))
    {

       mtext( TIT, side = 1, at = u[2] , line=1.5, adj=1)

    }

  invisible(list(n=nn, dy=dy,  minS=minS, maxS=maxS, meanS=meanS, DX=range(tim[tflag]), DY=DY, DIST=DIST ))
  
 

}

