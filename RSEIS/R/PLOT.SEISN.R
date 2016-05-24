`PLOT.SEISN` <-
    function(GH, tim=1, dt=1,  sel=c(1:4), WIN=c(1,0), labs=c("CE1"), notes="CE1.V", subnotes=NA, tags="CE1.V",
             sfact=1, LOG="", COL='red', add=1, pts=FALSE, YAX=1, TIT=NULL, SHIFT=NULL , COLLAPSE=FALSE, 
             rm.mean=TRUE, UNITS="volts", MARK=TRUE, xtickfactor = 1 )
{
  ### plot a matrix of seismograms on a simple panel display
  ###   GH = structure of traces

  ###  add = 1,2,3  if add=1 plot and add traces
   ###                  add =2 plot, but no traces
  ###                   add = 3 no plot, but add traces

  ###   sfact >= 2 = scale by window

####  COLLAPSE = TRUE means plot all traces on same panel
    
  if(missing(sel)) { sel = 1:length(GH$JSTR) }
  if(missing(sfact)) { sfact=1}
  
  if(missing(dt)) { dt=rep(GH$info$dt[1],length(GH$JSTR)) }
  
  if(missing(LOG)) { LOG=""  }
  
  if(missing(add)) { add=1 }
  if(missing(pts)) {  pts=FALSE  }
  if(missing(YAX)) {  YAX=1  }   ######  should be 1,2,3,4....different options
  if(missing(TIT)) { TIT=NULL }
  if(missing(SHIFT)) { SHIFT=NULL }
  if(missing(UNITS)) { UNITS="volts" }
  if(missing(rm.mean)) { rm.mean=TRUE }
  if(missing(MARK)) { MARK=TRUE  }

  if(missing(xtickfactor)) { xtickfactor = 1 }

  
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

   if(missing(subnotes))
    {
      subnote.flag = FALSE
    }
  else
    {
      subnote.flag = TRUE
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
  
  tflag = tim>=WIN[1]&tim<=WIN[2]
  
  tr1 = 0.05
  tr2 = .9

 ###  print("plot.seisn sel")
  ### print(sel)

  if(is.logical(sel)) { sel = which(sel)  } 
  
  nn = length(sel)
  
  if(missing(COL)) { COL=rep(1, nn)  }
if(length(COL)<nn) {  COL=c(COL, rep(1, nn-length(COL))) }

  if(length(COL)< length(GH$JSTR))
    {

      ncol = rep("black", times=length(GH$JSTR))
      ncol[sel] = COL
      COL = ncol

    }

  
  
    if(missing(labs)) {
      labs=rep(NA, nn)
      lab.flag = FALSE
    }
    else
      {
        lab.flag = TRUE
      }

  
  ##  print(paste("i am here",xtickfactor ))
  if(xtickfactor==1)
    {
      XLAB = "Time(s)"
      ttics = pretty(tim[tflag], n=10 )
      atics = ttics
    }
  else
    {
      
      ttics = pretty(  tim[tflag]/xtickfactor  , n=10 )
      atics = ttics
      ttics = ttics*xtickfactor
      if(xtickfactor==86400)   XLAB = "Time(Days)"
      if(xtickfactor==3600)    XLAB = "Time(Hours)"
      if(xtickfactor==60)      XLAB = "Time(Minutes)"
      if(xtickfactor>31535999  )      XLAB = "Time(Years)"
      
      
    }

  
  
  if(LOG=='x')
    {
      periods = c(30,20,10,5,2,1)
      hz = 1/periods
      at1 = c(pretty(1:10), pretty(tim))
      at2 = at1[at1>0&at1<max(tim)]
      ttics = c(hz, at2 )
      
      btics = c(periods, at2 )
      atics = as.character(btics)
      atics[length(atics)] = paste(sep=' ', atics[length(atics)],"Hz")
      
      atics[btics==1] = paste(sep=' ', atics[btics==1],"Hz")
           atics[1] = paste(sep=' ', atics[1],"s")
 
      
  }
  ########  here we must set the params for the individualized traces
  dy = (1/nn)
  if(COLLAPSE) dy = 1
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
      plot(range(tim[tflag]), c(0,1), type='n', axes=FALSE, xlab="", ylab="", log=LOG)
       
    }
   ##############################
   ##############################
   ##############################
   #########################  actual plotting  all traces ####################
   box(col=grey(0.8))

  upar = par("usr")
  
  for(i in 1:length(sel))
    {
      ii = sel[i]
      tim = GH$dt[ii]*seq(from=0,to=length(GH$JSTR[[ii]])-1)
      tflag = tim>=WIN[1]&tim<=WIN[2]
      if(!is.null(SHIFT)==TRUE)
        {
       ###   print(paste(sep=" ", i, ii, SHIFT[ii]))
          tim = tim-SHIFT[ii]
        }
      tflag = tim>=WIN[1]&tim<=WIN[2]      
      amp = GH$JSTR[[ii]][tflag]
      if(rm.mean==TRUE) amp = amp - meanS[i]
      tcol = COL[ii]
      lamp = length(amp[!is.na(amp)])
     ##  print(paste(sep=' ',i, ii, lamp))
      y3 = 1-(dy*i)
      if(COLLAPSE) y3 = 0
         if(note.flag==TRUE)
        {
                                      ##   print( paste(sep=' ', "IN PLOT.MATN", i, notes[i]))
          
          if(add!=3)text(max(tim[tflag]), y3+dy-dy*0.1, notes[i], adj=1)
          
        }

      if(lamp<1)
        {
            ##  if the trace is empty we have a problem.
           ##   print( paste(sep=' ', "IN PLOT.SEISN", "PROBLEMS: no samples?" ))
          next;
        }

      ###  here I tried removing the mean value before plotting....this is wrong
      ###  amp = amp-mean(amp[!is.na(amp)])

      
     
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

      if(add!=3) addtix(side=3, pos=y3+dy,   tck=0.005, at=ttics, labels=FALSE, col=gray(0.8) )
      z = RPMG::RESCALE(amp, y3, y3+dy, minamp, maxamp )

      
      if(add!=3)abline(h=y3, lty=2, col=grey(0.8))
      if(add!=2)lines(tim[tflag], z, col=tcol)
      if(pts==TRUE)points(tim[tflag], z, col=tcol, pch=4)
      
   ###  print( paste(sep=' ', "IN PLOT.SEISN", y3, y3+dy, minamp, maxamp))
      
      cmm = c(minamp, maxamp)
      lcmm = length(cmm[!is.na(cmm)])
      dmm = maxamp-minamp
      if( lcmm < 2   | dmm<=0)
        {
            ##    print( paste(sep=' ', "IN PLOT.SEISN", "PROBLEMS", lcmm ,dmm ))
            cmm = c(0, 1)
            yy = c(minamp, maxamp)
      
            
            yt = yy 
      
            yts = RPMG::RESCALE(yt, y3, y3+dy, minamp, maxamp )
          ## next;
        }
      else
          {
      yy = pretty(cmm, n = 5)
      
      flg = yy>minamp & yy<maxamp
      yt = yy[flg]
      
      yts = RPMG::RESCALE(yt, y3, y3+dy, minamp, maxamp )
  }
       ### print(paste(sep =  ' ' ,minamp,maxamp,  paste(collapse=" ", yt) ))
                                        #

      #########  plot all the mini axes on the y-axis if YAX == 2
      axis.side = 2
      axis.pos = upar[1]
      
      ylab = labs[i]
      vfonts=c("serif", "plain" )
      lab.pos = 1.2
      if(YAX == 2)
        {

          axis(axis.side, pos=axis.pos   ,tck=-0.005 , at=yts, labels=yt, las=2 , line=0.1 )
          if( any( !is.na(UNITS) ) )
            {
              ##   mtext(side=axis.side, at=y3+dy/2, text=ylab , line=-1)
             ## text(x=axis.pos, y=y3+dy/2, labels=ylab , vfont=vfonts , srt=90,  pos=4)
              text(x=axis.pos, y=y3+dy/2, labels=ylab , vfont=vfonts, srt=90  , adj=c(.5, lab.pos  ), xpd=TRUE )  
            }
        }
      #####################  alternate axis side if YAX == 3
      if(YAX == 3) {


         if( (i%%2)==0 )
        {
          axis.side = 2
          axis.pos = upar[1]
          lab.pos = 1.2
        }
      else
        {
          axis.side = 4
          axis.pos = upar[2]
          lab.pos = -1
        }

         
      
        axis(axis.side, pos=axis.pos   ,tck=-0.005 , at=yts, labels=yt, las=2 , line=0.1 )
         if( any( !is.na(UNITS) ) )
           {
            
             ##    mtext(side=axis.side, at=y3+dy/2, text=ylab , line=-1, font=vfonts)
             text(x=axis.pos, y=y3+dy/2, labels=ylab , vfont=vfonts, srt=90  , adj=c(.5,  lab.pos ), xpd=TRUE )

             
           }
         

      }
      
      
      if(i==KDIFF & YAX <2 )
        {
          if(add!=3)
           if(!is.na(UNITS) )  axis(2, pos= upar[1] ,tck=-0.005 , at=yts, labels=yt, las=2 , line=0.1 )
        }
      else
        {

           labstring= paste(sep="", "\\mu",format.default(diffS[KDIFF]/diffS[i], digits=4) )

          ##  bnum = paste(sep='', "X", format.default(diffS[KDIFF]/diffS[i], digits=4))
          ##  blab=bnum 
          if(add!=3)text(min(tim[tflag]), y3+0.75*dy, labels=labstring , adj=0, vfont=c("sans serif", "plain"))
        }

      
      
      # axis(side=3, pos=y3+dy,   tck=0.005, at=ttics, labels=FALSE, col=2 )
      
      
      ##  mtext(side=2, at=y3+dy/2, text=ylab , line=1)

      
      
                                        #    print( paste(sep=' ', "IN PLOT.SEISN",note.flag))
      
   
      if(subnote.flag==TRUE)
        {
                                        #  print( paste(sep=' ', "IN PLOT.MATN", notes[i]))
          
          if(add!=3)text(max(tim[tflag]), y3+dy*0.1, subnotes[i], adj=0)
          
        }

      if(add!=3 & tag.flag==TRUE)
        {
          text(max(tim[tflag]), y3+dy/2, labels=tags[i], pos=4 , col=tag.col)
        }
      
      if(add!=3 & tag.flag==FALSE)
        text(max(tim[tflag]), y3+dy/2, labels=i, pos=4, col=gray(0.8))
      
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
      
      title(xlab=XLAB, line=1.4, cex=1.2)


      
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

  invisible(list(n=nn, dy=dy,  minS=minS, maxS=maxS, meanS=meanS, DX=range(tim[tflag]) ))
  
 

}

