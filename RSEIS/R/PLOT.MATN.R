`PLOT.MATN` <-
function(ascd, tim=1, dt=1,  WIN=c(0,1), labs="", notes=notes, sfact=1,ampboost=0,  shift=NULL, LOG="", COL='red', add=1, AXES=1, units=NULL, VS=FALSE)
{
###  source("/home/lees/R_PAX/RSEIS/R/PLOT.MATN.R")
  
### plot a matrix of seismograms on a simple panel display
###   ascd = matrix(time, trace)

###  add = 1,2,3  if add=1 plot and add traces
###                  add =2 plot, but no traces
###                   add = 3 no plot, but add traces
  if(missing(sfact)) { sfact=1}
   if(missing(ampboost)) {ampboost=0 }
 
  if(missing(dt)) { dt=1}
  if(missing(LOG)) { LOG=""  }
  
  if(missing(add)) { add=1 }

##########  if AXES=0  no Y axes
##########  if AXES=1  plot scale for largest amplitude band and a multiplier for all others
##########  if AXES=2  plot a y-axis for each band, add units for scale, left side
##########  if AXES=3  plot a y-axis for each band, add units for scale, right side
##########  if AXES=4  plot a y-axis for each band, add units for scale, alternate sides
  
  
  if(missing(AXES)) {  AXES=1 } 
  if(missing(units)) {  units=colnames(ascd) } 
  if(missing(VS)) { VS=FALSE } 

  
  if(missing(tim))
    {
      tim = dt*seq(from=0,to=length(ascd[,1])-1)
    }

  if(missing(WIN))
    {
      WIN =range(tim)
    }
  if(missing(notes))
    {
      note.flag = FALSE
      notes = NULL
    }
  else
    {
      note.flag = TRUE
    }
  if(is.null(WIN)==TRUE){ WIN =range(tim) }

 
  
  tflag = tim>=WIN[1]&tim<=WIN[2]
  
  tr1 = 0.05
  tr2 = .9
  
  matsiz = dim(ascd)
  nn = matsiz[2]


    if(missing(shift)) { shift=rep(0, length=nn)  } 

  
    if(is.null(units)) { units=rep("",nn ) }


  
  if(missing(COL)) { COL=rep(1, nn)  }

  
  opar <- par(no.readonly = TRUE)


  if(AXES==3 | AXES==4)
    {
      mai = par("mai")
      mai[4] = mai[2]

      par(mai=mai )
    }
       
  

  ncol = length(COL)
  KX = length(COL)
  
  if(ncol<nn) {  COL=c(COL, rep(1, nn-ncol)) }

  if(missing(labs)) { labs=rep(" ", nn) }
  ttics = pretty(tim[tflag], n=10 )
  atics = ttics
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
  dy = (1/nn)
  maxS = rep(0,nn)
  minS = rep(0,nn)
  diffS = rep(0,nn)

  windiv = matrix(nrow=nn, ncol=4)
  
  for(i in 1:nn)
    {
      amp = ascd[tflag,i]
      lamp = length(amp[!is.na(amp)])
      if(lamp<1)
        {
          maxS[i] = 0
          minS[i] = 0
          diffS[i] = 0

        }
      else
        {
          
          maxS[i] = max(amp[!is.na(amp)])
          minS[i] = min(amp[!is.na(amp)])
          diffS[i] = maxS[i]-minS[i]
        }
    }
  ##  abs waiting using only COMP
  KDIFF = which.max(diffS)
  
  if(sfact>=2)
    {
      MAXy = max(maxS)
      MINy = min(minS)
      
      maxS =rep(MAXy, nn)
      minS =rep(MINy, nn)
    }

  if(add==1)
    {
      plot(range(tim[tflag], na.rm=TRUE ), c(0,1), type='n', axes=FALSE, xlab="", ylab="", log=LOG)
      
    }
  if(add==2)
    {
      plot(range(tim[tflag], na.rm=TRUE ), c(0,1), type='n', axes=FALSE, xlab="", ylab="", log=LOG)

    }

  
  box(col=grey(0.8))

  upar = par("usr")
  
  nump = length(tim[tflag])
  if(ncol>nn)
    {
      KR = nump/(ncol-1)
      KX = floor(seq(from=1, length=nump)/(KR))+1
    }
  
  for(i in 1:nn)
    {
      amp = ascd[tflag,i]
      lamp = length(amp[!is.na(amp)])
      if(lamp<1)
        {
          
          next;
        }
      
      amp = amp-mean(amp[!is.na(amp)])
      y3 = 1-(dy*i)
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

    

      ylow= y3-ampboost
      
      yhigh= y3+dy+ampboost

      
      z = RPMG::RESCALE(amp, ylow, yhigh, minamp, maxamp )


      
      windiv[i, ] = c(y3, y3+dy, minamp, maxamp)
      
      
      if(add!=3)
        {

          addtix(side=3, pos=y3+dy,   tck=0.005, at=ttics, labels=FALSE, col=gray(0.8) )
 
          abline(h=y3, lty=2, col=grey(0.8))
          
          axis(side=1, tck=0.01, at=ttics, labels=FALSE)
          axis(side=1, tick=FALSE,  at=ttics, labels=atics, line=-1)
          
          
          moretics = seq(from=min(ttics), to=max(ttics), by=1)
          if(length(moretics)<500)
            {
              axis(side=3, tck=0.01, at=moretics, labels=FALSE)
            }
          title(xlab='Time (s)', line=1.4, cex=1.2) 
        }
      

      if(add!=2)
        {
          if(ncol<=nn)
            {
              ex= tim[tflag]-shift[i]
              if(VS)
                {

                  polygon(c(ex, ex[1]) , c(z, mean(z)) , col=COL[i])

                }
                ex= tim[tflag]-shift[i]
              lines(ex, z, col=COL[i])
            }
          else
            {
              ni = length(z)
              E = tim[tflag]
              ex= tim[tflag]-shift[i]
                
              lines(ex, z, col=1)
              segments(E[1:(ni-1)],z[1:(ni-1)],E[2:ni],z[2:ni], col=COL[KX[i]])
              

            }
        }
                     #   print( paste(sep=' ', "IN PLOT.MATN", y3, y3+dy, minamp, maxamp))
      
      cmm = c(minamp, maxamp)
      lcmm = length(cmm[!is.na(cmm)])
      dmm = maxamp-minamp
      if( lcmm < 2   | dmm<=0)
        {
                                        #   print( paste(sep=' ', "IN PLOT.MATN", "PROBLEMS", lcmm ,dmm ))
          next;
        }
      yy = pretty(cmm, n = 5)
      
      flg = yy>minamp & yy<maxamp
      yt = yy[flg]
      yts = RPMG::RESCALE(yt, y3, y3+dy, minamp, maxamp )
      


                                        #   axis(2, tck=0.01 , at=yts, labels=yt, las=2 , line=0.1 )

      if(AXES==1)
        {
          ############  this is the default
          if(i==KDIFF)
            {
              if(add!=3)axis(2, pos= upar[1] ,tck=-0.005 , at=yts, labels=yt, las=2 , line=0.1 )
            }
          else
            {
              bnum = paste(sep='', "X", format.default(diffS[KDIFF]/diffS[i], digits=4))
              blab=bnum 
              if(add!=3) text(min(tim[tflag]), y3+0.75*dy, labels=blab, adj=0)
            }
        }


     if(AXES==2) 
        {
          #######  put a Y axis on each trace, left side
          side = 2; pos= upar[1]
          axis(side, pos=pos ,tck=-0.005 , at=yts, labels=yt, las=2 , line=0.1, cex=.7 )
          mtext(units[i], side =side , line = -1, at=mean(yts))
          
        }

     if(AXES==3) 
        {
          #######  put a Y axis on each trace, right side
          side = 4; pos= upar[2]
          axis(side, pos=pos ,tck=-0.005 , at=yts, labels=yt, las=2 , line=0.1, cex=.7 )
          mtext(units[i], side =side , line = -1, at=mean(yts))
          
        }

      if(AXES==4) 
        {
          #######  put a Y axis on each trace, alternate left and right
          if((i%%2) == 1) { side = 2; pos= upar[1] } else { side = 4; pos= upar[2]  }
          axis(side, pos=pos ,tck=-0.005 , at=yts, labels=yt, las=2 , line=0.1, cex=.7 )
          mtext(units[i], side =side , line = -1, at=mean(yts))
          
        }
 


      
      
                                        # axis(side=3, pos=y3+dy,   tck=0.005, at=ttics, labels=FALSE, col=2 )
      
      ylab = labs[i]
      mtext(side=2, at=y3+dy/2, text=ylab , line=1)
                                        #   print( paste(sep=' ', "IN PLOT.MATN",note.flag))
      
      if(note.flag==TRUE)
        {
                                        #  print( paste(sep=' ', "IN PLOT.MATN", notes[i]))
          
          if(add!=3)text(max(tim[tflag]), y3+dy-dy*0.1, notes[i], adj=1)
          
        }
      else
        {
          if(add!=3 & is.null(notes[i]) )text(max(tim[tflag]), y3+dy/2, labels=i, pos=4, col=gray(0.8))
        }
      
    }

  u = par("usr")
  
  colnames(windiv) = c("winymin", "winymax", "usrymin", "usrymax")
  
    if(AXES==3 | AXES==4)
    {
      par(opar)
    }
       
  
  invisible(list(n=nn, windiv=windiv) )
  
  

}

