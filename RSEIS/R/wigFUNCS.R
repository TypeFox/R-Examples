#### BUTTONS

##########################################  BUttons

###  first button simply returns the name of the button pushed
######    and the clicks and picks

NEXT<-function(nh, g)
{
  #####  BUTTONDOC:NEXT:'Next BATCH of FILES'

  if(g$zenclick>1)
    {
      rd = getrdpix(g$zloc, g$zenclick, g$sel, nh)
    }
  else
    {
      rd=list(PUSHED="NEXT")
    }
  g$action = "break"
  g$rd = rd
  g$zloc = list(x=NULL, y=NULL)
  invisible(list(global.vars=g) )	

}

PREV<-function(nh, g)
{
#####  BUTTONDOC:PREV:'Previous BATCH of FILES'
  if(g$zenclick>1)
    {
      rd = getrdpix(g$zloc, g$zenclick, g$sel, nh)
    }
  else
    {
      rd=list(PUSHED="PREV")
    }
  g$action = "break"
  g$rd = rd
  g$zloc = list(x=NULL, y=NULL)
  invisible(list(global.vars=g))	

}


HALF<-function(nh, g)
{
#####  BUTTONDOC:HALF:'Shift Half a window'
  if(g$zenclick>1)
    {
      rd = getrdpix(g$zloc, g$zenclick, g$sel, nh)
    }
  else
    {
      rd=list(PUSHED="HALF")
    }
  g$action = "break"
  g$zloc = list(x=NULL, y=NULL)
  g$rd = rd
  invisible(list(global.vars=g))	

}

CENTER<-function(nh, g)
  {
#####  BUTTONDOC:CENTER:'Center a window'
    if (g$zenclick > 1) {
        rd = getrdpix(g$zloc, g$zenclick, g$sel, nh)
    }
    else {
        rd = list(PUSHED = "CENTER")
    }
    g$action = "break"
    g$rd = rd
    g$zloc = list(x = NULL, y = NULL)
    invisible(list(global.vars = g))


  }


MARK<-function(nh, g)
{
#####  BUTTONDOC:MARK:'Mark a trace' 
  if(g$zenclick>1)
    {
      rd = getrdpix(g$zloc, g$zenclick, g$sel, nh)
    }
  else
    {
      rd=list(PUSHED="MARK")
    }
  g$action = "break"
  g$rd = rd
  g$zloc = list(x=NULL, y=NULL)
  invisible(list(global.vars=g))	

}

##########################################
DOC<-function(nh, g)
{
  #####  BUTTONDOC:DOC:'Show documentation' 
  PICK.DOC(g$BLABS)
  g$zloc = list(x=NULL, y=NULL)
  g$action = "replot"
  invisible(list(global.vars=g))	
}
##########################################
REFRESH<-function(nh, g)
  {
    #####  BUTTONDOC:REFRESH:'Refresh screen' 
    u = par("usr")
    L = length(g$sloc$x)
    if(L>1)
      {
        abline(v=g$sloc$x[c(L-1,L)], col=gray(0.8), lty=2)
      }
    g$sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
    g$zloc = list(x=NULL, y=NULL)
    g$action = "replot"
    invisible(list(global.vars=g))	

  }
#######
RESTORE<-function(nh, g)
  {
    #####  BUTTONDOC:RESTORE:'Restore from zoom' 

    u = par("usr")
    L = length(g$sloc$x)
    if(L>1)
      {
        abline(v=g$sloc$x[c(L-1,L)], col=gray(0.8), lty=2)
      }
    g$sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
    g$zloc = list(x=NULL, y=NULL)
    g$WIN = NULL
    g$action = "replot"
    invisible(list(global.vars=g))	
  }
#######
ZOOM.out<-function(nh, g)
  {
#####  BUTTONDOC:ZOOM.out:'Zoom out' 
    
    u = par("usr")
    DX = (u[2]-u[1])*0.3
    zloc = list(x= c(u[1]-DX, u[2]+DX))
    g$WIN = zloc$x
    
    g$zloc = list(x=NULL, y=NULL)
      
    g$action = "replot"
    invisible(list(global.vars=g))	
  }

ZOOM.in<-function(nh, g)
  {
    #####  BUTTONDOC:ZOOM.in:'Zoom in' 
    zenclick = length(g$zloc$x)
    if(zenclick>=3)
      {
        n1=g$zenclick-2
        pwin = sort(g$zloc$x[c(n1,n1+1)])
        g$WIN = pwin
      }
    else
      {
        u = par("usr")
        DX = (u[2]-u[1])*0.3
        zloc = list(x= c(u[1]+DX, u[2]-DX))
        g$WIN = zloc$x
      }

    
    g$zloc = list(x=NULL, y=NULL)
    
    g$action = "replot"
    invisible(list(global.vars=g))	
  }
#######
LEFT<-function(nh, g)
  {
    #####  BUTTONDOC:LEFT:'Shift Left'
    u = par("usr")
    DX = (u[2]-u[1])*0.3
####  zloc = list(x= c(u[1]+DX, u[2]+DX))
    g$WIN  =c(u[1]-DX, u[2]-DX)
    
    g$zloc = list(x=NULL, y=NULL)
    
    g$action = "replot"
    invisible(list(global.vars=g))	

  }

RIGHT<-function(nh, g)
  {
    #####  BUTTONDOC:RIGHT:'Shift Right'
    u = par("usr")
    DX = (u[2]-u[1])*0.3
####  zloc = list(x= c(u[1]+DX, u[2]+DX))
    g$WIN  =c(u[1]+DX, u[2]+DX)
    
    g$zloc = list(x=NULL, y=NULL)
    
    g$action = "replot"
    invisible(list(global.vars=g))	

  }







#######
SCALE<-function(nh, g)
  {
    #####  BUTTONDOC:SCALE:'Toggle Scale by trace/window'
             if(g$ScaleFACT==1)
            {

              g$ScaleFACT=2
            }
          else
            {
              g$ScaleFACT=1

            }

    g$action = "replot"
     g$zloc = list(x=NULL, y=NULL) 
          
    invisible(list(global.vars=g))	
  }
########################################
Xwin<-function(nh, g)
  {
#####  BUTTONDOC:Xwin:'Delete all windows except main'
          ALLdevs =     dev.list()


          ww  = ALLdevs[ which(g$MAINdev != ALLdevs)]

              for(i in 1:length(ww))
                {
                  dev.off(which = ww[i])

                }
              
         
          dev.set(g$MAINdev)
            g$zloc = list(x=NULL, y=NULL) 
     
    g$action="donothing"
    invisible(list(global.vars=g))
 
  }



########################################
PSEL<-function(nh, g)
        {
          #####  BUTTONDOC:PSEL:'Pick trace Sta/COMP to show' 

          sel = SELSTA(nh, sel=g$sel, newdev=TRUE, STAY=FALSE)
          
          NSEL = length(nh$dt[g$sel])

          g$du = 1/NSEL
         
          isel = sel[1]
          
          Torigin = list(jd=nh$info$jd[isel], hr=nh$info$hr[isel],
            mi=nh$info$mi[isel],
            sec=(nh$info$sec[isel]+nh$info$msec[isel]/1000+nh$info$t1[isel]-nh$info$off[isel]))

          g$Torigin=Torigin
          g$sel = sel

          
          g$STNS = nh$STNS[sel]
          g$COMPS = nh$COMPS[sel]
           g$zloc = list(x=NULL, y=NULL) 
 
          g$action = "replot"
           invisible(list(global.vars=g))
        }
#######################################
####  this needs work
FLIP<-function(nh, g)
  {
    #####  BUTTONDOC:FLIP:'Flip selected trace' 
    zenclick = length(g$zloc$x)

    if(zenclick>1)
      {
        nc = 1:(zenclick-1)
        lnc = length(nc)
        
        ypick = length(g$sel)-floor(length(g$sel)*g$zloc$y[nc])
        ipick = unique( g$sel[ypick] )

        cat("FLIP: pwig POLARITY REVERSED: "); cat(ipick, sep=" " ); cat("\n")
        
        for(JJ in 1:length(ipick) )
          {
            jtr  = ipick[JJ]
            nh$JSTR[[jtr]] = (-1)*nh$JSTR[[jtr]]
          }
      }
    else
      {
        cat("FLIP: No traces selected: Try Again"); cat("\n")

      }
     g$zloc = list(x=NULL, y=NULL) 
 
    g$action = "replace"
    invisible(list(NH=nh, global.vars=g))

    
  }        
########################
PTS<-function(nh, g)
  {
#####  BUTTONDOC:PTS:'Show sample points' 

    
    g$pts=!g$pts
    g$action = "replot"
     g$zloc = list(x=NULL, y=NULL) 
 
    invisible(list(global.vars=g))
    
  }

FILT<-function(nh, g)
  {
#####  BUTTONDOC:FILT:'Filter trace'
    ### print( data.frame(g$filters) )

     
    Fdef = choosfilt(thefilts=g$filters, ncol=5)

    if(!is.null(Fdef))
      {

        if(Fdef$type=="None")
          {
            dev.set( g$MAINdev)
            g$SUBTIT = NA
            
            g$action = "revert"
             KF = nh
            return(list(global.vars=g))
          }
        else
          {
           ###  g$SUBTIT = paste(Fdef$type,Fdef$fl, Fdef$fh, sep=" ")
           g$SUBTIT =   filterstamp(Fdef$fl, Fdef$fh, Fdef$type)
             g$action = "replace"
              
            KF = FILT.SEISN(nh, sel = g$sel, FILT=Fdef)
          }
###  X11()
      }
    else
      {

        
        g$action = "replot"
        KF = nh
        

      }
 g$zloc = list(x=NULL, y=NULL) 
 
    dev.set( g$MAINdev)
   
    invisible(list(NH=KF, global.vars=g))
    
    
  }


UNFILT<-function(nh, g)
  {
    #####  BUTTONDOC:UNFILT:'Unfilter traces'
    dev.set( g$MAINdev)
    g$SUBTIT = NA
    g$action = "revert"
     g$zloc = list(x=NULL, y=NULL) 
 
    invisible(list(global.vars=g))
  }
#########################



fspread<-function(nh, g)
  {
    #####  BUTTONDOC:fspread:'do a filter spread on selection' 

    ###  click on a trace panel and do a filter spread
    

 zenclick = length(g$zloc$x)


    if(zenclick>=3)
      {
        ypick = length(g$sel)-floor(length(g$sel)*g$zloc$y[zenclick-1])
        ipick = g$sel[ypick]
        print(paste(sep=' ',"fspread", ypick, nh$info$name[ ipick]))

        famp = nh$JSTR[[ipick]]

        pwin = sort(c(g$zloc$x[zenclick-2], g$zloc$x[zenclick-1]))

        ex = seq(from=nh$info$t1[ipick], by=nh$info$dt[ipick], length.out=length(famp) )
        temp =  famp[ ex > pwin[1] & ex <pwin[2]]

        
#### Xamp =  -1*temp
        smallex = ex[ ex > pwin[1] & ex <pwin[2]]

        asec = nh$info$sec[ipick]+nh$info$msec[ipick]/1000+nh$info$t1[ipick]-nh$info$off[ipick]+pwin[1]
        
        spaz = recdate( nh$info$jd[ipick], nh$info$hr[ipick], nh$info$mi[ipick], asec,  nh$info$yr[ipick] )
        
        spaz$yr =   as.integer(nh$info$yr[ipick])
        
        MODAY = getmoday(spaz$jd,  spaz$yr)
        
        TP = list(yr=spaz$yr[1], jd=spaz$jd, mo=MODAY$mo,
          dom= MODAY$dom  ,hr=spaz$hr, mi=spaz$mi, sec=spaz$sec )

       dst = dateStamp(TP)

        titl = paste(nh$STNS[ipick], nh$COMPS[ipick], dst)
        
        fh=c(1/20, 1/10, 1/5, .5, 1, 2, 3)
        fl=rep(1/100, times=length(fh) )

        
         dev.new(width=14, height=10)

        jex = range(smallex)
        jr =  jex[2] - jex[1]
        j10 = jr*0.2

        jwin = c(jex[1]+j10, jex[2]-j10)
       #  jwin = NULL
            
        FILT.spread(smallex, temp, nh$dt[ipick], fl = fl, fh = fh, sfact = 1, WIN = jwin, PLOT = TRUE, TIT =titl , TAPER = 0.1, POSTTAPER=NULL)

        dev.set(g$MAINdev)
        
        g$zloc = list(x=NULL, y=NULL) 
        g$action="donothing"
        invisible(list(global.vars=g))
        
        
      }
    else
      {
        cat("XTR WARNING: no window or trace has been selected:", sep="\n")
        RETX=NULL
        g$zloc = list(x=NULL, y=NULL) 
        
        g$action="donothing"
        invisible(list(global.vars=g))
        
        
      }

  }
















SPEC<-function(nh, g)
  {
    #####  BUTTONDOC:SPEC:'Display Spectrum' 

    nclick = length(g$zloc$x)

    if(nclick>=3)
      {
        nc = 1:(nclick-1)
        lnc = length(nc)
        
        ypick = length(g$sel)-floor(length(g$sel)*g$zloc$y[nc])
        ipick = g$sel[ypick]
### print(paste(sep=' ',ypick, NH$info$name[ ipick]))

        print(ipick)
        
        i1 = seq(from=1, to=max(nc), by=2)
        i1 = i1[i1<max(nc)]
        amp = list()
        dees = list()
        stamps =  list()
        speccol = vector()
        ni = 0

        for(ipix in i1)
          {
            pwin = sort(c(g$zloc$x[ipix], g$zloc$x[ipix+1]))
            print(c(ipix, pwin))

            kpix = ipick[ipix]

            famp = nh$JSTR[[kpix]]


            ex = seq(from=nh$info$t1[kpix], by=nh$info$dt[kpix], length.out=length(famp) )
            temp =  famp[ ex > pwin[1] & ex <pwin[2]]


            if(any(is.na(temp)))
              {
                print(paste("getting NA in trace",kpix, nh$STNS[kpix],nh$COMPS[kpix],pwin[1],  pwin[2]  ))
                next
              }
            
            ni = ni +1

            amp[[ni]] = temp-mean(temp)
            dees[ni] = nh$dt[kpix]

            speccol[ni] = g$pcols[kpix]

            ftime = Zdate(nh$info, kpix, pwin[1])
            psta = nh$STNS[kpix]
            pcomp =  nh$COMPS[kpix]
            STAMP = paste(sep=" ", psta, pcomp, ftime)
            stamps[ni] = STAMP
            
          }

        print(stamps)
        a = list(y=amp, dt=dees, stamps=stamps)

        if(length(a$y)>0)
          {
            dev.new(width=10, height=10)
            
            f1 = 0.1
            f2 = floor(0.33*(1/nh$dt[ipick]))
            
###  oop=par(no.readonly = TRUE)
###  par(mfrow=c(length(a$y), 1) )
###  for(io in 1:length(a$y)) plot(a$y[[io]], type='l')
###  par(oop)
###  readline("type in something")
            
            MTM.drive(a, f1, f2[1], COL=speccol, PLOT=TRUE)
          }
        dev.set(g$MAINdev)
      }
    else
      {
        cat("SPEC WARNING: no window or trace has been selected:", sep="\n")
      }
    

    g$zloc = list(x=NULL, y=NULL) 

    g$action="donothing"
    invisible(list(global.vars=g))
  }


WWIN<-function(nh, g)
  {
    #####  BUTTONDOC:WWIN:'Window' 
    nclick = length(g$zloc$x)

    if(nclick>=3)
      {
        nc = 1:(nclick-1)
        lnc = length(nc)
        
        ypick = length(g$sel)-floor(length(g$sel)*g$zloc$y[nc])
        ipick = g$sel[ypick]
### print(paste(sep=' ',ypick, NH$info$name[ ipick]))

        print(ipick)
        
        i1 = seq(from=1, to=max(nc), by=2)
        i1 = i1[i1<max(nc)]
        amp = list()
        dees = list()
        stamps =  list()
        speccol = vector()
        ni = 0

        for(ipix in i1)
          {
            pwin = sort(c(g$zloc$x[ipix], g$zloc$x[ipix+1]))
            print(c(ipix, pwin))

            kpix = ipick[ipix]

            famp = nh$JSTR[[kpix]]


            ex = seq(from=nh$info$t1[kpix], by=nh$info$dt[kpix], length.out=length(famp) )
            temp =  famp[ ex > pwin[1] & ex <pwin[2]]


            if(any(is.na(temp)))
              {
                print(paste("getting NA in trace",kpix, nh$STNS[kpix],nh$COMPS[kpix],pwin[1],  pwin[2]  ))
                next
              }
            
            ni = ni +1

            amp[[ni]] = temp-mean(temp)
            dees[ni] = nh$dt[kpix]

            speccol[ni] = g$pcols[kpix]

            ftime = Zdate(nh$info, kpix, pwin[1])
            psta = nh$STNS[kpix]
            pcomp =  nh$COMPS[kpix]
            STAMP = paste(sep=" ", psta, pcomp, ftime)
            stamps[ni] = STAMP
            
          }

        dev.new(width=10, height=10)

        for(i in 1:ni) {
          y = amp[[i]]
          len = length(amp[[i]])
          print(c(i, len, dees[[i]]))
          xt = seq(from=0, by=dees[[i]], length=len)
          plot(xt , amp[[i]],main=stamps[[i]], type='l');
          
          locator(1) }

         dev.set(g$MAINdev)
        
      }


    g$zloc = list(x=NULL, y=NULL) 

    g$action="donothing"
    invisible(list(global.vars=g))
 
    
  }
##########################


SGRAM<-function(nh, g)
  {
    #####  BUTTONDOC:SGRAM:'Spectrogram' 

 zenclick = length(g$zloc$x)


          if(zenclick>=2)
            {
              if(zenclick==2)
                {
                  ypick = length(g$sel)-floor(length(g$sel)*g$zloc$y[zenclick-1])
                  ipick = g$sel[ypick]
###  print(paste(sep=' ',ypick, NH$info$name[ipick]))
                  pwin = g$WIN
                }
              else
                {
                  
                  ypick = length(g$sel)-floor(length(g$sel)*g$zloc$y[zenclick-1])
                  ipick = g$sel[ypick]
### print(paste(sep=' ',ypick, NH$info$name[ ipick]))
                  pwin = sort(c(g$zloc$x[zenclick-2], g$zloc$x[zenclick-1]))
                }
           
          g$LASTwin = pwin
   
          ### print(paste(sep=" ", "DOING SGRAM  Nclick, ipick, pwin", Nclick, ipick, pwin))
          
          famp = nh$JSTR[[ipick]]
         
          ex = seq(from=nh$info$t1[ipick], by=nh$info$dt[ipick], length.out=length(famp))
         
          temp =  famp[ ex > pwin[1] & ex <pwin[2]]

          Xamp =   temp-mean(temp)

       #   ftime = Zdate(nh$info, g$sel[ypick], pwin[1])
            ftime =  ghstamp(nh, sel=g$sel[ypick], WIN=pwin )

              
	print(paste(sep=" ",min(ex), max(ex)))
	
	print(paste(sep=" ",pwin[1], pwin[2]))
	
          print(paste(sep=" ", ipick, length(famp),length(temp),length(Xamp), nh$dt[ipick],ftime)) 

          SPECT.drive(Xamp, DT=nh$dt[ipick], STAMP=ftime)

        ###   plotevol(DEV, log=1, fl=0, fh=15, col=rainbow(50))
           }
           else
             {
              pwin = g$LASTwin
              ypick = 1
              ipick = g$sel[1]
              cat("SGRAM WARNING: no window or trace has been selected:" , sep="\n")
            }

          dev.set(g$MAINdev)
            g$zloc = list(x=NULL, y=NULL) 
     
    g$action="donothing"
    invisible(list(global.vars=g))
 
    

  }


WLET<-function(nh, g)
  {
    #####  BUTTONDOC:WLET:'Wavelet Transform'
 zenclick = length(g$zloc$x)


          if(zenclick>=2)
            {
              if(zenclick==2)
                {
                  ypick = length(g$sel)-floor(length(g$sel)*g$zloc$y[zenclick-1])
                  ipick = g$sel[ypick]
###  print(paste(sep=' ',ypick, NH$info$name[ipick]))
                  pwin = g$WIN
                }
              else
                {
                  
                  ypick = length(g$sel)-floor(length(g$sel)*g$zloc$y[zenclick-1])
                  ipick = g$sel[ypick]
### print(paste(sep=' ',ypick, NH$info$name[ ipick]))
                  pwin = sort(c(g$zloc$x[zenclick-2], g$zloc$x[zenclick-1]))
                }
           
          g$LASTwin = pwin
   
          ### print(paste(sep=" ", "DOING SGRAM  Nclick, ipick, pwin", Nclick, ipick, pwin))
          
          famp = nh$JSTR[[ipick]]
         
          ex = seq(from=nh$info$t1[ipick], by=nh$info$dt[ipick], length.out=length(famp))
         
          temp =  famp[ ex > pwin[1] & ex <pwin[2]]

          Xamp =   temp-mean(temp)

          ## ftime = Zdate(nh$info, g$sel[ypick], pwin[1])

              ftime =  ghstamp(nh, sel=g$sel[ypick], WIN=pwin )


  
              wlet.drive(Xamp, nh$dt[ipick], STAMP=ftime)
            

              
        ###   plotevol(DEV, log=1, fl=0, fh=15, col=rainbow(50))
           }
           else
             {
              pwin = g$LASTwin
              ypick = 1
              ipick = g$sel[1]
              cat("WLET WARNING: no window or trace has been selected:" , sep="\n")
            }

          dev.set(g$MAINdev)
            g$zloc = list(x=NULL, y=NULL) 
     
    g$action="donothing"
    invisible(list(global.vars=g))
 
    

  }


XTR<-function(nh, g)
  {
    #####  BUTTONDOC:XTR:'Extract single trace' 
    zenclick = length(g$zloc$x)


    if(zenclick>=3)
      {
        ypick = length(g$sel)-floor(length(g$sel)*g$zloc$y[zenclick-1])
        ipick = g$sel[ypick]
        print(paste(sep=' ',"EXTRACT", ypick, nh$info$name[ ipick]))

        famp = nh$JSTR[[ipick]]

        pwin = sort(c(g$zloc$x[zenclick-2], g$zloc$x[zenclick-1]))

        ex = seq(from=nh$info$t1[ipick], by=nh$info$dt[ipick], length.out=length(famp) )
        temp =  famp[ ex > pwin[1] & ex <pwin[2]]

        
#### Xamp =  -1*temp
        smallex = ex[ ex > pwin[1] & ex <pwin[2]]

        asec = nh$info$sec[ipick]+nh$info$msec[ipick]/1000+nh$info$t1[ipick]-nh$info$off[ipick]+pwin[1]
        
        spaz = recdate( nh$info$jd[ipick], nh$info$hr[ipick], nh$info$mi[ipick], asec,  nh$info$yr[ipick] )
        
        spaz$yr =   as.integer(nh$info$yr[ipick])
        
        MODAY = getmoday(spaz$jd,  spaz$yr)
        
        TP = list(yr=spaz$yr, jd=spaz$jd, mo=MODAY$mo,
          dom= MODAY$dom  ,hr=spaz$hr, mi=spaz$mi, sec=spaz$sec )

        RETX = list(but="RET", x=smallex, y=temp, dt=nh$dt[ipick], STNS=nh$STNS[ipick],
          COMPS=nh$COMPS[ipick],  fname=nh$info$name[ipick] , TIMEpick=TP, mark=TRUE, deltat=nh$dt[ipick] )
        g$zloc = list(x=NULL, y=NULL) 
        
        g$action="exit"
        invisible(list(RETX = RETX, global.vars=g))
        
        
      }
    else
      {
        cat("XTR WARNING: no window or trace has been selected:", sep="\n")
        RETX=NULL
        g$zloc = list(x=NULL, y=NULL) 
        
        g$action="donothing"
        invisible(list(global.vars=g))
        
        
      }


    
  }

########################################

Pinfo<-function(nh, g)
  {
    #####  BUTTONDOC:Pinfo:'Pick information' 
    zenclick = length(g$zloc$x)

    if(zenclick>=2)
      {
       ### NSEL = length(nh$dt[g$sel])

       ### du = 1/NSEL
          
        kix = legitpix(g$sel, g$zloc, zenclick)
        ypick =  kix$ypick
        ppick = kix$ppick
        
        dpick = c(0, diff(ppick))
        ipick = g$sel[ypick]

        m = match(g$STNS[ipick],g$UNIsta)
      ###  jj = floor(( g$zloc$y[zenclick-1])/du)
        asec = nh$info$sec[ipick]+nh$info$msec[ipick]/1000+nh$info$t1[ipick]-nh$info$off[ipick]+ppick[zenclick-1]

        print(paste(sep=" ", "PICK=",
                    nh$info$yr[ipick], nh$info$jd[ipick], nh$info$hr[ipick],
                    nh$info$mi[ipick], asec, "sta=", nh$STNS[ipick], "comp=", nh$COMPS[ipick] ))
        print(ppick)

        ##  pstas = paste(nh$STNS[ipick], nh$COMPS[ipick], sep=".")

        rd = getrdpix(g$zloc, zenclick, g$sel, nh)
        
        RDtmes = rd$yr+rd$jd/366+rd$hr/(366*24)+rd$mi/(366*24*60)+rd$sec/(366*24*3600)
        
        wearliest = which.min(RDtmes)
        PAS = paste(sep="_", "Jtim(", rd$jd[wearliest], ", hr=" , rd$hr[wearliest] ,
          ", mi=", rd$mi[wearliest], ",sec=", rd$sec[wearliest], ")")

        DEEtimes = YRsecdif(
          rd$jd[wearliest],rd$hr[wearliest],rd$mi[wearliest], rd$sec[wearliest],
          rd$jd,  rd$hr, rd$mi, rd$sec, rd$yr[wearliest],  rd$yr) 

        apickorg = paste(sep=",", rd$yr[wearliest], rd$jd[wearliest],rd$hr[wearliest],rd$mi[wearliest], rd$sec[wearliest])
        
        ##  pstas =  nh$STNS[ipick]

        apstas = paste(sep="", '"', paste(rd$stn, collapse='","'), '"')


        ##    pcomps =nh$COMPS[ipick]

        apcomps = paste(sep="", '"', paste(rd$comp, collapse='","'), '"')

        cat("", sep="\n")
        cat("", sep="\n")
        cat("##################", sep="\n")
        cat( paste(sep=" ", "orgtim=c(", apickorg , ")") , sep="\n")
        
        cat("", sep="\n")
        cat( paste(sep=" ", "stns=c(", apstas, ")") , sep="\n")
        cat( paste(sep=" ", "comps=c(", apcomps, ")") , sep="\n")

        cat( paste(sep=" ", "tims=c(", paste(DEEtimes, collapse=","), ")") , sep="\n")

        cat("", sep="\n")
        cat("##################", sep="\n")
        cat("", sep="\n")
        cat("Time Differences between picks:", sep="\n")
        
        cat(paste(dpick), sep="\n")

        cat("", sep="\n")
####  print(zloc$y[1:(zenclick-1)])  
####  print(ypick)     
####  print(ipick)
         cat("##################", sep="\n")
        cat("rd = scan(file='', what=list(jd=0,hr=0,mi=0,sec=0,yr=0,stn='',comp=''))" , sep="\n")
        write.table(file="", data.frame(rd), row.names =FALSE, col.names =FALSE )
        cat(" ", sep="\n")

        
        cat("GMT TIME: ", sep="\n")
        showdatetime(rd)

        cat(" ", sep="\n")
        
        PAS = paste(sep=" ", "Jtim(", rd$jd, ", hr=" , rd$hr , ", mi=", rd$mi, ",sec=", rd$sec, ")")
        cat("", sep="\n")
        cat(PAS, sep="\n")


        if(!is.null(nh$TZ))
          {
            rdlocal = recdate(jd=rd$jd, hr=rd$hr+nh$TZ, mi=rd$mi, sec=rd$sec , yr=rd$yr)
            cat(" ", sep="\n")
            
            cat(paste(sep=" ", "LOCAL TIMES, SHIFT=", nh$TZ) , sep="\n")
            showdatetime(rdlocal, AMPM=TRUE)
            
          }

      }
    else
      {
        cat("Pinfo WARNING: no pick or trace has been selected:", sep="\n")
        
      }
  
    
      g$zloc = list(x=NULL, y=NULL) 

    g$action="donothing"
    invisible(list(global.vars=g))
 
  }
#################################
#################################
TSHIFT<-function(nh, g)
  {
    #####  BUTTONDOC:TSHIFT:'Shift traces to line up with first pick'
    zenclick = length(g$zloc$x)
          if(zenclick>=2)
            {
              
              kix = legitpix(g$sel, g$zloc, zenclick)

              
              ypick =  kix$ypick
              ppick = kix$ppick
      
              dpick = c(0, diff(ppick))
              ipick = g$sel[ypick]

              print(paste(nh$STNS[ipick], nh$COMPS[ipick]))

              
              tshft = rep(0,times=length(nh$STNS))

              
              tshft[ipick] = ppick-ppick[1]
              
              ## print(data.frame(list(sta=nh$STNS, comp=nh$COMPS, tshft=tshft)))

              print(data.frame(list(sta=nh$STNS[ipick] , comp=nh$COMPS[ipick] , tshft=tshft[ipick] )))


              
              Tshift = list(name = nh$STNS[ipick], t=tshft[ipick])

              cat(file = "","\n\n")
              
              nam = "kshift"
              cat(file = "", paste(sep = "", nam, "=list()")   )
              cat(file = "","\n")
              cat(file = "", paste(sep = "", nam, "$name=c(\"", paste(format(Tshift$name), collapse = "\",\""), "\")"), fill = TRUE)
              cat(file = "", paste(sep = "", nam, "$t=c(", paste(format(Tshift$t), collapse = ","), ")"), fill = TRUE)
              cat(file = "","\n")
              
              g$ASHIFT = tshft
               g$BLAHSHIFT = tshft
            }
          else
            {

              g$ASHIFT = g$SHIFT.ORIG

            }

    g$zloc = list(x=NULL, y=NULL)
      g$action = "replot"
    invisible(list(NH=nh, global.vars=g))


}
#################################
#################################
#################################
RMS<-function(nh, g)
  {
    #####  BUTTONDOC:RMS:'Root Mean Square of selection'
    zenclick = length(g$zloc$x)
    sel = g$sel
    if(zenclick>=2)
      {
        kix = legitpix(g$sel, g$zloc, zenclick)
        ypick =  kix$ypick
        ppick = kix$ppick
        
        myinfo = list(yr=nh$info$yr, jd=nh$info$jd, hr=nh$info$hr, mi=nh$info$mi, sec=rep(0, times=length(nh$info$mi)))
        
        if(length(ypick)>0)
          {   ############   length(ypick) proceed only if have legitimate picks
            
            ipick = sel[ypick]
            npick = length(ypick)
            
            pairseq = seq(from=1, to=npick-1, by=2)
            
#####    Output1 = vector(length=length(pairseq))
            Output2 = vector(length=length(pairseq))
            for(iz in  pairseq)
              {   ###############   loop over pairs of picks
                i1 = ipick[iz]
################  this is the time in sec from the beginning of the trace
                asec = nh$info$sec[i1]+nh$info$msec[i1]/1000+nh$info$t1[i1]-nh$info$off[i1]+ppick[iz]
                if(npick<2)
                  {
                    bsec = asec+5
                  }
                else
                  {
                    iz1 = ipick[iz+1]
                    bsec = nh$info$sec[iz1]+nh$info$msec[iz1]/1000+nh$info$t1[iz1]-nh$info$off[iz1]+ppick[iz+1]
                  }
                
                rsig1 = nh$JSTR[[i1]]
                
                t1 = seq(from=0, length=length(rsig1), by=nh$dt[i1])
                
                which.time = which( t1>ppick[iz]  & t1< ppick[iz+1] )
                rwhich = range(which.time)
                
                rsig  = rsig1[ which.time ]
                
                rsig = rsig-mean(rsig)
                
                rms = sqrt( mean( rsig^2 ))
                
                cat(paste(sep=" ", "#########", iz, i1, format(ppick[iz]), format(ppick[iz+1]),
                          format(asec) , format(bsec), length(rsig), format(rms) ), sep="\n" )

                Output2[iz] = paste(sep=" ",
                         nh$STNS[i1],
                         nh$COMPS[i1] ,
                         myinfo$yr[i1],
                         myinfo$jd[i1],
                         myinfo$hr[i1],
                         myinfo$mi[i1],
                         format(asec),
                         format(bsec),
                         format(rms))
                
                dur = diff(c(asec, bsec) )
             
                if(is.null(dur)) dur = 0



               #### g$WPX =  pickhandler(i1=i1, ppick=ppick[iz], kzap=kzap, err=NA, ycol=ycol, NPX=g$NPX, g$WPX, nh)
               #### g$NADDPIX = g$NADDPIX+1
                
               #### g$NPX = g$NPX+1

               #### Nn = names(g$WPX)
               #### g$WPX =rbind(g$WPX, rep(NA, length(Nn)))
              }

            cat("############", sep="\n")
            cat( "OUTrms =scan(file=\"\", what=list(stn=\"\", comp=\"\", yr=0, jd=0, hr=0, mi=0, t1=0, t2=0, rms=0))", sep="\n" ) 
            for(iz in  pairseq)
              {
                cat(Output2[iz], sep="\n")
              }
            cat("\n" )
            cat("######", sep="\n")
          }
        else
          {
            print("not enough legitimate picks, need at least 2 or more")
            
          }
        
      }
    else
      {
        
       
         print("not enough legitimate picks, need at least 2 or more")
      }

     g$zloc = list(x=NULL, y=NULL) 
    g$action = "donothing"
     
    invisible(list(global.vars=g))

}

###############################
LocStyle<-function(nh, g)
  {
    #####  BUTTONDOC:LocStyle:'choose the locator style for picking in swig' 
    ###  choose the locator style for picking in swig

    g$ilocstyle = -1
    inum= c(-1, 0, 1, 2, 3)
    achoice = c("points", "abline", "segs(default)", "segs+abline", "segs+long-abline")


    P2 = RPMG::chooser(achoice, ncol=5, nsel=1, newdev=TRUE, STAY=FALSE,
      cols =rgb(1, .7, .7) , main="" , pch=21, cex=3,  col='red' , bg='blue' )

    i = which(P2==achoice)
    g$ilocstyle = inum[i]
    g$iloc
    g$action = "donothing"
    invisible(list(NH=nh, global.vars=g))


  }

