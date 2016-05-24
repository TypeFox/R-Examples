`NEWPLOT.WPX` <-
function(t0, STNS, COMPS, YPX, FILL=FALSE, FORCE=TRUE, cex=cex, srt=srt)
  {
    
    ##     YPX[[NPX]] =  list( pick=c(NH$info$yr[ipick], NH$info$jd[ipick], NH$info$hr[ipick], NH$info$mi[ipick], asec), kind="SHAPE", sta= STNS[ypick], comp=COMPS[ypick], PULS=PULS)

    ### match(JH$COMPS, "V")
    ###  match(JH$STN, "MAS")&match(JH$COMPS, "V")


    if(missing(FILL)) { FILL=FALSE }
    if(missing(FORCE)) { FORCE=TRUE }
    if(missing(cex))    { cex=1 }
    if(missing(srt))    { srt=0 }

    
#########  if FORCE==true then if there is no match, plot picks anyway on first trace
   if(length(STNS)<=0) { return(0) }
    
    du = 1/length(STNS)

    
 ###  print(YPX)
  ### print(paste(sep=' ', "PLOT.WPX NPX=", length(YPX), "du=", du))

    ### print(cbind(STNS, COMPS))

    for(i in 1:length(YPX))
      {

        qpix = YPX[[i]]

        if(is.null(qpix$dispcomp)) qpix$dispcomp = qpix$comp

        if(is.null(qpix$kind)) { qpix$kind = NA }

        
        if( any(as.character(qpix$kind)=="SHAPE") )  { next; }

        if(!is.null(STNS) )
          {
            msta = which(STNS %in% qpix$sta )
            
          }

          if(!is.null(COMPS) )
          {
            mcomp = which(COMPS %in% qpix$dispcomp)
          }
      
        
     ##    imatch = unique( c(msta , mcomp) )

        im = which(  msta  %in%  mcomp)

        imatch = msta[im]

        if(length(imatch)<1 )
          {
            if(FORCE==TRUE)
              {
                print(paste(sep=' ', "ERROR: No Match in PLOT.WPX", i, imatch, qpix$sta, qpix$comp))
                imatch=1
              }
            else
              {
                next
              }
          }
        
      ### print(paste(sep=' ', "PLOT.WPX", i, imatch, qpix$sta, qpix$comp))

              
        ypixA = (length(STNS)-imatch)*du
        ypixB =  ypixA+du

        
        x1 = secdif(   t0$jd, t0$hr, t0$mi, t0$sec, qpix$pick[2],  qpix$pick[3],  qpix$pick[4], qpix$pick[5])

        
      ###  print(paste(sep=' ', "PLOT.WPX", i, x1, ypixA, x1, ypixB, qpix$sta, qpix$comp, qpix$kind, imatch))

       ###  if(x1>0 & x1 <3600)
        #####   {
       #####      print(paste(sep=' ', "PLOT.WPX", i, x1, ypixA, x1, ypixB))
       #####    }

        
       if(FILL==TRUE)
         {
           #############  draw the pick down the whole window
             segments(x1, 0, x1, 1, col=grey(.8) )
         }

               #############   draw the pick  on the designated panel
        segments(x1, ypixA, x1, ypixB, col=qpix$col)
        text(x1, ypixB, labels=qpix$kind, col=qpix$col, pos=4, cex=cex, srt=srt)

 
         #############    if a duration is given, mark it with a bar and a hook
       if(length(qpix$dur)>0)
         {
             segments( x1, ypixB, x1+qpix$dur, ypixB, col=qpix$col  )
             segments( x1+qpix$dur, ypixB, x1+qpix$dur, ypixB-.2*du, col=qpix$col  )

             
         }
       #############    if an err is given, mark it with a small error bat
       if(length(qpix$err)>0)
         {
           if(is.null(qpix$ecol)) qpix$ecol = qpix$col
             segments( x1-qpix$err, ypixB-.1*du, x1+qpix$err, ypixB-.1*du, col=qpix$ecol  )
           
             segments( x1+qpix$err, ypixB-.05*du, x1+qpix$err, ypixB-.15*du, col=qpix$ecol  )
             segments( x1-qpix$err, ypixB-.05*du, x1-qpix$err, ypixB-.15*du, col=qpix$ecol  )

             
         }
       
        
      }
    
  }

