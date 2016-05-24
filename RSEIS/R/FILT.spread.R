`FILT.spread` <-
function(x,y, dt, fl=fl, fh=fh, sfact=1, WIN=NULL, PLOT=TRUE, TIT=NULL, TAPER=0.05, POSTTAPER=0.05  )
  {
    ####  filter sweep 
    if(missing(fl))
      { fl = c(.1, 1, 2, 3, 4) }
    if(missing(fh))
      {  fh = c(1,2,3,4,5) }
    if(missing(sfact))  { sfact = 1 }
    if(missing(WIN))  { WIN=NULL }
    if(missing(PLOT))  { PLOT=TRUE }
        if(missing(TIT))  { TIT=NULL }
        if(missing(TAPER))  { TAPER=NULL }
        if(missing(POSTTAPER))  { POSTTAPER=0.05 }


    
    n=length(fl)
    ##  graphics.off()
   ##  par(mfrow=c(n+1, 1))
   ##  par(mai=c(0.0, .7, 0.1, 0.5))


    oky = !is.na(y)
    yrng = range(y[oky]-mean(y[oky]))
    
    FR = matrix(nrow=length(y), ncol=n+1)
    Notes = as.vector(1:(n+1))
    
    
    for(i in 1:n)
      {

        gy = y

        if(fl[i]>fh[i])
          {
            print(paste(sep=' ', "Warning on Filter definition: FL=", fl[i], " FH=", fh[i], "HZ"))
          }

        if(fh[i]>1/(2*dt))
          {
            print(paste(sep=' ', "Warning on Filter definition: FL=", fl[i], " FH=", fh[i], "HZ", "NYQ=", 1/(2*dt)))
          }

        if(!is.null(TAPER))
          {
            ###  print("applying tape at ", TAPER)
           
             tapy = applytaper(y, p=TAPER)
            y = tapy
          }
        
        fy = butfilt(y[oky],fl[i], fh[i], dt, "BP", "BU" )

        
          if(!is.null(POSTTAPER))
          {
            ###  print("applying tape at ", TAPER)
            fy = fy - mean(fy)
             ftapy = applytaper(fy, p=POSTTAPER)
            fy  =  ftapy 
          }
        
        gy[oky] = fy
       
        FR[,i] = gy

        khigh = format.default(fh[i], digits=3)
        lhigh  = "Hz"
        if(fh[i]<1)
          {
            khigh = format.default(1/fh[i], digits=3)
              lhigh  = "s"
          }
        klow = format.default(fl[i], digits=3)
        llow  = "Hz"
        
        if(fl[i]<1)
          {
            klow = format.default(1/fl[i], digits=3)
             llow  = "s"
          }
        Notes[i] = paste(sep=' ', "BP FILTER",klow,llow , "to",khigh , lhigh)
      }
    FR[,n+1] = y-mean(y[oky])
    Notes[n+1] =paste(sep=' ',"Unfiltered")

    if(PLOT==TRUE)
      {
        if(is.null(WIN)==FALSE)
          {
            
            PLOT.MATN( FR, tim=x, WIN=WIN, dt=dt, sfact=sfact, notes=Notes, add=1)
          }
        else
          {
            PLOT.MATN( FR, tim=x, dt=dt, sfact=sfact, notes=Notes, add=1)
          }
        if(!is.null(TIT)) { title(main=TIT) }
      }
    
    invisible(list(FMAT=FR, Notes=Notes) )
    
}

