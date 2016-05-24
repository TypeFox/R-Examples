`DISPLACE.SEISN`<-
function(TH, sel=1:length(TH$JSTR), inst=1, Kal=Kal, waterlevel=1.e-8, FILT=list(ON=FALSE, fl=1/30, fh=7.0, type="HP", proto="BU") )
{
  ########  convert the seismic to displacement
  if(missing(Kal)) {   Kal = PreSet.Instr() }
  if(missing(inst)) {   inst = rep(1,length=length(TH$JSTR))  }

  if(missing(sel)) { sel = 1:length(TH$JSTR) }
    if(missing(FILT)) { FILT = list(ON=FALSE, fl=1/30, fh=7.0, type="HP", proto="BU")  }
 
  Calibnew = c(1,1.0, 0.0 )

  
if(is.logical(sel)) { sel = which(sel) }

  
 for(i in 1:length(sel))
   {
     ii = sel[i]
     dt = TH$dt[ii]
     ins = inst[ii]
     if(ins==0) next
     
     y = TH$JSTR[[ii]]

   
##############################    for now, do not analyse wigs  that have NA's
    if(any(is.na(y)))
       {

      ####  since some of the data does not
         ####    need to fill in and/delete

         TH$units[ii] = "NA"
         next

       }
     

     
     y = y-mean(y, na.rm =TRUE)
     y = detrend(y)
     y = applytaper(y)
     dy  = deconinst(y, dt, Kal, ins, Calibnew, waterlevel=waterlevel)
     
     ty = applytaper(dy-mean(dy, na.rm =TRUE), p=0.05)
     tapy = detrend(ty)
     fy = tapy-mean(tapy, na.rm =TRUE)

     #####  use this to get micro meters
     ## amp = fy*1000000
     amp = fy
     damp  = trapz(amp, dt)

     if(is.null(damp)) { TH$units[ii] = "NA"; next }

     if(FILT$ON==TRUE)
       {
         fy = butfilt(damp, FILT$fl, FILT$fh , dt, FILT$type , FILT$proto )
       }
     else
       {
         fy =damp
       }
     
     TH$JSTR[[ii]] = fy
      TH$units[ii] = "m"
   }

    proc = paste(sep=" ", "DISPLACE.SEISN", "FILT", FILT$ON, FILT$fl, FILT$fh, FILT$type, FILT$proto) 

  TH$process = c(TH$process, proc)
 invisible(TH)
}

