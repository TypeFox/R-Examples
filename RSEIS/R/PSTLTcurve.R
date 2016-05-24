`PSTLTcurve` <-
function(y, dt=0.008, fwlen =  125,  bwlen  = 125, perc=0.05, stretch=1000, MED=255, PLOT=FALSE)
{
###  do automatic picking on a trace
###  scale y so it can be used as an integer array

     if(missing(dt))  {  dt=0.008 }

  if(missing(fwlen))  {   fwlen =  1000 }
  if(missing(bwlen))  {   bwlen =  1000 }
  if(missing(stretch))  {  stretch =  1000 }

  if(missing(MED))  {  MED = 255 }
  if(missing(perc))  {  perc = 0.05  }
     
   if(missing(PLOT))  {  PLOT=FALSE }
 
  lx = length(y)

     if(lx<(fwlen+bwlen+MED))
       {

         return(list(ind=1, eye=1, rat=1))

       }

     
  r = range(y)

  rat = rep(0, length(y))

  logflg = 0


     #######  new preconditioning

     ey = envelope(y-mean(y))

     
     ##########apply a robust smoothing filter (running median)
     #### print(paste(sep=' ', "################ in PSTLT ", lx)) 
     s = runmed(ey, MED, algorithm ="Stuetzle", endrule ="constant")    
     rs = range(s)
     s = stretch*(0.5+(s-rs[1])/(rs[2]-rs[1]))


     quack = .C("CALL_DFBRAT",PACKAGE = "RSEIS",
       as.double(s),  as.double(rat),
       as.integer(lx),as.integer(fwlen), as.integer(bwlen), as.integer(logflg) )


####   plot.ts(quack[[2]], ylab="ratio")
      ###  this is the maximum of the ratio curve:
     therat = quack[[2]]
     ix = which.max(therat)
     mix = ix

     if(1 == 1 )
       {
     prat = peaks(therat)
     ####print(therat[prat])
     mrat = mean(c(mean(therat), therat[ix]), na.rm=TRUE)
      ####print(mrat)
     rx =   which(prat)
     #### print(rx)
     gx = rx[therat[rx]>=mrat]
     ####print(gx)
     mix = min(gx, na.rm=TRUE)
   }
     if(is.na(mix)) { mix = ix }
     
     ###  need to look for changes in s
     ###  s should be flat in the noise region and then start rising
     ###   where the signal starts - need to capture that point.
     
     s2 = runmed(ey, 4*MED+1, algorithm ="Stuetzle", endrule ="constant")    
     rs = range(s2)
     s2 = stretch*(0.5+(s2-rs[1])/(rs[2]-rs[1]))
     ex = seq(1,length(s2))

     M = fwlen
     el = s2[1:M]
     
     if(length(el)<1)
       {
         print(paste(sep=' ',"problem in PSTLTcurve", M) )
         return(NA)
       }

     if(length(which(is.na(el)))>1)
       {
         print(paste(sep=' ',"problem in PSTLTcurve", M) )
         return(NA)
       }
     
     m2 = sqrt(var(el))
     mel = jstats(el)
     ### abline(h=mel$mean, col=3)
     ### abline(h=mel$mean+6*mel$std, col=4)
 
     m1 = mean(el)
     m3 = max(s2)
                     ####  source("/home/lees/Progs/R_stuff/autopix.R")

     ###  slowly lower the perc value and calculate eye each time.
     ###   when done take median value of eye
    
     vim = pretty(c(.5,  perc), n=10)
     KAPS = (m1+vim*(m3-m1))
     
     jout = vector()
     
     for(kj in 1:length(vim))
       {

         perc1 = vim[kj]
         

         kap = (m1+perc1*(m3-m1))
         ef = s2>kap&ex>M

         if(length(ex[ef])>2)
           {
             eye1 = min( ex[ef])
                 ###############  if eye1 is greater than mix....look in a different way
             if(eye1 >(mix))
               {
                 M = max(M, eye1-2*M)
                 el = s2[1:M]
                 
                    ### abline(h=mel$mean, col=3)
                    ### abline(h=mel$mean+6*mel$std, col=4)
                 m1 = mean(el)
                 m3 = max(s2)
                 
                 kap = (m1+perc1*(m3-m1))

                 ### kap = max(el)
                  ### kap = kap+0.1*(kap-mel$mean)
                 ef = s2>kap&ex>M

               }
           }
  

         if(length(ex[ef])>2)
           {
             eye = min( ex[ef])
           }
         else
           {
             eye = ix
           }
         jout=c(jout, eye)
       }

     ###  here we want the smalles jout that is not an outlier

     ###  Key = stats(jout)

     ###  Key = quantile(jout, c(1,3)/4 )

     
     eye = round( quantile(jout, c(1)/5 )  )

     
     

     
     SNR = 0
     isig1 = ix+1
     isig2 = min(c((ix+fwlen), length(ex)))
     
     inois1 = ix-1
     inois2 = max(c((ix-fwlen), 1))

     nois = y[inois1:inois2]
     if(length(nois)>2)
       {
         nois2 = sum(nois^2)
         if(nois2>0)
           {
             sig2 = sum(y[isig1:isig2]^2)
             if(sig2>=0)
               {
                 SNR = sig2/nois2
               }
           }
       }

     

     if(PLOT==TRUE)
       {
         opar <- par(no.readonly = TRUE)
         x = seq(1:length(y))
         
         par(mfrow=c(3,1))
         par(mai=c(0.5, .7, 0.2, 0.5) )
         
         plot(x, y, type='n')
         abline(v=ix, col=2)
         abline(v=eye, col=4)
         abline(v=mix, col=3)
         lines(x,y)

         mtext(c("ix", "eye", "mix"), line=.2, col=c(2,4,3), at=c(ix, eye, mix))
         
         plot(x, s, type='n')
         abline(v=ix, col=2)
         abline(v=mix, col=3)        
         abline(h=KAPS, lty=2, col=rgb(.9,.9,.9)  )
         abline(h=kap, lty=2, col=3)
         
         lines(x,s)
         lines(x,s2, col=5)

         abline(v=eye, col=4)
         abline(v=eye1, col=rgb(1, .5, .5) )

         
         plot(x, therat, ylab="ratio", type='n')
         abline(v=ix, col=2)
         abline(v=eye, col=4)
         abline(v=mix, col=3)         
         lines(x, therat)
         title(paste(sep=' ', "Ratio Curve=", fwlen, bwlen, " SNR=", format.default(SNR, digits=5)))
       mtext(c("ix", "eye", "mix"), line=.5, col=c(2,4,3), at=c(ix, eye, mix))

      u = par("usr")
         L = 0
         z = u[3]+0.95*(u[4]-u[3])
         
     ######    print(paste(sep=' ', "psegments 1", L, z, L+bwlen, z))
         if(all(is.numeric(c(L, z, L+bwlen, z))))
           {
             segments(L, z, L+bwlen, z, col=2, lwd=3, xpd=TRUE)
           }
         
    ######    print(paste(sep=' ', "psegments 2",L+bwlen+100, z,  L+bwlen+100+fwlen, z))
         if(all(is.numeric(c(L+bwlen+10, z,  L+bwlen+10+fwlen, z))))
           {
             segments(L+bwlen+10, z,  L+bwlen+10+fwlen, z, col=4, lwd=3, xpd=TRUE)
           }
         par(opar)
 
       }


     if(is.na(ix)  | is.null(ix) )  ind=NA
     if(is.na(eye) | is.null(eye) ) eye=NA
     if(is.na(mix) | is.null(mix) ) mix=NA
     
  return(list(flag=1, ind=ix, eye=eye, mix=mix, SNR=SNR,  s2=s2, rat=therat))
}

