`STLTcurve` <-
function(y, dt=0.008, fwlen =  125,  bwlen  = 125, stretch=1000, MED=255, PLOT=FALSE)
{
###  do automatic picking on a trace
###  scale y so it can be used as an integer array

     if(missing(dt))  {  dt=0.008 }

  if(missing(fwlen))  {   fwlen =  1000 }
  if(missing(bwlen))  {   bwlen =  1000 }
  if(missing(stretch))  {  stretch =  1000 }

  if(missing(MED))  {  MED = 255 }
     
   if(missing(PLOT))  {  PLOT=FALSE }
 
  lx = length(y)
  r = range(y)

  rat = rep(0, length(y))

  logflg = 0

     ###   precondition the seismogram

##########################################################
############     ####  this is old way of doing things:
#############s =  abs(y-mean(y))
############
#############rs = range(s)
############
#############s = 10*(0.5+(s-rs[1])/(rs[2]-rs[1]))
##########################################################

##########################################################
############     ###  this is how it is done for picking in snaps:
############ ###  s = abs(10000*(y))

############ ###    s = 10*s
############ ###   s[s>0] = s[s>0] +0.5autopix.R
############ ###   s[s<0] = s[s<0] -0.5
##########################################################

     #######  new preconditioning

     ey = envelope(y-mean(y))
     ##########apply a robust smoothing filter (running median)
    ####       print(paste(sep=' ', "################ in STLT ", lx)) 
     s = runmed(ey, MED, algorithm ="Stuetzle")    
     rs = range(s)
     ess = stretch*(0.5+(s-rs[1])/(rs[2]-rs[1]))     

     s = ess
     ####  here we need to rectify problems associate with very large variations

   ####  kappa = mean(s)+2*sd(s)

  ####   s[s>kappa] = kappa

     ###  need to look for changes in s
     ###  s should be flat in the noise region and then start rising
     ###   where the signal starts - need to capture that point.
     

  quack = .C("CALL_DFBRAT",PACKAGE = "RSEIS",
    as.double(s),  as.double(rat),
    as.integer(lx),as.integer(fwlen), as.integer(bwlen), as.integer(logflg) )


     ####   plot.ts(quack[[2]], ylab="ratio")
      
  ix = which.max(quack[[2]])


     if(PLOT==TRUE)
       {
         opar <- par(no.readonly = TRUE)
         x = seq(1:length(y))
         
         par(mfrow=c(3,1))
         plot(x, y, type='n')
         abline(v=ix, col=2)
         lines(x,y)
         
         plot(x, s, type='n')
         abline(v=ix, col=2)
         lines(x,s)
         
         plot(x, quack[[2]], ylab="ratio", type='n')
         abline(v=ix, col=2)
         lines(x, quack[[2]])
         title(paste(sep=' ', "Ratio Curve=", fwlen, bwlen))

      u = par("usr")
         L = 0
         z = u[3]+0.95*(u[4]-u[3])
      segments(L, z, L+bwlen, z, col=2, lwd=3)
      segments(L+bwlen+1000, z,  L+bwlen+1000+fwlen, z, col=4, lwd=3)
         invisible( par(opar))
 
       }

  return(list(ind=ix, rat=quack[[2]] , ess=ess ))
}

