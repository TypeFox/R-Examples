apc.frame <-
function( a.lab,
         cp.lab,
          r.lab,
         rr.lab = r.lab / rr.ref,
         rr.ref = r.lab[length(r.lab)/2],
          a.tic = a.lab,
         cp.tic = cp.lab,
          r.tic = r.lab,
         rr.tic = r.tic / rr.ref,
        tic.fac = 1.3,
          a.txt = "Age",
         cp.txt = "Calendar time",
          r.txt = "Rate per 100,000 person-years",
         rr.txt = "Rate ratio",
       ref.line = TRUE,
            gap = diff(range(c(a.lab,a.tic)))/10,
       col.grid = gray( 0.85 ),
          sides = c(1,2,4) )
  {
cp  <- min( c(cp.tic,cp.lab) ) - max( c(a.lab,a.tic) ) - gap
xl  <- c(min( c(a.lab,a.tic) ),
         max( c(a.lab,a.tic) ) + gap + diff( range( c(cp.lab,cp.tic) ) ) )
yl  <- range( c(r.lab,r.tic) )
rrtck <- outer( c(0.5,1,1.5,2:9), 10^(-5:5), "*" )
# Empty plot frame
plot( NA,
      xlab="", xlim=xl, xaxt="n", xaxs="i",
      ylab="", ylim=yl, yaxt="n", yaxs="i", log="y" )
# Grid lines
abline( h=c(r.tic,outer( c(0.5,1,1.5,2:9), 10^(-5:5), "*" )),
        v=c(a.tic,cp.tic - cp), col=col.grid )
# Reference line for the RR=1
if ( ref.line )
segments( min(c(cp.lab,cp.tic))-cp, rr.ref,
          max(c(cp.lab,cp.tic))-cp, rr.ref )
# Close it nicely off:
box()
# Axis construction (tickmarks, labels and annotation)
if ( 1 %in% sides )
   {
   axis( side=1, at=a.lab )
   axis( side=1, at=a.tic, labels=NA, tcl=par("tcl")/tic.fac )
   axis( side=1, at=cp.lab - cp, labels=cp.lab )
   axis( side=1, at=cp.tic - cp, labels=NA, tcl=par("tcl")/tic.fac )
   axis( side=1, at=mean( a.lab ), labels=a.txt, line=1, tcl=0 )
   axis( side=1, at=mean( cp.lab - cp ), labels=cp.txt, line=1, tcl=0 )
   }
if ( 2 %in% sides )
   {
   axis( side=2, at=r.lab, labels=paste( r.lab ) )
   axis( side=2, at=r.tic, labels=NA, tcl=par("tcl")/tic.fac )
   mtext( side=2, r.txt, line=2.5, las=0 )
   }
if ( 3 %in% sides )
   {
   axis( side=3, at=a.lab )
   axis( side=3, at=a.tic, labels=NA, tcl=par("tcl")/tic.fac )
   axis( side=3, at=cp.lab - cp, labels=cp.lab )
   axis( side=3, at=cp.tic - cp, labels=NA, tcl=par("tcl")/tic.fac )
   axis( side=3, at=mean( a.lab ), labels=a.txt, line=1, tcl=0 )
   axis( side=3, at=mean( cp.lab - cp ), labels=cp.txt, line=1, tcl=0 )
   }
if ( 4 %in% sides )
   {
   axis( side=4, at=c(rr.ref,rr.lab*rr.ref), labels=paste( c(1,rr.lab) ) )
   axis( side=4, at=c(rr.ref,rr.tic*rr.ref), labels=NA, tcl=par("tcl")/tic.fac )
   mtext( side=4, rr.txt, line=2.5, las=0 )
   }
# Return the offset for the cohort/period and the RR-factor.
options( apc.frame.par = c("cp.offset"=cp,"RR.fac"=rr.ref) )
invisible( c("cp.offset"=cp,"RR.fac"=rr.ref) )
  }
