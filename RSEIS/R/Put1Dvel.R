`Put1Dvel` <-
function(vel, outfile)
  {
    #    v =  scan(file=infile, skip=2, list(zp=0, vp=0, ep=0, zs=0, vs=0, es=0), quiet =TRUE)
    #
    #  a1 = "# MODELK2      IASPEI 91 Velocity"
    #  a2 = "# P DEPTH   P VEL      PERR      S DEPTH    S VEL    SERR"

   write(file =  outfile, vel$descriptor)

   vel$vs = round(vel$vs*100)/100

   v = cbind(vel$zp, vel$vp, vel$ep, vel$zs, vel$vs, vel$es)
   
   write.table(file =  outfile, v , append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="     ")

  }

