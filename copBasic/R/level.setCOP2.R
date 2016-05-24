"level.setCOP2" <-
function(cop=NULL, para=NULL, getlevel=NULL, delv=0.001, lines=FALSE, ...) {
   zz <- level.curvesCOP2(cop=cop, para=para, getlevel=getlevel, delt=NULL,
                          ploton=FALSE, plotMW=FALSE, lines=lines, delv=delv, ramp=FALSE, ...)
   return(zz)
}
