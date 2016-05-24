"level.setCOP" <-
function(cop=NULL, para=NULL, getlevel=NULL, delu=0.001, lines=FALSE, ...) {
   zz <- level.curvesCOP(cop=cop, para=para, getlevel=getlevel, delt=NULL,
                         ploton=FALSE, plotMW=FALSE, lines=lines, delu=delu, ramp=FALSE, ...)
   return(zz)
}
