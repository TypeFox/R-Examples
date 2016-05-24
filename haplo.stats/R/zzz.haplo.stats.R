#$Author: sinnwell $
#$Date: 2013/01/14 19:33:16 $
#$Header: /projects/genetics/cvs/cvsroot/haplo.stats/R/zzz.haplo.stats.R,v 1.4 2013/01/14 19:33:16 sinnwell Exp $
#$Locker:  $
#$Log: zzz.haplo.stats.R,v $
#Revision 1.4  2013/01/14 19:33:16  sinnwell
#small changes for 1.5.9
#
#Revision 1.3  2012/12/28 21:39:36  sinnwell
#version with update zzz file
#
#Revision 1.2  2011/12/05 20:56:10  sinnwell
#final manual changes, updated test suite
#
#Revision 1.1  2004/07/02 14:18:11  sinnwell
#Initial revision
#
#Revision 1.5  2003/10/06 15:45:17  sinnwell
#dyn.load haplo.stats, not haplo.score
#
#Revision 1.4  2003/08/26 16:41:23  sinnwell
#change License statement
#
#Revision 1.3  2003/03/20 22:26:56  sinnwell
#remove source() command
#
#Revision 1.2  2003/03/17 16:45:41  sinnwell
#autoload data in hla.demo.q in .First.Lib()
#
#Revision 1.1  2003/03/06 23:25:55  sinnwell
#Initial revision
#

#.onLoad <- function(lib, pkg) {
#   library.dynam("haplo.stats", pkg, lib)
#}

#.onAttach <- function(lib, pkg) {
#   library.dynam("haplo.stats", pkg, lib)
#}

##.First.lib <- function(lib, pkg) {
##
##library.dynam("haplo.stats", pkg, lib)
##}


##.Last.lib <- function(libpath) {
##  library.dynam.unload("haplo.stats", libpath)
##}

.onUnload <- function(libpath) {
  library.dynam.unload("haplo.stats", libpath)
}
