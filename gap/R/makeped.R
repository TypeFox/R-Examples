# 18-10-03 start to implement
# 19-10-03 worked but could not figure out loops
# 20-10-03 polished and tested with Abbas Parsian's homozygosity mapping pedigrees, add documentation
#
makeped<-function(pifile="pedfile.pre",pofile="pedfile.ped",auto.select=1,
                  with.loop=0,loop.file=NA,auto.proband=1,proband.file=NA)
{
  z<-.C("makeped",pifile=as.character(pifile),pofile=as.character(pofile),autoselect=as.integer(auto.select),
        withloop=as.integer(with.loop),loopfile=as.character(loop.file),
        autoproband=as.integer(auto.proband),probandfile=as.character(proband.file),PACKAGE="gap")
}
