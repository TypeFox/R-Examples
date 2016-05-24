convert2Rseis<-function(FLS, NEWDIR=".", kind = 1, Iendian="little" , BIGLONG=FALSE, NEWsta="", NEWcomp="" )
  {

############   takes an array  of binary file
############    converts to native R format for
####    later use and storage
    if(missing(kind)) { kind=1 }
    if(missing(Iendian)) { Iendian=1 }
    if(missing(BIGLONG)) { BIGLONG=FALSE}
    if(missing(NEWDIR)) {NEWDIR="." }
    if(missing(NEWsta)) {NEWsta=NULL }
    if(missing(NEWcomp)) {NEWcomp=NULL }

###  the output file is names after the original file name.
    ####  the list returned is DAT
    ###  so if this file is loaded, use DAT to extract the components

    ###  uassumes files contain only one trace at a time

##########   JGET.seis returns a list of lists of:
     ###   list(fn = fn, sta = thesta, comp = thecomp,
      ###          dt = dt, DATTIM = tstart, N = N, units = aunits,
      ###          amp = x)


######################  convert data to R-format

    for(i in 1:length(FLS))
      {

        f1 = FLS[i]

        a = JGET.seis(f1, kind = kind, Iendian=Iendian, BIGLONG=BIGLONG , PLOT = -1)

        DAT = a[[1]]

        bn = basename(DAT$fn)
        d2 = dirname(DAT$fn)

        if(!is.null(NEWsta))
          {
            DAT$sta = NEWsta[i]
            bn = paste(sep=".", bn,NEWsta[i])
          }

        if(!is.null(NEWcomp))
          {
            DAT$comp = NEWcomp[i]
             bn = paste(sep=".", bn,NEWcomp[i])
          }


        
        
        NEWbn = paste(sep=".", bn, "RDATA")

        fileout = paste(sep="/", NEWDIR, NEWbn)
        DAT$fn = fileout
        save(file=fileout, DAT)



      }



  }
