convertseis2R<-function(fn, destdir=".", kind=1, Iendian=1, BIGLONG=FALSE)
  {
########   convert a data set to native R binary format
###  this code reads in seismic data in a variety of formats
###  and saves them as binary R data files for fast loading and processing
###  fpath = "/Users/lees/Site/CHAC/DATA"
    ##  lf =  list.files(path = fpath, pattern = "20", all.files = TRUE,
    ##                 full.names = TRUE, recursive = TRUE,
    ##                 ignore.case = FALSE, include.dirs = FALSE)

    for(i in 1:length(fn))
      {
        fn2 = fn[i]
        ###  change the suffix from SAC or SEGY to RDATA
        bn = basename(fn2)
        bn2 = unlist( strsplit(bn, split="\\.") )
        bn3 = paste(collapse=".", bn2[1:(length(bn2)-1)])
        bn4 = paste(sep=".", bn3, "RDATA")
        fndest = paste(sep="/", destdir, bn4)
        ####  read in data
        DAT = RSEIS::JGET.seis(fn2, kind = kind, Iendian = Iendian, BIGLONG = BIGLONG, 
          HEADONLY = FALSE, PLOT = -1)
        ####  save in R format
        save(file=fndest, DAT)
      }
  }

