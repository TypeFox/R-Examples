reload.constants <- function(xls.path,rda.path=xls.path)
  {

    
    if(!require("gdata"))
      stop("'gdata' package must be installed to run this function.")

    read.xls <- gdata::read.xls
    
    "critical" <- read.xls(
                           xls=file.path(xls.path,"SII_Constants.xls"),
                           sheet=1,
                           skip=5, 
                           header=TRUE, 
                           nrows=21, 
                           check.names=FALSE,
                           row.names=1
                           )
                           
    save(list="critical", file=file.path(rda.path, "critical.rda") )
    
    "equal" <- read.xls(
                                     xls=file.path(xls.path,"SII_Constants.xls"),
                                     sheet=2,
                                     skip=5, 
                                     header=TRUE, 
                                     nrows=17, 
                                     check.names=FALSE,
                                     row.names=1
                                     )
    save(list="equal", file=file.path(rda.path, "equal.rda") )

    "onethird" <- read.xls(
                                   xls=file.path(xls.path,"SII_Constants.xls"),
                                   sheet=3,
                                   skip=5, 
                                   header=TRUE, 
                                   nrows=18, 
                                   check.names=FALSE,
                                   row.names=1
                                   )
    save(list="onethird", file=file.path(rda.path, "onethird.rda") )

    "octave" <- read.xls(
                         xls=file.path(xls.path,"SII_Constants.xls"),
                         sheet=4,
                         skip=5, 
                         header=TRUE, 
                         nrows=6, 
                         check.names=FALSE,
                         row.names=1
                         )

    save(list="octave", file=file.path(rda.path, "octave.rda") ) 

    "sic.critical" <- read.xls(
                               xls=file.path(xls.path,"SII_Constants.xls"),
                               sheet=5,
                               skip=3, 
                               header=TRUE, 
                               nrows=21, 
                               check.names=FALSE,
                               row.names=1
                               )
                           
    save(list="sic.critical", file=file.path(rda.path, "sic.critical.rda") )

    "sic.onethird" <- read.xls(
                               xls=file.path(xls.path,"SII_Constants.xls"),
                               sheet=6,
                               skip=3, 
                               header=TRUE, 
                               nrows=18, 
                               check.names=FALSE,
                               row.names=1
                               )
    save(list="sic.onethird", file=file.path(rda.path, "sic.onethird.rda") )

    "sic.octave" <- read.xls(
                         xls=file.path(xls.path,"SII_Constants.xls"),
                         sheet=7,
                         skip=3, 
                         header=TRUE, 
                         nrows=6, 
                         check.names=FALSE,
                         row.names=1
                         )

    save(list="sic.octave", file=file.path(rda.path, "sic.octave.rda") ) 

    ls()
  }
