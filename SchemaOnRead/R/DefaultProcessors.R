##
## File:   DefaultProcessors.R
## Author: Michael J. North
## Date:   November 25, 2015
##

##
## Define the default processors list getter.
##
defaultProcessors <- function() {

  ## Define the default processors list.
  list(
    processDirectory,
    processCSVFile,
    processTextFile,
    processRDSFile,
    processXMLFile,
    processXLSandXLSXFile,
    processDIFFile,
    processODSFile,
    processBitmapFile,
    processGIFFile,
    processTIFFFile,
    processNetCDandH5FFile,
    processPajekFile,
    processHTMLFile,
    processARFFFile,
    processEPIINFOFile,
    processSPSSFile,
    processSASFile,
    processStataFile,
    processDefaultFile
  )

}
