`readFrameFromFITS` <-
function (file, hdu = 1) 
{
### Use is analogous to readFITS, but it returns a data frame.
###
### Takes:
  ## FITS file name: file
  ## Index number of header and data unit in FITS file: hdu
### Returns:
  ## Data frame with bintable data
### Requires/Used by:
  ## Requires .nameFrameColumn.r
###
### Written by Eric H. Neilsen, Jr., FNAL
###   Modified set data frame properly at line 28, 2009-04-08
###
    ## Get data
    fitsContent <- readFITS(file, hdu)
    ## Make initial data frame, set column name
    dataFrame <- data.frame(fitsContent$col[[1]])
    dataFrame <- .nameFrameColumn(dataFrame, fitsContent$col[[1]], 
        fitsContent$colNames[1])
    ## Add subsequent columns of data and names
    for (i in 2:length(fitsContent$colNames)) {
        newFrame <- data.frame(fitsContent$col[[i]])
        newFrame <- .nameFrameColumn(newFrame, fitsContent$col[[i]], 
            fitsContent$colNames[i])
        dataFrame <- data.frame(c(dataFrame, newFrame))
    }
    ## Return
    dataFrame
}

`.nameFrameColumn` <-
function (nframe, datacol, colname) 
{
### Function makes column names for data frame
### (preceeded by . to hide it from users in package build)
### 
### Takes:
  ## Data frame: nframe
  ## New data for column: datacol
  ## New data column name: colname
### Returns:
  ## Data frame name vector for column
### Requires/Used by:
  ## Used by readFrameFromFITS.r
###
### Written by Eric H. Neilsen, Jr., FNAL
### 
    if (length(dim(datacol)) < 1) {  # flat file names
        names(nframe) <- colname
    } else {                         # sequential names for '3D' table entries
        names(nframe) <- sprintf("%s.%d", colname,
            c(1:length(datacol[1,])))
    }
    ## Return
    nframe
}

