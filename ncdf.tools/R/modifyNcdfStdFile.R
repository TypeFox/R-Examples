modifyNcdfStdFile <-  function(
    file.con  ##<< character string (i.e. file name) or file connection of the file to change.
)
    ##title<< Standardize file name and missing value attribute of a ncdf file.
    ##description<< Wrapper function around modifyNcdfStdFile and modifyNcdfSetMissing
    ## to standardize the name and the missing value attributes of a ncdf file.
    
{

    if (inherits(file.con,  'character')) {
        if (!file.exists(file.con))
            stop(paste('File ', file.con, ' not existent!'))
        file.con <- open.nc(file.con,  write =  TRUE)
        close.file = TRUE
    } else {
        close.file = FALSE
    }
    modifyNcdfStdNames(file.con)
    modifyNcdfSetMissing(file.con, readNcdfVarName(file.con))
    if (close.file)
        close.nc(file.con) 
}
