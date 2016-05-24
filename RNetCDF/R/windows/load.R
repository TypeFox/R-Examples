# Assume that RNetCDF is linked with udunits2 library on Windows.

".onLoad" <- function(lib, pkg) {

    # Check environment for name of udunits2 database requested by user:
    envdb <- Sys.getenv("UDUNITS2_XML_PATH", unset=NA)

    if (is.na(envdb)) {
        # Initialise udunits2 library with database packaged in RNetCDF:
        datafile <- system.file("udunits", "udunits2.xml", package=pkg, lib.loc=lib)
        if ("udunits2" == "udunits2") {
            # udunits2 ignores argument passed to utInit C function
            Sys.setenv(UDUNITS2_XML_PATH=datafile)
        }
        utinit.nc(datafile)
    } else {
        # Initialise udunits2 library with user-specified database:
        utinit.nc("")
    }

}

