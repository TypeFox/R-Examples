`getVersion` <-
function(data=NULL) {
    if (is.SGP(data) && .hasSlot(data, "Version")) {
        return(list(SGP_Package_Version=c(data@Version[["SGP_Package_Version"]], as.character(packageVersion("SGP"))),
            Date_Prepared=c(data@Version[["Date_Prepared"]], date())))
    } else {
        return(list(SGP_Package_Version=as.character(packageVersion("SGP")), Date_Prepared=date()))
    }
} ### END getVersion function
