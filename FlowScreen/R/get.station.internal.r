#' Get station information for USGS or WSC hydrometric stations
#' 
#' @param stnID Character string of station ID.
#' @return Returns a list of station information
#' @author Jennifer Dierauer

## uses HYDAT_list - internal dataset in sysdata.rda

get.station.internal <- function (stnID) {

    rhbn <- NULL
    for (istn in 1:length(WSC.site.info$Station)) {
    	if (WSC.site.info$Station[istn] == stnID) 
    		{st_info <-(WSC.site.info[istn,1:17])
             
             if (WSC.site.info$RHBN[istn]==TRUE) (rhbn="*")
             
    		 st_info[18] <-paste(WSC.site.info[istn,1]," - ", WSC.site.info[istn,2], 
                                 " - ", WSC.site.info[istn,4], rhbn, sep="")
    		}
    }
    names (st_info) [18] <-"Station_lname"
    return (st_info)

}
