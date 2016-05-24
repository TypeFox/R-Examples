#' populateIncucyteDRCSetMetadata
#'
#' Function to populate the metadata slot of an IncucyteDRCSet object.  Not exported, used at part of the
#' splitIncucyteDRCPlateData function.
#'
#' @param idrc_set IncucyteDRCSet object
#' @param group_columns Vector of columns that are present in the data frame to generate the groups.
#'
#' @return IncucyteDRCSet object


populateIncucyteDRCSetMetadata <- function(idrc_set, group_columns) {

    #create the metadata data frame which should be a single row
    metadata <- idrc_set$platemap %>%
        dplyr::select_(.dots=c('group_idx', group_columns)) %>%
        dplyr::mutate(metric=idrc_set$platedata$metric,
                      plateid=idrc_set$platedata$plateid) %>%
        dplyr::distinct() %>%
        as.data.frame()

    stopifnot(nrow(metadata)==1)

    outdata <- idrc_set
    outdata$metadata <- metadata

    return(outdata)

}
