#' exportDRCDataToDataFrame
#'
#' Exports data into a standard data frame format for an IncucyteDRCSet
#'
#' @param idrc_set IncucyteDRCSet object
#' @param include_control Whether to include control sample as zero conc control
#' @param add_metadata Whether or not to merge IncucyteDRCSet metadata into the output
#'
#' @return IncucyteDRCSet object
#' @export
#'
#' @examples
#' pm_file <- system.file(file='extdata/example.PlateMap', package='IncucyteDRC')
#' test_pm <- importPlatemapXML(pm_file)
#' data_file <- system.file(file='extdata/example_data.txt', package='IncucyteDRC')
#' test_data <- importIncucyteData(data_file, metric='pc')
#'
#' test_list <- splitIncucyteDRCPlateData(test_pm, test_data, group_columns='growthcondition')
#'
#' print(test_list)
#'
#' test_idrc_set <- fitGrowthCurvesGrouped(test_list[[2]])
#' test_idrc_set <- fitGrowthCurvesIndividual(test_idrc_set)
#' test_idrc_set <- calculateDRCData(test_idrc_set, cut_time=100)
#' exportDRCDataToDataFrame(test_idrc_set)
#' exportDRCDataToDataFrame(test_idrc_set, include_control=TRUE)
#' exportDRCDataToDataFrame(test_idrc_set, include_control=TRUE, add_metadata=TRUE)
exportDRCDataToDataFrame <- function(idrc_set, include_control=FALSE, add_metadata=FALSE) {

    if(is.null(idrc_set$drc_data)) {
        stop('Need to calculate the dose response curve data first using calculateDRCData')
    }

    #get the sample data
    out_df <- idrc_set$drc_data %>%
                    dplyr::filter(samptype=='S') %>%
                    dplyr::transmute(sampleid, conc, value=cut_val)

    if(include_control) {
        #get the control data and wrangle it to look like zero conc data for each sample
        ctrl_df <- idrc_set$drc_data %>%
            dplyr::filter(samptype=='C') %>%
            dplyr::transmute(conc, value=cut_val)

        if(nrow(ctrl_df) ==0 )  stop('No control data - check your platemap configuration or set the include_control parameter to FALSE')

        ctrl_df <- ctrl_df %>%
            merge(y=unique(out_df$sampleid)) %>%
            dplyr::transmute(sampleid=y, conc=0, value)

        out_df <- dplyr::bind_rows(out_df, ctrl_df)
    }

    #add a replicate number
    out_df <- out_df %>%
        dplyr::group_by(sampleid, conc) %>%
        dplyr::mutate(idx=row_number()) %>%
        dplyr::ungroup()  %>%
        as.data.frame()

    #add metadata if needed
    if(add_metadata & is.data.frame(idrc_set$metadata)) {
        out_df <- merge(out_df, idrc_set$metadata)
    }

    return(as.data.frame(out_df))

}
