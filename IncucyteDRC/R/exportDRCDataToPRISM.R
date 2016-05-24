#' exportDRCDataToPRISM
#'
#' Exports data in PRISM format for an IncucyteDRCSet
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
#' exportDRCDataToPRISM(test_idrc_set)
#' exportDRCDataToPRISM(test_idrc_set, include_control=TRUE)
#' exportDRCDataToPRISM(test_idrc_set, include_control=TRUE, add_metadata=TRUE)
#'
exportDRCDataToPRISM <- function(idrc_set, include_control=FALSE, add_metadata=FALSE) {

    out_df <- exportDRCDataToDataFrame(idrc_set, include_control, add_metadata=FALSE) %>%
                    dplyr::transmute(sampleid, col_id=paste(round(conc,4), idx, sep='_'), value) %>%
                    tidyr::spread(col_id, value) %>%
                    as.data.frame()

    if(add_metadata & is.data.frame(idrc_set$metadata)) {
        out_df <- merge(out_df, idrc_set$metadata)
    }

    return(as.data.frame(out_df))

}
