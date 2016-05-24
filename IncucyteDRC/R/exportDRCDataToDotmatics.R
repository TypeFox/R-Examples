#' exportDRCDataToDotmatics
#'
#' Exports data in Dotmatics format for an IncucyteDRCSet or IncucyteDRCSetList object that represents
#' a full plate of data.
#' The platemap must be provided as an IncucyteDRCPlatemap object and should correspond to the data
#' contained in the IncucyteDRCSet or IncucyteDRCSetList object.
#'
#' @param idrc_set_list IncucyteDRCSetList object
#' @param platemap IncucyteDRCPlatemap object from importPlatemap or importPlatemapXML
#'
#' @return list object with samplelist and data elements
#' @export
#'
#' @examples
#' pm_file <- system.file(file='extdata/example.PlateMap', package='IncucyteDRC')
#' test_pm <- importPlatemapXML(pm_file)
#' data_file <- system.file(file='extdata/example_data.txt', package='IncucyteDRC')
#' test_data <- importIncucyteData(data_file, metric='pc')
#'
#' test_drc <- splitIncucyteDRCPlateData(test_pm, test_data, group_columns='growthcondition')
#' test_drc <- lapply(test_drc, fitGrowthCurvesIndividual)
#' test_drc <- lapply(test_drc, fitGrowthCurvesGrouped)
#' test_drc <- lapply(test_drc, calculateDRCData, cut_time=100)
#' print(test_drc)
#' print(test_drc[[2]])
#' exportDRCDataToDotmatics(test_drc, test_pm)
#' exportDRCDataToDotmatics(test_drc[[2]], test_pm)
exportDRCDataToDotmatics <- function(idrc_set_list, platemap) {

    if(inherits(idrc_set_list, 'IncucyteDRCSet')) {
        message('Generating Dotmatics data for IncucyteDRCSet')
        idrc_set_list <- list(idrc_set_list)
    } else if(inherits(idrc_set_list, 'list') & all(sapply(test_drc, class) == 'IncucyteDRCSet')) {
        message('Generating Dotmatics data for IncucyteDRCSetList')
    } else {
        stop("Input should be an IncucyteDRCSet or IncucyteDRCSetList object")
    }

    if(is.null(attr(platemap, 'IncucyteDRCPlatemap'))) {
        warning('Recommended that platemap data frames are parsed through importPlatemap function to check formatting')
    }

    #combine the drc data from the list object
    data_df <- lapply(idrc_set_list, function(x) {
        x$drc_data %>% dplyr::transmute(platename=x$platedata$plateid, wellid, value=cut_val)
    }) %>% dplyr::bind_rows()

    #make rows for the blanks
    blanks_df <- platemap %>%
        dplyr::filter(samptype == 'B') %>%
        dplyr::transmute(row, col, wellid, value=0)

    #combine data with platemap
    combined_df <- data_df %>%
        dplyr::inner_join(platemap, by='wellid') %>%
        dplyr::arrange(row,col)

    #combine data and blanks, order, add platename
    out_df <- combined_df %>%
        dplyr::select(row, col, wellid, value) %>%
        dplyr::bind_rows(blanks_df) %>%
        dplyr::arrange(row,col) %>%
        dplyr::transmute(platename=unique(combined_df$platename), wellid, value) %>%
        as.data.frame()

    #check that the platemap and out_df are the same length
    if(nrow(out_df) != nrow(platemap)) {
        warning('Data is missing for some wells in the supplied platemap')
    }

    #try to derive a sensible sample list - this may have to be modified by the user
    samples_df <- combined_df %>%
        dplyr::filter(samptype=='S') %>%
        dplyr::group_by(platename, sampleid, celltype, growthcondition) %>%
        dplyr::summarise(max_conc = max(conc),
                  first_row = first(row),
                  first_col = first(col),
                  N = n()) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(first_row, first_col) %>%
        as.data.frame()

    #return a list object
    return(list(data=out_df,
                samplelist=samples_df,
                platemap=platemap))

}
