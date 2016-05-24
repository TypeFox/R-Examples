
#' rangeMapper browser
#'
#' @param path path to a rangeMapper project. If missing a demo project is created on the fly.
#' @export
#' @examples
#' \dontrun{
#' View_rmap()
#'}
#'
View_rmap <- function(path) {

    if(missing(path)) {

        br = rgdal::readOGR(system.file(package = "rangeMapper","extdata", "wrens", "vector_combined"), "wrens", verbose = FALSE)
        d = read.csv2(system.file(package = "rangeMapper","data", "wrens.csv")) %>%
             subset(., select = c('sci_name', 'body_size', 'body_mass', 'clutch_size') )
        con = ramp("wrens.sqlite", gridSize = 2, spdf = br, biotab = d, ID = "sci_name",
                    metadata = rangeTraits()[1], FUN = "median", overwrite = TRUE)
        path = paste(tempdir(), "wrens.sqlite", sep = .Platform$file.sep)
        RSQLite::dbDisconnect(con)
        }

     options(rangeMapper.path = path)

    if(interactive())
	   shiny::runApp(system.file('GUI', package = 'rangeMapper'))
	}


