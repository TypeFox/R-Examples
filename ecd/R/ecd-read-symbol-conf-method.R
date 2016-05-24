#' Read conf for sample data 
#' 
#' Read conf for sample data 
#'
#' @param symbol    Character. The symbol of sample data. Default: dji.
#' @param conf_file File name fof symbol config, default to conf/ecd-fit-conf.yml
#'
#' @return the conf object
#'
#' @keywords fit sample-data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' conf <- ecd.read_symbol_conf("dji")
#' }
### <======================================================================>
"ecd.read_symbol_conf" <- function(symbol, conf_file = "conf/ecd-fit-conf.yml") 
{
    # read the conf file
    conf_all <- yaml.load_file(conf_file)
    conf1 <- Filter(function(x) {x$symbol == symbol}, conf_all)
    if (length(conf1) != 1) {
        stop(paste("symbol", symbol, "doesn't exist in conf:", conf_file))
    }
    conf <- conf1[[1]]

    if (!("ecd" %in% names(conf))) {
        stop(paste("symbol", symbol, "doesn't have ecd configuration"))
    }
    if (!("breaks" %in% names(conf))) {
        stop(paste("symbol", symbol, "doesn't have breaks configuration"))
    }
    conf
}
### <---------------------------------------------------------------------->
