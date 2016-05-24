#' Sample data fit
#' 
#' Fitting sample data to ecd with a starting set of parameters.
#' This is the highest level wrapper of the fitting routine.
#'
#' @param symbol    Character. The symbol of sample data. Default: dji.
#' @param iter      A length-one numeric. Number of maximum iterations. Default: 1000.
#' @param FIT       Logical, indicating whether to call linear regression, default = FALSE
#' @param EPS       Logical, indicating whether to save the plot to EPS, default = FALSE
#' @param conf_file File name fof symbol config, default to conf/ecd-fit-conf.yml
#' @param eps_file  File name for eps output
#' @param qa.fit    Logical, qa the standardfit_fn once.
#'
#' @return Final ecd object
#'
#' @keywords fit sample-data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dji <- ecd.fit_data("dji", FIT=T)
#' }
### <======================================================================>
"ecd.fit_data" <- function(symbol = "dji", iter=1000,
                           FIT = FALSE, 
                           EPS = FALSE,
                           conf_file = "conf/ecd-fit-conf.yml", 
                           eps_file = NULL, 
                           qa.fit = FALSE) 
{
    ts <- ecd.data(symbol)
    conf <- ecd.read_symbol_conf(symbol, conf_file)    
    
    ecd.fit_ts_conf(ts, conf,
                    iter = iter,
                    FIT = FIT, 
                    EPS = EPS,
                    eps_file = eps_file, 
                    qa.fit = qa.fit)

}
### <---------------------------------------------------------------------->
