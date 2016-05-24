#' Timeseries fitting utility
#' 
#' Fitting timeseries with provided conf as starting set of parameters.
#'
#' @param ts        An xts object from either \code{ecd.data} or \code{ecd.df2ts}.
#' @param conf      A nested list object, the configuration.
#' @param iter      A length-one numeric. Number of maximum iterations. Default: 1000.
#' @param FIT       Logical, indicating whether to call linear regression, default = FALSE
#' @param EPS       Logical, indicating whether to save the plot to EPS, default = FALSE
#' @param eps_file  File name for eps output
#' @param qa.fit    Logical, qa the standardfit_fn once.
#'
#' @return Final ecd object
#'
#' @keywords fit timeseries
#'
#' @export
#'
#' @examples
#' \dontrun{
#' d <- ecd.fit_ts_conf(ts, conf)
#' }
### <======================================================================>
"ecd.fit_ts_conf" <- function(ts, conf, 
                           iter=1000,
                           FIT = FALSE, 
                           EPS = FALSE,
                           eps_file = NULL, 
                           qa.fit = FALSE) 
{
    symbol <- conf$symbol

    if (!("ecd" %in% names(conf))) {
        stop(paste("symbol", symbol, "doesn't have ecd configuration"))
    }
    if (!("breaks" %in% names(conf))) {
        stop(paste("symbol", symbol, "doesn't have breaks configuration"))
    }
            
    breaks <- conf$breaks
    
    # optional
    merge_tails <- if ("merge_tails" %in% names(conf)) {
        t <- conf$merge_tails
        if ("fit" %in% names(t)) t$fit
        else if ("plot" %in% names(t)) t$plot
        else t
    } else c(0, 0)

    weights <- list()
    if("weights" %in% names(conf)) { weights <- conf$weights }

    # ---------------------------------------
    # suggestion - first round, disable qq_df and pdf_df
    # ---------------------------------------
    
    print(paste(
        "Loaded symbol", symbol,";",
        "breaks=", breaks,
        "merge_tails= (", paste(merge_tails, collapse = ','),")",
        "weights= (", paste(weights, collapse = ','),")" # TODO this isn't good
        ))
    
	ts1 <- ecd.data_stats(ts, breaks=breaks, merge_tails=merge_tails)
    attr <- xtsAttributes(ts1)
    
    # construct ecd object from conf
    d <- as.numeric(conf$ecd)
    d[3] <- if(is.na(d[3])) attr$stdev else d[3] # sigma
    d[4] <- if(is.na(d[4])) 0 else d[4] # beta
    d[5] <- if(is.na(d[5])) attr$mean else d[5] # mu
    object <- ecd(d[1], d[2], d[3], d[4], d[5])

    print(paste("Constructed ecd:",
                "alpha=", object@alpha, "gamma=", object@gamma,
                "R=", sprintf("%.2f",object@R),
                "theta=", sprintf("%.2f",object@theta),
                "degree=", sprintf("%.1f",object@theta/pi*180)))

    # ready to go into fit API
    if ( qa.fit ) {
        fn.out <- ecd.standardfit_qa(object, ts1, 
                                     plotqq=1, weights=weights, debug=2)
        return (fn.out)
    }
    else if ( FIT ) {
        fit <- ecd.standardfit(object, ts1, 
                               trace=1, iter=iter, plotqq=1, 
                               weights=weights)
        object <- fit$dist
        #  opt.out <- fit$opt.out
    }

    # graphs
    # merge_tails may be different
    merge_tails <- if ("merge_tails" %in% names(conf)) {
        t <- conf$merge_tails
        if ("plot" %in% names(t)) t$plot else t
    } else c(0, 0)

    ts2 <- ecd.data_stats(ts, breaks=breaks,
                          merge_tails=merge_tails,
                          with.tail=TRUE)
    
    plot_2x2(object, ts2, EPS, eps_file = eps_file)
	invisible(object)
}
### <---------------------------------------------------------------------->
