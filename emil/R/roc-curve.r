#' Calculate ROC curves
#'
#' @param result Modeling results, as returned by \code{\link{evaluate}}.
#' @param y True response vector used to create \code{result}.
#' @param resample Resampling scheme used to create \code{result}.
#' @param class The class of interest to create ROC-curves for.
#' @param statistic The name of the statistic (as returned by the prediction
#'   function of the modeling procedure).
#' @return A data frame of class \dQuote{roc}.
#' @example examples/roc-curve.r
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @importFrom dplyr do_
#' @importFrom magrittr "%<>%"
#' @export
roc_curve <- function(result, y, resample, class=levels(y), statistic="probability"){
    stopifnot(inherits(result, "modeling_result"))
    stopifnot(is.factor(y))
    if(is.numeric(class)) class <- levels(y)[class]
    roc_fun <- function(stat){
        y_sorted <- y[stat$id][order(stat$statistic)] == as.character(stat$class[1])
        data.frame(
            sensitivity = c(rev(cumsum(rev(y_sorted)))/sum(y_sorted), 0),
            specificity = c(0, cumsum(!y_sorted)/sum(!y_sorted))
        )
    }
    if(is_multi_procedure(result)){
        stat_table <- select(result, fold=resample, method=TRUE, "prediction",
                             statistic, statistic = class)
        if(!is.factor(stat_table$method))
            stat_table$method <- factor(stat_table$method)
    } else {
        stat_table <- select(result, fold=resample, "prediction",
                             statistic, statistic = class)
    }
    if(length(class) > 1){
        stat_table %<>% gather_("class", "statistic", class)
    } else {
        stat_table %<>% mutate_(class = "factor(class)")
    }
    stat_table %<>%
        group_by_(.dots=intersect(c("fold", "method", "class"), colnames(stat_table))) %>%
        do_(~roc_fun(.))
    structure(stat_table %>% as.data.frame,
              class=c("roc_curve", "data.frame"))
}

#' @param x Roc curve object, as returned by \code{roc_curve}.
#' @param ... Sent to \code{\link{as.data.frame}} or \code{\link{as.data.table}}.
#' @method as.data.table roc_curve
#' @rdname roc_curve
#' @importFrom data.table as.data.table
#' @export
as.data.table.roc_curve <- function(x, ...){
    if(inherits(x, "data.table")) return(x)
    class(x) <- class(x)[-1]
    structure(as.data.table(x, ...), class=c("roc_curve", "data.table"))
}

#' @method as.data.frame roc_curve
#' @rdname roc_curve
#' @export
as.data.frame.roc_curve <- function(x, ...){
    if(inherits(x, "data.frame")) return(x)
    class(x) <- class(x)[-1]
    structure(as.data.frame(x, ...), class=c("roc_curve", "data.frame"))
}

#' @importFrom ggplot2 ggplot aes_string geom_segment geom_path facet_wrap
#' @method plot roc_curve
#' @rdname roc_curve
#' @export
plot.roc_curve <- function(x, ...){
    p <- ggplot(x, aes_string(x="1-specificity", y="sensitivity",
                              group="paste(class, fold, method)")) + 
        geom_segment(x=0, y=0, xend=1, yend=1, linetype="dashed", size=.3) +
        geom_path(aes_string(colour="class"))
    if("method" %in% colnames(x)){
        p + facet_wrap(~method)
    } else {
        p
    }
}
