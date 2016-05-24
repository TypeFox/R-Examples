#' Plot data representing Glycans in boxplot or violin mode using ggplot2
#'
#' This function constructs standard plots in exploratory analysis 
#' of N-Glycans. 
#'
#' @author Ivo Ugrina
#' @export
#' @importFrom grDevices boxplot.stats
#' @importFrom stats IQR kruskal.test p.adjust
#' @param data data frame which holds columns representing Glycans.
#'   These column names must start with 'GP'.
#' @param collapse should Glycans be presented in one facet (default)
#'   or with more facets (one per Glycan).
#' @param violin should Glycans be presented in a boxplot (default)
#'   or violin format.
#' @param group this a possible grouping parameter on which
#'   stratification of \code{data} should be conducted. It should be
#'   a name of one of the columns in dataframe \code{data}
#'   and of type \code{factor}.
#' @param all should all of the variables (default) be presented in the plot
#'   or only those that have significant p-values. This variable is
#'   meaningful only when \code{group} is not \code{NULL} since the testing
#'   of differences is conducted  between different groups represented
#'   by \code{group} variable. If \code{group} has only 2 levels then
#'   Mann-Whitney-Wilcoxon (\code{\link{wilcox.test}}) test is conducted.
#'   Otherwise, Kruskal-Wallis test is conducted (\code{\link{kruskal.test}}).
#'   Obtained p-values are adjusted to multiple testing with \code{\link{p.adjust}}.
#' @param p.adjust.method method used for adjustment of p-values to multiple
#'   testing. Variable p.adjust.method must be an element of
#'   \code{\link{p.adjust.methods}}.
#' @param print.p.values should p-values be printed on plots
#' @param log.transform should Glycans be log transform prior to plotting.
#' @param glyco.names names of columns that represent glycan data. If \code{NULL}
#'   all columns starting with 'GP' in their names will be used
#' @return Returns a list consisting of p-values, adjusted p-vales and the plot.
#' @examples
#' devAskNewPage(TRUE)
#' exampleData <- data.frame(ID=1:100, GP1=runif(100),
#'   GP2=rexp(100,0.2), GP3=rgamma(100, 3),
#'   Plate=factor(sample(1:2,100,replace=TRUE)))
#' glyco.plot(exampleData)
#' glyco.plot(exampleData, group='Plate', collapse=FALSE, log=TRUE)
glyco.plot <- function(data, collapse = TRUE, violin = FALSE, group = NULL, all = TRUE, 
    p.adjust.method = "holm", print.p.values = TRUE, log.transform = FALSE, glyco.names = NULL) {
    
    # basic tests
    
    # test: if data is data.frame if group variable is indeed a factor variable and a
    # part of the data if collapse, all, violin, log.transform are logical if
    # p.adjust.method is in stats::p.adjust.methods
    stopifnot(is.data.frame(data))
    if (!is.null(group)) {
        stopifnot(group %in% names(data))
        stopifnot(is.factor(data[[group]]))
    }
    stopifnot(p.adjust.method %in% stats::p.adjust.methods, is.logical(collapse), 
        is.logical(all), is.logical(violin), is.logical(log.transform))
    
    # other
    if (is.null(glyco.names)) {
        tmp <- grep("^GP", names(data))
        gps <- names(data)[tmp]
        not.gps <- names(data[-tmp])
    } else {
        gps <- glyco.names
        not.gps <- names(data)[!names(data) %in% gps]
    }
    
    stopifnot(length(gps) > 0)
    
    cols <- c(group, gps)
    
    newdata <- tidyr::gather_(data[,cols], "variable", "value", gps)
    if (TRUE == log.transform) {
        newdata$value <- log(newdata$value)
    }
    
    if (is.null(group)) {
        if (collapse) {
            p <- ggplot2::ggplot(newdata, ggplot2::aes(variable, value))
            
            if (violin) 
                p2 <- p + ggplot2::geom_violin() else p2 <- p + ggplot2::geom_boxplot()
            
            p2 + ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
        } else {
            p <- ggplot2::ggplot(newdata, ggplot2::aes(variable, value))
            
            if (violin) 
                p2 <- p + ggplot2::geom_violin() else p2 <- p + ggplot2::geom_boxplot()
            
            p2 + ggplot2::facet_wrap(~variable, ncol = 6, scales = "free") + ggplot2::xlab("") + 
                ggplot2::ylab("") + ggplot2::theme(axis.text.x = ggplot2::element_blank())
        }
        
    } else {
        if (2 == length(levels(data[[group]]))) {
            wilcox.test.multinumeric <- function(Y, x) {
                apply(Y, 2, function(column) {
                  tmp <- coin::wilcox_test(column ~ x)
                  coin::pvalue(tmp)
                })
            }
            
            X <- data
            p.val.unadj <- wilcox.test.multinumeric(X[, gps], X[[group]])
            p.val <- p.adjust(p.val.unadj, method = p.adjust.method)
        } else {
            kw.test.multinumeric <- function(Y, x) {
                apply(Y, 2, function(column) {
                  tmp <- kruskal.test(column, x)
                  tmp$p.value
                })
            }
            
            X <- data
            p.val.unadj <- kw.test.multinumeric(X[, gps], X[[group]])
            p.val <- p.adjust(p.val.unadj, method = p.adjust.method)
        }
        
        if (TRUE == all) {
            significant.gps <- gps
            s.pval <- p.val
        } else {
            s.ind <- which(p.val < 0.05)
            significant.gps <- names(s.ind)
            s.pval <- p.val[s.ind]
        }
        
        if (0 == length(s.pval)) {
            return(list(p.val.unadj = p.val.unadj, p.val = p.val, plot = NULL))
        }
        
        tmp <- data.frame(variable = significant.gps, p = round(s.pval, 3), stringsAsFactors = FALSE)
        tmp2 <- paste0(tmp$variable, ", p=", tmp$p)
        
        newdata <- merge(tmp, newdata, by = "variable")
        
        if (print.p.values) {
            newdata$desc <- factor(paste0(newdata$variable, ", p=", newdata$p), levels = tmp2)
        } else {
            newdata$desc <- factor(newdata$variable)
        }
        
        p <- ggplot2::ggplot(newdata, ggplot2::aes(desc, value))
        
        if (TRUE == violin) {
            p2 <- p + ggplot2::geom_violin(ggplot2::aes_string(fill = group))
        } else {
            p2 <- p + ggplot2::geom_boxplot(ggplot2::aes_string(fill = group))
        }
        
        if (collapse) {
            plot <- p2 + ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
        } else {
            plot <- p2 + ggplot2::facet_wrap(~desc, ncol = 6, scales = "free") + 
                ggplot2::xlab("") + ggplot2::ylab("") + ggplot2::theme(strip.text.x = ggplot2::element_blank())
        }
        
        return(list(p.val.unadj = p.val.unadj, p.val = p.val, plot = plot))
    }
}
