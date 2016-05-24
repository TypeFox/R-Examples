#####################################################
## report_aggregator.r
##
## tools for aggregating respondent reports, taking
## account of weights
##

#####################################################
##' aggregate a reported quantity by groups
##' 
##' This function takes a quantity and aggregates it by groups,
##' using the design weights.
##'
##' @param resp.data the data 
##' @param attribute.names the names of the variables that define
##'          the groups for which the qoi should be aggregated
##' @param qoi the variable with quantity to aggregate
##' @param weights analysis weights
##' @param qoi.name the name of the qoi
##' @param dropmiss NOT YET IMPLEMENTED
##' @return the estimated average degree for respondents in each
##'         of the categories given by \code{attribute.names}
##' @rdname report.aggregator
report.aggregator_ <- function(resp.data,
                               attribute.names,
                               qoi,
                               weights,
                               qoi.name,
                               dropmiss=FALSE) {
  
  resp.data <- resp.data

  wdat <- select_(resp.data, .dots=weights)
  qdat <- select_(resp.data, .dots=qoi)
  adat <- select_(resp.data, .dots=attribute.names)

  df <- bind_cols(wdat, qdat, adat)

  # see
  # http://stackoverflow.com/questions/21208801/group-by-multiple-columns-in-dplyr-using-string-vector-input

  wgt.col <- as.symbol(names(df)[1])
  qoi.col <- as.symbol(names(df)[2])

  grouping.cols <- names(df)[-1:-2]
  dots <- lapply(grouping.cols, as.symbol)

  df.summ <- df %>% group_by_(.dots=dots) %>%
             dplyr::summarise_(mean.qoi = interp(~weighted.mean(a, w=b), a=qoi.col, b=wgt.col),
                        sum.qoi = interp(~sum(a*b), a=qoi.col, b=wgt.col),
                        wgt.total = interp(~sum(b), b=wgt.col),
                        wgt.inv.total = interp(~sum(1/b), b=wgt.col),
                        ## this is a hack because n() generates errors due to
                        ## plyr/dplyr import conflict ( and it is hard to regulate
                        ## import order with package infrastructure)
                        num.obs = interp(~length(a), a=qoi.col))

  toren <- list(~mean.qoi, ~sum.qoi, ~wgt.total, ~wgt.inv.total, ~num.obs)
  newnames <- paste0(c("mean.", "sum.", "wgt.total.", 
                       "wgt.inv.total.", "num.obs."), qoi.name)

  df.summ <- dplyr::rename_(df.summ,
                     .dots=setNames(toren, newnames))

  return(df.summ)

}

#####################################################
##' @rdname report.aggregator
report.aggregator <- function(resp.data,
                              attribute.names=NULL, 
                              qoi,
                              weights,
                              qoi.name=NULL,
                              dropmiss=FALSE) {

    report.aggregator_(resp.data,
                       attribute.names=lazy(attribute.names, env=resp.data),
                       qoi=lazy(qoi, resp.data),
                       weights=lazy(weights, resp.data),
                       ifelse(is.null(qoi.name), "qoi", lazy(qoi.name)),
                       dropmiss)

}

