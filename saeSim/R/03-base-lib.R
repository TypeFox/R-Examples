#' Construct data with id-variables
#' 
#' @description This function constructs a data frame with grouping/\code{id} variables.
#' 
#' @param nDomains The number of domains.
#' @param nUnits The number of units in each domain. Can have \code{length(nUnits) > 1}.
#' @param nTime The number of time points for each units.
#'  
#' @return
#' Return a \code{data.frame} with variables \code{idD} as ID-variable for
#' domains, and \code{idU} as ID-variable for units.
#' 
#' @rdname base_id
#' @export
#' @examples
#' base_id(2, 2)
#' base_id(2, c(2, 3))
base_id <- function(nDomains = 10, nUnits = 10) {
  
  stopifnot(
    length(nDomains) == 1, 
    if (length(nUnits) > 1) length(nUnits) == nDomains else TRUE
  )
  
  out <- data.frame(idD = rep(1:nDomains, times = nUnits)) %>% 
    group_by_("idD") %>% mutate(idU = 1:n()) %>% arrange_("idD", "idU")
  
  out <- if (length(nUnits) == 1 && nUnits == 1) out["idD"] else out
  as.data.frame(out)
  
}

#' @export
#' @rdname base_id
base_id_temporal <- function(nDomains = 10, nUnits = 10, nTime = 10) {
  stopifnot(length(nTime) == 1)
  dat <- base_id(nDomains, nUnits)
  dat <- dat[rep(1:nrow(dat), nTime), , drop = FALSE]
  dat <- split(dat, dat, drop = TRUE) %>% 
    lapply(function(df) {df$idT <- 1:nrow(df); df}) %>%
    do.call(what = rbind)
  rownames(dat) <- NULL
  dat[order(dat$idD), ]
}

#' Add id-variables to data
#' 
#' Use this function to add id-variables to your data.
#' 
#' @param data a data.frame.
#' @param domainId variable names in \code{data} as character which will identify the areas/domains/groups/cluster in the data.
#' 
#' @rdname base_add_id
#' @export
base_add_id <- function(data, domainId) {
  dataList <- split(data, data[domainId])
  nUnits <- sapply(dataList, nrow)
  nDomains <- length(nUnits)
  cbind(base_id(nDomains, nUnits), rbind_all(dataList))
}