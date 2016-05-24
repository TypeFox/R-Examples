#' Aggregation function
#' 
#' This function is intended to be used with \code{\link{sim_agg}} and not interactively. This is one implementation for aggregating data in a simulation set-up.
#' 
#' @param groupVars variable names as character identifying groups to be aggregated.
#' 
#' @details This function follows the split-apply-combine idiom. Each data set is split by the defined variables. Then the variables within each subset are aggregated (reduced to one row). Logical variables are reduced by \code{\link{any}}; for characters and factors dummy variables are created and the aggregate is the mean of each dummy; and for numerics the mean (removing NAs).
#' 
#' @seealso \code{\link{sim_agg}}
#' @export
#' 
#' @examples
#' sim_base() %>% sim_gen_x() %>% sim_gen_e() %>% sim_agg(agg_all())
agg_all <- function(groupVars = "idD") {
  force(groupVars)
  aggFun <- function(dat) {
    # character > factor > dummies (numeric)
    dat[sapply(dat, is.character)] <- dat[sapply(dat, is.character)] %>% lapply(as.factor)
    dummies <- names(dat)[sapply(dat, is.factor)] %>% 
      lapply(function(fac) {
        paste("~ -1", fac, sep = " + ") %>% as.formula %>% 
          model.matrix(data=dat) %>% as.data.frame
    }) %>% do.call(what=cbind)
    dat <- c(dat[!sapply(dat, is.factor)], dummies) %>% as.data.frame
    # Delete vars:
    dat <- dat[!(names(dat) %in% c("idU"))]
    # split:
    datList <- split(dat, as.list(dat[groupVars]))
    # apply:
    datList <- lapply(datList, 
           function(df) {
             dfOut <- df[1,]
             # logicals
             dfOut[sapply(df, is.logical)] <- lapply(df[sapply(df, is.logical)], any)
             # numeric
             dfOut[sapply(df, function(var) is.numeric(var) | is.integer(var))] <- 
               lapply(df[sapply(df, is.numeric)], mean, na.rm = TRUE)
             
             dfOut
           })
    # combine:
    out <- rbind_all(datList) %>% as.data.frame
    out
  }
  preserve_attributes(aggFun)
} 
