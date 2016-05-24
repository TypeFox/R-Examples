globalVariables(".")

#' An implementation of a stack
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' 
#' @section Fields:
#' \describe{
#'   \item{\code{stack}}{A list which can directly accessed or accessed through the functions below. The bottom (original) item, set to NULL, should not be modified.}}
#' 
#' @section Methods:
#' \describe{
#'   \item{\code{push(item, name = "")}}{Append an \code{item} to the \code{stack} with an optional \code{name} and return the item}}
#' 
#' @section Active Bindings:
#' \describe{
#'   \item{\code{height}}{This will return the length of the \code{stack}}}
#' 
#' \describe{
#'   \item{\code{pop}}{This will remove the last (most recent) item from the \code{stack} and return it}}
#' 
#' \describe{
#'   \item{\code{pop}}{This will return, but not remove, the last (most recent) item from the \code{stack}}}

stackClass = R6Class(
  "stackClass",
  public = list(
    stack = list("bottom" = NA),
    push = function(item, name = "") {
      self$stack[[self$height + 1]] = item
      names(self$stack)[self$height] = name
      self$peek}),
  active = list(
    height = function() length(self$stack),
    pop = function() {
      result = self$stack[[self$height]]
      self$stack[[self$height]] = NULL
      result},
    peek = function()
      self$stack[[self$height]]))


#' An implementation of an loop
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' 
#' @format An \code{\link{R6Class}} generator object
#' @keywords data
#' 
#' @section Inherits:
#' \describe{\item{\code{\link{stackClass}}}{}}
#' 
#' @section Methods:
#' \describe{\item{\code{
#'   begin(item, name = "")}}{
#'     Alias for \code{\link{stackClass}$push}}}
#' 
#' \describe{\item{\code{
#'   end(endData, FUN, ...)}}{
#'     Will return \code{FUN(\link{stackClass}$pop, endData, ...)}}}
#' 
#' \describe{\item{\code{
#'   cross(crossData, FUN, ...)}}{
#'     Will return \code{FUN(crossData, \link{stackClass}$pop, ...)}}}
loopClass = R6Class(
  inherit = stackClass,
  public = list(
    begin = function(item, name = "")
      self$push(item, name),
    end = function(endData, FUN, ...)
      FUN(self$pop, endData, ...),
    cross = function(crossData, FUN, ...)
      FUN(crossData, self$pop, ...)))

#' Amend variables with new information
#' 
#' Replace all non-NA values in one set of columns with values from another matching set
#' @importFrom magrittr %>%
#' @export
#' 
#' @param data A data frame
#' @param originalNames A vector of column names with out-of-date information
#' @param amendNames A vector of column names with amended information. They will be removed at the end of processing.
#' @return An amended \code{\link{tbl_df}}

amendColumns = function(data, originalNames, amendNames) {
  dataNames = dplyr::data_frame(originalNames,
                                amendNames) %>%
    dplyr::mutate(index = 1:length(originalNames))
  #build calls using ifelse
  calls = plyr::dlply(.data = dataNames,
                      .fun = function(row) lazyeval::lazy(
                        ifelse(is.na(amendNames), 
                               originalNames, 
                               amendNames)) %>%
                        lazyeval::interp(amendNames = as.name(row$amendNames),
                                         originalNames = as.name(row$originalNames)),
                      .variables = "index") %>%
    setNames(originalNames)
  
  data %>% 
    dplyr::mutate_(.dots = calls) %>%
    dplyr::select_(.dots = sprintf("-`%s`", amendNames))}

#' Fill variables with new information
#' 
#' Replace all NA values in one set of columns with values from another matching set
#' @importFrom magrittr %>%
#' @export
#' 
#' @param data A data frame
#' @param originalNames A vector of column names with out-of-date information
#' @param fillNames A vector of column names with new information. They will be removed at the end of processing.
#' @return A \code{\link{tbl_df}}
fillColumns = function(data, originalNames, fillNames) {
  dataNames = dplyr::data_frame(originalNames,
                                fillNames) %>%
    dplyr::mutate(index = 1:length(originalNames))
  
  #build calls using ifelse
  calls = 
    plyr::dlply(.data = dataNames,
                .fun = function(row) 
                  lazyeval::lazy(
                    ifelse(is.na(originalNames), 
                           fillNames, 
                           originalNames)) %>%
                  lazyeval::interp(fillNames = as.name(row$fillNames),
                                   originalNames = as.name(row$originalNames)),
                .variables = "index") %>%
    setNames(originalNames)
  
  data %>% 
    dplyr::mutate_(.dots = calls) %>%
    dplyr::select_(.dots = paste0("-", fillNames))}

#' Amend a dataframe with new information
#' 
#' \code{\link{full_join}} two dataframes. If there are matching columns, 
#' amend each \code{data} column with the corresponding \code{amendData} column using \code{\link{amendColumns}}.
#' 
#' @importFrom magrittr %>%
#' @export
#' 
#' @param data A data frame
#' @param amendData A data frame
#' @param by A quoted vector of column names to join by. If set to NULL or unspecified, will default to the grouping columns in data
#' @param suffix A suffix used internally. No existing column names should use this suffix.
#' @return An amended \code{\link{tbl_df}}
amend = function(data, amendData, by = NULL, suffix = "toFix") {
  
  #default by variables to groups
  if (is.null(by)) by = 
    data %>% 
    dplyr::groups() %>% 
    lapply(deparse) %>% 
    unlist %>% 
    as.vector 
  
  if (is.null(by)) stop("Defaulted to merging by data grouping variables. However, no grouping variables found")
  
  #figure out which columns need to be merged.
  commonNames = intersect(names(data), names(amendData)) %>% 
    setdiff(by)
  if (length(commonNames) != 0) message("Amending columns: ", paste(commonNames, collapse = ", "))
  
  #if no columns need to be merge, a simple full join
  if (length(commonNames) == 0) dplyr::full_join(data, amendData) else {
    
    #else update columns then join
    toFix = paste0(commonNames, suffix = suffix)
    
    if (sum(toFix %in% names(amendData)) > 0) stop ("suffix conflict. Please choose another suffix.")
    
    names(toFix) = commonNames
    
    byLiteral = by %>% sprintf("`%s`", .)
    
    amendData %>%
      plyr::rename(toFix) %>%
      dplyr::full_join(data, by) %>% 
      amendColumns(commonNames, unname(toFix)) %>%
      dplyr::arrange_(.dots = byLiteral)}}

#' Insert new information into a dataframe.
#' 
#' \code{\link{anti_join}} data with insertData, then \code{\link{bind_cols}} of insertData, then arrange by \code{by} variables.
#' @importFrom magrittr %>%
#' @export
#' 
#' @param data A data frame
#' @param insertData A data frame
#' @param by A quoted vector of column names to join by.
#' @return An inserted \code{\link{tbl_df}}
insert = function(data, insertData, by)
  data %>%
  dplyr::anti_join(insertData, by = by) %>%
  dplyr::bind_rows(insertData) %>%
  dplyr::arrange_(.dots = by)