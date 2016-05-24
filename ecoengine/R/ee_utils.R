#' @noRd
ee_base_url <- function() {
    "https://ecoengine.berkeley.edu/api/" 
    # "https://54.215.216.86/" # This is the current URL
}

#' ee_cbind
#'
#' Allows for combining split ecoengine calls (e.g. paginated calls) back into one single result object
#' @param results A list of objects of class \code{ecoengine}
#' @export
#' @importFrom assertthat assert_that
#' @examples \dontrun{
#' x1 <- ee_observations(genus = "Lynx", page = 1)
#' x2 <- ee_observations(genus = "Lynx", page = 2)
#' x12 <- ee_cbind(list(x1, x2))
#'}
ee_cbind <- function(results) {
    assert_that(class(results) == "list")
    # Fix to remove ldply
    res2 <- ldply(results, function(x) {
        y <- x[-4]
        data.frame(LinearizeNestedList(y))
    })

    page_loc <- which(names(res2) == "call.page")
    if(nrow(res2[!duplicated(res2[, -page_loc]), ]) == 1) {
        data <- ldply(results, function(x) x[[4]])
        res <- results[[1]]
        res$data <- data
        # Next two lines add pages back. Combines all pages, sorts, then adds.
        pages <- res2$call.page
        res$call$page <- paste0(pages[order(pages)], collapse = " ")
        res
    } else {
        stop("Cannot combine results from unidentical calls")
    }

}



#' Print a summary for an ecoengine object
#' @method print ecoengine
#' @export  
#' @param x An object of class \code{ecoengine}
#'   
#' @param ... additional arguments
print.ecoengine <- function(x, ...) {
cat(sprintf("[Total results on the server]: %s \n", x$results))
cat("[Args]: \n")
suppressWarnings(pretty_lists(x$call))
cat(sprintf("[Type]: %s \n", x$type))
cat(sprintf("[Number of results retrieved]: %s \n", nrow(x$data)))
}


#' ee_pages - Returns total number of pages for any ecoengine request
#'
#' @param ee Object of class \code{ecoengine}
#' @param  page_size Default page size. Currently set to \code{1000} package wide.
#' @export
#' @return integer
ee_pages <- function(ee, page_size = 1000) {
    if(!identical(class(ee), "ecoengine"))
        stop("Object must be of class ecoengine")

   ceiling(ee$results/page_size)
}


#'ecoengine paginator
#'
#'Takes a page range and total number of observations to return the right sequence of pages that need to be crawled.
#' @param page requested page number or page range. Can also be "all"
#' @param  total_obs Total number of records available for any search query
#' @param  page_size Default is \code{25}. Set higher if needed.
#' @export
#' @examples \dontrun{
#' ee_paginator(1, 100)
#' ee_paginator("all", 100)
#' ee_paginator(1:2, 1000)
#' ee_paginator(1:4, 3800)
#' # This will return an error since there are only 4 pages per 100 observations
#' ee_paginator(1:5, 10000)
#' }
ee_paginator <- function(page, total_obs, page_size = 1000) {
        all_pages <- ceiling(total_obs/page_size)
        if(total_obs < page_size) { req_pages <- 1 }
        if(identical(page, "all")) { req_pages <- seq_along(1: all_pages)}
        if(length(page) == 1 & identical(class(page), "numeric")) { req_pages <- page }
        if(identical(class(page), "integer")) {
            if(max(page) > all_pages) {
                stop("Pages requested outside the range")
            } else {
                req_pages <- seq_along(page)
            }
        }

        req_pages
    }


#' @noRd
ee_compact <- function(l) Filter(Negate(is.null), l)

#' @noRd
# Internal function to convert list to data.frame when it contains NULL
rbindfillnull <- function(x) {
  x <- unlist(x)
  x[is.null(x)] <- "none"
  data.frame(as.list(x), stringsAsFactors = FALSE)
}





#' @noRd
# This is an internal function to linearize lists
# Source: https://gist.github.com/mrdwab/4205477
# Author page (currently unreachable):  https://sites.google.com/site/akhilsbehl/geekspace/articles/r/linearize_nested_lists_in
# Original Author: Akhil S Bhel
# Notes: Current author could not be reached and original site () appears defunct. Copyright remains with original author
LinearizeNestedList <- function(NList, LinearizeDataFrames=FALSE,
                                NameSep="/", ForceNames=FALSE) {
    # LinearizeNestedList:
    #
    # https://sites.google.com/site/akhilsbehl/geekspace/
    #         articles/r/linearize_nested_lists_in_r
    #
    # Akhil S Bhel
    # 
    # Implements a recursive algorithm to linearize nested lists upto any
    # arbitrary level of nesting (limited by R's allowance for recursion-depth).
    # By linearization, it is meant to bring all list branches emanating from
    # any nth-nested trunk upto the top-level trunk s.t. the return value is a
    # simple non-nested list having all branches emanating from this top-level
    # branch.
    #
    # Since dataframes are essentially lists a boolean option is provided to
    # switch on/off the linearization of dataframes. This has been found
    # desirable in the author's experience.
    #
    # Also, one'd typically want to preserve names in the lists in a way as to
    # clearly denote the association of any list element to it's nth-level
    # history. As such we provide a clean and simple method of preserving names
    # information of list elements. The names at any level of nesting are
    # appended to the names of all preceding trunks using the `NameSep` option
    # string as the seperator. The default `/` has been chosen to mimic the unix
    # tradition of filesystem hierarchies. The default behavior works with
    # existing names at any n-th level trunk, if found; otherwise, coerces simple
    # numeric names corresponding to the position of a list element on the
    # nth-trunk. Note, however, that this naming pattern does not ensure unique
    # names for all elements in the resulting list. If the nested lists had
    # non-unique names in a trunk the same would be reflected in the final list.
    # Also, note that the function does not at all handle cases where `some`
    # names are missing and some are not.
    #
    # Clearly, preserving the n-level hierarchy of branches in the element names
    # may lead to names that are too long. Often, only the depth of a list
    # element may only be important. To deal with this possibility a boolean
    # option called `ForceNames` has been provided. ForceNames shall drop all
    # original names in the lists and coerce simple numeric names which simply
    # indicate the position of an element at the nth-level trunk as well as all
    # preceding trunk numbers.
    #
    # Returns:
    # LinearList: Named list.
    #
    # Sanity checks:
    #
    stopifnot(is.character(NameSep), length(NameSep) == 1)
    stopifnot(is.logical(LinearizeDataFrames), length(LinearizeDataFrames) == 1)
    stopifnot(is.logical(ForceNames), length(ForceNames) == 1)
    if (! is.list(NList)) return(NList)
    #
    # If no names on the top-level list coerce names. Recursion shall handle
    # naming at all levels.
    #
    if (is.null(names(NList)) | ForceNames == TRUE)
        names(NList) <- as.character(1:length(NList))
    #
    # If simply a dataframe deal promptly.
    #
    if (is.data.frame(NList) & LinearizeDataFrames == FALSE)
        return(NList)
    if (is.data.frame(NList) & LinearizeDataFrames == TRUE)
        return(as.list(NList))
    #
    # Book-keeping code to employ a while loop.
    #
    A <- 1
    B <- length(NList)
    #
    # We use a while loop to deal with the fact that the length of the nested
    # list grows dynamically in the process of linearization.
    #
    while (A <= B) {
        Element <- NList[[A]]
        EName <- names(NList)[A]
        if (is.list(Element)) {
            #
            # Before and After to keep track of the status of the top-level trunk
            # below and above the current element.
            #
            if (A == 1) {
                Before <- NULL
            } else {
                Before <- NList[1:(A - 1)]
            }
            if (A == B) {
                After <- NULL
            } else {
                After <- NList[(A + 1):B]
            }
            #
            # Treat dataframes specially.
            #
            if (is.data.frame(Element)) {
                if (LinearizeDataFrames == TRUE) {
                    #
                    # `Jump` takes care of how much the list shall grow in this step.
                    #
                    Jump <- length(Element)
                    NList[[A]] <- NULL
                    #
                    # Generate or coerce names as need be.
                    #
                    if (is.null(names(Element)) | ForceNames == TRUE)
                        names(Element) <- as.character(1:length(Element))
                    #
                    # Just throw back as list since dataframes have no nesting.
                    #
                    Element <- as.list(Element)
                    #
                    # Update names
                    #
                    names(Element) <- paste(EName, names(Element), sep=NameSep)
                    #
                    # Plug the branch back into the top-level trunk.
                    #
                    NList <- c(Before, Element, After)
                }
                Jump <- 1
            } else {
                NList[[A]] <- NULL
                #
                # Go recursive! :)
                #
                if (is.null(names(Element)) | ForceNames == TRUE)
                    names(Element) <- as.character(1:length(Element))
                Element <- LinearizeNestedList(Element, LinearizeDataFrames,
                                               NameSep, ForceNames)
                names(Element) <- paste(EName, names(Element), sep=NameSep)
                Jump <- length(Element)
                NList <- c(Before, Element, After)
            }
        } else {
            Jump <- 1
        }
        #
        # Update book-keeping variables.
        #
        A <- A + Jump
        B <- length(NList)
    }
    return(NList)
}


#' @noRd
pretty_lists <- function(x)
{
   for(key in names(x)){
      value <- format(x[[key]])
      if(value == "") next
      cat(key, "=", value, "\n")
   }
   invisible(x)
}
