##' Retrieving or updating tip and node data in phylo4d objects
##'
##' Methods to retrieve or update tip, node or all data associated with a
##' phylogenetic tree stored as a phylo4d object
##'
##' @param x A \code{phylo4d} object
##' @param type The type of data to retrieve or update: \dQuote{\code{all}}
##' (default) for data associated with both tip and internal nodes,
##' \dQuote{\code{tip}} for data associated with tips only,
##' \dQuote{\code{internal}} for data associated with internal nodes only.
##' @param label.type How should the tip/node labels from the tree be returned?
##' \dQuote{\code{row.names}} returns them as row names of the data frame,
##' \dQuote{\code{column}} returns them in the first column of the data frame.
##' This options is useful in the case of missing (\code{NA}) or non-unique
##' labels.
##' @param empty.columns Should columns filled with \code{NA} be returned?
##' @param merge.data if tip or internal node data are provided and data already
##' exists for the other type, this determines whether columns with common names
##' will be merged together (default TRUE). If FALSE, columns with common names
##' will be preserved separately, with \dQuote{.tip} and \dQuote{.node} appended
##' to the names. This argument has no effect if tip and node data have no
##' column names in common, or if type=\dQuote{all}.
##' @param clear.all If only tip or internal node data are to be replaced,
##' should data of the other type be dropped?
##' @param \dots For the \code{tipData} and \code{nodeData} accessors,
##' further arguments to be used by \code{tdata}. For the replacement
##' forms, further arguments to be used to control matching between
##' tree and data (see Details section of \code{\link{phylo4d-methods}}).
##' @param value a data frame (or object to be coerced to one) to replace the
##' values associated with the nodes specified by the argument \code{type}
##' @return \code{tdata} returns a data frame
##' @section Methods: \describe{
##' \item{tdata}{\code{signature(object="phylo4d")}: retrieve or update data
##' associated with a tree in a \code{phylo4d} object} }
##' @author Ben Bolker, Thibaut Jombart, Francois Michonneau
##' @seealso \code{\link{phylo4d-methods}}, \code{\linkS4class{phylo4d}}
##' @export
##' @keywords methods
##' @include phylo4d-methods.R
##' @rdname tdata-methods
##' @examples
##'    data(geospiza)
##'    tdata(geospiza)
##'    tipData(geospiza) <- 1:nTips(geospiza)
##'    tdata(geospiza)
setGeneric("tdata", function(x, ...) {
    standardGeneric("tdata")
})

##' @rdname tdata-methods
##' @aliases tdata,phylo4d-method
setMethod("tdata", signature(x="phylo4d"),
  function(x, type=c("all", "tip", "internal"),
           label.type=c("row.names","column"),
           empty.columns=TRUE) {

      ## Returns data associated with the tree
      ## Note: the function checks for unique labels. It's currently unecessary
      ## but could be useful in the future if non-unique labels are allowed.

      type <- match.arg(type)
      label.type <- match.arg(label.type)

      ids <- nodeId(x, type)
      labs <- labels(x, type)
      ## replace any missing labels with node numbers
      labs[is.na(labs)] <- names(labs)[is.na(labs)]
      tdata <- x@data[match(ids, row.names(x@data)), , drop=FALSE]
      row.names(tdata) <- ids
      data.names <- labs[match(names(labs), rownames(tdata))]

      if (label.type == "row.names") {
          if (!any(duplicated(data.names)) &&
              ## length(data.names) > 0 &&
              !any(is.na(data.names)) ) {
              row.names(tdata) <- data.names
          }
          else {
              warning("Non-unique or missing labels found, ",
                      "labels cannot be coerced to tdata row.names. ",
                      "Use the label.type argument to include labels ",
                      "as first column of data.")
          }
      }
      if (identical(label.type,"column")) {
          tdata <- data.frame(label=data.names, tdata)
      }

      ## remove empty columns (filled with NAs)
      if(!empty.columns) {
          emptyCol <- apply(tdata, 2, function(x) all(is.na(x)))
          tdata <- tdata[, !emptyCol, drop=FALSE]
      }

      tdata
  })

##' @rdname tdata-methods
##' @aliases tdata<-
##' @export
setGeneric("tdata<-", function(x, ..., value) {
    standardGeneric("tdata<-")
})

##' @name tdata<-
##' @rdname tdata-methods
##' @aliases tdata<-,phylo4d-method tdata<-,phylo4d,ANY-method
setReplaceMethod("tdata", signature(x="phylo4d", value="ANY"),
    function(x, type = c("all", "tip", "internal"), merge.data = TRUE,
        clear.all = FALSE, ..., value) {

    type <- match.arg(type)

    ## format new data
    value <- formatData(x, value, type, keep.all=TRUE, ...)

    ## get old data to keep (if any)
    if (clear.all || type=="all") {
        keep <- NULL
    } else {
        if (type=="tip") {
            keep <- tdata(x, type="internal", empty.column=FALSE)
            keep <- formatData(x, keep, "internal", match.data=FALSE)
        } else if (type=="internal") {
            keep <- tdata(x, type="tip", empty.column=FALSE)
            keep <- formatData(x, keep, "tip", match.data=FALSE)
        }
    }

    ## create updated data
    updated.data <- switch(type,
        tip = .phylo4Data(x, tip.data=value, node.data=keep,
            merge.data=merge.data),
        internal = .phylo4Data(x, tip.data=keep, node.data=value,
            merge.data=merge.data),
        all = .phylo4Data(x, all.data=value, merge.data=merge.data))

    ## try to arrange new columns after old columns
    kept <- names(updated.data) %in% names(keep)
    old.cols <- names(updated.data)[kept]
    new.cols <- names(updated.data)[!kept]
    x@data <- updated.data[c(old.cols, new.cols)]

    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    return(x)
})

### Tip data wrappers
##' @rdname tdata-methods
##' @aliases tipData tipData-method
##' @export
setGeneric("tipData", function(x, ...) {
    standardGeneric("tipData")
})

##' @name tipData
##' @rdname tdata-methods
##' @aliases tipData,phylo4d-method
setMethod("tipData", signature(x="phylo4d"), function(x, ...) {
    tdata(x, type="tip", ...)
})

## tipData<-
##' @rdname tdata-methods
##' @aliases tipData<-
##' @export
setGeneric("tipData<-", function(x, ..., value) {
    standardGeneric("tipData<-")
})

##' @name tipData<-
##' @rdname tdata-methods
##' @aliases tipData<-,phylo4d-method tipData<-,phylo4d,ANY-method
setReplaceMethod("tipData", signature(x="phylo4d", value="ANY"),
    function(x, ...,  value) {
    tdata(x, type="tip", ...) <- value
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    return(x)
})

### Node data wrappers
##' @rdname tdata-methods
##' @aliases  nodeData nodeData-method
##' @export
setGeneric("nodeData", function(x, ...) {
    standardGeneric("nodeData")
})

##' @name nodeData
##' @rdname tdata-methods
##' @aliases nodeData,phylo4d-method
setMethod("nodeData", signature(x="phylo4d"), function(x, ...) {
    tdata(x, type="internal", ...)
})

## nodeData<-
##' @rdname tdata-methods
##' @aliases nodeData<-
##' @export
setGeneric("nodeData<-", function(x, ..., value) {
    standardGeneric("nodeData<-")
})

##' @name nodeData<-
##' @rdname tdata-methods
##' @aliases  nodeData<-,phylo4d-method nodeData<-,phylo4d,ANY-method
setReplaceMethod("nodeData", signature(x="phylo4d", value="ANY"),
    function(x, ...,  value) {
    tdata(x, type="internal", ...) <- value
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    return(x)
})
