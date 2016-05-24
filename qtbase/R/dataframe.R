## Implement the QItemModel in QtCore using an R data.frame

##' The \code{qdataFrameModel} function creates a
##' \code{DataFrameModel}, an implementation of
##' \code{QAbstractItemModel} using a \code{data.frame}. This makes it
##' easy and fast to display and edit a \code{data.frame} in a
##' \code{QTableView} or any other derivative of
##' \code{QAbstractItemView}. The \code{qdataFrame} and
##' \code{qdataFrame<-} functions allow one to get and set the
##' \code{data.frame} underlying the model after construction.
##'
##' While a simple data.frame can be displayed as a textual table,
##' fancier tables require multiple data columns mapped to a single
##' model column, each playing a separate 'role'. To specify
##' additional roles, pass \code{useRoles = TRUE}. A role may be any
##' string; those used by Qt are listed in the \code{Qt::ItemDataRole}
##' enumeration. The \code{display} and \code{edit} roles are reserved
##' (see below). See the documentation of the
##' \code{QStyledItemDelegate} class for its expected data types for
##' each role.
##' 
##' A simple way to encode this is in the column name, syntax:
##' \code{[.headerName1][.headerName2][.etc].role}.
##' Examples:
##' \itemize{
##' \item{\code{.carColor.background} (background color for carColor column)}
##' \item{\code{.foreground} (foreground color for all columns)}
##' \item{\code{.firstName.lastName.font} (special font for first and last
##' name columns)}
##' }
##'
##' The set of model columns is derived from the unique header names.
##' Display-role columns are those not prefixed by a period. If the
##' column name in the data matches a string in the \code{editable}
##' argument, the data is used for both the edit and display roles.
##'
##' @note Calling the \code{headerData} method on
##' \code{DataFrameModel} from R will not yield the expected result,
##' because Smoke does not know of DataFrameModel and thus misses the
##' override of the non-pure virtual.  We can special case this if
##' need-be.
##' @title DataFrameModel
##' @param df The \code{data.frame} that provides the data of the model
##' @param parent The parent \code{QObject} for the model. Important
##' for preventing garbage collection of the model if the only
##' reference to it is through a view.
##' @param useRoles Whether to interpret column names as indicating
##' alternative roles; see details.
##' @param editable Character vector of column names in the
##' \code{data.frame} that should be editable
##' @param ... Extra arguments passed to \code{qdataFrame<-},
##' which actually loads the \code{data.frame} into the model.
##' @return An instance of C++ \code{DataFrameModel}
##' @author Michael Lawrence
##' @rdname DataFrameModel
qdataFrameModel <- function(df, parent = NULL, useRoles = FALSE,
                            editable = character(), ...)
{
  model <- .Call("qt_qdataFrameModel", parent, useRoles, editable,
                 PACKAGE="qtbase")
  qdataFrame(model, ...) <- df
  model
}

##' @param model \code{DataFrameModel} instance
##' @param value A \code{data.frame} that provides the data of the model
##' @rdname DataFrameModel
`qdataFrame<-` <- function(model, value)
{
  stopifnot(inherits(model, "DataFrameModel"))
  useRoles <- quseRoles(model)
  editable <- qeditable(model)
  df <- as.data.frame(value)
  ## this order must match the order of the Qt::ItemDataRole enumeration
  roleNames <- c("display", "decoration", "edit", "toolTip", "statusTip",
                 "whatsThis", "font", "textAlignment", "background",
                 "foreground", "checkState", "accessibleText",
                 "accessibleDescription", "sizeHint")
  createRoleList <- function()
    structure(vector("list", length(roleNames)), names = roleNames)
  roles <- createRoleList()
  if (useRoles) {
    getHasRole <- function(cn) structure(grepl("^\\.", cn), names = cn)
    getHeaderNames <- function(x) {
      headerNames <- sub("^\\.", "", sub("\\.[^.]*$", "", x))
      hasRole <- hasRole[x]
      headerNameList <- as.list(headerNames)
      headerNameList[hasRole] <- strsplit(headerNames[hasRole], "\\.")
      headerNameList
    }
    cn <- colnames(df)
    hasRole <- getHasRole(cn)
    header <- as.list(cn)
    header[hasRole] <- getHeaderNames(cn[hasRole])
    header <- unique(unlist(header))
    editOrDisplay <- ifelse(cn[!hasRole] %in% editable, "edit", "display")
    cn[!hasRole] <- paste(sub("(^[^\\.].*)", ".\\1.", cn[!hasRole]),
                          editOrDisplay, sep = "")
    names(hasRole) <- cn
    getRoleNames <- function(x) gsub(".*\\.", "", x)
    dataRoles <- getRoleNames(cn)
    resolveRole <- function(role) {
      headerNames <- getHeaderNames(role)
      nheaders <- sapply(headerNames, length)
      otherNames <- setdiff(header, unlist(headerNames))
      headerNames[nheaders == 0L] <- list(otherNames)
      nheaders <- sapply(headerNames, length)
      headerNames <- unlist(headerNames)
      if (anyDuplicated(headerNames))
        stop("Redundant role information: ", paste(role, collapse=", "))
      role <- rep(role, nheaders)
      headerInd <- match(headerNames, header)
      map <- integer(length(header))
      map[headerInd] <- match(role, cn)
      map - 1L
    }
    resolvedRoles <- tapply(cn, factor(dataRoles, unique(dataRoles)),
                            resolveRole, simplify=FALSE)
    roles[names(resolvedRoles)] <- resolvedRoles
  } else {
    isEditable <- colnames(df) %in% editable
    inds <- seq_len(ncol(df)) - 1L
    roles$display <- ifelse(isEditable, -1L, inds)
    roles$edit <- ifelse(isEditable, inds, -1L)
    header <- colnames(df)
  }
  attrs <- attributes(df)
  getHeaderRoles <- function(attrName, display = attrs[[attrName]]) {
    headerRoles <- createRoleList()
    headerRoles$display <- display
    if (useRoles) {
      attrPrefix <- paste(attrName, ".", sep = "")
      headerAttrs <- grep(attrPrefix, names(attrs), fixed=TRUE, value=TRUE)
      headerRoles[getRoleNames(headerAttrs)] <- attrs[headerAttrs]
    }
    headerRoles
  }
  rowRoles <- getHeaderRoles("row.names")
  colRoles <- getHeaderRoles("names", header)    
  .Call("qt_qsetDataFrame", model, df, roles, rowRoles, colRoles,
        PACKAGE="qtbase")
  model
}

##' @rdname DataFrameModel
qdataFrame <- function(model) {
  stopifnot(inherits(model, "DataFrameModel"))
  .Call("qt_qdataFrame", model, PACKAGE="qtbase")
}

### TODO: add setters

quseRoles <- function(model) .Call("qt_quseRoles", model, PACKAGE="qtbase")
qeditable <- function(model) .Call("qt_qeditable", model, PACKAGE="qtbase")
