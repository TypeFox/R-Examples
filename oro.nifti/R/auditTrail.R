#' @name Audit Trails
#' @title Facilitate the Creation and Modification of Audit Trails
#' 
#' @description Facilitate the creation and modification of audit trails for NIfTI class
#' objects.
#' @details The function \code{oro.nifti.info} is used to find the \code{ecode} or the
#' XML \code{namespace} relevant to the audit trail.
#' 
#' The function \code{enableAuditTrail} is turned \dQuote{off} by default to
#' minimize package dependencies.  Should one wish to turn \dQuote{on} the
#' audit trail functionality, then one should set the option
#' \code{NIfTI.audit.trail} to \code{TRUE} and call the function
#' \code{enableAuditTrail}.  Setting the option \code{NIfTI.audit.trail} to
#' \code{FALSE} will disable the audit trail.
#' 
#' The function \code{newAuditTrail} returns an \code{XMLAbstractNode}
#' representing the root node of an audit trail.  This is mostly intended as an
#' internal function.
#' 
#' The function \code{niftiExtensionToAuditTrail} takes an object representing
#' a NIfTI object, casts it as a \code{niftiAuditTrail} and checks if there is
#' an extension (a \code{niftiExtensionSection}) with \code{ecode} equal to
#' \code{oro.nifti.info("ecode")}; i.e. has a extension with data representing
#' a serialized audit trail.  The function will then strip the object of this
#' extension parsing the serialized \code{edata} into an audit trail and adding
#' a \sQuote{read} event to the trail.
#' 
#' The function \code{niftiAuditTrailToExtension} takes a
#' \code{niftiAuditTrail} and returns a \code{niftiExtensionSection} with
#' \code{edata} containing the serialized form of the audit trail after adding
#' a \sQuote{saved} event to the trail.
#' 
#' The function \code{niftiAuditTrailSystemNodeEvent} adds an element with name
#' equal to \code{type} to the \code{trail}.  It uses the
#' \code{niftiAuditTrailSystemNode} function to create the node.
#' 
#' The function \code{niftiAuditTrailSystemNode} is an internal function
#' creating an \code{XMLAbstractNode} element with name \code{type} and
#' attributes giving information about the R system and library.  The
#' \code{filename} and \code{call} will also be added as attributes if
#' available.
#' 
#' The function \code{niftiAuditTrailEvent} adds an element with name
#' \code{event} to the \code{trail}.  The arguments \code{type},
#' \code{filename}, \code{call} are added as attributes and the \code{comment}
#' is the text value of the element.
#' 
#' The function \code{niftiAuditTrailCreated} will create a new audit trail
#' containing a system node element \code{created} with the child
#' \code{history} with the contents \code{history}.  If the last element of the
#' \code{history} given is an \code{event} with \code{type="processing"}, then
#' this node will be removed from the \code{history} and its \code{call}
#' attribute will be used as the value of the \code{call} attribute on the
#' \code{created} node.
#' 
#' The function \code{getLastCallWithName} will search the call stack for a
#' call of the function \code{functionName}, returning last call to that
#' function if possible.  It will default to the call of the function which
#' called the function which called getLastCallWithName if there was no such
#' call (and if there was no such call it will return the call of itself).
#' 
#' @aliases oro.nifti.info enableAuditTrail newAuditTrail
#' niftiExtensionToAuditTrail niftiAuditTrailToExtension
#' niftiAuditTrailSystemNode niftiAuditTrailSystemNodeEvent
#' niftiAuditTrailCreated niftiAuditTrailEvent getLastCallWithName
#' @param nim is an object of class \code{niftiAuditTrail} or can be converted
#' to such.
#' @param workingDirectory The working directory associated with the
#' \sQuote{filename}.
#' @param filename The filename associated with the nifti object.
#' @param call A \code{call}, function name in the call-stack or a string.
#' @param type An identifier to add some meaning to the event.
#' @param trail The \code{XMLAbstractNode} representing the audit trail or the
#' \code{niftiAuditTrail} object with a trail that will be amended.
#' @param comment Some textual comment
#' @param history An \code{XMLAbstractNode} to store historical events for
#' inclusion in the \sQuote{trail}.
#' @param functionName The name of a function on the call stack.
#' @note These functions are mostly intended to be used internally in order to
#' document the changes that occur to NIfTI objects due to functions that are
#' audit-trail aware.  However, as the precise manner in which these functions
#' are used is not documented anywhere else, we shall proceed to describe which
#' functions are audit-trail aware and how they interact with the audit trail.
#' 
#' \code{as.nifti} and its S4 alias \code{as(nim, "nifti")} will always produce
#' \code{niftiAuditTrail} objects if the functionality is turned on.  The
#' function \code{niftiAuditTrailCreated} will be used and if an exemplar
#' object is provided (e.g., \code{as.nifti(array, niftiExemplar)}) then the
#' trail of the exemplar will be used as the \code{history}.
#' 
#' \code{readNIfTI} and \code{writeNIfTI} also always produce
#' \code{niftiAuditTrail} objects if the functionality is turned on.  The
#' functions \code{niftiExtensionToAuditTrail} and
#' \code{niftiAuditTrailToExtension} are used internally by these functions to
#' facilitate this behaviour.
#' @author Andrew Thornton \email{zeripath@@users.sourceforge.net} and Brandon
#' Whitcher \email{bwhitcher@@gmail.com}
#' 
#' @examples 
#' ## A good example of the use of these functions is shown by this
#' ## wrapper function which takes a function fun(nim, ...) returning
#' ## lists of arrays which are nifti-ized using as(...)
#' options("niftiAuditTrail"=TRUE)
#' enableAuditTrail()
#' 
#' wrapper <- function(functionToWrap, nameOfCallingFunction, nim, ...) {
#'   if (!is(nim, "nifti")) 
#'     nim <- as(nim, "nifti")
#'   
#'   if (is(nim, "niftiAuditTrail")) {
#'     ## This will force as(...) to set the call which created the
#'     ## results to the calling function's call rather than
#'     ## as(result, nifti) as it would otherwise do
#'     slot(nim, "trail") <- niftiAuditTrailEvent(slot(nim, "trail"), "processing",
#'                                       nameOfCallingFunction)
#'   }
#'   
#'   result <- functionToWrap(nim, ...)
#'   as(result, "nifti") <- nim
#'   return(result)
#' }
#' 
#' ## An example of how wrapper is used follows:
#' functionToWrap <- function(ignored, x, y) {
#'   return (array(1, dim=c(x,y)))
#' }
#' 
#' ## The nifti-ized form
#' niftiizedForm <- function(nim,...) {
#'   return(wrapper(functionToWrap, "niftiizedForm", nim, ...))
#' }
#'  
#' \dontrun{
#'   if (isTRUE(getOption("niftiAuditTrail"))) {
#'     print(slot(as.nifti(functionToWrap(nifti(), 4, 4), nifti()), "trail"))
#'     print(slot(niftiizedForm(nifti(), 4, 4), "trail"))
#'   }
#' }
#' @rdname audit_trail
#' @export
oro.nifti.info <- function(type) {
  switch(type,
         ecode = 6,
         namespace = "http://www.dcemri.org/namespaces/audit-trail/1.0")
}
#' @rdname audit_trail
#' @export
enableAuditTrail <- function() {
  if (requireNamespace("XML", quietly=TRUE)) {
    options("niftiAuditTrail" = TRUE)
  }
}

## Look back through call stack for the last call with the name functionName
## otherwise return the function that called the function that called us
#' @rdname audit_trail
#' @export
getLastCallWithName <- function(functionName) {
  theCalls <- sys.calls()
  correctCalls <- which(sapply(theCalls, function(x) x[[1]] == functionName))
  if (length(correctCalls) == 0) {
    return(theCalls[max(1, length(theCalls) - 2)])
  }
  return(theCalls[[max(correctCalls)]])
}
#' @rdname audit_trail
#' @export
newAuditTrail <- function() {
  if (isTRUE(getOption("niftiAuditTrail")) && requireNamespace("XML", quietly=TRUE)) {
    trail <- XML::newXMLNode("audit-trail",
                     namespaceDefinitions=oro.nifti.info("namespace"))
    return(trail)
  }
}

.listToNodes <- function(thelist) {
  return(lapply(names(thelist), function(x) { XML::newXMLNode(x, thelist[x]) }))
}
#' @rdname audit_trail
#' @export
niftiExtensionToAuditTrail <- function(nim, workingDirectory=NULL,
                                       filename=NULL, call=NULL) {
  if (isTRUE(getOption("niftiAuditTrail")) && requireNamespace("XML", quietly=TRUE)) {
    if (!is(nim, "niftiAuditTrail")) {
      nim <- as(nim, "niftiAuditTrail")
    }
    ## We enforce that there is a single extension with ecode == oro.nifti.ecode
    ecodes <- lapply(nim@extensions, function(x) x@ecode)
    oei <- which(ecodes == oro.nifti.info("ecode"))
    
    if (length(oei) == 0) {
      audit.trail(nim) <-
        niftiAuditTrailCreated(workingDirectory=workingDirectory,
                               filename=filename, call=call)
    } else {
      ## One Trail
      if (length(oei) > 1) {
        warning("Found more than one extension with ecode ==",
                oro.nifti.info("ecode"), ", Appending to last trail only")
        oei <- oei[length(oei)]
      }
      oe <- nim@extensions[[oei]]@edata
      nim@extensions[[oei]] <- NULL
      audit.trail(nim) <-
        niftiAuditTrailSystemNodeEvent(XML::xmlRoot(XML::xmlParse(iconv(oe, to="UTF-8"),
                                                                  asText=TRUE)),
                                       type="read",
                                       workingDirectory=workingDirectory,
                                       filename=filename, call=call)
    }
  }
  return(nim)
}

niftiAuditTrailToExtension <- function(nim, workingDirectory=NULL,
                                       filename=NULL, call=NULL) {
  if (isTRUE(getOption("niftiAuditTrail")) && is(nim, "niftiAuditTrail") &&
        requireNamespace("XML", quietly=TRUE)) {
    sec <- new("niftiExtensionSection")
    sec@ecode <- oro.nifti.info("ecode")
    audit.trail(nim) <-
      niftiAuditTrailSystemNodeEvent(audit.trail(nim), "saved",
                                     workingDirectory=workingDirectory,
                                     filename=filename, call=call)
    ## Serialize the XML to sec@edata
    sec@edata <- XML::saveXML(audit.trail(nim))
    ## Fix the esize to be congruent to 0 mod 16
    sec@esize <- nchar(sec@edata, type="bytes") + 8
    sec@esize <- (-sec@esize %% 16) + sec@esize
    return(sec)
  }
}
#' @rdname audit_trail
#' @export
niftiAuditTrailSystemNode <- function(type="system-info",
                                      workingDirectory=NULL, filename=NULL,
                                      call=NULL) {
  if (isTRUE(getOption("niftiAuditTrail")) && requireNamespace("XML", quietly=TRUE)) {
    if (is(call, "character") && is(try(get(call, mode="function"),
                                        silent=TRUE), "function")) {
      call <- as.character(as.expression(getLastCallWithName(call)))
    }
    if (is(call, "call")) {
      call <- as.character(as.expression(call))
    }
    currentDateTime <- format(Sys.time(), "%a %b %d %X %Y %Z")
    children <- .listToNodes(c("workingDirectory"=workingDirectory,
                               "filename"=filename, "call"=call))
    sysinfo <-
      .listToNodes(c("r-version"=R.version$version.string,
                     "date"=currentDateTime,
                     "user"=Sys.getenv("LOGNAME"),
                     "package-version"=packageDescription("oro.nifti")$Version))
    if (is.null(children)) {
      children <- sysinfo
    } else {
      children <- c(children, XML::newXMLNode("system", sysinfo))
    }
    system <- XML::newXMLNode(type, children)
    return(system)
  }
}
#' @rdname audit_trail
#' @export
niftiAuditTrailSystemNodeEvent <- function(trail, type=NULL, call=NULL,
                                           workingDirectory=NULL,
                                           filename=NULL, comment=NULL) {
  if (isTRUE(getOption("niftiAuditTrail")) && requireNamespace("XML", quietly=TRUE)) {
    if (is(trail, "niftiAuditTrail")) {
      return(niftiAuditTrailSystemNodeEvent(audit.trail(trail), type, call,
                                            workingDirectory, filename,
                                            comment))
    }
    eventNode <- niftiAuditTrailSystemNode(type=type, call=call,
                                           workingDirectory=workingDirectory,
                                           filename=filename)
    if (!is.null(comment))
      eventNode <- XML::addChildren(eventNode, XML::newXMLTextNode(comment))
    trail <- XML::addChildren(trail, eventNode)
    return(trail)
  }
}
#' @rdname audit_trail
#' @export
niftiAuditTrailCreated <- function(history=NULL, call=NULL,
                                   workingDirectory=NULL, filename=NULL) {
  if (isTRUE(getOption("niftiAuditTrail")) && requireNamespace("XML", quietly=TRUE)) {
    if (is(history, "niftiAuditTrail")) {
      return(niftiAuditTrailCreated(audit.trail(history), call,
                                    workingDirectory, filename))
    } else {
      trail <- newAuditTrail()
      if (is.null(history) || length(XML::xmlChildren(history)) == 0) {
        created <-
          niftiAuditTrailSystemNode("created",
                                    "workingDirectory"=workingDirectory,
                                    "filename"=filename, "call"=call)
      } else {
        historyChildren <- XML::xmlChildren(history)
        
        lastEvent <- historyChildren[[length(historyChildren)]]
        
        if (XML::xmlName(lastEvent) == "event" &&
              XML::xmlValue(lastEvent[["type"]]) == "processing") {
          ## We are in some processing history;
          ## the given call is not the correct call
          call <- XML::xmlValue(lastEvent[["call"]])
          historyChildren <- historyChildren[1:(length(historyChildren) - 1)]
        } 
        
        created <-
          niftiAuditTrailSystemNode("created",
                                    "workingDirectory"=workingDirectory,
                                    "filename"=filename, "call"=call)
        historyNode <- XML::newXMLNode("history")
        ## OK, serialize and reParse the history
        historyChildren <-
          lapply(historyChildren,
                 function(x) {
                   XML::xmlRoot(XML::xmlParse(iconv(XML::saveXML(x), to="UTF-8"), 
                                              asText=TRUE))
                 })
        historyNode <- XML::addChildren(historyNode, historyChildren)
        created <- XML::addChildren(created, historyNode)
      }
      trail <- XML::addChildren(trail, created) 
      return(trail)
    }
  }
}
#' @rdname audit_trail
#' @export
niftiAuditTrailEvent <- function(trail, type=NULL, call=NULL, comment=NULL) {
  if (isTRUE(getOption("niftiAuditTrail")) && requireNamespace("XML", quietly=TRUE)) {
    if (is(trail,"niftiAuditTrail")) {
      return(niftiAuditTrailEvent(audit.trail(trail), type, call, comment))
    }
    if (is(call, "character") &&
          is(try(get(call, mode="function"), silent=TRUE), "function")) {
      call <- as.character(as.expression(getLastCallWithName(call)))
    }
    if (is(call, "call")) {
      call <- as.character(as.expression(call))
    }
    currentDateTime <- format(Sys.time(), "%a %b %d %X %Y %Z")
    eventNode <- XML::newXMLNode("event",
                                 .listToNodes(c("type"=type, "call"=call,
                                                "date"=currentDateTime, 
                                                "comment"=comment, 
                                                "user"=Sys.getenv("LOGNAME"))))
    trail <- XML::addChildren(trail, eventNode)
    return(trail)
  }
}
