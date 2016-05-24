#' The SerialisableObject class
#' 
#' This reference class extends the standard \code{\linkS4class{envRefClass}}
#' class, adding a function for simple serialisation of the data fields of an
#' object, and one for finding all of the methods available for an object. A
#' serialised object may be deserialised using the
#' \code{\link{deserialiseReferenceObject}} function.
#' 
#' @export
SerialisableObject <- setRefClass("SerialisableObject", methods=list(
    fields = function ()
    {
        "Retrieve a list of all field names"
        allFieldNames <- names(.self$getRefClass()$fields())
        if (is.null(allFieldNames))
            return (NULL)
        else
            return (allFieldNames[!(allFieldNames %~% "\\.$")])
    },
    
    methods = function () { return (.self$getRefClass()$methods()) },
    
    serialise = function (file = NULL)
    {
        "Serialise the object to a list or file"
        originalClass <- class(.self)
        originalPackage <- attr(originalClass,"package")
        attributes(originalClass) <- NULL
        
        # Fields with names ending in "." will not be returned, and therefore not be serialised
        fields <- .self$fields()
        serialisedObject <- list()

        for (field in fields)
        {
            fieldValue <- get(field)

            if (is(fieldValue, "SerialisableObject"))
                serialisedObject <- c(serialisedObject, list(fieldValue$serialise(NULL)))
            else
                serialisedObject <- c(serialisedObject, list(fieldValue))
        }

        names(serialisedObject) <- fields
        attr(serialisedObject, "originalClass") <- originalClass
        attr(serialisedObject, "originalPackage") <- originalPackage

        if (!is.null(file))
            save(serialisedObject, file=ensureFileSuffix(file,"Rdata"))
        
        invisible(serialisedObject)
    }
))

setMethod("show", "SerialisableObject", function (object)
{
    if ("summarise" %in% object$methods())
    {
        summaryList <- object$summarise()
        if (is.list(summaryList) && all(c("labels","values") %in% names(summaryList)))
            printLabelledValues(summaryList$labels, summaryList$values)
        else if (!is.null(names(summaryList)))
            printLabelledValues(names(summaryList), as.character(summaryList))
    }
    else
        cat(paste("An object of class \"", class(object)[1], "\"\n", sep=""))
})

.NilObject <- SerialisableObject$new()

#' The nil object
#' 
#' The nil object is an empty object of class \code{\link{SerialisableObject}}.
#' It can be used as a placeholder where such an object of this class, or one
#' of its subclasses, is required. It serialises to the empty list.
#' 
#' @param object Any object.
#' @return \code{nilObject} returns the nil object. \code{is.nilObject} returns
#'   \code{TRUE} if its argument is identical to the nil object, or if it is
#'   equivalent in the sense of serialising to an identical result.
#' 
#' @author Jon Clayden
#' @seealso \code{\link{SerialisableObject}}
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @export
nilObject <- function ()
{
    return (.NilObject)
}

#' @rdname nilObject
#' @export
is.nilObject <- function (object)
{
    if (identical(object, .NilObject))
        return (TRUE)
    else if (identical(object$serialise(), .NilObject$serialise()))
        return (TRUE)
    else
        return (FALSE)
}

#' Reference object serialisation and deserialisation
#' 
#' Rather than using R's \code{\link{save}} and \code{\link{load}} functions
#' directly for reference objects, TractoR uses the
#' \code{\linkS4class{SerialisableObject}} class and these functions to save
#' and load objects. The main difference is that this approach stores only the
#' data in the object, and not the functions which operate on them. This helps
#' backward compatibility when new member functions are added.
#' 
#' The \code{serialiseReferenceObject} function, or the \code{serialise} member
#' function of the \code{\link{SerialisableObject}} class can be used to create
#' and/or \code{\link{save}} a version of an object which contains a
#' hierarchical representation of the data embedded in it. These serialised
#' objects are standard R lists, with an \code{"originalClass"} attribute
#' describing the class of the original object. The
#' \code{deserialiseReferenceObject} function can be used to deserialise them.
#' Custom deserialisers can be specified using \code{registerDeserialiser},
#' typically for legacy classes.
#' 
#' Note that this should generally NOT be used as the primary mechanism for
#' saving and loading \code{\link{MriImage}} objects. Saving to standard
#' NIfTI/Analyze format is usually preferable, and can be done using
#' \code{\link{writeImageFile}}.
#' 
#' @param object For \code{serialiseReferenceObject}, a list or object
#'   inheriting from \code{\linkS4class{SerialisableObject}}. For other
#'   functions, an object in (raw) serialised form. See Details.
#' @param expectedClass A class name which the object is expected to inherit.
#'   Any class is acceptable if this parameter is \code{NULL}.
#' @param file A file name to deserialise from.
#' @param raw If \code{TRUE}, the raw serialised object is returned; otherwise
#'   the object is converted back to its original class.
#' @param className A string naming a class to be handled by the specified
#'   deserialiser.
#' @param deserialiser A function taking as its argument a list of serialised
#'   fields, and returning a suitable deserialised object.
#' @return \code{isDeserialisable} returns \code{TRUE} if the \code{object} is
#'   deserialisable and inherits from the specified class.
#'   \code{deserialiseReferenceObject} returns a raw or reconstituted object
#'   after deserialisation.
#' 
#' @author Jon Clayden
#' @seealso \code{\linkS4class{SerialisableObject}}, \code{\link{save}},
#' \code{\link{load}}, \code{\link{writeImageFile}}.
#' @references Please cite the following reference when using TractoR in your
#' work:
#' 
#' J.D. Clayden, S. Muñoz Maniega, A.J. Storkey, M.D. King, M.E. Bastin & C.A.
#' Clark (2011). TractoR: Magnetic resonance imaging and tractography with R.
#' Journal of Statistical Software 44(8):1-18.
#' \url{http://www.jstatsoft.org/v44/i08/}.
#' @aliases serialisation
#' @rdname serialisation
#' @export
isDeserialisable <- function (object, expectedClass = NULL)
{
    if (is.null(object) || is.null(attr(object,"originalClass")))
        return (FALSE)
    else if (!is.null(expectedClass) && !(expectedClass %in% attr(object,"originalClass")))
        return (FALSE)
    else
        return (TRUE)
}

#' @rdname serialisation
#' @export
serialiseReferenceObject <- function (object, file = NULL)
{
    if (is(object, "SerialisableObject"))
        serialisedObject <- object$serialise()
    else if (is.list(object))
        serialisedObject <- lapply(object, serialiseReferenceObject)
    else
        report(OL$Error, "Object to serialise must be a list or a SerialisableObject")
    
    if (!is.null(file))
        save(serialisedObject, file=ensureFileSuffix(file,"Rdata"))
    
    invisible(serialisedObject)
}

#' @rdname serialisation
#' @export
deserialiseReferenceObject <- function (file = NULL, object = NULL, raw = FALSE)
{
    if (is.null(object))
    {
        if (is.null(file))
            report(OL$Error, "Either a file or raw deserialised object must be specified")
        object <- get(load(ensureFileSuffix(file,"Rdata")))
    }
    
    if (!isDeserialisable(object))
    {
        if (is.list(object))
            return (invisible(lapply(object, function(x) deserialiseReferenceObject(object=x,raw=raw))))
        else
            report(OL$Error, "The specified object or file is not deserialisable")
    }
    else if (raw)
        return (invisible(object))
    
    fields <- lapply(object, function (field) {
        if (isDeserialisable(field))
            return (deserialiseReferenceObject(object=field))
        else
            return (field)
    })
    names(fields) <- names(object)
    
    packageName <- attr(object, "originalPackage")
    if (!is.null(packageName) && !(paste("package",packageName,sep=":") %in% search()))
        require(packageName, character.only=TRUE, quietly=TRUE)
    
    className <- attr(object, "originalClass")
    if (className %in% names(.Workspace$deserialisers))
        finalObject <- .Workspace$deserialisers[[className]](fields)
    else
    {
        if (!is.null(packageName))
            class <- getRefClass(className, where=as.environment(paste("package",packageName,sep=":")))
        else
            class <- getRefClass(className)
        finalObject <- do.call(class$new, fields)
    }
    
    invisible(finalObject)
}

#' @rdname serialisation
#' @export
registerDeserialiser <- function (className, deserialiser)
{
    if (!is.character(className) || length(className) != 1)
        report(OL$Error, "Class name should be specified as a character string")
    
    deserialiser <- match.fun(deserialiser)
    .Workspace$deserialisers[[className]] <- deserialiser
    
    invisible(NULL)
}
