##' @rdname catalog-extract
##' @export
setMethod("[[", c("ProjectCatalog", "character"), function (x, i, ...) {
    w <- whichNameOrURL(x, i)
    x[[w]]
})

##' @rdname catalog-extract
##' @export
setMethod("[[", c("ProjectCatalog", "ANY"), function (x, i, ...) {
    b <- callNextMethod(x, i, ...)
    if (is.null(b)) return(NULL)
    CrunchProject(index_url=self(x), entity_url=urls(x)[i],
        body=b)
})

##' @rdname catalog-extract
##' @export
setMethod("[[<-", c("ProjectCatalog", "character", "missing", "list"),
    function (x, i, j, value) {
        if (i %in% names(x)) {
            ## TODO: update team attributes
            halt("Cannot (yet) modify project attributes")
        } else {
            ## Creating a new team
            u <- crPOST(self(x), body=toJSON(list(name=i)))
            x <- refresh(x)
            ## Add members to team, if given
            if (!is.null(value[["members"]]))
            members(x[[i]]) <- value[["members"]]
            return(x)
        }
    })

##' @rdname catalog-extract
##' @export
setMethod("[[<-", c("ProjectCatalog", "character", "missing", "CrunchProject"),
    function (x, i, j, value) {
        ## TODO: something
        ## For now, assuming that modifications have already been persisted
        ## by other operations on the team entity (like members<-)
        return(x)
    })

##' @rdname teams
##' @export
setMethod("members", "CrunchProject", function (x) {
    MemberCatalog(crGET(shojiURL(x, "catalogs", "members")))
})

##' @rdname teams
##' @export
setMethod("members<-", c("CrunchProject", "MemberCatalog"), function (x, value) {
    ## TODO: something
    ## For now, assume action already done in other methods, like NULL
    ## assignment above.
    return(x)
})

##' @rdname teams
##' @export
setMethod("members<-", c("CrunchProject", "character"), function (x, value) {
    value <- setdiff(value, emails(members(x)))
    if (length(value)) {
        payload <- sapply(value,
            function (z) emptyObject(),
            simplify=FALSE)
        crPATCH(self(members(x)), body=toJSON(payload))
    }
    return(x)
})

##' @rdname tuple-methods
##' @export
setMethod("entity", "CrunchProject", function (x) {
    return(ProjectEntity(crGET(x@entity_url)))
})

##' @rdname delete
##' @export
setMethod("delete", "CrunchProject", function (x, confirm=requireConsent(), ...) {
    prompt <- paste0("Really delete project ", dQuote(name(x)), "? ",
        "This cannot be undone.")
    if (confirm && !askForPermission(prompt)) {
        halt("Must confirm deleting project")
    }
    u <- self(x)
    out <- crDELETE(u)
    dropCache(absoluteURL("../", u))
    invisible(out)
})

##' A project's datasets
##' @param x a \code{CrunchProject}
##' @param value \code{CrunchDataset} for the setter
##' @return An obect of class \code{DatasetCatalog}. The setter returns the
##' catalog with the given dataset added to it (via changing its owner to be
##' the specified project, \code{x}).
##' @name project-datasets
##' @export
datasets <- function (x) DatasetCatalog(crGET(shojiURL(x, "catalogs", "datasets")))

##' @rdname project-datasets
##' @export
`datasets<-` <- function (x, value) {
    stopifnot(inherits(x, "CrunchProject"))
    value <- setEntitySlot(value, "owner", self(x))
    return(refresh(x))
}
