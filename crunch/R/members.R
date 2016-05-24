.backstopUpdate <- function (x, i, j, value) {
    ## Backstop error so you don't get "Object of class S4 is not subsettable"
    halt(paste("Cannot update", class(x), "with type", class(value)))
}

##' @rdname catalog-extract
##' @export
setMethod("[[<-", c("MemberCatalog", "ANY", "missing", "ANY"), .backstopUpdate)

##' @rdname catalog-extract
##' @export
setMethod("[[<-", c("MemberCatalog", "character", "missing", "NULL"),
    function (x, i, j, value) {
        ## Remove the specified user from the catalog
        payload <- sapply(i, function (z) NULL, simplify=FALSE)
        crPATCH(self(x), body=toJSON(payload))
        return(refresh(x))
    })

##' @rdname teams
##' @export
setMethod("members<-", c("CrunchTeam", "MemberCatalog"), function (x, value) {
    ## TODO: something
    ## For now, assume action already done in other methods, like NULL
    ## assignment above.
    return(x)
})

##' @rdname teams
##' @export
setMethod("members<-", c("CrunchTeam", "character"), function (x, value) {
    payload <- sapply(value,
        function (z) emptyObject(),
        simplify=FALSE)
    crPATCH(self(members(x)), body=toJSON(payload))
    return(refresh(x))
})

##' @rdname describe-catalog
##' @export
setMethod("emails", "ShojiCatalog", function (x) getIndexSlot(x, "email"))
