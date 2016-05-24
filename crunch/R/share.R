setMethod("permissions", "CrunchDataset", function (x) {
    perm_url <- shojiURL(x, "catalogs", "permissions")
    return(PermissionCatalog(crGET(perm_url)))
})

##' @rdname describe-catalog
##' @export
setMethod("emails", "PermissionCatalog", function (x) getIndexSlot(x, "email"))

is.editor <- function (x) {
    out <- vapply(index(x), function (a) {
            isTRUE(a[["dataset_permissions"]][["edit"]])
        }, logical(1), USE.NAMES=FALSE)
    names(out) <- emails(x)
    structure(return(out))
}

userCanEdit <- function (email, dataset) {
    e <- is.editor(permissions(dataset))
    return(ifelse(email %in% names(e), e[email], FALSE))
}

iCanEdit <- function (dataset) {
    is.editor(permissions(dataset)[userURL()])
}

userCanView <- function (email, dataset) {
    email %in% emails(permissions(dataset))
}

##' Share a dataset
##'
##' @param dataset a CrunchDataset
##' @param users character: email address(es) or URLs of the users or teams with
##' whom to share the dataset. If there is no Crunch user associated with an
##' email, an invitation will be sent.
##' @param edit logical: should the specified user(s) be given edit privileges
##' on the dataset? Default is \code{FALSE}. \code{edit} can be a single value
##' or, if inviting multiple users, a vector of logical values of equal length
##' of the number of emails given.
##' @param notify logical: should users who are getting new privileges on this
##' dataset be sent an email informing them of this fact? Default is
##' \code{TRUE}.
##' @return Invisibly, the dataset.
##' @export
share <- function (dataset, users, edit=FALSE, notify=TRUE) {
    perms <- permissions(dataset)
    if (length(edit) == 1) {
        edit <- rep(edit, length(users))
    }
    if (length(edit) != length(users)) {
        halt("Must supply `edit` permissions of equal length as the number of `emails` supplied")
    }
    if (!any(edit) && all(emails(perms)[is.editor(perms)] %in% users)) {
        halt("Cannot remove editor from the dataset without specifying another")
    }
    payload <- lapply(edit,
        function (e) list(dataset_permissions=list(edit=e, view=TRUE)))
    names(payload) <- users
    payload$send_notification <- notify
    if (notify) {
        payload$url_base <- passwordSetURLTemplate()
        payload$dataset_url <- webURL(dataset)
    }
    payload <- toJSON(payload)
    crPATCH(self(perms), body=payload)
    invisible(dataset)
}

passwordSetURLTemplate <- function () {
    absoluteURL("/password/change/${token}/", getOption("crunch.api"))
}

## TODO: test and release this
# shareDataset <- function (x, emails, notify=TRUE) {
#     ## Share one more more datasets without loading them
#     dscat <- active(datasetCatalog())
#     if (!is.numeric(x)) {
#         x <- selectDatasetFromCatalog(x, dscat, strict=TRUE)
#     }
#     dsurls <- urls(dscat)[x]
#     perm_urls <- vapply(dsurls, function (x) absoluteURL("permissions/", x),
#         character(1))
#
#     payload <- sapply(emails,
#             function (x) list(dataset_permissions=list(edit=FALSE, view=TRUE)),
#             simplify=FALSE)
#     payload$send_notification <- notify
#     payload <- toJSON(payload)
#
#     out <- lapply(perm_urls, function (x) {
#         message("Sharing ", x)
#         crPATCH(x, body=payload)
#     })
#     invisible(out)
# }
