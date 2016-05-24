userURL <- function () rootURL("user")

getUser <- function (x=userURL()) {
    ShojiObject(crGET(x))
}

getUserCatalog <- function (x=sessionURL("users")) {
    UserCatalog(crGET(x))
}

getAccount <- function (x=rootURL("account", getUser())) {
    ShojiObject(crGET(x))
}

##' Find all users on your account
##'
##' @param x URL of the user catalog. Default is the right thing; you shouldn't
##' specify one
##' @return a \code{UserCatalog}
##' @export
getAccountUserCatalog <- function (x=shojiURL(getAccount(), "catalogs", "users")) {
    UserCatalog(crGET(x))
}

##' @rdname describe-catalog
##' @export
setMethod("emails", "UserCatalog", function (x) getIndexSlot(x, "email"))

invite <- function (email, name=NULL, notify=TRUE, id_method="pwhash",
                    advanced=FALSE, admin=FALSE, ...) {
    payload <- list(
        email=email,
        send_invite=notify,
        id_method=id_method,
        account_permissions=list(
            alter_users=isTRUE(admin),
            create_datasets=isTRUE(advanced)),
        ...)
    if (!is.null(name)) {
        payload$first_name <- name
    }
    if (id_method == "pwhash") {
        payload$url_base <- "/password/change/${token}/"
    }

    url <- shojiURL(getAccount(), "catalogs", "users")
    return(crPOST(url, body=toJSON(list(element="shoji:entity", body=payload))))
}
