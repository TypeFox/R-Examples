chooseTaxa <- function(object, ...) {
    UseMethod("chooseTaxa")
}

chooseTaxa.default <- function(object, n.occ = 1, max.abun = 0,
                               type = c("AND","OR"), value = TRUE,
                               na.rm = FALSE, ...) {
    if(missing(type))
        type <- "AND"
    type <- match.arg(type)
    ## issue a warning if any values in object are NA
    if(!na.rm && any(is.na(object)))
        warning("'NA's present in data; results may not be what you expect.\nConsider using 'na.rm = TRUE'.")
    occ.want <- colSums(object > 0, na.rm = na.rm) >= n.occ
    abun.want <- apply(object, 2, max, na.rm = na.rm) >= max.abun
    want <- if(isTRUE(all.equal(type, "AND"))) {
        occ.want & abun.want
    } else {
        occ.want | abun.want
    }
    if(value) {
        rname <- rownames(object)
        cname <- colnames(object)
        object <- object[, want, drop = FALSE]
        rownames(object) <- rname
        colnames(object) <- cname[want]
    } else {
        object <- want
    }
    object
}
