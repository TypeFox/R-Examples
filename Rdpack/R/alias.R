                                                      # todo: are aliases updated for classes?

update_aliases_tmp <- function(rdo, package = NULL){#only for methods, currently; only appends
    fname <- .get.name_content(rdo)$short # name of the generic

    rdoaliases <- toolsdotdotdot.Rd_get_metadata(rdo, "alias")  # aliases currently in rdo

    curtxt <- get_sig_text(rdo, package = package)       # process signatures described in rdo
    ucur <- parse_usage_text(curtxt)

    methaliases <- sapply(ucur, function(x) .methsig2alias(fname, x$defaults))

    new_aliases <- methaliases[ !(methaliases %in% rdoaliases) ]
    if(length(new_aliases) > 0)                                           # update the aliases
        for(alias in new_aliases){
            rdo <- Rdo_insert(rdo, char2Rdpiece(alias, "alias"))
        }
    rdo
}

.methsig2alias <- function(name, sig){
    res <- gsub("\"", "", paste(c(name, sig), collapse=","))
    paste(res, "-method", sep="")   # no comma after the last arg
}
