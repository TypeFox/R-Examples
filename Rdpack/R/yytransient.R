                                        # 2013-12-01 new; need testing
.ascharRd <- function(rdo){
    if(is.character(rdo)) # note: applying as.character.Rd, as below gives error in this case
        res <- as.vector(rdo)   # to drop attributes;
    else{
        class(rdo) <- "Rd"
        res <- as.character(rdo)
    }

    paste(res, collapse="")  # todo: argument to make this optional?
}

.safeRdtag <- function(x, nulltag = ""){         # replaces NULL with nulltag
    res <- attr(x, "Rd_tag")
    if(is.null(res))
        res <- nulltag
    res
}

                                        # a replacement for tools:RdTags; takes
                                        # 2013-12-08 renamed from .top_RdTags
Rdo_tags <- function(rdo, nulltag = ""){            # note: absent Rd_tag's are returned as ""
    locf <- function(x)
        .safeRdtag(x, nulltag)

    res <- sapply(rdo, locf)
    if (length(res)==0)
        res <- character()
    res
}

Rdo_which <- function(rdo, fun){  # fun - predicate
    wrk <- sapply(rdo, fun)
    which(wrk)
}

                                        # 2013-12-08 renamed from .whichtageq
Rdo_which_tag_eq <- function(rdo, tag){
    tags <- Rdo_tags(rdo)
    pos <- which(tags == tag)

    pos
}

                                        # 2013-12-08 renamed from .whichtagin
Rdo_which_tag_in <- function(rdo, tags){                 # TODO: needs testing!
    alltags <- Rdo_tags(rdo)
    pos <- which(sapply(alltags, function(x) .tag_in(x,tags) ))
    pos
}


Rdo_get_item_labels <- function(rdo){
    wrk <- Rdo_locate(rdo,
                         f = function(x){ # attr(x,"Rd_tag") == "\\item"
                             .tag_eq(x, "\\item")
                             },
                         lists = TRUE,
                         nested = FALSE,
                         pos_only = function(x) if(length(x) > 0) .ascharRd(x[[1]])
                                                else "" #in \itemize items do not have labels
                         )

    sapply(wrk, function(x) x$value)
}

Rdo_get_argument_names <- function(rdo){
               # 2013-12-08 was: indx <- Rdo_locate_predefined_section(rdo, "\\arguments")
               # todo: shouldn't this be using some "locate"-type function?  I removed
               #       Rdo_locate_predefined_section() since it simply checks the tags as
               #       below but probably should vreate a ne function based on a "locate"
    indx <- Rdo_which_tag_eq(rdo, "\\arguments")

    if(length(indx) == 0)
        return(character(0))    # no "arguments" section in rdo; initially was: NA_character_

    txt <- Rdo_get_item_labels(rdo[[indx]])

    if (length(txt)==0)               # this chunk: copied from tools:::.Rd_get_argument_names
        return(character())
    txt <- unlist(strsplit(txt, ", *"))
    txt <- gsub("\\\\l?dots", "...", txt)
    txt <- gsub("\\\\_", "_", txt)
    txt <- .strip_whitespace(txt)


    txt <- unique(txt) # cater for duplication due to #ifdef unix or windows;
                       # see e.g. system.Rd or system2.Rd in base package
                       #
                       # todo: more thought on this?
    txt
}

                                        # 2013-10-23 new
                                        # pattern = ".*[.]Rd$" would be too dangerous.
                                        # but may be exclude could exclude foo-package
Rdreplace_section <- function(text, sec, pattern, path = "./man", exclude = NULL, ...){
    rdnames <- dir(path, pattern, full.names=TRUE)
                                           # todo: allow 'exclude' to have length more than 1
    if(!is.null(exclude)){
        rdexcl <- dir(path, exclude, full.names=TRUE)
        rdnames <- setdiff(rdnames, rdexcl)
    }

    for(nam in rdnames){
        try(Rdo_set_section(text, sec, nam))
    }

    rdnames # todo: return a logical vector of success/failures?
}

Rdo_locate_core_section <- function(rdo, sec){                       # 2013-12-08
    # Rdo_which_tag_eq(rdo, sec)
    secall <- Rdo_sections(rdo)
    sec.names <- sapply(secall, function(x) x$title)
    indx <- which(sec.names == sec)

    if(length(indx)==0)
        list()
    else
        secall[[indx]]
}
