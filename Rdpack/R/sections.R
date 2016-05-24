.sec_name <- function(x){           # this assumes that x is a section Rd piece.
    a <- attr(x, "Rd_tag")  #TODO: paste tuk e prapka, dobavi treti element!
    if(a != "\\section")  a  else  paste(as.character(x[[1]]), collapse="") # 'as' is to drop attributes
}

.sec_content <- function(x){           # this assumes that x is a section Rd piece.
    a <- attr(x, "Rd_tag")
    if(a != "\\section")  x  else  x[[2]]
}

.order_sections <- function(rdo){     # todo: take care of empty lines at top level.
    old <- Rdo_tags(rdo)

    n <- length(rdo)
    newpos <- numeric(n)
    curpos <- 1
    for(i in 1:length(.rd_sections)){
        ind <- which(old == .rd_sections[i])
        len <- length(ind)
        if(len > 0){
            newpos[curpos:(curpos+len-1)] <- ind
            curpos <- curpos + len
        }
    }
    if(curpos < n){  # note everything in "old" accounted for. Can there be things
                               # other than newlines?
        missedpos <- seq_len(n)[ !(seq_len(n) %in% newpos) ]

        # browser()   # for testing

        newpos[curpos:n] <- missedpos   # to make sure that nothing is omitted
    }

    if(identical(newpos, seq_along(newpos))){
        rdo  # nothing to reshuffle here
    }else{
        res <- rdo[newpos]
        attributes(res) <- attributes(rdo)
        res
    }
}

.replace_dummy <- function(rdo, unchanged_titles){
    secall <- Rdo_sections(rdo)
    for(sec in secall)
        if(sec$title %in% unchanged_titles)
            rdo[[ sec$pos ]] <- char2Rdpiece("dummy", sub("^\\\\", "", sec$title),
                                             force.sec = TRUE)
    rdo
}

.srcrefpos <- function(rdo, rdoorig, pos_list){
    secall <- Rdo_sections(rdo)
    res <- vector(mode="list", length(pos_list))
    ind <- 0
    for(sec in secall){
        pos <- .tmp_pos(sec$title, pos_list)
        if(!is.null(pos)){
            ind <- ind + 1
            res[[ind]] <- list(attr(rdo[[sec$pos]], "srcref"),
                               attr(rdoorig[[pos]], "srcref"))
        }
    }
    res
}

.replace_from <- function(rdo, rdoorig, pos_list){
    secall <- Rdo_sections(rdo)
    for(sec in secall){
        pos <- .tmp_pos(sec$title, pos_list)
        if(!is.null(pos))
            rdo[[sec$pos]] <- rdoorig[[pos]]
    }
    rdo
}

Rdo_sections <- function(rdo){
    res <- Rdo_locate(rdo, function(x) .tag_in(x, .rd_sections), .sec_name,
                      lists = TRUE, nested = FALSE)          # 2012-10-20 new arg. nested used
    lapply(res, function(x){ names(x)[names(x) == "value"] <- "title"; x})
}

Rdo_modify_simple <- function(rdo, text, section, ...){
    val <- char2Rdpiece(text, section)
    Rdo_modify(rdo, val, ...)
}

Rdo_empty_sections <- function(rdo, with_bs = FALSE){           # rdo is Rd object or filename
    if(is.character(rdo))
        rdo <- parse_Rd(rdo)

    pat <- ".*Dropping empty section[ ]+[\\](.*)$"
    wrk <- checkRd(rdo)
    ind <- grepl(pat, wrk)
    res <- gsub(".*Dropping empty section[ ]+[\\]", "", wrk[ind])
    if(with_bs)
        res <- sapply(res, function(x) paste("\\",x,sep="")) # prepend backslash
    res
}

Rdo_drop_empty <- function(rdo, sec = TRUE){ # todo: drop other empty things?
    emptysec <- Rdo_empty_sections(rdo, with_bs = TRUE)

    rdotags <- Rdo_tags(rdo)

    stopifnot(length(rdotags) == length(rdo)) # paranoic; but tools:::RdTags above is
                                              # non-exported function
    res <- rdo[ !(rdotags %in% emptysec) ]
    class(res) <- class(rdo)

    res
}
                                                       # rd is a character vector of filenames
Rdo_replace_section <- function(rdo, val, create = FALSE, replace = TRUE){
    Rdo_modify(rdo, val, create=create, replace = replace)
}

Rdo_set_section <- function(text, sec, file, ...){
    piece <- char2Rdpiece(text, sec)
    rdo <- Rdo_replace_section(parse_Rd(file), piece, ...)
    Rdo2Rdf(rdo, file = file)
}
                                     # todo: this does not consider the possibility for #ifdef
                          # return the location (full index) of the content of a section.  for
                          # standard sections this is at the top level but for user defined
                          # sections it is deeper.  see also tools:::.Rd_get_section
.locate_sec_content <- function(rdo, sec,  onlyone = TRUE){  # 2011-12-12 add 'more' arg.
    tags <- Rdo_tags(rdo)
    indx_sec <- which(tags == sec) # in case this is a standard section, todo: do it properly!

    if(length(indx_sec) > 0){  # standard section
        if(length(indx_sec) > 1  &&  onlyone){
            warning("More than one section with this name, using the first one.")
            indx_sec <- indx_sec[[1]]
        }
        res <- indx_sec        # for standard sections the content is simply that tagged list.
    }else{  # user section;  todo: arg. 'onlyone' silently ignored here
        indx_sec <- which(tags == "\\section")
        if(length(indx_sec) > 0 ){             # wraps in c() here to get rid of attributes;
                                               #  note also the two 1's in the index (seems
                                               #     that each arg of "\\section" is a list).

            methindx <- sapply(indx_sec, function(x) identical(c(rdo[[c(x,1,1)]]), sec) )
            indx_sec <- indx_sec[methindx] # take only those that match.
        }
        if(length(indx_sec) == 0)
            res <- NA # stop("No section `??'. Todo: create it instead of throwing an error.")
        else{
            if(length(indx_sec)>1){
                warning("More than 1 sections `Methods' found, using the first one.")
                indx_sec <- indx_sec[1]
            }
            res <- c(indx_sec,2)   # for user sections the content is deeper.
        }
    }
    res
}

.rdo_equal_sections <- function(rdo1, rdo2, sections = NULL){ # todo: arg. sections ne se
                                                              #       izpolzva veche!
    wrk1 <- Rdo_remove_srcref(rdo1)
    wrk2 <- Rdo_remove_srcref(rdo2)

    sec1 <- Rdo_sections(wrk1)
    sec2 <- Rdo_sections(wrk2)

    nams1 <- sapply(sec1, function(x) x$title)
    nams2 <- sapply(sec2, function(x) x$title)

    u1nam <- names(which(table(nams1) ==1))   # u1nam <- .without_duplicates(nams1)
    u2nam <- names(which(table(nams2) ==1))   # u2nam <- .without_duplicates(nams2)

    unam <- intersect(u1nam, u2nam)

    i1 <- which(nams1 %in% unam)
                                    # 2012-10-22 dobavyam all.equal, for some reason sometimes
                                    #         seeemingly identical things are not "identical",
                                    #         e.g. Special.Rd dava takova nesto za details.
                                    # todo: tova e krapka razibirase.
    isame <- sapply(i1, function(x){
                            i2 <- which(nams2 == nams1[x]);
                            if(length(i2)==1 &&
                               ( identical(wrk1[[sec1[[x ]]$pos]], wrk2[[sec2[[i2]]$pos]]) ||
                                 isTRUE(all.equal(wrk1[[sec1[[x ]]$pos]], wrk2[[sec2[[i2]]$pos]]))
                               ))
                                i2
                            else
                                NA
                        })
    # browser()

    isame <- isame[!is.na(isame)]
    res <- sec2[isame]

    res
}

# tools:::.Rd_get_name
# tools:::.Rd_deparse
#
# . function (x, sec)    # get the Rd
# {
#     x <- .Rd_get_section(x, "section")
#     if (length(x))
#         .strip_whitespace(.Rd_deparse(x, tag = FALSE))
#     else character()
# }


# list_files_with_exts("./man","Rd")
# list_files_with_type("./man","docs")
# Rdo_replace_section(a, char2Rdpiece("documentation", "keyword"), replace = "Rd")

# files <- list_files_with_exts(dir, "Rd")

