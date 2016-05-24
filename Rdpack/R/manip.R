                                                          # generic rdo manipulation functions
append_to_Rd_list <- function(rdo, x, pos){  # pos is a numeric vector
    rdo[[pos]] <- structure( c(rdo[[pos]], x), Rd_tag = attr(rdo[[pos]],"Rd_tag") )
    rdo
}
                                   # todo: allow vector `pos' to insert deeper into the object
               # todo: character `pos' to insert at a position specified by "tag" for example?
Rdo_insert_element <- function(rdo, val, pos){
    n <- length(rdo)
    if(missing(pos))
        pos <- n+1
    stopifnot(1 <= pos, pos <= n+1)

    res <- if(pos==1)        c(list(val), rdo)
           else if(pos==n+1) c(rdo, list(val))
           else              c( rdo[1:(pos-1)], list(val), rdo[pos:n])
    attributes(res) <- attributes(rdo)             # todo: more guarded copying of attributes?
    res
}

Rdo_get_insert_pos <- function(rdo, tag){
    tags <- Rdo_tags(rdo)
    pos <- which(tags == tag)

    pos <- if(length(pos)>0)
               max(pos) + 1     # put after the last occurence of the same section
           else if(tag %in% .rd_sections){
               indx <- which(.rd_sections == tag)
               iafter <- which( tags %in% .rd_sections[-seq_len(indx)] )
               iafter[1]                    # put in front of the next section in .rd_sections
           }else{                                            # otherwise put at the end of rdo
               length(rdo)+1
           }
    pos
}


Rdo_insert <- function(rdo, val, newline = TRUE){
    tag <- attr(val,"Rd_tag")
    if(is.null(tag))
        stop("val has no Rd_tag attribute.")

    pos <- Rdo_get_insert_pos(rdo, tag)

                                                       # 2013-03-29 - additional check for pos
                                                       #        todo: this is incomplete!
    if(newline && (   (length(pos) == 1 && pos <= length(rdo) && !Rdo_is_newline(rdo[[pos]]))
                   || (length(pos) == 1 && pos >  length(rdo))
                   || length(pos) > 1   # todo: this is incomplete
                   )
       )
        rdo <- Rdo_insert_element(rdo, Rdo_newline(), pos)

    rdo <- Rdo_insert_element(rdo, val, pos)
    if(newline   &&  pos > 1  &&  !Rdo_is_newline(rdo[[pos-1]]))
        rdo <- Rdo_insert_element(rdo, Rdo_newline(), pos)

    rdo
}
                             # this is not specific to "sections", any Rd_tag in 'val' will do
Rdo_modify <- function(rdo, val, create=FALSE, replace=FALSE, top = TRUE){
    sec <- attr(val,"Rd_tag")
    if(is.null(sec))
        stop("val has no Rd_tag attribute.")

    rdotag <- attr(rdo,"Rd_tag")
    if(top && !is.null(rdotag)  && rdotag == sec){                 # rdo is the section `sec';
        if(replace)   # todo: character 'replace' not handled yet.         # append to the end
            return(val)
        else{
            rdo <- append_to_Rd_list( list(rdo), val, 1)
            return( rdo[[1]] )
        }
    }

    pos <- Rdo_which_tag_eq(rdo, sec)

    if(length(pos)>0){
        if(is.character(replace)){    # the value to replace
           goodpos <- sapply(pos,                        # todo: the check below is simplistic
                             function(x){
                                 locval <- rdo[[x]]
                                 is.list(locval) && length(locval) == 1  &&
                                                    locval[[1]] == replace
                             })
           goodpos <- pos[goodpos]
           if(length(goodpos) > 0)
               rdo[[ goodpos[1] ]] <-  val
           else if(create)
               rdo <- Rdo_insert(rdo, val)

       }else if(replace)
           rdo[[ pos[1] ]] <- val                     # pos[1] in case there are more
       else
           rdo <- append_to_Rd_list(rdo, val, pos[1])
    }else if(create)                                # argument section not present, create one
        rdo <- Rdo_insert(rdo, val)    # rdo[[ length(rdo)+1 ]] <- val
    else
        stop("No section ", sec, " to replace or modify.")
    rdo
}


## The 'locate_' family of functions finds (recursively) positions of elements of rdo objects.
## A position is a vector suitable for use in '[[' or a list of such vectors.

          # find the index the element enclosing rdo[[ind]] that has Rd_tag 'tag'.
          # (e.g. if rdo[[ind]] is an '\item', we may wish to grab the enclosing '\describe'.)
.locate_enclosing_tag <- function(rdo, ind, tag, baseind = numeric(0), mindepth=0){
    repeat{
        fullind <- c(baseind,ind)
        rdtag <- if(length(fullind) > 0) attr(rdo[[fullind]], "Rd_tag")
                 else                    attr(rdo           , "Rd_tag")

        if(is.character(rdtag) && rdtag == tag)
            return(ind)
        else if(length(ind)==mindepth)
            return(as.numeric(NA))

        ind <- head(ind,-1) # drop the last level
    }
}

## The '_top_' family of functions look only at the top level of rdo, in the sense that they
## examine the elements of rdo at the top level and do not search recursively its structure.
## Some of these functions examine also rdo as a whole (todo: maybe all of them should do
## this for consistency).
.locate_top_tag <- function(rdo, sec, tag = "\\describe"){
    if(sec == ""){     # search at top level of rdo
        indx1 <- integer(0)
    }else{
        indx1 <- .locate_sec_content(rdo, sec)
        if(sec %in% c("\\arguments", "\\value")) # these sections are themselves lists
            return(indx1)
        rdo <- rdo[[indx1]]
    }

    indx3 <- Rdo_which_tag_eq(rdo, tag)
    if(length(indx3)==0)
        return(integer(0))

    indx3 <- indx3[[1]]  # playing safe, in case more than one
    c(indx1,indx3)                             # todo: option to return a list of all matches?
}

.locate_top_items <- function(rdo, sec, labels=FALSE){
    indx1 <- .locate_top_tag(rdo, sec)
    if(length(indx1) == 0)
        return(list())

    indx2 <- Rdo_which_tag_eq(rdo[[indx1]], "\\item")

    res <- lapply( indx2, function(x) c(indx1, x))               # returns a list of positions

             # todo: item labels (i.e., the 1st arg. of \item) may be long, is this a problem?
    if(isTRUE(labels) || is.character(labels))          # set names attribute if requested
        names(res) <- sapply(res, function(x) .get_item_label(rdo,x) )

    if(is.character(labels))   # take only the requested items
        res <- res[labels]
    res
}

.nl_and_indent <- function(rdo, newlines = TRUE, indent = "    "){
    nl <- Rdo_newline()
    tab <- indent # todo: allow it to be dropped
    wrk <- lapply(rdo, function(x) list(nl, tab, x) )
    res <- do.call("c", wrk)       # todo: is this the correct arrangement - ne, gubi Rd_tag
    attr(res, "Rd_tag") <- attr(rdo, "Rd_tag")  # 2012-10-16 todo: "Rd" attribute? - use c_Rd?
    res
}

.get_subset <- function(rdo, pos, rdtag = FALSE, keep.class = FALSE, newlines = FALSE){
    if(newlines){
        nl <- Rdo_newline()
        wrk <- lapply(pos, function(x) list(nl, rdo[[x]]) )
        res <- do.call("c", wrk)       # todo: check if this is the correct arrangement
    }else
        res <- lapply(pos, function(x) rdo[[x]])
    if(isTRUE(rdtag))
        attr(res, "Rd_tag") <- attr(rdo, "Rd_tag")
    else if(is.character(rdtag))
        attr(res, "Rd_tag") <- rdtag

    if(keep.class)
        class(res) <- class(rdo)

    res
}

.get_item_label <- function(rdo, pos){
                 # 2013-12-01 no need of looking for leaves after introducing .ascharRd()
                 #         indx <- c(pos, .locate_item_label(rdo, pos))
                 #         res <- rdo[[indx]]
    indx <- c(pos, 1)
    res <- .ascharRd(rdo[[indx]])

    res
}

                                       # 'labels' is a logical or a character vector of labels
.get_top_items <- function(rdo, sec, labels){               # e.g. sec = "Slots", "arguments"
    itemind <- .locate_top_items(rdo, sec, labels=labels)
    lapply(itemind, function(x) rdo[[ x ]] )           # see also .get_subset()
}

.get_top_labels <- function(rdo, sec){
    indx <- .locate_top_items(rdo, sec)
    sapply(indx, function(x) .get_item_label(rdo,x))
}
                                                           # manipulation of specific sections

                        # todo: this is very similar to parse_Rdname
                        #       2012-11-04 partly done, introduced function .parse_long_name()
.get.name_content <- function(rdo){      # todo: error processing             Rd section \name
    name <- rdo[[ Rdo_which_tag_eq(rdo, "\\name") ]]
                                   # 2012-11-04   # dropping "-methods" or similar, if present
                                                  # short <- gsub("^([^-]+)-.*", "\\1", name)
    short <- .parse_long_name(name)["name"]
    c(name = name, short = short)
}
                                                                       # Rd section \arguments
Rdo_append_argument <- function(rdo, argname, description = NA, indent = "  ", create=FALSE){
    if(is.na(description))
        description <- "~~ TODO: describe this argument. ~~"

    if(length(description) == 1  && length(argname) > 1)
        description <- rep(description, length(argname))

    wrk1 <- lapply(1:length(argname),
                   function(i)
                      list(indent, Rdo_item(argname[i], description[i]), Rdo_newline())
                   )

    wrk <- Rdo_macro( do.call("c", wrk1), "\\arguments")
    Rdo_modify(rdo, wrk, create=TRUE)
}
                                                                           # Rd section \usage

     # todo: get_usage_text can be generalised to any Rd section but it is better to use a
     #       different approach since print.Rd() does not take care for some details (escaping
     #       %, for example). Also, the functions that use this one assume that it returns R
     #       code which may not be the case if the usage section contains Rd comments.

                                        #          get_usage_text("../gbRd/man/Rd_title.Rd")
get_usage_text <- function(rdo){        #    ut <- get_usage_text("../gbRd/man/Rd_fun.Rd")
    if(is.character(rdo) && length(rdo)==1)
        rdo <- parse_Rd(rdo)
                                                # Alternatives:
    pos <- .locate_sec_content(rdo, "\\usage")  # wrk <- Rdo_section(rdo, "\\usage")
                                                # wrk <- tools:::.Rd_get_section(rdo, "usage")
                    # todo: maybe the following 3 lines would be more generally useful?

                                        #2012-10-04 drop comments and specials
    rdo_u <- toolsdotdotdot.Rd_drop_nodes_with_tags(rdo[[pos]], c("COMMENT","\\special"))

    wrk <- structure(rdo_u, Rd_tag = "Rd")
    class(wrk) <- "Rd"
    paste(capture.output(print(wrk)), collapse="\n")
}

Rdo_show <- function(rdo){                                   # 2012-09-22 renamed, was Rd_show
    outfile <- tempfile(fileext = ".txt")
    file.show(Rd2txt(rdo, outfile))
    unlink(outfile)
    invisible(NULL)
}

Rdo_reparse <- function(rdo){
    outfile <- tempfile(fileext = "Rd")
    Rdo2Rdf(rdo, file = outfile)
    rdo <- parse_Rd(outfile)
    unlink(outfile)
    rdo
}

Rdo_collect_aliases <- function(rdo){
    pos <- Rdo_locate(rdo, function(x) .tag_eq(x,"\\alias"), lists=TRUE)
    res <- sapply(pos, function(x) as.character(rdo[[ x ]]))
    nams <- sapply(pos, function(x){
                           if( length(x) == 1) ""
                           else    # 2012-10-13 as.character(rdo[[ c(x[1:(length(x)-2)],1) ]])
                               as.character(rdo[[ c(x[seq_len(length(x)-2)],1) ]])
                        })
    names(res) <- gsub("[ \n]" , "", nams)

    res
}

                                               # 2013-12-02 derived from Rdo_collect_aliases()
Rdo_collect_metadata <- function(rdo, sec){   #sec is "alias", "keyword" or similar simple one
    secname <- paste0("\\", sec)
    pos <- Rdo_locate(rdo, function(x) .tag_eq(x, secname), lists=TRUE)
    res <- sapply(pos, function(x) as.character(rdo[[ x ]]))
    nams <- sapply(pos, function(x){
                           if( length(x) == 1) ""
                           else    # 2012-10-13 as.character(rdo[[ c(x[1:(length(x)-2)],1) ]])
                               as.character(rdo[[ c(x[seq_len(length(x)-2)],1) ]])
                        })
    names(res) <- gsub("[ \n]" , "", nams)

    res
}

