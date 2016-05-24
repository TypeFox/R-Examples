                                          # backslashes need to be escaped before braces
                                          # otherwise backslahes in \{ and \} may get escaped.
.bs2 <- function(x) gsub("([\\]+)", "\\1\\1", x)                       # escape backslash (bs)
.esc_br <- function(x) gsub("([{}])", "\\\\\\1", x)                    # escape {, }

.bspercent <- function(x) gsub('([^\\%]?([\\][\\])*)%', '\\1\\\\%', x) # escape percents
.anypercent <- function(x){
    tag <- attr(x, "Rd_tag")
    if(is.null(tag)  ||  tag != "COMMENT" )
        .bspercent(x)   # expects, correctly that `tag' does not disappear
    else                # but see the comment above about usage with Rdapply.
        x
}

           # this seems incomplete, since \v and \l should be doubled only in R strings when
           #         in RCODE but in contexts where parse_Rd considers them (and other
           #         escaped sequences) to be markup macros, they will be in the Rd_tag
           #         attribute, not in the string. So, it seems that more complicated
           #         processing is not needed.

           # 2012-10-14 dobavyam obrabotka na poveche ot edna cherta. Togava vsichki osven
           # poslednata tryabva da se udvoyat (may), a poslednata se udvoyava ako e
           # posledvana ot v, l ili kray na string.
.escape_bs_in_Rcode <- function(rdo){
    f <- function(x) if(grepl("\\\\", x)){  # if x contains any backslashes
                                     # before 2012-10-14: gsub("(\\\\+)(v|l)", "\\1\\1\\2", x)
                         res <- x
                         res <- gsub("(\\\\+)(\\\\)", "\\1\\1\\2", res)
                         res <- gsub("(\\\\)(v|l)", "\\1\\1\\2", res)
                         res <- gsub("(\\\\)$", "\\1\\1", res)
                         res <- gsub("(\\\\)('|`|\")", "\\1\\1\\2", res)
                         # browser()
                         res
                     }else
                         x

    Rdtagapply(rdo, f, "RCODE")
}
                                                    # 2012-09-29 included  srcref  processing
                                   # todo: this function was patched many times, needs rewrite
Rdo2Rdf <- function(rdo, deparse = FALSE, ex_restore = FALSE, file = NULL, rcode = TRUE,
                    srcfile = NULL){
    if(is.character(rdo))                # otherwise print(rdo) may go into infinite recursion
        rdo <- list(rdo)

    if(class(rdo) != "Rd")                  # to force using print.Rd()
        class(rdo) <- "Rd"

    if(is.character(srcfile)){          # remember which sections have not changed
        rdoorig <- parse_Rd(srcfile)
        unchanged_sec <- .rdo_equal_sections(rdo, rdoorig)
    }

    # browser()

    if(rcode){
        rdo <- .escape_bs_in_Rcode(rdo)  # this also does the examples
    }else if(ex_restore){    # 2012-09-27 promenyam da izpolzva .escape_bs_in_Rcode
                   # There should be no more than one `examples' section in a proper Rd object
                   # and no `NULL' Rd_tag's at the top level.  Allow them here since `rdo' may
                   # be an `Rd' piece, not a whole Rd object or fragment.  In this way this
                   # function can be used for intermediate transformations.  Note though that
                   # print.Rd used below to produce the text output may be more picky.
        indx <- Rdo_which_tag_eq(rdo, "\\examples")
        for(i in seq_along(indx)){
            if(!is.null(indx[[i]]))
                rdo[[ indx[[i]] ]] <- .escape_bs_in_Rcode(rdo[[ indx[[i]] ]])
        }
    }

    rdo <- Rdtagapply(rdo, .esc_br, "nonmathVERB")                             # escape \{, \}
    rdo <- Rdtagapply(rdo, .anypercent, "nonmath")                             # escape %

    # 2012-10-14
    # rdo<-Rdtagapply(rdo, function(x)  gsub("((^|[^\\])([\\][\\])*)[\\]$", "\\1\\\\\\\\", x),
    #                   "VERB")
    rdo <- Rdtagapply(rdo, function(x)  gsub("(^|[^\\])(\\\\+)$", "\\1\\2\\2", x), "VERB")

                                          # pos_only = function(x){ res <- 1; browser(); res }
    pos_filecmd <- Rdo_locate(rdo, function(x) .tag_eq(x,"\\file"), lists = TRUE)
    for(pos in pos_filecmd)
        rdo[[pos]] <- Rdtagapply(rdo[[pos]], function(x) gsub("(\\\\)", "\\1\\1", x), "TEXT")

    # pos_filecmd <- Rdo_locate(rdo, function(x) .tag_eq(x,"\\samp"), lists = TRUE)
    # for(pos in pos_filecmd)
    #     rdo[[pos]] <- Rdtagapply(rdo[[pos]],
    #                              function(x) gsub("(\\\\+)([^\\])", "\\1\\1\\2", x),
    #                              "VERB")

       # krapka, za nesta kato \code{\\} (see in base R: basename.Rd, body.Rd, Arithmetic.Rd)
       #   2012-10-14 promenyam tozi (cautious) variant, koyto samo udvoyava edna cherta
       #   kogato sa necheten broy s variant koyto udvoyava vsichki v kraya na \code{} string.
       #   todo: da ne go pravi v examples section (no tam to trudno ste se sluchi)
       #            rdo <- Rdtagapply(rdo, function(x) gsub("((^|[^\\])([\\][\\])*)[\\]$",
       #             "\\1\\\\\\\\", x),  "RCODE")
       # 2012-10-14 otkomentiram, promenich v escape_bs po-gore.
    # rdo <- Rdtagapply(rdo, function(x)  gsub("([\\]+)$", "\\1\\1", x), "RCODE")

    if(is.character(srcfile)){            # replace unchanged sections with dummy contents
        unchanged_titles <- sapply(unchanged_sec, function(x) x$title)
        rdo <- .replace_dummy(rdo, unchanged_titles)
    }
                           # this was used for saving before introducing tfn, etc.
                           # res <- capture.output(print(rdo, deparse = deparse), file = file)

                           # as.character.Rd ima bug svarzana s newcommand i/ili argumentite
                           # na href, vzh. naprimer which.min.Rd v base, tam as.character dava
                           #     \href{{http://CRAN.R-project.org/package=nnet}{\pkg{nnet}}}
                           # (nay-vanshnite skobi {} sa izlishni). todo: tova e krapka!
    rdotxt <- paste0(as.character(rdo, deparse = deparse), collapse="")
    rdotxt <- gsub("\\\\href\\{(\\{[^{}]*\\}\\{[^}]*\\}+)\\}", "\\\\href\\1", rdotxt)

    # todo: krapka; \code{\{} becomes \code{{} which wrecks havoc for parse_Rd; \code{\}}
    rdotxt <- gsub("(\\\\code\\{)(\\{|\\})(\\})", "\\1\\\\\\2\\3", rdotxt)

    tfn <- tempfile()             # use a temporary file in case file and srcfile are the same
    on.exit(unlink(tfn))
    if(is.character(srcfile)){
        tfn <- tempfile()
        res <- capture.output(cat(rdotxt, sep = "", collapse = ""), file = tfn)# writes to tfn

        rdocur <- parse_Rd(tfn)  # to set srcref
        srcrefpos <- .srcrefpos(rdocur, rdoorig, unchanged_sec)

        rdotxt <- rdo_text_restore(rdocur, rdoorig, srcrefpos, file=tfn)

        nc_ind <- Rdo_which_tag_in(rdoorig,  c("\\newcommand", "\\renewcommand"))
        if(length(nc_ind) > 0){
            nclines <- sapply(nc_ind, function(x) as.character(attr(rdoorig[[x]],"srcref")))
            rdotxt <- c(nclines, rdotxt) # put before anything else todo: could try to put at
        }                                # original place or at least after any comments?

        writeLines(rdotxt, tfn)  #overwrites tfn

        res <- if(is.null(file))
                   paste0(rdotxt, collapse = "\n")
               else{
                   file.copy(tfn, file, overwrite = TRUE) # todo: check success
                   NULL  # for clarity; capture.output above set it to NULL as tfn is not NULL
               }
    }else
        # 2012-10-14 res <- capture.output(cat(rdotxt, sep = "", collapse = ""),  file = file)
        res <- capture.output(cat(rdotxt, sep = "", collapse = "",  file = file))

    if(is.null(file))
        res
    else{
        cat("\tThe Rd content was written to file ", file, "\n")
        invisible(res)    # res is NULL here
    }
}

.tmp_pos <- function(name, pos_list){
    for(elem in pos_list)
        if(elem$title == name)
            return(elem$pos)
    NULL
}

.rdo_srcref <- function(rdo, tag){  # todo: special cases!
    pos <- Rdo_which_tag_eq(rdo, tag)
    attr(rdo[[ pos[1] ]], "srcref")
}

.rdo_replace_at <- function(text, pospair){
    newtext <- as.character(pospair[[2]])
    m <- length(pospair[[1]])

    beg_line <- pospair[[1]][1]
    beg_col  <- if(m > 4) pospair[[1]][5] else pospair[[1]][2]

    end_line <- pospair[[1]][3]
    end_col  <- if(m > 4) pospair[[1]][6] else pospair[[1]][4]

    res <- c(text[seq_len(beg_line - 1)],
             if(beg_col > 1) paste0(substr(text[beg_line], 1, beg_col - 1),  newtext[1])
             else            newtext[1],
             newtext[-1]  )

    le <- nchar(text[end_line]) # 2012-10-13 was: length(text[end_line])  (!?)
    if(end_col < le)
        res[length(res)] <- paste0(res[length(res)], substr(text[end_line], end_col + 1, le))

    c(res, text[-(1:end_line)])
}

         # 2012-10-13 dobavyam 'ends' po-dolu, sluchva se sections da ne zapochvat na nov red!
         #            vzh. NumericConstants.Rd v base ( stava: \note{dummy} \seealso{dummy} )
rdo_text_restore <- function(cur, orig, pos_list, file){
    res <- readLines(file)
    if(length(pos_list) == 1){
        res <- .rdo_replace_at(res, pos_list[[1]])
    }else{
        starts <- sapply(pos_list, function(x) x[[1]][[1]])
        ends <- sapply(pos_list,
                       function(x) if(length(x[[1]]) > 4) x[[1]][5] else x[[1]][2]  )
        p <- order(starts, ends, decreasing = TRUE)

        dec_pos_list <- pos_list[p]
        for(pos in dec_pos_list)
            res <- .rdo_replace_at(res, pos)
    }
    res
}

.without_duplicates <- function(x){
    x[!(x %in% unique(x[duplicated(x)]))]
}
