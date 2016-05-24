# TODO: krapka!
.patch_latex <- function(txt){   # print(bibentry,"latex") inserts \bsl macros.
    gsub("\\bsl{}", "", txt, fixed=TRUE)
}

bibentry_key <- function(x){                                                     # 2013-03-29
    attr(unclass(x[[1]][[1]])[[1]], "key")
}

get_bibentries <- function(..., package = NULL, bibfile = "REFERENCES.bib"){     # 2013-03-29
    fn <- if(is.null(package))
              file.path(..., bibfile)
          else
              system.file(..., bibfile, package = package, mustWork=TRUE)

    if(length(fn) > 1){
        warning("More than one file found, using the first one only.")
        fn <- fn[1]
    }
    res <- read.bib(file=fn, package=package)
    names(res) <- sapply(1:length(res), function(x) bibentry_key(res[[x]][[1]]))

    res
}


rebib <- function(infile, outfile, ...){                     # 2013-03-29
    rdo <- parse_Rd(infile)

    if(missing(outfile))
        outfile <- basename(infile)
    else if(identical(outfile, ""))  # 2013-10-23 else clause is new
        outfile <- infile

    rdo <- inspect_Rdbib(rdo, ...)

    Rdo2Rdf(rdo, file=outfile, srcfile=infile)

    rdo
}


inspect_Rdbib <- function(rdo, force = FALSE, ...){               # 2013-03-29
                   # 2013-12-08 was: pos <- Rdo_locate_predefined_section(rdo, "\\references")
    pos <- Rdo_which_tag_eq(rdo, "\\references")

    if(length(pos) > 1)
        stop(paste("Found", length(pos), "sections `references'.\n",
                   "There should be only one."
                   ))
    else if(length(pos) == 0)  # no section "refeences".
        return(rdo)

    bibs <- get_bibentries(...)

    fkey <- function(x){
                 m <- gregexpr("[ ]+", x)
                 rm <- regmatches(x, m, invert = TRUE)[[1]]
                 if(length(rm) >= 2 && rm[2] != "bibentry:")
                     rm[2]   # e.g. bibentry:all
                 else if(length(rm) < 3)     # % bibentry: xxx_key_xxx
                     ""   # NA_character_
                 else
                     rm[3]
             }

    fbib <- function(x) grepl("[ ]+bibentry:", x)
    posbibs <-  Rdo_locate(rdo[[pos]], f = fbib, pos_only = fkey)
    poskeys <- sapply(posbibs, function(x) x$value)

    print(posbibs)

    fendkey <- function(x){
                 m <- gregexpr("[ ]+", x)
                 rm <- regmatches(x, m, invert = TRUE)[[1]]
                 if(length(rm) >= 2 && rm[2] != "end:bibentry:")
                     rm[2]   # e.g. end:bibentry:all
                 else if(length(rm) < 3)     # % end:bibentry: xxx_key_xxx
                     ""   # NA_character_
                 else
                     rm[3]
             }

    fendbib <- function(x) grepl("end:bibentry:", x)
    posendbibs <-  Rdo_locate(rdo[[pos]], f = fendbib, pos_only = fendkey)
    posendkeys <- sapply(posendbibs, function(x) x$value)

    toomit <- which(poskeys %in% posendkeys)  # note: en@bibkeys:all is different! todo:
    if(length(toomit) > 0  && !force){
        poskeys <- poskeys[-toomit]
        posbibs <- posbibs[-toomit]
    }

    if(length(poskeys)==0)
        "nothing to do."
    else if(any(poskeys == "bibentry:all")){
        poskey <- posbibs[[ which(poskeys == "bibentry:all") ]]$pos

        bibstxt <- capture.output(print(bibs, "latex"))

        bibstxt <- .patch_latex(bibstxt)  # TODO: krapka!

        bibstxt <- paste(c("", bibstxt), "\n", sep="")
        endbibline <- Rdo_comment("% end:bibentry:all")

        keyflag <- "end:bibentry:all" %in% posendkeys
        if(keyflag && force){              #todo: more careful!
            endposkey <- posendbibs[[ which(posendkeys == "end:bibentry:all") ]]$pos
            rdo[[pos]] <- Rdo_flatremove(rdo[[pos]], poskey+1, endposkey)
        }

        if(!keyflag || force){
            rdo[[pos]] <- Rdo_flatinsert(rdo[[pos]], list(endbibline), poskey,
                                         before = FALSE)
            rdo[[pos]] <- Rdo_flatinsert(rdo[[pos]], bibstxt, poskey,
                                         before = FALSE)
        }
    }else{
        for(i in length(poskeys):1){
            bibkey <- posbibs[[i]]$value
            poskey <- posbibs[[i]]$pos

            bibstxt <- capture.output(print(bibs[poskeys[i]],"latex"))

            bibstxt <- .patch_latex(bibstxt)  # TODO: krapka!

            bibstxt <- list( paste( c("", bibstxt), "\n", sep="") )
            endbibline <- Rdo_comment(paste("% end:bibentry: ", bibkey))

            keyflag <- bibkey %in% posendkeys
            if(keyflag && force){                                       #todo: more careful!
                endposkey <- posendbibs[[ which(posendkeys == bibkey) ]]$pos
                rdo[[pos]] <- Rdo_flatremove(rdo[[pos]], poskey+1, endposkey)
            }

            if(!keyflag || force){ # this is always TRUE here but is left for cmmong look
                                   # with "all". todo: needs consolodation
                rdo[[pos]] <- Rdo_flatinsert(rdo[[pos]], list(endbibline), poskey,
                                             before = FALSE)
                rdo[[pos]] <- Rdo_flatinsert(rdo[[pos]], bibstxt, poskey,
                                             before = FALSE)
            }
        }
    }

    rdo
}

Rdo_flatremove <- function(rdo, from, to){  # 2013-03-30 todo: more careful!
    res <- rdo[-(from:to)]
    attributes(res) <- attributes(rdo)             # todo: more guarded copying of attributes?
    res
}

                                        # todo: move to another file later
Rdo_flatinsert <- function(rdo, val, pos, before = TRUE){                        # 2013-03-29
    depth <- length(pos)
    if(depth > 1){
        rdo[[pos]] <- Recall(rdo[[ pos[-depth] ]], val, pos[-depth])
        # todo: dali zapazva attributite na rdo?
        return(rdo)
    }

    n <- length(rdo)
    if(!before)
        pos <- pos + 1

    res <- if(pos==1)        c(val, rdo)
           else if(pos==n+1) c(rdo, val)
           else              c( rdo[1:(pos-1)], val, rdo[pos:n])
    attributes(res) <- attributes(rdo)             # todo: more guarded copying of attributes?
    res
}

