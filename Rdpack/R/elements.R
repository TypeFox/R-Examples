Rdo_newline <- function(){ structure("\n", Rd_tag = "TEXT") }

Rdo_is_newline <- function(rdo){
    is.character(rdo) && length(rdo) == 1 && rdo == "\n"
}

Rdo_comment <- function(x = "%%"){ structure(x, Rd_tag = "COMMENT") }

Rdo_text  <- function(x){ structure(x, Rd_tag = "TEXT") }         # todo: check for % and \\ ?
Rdo_verb  <- function(x){ structure(x, Rd_tag = "VERB") }
Rdo_Rcode <- function(x){ structure(x, Rd_tag = "RCODE") }

Rdo_tag    <- function(x, name){    structure(x, Rd_tag = name) }

Rdo_macro <- function(x, name){
    if(!is.character(x))                           # the argument has been prepared in advance
        return(Rdo_tag(x,name))                    # (bar the Rd_tag)

    cnt <- Rdo_text(x)
    if(name %in% c("\\eqn", "\\deqn"))   # these are like two argument macros.
        cnt <- list(cnt)

    Rdo_macro1(cnt, name)
}

.Rdo_elem <- function(x){
    if(is.null(attr(x, "Rd_tag"))){
        if(is.character(x))
            list(Rdo_text(x))
        else
            x
    }else
        list(x)
}

Rdo_macro1 <- function(x, name){    structure(list(x), Rd_tag = name) }
Rdo_macro2 <- function(x, y, name){
    structure(list(.Rdo_elem(x), .Rdo_elem(y)), Rd_tag = name)
}

Rdo_item <- function(x,y){ Rdo_macro2(x, y, "\\item") }

Rdo_sigitem <- function(x,y, newline=TRUE){                 # todo: newline not currently used
    if(missing(y))                   # default as in promptMethods()
        y <- list(  Rdo_newline()
                  , Rdo_comment("%%  ~~describe this method here~~")
                  , Rdo_newline() )

    Rdo_item(  list(Rdo_macro1(Rdo_Rcode(c(x)), "\\code"))   # c(x) to remove attributes.
             , list(y, Rdo_newline()) )                      # (todo: does it really do it?)
}
                                                            # as.vector() to remove name attr.
Rdo_piecetag <- function(name)    as.vector(Rdo_piece_types[name])
Rdo_sectype <- function(x)  as.vector( Rdo_predefined_sections[x] )

is_Rdsecname <- function(name)
    name %in% names(Rdo_predefined_sections) # 2012-10-08 todo: wished to change it to
                                             # .rd_sections but those are with backslashes.
                                             # need to consolidate.

.char2Rdtag <- function(content, name = NULL){
    if(is.null(name) || name == "")
        return(Rdo_text(content))

    tag <- Rdo_piecetag(name)
    if(is.na(tag))
        tag <- "TEXT"    # todo: check the content and mark it RawText if it contains % or \\

    Rdo_tag(content, tag)
}

char2Rdpiece <- function(content, name, force.sec = FALSE){
    content <- .char2Rdtag(content, name = name)

    if(!force.sec || is_Rdsecname(name))
        Rdo_macro1(content, paste("\\", name, sep=""))
    else
        Rdo_macro2(Rdo_text(name), content, "\\section")  #todo: dali name e "TEXT" tuk
}

list_Rd <- function(..., Rd_tag=NULL, Rd_class = FALSE){
    dots <- list(...)
    nams <- allNames(dots)

    mflags <-  grepl("^\\\\", nams)                                             # LaTeX macros
    wrk <- mapply(function(x, i) if(mflags[i]) Rdo_macro(x, nams[i])
                                 else          x
                  , dots, seq_along(dots), SIMPLIFY = FALSE, USE.NAMES=FALSE)

                        # process other names - title, name, etc., not starting with backslash
    sflags <-  grepl("^([^\\\\]|$)", nams)    # other named  elements, including ""

    wrk <- mapply(function(x, i)
                      if(sflags[i]){
                          if(is.character(x)) #  Rdo_piece(x, nams[i])
                              char2Rdpiece(x, nams[i])
                          else           # structure(x, Rd_tag = paste("\\", nams[i], sep=""))
                              Rdo_tag(x, paste("\\", nams[i], sep=""))
                      }else
                          x
                  , wrk, seq_along(dots), SIMPLIFY = FALSE, USE.NAMES=FALSE)

    res <- wrk
    if(!is.null(Rd_tag))
        attr(res, "Rd_tag") <- Rd_tag

    if(Rd_class)
        class(res) <- "Rd"

    res
}

c_Rd <- function(...){
    dots <- list(...)
    rdtags <- sapply(dots, function(x) attr(x,"Rd_tag") )

    if(!is.null(names(rdtags)))
        names(rdtags) <- NULL

    for(i in seq_along(dots)){          # process "character" elements
        elem <- dots[[i]]
        if(!is.character(elem))
            next

        if(is.null(rdtags[[i]])){             # todo: convert the string to Rd object (parse
            attr(elem, "Rd_tag") <- "TEXT"    # it); for now, simply attach a TEXT Rd_tag (not
        }                                     # necessarilly correct!)  and embed in a list.
        dots[[i]] <- list(elem)
    }

                                    # now everything should be a list; # tags may have changed
    res <- do.call("c", dots)
    if(any(sapply(dots, function(x) inherits(x, "Rd"))))
        class(res) <- "Rd"                              # if any arg. is "Rd" res is also "Rd"
    else{                                               # set Rd_tag to the 1st non-null rdtag
        indx <- which(sapply(dots, function(x) !is.null(attr(x,"Rd_tag"))))
        if(length(indx)>0)
            attr(res, "Rd_tag") <- rdtags[[ indx[1] ]]
    }
    res
}

