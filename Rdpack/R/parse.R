.parse_Rdlines <- function(lines){
    tmpfile <- tempfile("Rdlines", fileext = ".Rd")
                                 # unlist is harmless if `lines' is already a character vector
    cat(unlist(lines), file = tmpfile, sep = "\n")          # todo: catch errors
    res <- parse_Rd(tmpfile)                              # todo: catch errors
    unlink(tmpfile)
    res
}

parse_Rdtext <- function(text, section = NA){
    if(!is.na(section))
        text <- c(paste(section,"{",sep=""), text, "}")

    res <- .parse_Rdlines(text)    # write(text, sep="\n", file = file); res <- parse_Rd(file)

    rdtag <- attr(res[[1]], "Rd_tag")
    if(!is.na(section) &&  !is.null(rdtag) && rdtag==section)     # todo: tova e krapka!
        res <- res[[1]]
    if(identical( c(res[[1]]), "\n" ))   # c( ) to strip attributes
       res[[1]] <- NULL

    res
}

parse_Rdpiece <- function(x, result=""){
    maket <- list_Rd(title="Dummy title"
                     , name="dummyname"
                     , "\\description"="Dummy description"
                     , note = x
                     , Rd_class=TRUE
                     )

    fn <- tempfile(pattern = "Rdpiece", fileext = "Rd")
    cat(as.character(maket), file=fn, sep="")    # todo: error processing
    wrk <- parse_Rd(fn)

    unlink(fn)

    if(result=="text"){                                             # Rd_help2txt is from gbRd
        res <- Rd_help2txt(wrk, keep_section="\\note", omit_sec_header=TRUE)
    }else{
        indx <- Rdo_which_tag_eq(wrk, "\\note")
        res <- wrk[[indx]]
        if(!is.null(attr(res,"Rd_tag")))
            attr(res,"Rd_tag") <- NULL
    }

    res
}
