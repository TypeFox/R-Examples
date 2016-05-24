                                        # See the help of Rd_db and the examples there!
                                        #  ""  also tools/R/Rd.R for other internal functions!
                                        #  ""  and RdConv2.R for utility functions.

Rdo_empty <- function(rdtag){     # Create an empty Rdo; na 2011-11-22 dobavyam argument rdtag
    if(missing(rdtag)){
        res <- list()
        class(res) <- "Rd"
        res
    }else
        structure(list(), Rd_tag = rdtag)
}

                                                                       # Create a minimal Rdo.
                       # todo: this is very basic, among other things, needs also empty lines.
Rdo_create <- function(arguments, title="Dummy title", name="dummy name"){
    res <- Rdo_empty()
    res[[1]] <- Rd_name(name)
    res[[2]] <- Rd_title(title)
    res[[3]] <- arguments
    res
}

   # 2011-10-21 - bug  fix - "Rd_tag is not set for Rd objects, their class attribute is "Rd"
Rdo_section <- function(rdo, sec){                                # 2011-10-30 - some clean up
    type <- attr( rdo, "Rd_tag")
    if(is.null(type))
        type <- ""

    if(identical(class(rdo),"Rd") || type == "Rd"){
        tags <- tools:::RdTags(rdo)
        indx <- which(tags==sec)             # todo: What to return when there are no matches?
        if(length(indx) > 1){
            warning(paste("Found more than one section named", sec,
                          ". Using only the first match.", sep=""))
            indx <- indx[1]
        }
        rdo[[which(tags==sec)]]
    }else if (type == sec)
        rdo
    else                  # assumes elements of rdo are elements of sec and wraps accordingly.
        structure(rdo, Rd_tag = sec)
}

Rdo_set_sectag <- function(s,sectag,eltag){
    attr( s, "Rd_tag") <- eltag           # using `structure' would be more elegant...
    res <- list(s)
    attr( res, "Rd_tag") <- sectag

    res
}

Rd_title <- function(s) Rdo_set_sectag(s, sectag="\\title"    , eltag="TEXT")
Rd_name <- function(s)  Rdo_set_sectag(s, sectag="\\name"     , eltag="VERB")
Rd_args <- function(s)  Rdo_set_sectag(s, sectag="\\arguments", eltag="VERB")

Rdo_get_args <- function(rd,args,...){     # tools:::RdTags(rd[[which(tags=="\\arguments")]])
    rdo <- Rd_fun(rd,                      # todo: argument ... is not used.
                  , keep_section  = "\\arguments"
                  )

    tags <- tools:::RdTags(rdo)
    rdargs <- rdo[[which(tags=="\\arguments")]]   # use of [[]] assumes only one element here
                                                  # not a problem for installed documentation.
              # to do: Ako iskam da vklyucha obrabotka na prazni redove i drug text mezhdu
              #        item-ite, obrabotkata ste se uslozhni. Tryabva i dopalnitelen argument!
    if(missing(args))
        return(rdargs)
                                        # uslozhnyavam f, za da obrabotva i sluchai, kogato
                                        # nyakolko argumenta sa opisani v edin item.
                                        # f <- function(x){x[[1]] %in% args}
    f <- function(x){           # tozi code tryabva da se ischisti; dokato otkriya tazi rabota
                                # (t.e. che e neobchodimo as.character) se omotach.
        wrk0 <- as.character(x[[1]])                # x[[1]] is tagged with Rd_tag or similar!
        if(wrk0 %in% args)
            return(TRUE)

        wrk <- strsplit(wrk0,",[ ]*")
        if(!is.character(wrk[[1]])){
            warning("wrk[[1]] is not a character vector! ", wrk)
            return(FALSE)
        }
        wrk <- any( wrk[[1]] %in% args )                    # ima li nuzhda ot trim na blanks?
        wrk
    }
    sel <- !sapply(rdargs, f)

    ## deal with "..." arg
    if("..." %in% args || "\\dots" %in% args){  # since formals() represents ... by "..."
        f2 <- function(x){
            if(is.list(x[[1]]) && length(x[[1]])>0 &&
               attr(x[[1]][[1]],"Rd_tag") == "\\dots")
                TRUE
            else
                FALSE
        }
        i2 <- sapply(rdargs, f2)
        sel[i2] <- FALSE
    }

    rdargs[sel] <- NULL   # keeps attributes (even if 0 or 1 elem remain).
    rdargs
}

Rdo_get_arg <- function(rd,arg){
    wrk <- Rdo_get_args(rd,arg)
    wrk[[1]]
}


Rdo_args2txt_list <- function(x,arg,...){
    rdo <- Rd_fun(x)
    if(missing(arg)){
        tmparg <- tools:::.Rd_get_argument_names(rdo)

        # correct for merged descriptions...
        arg <- character(0)
        for(s in tmparg){
            arg <- c(arg, strsplit( as.character(s), ",[ ]*")[[1]])
        }
    }

    res <- list()
    for(a in arg)
        res[[a]] <- Rdo_args2txt(rdo,a,...)
    res
}

Rdo_args2txt <- function(rdo,arg,title="Hhh",name="Aa",type="text"){
    wrk <- Rdo_get_args(rdo,arg)
    wrk2 <- Rdo_create(arguments=wrk,title=title,name=name)

    res <- Rd_help2txt(wrk2, keep_section = "\\arguments"
                           , omit_sec_header = TRUE
                       )
                              # nay-dobre e da ima programka, koyato da macha izlishni poleta!
    res <- paste(res,collapse="\n")
    res
}


                          # based on print.help_files_with_topic() in the sources of R-2.10.0.
Rd_help2txt <- function(x, topic, pkgname=""
                        , help_type="text"
                        , verbose=FALSE
                        , try.all.packages=FALSE
                        , keep_section = TRUE
                        , omit_sec_header = FALSE
                        ){
    rdo <- Rd_fun(x, topic=topic, pkgname=pkgname
                  , help_type        = help_type
                  , verbose          = verbose
                  , try.all.packages = try.all.packages
                  , keep_section     = keep_section
                  )

    temp <- tools::Rd2txt(rdo, out=tempfile("Rtxt"), package=pkgname)

    res <- readLines(temp) # note: temp is a (temporary) file name.
    unlink(temp)

                                         # krapka, iztrii title i/ili name ako ne sa poiskani.
    iomit <- numeric(0)
             # the code below assumes that each item is on one line followed by  a blank line.
    if(!("\\title" %in% keep_section))
        iomit <- c(iomit,1:2)

                           # !!! Poleto "name" ne vliza v teksta, zatoca nyama nuzhda ot tova!
                           # if(!("\\name" %in% keep_section))
                           #     iomit <- c(iomit,3:4)
    if(isTRUE(omit_sec_header))
        iomit <- c(iomit,3:4)

                                  # file.show(temp, title = gettextf("R Help on '%s'", topic),
                                  #                 delete.file = TRUE)
    if(length(iomit)>0)
        res <- res[-iomit]          # !!! ??? MNOGO GRUBA KRAPKA za omit-vane na title/name
    res
}


# based on print.help_files_with_topic() in the sources of R-2.10.0.
Rd_fun <- function(x, topic, pkgname   = ""
                    , help_type        = "text"
                    , verbose          = FALSE
                    , try.all.packages = FALSE
                    , keep_section     = TRUE
                   ){
    rdo <- NULL         # prepare the "Rd" object rdo; # is it better to check with "inherit"?
    if(class(x) == "Rd"){  # if(inherits(file, "Rd")) ...
        rdo <- x
    }else{
        if(class(x) != "help_files_with_topic" ){
               # The following comments baffle me now. Does `do.call' resolve the issues?
               #
               # help returns an object of class "help_files_with_topic" the
               #  eval(substitute()) wrapper (I saw it in tkGUI, vzh sasto help.R, sasto:
               #  .tryHelp in question.R) is needed to cover the case when x is a
               #  function. Without this wrapper the result is not correct.

               # Izglezhda, che bez substitute() argumentat se evvaluate-va na nepodochodyasto
               #  myasto.  If x is a name of a function, then the wrapper is not needed.

                       # wrk <- eval(substitute(help(x, help_type=help_type
                       #            , verbose=verbose
                       #            , try.all.packages=try.all.packages)))

                       # cat("KUKUKUUUU: ", substitute(x), "   class(x): ", class(x), "\n\n" )

            wrk <- do.call("help",list(x, help_type=help_type
                                        , verbose=verbose
                                        , try.all.packages=try.all.packages))
            x <- wrk
        }
        ## Check for errors! ???

        if(class(x) == "help_files_with_topic"){
                                                  # from print.help_files_with_topic in help.R
                                                  #
                                                  # browser <- getOption("browser")
            topic <- attr(x, "topic")
            type <- attr(x, "type")
            paths <- as.character(x) # removes attributes of x.
                                     # If more matches are found will `paths' have length > 1?
            file <- paths

                                             # !!! check for lenght(paths)==0  !!!! ??
                                             # but no error is raized, rdo simply remain NULL.
                        # the following commands are probably copied from utils:::.getHelpFile
            path <- dirname(file)
            dirpath <- dirname(path)
            pkgname <- basename(dirpath)
            RdDB <- file.path(path, pkgname)

                             # cat("\n\nx is: "    ,unclass(x)                      ,"\n\n\n")
                             # cat("paths is: ",paths                      ,"\n")
                             # cat("file is: ", file                       ,"\n")
                             # cat("path is: ", path                       ,"\n")
                             # cat("RdDB is: ", paste(RdDB, "rdx", sep="."),"\n")

            if(file.exists(paste(RdDB, "rdx", sep="."))) {
                rdo <- tools:::fetchRdDB(RdDB, basename(file))
                                                          # a debugging message, remove later!
                   # cat("Class of object returned by \"tools:::fetchRdDB: ", class(rdo),"\n")
                   # really returns "Rd".
            }
        }
    }
    if(is.null(rdo))                             # todo: should someting less radical be done?
        stop("rdo object is NULL!")

    if(is.character(keep_section) && length(keep_section)>0){
        tags <- tools:::RdTags(rdo)
        keep_tags <- unique(c("\\title","\\name",keep_section))
        rdo[which(!(tags %in% keep_tags))] <-  NULL
    }

    rdo
}


# wrk[[1]] <- Rd_name("Random name")
# wrk[[2]] <- Rd_title("Kukurigu")
# Rd2txt(wrk)
# Rd2HTML(wrk)

# 2011-10-30 premestvam v Rdpack funktsiite napisani sled publikuvane na originalnata
# versiya.
