### default settings
.CacheLength <- 0
.CacheFiles <- NULL
.alreadyDefinedPkgs <- NULL
.tobeDefinedPkgs <- NULL
.numberofRstyleDefs <- NULL

.SweaveListingOptions <- list(
Rset = list("fancyvrb" = "true", "escapechar" = "`",
        "extendedchars" = "false",
        "language" = "R", "basicstyle" = "{\\color{Rcolor}\\small}",
        "keywordstyle" = "{\\bf\\color{Rcolor}}",
        "commentstyle" = "{\\color{Rcommentcolor}\\ttfamily\\itshape}",
        "literate" = paste("{\"}{\\texttt{\"}}1",
                           "{<-}{{$\\leftarrow$}}2",
                           "{<<-}{{$\\twoheadleftarrow$}}2",#"%\n",
                           "{~}{{$\\sim$}}1", "{<=}{{$\\leq$}}2",#"%\n",
                           "{>=}{{$\\geq$}}2", "{^}{{$\\scriptstyle\\wedge$}}1", sep=""),
                     ## ~,^,<=, >= as suggested by Frank Harrell
        "alsoother" = "{$}", "alsoletter" = "{.<-}",
        "otherkeywords" = "{!,!=,~,$,*,\\&,\\%/\\%,\\%*\\%,\\%\\%,<-,<<-,/}",
        "escapeinside" = "{(*}{*)}" ## as suggested by Frank Harrell
        ),
Rdset = list("fancyvrb" = "true", "language" = "Rd", 
             "keywordstyle" = "{\\bf}",
             "basicstyle" = "{\\color{black}\\footnotesize}",
               "commentstyle" = "{\\ttfamily\\itshape}",
               "alsolanguage" = "R"),
Rin = list("style" = "Rstyle", "fancyvrb" = "true",
           "basicstyle" = "\\color{Rcolor}\\small"),
Rout = list("fancyvrb" = "false", "basicstyle" = "\\color{Routcolor}\\small"),
Rcode = list("style" = "Rstyle", "fancyvrb" = "true",
             "fontshape"= "sl", "basicstyle" = "\\color{Rcolor}"),
Rcolor  = c(0,0.5,0.5),
RRecomdcolor  = c(0,0.6,0.4),
Rbcolor       = c(0,0.6,0.6),
Routcolor     = c(0.461,0.039,0.102),
Rcommentcolor = c(0.101,0.043,0.432),
pkv = "2.5",
pkg = "distr",
Keywordstyle = "{\\bf}",
Recomd.Keywordstyle = "{\\bf\\color{RRecomdcolor}}",
interm.Keywordstyle = "{\\bf\\color{Rbcolor}}",
overwrite = FALSE,
intermediate = TRUE,
inSweave = FALSE,
fromRForge = TRUE,
base.url = paste("http://r-forge.r-project.org/scm/viewvc.php/",
                            "*checkout*/pkg/", 
                            sep = ""),
addRset = TRUE,
addRdset = TRUE,
addRinset = TRUE,
addRoutset = TRUE,
addRcodeset = TRUE,
fileCommand = "\\def\\file#1{{\\tt #1}}",
pkgCommand = "\\def\\pkg#1{{\\tt \"#1\"}}"
)

SweaveListingOptions <- function(...) {
    if (nargs() == 0) return(.SweaveListingOptions)
    current <- .SweaveListingOptions
    temp <- list(...)
    if (length(temp) == 1 && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
            list = temp <- arg,
            character = return(.SweaveListingOptions[arg]),
            stop("invalid argument: ", sQuote(arg)))
    }
    if (length(temp) == 0) return(current)
    n <- names(temp)
    if (is.null(n)) stop("options must be given by name")
    changed <- current[n]
    current[n] <- temp
#    if (sys.parent() == 0)
        env <- asNamespace("SweaveListingUtils")
#    else
#        env <- parent.frame()
    assign(".SweaveListingOptions", current, envir = env)

    invisible(current)
}

getSweaveListingOption <- function(x) SweaveListingOptions(x)[[1]]
SweaveListingoptions <- SweaveListingOptions
