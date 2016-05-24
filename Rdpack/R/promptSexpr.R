                     # like promptPackage but uses \Sexpr{} for automatic information updating
promptPackageSexpr <- function(package, filename= NULL, final = TRUE, overview = FALSE,
                               bib = TRUE){   # 2013-03-30 neew arg. `bib'
    paste0 <- function(...) paste(..., sep="")

    overv <- if(overview | !final)  ""  else  "% "
    hide <- "\\Sexpr[stage=build,results=hide]"
    pd <- paste0("{pd <- packageDescription(\"", package, "\")}")
    lb <- paste0("{lb <- library(help=\"", package, "\", character.only=TRUE)}")
    lbflag <- "lbflag <- !is.null(lb$info[[2]])"

    Rdtxt <-
        list(  paste0("\\name{", package, "-package}")
             , paste0("\\alias{", package, "-package}")
             , paste0("\\alias{", package, "}")
             , paste0("\\docType{", "package", "}")
             , "\\title{"
             , paste0("  ", hide, c(pd, lb))     # init. commands, cannot be outside sections.
             , "  \\Sexpr[stage=build]{pd$Title}"
             , "}"
             , "\\description{"
             , "  \\Sexpr[stage=build]{pd$Description}"
             , "}"
             , "\\details{"
             , "   \\tabular{ll}{"
             , "   Package:  \\tab \\Sexpr[stage=build]{pd$Package}\\cr"
             , "   Type:     \\tab \\Sexpr[stage=build]{pd$Type}\\cr"
             , "   Version:  \\tab \\Sexpr[stage=build]{pd$Version} \\cr"
             , "   Date:     \\tab \\Sexpr[stage=build]{pd$Date}\\cr"
             , "   License:  \\tab \\Sexpr[stage=build]{pd$License}\\cr"
             , "   LazyLoad: \\tab \\Sexpr[stage=build]{pd$LazyLoad}\\cr"
             , "   Built:    \\tab \\Sexpr[stage=build]{pd$Built}\\cr"
             , "   }"
             , ""
             , "   Index:"
             , paste0("  \\Sexpr[stage=build,results=rd]{paste(\"\\\\\\\\preformatted{\""
                      , ", paste(if(!is.null(lb$info[[2]])) lb$info[[2]] else \"No entries\",collapse=\"\\n\"), \"}\", sep=\"\")}")
             , ""
             , paste0(overv, c("~~ An overview of how to use the package, including the most important ~~"
                               , "~~ functions ~~"))
             , "}"
             , "\\author{"
             , "  \\Sexpr[stage=build]{pd$Author}"
             , ""
             , paste0("Maintainer: ", "\\Sexpr[stage=build]{pd$Maintainer}")
             , "}"
             , "\\references{"
             , if(bib) "% bibentry:all"
                 else  "% ~~ Literature or other references for background information ~~"
             , "}"
             , "% ~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~"
             , "% ~~ the R documentation directory ~~"
             , "\\keyword{ package }"
             , "% \\seealso{"
             , "% ~~ Optional links to other man pages, e.g. ~~"
             , "% ~~ \\code{\\link[<pkg>:<pkg>-package]{<pkg>}} ~~"
             , "% }"
             , "% \\examples{"
             , "% ~~ simple examples of the most important functions ~~"
             , "% }"
             )
    if(is.null(filename))
        filename <- paste0(package, "-package.Rd")

    if(is.na(filename))
        return(Rdtxt)

    cat(unlist(Rdtxt), file = filename, sep = "\n")

    message(gettextf("Created file named %s.", sQuote(filename)),
            "\n",
            gettext("Edit the file and move it to the appropriate directory."),
            domain = NA)

    invisible(filename)
}
