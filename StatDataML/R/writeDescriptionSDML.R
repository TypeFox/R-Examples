writeDescriptionSDML <- function(title = "RDataset",
                                 source = "R",
                                 version = " ",
                                 date = NULL,
                                 comment = "",
                                 properties = NULL,
                                 file = "",
                                 textdata,
                                 sep,
                                 na.string,
                                 null.string,
                                 posinf.string,
                                 neginf.string,
                                 nan.string,
                                 true,
                                 false)
{
    ## writes the description tag to file

    catSDML("<description>\n", file = file)
    catSDML("<title>", markup(title), "</title>\n", file = file)
    catSDML("<source>", markup(source), "</source>\n", file = file)
    cat("<date>", markup(date), "</date>\n", file = file, append = TRUE, sep = " ")
    catSDML("<version>", markup(version), "</version>\n", file = file)
    catSDML("<comment>", markup(comment), "</comment>\n", file = file)
    sdmlib <- system.file(package = "StatDataML")
    sdmlib <- substr(sdmlib, 1, nchar(sdmlib)-10)
    pkgversion <- packageDescription("StatDataML", lib.loc=sdmlib)$Version

    catSDML("<creator>R-", R.version$major, ".",  R.version$minor,
            ":StatDataML_", pkgversion, "</creator>\n", file = file)
    if (!is.null(properties)) {
        catSDML("<properties>\n", file = file)

        if (!is.list(properties))
            properties <- list(properties)
        writeListArraySDML(properties, file = file, textdata = textdata,
                           sep = sep, na.string = na.string,
                           null.string = null.string, posinf.string = posinf.string,
                           neginf.string = neginf.string, nan.string = nan.string,
                           true = true, false = false)

        catSDML("</properties>\n", file = file)
    }
    catSDML("</description>\n", file = file)
}
