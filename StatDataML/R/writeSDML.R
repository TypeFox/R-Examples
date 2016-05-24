writeSDML <- function(x,
                      file = "",
                      textdata = NULL,
                      dtd = NULL,
                      sep = " &#x000A;&#x000D;",
                      na.string = "NA",
                      null.string = "NULL",
                      posinf.string = "+Inf",
                      neginf.string = "-Inf",
                      nan.string = "NaN",
                      true = "1",
                      false = "0",
                      title = deparse(substitute(x)),
                      source = "R",
                      version = " ",
                      date = NULL,
                      comment = " ",
                      properties = NULL)
{
    if (is.null(date)) date <- date()
    if (is.null(dtd))
        dtd <- system.file("dtd/StatDataML.dtd", package = "StatDataML")[1]

    cat("<?xml version=\"1.0\"?>\n", file = file, sep = "")
    catSDML("<!DOCTYPE StatDataML PUBLIC \"StatDataML.dtd\" \"", dtd,
            "\" >\n", file = file)
    catSDML("<StatDataML xmlns=\"http://www.omegahat.org/StatDataML/\">\n",
            file = file)

    writeDescriptionSDML(title = title, source = source,
                         version = version, date = date,
                         comment = comment, file = file,
                         properties = properties,
                         textdata = textdata, sep = sep,
                         na.string = na.string, null.string = null.string,
                         posinf.string = posinf.string,
                         neginf.string = neginf.string,
                         true = true, false = false, nan.string = nan.string)

    writeDatasetSDML(x, file = file, textdata = textdata, sep = sep,
                     na.string = na.string, null.string = null.string,
                     posinf.string = posinf.string, neginf.string = neginf.string,
                     true = true, false = false, nan.string = nan.string)

    catSDML("</StatDataML>\n", file = file)
}
