## Author: Ingo Feinerer
## Readers

FunctionGenerator <-
function(x)
{
    class(x) <- c("FunctionGenerator", "function")
    x
}

getReaders <-
function()
    c("readDOC", "readPDF", "readPlain", "readRCV1", "readRCV1asPlain",
      "readReut21578XML", "readReut21578XMLasPlain", "readTabular",
      "readTagged", "readXML")

prepareReader <-
function(readerControl, reader = NULL, ...)
{
    if (is.null(readerControl$reader))
        readerControl$reader <- reader
    if (inherits(readerControl$reader, "FunctionGenerator"))
        readerControl$reader <- readerControl$reader(...)
    if (is.null(readerControl$language))
        readerControl$language <- "en"
    readerControl
}

processURI <-
function(uri)
{
    uri <- as.character(uri)
    if (identical(substr(uri, 1, 7), "file://"))
        uri <- substr(uri, 8, nchar(uri))
    uri
}

# readDOC needs antiword installed to be able to extract the text
readDOC <-
function(AntiwordOptions = "")
{
    stopifnot(is.character(AntiwordOptions))

    function(elem, language, id) {
        uri <- processURI(elem$uri)
        content <- system2("antiword",
                           c(AntiwordOptions, shQuote(normalizePath(uri))),
                           stdout = TRUE)
        PlainTextDocument(content, id = basename(elem$uri), language = language)
    }
}
class(readDOC) <- c("FunctionGenerator", "function")

readPDF <-
function(engine = c("xpdf", "Rpoppler", "ghostscript", "Rcampdf", "custom"),
         control = list(info = NULL, text = NULL))
{
    stopifnot(is.character(engine), is.list(control))

    engine <- match.arg(engine)

    pdf_info <-
        switch(engine,
               xpdf = function(x) pdf_info_via_xpdf(x, control$info),
               Rpoppler = Rpoppler::PDF_info,
               ghostscript = pdf_info_via_gs,
               Rcampdf = Rcampdf::pdf_info,
               custom = control$info)

    pdf_text <-
        switch(engine,
               xpdf = function(x) system2("pdftotext",
                                          c(control$text, shQuote(x), "-"),
                                          stdout = TRUE),
               Rpoppler = Rpoppler::PDF_text,
               ghostscript = pdf_text_via_gs,
               Rcampdf = Rcampdf::pdf_text,
               custom = control$text)

    if (!is.function(pdf_info) || !is.function(pdf_text))
        stop("invalid function for PDF extraction")

    function(elem, language, id) {
        uri <- processURI(elem$uri)
        meta <- pdf_info(uri)
        content <- pdf_text(uri)
        PlainTextDocument(content, meta$Author, meta$CreationDate, meta$Subject,
                          meta$Title, basename(elem$uri), language, meta$Creator)
     }
}
class(readPDF) <- c("FunctionGenerator", "function")

readPlain <-
function(elem, language, id) {
    if (!is.null(elem$uri))
        id <- basename(elem$uri)
    PlainTextDocument(elem$content, id = id, language = language)
}

readXML <-
function(spec, doc)
{
    stopifnot(is.list(spec), inherits(doc, "TextDocument"))

    function(elem, language, id) {
        tree <- XML::xmlParse(elem$content, asText = TRUE)
        content(doc) <- if ("content" %in% names(spec))
            .xml_content(tree, spec[["content"]])
        else
            XML::xmlTreeParse(elem$content, asText = TRUE)
        for (n in setdiff(names(spec), "content"))
            meta(doc, n) <- .xml_content(tree, spec[[n]])
        XML::free(tree)
        if (!is.null(elem$uri))
            id <- basename(elem$uri)
        if (!length(meta(doc, "id")))
            meta(doc, "id") <- as.character(id)
        if (!length(meta(doc, "language")))
            meta(doc, "language") <- as.character(language)
        doc
    }
}
class(readXML) <- c("FunctionGenerator", "function")

RCV1Spec <-
    list(author = list("unevaluated", ""),
         datetimestamp = list("function", function(node)
           as.POSIXlt(as.character(XML::getNodeSet(node, "/newsitem/@date")),
                      tz = "GMT")),
         description = list("unevaluated", ""),
         heading = list("node", "/newsitem/title"),
         id = list("attribute", "/newsitem/@itemid"),
         origin = list("unevaluated", "Reuters Corpus Volume 1"),
         publisher = list("attribute",
           "/newsitem/metadata/dc[@element='dc.publisher']/@value"),
         topics = list("attribute",
           "/newsitem/metadata/codes[@class='bip:topics:1.0']/code/@code"),
         industries = list("attribute",
           "/newsitem/metadata/codes[@class='bip:industries:1.0']/code/@code"),
         countries = list("attribute",
           "/newsitem/metadata/codes[@class='bip:countries:1.0']/code/@code"))
readRCV1 <- readXML(spec = RCV1Spec, doc = XMLTextDocument())
readRCV1asPlain <-
readXML(spec = c(RCV1Spec, list(content = list("node", "/newsitem/text"))),
        doc = PlainTextDocument())

Reut21578XMLSpec <-
    list(author = list("node", "/REUTERS/TEXT/AUTHOR"),
         datetimestamp = list("function", function(node)
           strptime(sapply(XML::getNodeSet(node, "/REUTERS/DATE"),
                           XML::xmlValue),
                    format = "%d-%B-%Y %H:%M:%S",
                    tz = "GMT")),
         description = list("unevaluated", ""),
         heading = list("node", "/REUTERS/TEXT/TITLE"),
         id = list("attribute", "/REUTERS/@NEWID"),
         topics = list("attribute", "/REUTERS/@TOPICS"),
         lewissplit = list("attribute", "/REUTERS/@LEWISSPLIT"),
         cgisplit = list("attribute", "/REUTERS/@CGISPLIT"),
         oldid = list("attribute", "/REUTERS/@OLDID"),
         origin = list("unevaluated", "Reuters-21578 XML"),
         topics_cat = list("node", "/REUTERS/TOPICS/D"),
         places = list("node", "/REUTERS/PLACES/D"),
         people = list("node", "/REUTERS/PEOPLE/D"),
         orgs = list("node", "/REUTERS/ORGS/D"),
         exchanges = list("node", "/REUTERS/EXCHANGES/D"))
readReut21578XML <- readXML(spec = Reut21578XMLSpec, doc = XMLTextDocument())
readReut21578XMLasPlain <-
readXML(spec = c(Reut21578XMLSpec,
                 list(content = list("node", "/REUTERS/TEXT/BODY"))),
        doc = PlainTextDocument())

readTabular <-
function(mapping)
{
    stopifnot(is.list(mapping))
    function(elem, language, id) {
        meta <- lapply(mapping[setdiff(names(mapping), "content")],
                       function(m) elem$content[, m])
        if (is.null(meta$id))
            meta$id <- as.character(id)
        if (is.null(meta$language))
            meta$language <- as.character(language)
        PlainTextDocument(elem$content[, mapping$content], meta = meta)
    }
}
class(readTabular) <- c("FunctionGenerator", "function")

readTagged <-
function(...)
{
    args <- list(...)
    function(elem, language, id) {
        if (!is.null(elem$content)) {
            con <- textConnection(elem$content)
            on.exit(close(con))
        } else
            con <- elem$uri

        if (!is.null(elem$uri))
            id <- basename(elem$uri)

        a <- c(list(con = con, meta = list(id = id, language = language)), args)

        do.call(NLP::TaggedTextDocument, a)
    }
}
class(readTagged) <- c("FunctionGenerator", "function")
