### Set up new environment for equation list and equation number counter
.hwriterGlobalEnv <- new.env()


newPage <- function(filename, dirname = NULL, title = filename,
                    doctype = "<!DOCTYPE html>\n",
                    link.javascript = NULL, link.css = NULL, css = NULL,
                    head = NULL, head.attributes = NULL,
                    body.attributes = NULL)
{
    ## today <- format(strptime(date(), "%a %b %d %H:%M:%S %Y"), "%B %d, %Y")

    ## Create equation numbers and list for labels
    ## Original version used global assignment
    ## hwriterEquation <<- 0
    ## hwriterEquationList <<- character(0)
    assign("hwriterEquation", 0, .hwriterGlobalEnv)
    assign("hwriterEquationList", character(0), .hwriterGlobalEnv)

    if (!is.null(dirname)) {
        if (!file.exists(dirname))
            dir.create(dirname, recursive = TRUE, showWarnings = FALSE)
        filename <- file.path(dirname, filename)
    }
    page <- file(filename, "wt")
    if (!is.null(link.javascript))
        link.javascript <- paste(hmakeTag("script", type = "text/javascript",
                                          src = link.javascript),
                                 collapse = "\n")
    if (!is.null(link.css))
        link.css <- paste(hmakeTag("link", rel = "stylesheet",
                                   type = "text/css", href = link.css),
                          collapse = "\n")
    if (!is.null(css))
        css <- paste(hmakeTag("style", css), collapse = "\n")
    head <- paste(hmakeTag("title", title), head, link.javascript,
                  link.css, css, sep = "\n")
    head <- do.call(hmakeTag, c(list("head", head, newline = TRUE),
                                head.attributes))
    bodyStart <- do.call(hmakeTag, c(list("body", NULL), body.attributes))
    bodyStart <- substr(bodyStart, 1, regexpr("</body>", bodyStart) - 1)
    hwrite(paste(doctype, "<html ", ">", head, bodyStart, sep = ""), page)
    page
}
