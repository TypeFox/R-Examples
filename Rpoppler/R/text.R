PDF_text <-
function(file)
{
    x <- PDF_doc_from_file(file)

    y <- .Call("PDF_text", x$doc)

    y

}
