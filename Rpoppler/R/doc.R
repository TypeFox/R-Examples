## <NOTE>
## No support for password protected PDF files for now.
## </NOTE>

PDF_doc <-
function(file)
{
    file <- file_path_to_URI(file)[1L]
    
    doc <- .Call("PDF_doc", file)
    y <- list(doc = doc, file = file)
    class(y) <- "PDF_doc"
    
    y
}

print.PDF_doc <-
function(x, ...)
{
    writeLines(sprintf("A reference to PDF file '%s'", x$file))
    invisible(x)
}


PDF_doc_from_file <-
function(file)
{
    if(inherits(file, "PDF_doc"))
        file
    else
        PDF_doc(file)
}

file_path_to_URI <-
function(path)
{
    sprintf("%s%s",
            if(.Platform$OS == "windows") "file:///" else "file://",
            normalizePath(path))
}
