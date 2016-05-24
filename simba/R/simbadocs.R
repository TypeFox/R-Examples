"simbadocs" <- 
function (doc = c("NEWS", "ChangeLog", "mps-coefficients.pdf", "simba_manual.pdf")) 
{
    doc <- match.arg(doc)
    if (length(grep(".pdf", doc)) > 0) {
        doc <- file.path(system.file(package = "simba"), "doc", 
            doc)
        if (.Platform$OS.type == "windows") 
            shell.exec(doc)
        else system(paste(getOption("pdfviewer"), doc, "&"))
    }
    else {
        file.show(system.file(package = "simba", doc))
    }
}