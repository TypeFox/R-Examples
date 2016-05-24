TANOVAmanual <- function(view = TRUE) 
{
    f <- system.file("doc", "TANOVAmanual.pdf", package = "TANOVA")
    if (view) {
        if (.Platform$OS.type == "windows") 
            shell.exec(f)
        else system(paste(Sys.getenv("R_PDFVIEWER"), f, "&"))
    }
    return(f)
}

