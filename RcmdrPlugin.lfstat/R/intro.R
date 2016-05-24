lfstat_intro <- function(){
    pdfviewer <- getOption("pdfviewer")
    introfile <- system.file("doc","intro_lfstat.pdf",package = "RcmdrPlugin.lfstat")
    if (identical(pdfviewer, "false")) {
        }
    else if (.Platform$OS.type == "windows" && identical(pdfviewer, 
        file.path(R.home("bin"), "open.exe"))) 
        shell.exec(introfile)
    else system2(pdfviewer, shQuote(introfile), wait = FALSE)
}

