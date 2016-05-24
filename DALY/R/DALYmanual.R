## Open DALY Manual PDF

DALYmanual <-
function(){
  pdf <- system.file("doc/DALYmanual.pdf", package = "DALY")

  if (.Platform$OS.type == "windows"){
    shell.exec(pdf)
  } else {
    system(paste(shQuote(getOption("pdfviewer")), shQuote(pdf)), wait = FALSE)
  }
}