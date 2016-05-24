hsdardocs <- function(doc)
{
if (doc == "References.pdf")
{
  doc <- file.path(system.file(package = "hsdar"), "doc", doc) 
  if (.Platform$OS.type == "windows")
  {
    shell.exec(doc)
  } else {
    system(paste(getOption("pdfviewer"), doc, "&"))
  }
}
if (toupper(doc) == "COPYRIGHT")
{
  doc <- file.path(system.file(package = "hsdar"), "COPYRIGHTS") 
  if (.Platform$OS.type == "windows")
  {
    file.show(doc)
  } else {
    system(paste(getOption("pager"), doc, "&"))
  }
}
if (doc == "hsdar-intro.pdf")
{
  doc <- file.path(system.file(package = "hsdar"), "doc", doc) 
  if (.Platform$OS.type == "windows")
  {
    shell.exec(doc)
  } else {
    system(paste(getOption("pdfviewer"), doc, "&"))
  }
}
}
