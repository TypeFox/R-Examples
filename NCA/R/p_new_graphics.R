p_new_pdf <- 
function (out.prefix, name1, name2, name3="", paper="a4") {
  if (name3 == "") {
    file.name <- paste(out.prefix, name1, name2, "pdf", sep=".")
  } else {
    file.name <- paste(out.prefix, name1, name2, name3, "pdf", sep=".") 
  }
  if (paper == "A4r") {
    pdf(file.name, paper=paper, width = 0, height = 0)
  } else {
    pdf(file.name, paper=paper)
  } 
}

p_new_window <-
function (title="", width=7, height=7) {
  # RStudio can't handle multiple windows (yet?)
  cmd <- commandArgs(trailingOnly=FALSE)[1]
  if (substr(tolower(cmd), nchar(cmd)-6, nchar(cmd)) == "rstudio" ||
      substr(tolower(cmd), nchar(cmd)-10, nchar(cmd)) == "rstudio.exe") {
    return ()
  }

  # Don't make windows smaller than 7x7
  width <- max(7, width)
  height <- max(7, height)
  dev.new(title=title, width=width, height=height, rescale='fixed')
}