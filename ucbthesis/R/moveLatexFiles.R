moveLatexFiles <- function(tmpExtensions = c("aux", "bbl", "bcf", "blg", "lof",
                               "log", "lot", "run.xml", "toc")) {
  
  filesInDir <- 
    unlist(str_split(system("ls", intern=TRUE), pattern = "[[:space:]]+"))
  
  sapply(filesInDir, function(file) {
    extension <- str_split(file, pattern="[.]", n=2)[[1]][2]
    if (extension %in% tmpExtensions) {
      cat("Moving", file, "to ./tmp/ ...\n")
      system(str_c("mv ", file, " ./tmp/", file))
    }
    invisible(file)
  })
  
  invisible()
}
