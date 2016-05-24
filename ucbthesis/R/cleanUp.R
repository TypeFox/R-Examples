# ls's through the current working directory and removes all files/folders that
# don't match the extensions in safeExtensions. In particular, this will delete
# the folder ./tmp/
cleanUp <- function(safeExtensions) {
  filesInDir <- 
    unlist(str_split(system("ls", intern=TRUE), pattern = "[[:space:]]+"))
  
  sapply(filesInDir, function(file) {
    extension <- str_split(file, pattern="[.]", n=2)[[1]][2]
    if (!(extension %in% safeExtensions)) {
      cat("Removing", file, "...\n")
      system(str_c("rm -rf ", file))
    }
    invisible(file)
  })
  
  invisible()
}
