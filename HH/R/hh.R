hh.old <-
  function(file) file.path(options()$HH.ROOT.DIR, file)


hh <- function(file) {
  ## file <- "datasets/abcde.dat"
  ## file <- "conc/code/normpdf.s"
  dataset.name <- sub("datasets/(.*)\\.dat", file, replacement="\\1")
  
  if (dataset.name != file) { ## dataset name
    cat("Access to the HH book files has been changed.\n")
    cat(paste0("Use\n   data(",
               dataset.name,
               ")\n"))
    cat("instead of\n")
    cat(paste0("   ", dataset.name, " <- read.table(hh('",
               file,
               "'))\n"))
  }
  else { ## R source or transcript
    cat("File\n")
    cat(paste0("   ", hh.old(file), "\n"))
    cat("is no longer distributed.  Instead, open file\n")
    
    chapter.name <- sub("(.*)/code/.*", file, replacement="\\1")
    chapter.scripts.files <- list.files(hh.old("scripts"))
    script.name <- grep(chapter.name, chapter.scripts.files, value=TRUE)
    cat(paste0("   ", hh.old(""), "scripts/", script.name, "\n"))
    
    cat("in a text editor.\n")
  }
  stop("")
}

## read.table(hh("datasets/abcde.dat"))
## hh("conc/code/normpdf.s")
## hh("conc/code/normpdf.r")
