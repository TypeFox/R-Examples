# Identify functions that should not be tested, including quit(), q() and install.packages().
# Replace all occurences with commented function
clean_script <- function(file, message=FALSE, target=c("rmd", "html")){
  target <- match.arg(target)
  noRunFunctions <- c("quit", "readClipboard", "RSiteSearch", "install.packages", "update.packages", "savehistory", "loadhistory", "findFn", "vignette", "png", "trellis.device", "dev.off", "save", "load")
  ptn <- sprintf("^!( *#* *).*%s", paste0("(", noRunFunctions, "\\(", ")", collapse="|"))
  
  ptn <- sprintf("^[^ #]*.*(%s)\\(", paste(noRunFunctions, collapse="|"))
  if(!file.exists(file)) stop("file doesn't exist")
  txt <- scan(file, what="character", sep="\n", quiet=TRUE, blank.lines.skip = FALSE)
  
  # Remove fancy quotes
  txt <- gsub("â€™", "'", txt, useBytes=TRUE)  ## fancy single quote
  txt <- gsub("’", "'", txt, useBytes=FALSE)  ## fancy single quote
  
  
  # Remove hard space
  txt <- gsub("ï»¿", "", txt, useBytes=TRUE)  ## weird hard space
  
  nonAsciiSpace <- "\xef\xbb\xbf"
  Encoding(nonAsciiSpace) <- "latin1"
  Encoding(txt) <- "latin1"
  txt <- gsub(nonAsciiSpace, " ", txt)
  Encoding(txt) <- "ASCII"
  
  
  
  
  # Comment lines containing noRunFunctions
  idx <- grep(ptn, txt, perl=TRUE)
  if(length(idx) > 0) {
    if(target == "html"){
      txt[idx] <- paste("##", txt[idx], sep=" ")
    } else {
      txt[idx] <- paste("\\dontrun{", txt[idx], "}", sep="\n")
    }
    if(message){
      cat(txt[idx], sep="\n")
      cat(txt[idx], sep="\n")
    }
  }
  
  Encoding(txt) <- "ascii"
  
  txt
}



# Run clean_script for all scripts in inst/scripts/1-orig, putting cleaned version in inst/scripts/2-clean
clean_all_scripts <- function(scriptPath, outPath, comment=TRUE, message=TRUE){
  if(!file.exists(scriptPath)) stop("script path doesn't exist")
  scripts <- normalizePath(list.files(scriptPath, pattern="\\.R$", full.names = TRUE))
  clean <- lapply(scripts, clean_script, message=message)
  ret <- sapply(seq_along(scripts), function(i){
    newFile <- file.path(outPath, basename(scripts[i]))
    if(message) message(newFile)
    txt <- paste(clean[[i]], collapse="\n")
    writeLines(txt, con=newFile)
    newFile
  })
  invisible(ret)
}

clean_all_scripts(scriptPath = "rfordummies/inst/scripts/1-orig", outPath = "rfordummies/inst/scripts/2-clean")



