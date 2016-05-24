essQuery <- function(essentia, aq = "", flags = "") 
{
  if (paste0(essentia,aq,flags) != "") {
    commandcount <- 1
    if (substr(essentia, 1, 3) != "ess") {
      flags <- aq
      aq <- essentia
      essentia <- "ess exec"
    }
    if (((substr(essentia, 1, 8) == "ess exec") && (!grepl("#Rignore", 
                                                           flags))) || ((substr(essentia, 1, 10) == "ess stream" || 
                                                                         substr(essentia, 1, 9) == "ess query") && (grepl("#Rinclude", 
                                                                                                                          flags)))) {
      if ((substr(essentia, 1, 10) == "ess stream") || (substr(essentia, 
                                                               1, 8) == "ess exec")) {
        aq <- paste(aq, "; echo 'RSTOPHERE'", sep = "")
      }
      else {
        flags <- paste("; echo 'RSTOPHERE'", flags, sep = " ")
      }
      line <- paste(essentia, "\"", aq, "\"", flags, sep = " ")
      colspec <- TRUE
      if (grepl("-notitle", line)) {
        colspec <- FALSE
      }
      varname <- "command"
      if (grepl("#R#", line)) {
        titleindex <- grepRaw("#R#", line, all = TRUE)
        varname <- substr(line, titleindex[[1]] + 3, titleindex[[2]] - 
                            1)
        remove(titleindex)
      }
      pipedline <- pipe(line, open = "r")
      t2 <- read.csv(pipedline, header = colspec, sep = ",", 
                     quote = "\"'", comment.char = "#", blank.lines.skip = FALSE, 
                     allowEscapes = TRUE, skip = 0)
      close(pipedline)
      index <- 1
      t3 <- NULL
      separate <- grepl(" #Rseparate", line)
      for (file in seq(1, length(which(t2[, 1] == "RSTOPHERE")), 
                       1)) {
        if (separate) {
          assign(sprintf("%s%i", varname, commandcount), 
                 t2[index:(which(t2[, 1] == "RSTOPHERE")[[file]] - 
                             1), 1:ncol(t2)], inherits = TRUE)
          index <- which(t2[, 1] == "RSTOPHERE")[[file]] + 
            1
          print(get(sprintf("%s%i", varname, commandcount)))
          commandcount <- commandcount + 1
        }
        else {
          t3 <- rbind(t3, t2[index:(which(t2[, 1] == "RSTOPHERE")[[file]] - 
                                      1), 1:ncol(t2)])
          index <- which(t2[, 1] == "RSTOPHERE")[[file]] + 
            1
          if ((file == length(which(t2[, 1] == "RSTOPHERE")))) {
            commandcount <- commandcount + 1
            remove(t2)
            return(t3)
          }
        }
      }
      remove(t2)
      if ((file > 1) && (separate)) {
        print(sprintf("---------------- Stream Completed: %i files stored in %i commands: %s%i to %s%i ----------------", 
                      file, file, varname, commandcount - file, varname, 
                      commandcount - 1))
        if (grepl("#filelist", line)) {
          linepart1 <- unlist(strsplit(line, split = " "))
          t1 <- pipe(paste(linepart1[[1]], linepart1[[2]], 
                           linepart1[[3]], linepart1[[4]], linepart1[[5]], 
                           linepart1[[6]], "\"aq_pp -f,eok - -d %cols 2> /dev/null | echo %file \"", 
                           sep = " "), open = "r")
          t2 <- read.csv(t1, header = FALSE, sep = ",", 
                         quote = "\"'", comment.char = "#", blank.lines.skip = FALSE, 
                         allowEscapes = TRUE, skip = 0)
          assign(sprintf("%s%i", varname, commandcount), 
                 t2[1:file, 1:ncol(t2)], inherits = TRUE)
          print(get(sprintf("%s%i", varname, commandcount)))
          print(sprintf("---------------- Filenames stored in %s%i ----------------", 
                        varname, commandcount))
          commandcount <- commandcount + 1
          close(t1)
          remove(t1)
          remove(t2)
        }
      }
      remove(colspec, varname, index, separate, t3)
    }
    else {
      if (((substr(essentia, 1, 8) == "ess exec") && (grepl("#Rignore", 
                                                            flags))) || ((substr(essentia, 1, 10) == "ess stream" || 
                                                                          substr(essentia, 1, 9) == "ess query") && (!grepl("#Rinclude", 
                                                                                                                            flags)))) {
        line <- paste(essentia, "\"", aq, "\"", flags, sep = " ")
        system(line)
      }
      else {
        line <- paste(essentia, aq, flags, sep = " ")
        system(line)
      }
    }
    remove(commandcount, line, essentia, aq, flags)
  }
}