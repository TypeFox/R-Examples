capture.essentia <- function(scriptcall = "", linenumber = "all", separator = " ")
{
  argumentlist <- read.delim(header=FALSE,sep=separator,quote="\"'",comment.char="#",text=scriptcall)[1,]
  print(argumentlist)
  tmpcommand1 <- FALSE
  ignoreargs <- FALSE
  if (length(argumentlist)==1) {
    argumentlist <- matrix(argumentlist,1,2)
    argumentlist[,2] <- "ignore"
    ignoreargs <- TRUE
  }
  file <- sprintf("%s",argumentlist[1,1])
  print(file)
  commandcount <- 1
  lineold <- ""
  if (linenumber == "all") {
    lines <- readLines(file)
  }
  else {
    lines <- readLines(file)[linenumber]
  }
  for (line in lines) {
    if (line != "") {
      if (grepl(".R", file)) {
        line <- substr(line, 2, nchar(line))
        if (line == "THIS IS THE END OF MY ESSENTIA COMMANDS") {
          break
        }
      }
      line <- paste(lineold, line, sep = "")
      if (substr(line, nchar(line), nchar(line)) == "\\") {
        lineold <- substr(line, 1, nchar(line) - 1)
        next
      }
      else {
        lineold <- ""
      }
      if ((length(argumentlist) > 1) && (ignoreargs == FALSE)) {
        for (index in seq(2,length(argumentlist))) {
          bashsub_type <- 0
          while (bashsub_type < 2) {
            if (bashsub_type == 0) {
              line <- strsplit(line,sprintf("\\$%i",index - 1))[[1]]
            }
            else {
              line <- strsplit(line,sprintf("\\$\\{%i\\}",index - 1))[[1]]
            }
            if (length(line) > 1) {
              for (i in seq(1,length(line)-1)) {
                print(line[i])
                print(argumentlist[1,index])
                line[i] <- paste0(line[i],sprintf("%s",argumentlist[1,index]))
              }
            }
            line <- paste(line,collapse="")
            bashsub_type <- bashsub_type + 1
          }
        }
      }
      print(line)
      if (((substr(line, 1, 8) == "ess exec") && (!grepl("#Rignore", 
                                                         line))) || ((substr(line, 1, 10) == "ess stream" || 
                                                                      substr(line, 1, 9) == "ess query") && (grepl("#Rinclude", 
                                                                                                                   line)))) {
        quotes <- grepRaw("\"", line, all = TRUE)
        aq <- quotes[[length(quotes)]]
        if ((substr(line, 1, 10) == "ess stream") || (substr(line, 
                                                             1, 8) == "ess exec")) {
          line <- paste(substr(line, 1, aq - 1), "; echo 'RSTOPHERE'", 
                        substr(line, aq, nchar(line)), sep = "")
        }
        else {
          line <- paste(substr(line, 1, aq), "; echo 'RSTOPHERE'", 
                        substr(line, aq + 1, nchar(line)), sep = "")
        }
        colspec <- TRUE
        if (grepl("-notitle", line)) {
          colspec <- FALSE
        }
        varname <- "command"
        if (grepl("#R#", line)) {
          titleindex <- grepRaw("#R#", line, all = TRUE)
          varname <- substr(line, titleindex[[1]] + 3, 
                            titleindex[[2]] - 1)
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
            if (varname == "command") {
              assign(sprintf("%s%i", varname, commandcount), 
                     t2[index:(which(t2[, 1] == "RSTOPHERE")[[file]] - 
                                 1), 1:ncol(t2)], inherits = TRUE)
              index <- which(t2[, 1] == "RSTOPHERE")[[file]] + 
                1
              numrows <- length(get(sprintf("%s%i", varname, commandcount))[,1])
              print(get(sprintf("%s%i", varname, commandcount))[1:min(10,numrows),])
              print(sprintf("---------------- There are a total of %i Rows ----------------", numrows))
              commandcount <- commandcount + 1
            }
            else {
              assign(sprintf("%s%i", varname, file), t2[index:(which(t2[, 
                                                                        1] == "RSTOPHERE")[[file]] - 1), 1:ncol(t2)], 
                     inherits = TRUE)
              index <- which(t2[, 1] == "RSTOPHERE")[[file]] + 
                1
              numrows <- length(get(sprintf("%s%i", varname, file))[,1])
              print(get(sprintf("%s%i", varname, file))[1:min(10,numrows),])
              print(sprintf("---------------- There are a total of %i Rows ----------------", numrows))
              commandcount <- commandcount + 1
            }
          }
          else {
            t3 <- rbind(t3, t2[index:(which(t2[, 1] == 
                                              "RSTOPHERE")[[file]] - 1), 1:ncol(t2)])
            index <- which(t2[, 1] == "RSTOPHERE")[[file]] + 
              1
            if ((file == length(which(t2[, 1] == "RSTOPHERE")))) {
              if (varname == "command") {
                if (commandcount == 1) {
                  tmpcommand1 <- t3
                  numrows <- length(tmpcommand1[,1])
                  print(tmpcommand1[1:min(10,numrows),])
                  print(sprintf("---------------- There are a total of %i Rows ----------------", numrows))
                }
                else {
                  assign(sprintf("command%i", commandcount), 
                         t3, inherits = TRUE)
                  numrows <- length(get(sprintf("%s%i", varname, commandcount))[,1])
                  print(get(sprintf("command%i", commandcount))[1:min(10,numrows),])
                  print(sprintf("---------------- There are a total of %i Rows ----------------", numrows))
                  print(sprintf("---------------- Output Stored in command%i ----------------", 
                                commandcount))
                }
              }
              else {
                assign(sprintf("%s", varname), t3, inherits = TRUE)
                numrows <- length(get(sprintf("%s", varname))[,1])
                print(get(sprintf("%s", varname))[1:min(10,numrows),])
                print(sprintf("---------------- There are a total of %i Rows ----------------", numrows))
                print(sprintf("---------------- Output Stored in %s ----------------", 
                              varname))
              }
              commandcount <- commandcount + 1
            }
          }
        }
        remove(t2)
        if ((file > 1) && (separate)) {
          if (varname == "command") {
            print(sprintf("---------------- Stream Completed: %i files stored in %i commands: command%i to command%i ----------------", 
                          file, file, commandcount - file, commandcount - 
                            1))
            if (grepl("#filelist", line)) {
              linepart1 <- unlist(strsplit(line, split = " "))
              t1 <- pipe(paste(linepart1[[1]], linepart1[[2]], 
                               linepart1[[3]], linepart1[[4]], linepart1[[5]], 
                               linepart1[[6]], "\"aq_pp -f,eok - -d %cols 2> /dev/null | echo %file \"", 
                               sep = " "), open = "r")
              t2 <- read.csv(t1, header = FALSE, sep = ",", 
                             quote = "\"'", comment.char = "#", blank.lines.skip = FALSE, 
                             allowEscapes = TRUE, skip = 0)
              assign(sprintf("command%i", commandcount), 
                     t2[1:file, 1:ncol(t2)], inherits = TRUE)
              numrows <- length(get(sprintf("command%i", commandcount))[,1])
              print(get(sprintf("command%i", commandcount))[1:min(10,numrows),])
              print(sprintf("---------------- There are a total of %i Files ----------------", numrows))
              print(sprintf("---------------- Filenames stored in command%i ----------------", 
                            commandcount))
              commandcount <- commandcount + 1
              close(t1)
              remove(t1)
              remove(t2)
            }
          }
          else {
            print(sprintf("---------------- Stream Completed: %i files stored in %i commands: %s%i to %s%i ----------------", 
                          file, file, varname, 1, varname, file))
            if (grepl("#filelist", line)) {
              linepart1 <- unlist(strsplit(line, split = " "))
              t1 <- pipe(paste(linepart1[[1]], linepart1[[2]], 
                               linepart1[[3]], linepart1[[4]], linepart1[[5]], 
                               linepart1[[6]], "\"aq_pp -f,eok - -d %cols 2> /dev/null | echo %file \"", 
                               sep = " "), open = "r")
              t2 <- read.csv(t1, header = FALSE, sep = ",", 
                             quote = "\"'", comment.char = "#", blank.lines.skip = FALSE, 
                             allowEscapes = TRUE, skip = 0)
              assign(sprintf("%s%i", varname, file + 1), 
                     t2[1:file, 1:ncol(t2)], inherits = TRUE)
              numrows <- length(get(sprintf("%s%i", varname, file + 1))[,1])
              print(get(sprintf("%s%i", varname, file + 
                                  1))[1:min(10,numrows),])
              print(sprintf("---------------- There are a total of %i Files ----------------", numrows))
              print(sprintf("---------------- Filenames stored in %s%i ----------------", 
                            varname, file + 1))
              close(t1)
              remove(t1)
              remove(t2)
            }
          }
        }
      }
    else {
      if (line != "") {
        system(line)
      }
    }
    }
  }
  if ((commandcount == 2) && (tmpcommand1 != FALSE)) {
    return(tmpcommand1)
  }
  else {
    if ((commandcount > 2) && (tmpcommand1 != FALSE)) {
    assign("command1",tmpcommand1,inherits=TRUE)
    print("---------------- First Command's Output Stored in command1 ----------------")
    }
  }
  print(sprintf("---------------- There are a total of %i commands ----------------", 
                commandcount - 1))
  remove(colspec, commandcount, index, line, lineold, lines, 
         separate, t3, varname, tmpcommand1)
}