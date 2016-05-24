# $Id: NMBU.GUI.File.R 35 2014-01-10 21:17:26Z khliland $

##
## GUI functions for the File menu
##


#    item	fileMenu			command			  "Save project ..."						      				saveAll						""		  					""
#    item	fileMenu			command			  "Save project as ..."				    						saveAllAs					""		  					""
#    item	fileMenu			command			  "Open project ..."      										loadAll						""		  					""
#    item    fileMenu        	separator         ""                                              				""                          ""          				""
#####################################
# Safe log, output and data
saveAllWorker <- function(allFile){ # Denne fungerer ikke med nyere R Commander
  # Collect, compress and save files
  wd <- getwd()
  dir.create(paste(wd,"/tmpRCE",sep=""))
  setwd(paste(wd,"/tmpRCE",sep=""))
  tmp.script <- paste(wd,"/tmpRCE/script.R",sep="")
  tmp.output <- paste(wd,"/tmpRCE/output.R",sep="")
  tmp.data   <- paste(wd,"/tmpRCE/data.RData",sep="")
  tmp.Rcmdr  <- paste(wd,"/tmpRCE/Rcmdr.RData",sep="")
  
  # Log-file
  log <- tclvalue(tkget(LogWindow(), "1.0", "end"))
  fileCon <- file(tmp.script, "w")
  cat(log, file = fileCon)
  close(fileCon)

  # Output-file
  output <- tclvalue(tkget(OutputWindow(), "1.0", "end"))
  fileCon <- file(tmp.output, "w")
  cat(output, file = fileCon)
  close(fileCon)
  
  # Data-file
  save.image(file = tmp.data)
  
  # Rcmdr variables
  pos <-  match("RcmdrEnv", search())
  vars <- mget(c("outputStack","commandStack","length.command.stack",
                 "dialog.values","length.output.stack","modelNumber"), pos.to.env(pos))
  active <- list(DataSet=ActiveDataSet(), Model=ActiveModel())
  save("vars","active", file="Rcmdr.RData")
  
  # Compress files, save and clean up
  zip(paste(wd,"/tmpRCE/tmp.zip",sep=""),c("script.R","output.R","data.RData","Rcmdr.RData"))
  file.copy(paste(wd,"/tmpRCE/tmp.zip", sep=""), allFile, overwrite=TRUE)
  setwd(wd)
  unlink("tmpRCE", recursive = TRUE)

  putRcmdr("allFileName", allFile)
  Message(paste(gettextRcmdr("R Commander environment saved to"), allFile), type="note")
}

saveAllAs <- function() {
  allFile <- tclvalue(tkgetSaveFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"R Commander project" {".RCP"}}'),
                                    defaultextension="RCP", initialfile="RCommander.RCP"))
  if (allFile == "") return()
  saveAllWorker(allFile)
}

saveAll <- function() {
  .allFileName <- getRcmdr("allFileName")
  if (is.null(.allFileName)) {
    saveAllAs()
    return()
  }
  saveAllWorker(.allFileName)
}

#####################################
# Load log, output and data
loadAll <- function(){
  allFile <- tclvalue(tkgetOpenFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"R Commander project" {".RCP"}}'),
                                    defaultextension="RCP"))
  if (allFile == "") return()
  
  wd <- getwd()
  dir.create(paste(wd,"/tmpRCE",sep=""))
  setwd(paste(wd,"/tmpRCE",sep=""))
  tmp.script <- paste(wd,"/tmpRCE/script.R",sep="")
  tmp.output <- paste(wd,"/tmpRCE/output.R",sep="")
  tmp.data   <- paste(wd,"/tmpRCE/data.RData",sep="")
  tmp.Rcmdr  <- paste(wd,"/tmpRCE/Rcmdr.RData",sep="")
  
  unzip(allFile, junkpaths = FALSE, exdir = ".", unzip = "internal", setTimes = FALSE)
  
  # Restore log (script)
  fileCon <- file(tmp.script, "r")
  contents <- readLines(fileCon)
  close(fileCon)
  .log <- LogWindow()
  if (tclvalue(tkget(.log, "1.0", "end")) != "\n"){
    response2 <- RcmdrTkmessageBox(message=gettextRcmdr("Replace current script with saved script"),
                                   icon="question", type="yesno", default="yes")
    if ("no" == tclvalue(response2)) return()
  }
  tkdelete(.log, "1.0", "end")
  tkinsert(.log, "end", paste(contents, collapse="\n"))
  
  # Restore output
  fileCon <- file(tmp.output, "r")
  contents <- readLines(fileCon)
  close(fileCon)
  .output <- OutputWindow()
  tkdelete(.output, "1.0", "end")
  for (line in contents) {
    red <- ifelse(grepl("^>",line),TRUE, grepl("^[\\+]",line))
    tkinsert(.output, "end", paste(line, "\n", sep = ""))
    if(red){
      tktag.add(.output, "currentLine", "end - 2 lines linestart", 
                "end - 2 lines lineend")
      tktag.configure(.output, "currentLine", foreground = getRcmdr("command.text.color"))
    }
    tkyview.moveto(.output, 1)
  }  
  
  # Restore data
  load(tmp.data,.GlobalEnv)
  
  # Restore Rcmdr variables
  .tmp.env <- new.env()
  load(tmp.Rcmdr, envir = .tmp.env)
  pos <-  match("RcmdrEnv", search())
  vars <- get("vars",envir=.tmp.env)
  putRcmdr("outputStack",vars$outputStack)
  putRcmdr("commandStack",vars$commandStack)
  putRcmdr("length.command.stack",vars$length.command.stack)
  putRcmdr("dialog.values",vars$dialog.values)
  putRcmdr("length.output.stack",vars$length.output.stack)
  putRcmdr("modelNumber",vars$modelNumber)
  active <- get("active",envir=.tmp.env)
  activeDataSet(active$DataSet)
  activeModel(active$Model)
  rm(.tmp.env)

  setwd(wd)
  unlink("tmpRCE", recursive = TRUE)
}
