# last modified 2016-03-20 by J. Fox

# File (and Edit) menu dialogs

loadLog <- function(){
	logFile <- tclvalue(tkgetOpenFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"Script Files" {".R"}}'),
					defaultextension="R",
					parent=CommanderWindow()))
	if (logFile == "") return()
	fileCon <- file(logFile, "r")
	contents <- readLines(fileCon)
	close(fileCon)
	currentLogFileName <- getRcmdr("logFileName")
	putRcmdr("logFileName", logFile)
	.log <- LogWindow()
	if (tclvalue(tkget(.log, "1.0", "end")) != "\n"){
		response2 <- RcmdrTkmessageBox(message=gettextRcmdr("Save current log file?"),
				icon="question", type="yesno", default="yes")
		if ("yes" == tclvalue(response2)) saveLog(currentLogFileName)
	}
	tkdelete(.log, "1.0", "end")
	tkinsert(.log, "end", paste(contents, collapse="\n"))
}

saveLog <- function(logfilename) {
  .logFileName <- if (missing(logfilename) || (logfilename == "%logfilename")) 
    getRcmdr("logFileName") else logfilename
  if (is.null(.logFileName)) {
    saveLogAs()
    return()
  }
  log <- tclvalue(tkget(LogWindow(), "1.0", "end"))
  fileCon <- file(.logFileName, "w")
  cat(log, file = fileCon)
  close(fileCon)
  Message(paste(gettextRcmdr("Script saved to"), .logFileName), type="note")
}

saveLogAs <- function() {
	logFile <- tclvalue(tkgetSaveFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"Script Files" {".R"}}'),
					defaultextension="R",
					initialfile="RCommander.R",
					parent=CommanderWindow()))
	if (logFile == "") return()
	log <- tclvalue(tkget(LogWindow(), "1.0", "end"))
	fileCon <- file(logFile, "w")
	cat(log, file = fileCon)
	close(fileCon)
	putRcmdr("logFileName", logFile)
	Message(paste(gettextRcmdr("Script saved to"), logFile), type="note")
}

loadRmd <- function(){
    RmdFile <- tclvalue(tkgetOpenFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"R Markdown Files" {".Rmd" ".rmd"}}'),
        defaultextension="Rmd",
        parent=CommanderWindow()))
    if (RmdFile == "") return()
    fileCon <- file(RmdFile, "r")
    contents <- readLines(fileCon)
    close(fileCon)
    currentRmdFileName <- getRcmdr("RmdFileName")
    putRcmdr("RmdFileName", RmdFile)
    .rmd <- RmdWindow()
    if (tclvalue(tkget(.rmd, "1.0", "end")) != "\n"){
        response2 <- RcmdrTkmessageBox(message=gettextRcmdr("Save current Rmd file?"),
            icon="question", type="yesno", default="yes")
        if ("yes" == tclvalue(response2)) saveLog(currentRmdFileName)
    }
    tkdelete(.rmd, "1.0", "end")
    tkinsert(.rmd, "end", paste(contents, collapse="\n"))
}

loadRnw <- function(){
    RnwFile <- tclvalue(tkgetOpenFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"knitr Files" {".Rnw" ".rnw" ".Snw" ".snw"}}'),
                                      defaultextension="Rnw",
                                      parent=CommanderWindow()))
    if (RnwFile == "") return()
    fileCon <- file(RnwFile, "r")
    contents <- readLines(fileCon)
    close(fileCon)
    currentRnwFileName <- getRcmdr("RnwFileName")
    putRcmdr("RnwFileName", RnwFile)
    .rnw <- RnwWindow()
    if (tclvalue(tkget(.rnw, "1.0", "end")) != "\n"){
        response2 <- RcmdrTkmessageBox(message=gettextRcmdr("Save current Rnw file?"),
                                       icon="question", type="yesno", default="yes")
        if ("yes" == tclvalue(response2)) saveLog(currentRnwFileName)
    }
    tkdelete(.rnw, "1.0", "end")
    contents <- paste(contents, collapse="\n")
    contents <- sub("\n\\\\end\\{document\\}\n", "", contents)
    tkinsert(.rnw, "end", contents)
}

saveRmd <- function(Rmdfilename) {
    .RmdFileName <- if (missing(Rmdfilename) || (Rmdfilename == "%Rmdfilename")) 
                        getRcmdr("RmdFileName") else Rmdfilename
    if ((.RmdFileName == "RcmdrMarkdown.Rmd") || (.RmdFileName == "RcmdrRMarkdown.Rmd") || is.null(.RmdFileName)) {
        saveRmdAs()
        return()
    }
    .rmd <- tclvalue(tkget(RmdWindow(), "1.0", "end"))
    fileCon <- file(.RmdFileName, "w")
    cat(.rmd, file = fileCon)
    close(fileCon)
    Message(paste(gettextRcmdr("R Markdown file saved to"), .RmdFileName), type="note")
}

saveRnw <- function(Rnwfilename) {
    .RnwFileName <- if (missing(Rnwfilename) || (Rnwfilename == "%Rnwfilename")) 
      getRcmdr("RnwFileName") else Rnwfilename
    if ((.RnwFileName == "RcmdrKnitr.Rnw") || is.null(.RnwFileName)) {
        saveRnwAs()
        return()
    }
    .rnw <- tclvalue(tkget(RnwWindow(), "1.0", "end"))
    .rnw <- paste(.rnw, "\n\\end{document}\n")
    fileCon <- file(.RnwFileName, "w")
    cat(.rnw, file = fileCon)
    close(fileCon)
    Message(paste(gettextRcmdr("knitr file saved to"), .RnwFileName), type="note")
}


saveRmdAs <- function() {
    RmdFile <- tclvalue(tkgetSaveFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"R Markdown Files" {".Rmd" ".rmd"}}'),
        defaultextension="Rmd",
        initialfile="RCommanderMarkdown.Rmd",
        parent=CommanderWindow()))
    if (RmdFile == "") return()
    .rmd <- tclvalue(tkget(RmdWindow(), "1.0", "end"))
    fileCon <- file(RmdFile, "w")
    cat(.rmd, file = fileCon)
    close(fileCon)
    putRcmdr("RmdFileName", RmdFile)
    Message(paste(gettextRcmdr("R Markdown file saved to"), RmdFile), type="note")
}

saveRnwAs <- function() {
    RnwFile <- tclvalue(tkgetSaveFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"knitr Files" {".Rnw" ".rnw" ".Snw" ".snw"}}'),
                                      defaultextension="Rnw",
                                      initialfile="RCommanderKnitr.Rnw",
                                      parent=CommanderWindow()))
    if (RnwFile == "") return()
    .rnw <- tclvalue(tkget(RnwWindow(), "1.0", "end"))
    .rnw <- paste(.rnw, "\n\\end{document}\n")
    fileCon <- file(RnwFile, "w")
    cat(.rnw, file = fileCon)
    close(fileCon)
    putRcmdr("RnwFileName", RnwFile)
    Message(paste(gettextRcmdr("knitr file saved to"), RnwFile), type="note")
}


saveOutput <- function() {
	.outputFileName <- getRcmdr("outputFileName")
	if (is.null(.outputFileName)) {
		saveOutputAs()
		return()
	}
	output <- tclvalue(tkget(OutputWindow(), "1.0", "end"))
	fileCon <- file(.outputFileName, "w")
	cat(output, file = fileCon)
	close(fileCon)
	Message(paste(gettextRcmdr("Output saved to"), .outputFileName), type="note")
}

saveOutputAs <- function() {
	outputFile <- tclvalue(tkgetSaveFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"Output Files" {".txt"}}'),
					defaultextension="txt",
					initialfile="RCommander.txt",
					parent=CommanderWindow()))
	if (outputFile == "") return()
	output <- tclvalue(tkget(OutputWindow(), "1.0", "end"))
	fileCon <- file(outputFile, "w")
	cat(output, file = fileCon)
	close(fileCon)
	putRcmdr("outputFileName", outputFile)
	Message(paste(gettextRcmdr("Output saved to"), outputFile), type="note")
}

saveWorkspaceAs <- function(){
	saveFile <- tclvalue(tkgetSaveFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"R Data Files" {".RData" ".rda" ".Rda" ".RDA"}}'),
					defaultextension="",
					initialfile=".RData",
					parent=CommanderWindow()))
	if (saveFile == "") return()
	save(list=ls(envir=.GlobalEnv), file=saveFile)
	putRcmdr("saveFileName", saveFile)
	Message(paste(gettextRcmdr("R workspace saved to"), saveFile), type="note")
}

saveWorkspace <- function() {
	.saveFileName <- getRcmdr("saveFileName")
	if (is.null(.saveFileName)) {
		saveWorkspaceAs()
		return()
	}
	else save(list=ls(envir=.GlobalEnv), file=.saveFileName)
	Message(paste(gettextRcmdr("R workspace saved to"), .saveFileName), type="note")
}

CloseCommander <- function() closeCommander(ask=getRcmdr("ask.to.exit"), ask.save=getRcmdr("ask.on.exit"))

closeCommander <- function(ask=TRUE, ask.save=ask){
	if (ask){
		response <- tclvalue(RcmdrTkmessageBox(message=gettextRcmdr("Exit?"),
						icon="question", type="okcancel", default="cancel"))
		if (response == "cancel") return(invisible(response))
	}
	else {
		ask.save=FALSE
		response <- "ok"
	}
	sink(type="message")
	if (!is.null(ActiveDataSet()) && getRcmdr("attach.data.set"))
		justDoIt(logger(paste("detach(", ActiveDataSet(), ")", sep="")))
	putRcmdr(".activeDataSet", NULL)
	putRcmdr(".activeModel", NULL)
	if (ask.save && getRcmdr("log.commands") && tclvalue(tkget(LogWindow(), "1.0", "end")) != "\n"){
		response2 <- RcmdrTkmessageBox(message=gettextRcmdr("Save script file?"),
				icon="question", type="yesno", default="yes")
		if ("yes" == tclvalue(response2)) saveLog()
	}
    
	if (ask.save && getRcmdr("markdown.output") && getRcmdr("log.commands") && tclvalue(tkget(RmdWindow(), "1.0", "end")) != "\n"){
	    response2 <- RcmdrTkmessageBox(message=gettextRcmdr("Save R Markdown file?"),
	                                   icon="question", type="yesno", default="yes")
	    if ("yes" == tclvalue(response2)) saveRmd()
	}
    
	if (ask.save && getRcmdr("knitr.output") && getRcmdr("log.commands") && tclvalue(tkget(RnwWindow(), "1.0", "end")) != "\n"){
	    response2 <- RcmdrTkmessageBox(message=gettextRcmdr("Save knitr file?"),
	                                   icon="question", type="yesno", default="yes")
	    if ("yes" == tclvalue(response2)) saveRnw()
	}
    
	if (ask.save && !getRcmdr("console.output") && tclvalue(tkget(OutputWindow(), "1.0", "end")) != "\n"){
		response3 <- RcmdrTkmessageBox(message=gettextRcmdr("Save output file?"),
				icon="question", type="yesno", default="yes")
		if ("yes" == tclvalue(response3)) saveOutput()
	}
  if (MacOSXP()){
    Sys.setenv(PATH=getRcmdr("PATH"))
  }
	if (!WindowsP()) options(getRcmdr("oldPager"))
	if (getRcmdr("suppress.X11.warnings")) {
		sink(type = "message")
		close(getRcmdr("messages.connection"))
	}
	options(getRcmdr("saveOptions"))
  options(help_type = getRcmdr("restore.help_type"))
 	options(device = getRcmdr("restore.device"))
#   if (getRcmdr("restore.use.external.help")) 
#     system("defaults delete org.R-project.R use.external.help")
	tkdestroy(CommanderWindow())
	putRcmdr("commanderWindow", NULL)
	putRcmdr("logWindow", NULL)
    putRcmdr("RmdWindow", NULL)
	putRcmdr("messagesWindow", NULL)
	putRcmdr("outputWindow", NULL)
	open.showData.windows <- getRcmdr("open.showData.windows")
	if (length(open.showData.windows) > 0){
	  for (window in open.showData.windows){
	    tkdestroy(window)
	  }
	  putRcmdr("open.showData.windows", list())
	}
	options(getRcmdr("quotes"))
	tkwait <- options("Rcmdr")[[1]]$tkwait  # to address problem in Debian Linux
	if ((!is.null(tkwait)) && tkwait) putRcmdr(".commander.done", tclVar("1"))
	return(invisible(response))
}

closeCommanderAndR <- function(){
	response <- CloseCommander()
	if (response == "cancel") return()
	cat("\n")
	quit(save="no")
}

Options <- function(){
  setOption <- function(option, default) {
    if (!is.null(current[[option]])) return(current[[option]])
    else if (!is.null(getRcmdr(option, fail=FALSE))) return(getRcmdr(option))
    return(default)
  }
  asLogical <- function(x) as.logical(as.numeric(x))
  initializeDialog(title=gettextRcmdr("Commander Options"))
  notebook <- ttknotebook(top)
  closeTab <- tkframe(top)
  fontTab <- tkframe(top)
  outputTab <- tkframe(top)
  otherTab <- tkframe(top)
  current <- getOption("Rcmdr")
  console.output <- getRcmdr("console.output")
  default.font.size <- getRcmdr("default.font.size")
  default.font.family <- getRcmdr("default.font.family")
  log.commands <- getRcmdr("log.commands")
  log.font.size <- getRcmdr("log.font.size")
  log.font.family <- getRcmdr("log.font.family")
  log.width <- setOption("log.width", 80)
  log.height <- if (!is.null(current$log.height)) current$log.height
  else if (!log.commands) 0 else 10
  output.height <- if (!is.null(current$output.height)) current$output.height
  else if (console.output) 0 else 2*log.height
  contrasts <- setOption("default.contrasts", c("contr.Treatment", "contr.poly"))
  grab.focus <- getRcmdr("grab.focus")
  double.click <- getRcmdr("double.click")
  sort.names <- getRcmdr("sort.names")
  show.edit.button <- setOption("show.edit.button", TRUE)
  scale.factor <- current$scale.factor
  suppress.icon.images <- getRcmdr("suppress.icon.images")
  number.messages <- getRcmdr("number.messages")
  etc <- getRcmdr("etc")
  etcMenus <- getRcmdr("etcMenus")
  log.font <- tclvalue(tkfont.actual("RcmdrLogFont"))
  log.font.family <- tclvalue(.Tcl("font actual RcmdrLogFont -family"))
  if (length(grep(" ", log.font.family)) > 1) log.font.family <- paste("{", log.font.family, "}", sep="")
  title.color <- getRcmdr("title.color")
  # On Windows 7 with Classic theme, a color name is returned instead of RGB values
  if(substr(title.color, 1, 1) != "#") {
    title.color.rgb <- as.numeric(tkwinfo("rgb", top, title.color))
    title.color <- rgb(rbind(title.color.rgb), maxColorValue=65535)
  }
  use.markdown<- getRcmdr("use.markdown")
  use.knitr<- getRcmdr("use.knitr")
  retain.selections <- getRcmdr("retain.selections")
  messages.height <- as.character(getRcmdr("messages.height"))
  ask.to.exit <- getRcmdr("ask.to.exit")
  ask.on.exit <- getRcmdr("ask.on.exit")
  attach.data.set <- getRcmdr("attach.data.set")
  log.text.color <- getRcmdr("log.text.color")
  command.text.color <- getRcmdr("command.text.color")
  output.text.color <- getRcmdr("output.text.color")
  error.text.color <- getRcmdr("error.text.color")
  warning.text.color <- getRcmdr("warning.text.color")
  prefixes <- getRcmdr("prefixes")
  multiple.select.mode <- getRcmdr("multiple.select.mode")
  suppress.X11.warnings <- getRcmdr("suppress.X11.warnings")
  showData.threshold <- getRcmdr("showData.threshold")
  retain.messages <- getRcmdr("retain.messages")
  crisp.dialogs <- getRcmdr("crisp.dialogs")
  length.output.stack <- getRcmdr("length.output.stack")
  length.command.stack <- getRcmdr("length.command.stack")
  quit.R.on.close <- getRcmdr("quit.R.on.close")
  variable.list.height <- getRcmdr("variable.list.height")
  variable.list.width <- getRcmdr("variable.list.width")
  placement <- setOption("placement", "")
  suppress.menus <- getRcmdr("suppress.menus")
  rmd.template <- getRcmdr("rmd.template")
  rmd.standard <- system.file("etc", "Rcmdr-Markdown-Template.Rmd", package="Rcmdr")
  rnw.template <- getRcmdr("rnw.template")
  rnw.standard <- system.file("etc", "Rcmdr-knitr-Template.Rnw", package="Rcmdr")
  use.rgl <- setOption("use.rgl", TRUE)
  checkBoxes(closeTab, frame="closeOptionsFrame", boxes=c("askToExit", "askOnExit", "quitR"),
             initialValues=c(ask.to.exit, ask.on.exit, quit.R.on.close),
             labels=gettextRcmdr("Ask to exit Commander", "Ask to save documents on exit", "Quit R on exit"))
  checkBoxes(outputTab, frame="outputOptionsFrame", 
             boxes=c("consoleOutput", "logCommands", "numberMessages", "retainMessages", "useMarkdown", "useKnitr"),
             initialValues=c(console.output, log.commands, number.messages, retain.messages, use.markdown, use.knitr),
             labels=gettextRcmdr("Send output to R Console", "Log commands to script window", "Number messages", 
                                 "Retain messages", "Use R Markown", "Use knitr"))
  cval <- function(x,y) -sum((x-y)^2)
  contrasting <- function(x)
    optim(rep(127, 3),cval,lower=0,upper=255,method="L-BFGS-B",y=x)$par
  # the following local function from Thomas Lumley via r-help
  convert <- function (color){
    rgb <- col2rgb(color)/255
    L <- c(0.2, 0.6, 0) %*% rgb
    ifelse(L >= 0.2, "#000060", "#FFFFA0")
  }
  env <- environment()
  pal <- c(log.text.color, command.text.color, output.text.color, error.text.color, warning.text.color, title.color)
  pickColor <- function(initialcolor, parent){
    newcolor <- tclvalue(.Tcl(paste("tk_chooseColor", .Tcl.args(title = "Select a Color",
                                                                initialcolor=initialcolor, parent=parent))))
    if (newcolor == "") initialcolor else newcolor
  }
  hexcolor <- colorConverter(toXYZ = function(hex,...) {
    rgb <- t(col2rgb(hex))/255
    colorspaces$sRGB$toXYZ(rgb,...) },
    fromXYZ = function(xyz,...) {
      rgb <- colorspaces$sRGB$fromXYZ(xyz,..)
      rgb <- round(rgb,5)
      if (min(rgb) < 0 || max(rgb) > 1) as.character(NA)
      else rgb(rgb[1],rgb[2],rgb[3])},
    white = "D65", name = "#rrggbb")
  cols <- t(col2rgb(pal))
  hex <- convertColor(cols, from="sRGB", to=hexcolor, scale.in=255, scale.out=NULL)
  for (i in 1:8) assign(paste("hex", i, sep="."), hex[i], envir=env)
  fontColorsFrame <- tkframe(fontTab)
  colorField1 <- labelRcmdr(fontColorsFrame, text=rgb2col(hex[1]), fg=hex[1])
  button1 <- tkbutton(fontColorsFrame, text=hex[1], bg = hex[1], width="10",
                      fg=convert(hex[1]),
                      command=function() {
                        color <- pickColor(hex[1], parent=button1)
                        fg <- convert(color)
                        tkconfigure(button1, bg=color, fg=fg, text=toupper(color))
                        tkconfigure(colorField1, text=rgb2col(color), foreground=color)
                        assign("hex.1", color, envir=env)
                      }
  )
  colorField2 <- labelRcmdr(fontColorsFrame, text=rgb2col(hex[2]), fg=hex[2])
  button2 <- tkbutton(fontColorsFrame, text=hex[2], bg = hex[2], width="10",
                      fg=convert(hex[2]),
                      command=function() {
                        color <- pickColor(hex[2], parent=button2)
                        fg <- convert(color)
                        tkconfigure(button2, bg=color, fg=fg, text=toupper(color))
                        tkconfigure(colorField2, text=rgb2col(color), foreground=color)
                        assign("hex.2", color, envir=env)
                      }
  )
  colorField3 <- labelRcmdr(fontColorsFrame, text=rgb2col(hex[3]), fg=hex[3])
  button3 <- tkbutton(fontColorsFrame, text=hex[3], bg = hex[3], width="10",
                      fg=convert(hex[3]),
                      command=function() {
                        color <- pickColor(hex[3], parent=button3)
                        fg <- convert(color)
                        tkconfigure(button3, bg=color, fg=fg, text=toupper(color))
                        tkconfigure(colorField3, text=rgb2col(color), foreground=color)
                        assign("hex.3", color, envir=env)
                      }
  )
  colorField4 <- labelRcmdr(fontColorsFrame, text=rgb2col(hex[4]), fg=hex[4])
  button4 <- tkbutton(fontColorsFrame, text=hex[4], bg = hex[4], width="10",
                      fg=convert(hex[4]),
                      command=function() {
                        color <- pickColor(hex[4], parent=button4)
                        fg <- convert(color)
                        tkconfigure(button4, bg=color, fg=fg, text=toupper(color))
                        tkconfigure(colorField4, text=rgb2col(color), foreground=color)
                        assign("hex.4", color, envir=env)
                      }
  )
  colorField5 <- labelRcmdr(fontColorsFrame, text=rgb2col(hex[5]), fg=hex[5])
  button5 <- tkbutton(fontColorsFrame, text=hex[5], bg = hex[5], width="10",
                      fg=convert(hex[5]),
                      command=function() {
                        color <- pickColor(hex[5], parent=button5)
                        fg <- convert(color)
                        tkconfigure(button5, bg=color, fg=fg, text=toupper(color))
                        tkconfigure(colorField5, text=rgb2col(color), foreground=color)
                        assign("hex.5", color, envir=env)
                      }
  )
  colorField6 <- labelRcmdr(fontColorsFrame, text=rgb2col(hex[6]), fg=hex[6])
  button6 <- tkbutton(fontColorsFrame, text=hex[6], bg = hex[6], width="10",
                      fg=convert(hex[6]),
                      command=function() {
                        color <- pickColor(hex[6], parent=button6)
                        fg <- convert(color)
                        tkconfigure(button6, bg=color, fg=fg, text=toupper(color))
                        tkconfigure(colorField6, text=rgb2col(color), foreground=color)
                        assign("hex.6", color, envir=env)
                      }
  )
  logFontSizeVar <- tclVar(log.font.size)
  fontFrame <- tkframe(fontTab)
  logFontSizeSlider <- tkscale(fontFrame, from=6, to=20, showvalue=TRUE, variable=logFontSizeVar,
                               resolution=1, orient="horizontal")
  logWidthVar <- tclVar(log.width)
  outputSliderFrame <- tkframe(outputTab)
  logWidthSlider <- tkscale(outputSliderFrame, from=30, to=120, showvalue=TRUE, variable=logWidthVar,
                            resolution=5, orient="horizontal")
  logHeightVar <- tclVar(log.height)
  logHeightSlider <- tkscale(outputSliderFrame, from=0, to=25, showvalue=TRUE, variable=logHeightVar,
                             resolution=1, orient="horizontal")
  outputHeightVar <- tclVar(output.height)
  outputHeightSlider <- tkscale(outputSliderFrame, from=0, to=50, showvalue=TRUE, variable=outputHeightVar,
                                resolution=5, orient="horizontal")
  messagesHeightVar <- tclVar(messages.height)
  messagesHeightSlider <- tkscale(outputSliderFrame, from=0, to=10, showvalue=TRUE, variable=messagesHeightVar,
                                  resolution=1, orient="horizontal")       
  contrasts1 <- tclVar(contrasts[1])
  contrasts2 <- tclVar(contrasts[2])
  contrastsFrame <- tkframe(otherTab)
  contrasts1Entry <- ttkentry(contrastsFrame, width="15", textvariable=contrasts1)
  contrasts2Entry <- ttkentry(contrastsFrame, width="15", textvariable=contrasts2)
  checkBoxes(otherTab, frame="otherOptionsFrame", 
             boxes=c("grabFocus", "doubleClick", "sortNames", "showEditButton", "SuppressIconImages",
                     "retainSelections", "useRgl"),
             initialValues=c(grab.focus, double.click, sort.names, show.edit.button, suppress.icon.images,
                             retain.selections, use.rgl),
             labels=gettextRcmdr("Active window grabs focus", "Double-click presses OK button", 
                                 "Sort variable names alphabetically", "Show edit button",
                                 "Suppress icon images", "Retain dialog selections", "Use rgl package")
  )
  scaleFactorFrame <- tkframe(otherTab)
  scaleFactorVar <- tclVar(if (is.null(scale.factor)) 1.0 else scale.factor)
  scaleFactorSlider <- tkscale(scaleFactorFrame, from=0.2, to=3.0, showvalue=TRUE, variable=scaleFactorVar,
                               resolution=0.2, orient="horizontal")
  defaultFontSizeVar <- tclVar(default.font.size)
  defaultFontSizeSlider <- tkscale(fontFrame, from=6, to=20, showvalue=TRUE, variable=defaultFontSizeVar,
                                   resolution=1, orient="horizontal")
  logFontFamilyVar <- tclVar(log.font.family)
  defaultFontFamilyVar <- tclVar(default.font.family)
  logFontEntry <- ttkentry(fontFrame, width="20", textvariable=logFontFamilyVar)
  defaultFontEntry <- ttkentry(fontFrame, width="20", textvariable=defaultFontFamilyVar)
  rmdTemplateVar <- tclVar(rmd.template)
  templateFrame <- tkframe(outputTab)
  rmdTemplateEntry <- ttkentry(templateFrame, width="75", textvariable=rmdTemplateVar)
  onSelectTemplate <- function(){
    templateFile <- tclvalue(tkgetOpenFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"R Markdown Files" {".Rmd" ".rmd"}}'),
                                           defaultextension="Rmd",
                                           parent=outputTab))
    if (templateFile == "") return()
    tclvalue(rmdTemplateVar) <- templateFile
    return(NULL)
  }
  templateButton <- buttonRcmdr(templateFrame, text=gettextRcmdr("Select file"), command=onSelectTemplate)
  rnwTemplateVar <- tclVar(rnw.template)
  rnwTemplateEntry <- ttkentry(templateFrame, width="75", textvariable=rnwTemplateVar)
  onSelectRnwTemplate <- function(){
    rnwTemplateFile <- tclvalue(tkgetOpenFile(filetypes=gettextRcmdr('{"All Files" {"*"}} {"knitr Files" {".Rnw" ".rnw" ".Snw" ".snw"}}'),
                                              defaultextension="Rnw",
                                              parent=outputTab))
    if (rnwTemplateFile == "") return()
    tclvalue(rnwTemplateVar) <- rnwTemplateFile
    return(NULL)
  }
  rnwTemplateButton <- buttonRcmdr(templateFrame, text=gettextRcmdr("Select file"), command=onSelectRnwTemplate)
  all.themes <- tk2theme.list()
  current.theme <- tk2theme()
  all.themes <- union(current.theme, all.themes)
  themesBox <- variableListBox(otherTab, all.themes, 
                               title = gettextRcmdr("Theme (pick one)"), 
                               initialSelection = varPosn(current.theme, vars=all.themes))
  onOK <- function(){
    theme <- getSelection(themesBox)
    closeDialog(top)
    ask.to.exit <- asLogical(tclvalue(askToExitVariable))
    ask.on.exit <- asLogical(tclvalue(askOnExitVariable))
    quit.R.on.close <- asLogical(tclvalue(quitRVariable))
    console.output <- asLogical(tclvalue(consoleOutputVariable))
    number.messages <- asLogical(tclvalue(numberMessagesVariable))
    retain.messages <- asLogical(tclvalue(retainMessagesVariable))
    use.markdown <- asLogical(tclvalue(useMarkdownVariable))
    use.knitr<- asLogical(tclvalue(useKnitrVariable))
    rmd.template <- tclvalue(rmdTemplateVar)
    if (rmd.template == rmd.standard) rmd.template <- NULL
    rnw.template <- tclvalue(rnwTemplateVar)
    if (rnw.template == rnw.standard) rnw.template <- NULL
    log.font.family <- tclvalue(logFontFamilyVar)
    default.font.family <- tclvalue(defaultFontFamilyVar)
    log.font.size <- round(as.numeric(tclvalue(logFontSizeVar)))
    default.font.size <- tclvalue(defaultFontSizeVar)
    scale.factor <- round(as.numeric(tclvalue(scaleFactorVar)), 1)
    if (scale.factor == 1) scale.factor <- NULL
    log.width <- round(as.numeric(tclvalue(logWidthVar)))
    log.height <- as.numeric(tclvalue(logHeightVar))
    log.commands <- asLogical(tclvalue(logCommandsVariable)) && (log.height != 0)
    output.height <- as.numeric(tclvalue(outputHeightVar))
    console.output <- asLogical(tclvalue(consoleOutputVariable)) || (output.height == 0)
    contrasts <- c(tclvalue(contrasts1), tclvalue(contrasts2))
    grab.focus <- asLogical(tclvalue(grabFocusVariable))
    double.click <- asLogical(tclvalue(doubleClickVariable))
    sort.names <- asLogical(tclvalue(sortNamesVariable))
    show.edit.button <- asLogical(tclvalue(showEditButtonVariable))
    suppress.icon.images <- asLogical(tclvalue(SuppressIconImagesVariable))
    retain.selections <- asLogical(tclvalue(retainSelectionsVariable))
    use.rgl <- asLogical(tclvalue(useRglVariable))
    options <- current
    options$ask.to.exit <- ask.to.exit
    options$ask.on.exit <- ask.on.exit
    options$quit.R.on.close <- quit.R.on.close
    options$number.messages <- number.messages
    options$retain.messages <- retain.messages
    options$use.markdown <- use.markdown
    options$rmd.template <- rmd.template
    options$use.knitr <- use.knitr
    options$rnw.template <- rnw.template
    options$log.font.family <- log.font.family
    options$default.font.family <- default.font.family
    options$log.font.size <- log.font.size
    options$default.font.size <- default.font.size
    options$scale.factor <- scale.factor
    options$log.width <- log.width
    options$log.height <- log.height
    options$log.commands <- log.commands
    options$output.height <- output.height
    options$console.output <- console.output
    options$default.contrasts <- contrasts
    options$grab.focus <- grab.focus
    options$double.click <- double.click
    options$sort.names <- sort.names
    options$show.edit.button <- show.edit.button
    options$suppress.icon.images <- suppress.icon.images
    options$retain.selections <- retain.selections
    options$use.rgl <- use.rgl
    colors <- c(hex.1, hex.2, hex.3, hex.4, hex.5, hex.6)
    colors <- rgb2col(colors)
    options$log.text.color <- colors[1]
    options$command.text.color <- colors[2]
    options$output.text.color <- colors[3]
    options$error.text.color <- colors[4]
    options$warning.text.color <- colors[5]
    options$title.color <- colors[6]
    options$theme <- theme
    options(Rcmdr=options)
    closeCommander()
    Commander()
  }
  OKCancelHelp(helpSubject="Commander")
  tkgrid(closeOptionsFrame, sticky="nw")
  tkgrid(labelRcmdr(fontFrame, text=gettextRcmdr("Dialog text font size (points)")), defaultFontSizeSlider, sticky="sw", padx=6)
  tkgrid(labelRcmdr(fontFrame, text=gettextRcmdr("Script and output font size (points)")), logFontSizeSlider, sticky="sw", padx=6)
  tkgrid(labelRcmdr(fontFrame, text=gettextRcmdr("Dialog font")), defaultFontEntry, sticky="w", padx=6)
  tkgrid(labelRcmdr(fontFrame, text=gettextRcmdr("Script and output font")), logFontEntry, sticky="w", padx=6)
  tkgrid(fontFrame, sticky="w")
  tkgrid(labelRcmdr(fontTab, text="")) 
  pal <- c(log.text.color, command.text.color, output.text.color, error.text.color, warning.text.color, title.color)
  tkgrid(labelRcmdr(fontColorsFrame, text=gettextRcmdr("Script text color ")), button1, colorField1, sticky="w", padx=6)
  tkgrid(labelRcmdr(fontColorsFrame, text=gettextRcmdr("Command text color ")), button2, colorField2, sticky="w", padx=6)
  tkgrid(labelRcmdr(fontColorsFrame, text=gettextRcmdr("Output text color ")), button3, colorField3, sticky="w", padx=6)
  tkgrid(labelRcmdr(fontColorsFrame, text=gettextRcmdr("Error text color ")), button4, colorField4, sticky="w", padx=6)
  tkgrid(labelRcmdr(fontColorsFrame, text=gettextRcmdr("Warning text color ")), button5, colorField5, sticky="w", padx=6)
  tkgrid(labelRcmdr(fontColorsFrame, text=gettextRcmdr("Dialog subtitles text color ")), button6, colorField6, sticky="w", padx=6)
  tkgrid(fontColorsFrame, sticky="w")
  tkgrid(labelRcmdr(outputSliderFrame, text=gettextRcmdr("Script window width (characters)")), logWidthSlider, sticky="sw", padx=6)
  tkgrid(labelRcmdr(outputSliderFrame, text=gettextRcmdr("Script window height (lines)")), logHeightSlider, sticky="sw", padx=6)
  tkgrid(labelRcmdr(outputSliderFrame, text=gettextRcmdr("Output window height (lines)")), outputHeightSlider, sticky="sw", padx=6)
  tkgrid(labelRcmdr(outputSliderFrame, text=gettextRcmdr("Messages window height (lines)")), messagesHeightSlider, sticky="sw", padx=6)
  tkgrid(outputSliderFrame, sticky="w")
  tkgrid(labelRcmdr(outputTab, text=" "), sticky="w")    
  tkgrid(outputOptionsFrame, sticky="nw", columnspan = 3)
  tkgrid(labelRcmdr(templateFrame, text="R Markdown template file"), rmdTemplateEntry, templateButton, sticky="w", padx=6)
  tkgrid(templateFrame, columnspan=2, sticky="w")
  tkgrid(labelRcmdr(templateFrame, text="R knitr template file"), rnwTemplateEntry, rnwTemplateButton, sticky="w", padx=6)
  tkgrid(labelRcmdr(scaleFactorFrame, text=gettextRcmdr("Scale factor for Tk elements")), scaleFactorSlider, sticky="sw")
  tkgrid(scaleFactorFrame, sticky="w")
  tkgrid(labelRcmdr(otherTab, text=""))
  tkgrid(labelRcmdr(contrastsFrame, text=""), labelRcmdr(contrastsFrame, text=gettextRcmdr("Unordered factors")), labelRcmdr(contrastsFrame, text="   "),
         labelRcmdr(contrastsFrame, text=gettextRcmdr("Ordered factors")), sticky="w")
  tkgrid(labelRcmdr(contrastsFrame, text=gettextRcmdr("Contrasts")), contrasts1Entry, labelRcmdr(contrastsFrame, text="   "), contrasts2Entry, sticky="sw")
  tkgrid(contrastsFrame, sticky="sw")
  tkgrid(labelRcmdr(otherTab, text=" "), sticky="w")    
  tkgrid(otherOptionsFrame, sticky="w", columnspan=2)
  tkgrid(labelRcmdr(otherTab, text=""))
  tkgrid(getFrame(themesBox), sticky="w")
  tkgrid(labelRcmdr(otherTab, text=""))
  tkadd(notebook, closeTab, text=gettextRcmdr("Exit"), padding=6)
  tkadd(notebook, fontTab, text=gettextRcmdr("Fonts"), padding=6)
  tkadd(notebook, outputTab, text=gettextRcmdr("Output"), padding=6)
  tkadd(notebook, otherTab, text=gettextRcmdr("Other Options"), padding=6)
  tkgrid(notebook)
  tkconfigure(OKbutton, text=gettextRcmdr("Restart R Commander"), width=20)
  tkgrid(buttonsFrame, columnspan=3, sticky="ew")
  dialogSuffix()
}

saveOptions <- function(){
    initializeDialog(title=gettextRcmdr("Save Commander Options"))
    onCopy <- function(){
        focused <- tkfocus()
        if (tclvalue(focused) != optionsWindow$ID) focused <- optionsWindow
        selection <- strsplit(tclvalue(tktag.ranges(focused, "sel")), " ")[[1]]
        if (is.na(selection[1])) return()
        text <- tclvalue(tkget(focused, selection[1], selection[2]))
        tkclipboard.clear()
        tkclipboard.append(text)
    }
    onDelete <- function(){
        focused <- tkfocus()
        if (tclvalue(focused) != optionsWindow$ID) focused <- optionsWindow
        selection <- strsplit(tclvalue(tktag.ranges(focused, "sel")), " ")[[1]]
        if (is.na(selection[1])) return()
        tkdelete(focused, selection[1], selection[2])
    }
    onCut <- function(){
        onCopy()
        onDelete()
    }
    onPaste <- function(){
        onDelete()
        focused <- tkfocus()
        if (tclvalue(focused) != optionsWindow$ID) focused <- optionsWindow
        text <- tclvalue(.Tcl("selection get -selection CLIPBOARD"))
        if (length(text) == 0) return()
        tkinsert(focused, "insert", text)
    }
    onFind <- function(){
        focused <- tkfocus()
        if (tclvalue(focused) != optionsWindow$ID) focused <- optionsWindow
        initializeDialog(title=gettextRcmdr("Find"))
        textFrame <- tkframe(top)
        textVar <- tclVar(getRcmdr("last.search"))
        textEntry <- ttkentry(textFrame, width="20", textvariable=textVar)
        checkBoxes(frame="optionsFrame", boxes=c("regexpr", "case"), initialValues=c("0", "1"),
                   labels=gettextRcmdr(c("Regular-expression search", "Case sensitive")))
        radioButtons(name="direction", buttons=c("foward", "backward"), labels=gettextRcmdr(c("Forward", "Backward")),
                     values=c("-forward", "-backward"), title=gettextRcmdr("Search Direction"))
        onOK <- function(){
            text <- tclvalue(textVar)
            putRcmdr("last.search", text)
            if (text == ""){
                errorCondition(recall=onFind, message=gettextRcmdr("No search text specified."))
                return()
            }
            type <- if (tclvalue(regexprVariable) == 1) "-regexp" else "-exact"
            case <- tclvalue(caseVariable) == 1
            direction <- tclvalue(directionVariable)
            stop <- if (direction == "-forward") "end" else "1.0"
            where <- if (case) tksearch(focused, type, direction, "--", text, "insert", stop)
            else tksearch(focused, type, direction, "-nocase", "--", text, "insert", stop)
            where <- tclvalue(where)
            if (where == "") {
                Message(message=gettextRcmdr("Text not found."),
                        type="note")
                if (GrabFocus()) tkgrab.release(top)
                tkdestroy(top)
                tkfocus(CommanderWindow())
                return()
            }
            if (GrabFocus()) tkgrab.release(top)
            tkfocus(focused)
            tkmark.set(focused, "insert", where)
            tksee(focused, where)
            tkdestroy(top)
        }
        .exit <- function(){
            text <- tclvalue(textVar)
            putRcmdr("last.search", text)
            return("")
        }
        OKCancelHelp()
        tkgrid(labelRcmdr(textFrame, text=gettextRcmdr("Search for:")), textEntry, sticky="w")
        tkgrid(textFrame, sticky="w")
        tkgrid(optionsFrame, sticky="w")
        tkgrid(directionFrame, sticky="w")
        tkgrid(buttonsFrame, sticky="w")
        dialogSuffix(focus=textEntry)
    }
    onSelectAll <- function() {
        focused <- tkfocus()
        if (tclvalue(focused) != optionsWindow$ID) focused <- optionsWindow
        tktag.add(focused, "sel", "1.0", "end")
        tkfocus(focused)
    }
    onClear <- function(){
        onSelectAll()
        onDelete()
    }
    onUndo <- function(){
        focused <- tkfocus()
        if (tclvalue(focused) != optionsWindow$ID) focused <- optionsWindow
        tcl(focused, "edit", "undo")
    }
    onRedo <- function(){
        focused <- tkfocus()
        if (tclvalue(focused) != optionsWindow$ID) focused <- optionsWindow
        tcl(focused, "edit", "redo")
    }
    onOK <- function(){
        optionsFile <- tclvalue(tkgetSaveFile(filetypes=gettextRcmdr('{"All Files" {"*"}}'),
                                              initialfile=".Rprofile",
                                              parent=top))
        if (optionsFile == "") return()
        options <- tclvalue(tkget(optionsWindow, "1.0", "end"))
        fileCon <- file(optionsFile, "w")
        cat(options, file = fileCon)
        close(fileCon)
        closeDialog(top)
        Message(paste(gettextRcmdr("R Profile saved to"), optionsFile), type="note")
    }
    contextMenu <- function(){
        contextMenu <- tkmenu(tkmenu(optionsWindow), tearoff=FALSE)
        tkadd(contextMenu, "command", label=gettextRcmdr("Cut"), command=onCut)
        tkadd(contextMenu, "command", label=gettextRcmdr("Copy"), command=onCopy)
        tkadd(contextMenu, "command", label=gettextRcmdr("Paste"), command=onPaste)
        tkadd(contextMenu, "command", label=gettextRcmdr("Delete"), command=onDelete)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Find..."), command=onFind)
        tkadd(contextMenu, "command", label=gettextRcmdr("Select all"), command=onSelectAll)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Undo"), command=onUndo)
        tkadd(contextMenu, "command", label=gettextRcmdr("Redo"), command=onRedo)
        tkadd(contextMenu, "separator")
        tkadd(contextMenu, "command", label=gettextRcmdr("Clear window"), command=onClear)
        tkpopup(contextMenu, tkwinfo("pointerx", optionsWindow), tkwinfo("pointery", optionsWindow))
    }
    optionsFrame <- tkframe(top)    
    optionsWindow <- tktext(optionsFrame, bg="white", foreground=getRcmdr("log.text.color"),
                            font=getRcmdr("logFont"), height=20, width=65, wrap="none", undo=TRUE)
    optionsXscroll <- ttkscrollbar(optionsFrame, orient="horizontal",
                                   command=function(...) tkxview(optionsWindow, ...))
    optionsYscroll <- ttkscrollbar(optionsFrame,
                                   command=function(...) tkyview(optionsWindow, ...))
    tkconfigure(optionsWindow, xscrollcommand=function(...) tkset(optionsXscroll, ...))
    tkconfigure(optionsWindow, yscrollcommand=function(...) tkset(optionsYscroll, ...))
    tkbind(top, "<Control-x>", onCut)
    tkbind(top, "<Control-X>", onCut)
    tkbind(top, "<Control-c>", onCopy)
    tkbind(top, "<Control-C>", onCopy)
    tkbind(top, "<Control-f>", onFind)
    tkbind(top, "<Control-F>", onFind)
    tkbind(top, "<F3>", onFind)
    tkbind(top, "<Control-a>", onSelectAll)
    tkbind(top, "<Control-A>", onSelectAll)
    tkbind(top, "<Control-w>", onRedo)
    tkbind(top, "<Control-W>", onRedo)
    if (MacOSXP()){
        tkbind(top, "<Meta-x>", onCut)
        tkbind(top, "<Meta-X>", onCut)
        tkbind(top, "<Meta-c>", onCopy)
        tkbind(top, "<Meta-C>", onCopy)
        tkbind(top, "<Meta-v>", onPaste)
        tkbind(top, "<Meta-V>", onPaste)
        tkbind(top, "<Meta-f>", onFind)
        tkbind(top, "<Meta-F>", onFind)
        tkbind(top, "<Meta-a>", onSelectAll)
        tkbind(top, "<Meta-A>", onSelectAll)
        tkbind(top, "<Meta-w>", onRedo)
        tkbind(top, "<Meta-W>", onRedo)
    }
    tkbind(top, "<Alt-BackSpace>", onUndo)
    tkbind(optionsWindow, "<ButtonPress-3>", contextMenu)
    OKCancelHelp(helpSubject="saveOptions")
    menu <- tkmenu(top)
    tkconfigure(top, menu=menu)
    editMenu <- tkmenu(menu, tearoff=FALSE)
    tkadd(editMenu, "command", label=gettextRcmdr("Cut"), command=onCut)
    tkadd(editMenu, "command", label=gettextRcmdr("Copy"), command=onCopy)
    tkadd(editMenu, "command", label=gettextRcmdr("Paste"), command=onPaste)
    tkadd(editMenu, "command", label=gettextRcmdr("Delete"), command=onDelete)
    tkadd(editMenu, "separator")
    tkadd(editMenu, "command", label=gettextRcmdr("Find..."), command=onFind)
    tkadd(editMenu, "command", label=gettextRcmdr("Select all"), command=onSelectAll)
    tkadd(editMenu, "separator")
    tkadd(editMenu, "command", label=gettextRcmdr("Undo"), command=onUndo)
    tkadd(editMenu, "command", label=gettextRcmdr("Redo"), command=onRedo)
    tkadd(editMenu, "separator")
    tkadd(editMenu, "command", label=gettextRcmdr("Clear window"), command=onClear)  
    tkadd(menu, "cascade", label=gettextRcmdr("Edit"), menu=editMenu)
    Rprofile <- if (file.exists(".Rprofile")) readLines(".Rprofile") else ""
    start <- grep("^###! Rcmdr Options Begin !###", Rprofile)
    end <- grep("^###! Rcmdr Options End !###", Rprofile)
    if (length(start) == 1 && length(end) == 1){
        Rprofile <- Rprofile[-(start:end)]
    }
    Rprofile <- sub("\\n*$", "", Rprofile)
    tkinsert(optionsWindow, "end", paste(Rprofile, collapse="\n"))
    options <- getOption("Rcmdr")
    con <- file(open="w+")
    dput(options, con)
    options <- readLines(con)
    close(con)
    options <- paste(c("", "",
                       "###! Rcmdr Options Begin !###",
                       "options(Rcmdr=",
                       options,
                       ")",
                       "",
                       "# Uncomment the following 4 lines (remove the #s)",
                       "#  to start the R Commander automatically when R starts:",
                       "",
                       "# local({",
                       "#    old <- getOption('defaultPackages')",
                       "#    options(defaultPackages = c(old, 'Rcmdr'))",
                       "# })",
                       "",
                       "###! Rcmdr Options End !###"),
                     collapse="\n"
    )
    tkinsert(optionsWindow, "end", options)
    tkgrid(labelRcmdr(top, text=
                          paste(gettextRcmdr("The following commands will be saved in the file .Rprofile.",
                                             "You may edit this file before saving it."), collapse="\n")),
           sticky="w")
    tkgrid(optionsWindow, optionsYscroll, sticky="news")
    tkgrid(optionsXscroll, sticky="ew", columnspan=2)
    tkgrid(optionsFrame, sticky="news", padx=10, pady=0)
    tkgrid(buttonsFrame, sticky="ew")
    dialogSuffix()
}

loadPackages <- function(){
	availablePackages <- sort(setdiff(.packages(all.available = TRUE), .packages()))
	if (length(availablePackages) == 0){
		errorCondition(message=gettextRcmdr("No packages available to load."))
		return()
	}
	initializeDialog(title=gettextRcmdr("Load Packages"))
	packagesBox <- variableListBox(top, availablePackages, title=gettextRcmdr("Packages (pick one or more)"),
			selectmode="multiple", listHeight=10)
	onOK <- function(){
		packages <- getSelection(packagesBox)
		closeDialog(top)
		if (length(packages) == 0){
			errorCondition(recall=loadPackages, message=gettextRcmdr("You must select at least one package."))
			return()
		}
		for (package in packages) {
			Library(package)
		}
		Message(paste(gettextRcmdr("Packages loaded:"), paste(packages, collapse=", ")), type="note")
	}
	OKCancelHelp(helpSubject="library")
	tkgrid(getFrame(packagesBox), sticky="nw")
	tkgrid(buttonsFrame, sticky="w")
	dialogSuffix()
}

Setwd <- function(){
	wd <- tclvalue(tkchooseDirectory(initialdir=getwd(), parent=CommanderWindow()))
	if (wd != "") doItAndPrint(paste('setwd("', wd, '")', sep=""))
}


editMarkdown <- function(){
  .rmd <- RmdWindow()
  buffer <- tclvalue(tkget(.rmd, "1.0", "end"))
  compile <- function() {
    .rmd <- RmdWindow()
    editor <- getRcmdr("editor.text")
    buffer <- tclvalue(tkget(editor, "1.0", "end"))
    tkdelete(.rmd, "1.0", "end")
    tkinsert(.rmd, "end", buffer)
    compileRmd()
  }
  removeLastBlock <- function(){
    .rmd <- RmdWindow()
    editor <- getRcmdr("editor.text")
    buffer <- tclvalue(tkget(editor, "1.0", "end"))
    tkdelete(.rmd, "1.0", "end")
    tkinsert(.rmd, "end", buffer)
    removeLastRmdBlock()
    buffer <- tclvalue(tkget(.rmd, "1.0", "end"))
    tkdelete(editor, "1.0", "end")
    tkinsert(editor, "end", buffer)
  }
  ok <- function(){
    .rmd <- RmdWindow()
    editor <- getRcmdr("editor.text")
    buffer <- tclvalue(tkget(editor, "1.0", "end"))
    tkdelete(.rmd, "1.0", "end")
    tkinsert(.rmd, "end", buffer)
  }
  RcmdrEditor(buffer,  title="Edit R Markdown document", ok=ok,
              help=list(label="Using R Markdown", command=browseRMarkdown),
              file.menu=list(list(label="Generate report", command=compile)), 
              toolbar.buttons=list(list(label="Generate report", command=compile, 
                                        image="::image::submitIcon")))
#   edited <- getRcmdr("buffer")
#   if (!is.null(edited)){
#     tkdelete(.rmd, "1.0", "end")
#     tkinsert(.rmd, "end", edited)
#     tkyview.moveto(.rmd, 1)
#   }
}

editKnitr <- function(){
    .rnw <- RnwWindow()
    buffer <- tclvalue(tkget(.rnw, "1.0", "end"))
    compile <- function() {
        .rnw <- RnwWindow()
        editor <- getRcmdr("editor.text")
        buffer <- tclvalue(tkget(editor, "1.0", "end"))
        tkdelete(.rnw, "1.0", "end")
        tkinsert(.rnw, "end", buffer)
        compileRnw()
    }
    removeLastBlock <- function(){
        .rnw <- RnwWindow()
        editor <- getRcmdr("editor.text")
        buffer <- tclvalue(tkget(editor, "1.0", "end"))
        tkdelete(.rnw, "1.0", "end")
        tkinsert(.rnw, "end", buffer)
        removeLastRnwBlock()
        buffer <- tclvalue(tkget(.rnw, "1.0", "end"))
        tkdelete(editor, "1.0", "end")
        tkinsert(editor, "end", buffer)
    }
    ok <- function(){
      .rnw <- RnwWindow()
      editor <- getRcmdr("editor.text")
      buffer <- tclvalue(tkget(editor, "1.0", "end"))
      tkdelete(.rnw, "1.0", "end")
      tkinsert(.rnw, "end", buffer)
    }
    RcmdrEditor(buffer,  title="Edit knitr document", ok=ok,
        file.menu=list(list(label="Generate PDF report", command=compile)), 
        toolbar.buttons=list(list(label="Generate PDF report", command=compile, image="::image::submitIcon")))
#     edited <- getRcmdr("buffer")
#     if (!is.null(edited)){
#         tkdelete(.rnw, "1.0", "end")
#         tkinsert(.rnw, "end", edited)
#         tkyview.moveto(.rnw, 1)
#     }
}

appNap <- function(){
  initializeDialog(title=gettextRcmdr("Mac OS X app nap for R.app"))
  radioButtons(name="appnap", buttons=c("off", "on"), labels=gettextRcmdr(c("off (recommended)", "on")),
               title=gettextRcmdr("Set app nap"), initialValue=appnap())
  onOK <- function(){
    setting <- tclvalue(appnapVariable)
    appnap(setting)
    closeDialog()
  }
  OKCancelHelp(helpSubject="Commander")
  tkgrid(appnapFrame, sticky="w")
  tkgrid(buttonsFrame, sticky="ew")
  dialogSuffix()
}

installSoftware <- function(){
    initializeDialog(title=gettextRcmdr("Install Auxiliary Software"))
    has <- unlist(getRcmdr("capabilities"))
    installed <- c("", gettextRcmdr("(already installed)"))[1 + has]
    checkBoxes(frame="selectSoftwareFrame", boxes=c("latex", "pandoc"),
        initialValues=!as.numeric(has),
        labels=paste(gettextRcmdr(c("LaTeX", "Pandoc")), installed),
        title=gettextRcmdr("Software to Install"))
    onOK <- function(){
        if (tclvalue(latexVariable) == "1"){
            if (WindowsP()) browseURL("http://miktex.org/download")
            else if (MacOSXP()) browseURL("http://www.tug.org/mactex/")
            else browseURL("http://latex-project.org/ftp.html")
        }
        if (tclvalue(pandocVariable) == "1") browseURL("http://johnmacfarlane.net/pandoc/installing.html")
        closeDialog()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="AuxiliarySoftware")
    tkgrid(labelRcmdr(top, text = paste(
        gettextRcmdr("Please read the help for this dialog\nbefore installing auxiliary software."), "\n")))
    tkgrid(selectSoftwareFrame, sticky="w")
    dialogSuffix(grid.buttons=TRUE)
}
