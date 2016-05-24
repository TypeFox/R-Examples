# These functions for Excel supportwritten by Erich Neuwirth
#  last modified: 20 March 2008 by J. Fox  (following instructions from Erich Neuwirth)

    RExcelSupported <- function(){
    	RExcelSupport <- getOption("Rcmdr")$RExcelSupport
        !is.null(RExcelSupport) && RExcelSupport && exists("RExcelEnv") &&
            exists("putRExcel", where="RExcelEnv")
    	}




    SubmitToCommander <- function(commands){
        .log <- LogWindow()
        lines<-commands
        lines <- strsplit(lines, "\n")[[1]]
        .console.output <- getRcmdr("console.output")
        .output <- OutputWindow()
        iline <- 1
        nlines <- length(lines)
        while (iline <= nlines){
            while (whitespaceonly(lines[iline])) iline <- iline + 1
            if (iline > nlines) break
            current.line <- lines[iline]
            if (.console.output) cat(paste("\nRcmdr> ", current.line,"\n", sep=""))
            else{
                tkinsert(.output, "end", paste("\n> ", current.line,"\n", sep=""))
                tktag.add(.output, "currentLine", "end - 2 lines linestart", "end - 2 lines lineend")
                tktag.configure(.output, "currentLine", foreground=getRcmdr("command.text.color"))
                }
            jline <- iline + 1
            while (jline <= nlines){
                if (class(try(parse(text=current.line),silent=TRUE))!="try-error") break
                if (.console.output)cat(paste("Rcmdr+ ", lines[jline], sep="\n"))
                else{
                    tkinsert(.output, "end", paste("+ ", lines[jline],"\n", sep=""))
                    tktag.add(.output, "currentLine", "end - 2 lines linestart", "end - 2 lines lineend")
                    tktag.configure(.output, "currentLine", foreground=getRcmdr("command.text.color"))
                    }
                current.line <- paste(current.line, lines[jline],sep="\n")
                jline <- jline + 1
                iline <- iline + 1
                }
            if (!(is.null(current.line) || is.na(current.line))){
            if (length(grep("<-", current.line)) > 0){
                justDoIt(current.line)
            	loggerForExcel(current.line)
                }
            else if (length(grep("^remove\\(", current.line)) > 0){
                current.line <- sub(")", ", envir=.GlobalEnv)", current.line)
                justDoIt(current.line)
            	loggerForExcel(current.print.line)
                }
##            else if (any(sapply(Commander.Input.exceptions,
##                    function(.x) length(grep(paste("^", .x, "\\(", sep=""), current.line)) > 0))){
##                justDoIt(current.line)
##            	loggerForExcel(current.line)
##                }
            else if (length(current.line)>0) {
		          doItAndPrint(current.line, log=FALSE)
            	loggerForExcel(current.line)
		          }
            }
            iline <- iline + 1
        }
    tkyview.moveto(.output, 1)
    }




loggerForExcel <- function(command){
    .log <- LogWindow()
    .output <- OutputWindow()
    if (getRcmdr("log.commands")) {
        tkinsert(.log, "end", paste(command,"\n", sep=""))
        tkyview.moveto(.log, 1)
        }
    }

