# Long Running Task - windows only
do.long.running.task <- function(func, values, func.name = deparse(substitute(func)), input.name=NULL, output.name = "output"){
  tf1 <- "process_pid.txt"
  unlink(tf1)
  infile.text <- paste(
    # write PID to system pid file
  'cat(Sys.getpid(), file = "process_pid.txt")',
  'load("input.Rdata")',
  'retval <- NULL',
  'tryCatch(retval <- do.call(func, values), error = function(e) .GlobalEnv$retval <- e)',
  'save(retval, file="outfile.Rdata")', sep="\n")

    # remove the output file
  outfilename = "outfile.Rdata"
  unlink(outfilename)
  unlink("input.Rdata")
  the.formals <- formals(func)
  the.names <- as.character(unlist(the.formals[sapply(the.formals, class)=="name"]))
  the.names <- c(the.names, as.character(unlist(values[sapply(values, class)=="name"])))
  the.names <- the.names[nchar(the.names) > 0]
  cat("Starting new R session...\n")
  save(list=c("func", "values", ls(envir=.GlobalEnv)), file="input.Rdata")
  cat(infile.text, file="infile.R")
    #
  check_for_file <- function(data){
    outfilename <- data$outfilename
    pb <- data$pb
    tf1 <- data$tf1
    pb$pulse() # update progress bar
    dialog <- data$dialog
    
      # looking for file containing PID from process we started
    if(file.exists(tf1)){
      process.pid <- scan(tf1, quiet=T)
      unlink(tf1)
      cat(paste("Process pid is", process.pid, "\n"))
      dialog$setResponseSensitive(1, TRUE)
      gSignalConnect(dialog, "response", function(dlg, response) {
        print("Killing..")
        system(paste("taskkill /F /PID", process.pid), wait=F)
        gSourceRemove(timeout.id)
        dlg$destroy()
      })
    }

    rv <- TRUE
    if(file.exists(outfilename)){
      tryCatch({
        load(outfilename)  # Error reported in reading connection here
        dialog$destroy()
        unlink(outfilename)      
        retval <- get("retval")
        if("error"%in%class((retval))){
          quick_message(paste("An error occurred in", func.name, "\n\n", as.character(retval)))
        } else if(!is.null(retval) && !is.null(output.name)){
   	      assign.string = paste(".GlobalEnv$", output.name, " <- ", "retval", sep="")
    	    eval(parse(text=assign.string))
        }
        rv <- FALSE
      }, error = function(e) {
         cat(paste("Trying to open output:", e$message, "... Continuing...\n"))
      })
    }
    return(rv)    
  }

  system("R CMD BATCH infile.R outfile.R", wait=F)

  dialog <- gtkDialog("Running", NULL, c("modal","destroy-with-parent"), "gtk-cancel", 1,
                      show = FALSE)
  pb = gtkProgressBar()
  dialog[["vbox"]]$packStart(pb, TRUE, FALSE, 10)
  dialog$setPosition(GtkWindowPosition["center"])
  dialog$setResizable(FALSE)
  dialog$setResponseSensitive(1, FALSE)
  
  timeout.id <- gTimeoutAdd(1000, check_for_file,
    data=list(outfilename=outfilename, pb=pb, dialog=dialog, tf1=tf1))

  dialog$showAll()

}
