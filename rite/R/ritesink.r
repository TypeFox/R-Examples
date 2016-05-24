sinkstart <- function(
    echo=TRUE, split=FALSE,
    fontFamily='Courier', fontSize=10,
    col.bg='white', col.call=c('black',col.bg), col.result=c('black',col.bg),
    col.err=c('red',col.bg), col.warn=c('purple',col.bg), col.msg=c('blue',col.bg))
{
    if('outsink' %in% getTaskCallbackNames())
        stop('ritesink is already active')
    
    riteenv$echo <- echo
    riteenv$split <- split
    
    # setup colors
    if(is.null(col.bg) || is.na(col.bg) || col.bg=='')
        col.bg <- 'white'
    if(length(col.bg)>1){
        col.bg <- col.bg[1]
        warning('More than one color specified for background. Only first is used.')
    }
    if(is.null(col.call))
        col.call <- c('black',col.bg)
    if(is.null(col.result))
        col.result <- col.call
    if(is.null(col.err))
        col.err <- c('red',col.bg)
    if(is.null(col.warn))
        col.warn <- c('purple',col.bg)
    if(is.null(col.msg))
        col.msg <- c('blue',col.bg)
    if(length(col.call)==1)
        col.call <- c(col.call,col.bg)
    if(length(col.result)==1)
        col.result <- c(col.result,col.bg)
    if(length(col.err)==1)
        col.err <- c(col.err,col.bg)
    if(length(col.warn)==1)
        col.warn <- c(col.warn,col.bg)
    if(length(col.msg)==1)
        col.msg <- c(col.msg,col.bg)
    if(col.bg==col.call[1])
        stop('Background and call foreground colors are the same')
    if(col.bg==col.result[1])
        stop('Background and result foreground colors are the same')
    
    # sinks
    riteenv$stdsink <- file(riteenv$otmp <- tempfile(),'w+')
    sink(riteenv$stdsink, split=split, type='output') # output
    writeLines('# rite output sink', riteenv$stdsink)
    riteenv$lengtho <- paste(scan(riteenv$stdsink, what='character',
                                sep='\n', quiet=TRUE),collapse='\n')
    
    riteenv$errsink <- file(riteenv$etmp <- tempfile(),'w+')
    sink(riteenv$errsink, type='message') # message
    writeLines('# rite error sink', riteenv$errsink)
    riteenv$lengthe <- paste(scan(riteenv$errsink, what='character',
                                sep='\n', quiet=TRUE),collapse='\n')
    
    # error handler
    riteenv$errhandler <- function(){
        tkinsert(riteenv$output,'insert',paste('\n', geterrmessage(), sep=''), ('error'))
        tksee(riteenv$output, 'insert')
        invisible(NULL)
    }
        
    # callback function
    outsink <- function() {
        function(expr, value, ok, visible) {
            #e <- as.character(as.expression(expr))
            e <- deparse(expr)
            if(ok){
                if(riteenv$echo)
                    tkinsert(riteenv$output, 'insert', paste('\n>',e), ('call'))
                # Output sink (for `cat` and `print`)
                osink <- paste(scan(riteenv$otmp, what='character',
                                sep='\n', quiet=TRUE),collapse='\n')
                if(!identical(osink, riteenv$lengtho) && length(osink)>0){
                    last <- substr(osink, nchar(riteenv$lengtho)+1, nchar(osink))
                    # handle `simpleError` etc. that trigger callback
                    if(grepl("simpleError", last))
                        tkinsert(riteenv$output, 'insert', last, ('error'))
                    else if(grepl("simpleWarning", last))
                        tkinsert(riteenv$output, 'insert', last, ('warning'))
                    else if(grepl("simpleMessage", last))
                        tkinsert(riteenv$output, 'insert', last, ('message'))
                    else if(grepl("simpleCondition", last))
                        tkinsert(riteenv$output, 'insert', last, ('message'))
                    else
                        tkinsert(riteenv$output, 'insert', last, ('result'))
                    riteenv$lengtho <- osink
                }
                #tkinsert(riteenv$output,'insert','ok\n', ('message')) # print confirm 'ok' on non-printing calls
                riteenv$lengtho <- paste(scan(riteenv$otmp, what='character',
                                        sep='\n', quiet=TRUE),collapse='\n')
            }
            else if(visible && !ok) # !ok doesn't happen
                tkinsert(riteenv$output,'insert','Error\n', ('error'))
            else if(!visible && !ok) # !ok doesn't happen
                tkinsert(riteenv$output,'insert','Non-printing error\n', ('error'))
            
            # Error sink (for `warning` and `message`)
            esink <- paste(scan(riteenv$etmp, what='character',
                            sep='\n', quiet=TRUE),collapse='\n')
            if(!identical(esink, riteenv$lengthe) && length(esink)>0){
                fromerr <- substr(esink,nchar(riteenv$lengthe)+1,nchar(esink))
                if(any(grepl('error', fromerr))){
                    tkinsert(riteenv$output, 'insert',
                        paste(fromerr, '\n', collapse='\n'), ('error'))
                }
                else if(any(grepl('Warning', fromerr))){
                    tkinsert(riteenv$output, 'insert',
                        paste(fromerr, collapse='\n'), ('warning'))
                }
                else{
                    tkinsert(riteenv$output, 'insert',
                        paste(fromerr, collapse='\n'), ('message'))
                }
                riteenv$lengthe <- esink
            }
            tksee(riteenv$output, 'insert')
            TRUE
        }
    }
    addTaskCallback(outsink(), name='outsink')
        
    if(!'thesink' %in% ls(riteenv)){
        riteenv$thesink <- tktoplevel(borderwidth=0)
        exitsink <- function() {
            tkdestroy(riteenv$thesink)
            if(!sink.number()==0)
                sink()
            if(!sink.number('message')==2)
                sink(type='message')
            rm(thesink, output, scr, stdsink, errsink, envir = riteenv)
            sinkstop()
            invisible()
        }

        # widget
        tkwm.protocol(riteenv$thesink, 'WM_DELETE_WINDOW', exitsink)
        tkwm.title(riteenv$thesink, 'rite sink')        # title
        riteenv$scr <- tkscrollbar(riteenv$thesink, 
                        repeatinterval=25,
                        command=function(...){ tkyview(riteenv$output,...) })
        riteenv$output <- tktext(riteenv$thesink, bg=col.bg, fg=col.result[1], undo='true',
                                   yscrollcommand=function(...) tkset(riteenv$scr,...),
                                   font=tkfont.create(family=fontFamily, size=fontSize))
        tcl('wm', 'attributes', riteenv$thesink, topmost=TRUE)
        tkgrid(riteenv$output, sticky='nsew', column=1, row=1)
        tkgrid(riteenv$scr, sticky='nsew', column=2, row=1)
        tkgrid.columnconfigure(riteenv$thesink,1,weight=1)
        tkgrid.columnconfigure(riteenv$thesink,2,weight=0)
        tkgrid.rowconfigure(riteenv$thesink,1,weight=1)
        
        # tags/fonts
        if(!exists('riteenv$defaultfont'))
            riteenv$defaultfont <- tkfont.create(family=fontFamily, size=fontSize)
        tktag.configure(riteenv$output, 'call',
            foreground=col.call[1],
            background=col.call[2],
            font=riteenv$defaultfont,
            underline=0)
        tktag.configure(riteenv$output, 'result',
            foreground=col.result[1],
            background=col.result[2],
            font=riteenv$defaultfont,
            underline=0)
        if(!exists('riteenv$boldfont'))
            riteenv$boldfont <- tkfont.create(family=fontFamily, size=fontSize, weight='bold')
        tktag.configure(riteenv$output, 'error',
            foreground=col.err[1],
            background=col.err[2],
            font=riteenv$boldfont,
            underline=0)
        tktag.configure(riteenv$output, 'warning',
            foreground=col.warn[1],
            background=col.warn[2],
            font=riteenv$boldfont,
            underline=0)
        tktag.configure(riteenv$output, 'message',
            foreground=col.msg[1],
            background=col.msg[2],
            font=riteenv$boldfont,
            underline=0)
        
        # bind option('width') to window resize
        resize <- function(){
            w <- tkwinfo('width',riteenv$output)
            m <- tkfont.measure(riteenv$defaultfont,'m')
            nw <- round((as.numeric(w)-20)/as.numeric(m))
            options(width=nw)
        }
        tkbind(riteenv$thesink, '<Configure>', resize)
        
        # context menu (and associated functions and bindings)
        saveSink <- function(){
            outfilename <- tclvalue(tkgetSaveFile(initialdir=getwd(),
                            title='Save Output',
                            filetypes=paste('{{Text file} {*.txt}} {{All files} {*.*}}'),
                            defaultextension='.txt'))
            if(!length(outfilename) || outfilename=="")
                invisible()
            chn <- tclopen(outfilename, 'w')
            tclputs(chn, tclvalue(tkget(riteenv$output,'0.0','end')))
            tclclose(chn)
        }
        tkbind(riteenv$output, '<Control-S>', expression(saveSink, break))
        tkbind(riteenv$output, '<Control-s>', expression(saveSink, break))
        selectAll <- function(){
            tktag.add(riteenv$output,'sel','0.0','end')
            tkmark.set(riteenv$output,'insert','end')
        }
        tkbind(riteenv$output, "<Control-A>", expression(selectAll, break))
        tkbind(riteenv$output, "<Control-a>", expression(selectAll, break))
        copyText <- function(docut=FALSE){
            selrange <- strsplit(tclvalue(tktag.ranges(riteenv$output,"sel"))," ")[[1]]
            if(!tclvalue(tktag.ranges(riteenv$output,"sel"))==""){
                tkclipboard.clear()
                tkclipboard.append(tclvalue(tkget(riteenv$output, selrange[1], selrange[2])))
                if(docut==TRUE)
                    tkdelete(riteenv$output, selrange[1], selrange[2])
            }
            else {
                selectAll()
                copyText()
            }
        }
        pasteText <- function(){
            if("windows"==.Platform$OS.type)
                cbcontents <- readLines("clipboard")
            else if("unix"==Sys.getenv("OS"))
                cbcontents <- readLines(pipe("pbpaste"))
            else
                cbcontents <- ""
            tkinsert(riteenv$output, "insert", paste(cbcontents,collapse="\n"))
        }
        clearSink <- function() tkdelete(riteenv$output, '1.0', 'end')
        tkbind(riteenv$output, "<Control-L>", expression(clearSink, break))
        tkbind(riteenv$output, "<Control-l>", expression(clearSink, break))
        
        contextMenu <- tkmenu(riteenv$output, tearoff = FALSE)
            tkadd(contextMenu, "command", label = "Save <Ctrl-S>",
                command = saveSink)
            tkadd(contextMenu, "command", label = "Clear All <Ctrl-L>",
                command = clearSink)
            tkadd(contextMenu, "separator")
            tkadd(contextMenu, "command", label = "Select All <Ctrl-A>",
                command = selectAll)
            tkadd(contextMenu, "command", label = "Copy <Ctrl-C>",
                command = copyText)
            tkadd(contextMenu, "command", label = "Cut <Ctrl-X>",
                command = function() copyText(docut=TRUE))
            tkadd(contextMenu, "command", label = "Paste <Ctrl-V>",
                command = pasteText)
            tkadd(contextMenu, "separator")
            tkadd(contextMenu, "command",
                label = paste("Toggle echo on/off"),
                command = function(){
                    if(riteenv$echo)
                        riteenv$echo <- FALSE
                    else
                        riteenv$echo <- TRUE
                })
        rightClick <- function(x, y) {
            rootx <- as.integer(tkwinfo("rootx", riteenv$output))
            rooty <- as.integer(tkwinfo("rooty", riteenv$output))
            xTxt <- as.integer(x) + rootx
            yTxt <- as.integer(y) + rooty
            tkmark.set(riteenv$output,"insert",paste("@",xTxt,",",yTxt,sep=""))
            .Tcl(paste("tk_popup", .Tcl.args(contextMenu, xTxt, yTxt)))
        }
        tkbind(riteenv$output, "<Button-3>", rightClick)
        tkbind(contextMenu, "<Button-3>", function() tkunpost(contextMenu))
    
    }
    
    options('show.error.messages'=FALSE) # default TRUE
    options('error'=riteenv$errhandler) # default NULL
    
    invisible(NULL)
}

sinkstop <- function(quiet = TRUE){
    # check for active sink
    if(!exists('riteenv'))
        stop('sink already closed and removed')
    # reset options defaults
    options('show.error.messages'=TRUE) # default TRUE
    options('error'=NULL) # default NULL
    options('width'=80) # default 80
    
    # stop sinks
    if(!sink.number()==0)
        sink()
    if(!sink.number('message')==2)
        sink(type='message')
    
    # remove call back
    if('outsink' %in% getTaskCallbackNames())
        removeTaskCallback('outsink')
    
    # close connections
    if(riteenv$otmp %in% showConnections()){
        thiscon <- rownames(showConnections())[which(riteenv$otmp==showConnections()[,1])]
        close(getConnection(thiscon))
    }
    if(riteenv$etmp %in% showConnections()){
        thiscon <- rownames(showConnections())[which(riteenv$etmp==showConnections()[,1])]
        close(getConnection(thiscon))
    }
    
    # remove temporary sink files
    unlink(riteenv$otmp)
    unlink(riteenv$etmp)
    
    # remove `riteenv`
    if(!"thesink" %in% ls(riteenv)){
        if(!quiet)
            message('rite sink closed and removed')
    }
    else if(!quiet)
        message('rite sink closed')
    
    invisible(NULL)
}

riteenv <- new.env()

if(getRversion() >= "2.15.1")
    utils::globalVariables(c('thesink', 'output', 'scr', 'stdsink', 'errsink'))
