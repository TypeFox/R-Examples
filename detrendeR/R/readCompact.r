readCompact = function ()
{
    done <- tclVar(0)
	path<-tclvalue(tkgetOpenFile())
    filenamevar <- tclVar(path)
    tabnamevar <- tclVar(basename(path))
    header.var <- tclVar(0)
    show.data <- tclVar(0)
    header.bl = FALSE
    choosefile <- function() {
        fictrt <- tkgetOpenFile()
        fpath <- tclvalue(fictrt)
        tkfocus(tff)
        if (fpath != "") {
            tkdelete(file.entry, 0, "end")
            tkinsert(file.entry, "end", fpath)
            tkdelete(tab.entry, 0, "end")
            tkinsert(tab.entry, "end", basename(fpath))
        }
    }
    choosefic <- function() {
        if (tclvalue(tabnamevar) != "") {
            tabname <- parse(text = tclvalue(tabnamevar))[[1]]
        }
        else tabname <- "untitled"
        if (tclvalue(filenamevar) != "") {
            filename <- tclvalue(filenamevar)
        }
        else return()
        tkdestroy(tff)
           
        rdcom <- paste(tabname, " <<- read.compact('", filename, 
            "')", sep = "")
        eval(parse(text = rdcom))
        show.data.flag <- as.logical(tclObj(show.data))
        if (show.data.flag) 
            eval(parse(text = paste("edit(", tabname, ")", sep = "")))
    }
    tff <- tktoplevel()
    tkwm.title(tff, "Open a compact file")
    tkwm.resizable(tff, 0, 0)
    tkwm.geometry(tff, paste("+0+", .heigth, sep = ""))
    tkwm.deiconify(tff)
    tkgrab.set(tff)
    tkfocus(tff)
	
    frame1 <- tkframe(tff, relief = "groove")
    frame2 <- tkframe(tff, relief = "groove")
    frame.preview <- tkframe(tff, relief = "groove")
    tab.entry <- tkentry(frame1, textvariable = tabnamevar)
    file.entry <- tkentry(frame1, textvariable = filenamevar)
    separator <- tklabel(frame1, text = "")
    choosefile.but <- tkbutton(frame1, text = "...", command = function() choosefile())
    tkgrid(tklabel(frame1, text = "Select a file to read: "), 
        file.entry, separator, choosefile.but, sticky = "w")
    tkgrid(tklabel(frame1, text = "Enter name for data set: "), 
        tab.entry, sticky = "w")
    tkpack(frame1, fill = "x")
    
    show.data.cbut <- tkcheckbutton(frame2, text = "Show data frame                       ", 
        variable = show.data)

    tkpack(show.data.cbut, side = "left")

    tkpack(frame2)
    frame.exit <- tkframe(tff, relief = "groove")
    fr.exit.space <- tklabel(frame.exit, text = " ")
    ok.but <- tkbutton(frame.exit, text = "      Ok      ", command = function() choosefic())
    cancel.but <- tkbutton(frame.exit, text = "    Cancel    ", 
        command = function() tkdestroy(tff))
    tkgrid(cancel.but, fr.exit.space, ok.but)
    tkpack(frame.exit)
    tkfocus(tff)
	tkbind(tff, "<Destroy>", function() tclvalue(done) <- 2)
    tkbind(tff, "<KeyPress-Return>", function() choosefic())
    tkbind(tff, "<KeyPress-Escape>", function() tkdestroy(tff))
    tkwait.variable(done)
    if (tclvalue(done) == "2") 
    tkdestroy(tff)
}

#readCompact()

