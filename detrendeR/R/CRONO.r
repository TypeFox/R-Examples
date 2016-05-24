CRONO =function (input = "", ...) 
{
    flag = FALSE
    if (substring(input, 1, 1) == "<") 
        flag <- TRUE
    if (flag) 
        input <- ""
    filenamevar <- tclVar(input)
    output = ""
    if (input != "") 
        output <- paste(input, "crn", sep = ".")
    tabnamevar <- tclVar(output)
    fictrt <- tclVar()
    arValue <- tclVar(arMAX)
    listDataSets = function(envir = .GlobalEnv, ...) {
        Vars <- ls(envir = envir, all.names = TRUE)
        if (length(Vars) == 0) 
            return(Vars)
        out = names(which(sapply(Vars, function(.x) is.data.frame(get(.x, 
            envir = envir)) || is.matrix(get(.x, envir = envir)))))
        out
    }
    tfCRONO <- tktoplevel()
    tkwm.title(tfCRONO, "Make chrono")
    tkwm.resizable(tfCRONO, 0, 0)
    tkwm.geometry(tfCRONO, paste("+0+", .heigth, sep = ""))
    tkwm.deiconify(tfCRONO)
    tkgrab.set(tfCRONO)
    tkfocus(tfCRONO)
    choose.data = function() {
        input <- tk_select.list(sort(listDataSets()), title = "Select one")
        output <- paste(input, "crn", sep = ".")
        tkgrab.set(tfCRONO)
        if (input != "") {
            tkdelete(file.entry, 0, "end")
            tkinsert(file.entry, "end", input)
            tkdelete(tab.entry, 0, "end")
            tkinsert(tab.entry, "end", output)
        }
    }
    done <- tclVar(0)
    frame1.parent <- tkframe(tfCRONO, relief = "groove")
    frame1.a <- tkframe(frame1.parent, relief = "groove")
    frame1 <- tkframe(frame1.parent, relief = "groove")
    tkgrid(tklabel(frame1.a, text = "Options:", foreground = "blue"))
    tkpack(frame1.a, fill = "x")
    tab.entry <- tkentry(frame1, textvariable = tabnamevar)
    file.entry <- tkentry(frame1, textvariable = filenamevar)
    choosefile.but <- tkbutton(frame1, text = "...", command = function() choose.data())
    tkgrid(tklabel(frame1, text = "Input name: "), file.entry, 
        tklabel(frame1, text = " "), choosefile.but, sticky = "w")
    tkgrid(tklabel(frame1, text = "Output name:"), tab.entry, 
        sticky = "w")
    tkpack(frame1, fill = "x")
    frame3 <- tkframe(frame1.parent, relief = "groove")
    tkgrid(tklabel(frame3, text = "\nPrewhitened chronology:", 
        foreground = "blue"), columnspan = 1)
    tkpack(frame3, fill = "x")
    frame3.1 <- tkframe(frame1.parent)
    makeAr.value <- tclVar(makeAr)
    arMAX.value <- tclVar(arMAX)
    makeAr.cbut <- tkcheckbutton(frame3.1, text = "AR model of max order:", 
        variable = makeAr.value)
    slider <- tkscale(frame3.1, from = 1, to = 10, showvalue = T, 
        variable = arValue, resolution = 1, orient = "horizontal")
    tkgrid(makeAr.cbut, slider)
    tkpack(frame3.1, fill = "x")
    frame4.0 <- tkframe(frame1.parent, relief = "groove")
    tkgrid(tklabel(frame4.0, text = "Mean:", foreground = "blue"))
    tkpack(frame4.0, fill = "x")
    frame4 <- tkframe(frame1.parent, relief = "groove")
    rb1 <- tkradiobutton(frame4)
    rb2 <- tkradiobutton(frame4)
    rbValue <- tclVar(biweightMean)
    tkconfigure(rb1, variable = rbValue, value = TRUE)
    tkconfigure(rb2, variable = rbValue, value = FALSE)
    tkgrid(tklabel(frame4, text = "Robust     "), rb1)
    tkgrid(tklabel(frame4, text = "Arithmetic "), rb2)
    tkpack(frame4, fill = "x")
    frame.exit <- tkframe(frame1.parent, relief = "groove")
    OnOk = function() {
        flag <- try(exists(tclvalue(filenamevar)), silent = T)
        if (flag == TRUE) {
            tclvalue(done) <- 2
            eval(parse(text = paste("temp<-", tclvalue(filenamevar))))
            if (as.logic(tclvalue(makeAr.value))) 
                eval(parse(text = paste("temp <- apply(temp,  2,ar.func, order.max=", 
                  as.numeric(tclvalue(arValue)), ")", sep = "")))
            .assign("arMAX", as.numeric(tclvalue(arValue)))
            .assign("makeAr", as.logic(tclvalue(makeAr.value)))
            .assign("biweightMean", as.logic(tclvalue(rbValue)))
            eval(parse(text = paste(tclvalue(tabnamevar), "<<-Chron(temp, stc=stc, biweight=", 
                as.logic(tclvalue(rbValue)), ")", sep = "")))
        }
    }
    fr.exit.space <- tklabel(frame.exit, text = " ")
    ok.but <- tkbutton(frame.exit, text = "      Ok      ", command = function() OnOk())
    cancel.but <- tkbutton(frame.exit, text = "    Cancel    ", 
        command = function() tkdestroy(tfCRONO))
    tkgrid(cancel.but, fr.exit.space, ok.but)
    tkpack(frame1.parent, side = "top")
    tkpack(frame.exit, side = "right")
    tkbind(tfCRONO, "<Destroy>", function() tclvalue(done) <- 2)
    tkbind(tfCRONO, "<KeyPress-Return>", function() OnOk())
    tkbind(tfCRONO, "<KeyPress-Escape>", function() tkdestroy(tfCRONO))
    tkwait.variable(done)
    tkgrab.release(tfCRONO)
    if (tclvalue(done) == "2") 
        tkdestroy(tfCRONO)
}
