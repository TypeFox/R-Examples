detrending = function (TwoSteps = TRUE, input = "", ...)
{

listDataSets = function (envir = .GlobalEnv, ...){
    Vars <- ls(envir = envir, all.names = TRUE)
    if (length(Vars) == 0)
        return(Vars)
    out=names(which(sapply(Vars, function(.x) is.data.frame(get(.x, envir = envir))||is.matrix(get(.x, envir = envir)))))
    out
    }
    
    flag = FALSE
    if (substring(input, 1, 1) == "<")
        flag <- TRUE
    if (flag)
        input <- ""
    filenamevar <- tclVar(input)
    output = ""
    if (input != "")
        output <- input
    tabnamevar <- tclVar(output)
    choose.data = function() {
        output <- input <- tk_select.list(sort(listDataSets()),
            preselect = NULL, multiple = FALSE, title = "Select one")
        tkgrab.set(top_detrending)
        if (input != "") {
            tkdelete(file.entry, 0, "end")
            tkinsert(file.entry, "end", input)
            tkdelete(tab.entry, 0, "end")
            tkinsert(tab.entry, "end", output)
        }
    }
    top_detrending <- tktoplevel()
    tkwm.geometry(top_detrending, paste("+0+", .heigth, sep = ""))
    tkwm.resizable(top_detrending, 0, 0)
    tkwm.title(top_detrending, "Detrending options")
    tkwm.deiconify(top_detrending)
    tkgrab.set(top_detrending)
    size = c(268, 180, 0, 132)
    done <- tclVar(0)
    frame1.b <- tkframe(top_detrending, relief = "groove", borderwidth = 2)
    tab.entry <- tkentry(frame1.b, textvariable = tabnamevar)
    file.entry <- tkentry(frame1.b, textvariable = filenamevar)
    choosefile.but <- tkbutton(frame1.b, text = "...", command = function() {
        choose.data()
    })
    tkgrid(tklabel(frame1.b, text = "Input name: ", foreground = "blue"),
        file.entry, tklabel(frame1.b, text = " "), choosefile.but,
        sticky = "w")
    tkgrid(tklabel(frame1.b, text = "Output name:", foreground = "blue"),
        tab.entry, sticky = "w")
    tkpack(frame1.b, fill = "x")
    top_detrending_frame2 <- tkframe(top_detrending, relief = "groove",
        borderwidth = 0)
    Det1frame <- tkwidget(top_detrending_frame2, "labelframe",
        foreground = "blue", text = "First detrend: ", relief = "groove",
        borderwidth = 2)
    Det2frame <- tkwidget(top_detrending_frame2, "labelframe",
        foreground = "blue", text = "Second detrend: ", relief = "groove",
        borderwidth = 2)
    Det1.1frame <- tkframe(Det1frame, relief = "groove", borderwidth = 0)
    Det1.2frame <- tkwidget(Det1frame, "labelframe", foreground = "blue",
        text = "Spline options: ", relief = "groove", borderwidth = 2)
    Det1.2.1frame <- tkframe(Det1.2frame, relief = "groove",
        borderwidth = 0)
    Det1.2.2frame <- tkframe(Det1.2frame, relief = "groove",
        borderwidth = 0)
    Det1.2.3frame <- tkframe(Det1.2frame, relief = "groove",
        borderwidth = 0)
    Det2.1frame <- tkframe(Det2frame, relief = "groove", borderwidth = 0)
    Det2.2frame <- tkframe(Det2frame, relief = "groove", borderwidth = 0)
    Det2.2frame <- tkwidget(Det2frame, "labelframe", foreground = "blue",
        text = "Spline options: ", relief = "groove", borderwidth = 2)
    Det2.2.1frame <- tkframe(Det2.2frame, relief = "groove",
        borderwidth = 0)
    Det2.2.2frame <- tkframe(Det2.2frame, relief = "groove",
        borderwidth = 0)
    Det2.2.3frame <- tkframe(Det2.2frame, relief = "groove",
        borderwidth = 0)
    RadioButton = function(FRAME, variable = NULL, BUTTON = c("b.r1",
        "b.r2"), VALUE = c(1, 2)) {
        BUTTON <- as.vector(BUTTON)
        for (i in 1:length(BUTTON)) {
            tkpack(tkradiobutton(FRAME, text = BUTTON[i], value = VALUE[i],
                variable = variable), anchor = "w")
        }
    }
    method1.value = tclVar(method1)
    method2.value = tclVar(method2)
    n1.value = tclVar(n1)
    nPerc1.value = tclVar(nPerc1)
    p1.value = tclVar(p1)
    n2.value = tclVar(n2)
    nPerc2.value = tclVar(nPerc2)
    p2.value = tclVar(p2)
    detrend.types = c("Neg Exp", "Spline", "Spline%", "Mean")
    detrend.values = c("ModNegExp", "Spline", "Spline%", "Mean")
    RadioButton(Det1.1frame, variable = method1.value, BUTTON = detrend.types,
        VALUE = detrend.values)
    RadioButton(Det2.1frame, variable = method2.value, BUTTON = detrend.types,
        VALUE = detrend.values)
    n1.entry <- tkentry(Det1.2.1frame, textvariable = n1.value,
        width = 5)
    Det1.2.1lab <- tklabel(Det1.2.1frame, text = "Spline length:")
    tkpack(Det1.2.1lab, n1.entry, side = "left")
    nPerc1.entry <- tkentry(Det1.2.2frame, textvariable = nPerc1.value,
        width = 5)
    Det1.2.2lab <- tklabel(Det1.2.2frame, text = "Spline ratio:  ")
    tkpack(Det1.2.2lab, nPerc1.entry, side = "left", anchor = "w")
    p1.entry <- tkentry(Det1.2.3frame, textvariable = p1.value,
        width = 5)
    Det1.2.3lab <- tklabel(Det1.2.3frame, text = "Value of p:    ")
    tkpack(Det1.2.3lab, p1.entry, side = "left", anchor = "w")
    tkpack(Det1.2.1frame, Det1.2.2frame, Det1.2.3frame, side = "top")
    n2.entry <- tkentry(Det2.2.1frame, textvariable = n2.value,
        width = 5)
    Det2.2.1lab <- tklabel(Det2.2.1frame, text = "Spline length:")
    if (TwoSteps)
        tkpack(Det2.2.1lab, n2.entry, side = "left")
    nPerc2.entry <- tkentry(Det2.2.2frame, textvariable = nPerc2.value,
        width = 5)
    Det2.2.2lab <- tklabel(Det2.2.2frame, text = "Spline ratio:  ")
    if (TwoSteps)
        tkpack(Det2.2.2lab, nPerc2.entry, side = "left", anchor = "w")
    p2.entry <- tkentry(Det2.2.3frame, textvariable = p2.value,
        width = 5)
    Det2.2.3lab <- tklabel(Det2.2.3frame, text = "Value of p:    ")
    if (TwoSteps)
        tkpack(Det2.2.3lab, p2.entry, side = "left", anchor = "w")
    if (TwoSteps)
        tkpack(Det2.2.1frame, Det2.2.2frame, Det2.2.3frame, side = "top")
    tkpack(Det1.1frame, Det1.2frame, side = "left", expand = 1,
        fill = "x")
    if (TwoSteps)
        tkpack(Det2.1frame, Det2.2frame, side = "left", expand = 1,
            fill = "x")
    if (TwoSteps) {
        tkpack(Det1frame, Det2frame, side = "left", expand = 1,
            fill = "x")
    }
    else {
        tkpack(Det1frame, side = "left", expand = 1, fill = "x")
    }
    tkpack(top_detrending_frame2, fill = "x")
    OnOk = function() {
        makeFirstDetrending <<- TRUE
        method1 <<- tclvalue(method1.value)
        n1 <<- toNumber(tclvalue(n1.value))
        nPerc1 <<- toNumber(tclvalue(nPerc1.value))
        p1 <<- toNumber(tclvalue(p1.value))
        first.detrending.method <<- GetDetrendMethod(method1,
            n1, nPerc1, p1)
        makeSecondDetrending <<- TRUE
        method2 <<- tclvalue(method2.value)
        n2 <<- toNumber(tclvalue(n2.value))
        nPerc2 <<- toNumber(tclvalue(nPerc2.value))
        p2 <<- toNumber(tclvalue(p2.value))
        second.detrending.method <<- GetDetrendMethod(method2,
            n2, nPerc2, p2)
        interactive.detrend <<- as.logic(tclvalue(interactive.detrend.value))
        tclvalue(done) <- 1
    }
    top_detrending_frame5 <- tkframe(top_detrending, relief = "groove",
        borderwidth = 2)
    interactive.detrend.value <- tclVar(interactive.detrend)
    interactive.detrend.cbut <- tkcheckbutton(top_detrending_frame5,
        text = "Interactive detrending", variable = interactive.detrend.value)
    tkpack(interactive.detrend.cbut, side = "left")
    Cancel.but <- tkbutton(top_detrending_frame5, text = "Cancel",
        command = function() tkdestroy(top_detrending))
    Ok.but <- tkbutton(top_detrending_frame5, text = "  Ok  ",
        command = OnOk)
    tkpack(Ok.but, Cancel.but, side = "right", expand = "FALSE",
        fill = "y")
    tkpack(top_detrending_frame5, fill = "x")
    tkbind(top_detrending, "<Destroy>", function() tclvalue(done) <- 2)
    tkbind(top_detrending, "<KeyPress-Return>", function() OnOk())
    tkbind(top_detrending, "<KeyPress-Escape>", function() tclvalue(done) <- 2)
    tkfocus(top_detrending)
    tkwait.variable(done)
    tkgrab.release(top_detrending)
    if (tclvalue(done) == "2") {
        tkdestroy(top_detrending)
    }
    if (tclvalue(done) == "1") {
        .input = tclvalue(filenamevar)
        .output = tclvalue(tabnamevar)
        tkdestroy(top_detrending)
        interactiveDETRENDING(input = .input, output = .output,
            TwoSteps = TwoSteps)
    }
}
