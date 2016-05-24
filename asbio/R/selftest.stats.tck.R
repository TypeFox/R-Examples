selftest.stats.tck1<-function(){

tclServiceMode(FALSE)
tt <- tktoplevel()
 have_ttk <- as.character(tcl("info", "tclversion")) >=
            "8.5"
 if (have_ttk) {
            tkbutton <- ttkbutton
            tkcheckbutton <- ttkcheckbutton
            tkentry <- ttkentry
            tkframe <- ttkframe
            tklabel <- ttklabel
            tkradiobutton <- ttkradiobutton
        }

rb1 <- tkradiobutton(tt)
rb2 <- tkradiobutton(tt)
rb3 <- tkradiobutton(tt)
rb4 <- tkradiobutton(tt)
rb5 <- tkradiobutton(tt)
rbValue <- tclVar("")
tkconfigure(rb1, variable=rbValue,value="a")
tkconfigure(rb2, variable=rbValue,value="b")
tkconfigure(rb3, variable=rbValue,value="c")
tkconfigure(rb4, variable=rbValue,value="d")
tkconfigure(rb5, variable=rbValue,value="e")
tkwm.title(tt, "Statistics and parameters")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "The efficacy of a statistical estimator can be be quantified using its..."), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    stress.          "),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    bias.              "),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    efficiency.        "),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    consistency.        "),rb4, sticky = "w")
tkgrid(tklabel(tt,text="  e.    all but (a) are correct.        "),rb5, sticky = "w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Incorrect.", icon = "error")
    if (rbVal=="b")
    	tkmessageBox(message="Partially true.", icon = "error")
    if (rbVal=="c")
    	tkmessageBox(message="Partially true.", icon = "error")
    if (rbVal=="d")
    	tkmessageBox(message="Partially true.", icon = "error")
	if (rbVal=="e")
    	tkmessageBox(message="Correct.")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.stats.tck2())),
tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}



selftest.stats.tck2<-function(){

tclServiceMode(FALSE)
tt <- tktoplevel()
 have_ttk <- as.character(tcl("info", "tclversion")) >=
            "8.5"
 if (have_ttk) {
            tkbutton <- ttkbutton
            tkcheckbutton <- ttkcheckbutton
            tkentry <- ttkentry
            tkframe <- ttkframe
            tklabel <- ttklabel
            tkradiobutton <- ttkradiobutton
        }

rb1 <- tkradiobutton(tt)
rb2 <- tkradiobutton(tt)
rb3 <- tkradiobutton(tt)
rb4 <- tkradiobutton(tt)
rbValue <- tclVar("")
tkconfigure(rb1,variable=rbValue,value="a")
tkconfigure(rb2,variable=rbValue,value="b")


tkwm.title(tt, "Statistics and parameters")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "(T or F) If a sample includes the entire population then statistics serve a descriptive instead of an inferential role."), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    TRUE."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    FALSE."),rb2, sticky = "w")
tkgrid(tklabel(tt,text=""))

OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Correct.")
    if (rbVal=="b")
    	tkmessageBox(message="Incorrect.", icon = "error")
 }

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=function()tkmessageBox(message="No further questions.")),tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}

