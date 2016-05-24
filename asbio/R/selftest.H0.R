selftest.H0.tck1 <-function(){

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
tkconfigure(rb3,variable=rbValue,value="c")
tkconfigure(rb4,variable=rbValue,value="d")
tkwm.title(tt, "Self-test, H\u2080 hypothesis tests")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "I reject H\u2080 but H\u2080 is true.  This is a...?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    type I error."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    type II error."),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    type III error."),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    none of the above."),rb4, sticky = "w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Correct.")
	if (rbVal=="b")
    	tkmessageBox(message="Incorrect.", icon = "error")
    if (rbVal=="c")
    	tkmessageBox(message="Incorrect.", icon = "error")
	if (rbVal=="d")
    	tkmessageBox(message="Incorrect.", icon = "error")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.H0.tck2())),
tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}



selftest.H0.tck2 <-function(){

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
tkconfigure(rb3,variable=rbValue,value="c")
tkconfigure(rb4,variable=rbValue,value="d")

tkwm.title(tt, "Self-test, H\u2080 hypothesis tests")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "We suspect that a soil N treatment will increase mean plant biomass compared to a control.  What would our alternative hypothesis be?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    \u03BCN < \u03BCControl."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    \u03BCN > \u03BCControl."),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    \u03BCN = \u03BCControl."),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    \u03BCN \u2260 \u03BCControl."),rb4, sticky = "w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()

{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Incorrect.", icon = "error")
	if (rbVal=="b")
    	tkmessageBox(message="Correct.")
    if (rbVal=="c")
    	tkmessageBox(message="Incorrect.", icon = "error")
	if (rbVal=="d")
    	tkmessageBox(message="Incorrect.", icon = "error")
}


OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=function()tkmessageBox(message="No further questions.")),tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}

