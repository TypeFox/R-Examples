

selftest.sampd.tck1<-function(){

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
tkwm.title(tt, "Self-test, Sampling Distributions")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "If the parent distribution is \u0059 ~ N(\u03bc, \u03c3\u00b2), \nthen the sampling distribution of \u0232 will be: "), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    \u0232 ~ N(\u0232, \u03c3\u00b2).          "),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    \u0232 ~ N(\u03bc, \u03c3/n).              "),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    \u0232 ~ N(\u03bc/n, \u03c3\u00b2).        "),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    \u0232 ~ N(\u03bc, \u03c3\u00b2/n).        "),rb4, sticky = "w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Incorrect.", icon = "error")
    if (rbVal=="b")
    	tkmessageBox(message="Inorrect.", icon = "error")
    if (rbVal=="c")
    	tkmessageBox(message="Incorrect.", icon = "error")
if (rbVal=="d")
    	tkmessageBox(message="Correct.")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.sampd.tck2())),
tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}



selftest.sampd.tck2<-function(){

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

tkwm.title(tt, "Self-test, Sampling Distributions")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "Which of the following describes the Central Limit Theorem?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    If Y ~ N(\u03bc, \u03c3\u00b2), then \u0232 ~ N(\u03bc/n, \u03c3\u00b2/n)."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    If the sample size from Y is large, then the variance for the sampling distribution of \u0232 will be small."),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    Even if Y is not normal, if the sample size from Y is large, then \u0232 will be approximately normal."),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    None of the above."),rb4, sticky = "w")
tkgrid(tklabel(tt,text=""))

OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Incorrect", icon = "error")
    if (rbVal=="b")
    	tkmessageBox(message="This does occur, but it can be more concisely attributed to the Law of Large Numbers.", icon = "error")
    if (rbVal=="c")
    	tkmessageBox(message="This is the most correct answer.")
    if (rbVal=="d")
    	tkmessageBox(message="Incorrect.", icon = "error")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=function()tkmessageBox(message="No further questions.")),tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}

