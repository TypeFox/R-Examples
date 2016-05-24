

selftest.pdfs.tck1<-function(){

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
tkwm.title(tt, "Self-test, the normal pdf")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "What does the normal pdf describe?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    A continuous random variable, e.g. height or weight.  "),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    The number of successes, with a lower limit of 0, and an unbounded upper limit.  "),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    The number of successes, with a lower limit of 0, and an upper bound of n.  "),rb3, sticky ="w")
tkgrid(tklabel(tt,text="  d.    A discrete random variable, e.g. height or weight.   "),rb4,sticky="w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Correct.")
    if (rbVal=="b")
    	tkmessageBox(message="Incorrect.  The normal distribution describes continuous random variables.",icon="error")
    if (rbVal=="c")
    	tkmessageBox(message="Incorrect.  The normal distribution describes continuous random variables.",icon="error")
    if (rbVal=="d")
    	tkmessageBox(message="Incorrect.  Height and weight are continuous.",icon="error")  	
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.pdfs.tck2())),
tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}



selftest.pdfs.tck2<-function(){

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
tkwm.title(tt, "Self-test, the binomial pdf")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "What does the binomial pdf describe?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    A continuous random variable, e.g. height or weight.  "),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    The number of successes, with a lower limit of 0, and an unbounded upper limit.  "),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    The number of successes, with a lower limit of 0, and an upper bound of n.  "),rb3, sticky ="w")
tkgrid(tklabel(tt,text="  d.    A discrete random variable, e.g. height or weight.   "),rb4,sticky="w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Incorrect.  The binomial pdf describes discrete random variables.",icon="error")
    if (rbVal=="b")
    	tkmessageBox(message="Incorrect.  The upper bound will equal the number of trials, n",icon="error")
    if (rbVal=="c")
    	tkmessageBox(message="Correct.")
    if (rbVal=="d")
    	tkmessageBox(message="Incorrect.  The binomial pdf describes discrete random variables (height and weight are not discrete).",icon="error")  	
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.pdfs.tck3())),
tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}

selftest.pdfs.tck3<-function(){

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
tkwm.title(tt, "Self-test, the binomial pdf")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "What does the Poisson pdf describe?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    A continuous random variable, e.g. height or weight.  "),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    The number of successes, with a lower limit of 0, and an unbounded upper limit.  "),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    The number of successes, with a lower limit of 0, and an upper bound of n.  "),rb3, sticky ="w")
tkgrid(tklabel(tt,text="  d.    A discrete random variable, e.g. height or weight.   "),rb4,sticky="w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Incorrect.  The Poisson pdf describes discrete random variables.",icon="error")
    if (rbVal=="b")
    	tkmessageBox(message="Correct.")
    if (rbVal=="c")
    	tkmessageBox(message="Inorrect.  This describes the binomial pdf", icon="error")
    if (rbVal=="d")
    	tkmessageBox(message="Incorrect.  The Poisson pdf describes discrete random variables (height and weight are not discrete).",icon="error")  	
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=function()tkmessageBox(message="No further questions.")),tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}