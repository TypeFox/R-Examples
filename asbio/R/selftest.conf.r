selftest.conf.tck1<-function(){

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
tkwm.title(tt, "Self-test, Confidence Intervals")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "A 95% confidence interval for \u03bc is calculated with the following confidence limits: [24, 28].  What is the correct interpretation of this result?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    There is a probability of 0.95 that \u03bc will be between 24 and 28."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    We have 95% confidence that \u03bc will be between 24 and 28."),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    We have 95% confidence that the sample mean will be between 24 and 28."),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    There is a probability of 0.95 that the sample mean will be between 24 and 28."),rb4, sticky = "w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Incorrect. The interval has been calculated.  The interval will either contain mu, or it will not, i.e. the probability that the interval contains mu is either 0 or 1.", icon = "error")
    if (rbVal=="b")
    	tkmessageBox(message="Correct.")
    if (rbVal=="c")
    	tkmessageBox(message="Incorrect. The confidence interval for mu is calculated using the sample mean, and as a result the sample mean will always be in the confidence interval interval for mu", icon = "error")
if (rbVal=="d")
    	tkmessageBox(message="Incorrect. This answer is incorrect for several reasons.  For one thing the confidence interval for mu is calculated using the sample mean, and as a result the sample mean will always be in the confidence interval interval for mu", icon = "error")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.conf.tck2())),
tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}



selftest.conf.tck2<-function(){

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

tkwm.title(tt, "Self-test, Confidence Intervals")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "If one hundred confidence intervals for \u03bc are calculated then, on average, a proportion of 1 - \u03b1 of these will contain \u03bc."), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    TRUE.  "),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    FALSE.  "),rb2, sticky = "w")
tkgrid(tklabel(tt,text=""))

OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
     if (rbVal=="a")
    	tkmessageBox(message="Correct.")
    if (rbVal=="b")
    	tkmessageBox(message="Inorrect.  The definition given is correct and is often used to help conceptualize confidence intervals.", icon = "error")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=function()tkmessageBox(message="No further questions.")),tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}

