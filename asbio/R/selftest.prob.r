

selftest.prob.tck1<-function(){

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
tkwm.title(tt, "Self-test, probability")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "Given the following situation: \n P(A) = 0.2, P(B) = 0.4, P(A\u2229B) = 0.08, \n\n are A and B disjoint?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    TRUE.  "),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    FALSE.  "),rb2, sticky = "w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Incorrect. The probability of the union of A and B \u2260 0", icon = "error")
    if (rbVal=="b")
    	tkmessageBox(message="Correct.")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.prob.tck2())),
tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}



selftest.prob.tck2<-function(){

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
tkwm.title(tt, "Self-test, probability")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "Given the following situation: \n P(A) = 0.2, P(B) = 0.4, P(A\u2229B) = 0.08, \n\n are A and B independent?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    TRUE.  "),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    FALSE.  "),rb2, sticky = "w")
tkgrid(tklabel(tt,text=""))

OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
     if (rbVal=="a")
    	tkmessageBox(message="Correct. The probability of the union of A and B = P(A)P(B)")
    if (rbVal=="b")
    	tkmessageBox(message="Inorrect.  A and B are independent because the probability of the union of A and B  = P(A)P(B)", icon = "error")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.prob.tck3())),tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}

selftest.prob.tck3<-function(){
invisible(tclServiceMode(FALSE))

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
tkwm.title(tt, "Self-test, probability")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "Given the following situation: \n P(A) = 0.2, P(B) = 0.4, P(A|B) = 0.2, \n\n are A and B independent?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    TRUE.  "),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    FALSE.  "),rb2, sticky = "w")
tkgrid(tklabel(tt,text=""))

OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
     if (rbVal=="a")
    	tkmessageBox(message="Correct. P(A|B) = P(A)")
    if (rbVal=="b")
    	tkmessageBox(message="Inorrect.  A and B are independent because P(A|B) = P(A)", icon = "error")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=function()tkmessageBox(message="No further questions.")),tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}
