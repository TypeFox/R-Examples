selftest.ANOVAsiminf.tck1 <-function(){

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

tkwm.title(tt, "Simultaneous inference")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "What is the set of comparisons for which Scheff\u00E9's procedure will control family-wise type I error?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    All possible comparisons."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    All possible pairwise comparisons."),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    All possible comparisons to a control."),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    None of the above."),rb4, sticky = "w")

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
tkgrid(tkbutton(tt,text="Next question",command=function()tkmessageBox(message="No further questions.")),tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}