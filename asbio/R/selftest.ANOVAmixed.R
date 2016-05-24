selftest.ANOVAmixed.tck1 <-function(){

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

tkwm.title(tt, "Mixed effects")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "Why do we use REML to estimate variances in random and mixed effect models?"), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    ML estimates are biased low."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    MOM estimates can be negative."),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    Lack of balance in multifactor designs prevents MOM estimates."),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    All of the above."),rb4, sticky = "w")

tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
      tkmessageBox(message="Incorrect.", icon = "error")
    if (rbVal=="b")
     	tkmessageBox(message="Incorrect.", icon = "error")
    if (rbVal=="c")
    	tkmessageBox(message="Incorrect.", icon = "error")
    if (rbVal=="d")
      tkmessageBox(message="Correct.")
    }


OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=function()tkmessageBox(message="No further questions.")),tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}