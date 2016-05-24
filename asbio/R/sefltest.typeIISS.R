selftest.typeIISS.tck1 <-function(){

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

tkwm.title(tt, "Type II and III SS")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "We would use type II or III SS whenever..."), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    We have an unbalanced one way ANOVA format."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    We have a quantitative predictor."),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    We have a balanced multiway ANOVA format."),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    We have either an unbalanced multiway ANOVA format, or some other multiple X format with at least one quantitative X."),rb4, sticky = "w")

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