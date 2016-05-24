selftest.se.tck1<-function(){

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
tkwm.title(tt, "Self-test, assignment of treatments")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "Random assignment of experimental units to treatments will provide..."), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    causal proof regarding the effect of the treatment.  "),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    causal evidence regarding the effect of the treatment.  "),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    inference to the sampled population.  "),rb3, sticky ="w")
tkgrid(tklabel(tt,text="  d.    inference to the sample.   "),rb4,sticky="w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Incorrect.  Proof of causality is not possible.",icon="error")
    if (rbVal=="b")
    	tkmessageBox(message="Correct")
    if (rbVal=="c")
    	tkmessageBox(message="Incorrect.  The method of sampling (not the method of treatment assignment) will determine inference to the population.",icon="error")
    if (rbVal=="d")
    	tkmessageBox(message="Incorrect.  The method of sampling (not the method of treatment assignment) will determine inference to the population.",icon="error")  	
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.se.tck2())),
tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}



selftest.se.tck2<-function(){

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
tkwm.title(tt, "Self-test, selection of samples")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "Random selection of samples from the population will provide..."), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    causal proof regarding the effect of the treatment.  "),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    causal evidence regarding the effect of the treatment.  "),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    inference to the sampled population.  "),rb3, sticky ="w")
tkgrid(tklabel(tt,text="  d.    inference to the sample.   "),rb4,sticky="w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
    	tkmessageBox(message="Incorrect.  First of all, proof of causality is not possible.  Secondly, The method of treatment assignment (not the method of sampling) will determine whether evidence will be causal. ",icon="error")
    if (rbVal=="b")
    	tkmessageBox(message="Incorrect.  The method of treatment assignment (not the method of sampling) will determine whether evidence will be causal.",icon="error")
    if (rbVal=="c")
    	tkmessageBox(message="Correct") 
    if (rbVal=="d")
    	tkmessageBox(message="Incorrect.  Random sampling from a population will allow inference to the population (not just the sample).",icon="error")  	
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=function()tkmessageBox(message="No further questions.")),tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}