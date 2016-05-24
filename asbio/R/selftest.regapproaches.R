selftest.regapproaches.tck1 <-function(){

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
tkconfigure(rb1,variable=rbValue,value="a")
tkconfigure(rb2,variable=rbValue,value="b")
tkconfigure(rb3,variable=rbValue,value="c")
tkconfigure(rb4,variable=rbValue,value="d")
tkconfigure(rb5,variable=rbValue,value="e")
tkwm.title(tt, "Regression approaches")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "We use model II regression..."), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    to address outliers in linear regression."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    to address curvilinear associations within the context of a linear model."),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    to address heteroscedasticity in a general linear model."),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    to address cases in which levels in X are not fixed, or are measured with error."),rb4, sticky = "w")
tkgrid(tklabel(tt,text="  e.    to address non-quantitative, discrete, and/or strictly bounded response variables."),rb5, sticky = "w")
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
    if (rbVal=="e")
    	tkmessageBox(message="Incorrect.", icon = "error")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.regapproaches.tck2())),
tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}



selftest.regapproaches.tck2 <-function(){

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
tkconfigure(rb1,variable=rbValue,value="a")
tkconfigure(rb2,variable=rbValue,value="b")
tkconfigure(rb3,variable=rbValue,value="c")
tkconfigure(rb4,variable=rbValue,value="d")
tkconfigure(rb5,variable=rbValue,value="e")
tkwm.title(tt, "Regression approaches")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "We use weighted least squares..."), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    to address outliers in linear regression."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    to address curvilinear associations within the context of a linear model."),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    to address heteroscedasticity in a general linear model."),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    to address cases in which levels in X are not fixed, or are measured with error."),rb4, sticky = "w")
tkgrid(tklabel(tt,text="  e.    to address non-quantitative, discrete, and/or strictly bounded response variables."),rb5, sticky = "w")
tkgrid(tklabel(tt,text=""))
OnOK <- function()
{
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="a")
     	tkmessageBox(message="Incorrect.", icon = "error")
    if (rbVal=="b")
    	tkmessageBox(message="Incorrect.", icon = "error")
    if (rbVal=="c")
    	tkmessageBox(message="Correct.")
    if (rbVal=="d")
      tkmessageBox(message="Incorrect.", icon = "error")
    if (rbVal=="e")
    	tkmessageBox(message="Incorrect.", icon = "error")
}



OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.regapproaches.tck3())),
tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}

selftest.regapproaches.tck3 <-function(){

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
tkconfigure(rb1,variable=rbValue,value="a")
tkconfigure(rb2,variable=rbValue,value="b")
tkconfigure(rb3,variable=rbValue,value="c")
tkconfigure(rb4,variable=rbValue,value="d")
tkconfigure(rb5,variable=rbValue,value="e")
tkwm.title(tt, "Regression approaches")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "We use polynomial regression..."), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    to address outliers in linear regression."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    to address curvilinear associations within the context of a linear model."),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    to address heteroscedasticity in a general linear model."),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    to address cases in which levels in X are not fixed, or are measured with error."),rb4, sticky = "w")
tkgrid(tklabel(tt,text="  e.    to address non-quantitative, discrete, and/or strictly bounded response variables."),rb5, sticky = "w")
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
    if (rbVal=="e")
    	tkmessageBox(message="Incorrect.", icon = "error")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=substitute(selftest.regapproaches.tck4())),
tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}


selftest.regapproaches.tck4 <-function(){

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
tkconfigure(rb1,variable=rbValue,value="a")
tkconfigure(rb2,variable=rbValue,value="b")
tkconfigure(rb3,variable=rbValue,value="c")
tkconfigure(rb4,variable=rbValue,value="d")
tkconfigure(rb5,variable=rbValue,value="e")
tkwm.title(tt, "Regression approaches")
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt, text = "We can use the framework of GLMs..."), columnspan = 2)
tkgrid(tklabel(tt,text=""))
tkgrid(tklabel(tt,text="  a.    to address outliers in linear regression."),rb1, sticky = "w")
tkgrid(tklabel(tt,text="  b.    to address curvilinear associations within the context of a linear model."),rb2, sticky = "w")
tkgrid(tklabel(tt,text="  c.    to address heteroscedasticity in a general linear model."),rb3, sticky = "w")
tkgrid(tklabel(tt,text="  d.    to address cases in which levels in X are not fixed, or are measured with error."),rb4, sticky = "w")
tkgrid(tklabel(tt,text="  e.    to address non-quantitative, discrete, and/or strictly bounded response variables."),rb5, sticky = "w")
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
      tkmessageBox(message="Incorrect.", icon = "error")
    if (rbVal=="e")
    	tkmessageBox(message="Correct.")
}

OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but,columnspan=2)
tkgrid(tklabel(tt,text=""))
tkgrid(tkbutton(tt,text="Next question",command=function()tkmessageBox(message="No further questions.")),tkbutton(tt,text="Exit",command=function()tkdestroy(tt)),sticky ="w")
invisible(tclServiceMode(TRUE))
}