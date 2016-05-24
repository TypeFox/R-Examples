choixvarfunc <- function(title, question, liste)
  {
    tt2 <- tktoplevel()
    scr <- tkscrollbar(tt2, repeatinterval = 5, command = function(...) tkyview(lstbox,...))
    lstbox <- tklistbox(tt2, height = 4, selectmode = "single",
    yscrollcommand = function(...) tkset(scr, ...), background = "white")
    tkgrid(tklabel(tt2, text = question))
    tkfocus(tt2)
    tkwm.title(tt2, title)
    varChoice <- ""
    tkgrid(lstbox, scr)
    tkgrid.configure(scr, rowspan = 4, sticky = "nsw")
    var <- liste

    for(i in 1:length(var))
     {tkinsert(lstbox, "end", var[i])
     }

    OnOK <- function()
     {
      varChoice <<- var[as.numeric(tkcurselection(lstbox))+1]
      tkdestroy(tt2)
     }

    OK.but <- tkbutton(tt2, text = "   OK   ", command = OnOK)
    tkgrid(OK.but)
    tkbind(tt2, "<Destroy>", function() {tkdestroy(tt2)})
    tkfocus(tt2)
    tkwait.window(tt2)
    return(varChoice)
 }
