"modalDialog" <-
function(title,itemlist,entryInit,entryWidth=20,
                        returnValOnCancel="ID_CANCEL")
{
  dlg <- tktoplevel()
  tkwm.deiconify(dlg)
  tkgrab.set(dlg)
  tkfocus(dlg)
  tkwm.title(dlg,title)
  
#  scr <- tkscrollbar(tt, repeatinterval=5,
#                     command=function(...)tkyview(tl,...))

  textEntryWidget <- rep(NA, times = length(itemlist))
  textEntryVarTcl <- rep(NA, times = length(itemlist))

  textEntryVarTcl <- list()
  textEntryWidget <- list()

  tkgrid(tklabel(dlg,text="       "))
  for (i in 1:length(itemlist)) {
    textEntryVarTcl[[i]] <- tclVar(paste(entryInit[i]))
    textEntryWidget[[i]] <- tkentry(dlg,width=paste(entryWidth),
                               textvariable=textEntryVarTcl[[i]])
    tkgrid(tklabel(dlg,text=itemlist[i]),textEntryWidget[[i]])
  }
  tkgrid(tklabel(dlg,text="       "))

  ReturnVal <- returnValOnCancel

  onOK <- function()
  {
    ReturnVal <<- rep(NA, times = length(itemlist))
    for (i in 1:length(itemlist)) {
      ReturnVal[i] <<- tclvalue(textEntryVarTcl[[i]])
    }
    tkgrab.release(dlg)
    tkdestroy(dlg)
#    tkfocus(tt)
  }
  onCancel <- function()
  {
    ReturnVal <<- returnValOnCancel
    tkgrab.release(dlg)
    tkdestroy(dlg)
#    tkfocus(tt)
  }
  OK.but     <-tkbutton(dlg,text="   OK   ",command=onOK)
  Cancel.but <-tkbutton(dlg,text=" Cancel ",command=onCancel)
  tkgrid(OK.but,Cancel.but)
  tkgrid(tklabel(dlg,text="    "))

  tkfocus(dlg)
  tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg)})
  tkwait.window(dlg)

  return(ReturnVal)

}






