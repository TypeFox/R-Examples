tklist.modal <- function(title,elements0,returnValOnCancel="ID_CANCEL", selectmode = "single")
{
  dlg <- tktoplevel()
  tkwm.title(dlg, "Question")
  tkwm.deiconify(dlg)
  tkgrab.set(dlg)
  tkfocus(dlg)

  scr <- tkscrollbar(dlg, repeatinterval=5,
                     command=function(...)tkyview(tl,...))

  width0 <- max(nchar(elements0)) + 10

  tl<-tklistbox(dlg,height=4,
                selectmode=selectmode,
                yscrollcommand=function(...)tkset(scr,...),
                background="white", width = width0)

  tkgrid(tklabel(dlg,text=title))
  tkgrid(tl,scr)

  tkgrid.configure(scr,rowspan=4,sticky="nsw")
  for (i in 1:length(elements0)) {
    tkinsert(tl,"end",elements0[i])
  }

  ReturnVal <- returnValOnCancel
  onOK <- function()
  {
    ReturnVal <<- elements0[as.numeric(tkcurselection(tl))+1]
    tkgrab.release(dlg)
    tkdestroy(dlg)
#    tkfocus(ttMain)
  }
  onCancel <- function()
  {
    ReturnVal <<- returnValOnCancel
    tkgrab.release(dlg)
    tkdestroy(dlg)
#    tkfocus(ttMain)
  }
  OK.but     <-tkbutton(dlg,text="OK",command=onOK)
  tkgrid(OK.but)
  tkgrid(tklabel(dlg,text="    "))

  tkfocus(dlg)
  tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg)})
  tkwait.window(dlg)

  return(ReturnVal)

}
