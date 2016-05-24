"selectDialog" <-
function (title = "Selection entry", question = "Select item", 
    itemNames, top = NULL, returnValOnCancel = "ID_CANCEL", do.grab = FALSE) 
{
    "subSelectDialog" <- function(popup, itemNames, title = title, 
        subtitle = question) {
        tkwm.title(popup, title)
        scr <- tkscrollbar(popup, repeatinterval = 5, 
                           command = function(...) tkyview(tl, ...))
        tl <- tklistbox(popup, height = 4, selectmode = "single", 
                        yscrollcommand = function(...) tkset(scr, ...), 
                        background = "white")
        tkgrid(tklabel(popup, text = subtitle))
        tkgrid(tl, scr)
        tkgrid.configure(scr, rowspan = 4, sticky = "nsw")
        for (i in (1:length(itemNames))) 
            tkinsert(tl, "end", itemNames[i])
        tkselection.set(tl, 0)
        return(tl)
    }
    popup <- tktoplevel()
    tkwm.deiconify(popup)
    if (do.grab) 
        tkgrab.set(popup)
    tkfocus(popup)
    ReturnVal <- returnValOnCancel
    tl <- subSelectDialog(popup, itemNames, title = title, subtitle = question)
    "OnOK" <- function() {
        ReturnVal <<- as.numeric(tkcurselection(tl)) + 1
        tkgrab.release(popup)
        tkdestroy(popup)
        if (!is.null(top)) 
            tkfocus(top)
    }
    OK.but <- tkbutton(popup, text = "   OK   ", command = OnOK)
    tkgrid(OK.but)
    tkwait.window(popup)
    return(ReturnVal)
}
