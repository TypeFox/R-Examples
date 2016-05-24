subsetBoxDoE <- defmacro(window = tab1, expr={
    ## default: all
    if (!exists("constraintVar", where="RcmdrEnv")) 
        putRcmdr("constraintVar", tclVar(gettextRcmdr("<all candidate data set rows eligible>")))
    
    subsetFrame <- ttklabelframe(window, text = gettextRcmdr("Constraints for design runs"))
    subsetEntry <- ttkentry(subsetFrame, width = "80", textvariable = constraintVar)
    subsetScroll <- ttkscrollbar(subsetFrame, orient = "horizontal",
        command = function(...) tkxview(subsetEntry, ...))
    tkconfigure(subsetEntry, xscrollcommand = function(...) tkset(subsetScroll,
        ...))
    tkgrid(tklabel(subsetFrame, text=gettextRcmdr("For restricting the usable rows from the candidate data set,")), sticky = "w")
    tkgrid(tklabel(subsetFrame, text=gettextRcmdr("type in a logical expression using variable names and constants.")), sticky = "w")
    tkgrid(subsetEntry, sticky = "w")
    tkgrid(subsetScroll, sticky = "ew")
    helpsyntaxButton <- buttonRcmdr(subsetFrame, text = gettextRcmdr("Syntax Help"), 
        foreground = "darkgreen", command = onHelpSyntax, 
        default = "normal", borderwidth = 3)

    tkgrid(tklabel(subsetFrame, text=gettextRcmdr("Syntax reminder: &, |, !, ==, >=, >, <=, <, %in%")), helpsyntaxButton, sticky = "w")
    putRcmdr("subsetFrame", subsetFrame)
}
)
onHelpSyntax <- function(){
     if (GrabFocus() && .Platform$OS.type != "windows") 
            tkgrab.release(topdes2)     
     print(help("Syntax"))

}