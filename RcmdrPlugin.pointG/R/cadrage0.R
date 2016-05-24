cadrage0<-function(){
defaults <- list(initial.variables = NULL)
    dialog.values <- getDialog("cadrage0", defaults)
    initializeDialog(title=gettextRcmdr("Donnees de cadrage"))
    variablesFrame <- tkframe(top)
    subsetBox()
    onOK <- function(){
        percents <- as.character(tclvalue(percentsVariable))
        chisq <- tclvalue(chisqTestVariable)
        chisqComp <- tclvalue(chisqComponentsVariable)
        expected <- tclvalue(expFreqVariable)
        fisher <- tclvalue(fisherTestVariable)
        subset <- tclvalue(subsetVariable)
        closeDialog()
        if (chisq == 1) {shell.exec("http://www.insee.fr/fr/themes/tableau.asp?reg_id=0&ref_id=ccc")}

                      if (expected == 1) {shell.exec("http://www.insee.fr/fr/themes/tableau.asp?reg_id=0&ref_id=NATCCV05109")}
            if (chisqComp == 1) {shell.exec("http://www.insee.fr/fr/publics/communication/recensement/particuliers/doc/fiche-PCS.pdf")}
                
        if (fisher == 1) {{shell.exec("http://www.acteursdusport.fr/491-pratique-sportive-et-consommation-sportive.htm")}}
        tkfocus(CommanderWindow())
        }
    OKCancelHelp(helpSubject="shell.exec")
    radioButtons(name="percents",
     buttons=c("rowPercents", "columnPercents", "totalPercents", 
	"nonePercents"),
       values=c("row", "column", "total", "none"), initialValue="none",
        labels=gettextRcmdr(c("Row percentages", "Column percentages", "Percentages of total", "No percentages")), title=gettextRcmdr("Compute Percentages"))
    
checkBoxes(frame="testsFrame", boxes=c("chisqTest", "chisqComponents", "expFreq", "fisherTest"), initialValues=c("0", "0", "0", "0"),
        labels=gettextRcmdr(c("Pyramide des ages", "Professions",
            "Consommation", "Sport")))
    tkgrid(variablesFrame, sticky="w")
    tkgrid(labelRcmdr(top, text=gettextRcmdr("Informations sur..."), fg="blue"), sticky="w")
    tkgrid(testsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=6, columns=1)
}
