twoWayTable0<-function(){
Library("ade4")
Library("abind")
Library("vcd")
        defaults <- list(initial.row=NULL, initial.column=NULL
, initial.tableau="ligne"
, initial.graf="mosaic"
, initial.test="x2"
, initial.effet="V"
, initial.local="chi2"
			)
        dialog.values <- getDialog("twoWayTable0", defaults)
        initializeDialog(title=gettextRcmdr("Tableau de contingence"))
        variablesFrame <- tkframe(top)
        .factors <- Factors()
        rowBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Row variable (pick one)"),
                        initialSelection=varPosn(dialog.values$initial.row, "factor"))
        columnBox <- variableListBox(variablesFrame, .factors, title=gettextRcmdr("Column variable (pick one)"),
                        initialSelection=varPosn(dialog.values$initial.column, "factor"))


        onOK <- function(){
                row <- getSelection(rowBox)
                column <- getSelection(columnBox)
           tableau <- as.character(tclvalue(tableauVariable))
                graf <- as.character(tclvalue(grafVariable))
                test <- as.character(tclvalue(testVariable))
                effet <- as.character(tclvalue(effetVariable))
                local <- as.character(tclvalue(localVariable))

                
               putDialog("twoWayTable0", list(initial.row=row, 
initial.column=column,
initial.tableau=tableau,
     initial.graf=graf,
     initial.test=test,
     initial.effet=effet,
     initial.local=local))

                if (length(row) == 0 || length(column) == 0){
                        errorCondition(recall=twoWayTable0, message=gettextRcmdr("You must select two variables."))
                        return()
                }
                if (row == column) {
                        errorCondition(recall=twoWayTable0, message=gettextRcmdr("Row and column variables are the same."))
                        return()
                }
                closeDialog()
                

                command <- paste("xtabs(~", row, "+", column, ", data=", ActiveDataSet(),	")", sep="")
        doItAndPrint(paste(".Table <- ", command, sep=""))
        doItAndPrint(".Table")
    
   
 if (tableau == "ligne") doItAndPrint("rowPercents(.Table)") 
 if (tableau == "colonne") doItAndPrint("colPercents(.Table)") 
 if (tableau == "total") doItAndPrint("totPercents(.Table)") 
                
                if (graf == "mosaic") {
command0<-paste("mosaicplot(.Table,main=\"\",col=brewer.pal(max(length(dimnames(.Table)[[2]]),3),\"Set1\"))",sep="")
doItAndPrint(command0)
}
                if (graf == "bar"){
command0<-paste("barplot(t(sweep(.Table,1,apply(.Table,1,sum),\"/\")),beside=TRUE,col=brewer.pal(max(length(dimnames(.Table)[[2]]),3),\"Set1\")[1:(length(dimnames(.Table)[[2]]))],xlab=names(dimnames(.Table))[1],legend.text=TRUE,ylab=\"",paste("pourcentage(",column,")",sep=""),"\")",sep="")
doItAndPrint(command0)}

                if (graf == "assoc"){
command0<-"assoc(.Table,shade=TRUE)"
doItAndPrint(command0)}

                if (graf == "factomap"){
command0<-"mapAFC(.Table)"
doItAndPrint(command0)}


                
 if (test == "x2") doItAndPrint("chisq.test(.Table)") 
 if (test == "fisher") doItAndPrint("fisher.test(.Table)") 

if (effet == "V") doItAndPrint("cramerV(.Table)") 
 if (effet == "w") doItAndPrint("cohenW(.Table)")
if (effet == "tau") doItAndPrint("tauK(.Table)") 
 if (effet == "afc") doItAndPrint("vpAFC(.Table)") 
 

if (local == "chi2") doItAndPrint("chisq.test(.Table)$stdres") 
if (local == "pem") doItAndPrint("pem(.Table)") 
 if (local == "Q") doItAndPrint("localYule(.Table)") 

               
        logger("remove(.Table)")
        remove(.Table, envir=.GlobalEnv)

                tkfocus(CommanderWindow())
        }
        OKCancelHelp(helpSubject="xtabs", reset="twoWayTable0")
        radioButtons(name="tableau",
                        buttons=c("tableau1", "tableau2", "tableau3"),
                        values=c("ligne", "colonne", "total"), initialValue=dialog.values$initial.tableau,
                        labels=gettextRcmdr(c("% en ligne", "% en colonne", "% du total")), title=gettextRcmdr("Affichage du tableau de..."))
       
 radioButtons(name="graf",
                        buttons=c("graf1", "graf2", "graf3", "graf4"),
                        values=c("mosaic", "bar", "assoc", "factomap"), initialValue=dialog.values$initial.graf,
                        labels=gettextRcmdr(c("en mosaique", "en barres", "d'association", "sur plan factoriel")), title=gettextRcmdr("Graphique"))
       
 radioButtons(name="test",
                        buttons=c("test1", "test2"),
                        values=c("x2", "fisher"), initialValue=dialog.values$initial.test,
                        labels=gettextRcmdr(c(paste("Chi-2 d'ind","\U00E9","pendance de Pearson",sep=""), "exact de Fisher")), title=gettextRcmdr(paste("Test d'hypoth","\U00E8","ses",sep="")))
       
 radioButtons(name="effet",
                        buttons=c("effet1", "effet2","effet3","effect4"),
                        values=c("V", "w","tau","afc"), initialValue=dialog.values$initial.effet,
                        labels=gettextRcmdr(c("V (Cramer)", "w (Cohen)","Tau (Goodman-Kruskal)","valeur propre 1 (AFC)")), title=gettextRcmdr("Taille d'effet"))

        radioButtons(name="local",
                        buttons=c("local1", "local2", "local3"),
                        values=c("chi2", "pem", "Q"), initialValue=dialog.values$initial.local,
                        labels=gettextRcmdr(c(paste("R", "\U00E9", "sidus standardis","\U00E9","s", 
            sep = ""), "pem (Cibois)", "Q (Yule)")), title=gettextRcmdr("Analyse par case"))

        
       tkgrid(getFrame(rowBox), labelRcmdr(variablesFrame, text="    "), getFrame(columnBox),
sticky="nw")
       
 tkgrid(variablesFrame, sticky="w")
        tkgrid(tableauFrame, sticky="w")
        tkgrid(grafFrame, sticky="w")
        tkgrid(testFrame, sticky="w")
        tkgrid(effetFrame, sticky="w")
        tkgrid(localFrame, sticky="w")

        tkgrid(buttonsFrame,sticky="w")
        dialogSuffix(rows=6, columns=1)
}



