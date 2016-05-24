nalfcheckcalc <- function(){
initializeDialog(title = gettextRcmdr("Missing Values"))
a <- eval(parse(text =paste('lfnacheck(',ActiveDataSet(),')')))
tt <- tkframe(top)
histFrame <- tkframe(top)


onOK <- function(){
  closeDialog()
  command <- paste('lfnacheck(',ActiveDataSet(),')',sep ="")
  doItAndPrint(command)
  tkfocus(CommanderWindow())
}
Pressedlenhist <- function()
{   barplot(a$duration,
            xlab = "Missing value duration in days",
            ylab = "Number of events",
            main="Missing Value duration")
}

Pressedyearhist <- function()
{   barplot(a$hydrologicalyear[,2],na = a$hydrologicalyear[,1],xlab = "Hydrological year", ylab = "Number of missing values", main = "Missing values per year")
}
lenhist <- buttonRcmdr(top,text=gettextRcmdr("Runlengths of missing values"),command=Pressedlenhist)
yearhist <- buttonRcmdr(top,text=gettextRcmdr("Missing values per year"),command=Pressedyearhist)

OKCancelHelp(helpSubject = "lfnacheck")
               
tkgrid(labelRcmdr(tt,text=gettextRcmdr(paste('The Dataset',
                       ActiveDataSet(),
                       'contains',
                       a$total,
                       'missing values.\nThis are', 100*round(a$percentage,4),
                       'percent of the values.'))))
tkgrid(tt,sticky = "w")
tkgrid(labelRcmdr(histFrame, text = gettextRcmdr("Barchart of:")),lenhist,yearhist, sticky = "w")
tkgrid(histFrame,sticky = "w")
tkgrid(buttonsFrame, sticky="w")
dialogSuffix(rows=7, columns=2)

               
}
