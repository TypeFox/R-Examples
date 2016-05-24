"concatenate"<-function(){
nwd<-0;nwd<<-nwd
if(length(txt3)==1 & length(txt4)==1) nwd<-levels(as.factor(c(txt1,txt2)))
if(length(txt3)!=1 & length(txt4)==1) nwd<-levels(as.factor(c(txt1,txt2,txt3)))
if(length(txt4)!=1)nwd<-txt4
nwd<<-nwd

tt1<-tktoplevel()
tkwm.title(tt1,"ManagementOligonucleotide database")
  tkgrid(tklabel(tt1,text="                                                                                         "))
  tkgrid(tklabel(tt1,text="Choose on which oligonucleotide database you want to replace with concatenate degenerate primer"))
  index<-c("Export","OligoNucleotide data base 1","OligoNucleotide data base 2","OligoNucleotide data base 3")
  index1<-tkwidget(tt1,"ComboBox",editable=FALSE,values=index,width=40,height=4)
  tkgrid(index1)
calc<-function(){
    dbdb <- unlist(as.numeric(tcl(index1, "getvalue")) + 1)
    newdb<-0;newdb<<-newdb
    if (dbdb == 1) newdb <<- nwd
    if (dbdb == 1) exportglobal(n=4)
    if (dbdb == 2) txt1 <<- nwd
    if (dbdb == 3) txt2 <<- nwd
    if (dbdb == 4) txt3 <<- nwd
    if (dbdb == 2) ONDB1name <<- "Concatenated degenerate primer"
    if (dbdb == 3) ONDB2name <<- "Concatenated degenerate primer"
    if (dbdb == 4) ONDB3name <<- "Concatenated degenerate primer"
    tkdestroy(tt1)}
tkgrid(tkbutton(tt1,text="Replace",command=calc))
tkgrid(tklabel(tt1,text="   "))
tkfocus(tt1)
}
