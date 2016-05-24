Factoprefpls <-
function(){
    top<-tktoplevel(borderwidth=5)
    tkwm.title(top,.Facto_gettext("Scatter plot with additional variables"))
    donnee<-get(getRcmdr(".activeDataSet"))
    nomdonnee<-activeDataSet()

######## FRAMES -----------------------------------------------------------------------------------------------------------------------


descFrame <- tkframe(top)
    listdesc<-tklistbox(descFrame,selectmode="extended",yscrollcommand=function(...) tkset(scrdesc,...))
    scrdesc <- tkscrollbar(descFrame,repeatinterval=5,command=function(...)tkyview(listdesc,...))
    tkgrid(listdesc, scrdesc,sticky = "nw")
    tkgrid.configure(scrdesc, sticky = "wns")
    tkgrid.configure(listdesc,sticky = "ew")

   vars<-colnames(donnee)
   vars.desc = NULL
   for (i in (1:ncol(donnee))){
     if (is.numeric(donnee[,i])){
       tkinsert(listdesc,"end",vars[i])
       vars.desc = c(vars.desc,vars[i])
     }
   }

interestFrame <- tkframe(top)
    col.lab<-tklabel(interestFrame,text=.Facto_gettext("'x' variable:"),fg="darkgreen")
    col.val<-tclVar("")
    colj <- tkentry(interestFrame,width=15,textvariable=col.val)
    col2.lab<-tklabel(interestFrame,text=.Facto_gettext("'y' variable:"),fg="darkgreen")
    col2.val<-tclVar("")
    col2j <- tkentry(interestFrame,width=15,textvariable=col2.val)
    
    tkgrid(col.lab, colj,sticky = "nw")
    tkgrid(col2.lab, col2j,sticky = "nw")

optionsFrame <- tkframe(top,relief="ridge",borderwidth=2)
    indiv.check <- tkcheckbutton(top)
    indiv.bool <- tclVar("1")
    tkconfigure(indiv.check,variable=indiv.bool)
    indiv.lab<-tklabel(optionsFrame,text=.Facto_gettext("Plot the graph of the individuals"))
    ortho.check <- tkcheckbutton(top)
    ortho.bool <- tclVar("0")
    tkconfigure(ortho.check,variable=ortho.bool)
    ortho.lab<-tklabel(optionsFrame,text=.Facto_gettext("Same scale for 'x' and 'y'"))
    variab.check <- tkcheckbutton(top)
    variab.bool <- tclVar("1")
    tkconfigure(variab.check,variable=variab.bool)
    variab.lab<-tklabel(optionsFrame,text=.Facto_gettext("Plot the graph of the additional variables"))

    tkgrid(indiv.lab,indiv.check,sticky="w")
    tkgrid(ortho.lab,ortho.check,sticky="w")
    tkgrid(tklabel(optionsFrame,text=""))
    tkgrid(variab.lab,variab.check,sticky="w")

doubleclick<-function(){
    col.nam<-tclvalue(col.val)
    col2.nam<-tclvalue(col2.val)
    if(col.nam==""){
      col.val <- vars.desc[as.numeric(tkcurselection(listdesc))+1]
      tkinsert(colj,"end",col.val)}
    else{
      if(col2.nam==""){
        col2.val <- vars.desc[as.numeric(tkcurselection(listdesc))+1]
        tkinsert(col2j,"end",col2.val)
      }
    }
}
tkbind(listdesc,"<Double-ButtonPress-1>",doubleclick)

###### Fonction principale qui lance prefpls sans fermer la fenêtre------------------------------------------------------------------------------------------------------------
App<-function(){
    variable<-vars.desc[as.numeric(tkcurselection(listdesc))+1]
    if (length(variable)<2) variable = vars.desc
    col.nam<-tclvalue(col.val)
    col2.nam<-tclvalue(col2.val)
    if (col.nam=="") tkmessageBox(message=.Facto_gettext("No variable selected for the 'x' axis"),icon="warning",type="ok")
    if (col2.nam=="") tkmessageBox(message=.Facto_gettext("No variable selected for the 'y' axis"),icon="warning",type="ok")
    done = 0
    if ((col.nam!="")&(col2.nam!="")){
      done = 1
      if (!(col2.nam%in%variable)) variable <- c(col2.nam,variable)
      if (!(col.nam%in%variable)) variable <- c(col.nam,variable)
      for (i in 1:length(variable)){
        if (variable[i]==col.nam) col.pos<-i
        if (variable[i]==col2.nam) col2.pos<-i
      }
      commande.data<-paste(nomdonnee,'.aux', '<-', nomdonnee,'[, c("', paste(variable, collapse='", "'), '")]',sep='')
      justDoIt(commande.data)
      logger(commande.data)
      if (tclvalue(indiv.bool)==1){
        if (tclvalue(ortho.bool)!=1) doItAndPrint(paste('prefpls(',nomdonnee,'.aux, var1=',col.pos,', var2=',col2.pos,', choix = "ind", asp=NA)',sep=''))
        else doItAndPrint(paste('prefpls(',nomdonnee,'.aux, var1=',col.pos,', var2=',col2.pos,', choix = "ind", asp=1)',sep=''))
      }
      if (tclvalue(variab.bool)==1) doItAndPrint(paste('prefpls(',nomdonnee,'.aux, var1=',col.pos,', var2=',col2.pos,', choix = "var")',sep=''))
      justDoIt(paste('remove(',nomdonnee,'.aux)',sep=""))
      logger(paste('remove(',nomdonnee,'.aux)',sep=""))
    }
    return(done)
}


###### Fonction principale qui lance prefpls et ferme la fenêtre------------------------------------------------------------------------------------------------------------
onOK <- function(){
  done = App()
  if (done >0) tkdestroy(top)
}


##### Positionnement des widgets et frames sur la fenêtre 'top' ------------------------------------------------------------------------------------------
App.but <- tkbutton(top,borderwidth=3,width=12,text=.Facto_gettext("Apply"),command=App,fg="blue")

OKCancelHelp(helpSubject="prefpls")
tkgrid(tklabel(top, text=.Facto_gettext("Variables (double-click to select the two variables)"), fg = "blue"), columnspan = 2, sticky = "w")
tkgrid(descFrame,interestFrame,sticky="w")
tkgrid(tklabel(top,text=""))
tkbind(listdesc,"<Double-ButtonPress-1>",doubleclick)
tkgrid(optionsFrame,sticky="w")
# tkgrid(App.but,sticky="w")
# tkgrid(tklabel(top,text=""))
# tkgrid(buttonsFrame, columnspan=2)
  tkgrid(buttonsFrame, App.but)
  tkgrid.configure(buttonsFrame, column=1,sticky="e")
  tkgrid.configure(App.but, column=2,sticky="w")
tkfocus(top)
}
