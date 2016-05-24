FactoHCPC <-
function(){
    top<-tktoplevel(borderwidth=5)
    tkwm.title(top,"HCPC")

# définition des polices
  font2<-tkfont.create(family="times",size=12,weight="bold")
  fontheading<-tkfont.create(family="times",size=13,weight="bold")

recupres<-function(envir = .GlobalEnv, ...){
    classres <- ls(envir = envir, all.names = TRUE)
    if (length(classres) == 0) 
        return(classres)
    names(which(sapply(classres, function(.x) inherits(get(.x, 
        envir = envir),c("PCA","CA","MCA","MFA","HMFA")))))
}
      
tousres<-c(listDataSets(),recupres())
  
if (length(which(sapply(ls(envir = .GlobalEnv, all.names = TRUE),function(.x) inherits(get(.x,envir= .GlobalEnv),"CA"))))>0)
  resca<-TRUE 
else resca<-FALSE


objintFrame <- tkframe(top)
    listobjint<-tklistbox(objintFrame,selectmode="single",exportselection=FALSE,yscrollcommand=function(...) tkset(scr,...))
    scr <- tkscrollbar(objintFrame,repeatinterval=5,command=function(...)tkyview(listobjint,...))
    for (i in (1:length(tousres))){
        tkinsert(listobjint,"end",tousres[i])
    }
    tkselection.set(listobjint,0)
    tkgrid(tklabel(objintFrame, text = .Facto_gettext("Select an object"), fg = "blue"), columnspan = 1, sticky = "w")    
    tkgrid(listobjint, scr,sticky = "nw")
    tkgrid.configure(scr, sticky = "wns")
    tkgrid.configure(listobjint,sticky = "ew")
    tkselection.set(listobjint,0)  


HcpcFrame<-tkframe(top,relief="ridge",borderwidth=2)   
 
      tkgrid(tklabel(HcpcFrame, text = .Facto_gettext("Main options"), fg = "blue"), column=1, sticky = "w")      

      if(resca==TRUE){
        lignes <- tkradiobutton (HcpcFrame)
        lignes.lab <- tklabel(HcpcFrame,text=.Facto_gettext("rows"))
        colonnes <- tkradiobutton (HcpcFrame)
        colonnes.lab <- tklabel(HcpcFrame,text=.Facto_gettext("columns"))
      }
      else{
        lignes <- tkradiobutton (HcpcFrame,state="disabled")
        lignes.lab <- tklabel(HcpcFrame,text=.Facto_gettext("rows"))
        colonnes <- tkradiobutton (HcpcFrame,state="disabled")
        colonnes.lab <- tklabel(HcpcFrame,text=.Facto_gettext("columns"))
      }      
      classifCAValue <- tclVar("rows")
      classifCA.lab <- tklabel(HcpcFrame,text=.Facto_gettext("For the result of a CA, perform clustering on:"))
      tkconfigure(lignes,variable=classifCAValue,value="rows")
      tkconfigure(colonnes,variable=classifCAValue,value="columns")      
      
      meth1 <- tkradiobutton (HcpcFrame)
      meth1.lab <- tklabel(HcpcFrame,text=.Facto_gettext("interactive"))
      meth2 <- tkradiobutton (HcpcFrame)
      meth2.lab <- tklabel(HcpcFrame,text=.Facto_gettext("automatic"))
      methValue <- tclVar("0")
      meth.lab <- tklabel(HcpcFrame,text=.Facto_gettext("Choice of the number of clusters:"))
      tkconfigure(meth1,variable=methValue,value="0")
      tkconfigure(meth2,variable=methValue,value="-1")
      
      minmaxhcpc.label<-tklabel(HcpcFrame,text=.Facto_gettext("The optimal number of clusters is chosen between:"))
      minhcpc<-tclVar("3")
      maxhcpc<-tclVar("10")
      minhcpc.entry <-tkentry(HcpcFrame,width="3",textvariable=minhcpc)
      maxhcpc.entry <-tkentry(HcpcFrame,width="3",textvariable=maxhcpc)

      consolid.lab <- tklabel(HcpcFrame,text=.Facto_gettext("Consolidate clusters"))
      consolid.check <- tkcheckbutton(HcpcFrame)
      consolidValue<-tclVar("0")    
      tkconfigure(consolid.check,variable=consolidValue)
      
      graphhcpc.lab <- tklabel(HcpcFrame,text=.Facto_gettext("Print graphs"))
      graphhcpc.check <- tkcheckbutton(HcpcFrame)
      graphhcpcValue <- tclVar("1")
      tkconfigure(graphhcpc.check,variable=graphhcpcValue) 
      
      resu.val<-tclVar("results")
      resu<-tkentry(HcpcFrame,width=6,textvariable=resu.val)
      resu.lab<-tklabel(HcpcFrame,text=.Facto_gettext("Name of the result object:"))                 

      tkgrid(classifCA.lab,lignes.lab,lignes)
      tkgrid(colonnes.lab,colonnes) 
      tkgrid(tklabel(HcpcFrame,text=""))           
      tkgrid(meth.lab,meth1.lab,meth1)
      tkgrid(meth2.lab,meth2)
      tkgrid(tklabel(HcpcFrame,text=""))
      tkgrid(minmaxhcpc.label,minhcpc.entry , maxhcpc.entry)
      tkgrid(tklabel(HcpcFrame,text=""))      
      tkgrid(consolid.lab,consolid.check)
      tkgrid(graphhcpc.lab,graphhcpc.check) 
      tkgrid(resu.lab,resu)
      tkgrid(tklabel(HcpcFrame,text=""))     
      
      tkgrid.configure(minmaxhcpc.label,meth.lab,classifCA.lab,consolid.lab,graphhcpc.lab,resu.lab,column=1,columnspan=4,sticky="w")
      tkgrid.configure(meth1,meth2,lignes,colonnes,consolid.check,graphhcpc.check,column=8,sticky="w")
      tkgrid.configure(meth1.lab,lignes.lab,column=7,columnspan=1,sticky="w")
      tkgrid.configure(meth2.lab,colonnes.lab,column=7,columnspan=1,sticky="w") 
      tkgrid.configure(minhcpc.entry,column=7,columnspan=1,sticky="e")
      tkgrid.configure(maxhcpc.entry,column=8,columnspan=1,sticky="w")
      tkgrid.configure(resu,column=8,columnspan=1,sticky="w")
      
      tkgrid.columnconfigure(HcpcFrame,0, minsize=5)
      tkgrid.columnconfigure(HcpcFrame,5, minsize=5)
      tkgrid.columnconfigure(HcpcFrame,8, minsize=5) 

#Fonction outputs
    env<-environment()
    Gclust <- TRUE
    Gdescvar <- TRUE
    Gdescaxes <- TRUE
    Gdescind <- TRUE

    OnHCPCSorties<-function()
    {
      SortiesWin<- tktoplevel()
      tkwm.title(SortiesWin, .Facto_gettext("Outputs options"))

      onOKSorties<-function()
      {
        if(tclvalue(clust.bool)=="1") assign("Gclust", TRUE, envir=env)
        else assign("Gclust", FALSE, envir=env)
        
        if(tclvalue(descvar.bool)=="1") assign("Gdescvar", TRUE, envir=env)
        else assign("Gdescvar", FALSE, envir=env)
        
        if(tclvalue(descaxes.bool)=="1") assign("Gdescaxes", TRUE, envir=env)
        else assign("Gdescaxes", FALSE, envir=env)
        
        if(tclvalue(descind.bool)=="1") assign("Gdescind", TRUE, envir=env)
        else assign("Gdescind", FALSE, envir=env)        

        tkdestroy(SortiesWin)
      }
      OKSorties.but<-tkbutton(SortiesWin, text="OK", width=8,command=onOKSorties)

      tkgrid(tklabel(SortiesWin, text=""))
      tkgrid(tklabel(SortiesWin, text = .Facto_gettext("Select outputs options"), fg = "red",font=font2), column=1, columnspan = 8, sticky = "w")

      clust.lab <- tklabel(SortiesWin,text=.Facto_gettext("Print clusters individuals belong to "))
      clust.check <- tkcheckbutton(SortiesWin)
      if(Gclust) clust.bool<-tclVar("1")
      else clust.bool<-tclVar("0")      
      tkconfigure(clust.check,variable=clust.bool)
      
      descvar.lab <- tklabel(SortiesWin,text=.Facto_gettext("Print description of clusters by continuous variables "))
      descvar.check <- tkcheckbutton(SortiesWin)
      if(Gdescvar) descvar.bool<-tclVar("1")
      else descvar.bool<-tclVar("0")      
      tkconfigure(descvar.check,variable=descvar.bool)
      
      descaxes.lab <- tklabel(SortiesWin,text=.Facto_gettext("Print description of clusters by factors"))
      descaxes.check <- tkcheckbutton(SortiesWin)
      if(Gdescaxes) descaxes.bool<-tclVar("1")
      else descaxes.bool<-tclVar("0")      
      tkconfigure(descaxes.check,variable=descaxes.bool)
      
      descind.lab <- tklabel(SortiesWin,text=.Facto_gettext("Print parangons and most typical individuals for each cluster"))
      descind.check <- tkcheckbutton(SortiesWin)
      if(Gdescind) descind.bool<-tclVar("1")
      else descind.bool<-tclVar("0")      
      tkconfigure(descind.check,variable=descind.bool)                  

      tkgrid(tklabel(SortiesWin,text=""))
      tkgrid(clust.lab,clust.check)
      tkgrid(descvar.lab,descvar.check)
      tkgrid(descaxes.lab,descaxes.check)
      tkgrid(descind.lab,descind.check)
      tkgrid(tklabel(SortiesWin,text=""))     
      tkgrid(OKSorties.but)
      tkgrid(tklabel(SortiesWin,text=""))
      tkfocus(SortiesWin)
      
      tkgrid.configure(clust.lab,descvar.lab,descaxes.lab,descind.lab,column=1,sticky="w")
      tkgrid.configure(clust.check,descvar.check,descaxes.check,descind.check,column=8,sticky="w")
      tkgrid.configure(OKSorties.but,column=4,columnspan=2)
    }
    
    #SortiesFrame<-tkframe(HcpcFrame)
    HcpcSorties.but<-tkbutton(HcpcFrame, text=.Facto_gettext("Outputs"), command=OnHCPCSorties, width=13,fg="darkred")
    tkgrid(HcpcSorties.but)
    tkgrid.configure(HcpcSorties.but, column=8, columnspan=1,sticky="w")
    

#Fonction qui lance l'HCPC
 App<-function(){
    objint<-c(tclvalue(tkcurselection(listobjint)))
    obj<-tousres[as.numeric(tkcurselection(listobjint))+1]
    resultat<-tclvalue(resu.val)
    done = 1 
    
    commande.hcpc<-paste(resultat, '<-HCPC(', obj, ' ,nb.clust=', tclvalue(methValue), ',consol=', tclvalue(consolidValue),',min=', as.numeric(tclvalue(minhcpc)),',max=',as.numeric(tclvalue(maxhcpc)),',cluster.CA="',tclvalue(classifCAValue),'",graph=', tclvalue(graphhcpcValue), ')', sep="")
    justDoIt(commande.hcpc)
    logger(commande.hcpc) 
    
    if(Gclust) doItAndPrint(paste(resultat,'$data.clust[,ncol(',resultat,'$data.clust),drop=F]', sep=""))     
    if(Gdescvar) doItAndPrint(paste(resultat,'$desc.var', sep=""))
    if(Gdescaxes) doItAndPrint(paste(resultat,'$desc.axes', sep=""))
    if(Gdescind) doItAndPrint(paste(resultat,'$desc.ind', sep=""))

    return(done)
 }
 
# Fonction principale qui lance HCPC et ferme la fenêtre------------------------------------------------------------------------------------------------------------
onOK <- function(){
  done = App()
  if (done >0) tkdestroy(top)
}

 
#Mise en place des frames
App.but <- tkbutton(top,borderwidth=3,width=12,text=.Facto_gettext("Apply"),command=App,fg="blue")
OKCancelHelp(helpSubject="HCPC")
titre<-tklabel(top, text=.Facto_gettext("Hierarchical Clustering on Principal Components (HCPC)"),font=fontheading)
tkgrid.configure(titre,column=1,columnspan=3)
tkgrid(tklabel(top,text=""))
tkgrid(objintFrame,HcpcFrame,sticky="w")
tkgrid.configure(objintFrame,column=1,columnspan=1,sticky="w")
tkgrid.configure(HcpcFrame,column=2,columnspan=2,sticky="w")
tkgrid(tklabel(top,text="")) 
# tkgrid(App.but)
# tkgrid.configure(App.but,column=1,columnspan=1,sticky="w")
# tkgrid(tklabel(top,text=""))
# tkgrid(buttonsFrame)
# tkgrid.configure(buttonsFrame, column=1,columnspan=2)
  tkgrid(buttonsFrame, App.but)
  tkgrid.configure(buttonsFrame, column=1,sticky="e")
  tkgrid.configure(App.but, column=2,sticky="w")
tkfocus(top)

tkgrid.columnconfigure(top,0, minsize=10)
tkgrid.columnconfigure(top,8, minsize=10) 
      
}
