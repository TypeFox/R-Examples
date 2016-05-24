FactoFAMD <-
function()
{

# Fonction pour la gestion des noms ############################################

  nom.correct<-function(text, liste=NULL)
  {
    text<-chartr("^\ ", "...", text)
    if(!is.null(liste)) {
      while(text %in% liste) text<-paste(text, ".bis", sep="")
    }
    return(text)  
  }
################################################################################

#    Création des fonctions pour les options via nouvelle fenêtre graphique   

  top<-tktoplevel(borderwidth=10)
  tkwm.title(top,.Facto_gettext("FAMD"))
  tkwm.geometry(top, "-50+50")
  
  # définition des polices
  font2<-tkfont.create(family="times",size=12,weight="bold")
  fontheading<-tkfont.create(family="times",size=11,weight="bold")

  # récupération du jeu de données actif
  donnee<-get(getRcmdr(".activeDataSet"))
  vars<-colnames(donnee)
  rows<-rownames(donnee)

    listFrame <- tkframe(top,borderwidth=2)                                                                                                          
    lab1 = tklabel(listFrame,text=.Facto_gettext("Select quantitative variables"),fg="blue")
    lab2 = tklabel(listFrame,text=.Facto_gettext("Select factors"),fg="blue")
    lab3 = tklabel(listFrame,text="      ")
    tkgrid(lab1,lab3,lab2)
    tkgrid.configure(lab1,column=1, columnspan=2, sticky = "nw")
    tkgrid.configure(lab2,column=4, columnspan=2, sticky = "ne")
    tkgrid.configure(lab3,column=3, columnspan=1, sticky = "n")

    listdesc<-tklistbox(listFrame,selectmode="extended",exportselection=FALSE,yscrollcommand=function(...) tkset(scr,...))
    scr <- tkscrollbar(listFrame,repeatinterval=5,command=function(...)tkyview(listdesc,...))
    tkselection.set(listdesc,0)
    listfact<-tklistbox(listFrame,selectmode="extended",exportselection=FALSE,yscrollcommand=function(...) tkset(scrfact,...))
    scrfact <- tkscrollbar(listFrame,repeatinterval=5,command=function(...)tkyview(listfact,...))
    vars<-colnames(donnee)
    vars.fact = NULL
    vars.desc = NULL
    for (i in (1:ncol(donnee))){
      if (is.numeric(donnee[,i])){
        tkinsert(listdesc,"end",vars[i])
        vars.desc = c(vars.desc,vars[i])
      }
      else {
        vars.fact = c(vars.fact,vars[i])
        tkinsert(listfact,"end",vars[i])
      }
    }

    tkgrid(listdesc, scr,tklabel(listFrame,text="                                  "),listfact,scrfact,sticky = "nw")
    tkgrid.configure(scr, sticky = "wns", column=2, columnspan=1)
    tkgrid.configure(listdesc,sticky = "ew", column=1, columnspan=1)
    tkgrid.configure(scrfact,sticky = "wns", column=5, columnspan=1)
    tkgrid.configure(listfact,sticky = "ew", column=4, columnspan=1)

 
  Fillu.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    variablefact<-NULL
    .FilluLabel<-tclVar(paste(firstLabel, "", sep=" "))
    .factors<-Factors()
    OnFillu<-function()
    {
      if(length(.factors)==0) errorCondition(recall=NULL, message=.Facto_gettext("No Factor available"))
          
      FilluWin<-tktoplevel()
      tkwm.title(FilluWin,.Facto_gettext("Choice of supplementary factors"))
      #création de la fonction FOK.funct
      FOK.funct<-function()
      {
        fact.select<-listfact.nom[as.numeric(tkcurselection(listfact))+1]
        if(length(fact.select)==0) {
          assign("variablefact", NULL, envir=env)
          tclvalue(.FilluLabel)<-paste(firstLabel, "", sep=" ")
          tkconfigure(Fillu.but, fg="black")
          tkdestroy(FilluWin)
          return()
        }
        assign("variablefact", fact.select, envir=env)
        tclvalue(.FilluLabel)<-paste(label, "", sep=" ")
        tkconfigure(Fillu.but, fg="blue")
        tkdestroy(FilluWin)
      }
      
      # création et mise en page de la fenetre Fillu
      listfact<-tklistbox(FilluWin,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrfact,...)) # Liste vide
      scrfact <-tkscrollbar(FilluWin,repeatinterval=5,command=function(...)tkyview(listfact,...))
      listfact.nom<-NULL
      indice<-0
      for (i in (1:ncol(donnee))) {
        if (is.factor(donnee[,i])) {
          tkinsert(listfact,"end",vars[i]) # On renseigne la liste
          listfact.nom<-c(listfact.nom,vars[i])
          if(vars[i] %in% variablefact) tkselection.set(listfact, indice)
          indice<-indice+1
        }
      }
  
      FOK.but<-tkbutton(FilluWin, text="OK", width=16,command=FOK.funct)

      tkgrid(tklabel(FilluWin, text=""))
      tkgrid(tklabel(FilluWin, text = .Facto_gettext("Select supplementary factor(s)"), fg = "blue"), column=1, columnspan = 1, sticky = "ew")
      tkgrid(listfact, scrfact, sticky = "nw")
      tkgrid.configure(scrfact, sticky = "ens", columnspan=1)
      tkgrid.configure(listfact, sticky = "ew", column=1, columnspan=1)
      tkgrid(tklabel(FilluWin, text=""))
      tkgrid(FOK.but, column=1,columnspan=1, sticky="ew")
      tkgrid(tklabel(FilluWin, text=""))
      tkgrid.columnconfigure(FilluWin,0, minsize=25)
      tkgrid.columnconfigure(FilluWin,2, minsize=25)
  }  

   FilluFrame<-tkframe(IlluFrame)
   if(length(.factors)==0){
     Fillu.but<-tkbutton(FilluFrame, text=.Facto_gettext("No factors available"), borderwidth=3)
     tkconfigure(Fillu.but, fg="grey")
   }
   else Fillu.but<-tkbutton(FilluFrame, textvariable=.FilluLabel, command=OnFillu, borderwidth=3)
   tkgrid(Fillu.but, sticky="ew")

##   Fillu.but<-tkbutton(FilluFrame, textvariable=.FilluLabel, command=OnFillu, borderwidth=3)
##   tkgrid(Fillu.but, sticky="ew")
  })

  #! fonction pour le choix des variables quantitatives supplémentaires 
  Dillu.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    variableillu<-NULL
    .DilluLabel<-tclVar(paste(firstLabel, "", sep=" "))
    OnDillu<-function()
    { 
      DilluWin<-tktoplevel()
      tkwm.title(DilluWin,.Facto_gettext("Select supplementary quantitative variables"))
      #création de la fonction DOK.funct
      DOK.funct<-function()
      {
        vsup.select<-listvar.nom[as.numeric(tkcurselection(listvar))+1]
        if(length(vsup.select)==0)
        {
          assign("variableillu", NULL, envir=env)
          tclvalue(.DilluLabel)<-paste(firstLabel, "", sep=" ")
          tkconfigure(Dillu.but, fg="black")
          tkdestroy(DilluWin)
          return()
        }
        assign("variableillu", vsup.select, envir=env)
        tclvalue(.DilluLabel)<-paste(label, "", sep=" ")
        tkconfigure(Dillu.but, fg="blue")
        tkdestroy(DilluWin)
      }
      
      # création et mise en page de la fenetre Dillu
      listvar<-tklistbox(DilluWin,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrvar,...)) # Liste vide
      scrvar <-tkscrollbar(DilluWin,repeatinterval=5,command=function(...)tkyview(listvar,...)) 
      listvar.nom<-NULL
      indice<-0
      for (i in (1:ncol(donnee))) {
          if (is.numeric(donnee[,i])) {
            tkinsert(listvar,"end",vars[i]) # On renseigne la liste
            listvar.nom<-c(listvar.nom,vars[i])
            if(vars[i] %in% variableillu) tkselection.set(listvar, indice)
            indice<-indice+1
          }
      }
  
      DOK.but<-tkbutton(DilluWin, text="OK", width=16,command=DOK.funct)

      tkgrid(tklabel(DilluWin, text=""))
        tkgrid(tklabel(DilluWin, text = .Facto_gettext("Select supplementary quantitative variables"), fg = "blue"), column=1, columnspan = 1, sticky = "ew")
        tkgrid(listvar, scrvar, sticky = "nw")
        tkgrid.configure(scrvar, sticky = "ens", columnspan=1)
        tkgrid.configure(listvar, sticky = "ew", column=1, columnspan=1)
        tkgrid(tklabel(DilluWin, text=""))
        tkgrid(DOK.but, column=1,columnspan=1, sticky="ew")
        tkgrid(tklabel(DilluWin, text=""))
        tkgrid.columnconfigure(DilluWin,0, minsize=25)
      tkgrid.columnconfigure(DilluWin,2, minsize=25)
  }  

   DilluFrame<-tkframe(IlluFrame)
   if(length(listNumeric())==0){
     Dillu.but<-tkbutton(DilluFrame, text=.Facto_gettext("No quantitative variable available"), borderwidth=3)
     tkconfigure(Dillu.but, fg="grey")
   }
   else Dillu.but<-tkbutton(DilluFrame, textvariable=.DilluLabel, command=OnDillu, borderwidth=3)
   tkgrid(Dillu.but, sticky="ew")
  })
  
  #! fonction pour le choix des individus supplémentaires 
  Iillu.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    individuillu<-NULL
    .IilluLabel<-tclVar(paste(firstLabel, "", sep=" "))
    OnIillu<-function()
    {   
      IilluWin<-tktoplevel()
      tkwm.title(IilluWin,.Facto_gettext("Select supplementary individual(s)"))
      #création de la fonction IOK.funct
      IOK.funct<-function()
      {
        ind.select<-rows[as.numeric(tkcurselection(listind))+1]
        if(length(ind.select)==0) {
          assign("individuillu", NULL, envir=env)
          tclvalue(.IilluLabel)<-paste(firstLabel, "", sep=" ")
          tkconfigure(Iillu.but, fg="black")
          tkdestroy(IilluWin)
          return()
        }
        assign("individuillu", ind.select, envir=env)
        tclvalue(.IilluLabel)<-paste(label, "", sep=" ")
        tkconfigure(Iillu.but, fg="blue")
        tkdestroy(IilluWin)
      }
      
      # création et mise en page de la fenetre Fillu
      listind<-tklistbox(IilluWin,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrind,...)) # Liste vide
      scrind <-tkscrollbar(IilluWin,repeatinterval=5,command=function(...)tkyview(listind,...)) 
      indice<-0
      for (i in (1:nrow(donnee))) {
        tkinsert(listind,"end",rows[i]) # On renseigne la liste
        if(rows[i] %in% individuillu) tkselection.set(listind, indice)
        indice<-indice+1
      }
  
      IOK.but<-tkbutton(IilluWin, text="OK", width=16,command=IOK.funct)

      tkgrid(tklabel(IilluWin, text=""))
      tkgrid(tklabel(IilluWin, text = .Facto_gettext("Select supplementary individual(s)"), fg = "blue"), column=1, columnspan = 1, sticky = "ew")
      tkgrid(listind, scrind, sticky = "nw")
      tkgrid.configure(scrind, sticky = "ens", columnspan=1)
      tkgrid.configure(listind, sticky = "ew", column=1, columnspan=1)
      tkgrid(tklabel(IilluWin, text=""))
      tkgrid(IOK.but, column=1,columnspan=1, sticky="ew")
      tkgrid(tklabel(IilluWin, text=""))
      tkgrid.columnconfigure(IilluWin,0, minsize=25)
      tkgrid.columnconfigure(IilluWin,2, minsize=25)
  }  

   IilluFrame<-tkframe(IlluFrame)
   Iillu.but<-tkbutton(IilluFrame, textvariable=.IilluLabel, command=OnIillu, borderwidth=3)
   tkgrid(Iillu.but, sticky="ew")
  })
    
  
    #! fonction pour la réinitialisation des paramètres
  Reinitializ.funct<-function()
  {
    tkdestroy(top)
    FactoFAMD()
  }


  #! fonction pour le choix des éléments de sortie
  Sortie.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.sortie<-0
    #déclaration des variables
    Rpropre<-FALSE
    RFichier <- ""
    Rgroupe<-FALSE
    Rindividu<-FALSE
    Rindsup<-FALSE
#    Rquantisummary<-FALSE
    Rquanti<-FALSE
    Rquantisup<-FALSE    
#    Rqualisummary<-FALSE
    Rquali<-FALSE
    Rqualisup<-FALSE
    Rdescdim<-FALSE

    .SortieLabel<-tclVar(paste(firstLabel, "", sep=" "))

    OnSortie<-function()
    {
      SortieWin<-tktoplevel()
      tkwm.title(SortieWin,.Facto_gettext("Output options"))

      #création de la fonction onOKsub
      onOK.sortie<-function()
      {
        assign("compteur.sortie", compteur.sortie+1, envir=env)
        if(compteur.sortie>0) tclvalue(.SortieLabel)<-paste(label, "", sep=" ")
        tkconfigure(Sortie.but, fg="blue")
        
        if(tclvalue(eigValue)=="1") assign("Rpropre", TRUE, envir=env)
        else assign("Rpropre", FALSE, envir=env)
        
        if(tclvalue(groupeValue)=="1") assign("Rgroupe", TRUE, envir=env)
        else assign("Rgroupe", FALSE, envir=env)
        
        if(tclvalue(indValue)=="1") assign("Rindividu", TRUE, envir=env)
        else assign("Rindividu", FALSE, envir=env)
        
        if(tclvalue(ind.sup.Value)=="1") assign("Rindsup", TRUE, envir=env)
        else assign("Rindsup", FALSE, envir=env)
        
#        if(tclvalue(quantiSummaryValue)=="1") assign("Rquantisummary", TRUE, envir=env)
#        else assign("Rquantisummary", FALSE, envir=env)
        
        if(tclvalue(quantiValue)=="1") assign("Rquanti", TRUE, envir=env)
        else assign("Rquanti", FALSE, envir=env)
        
        if(tclvalue(quantisupValue)=="1") assign("Rquantisup", TRUE, envir=env)
        else assign("Rquantisup", FALSE, envir=env)
        
#        if(tclvalue(qualiSummaryValue)=="1") assign("Rqualisummary", TRUE, envir=env)
#        else assign("Rqualisummary", FALSE, envir=env)
        
        if(tclvalue(qualiValue)=="1") assign("Rquali", TRUE, envir=env)
        else assign("Rquali", FALSE, envir=env)
        
        if(tclvalue(qualisupValue)=="1") assign("Rqualisup", TRUE, envir=env)
        else assign("Rqualisup", FALSE, envir=env)
        
        if(tclvalue(descdimValue)=="1") assign("Rdescdim", TRUE, envir=env)
        else assign("Rdescdim", FALSE, envir=env)
        
        if (tclvalue(Fichier)=="") assign("RFichier", NULL, envir=env)
        assign("RFichier", tclvalue(Fichier), envir=env)

        tkdestroy(SortieWin)
      
      }
      
        eig.lab <-tklabel(SortieWin, text=.Facto_gettext("Eigenvalues"))
        eig.check <- tkcheckbutton(SortieWin)
        if(Rpropre) eigValue <- tclVar("1")
        else eigValue <- tclVar("0")
        tkconfigure(eig.check,variable=eigValue)
        
        groupe.lab <-tklabel(SortieWin, text=.Facto_gettext("Results for the variables"))
        groupe.check <- tkcheckbutton(SortieWin)
        if(Rgroupe) groupeValue <- tclVar("1")
        else groupeValue <- tclVar("0")
        tkconfigure(groupe.check,variable=groupeValue)
        
        ind.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for the active individuals"))
        ind.check <- tkcheckbutton(SortieWin)
        if(Rindividu) indValue <- tclVar("1")
        else indValue <- tclVar("0")
        tkconfigure(ind.check,variable=indValue)

        ind.sup.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for the supplementary individuals"))
        ind.sup.check <- tkcheckbutton(SortieWin)
        if(Rindsup) ind.sup.Value <- tclVar("1")
        else ind.sup.Value <- tclVar("0")
        tkconfigure(ind.sup.check,variable=ind.sup.Value)
        
#        quantiSummary.lab<-tklabel(SortieWin,text=.Facto_gettext("Summary of the quantitative variables"))
#        quantiSummary.check <- tkcheckbutton(SortieWin)
#        if(Rquantisummary) quantiSummaryValue <- tclVar("1")
#        else quantiSummaryValue <- tclVar("0")
#        tkconfigure(quantiSummary.check,variable=quantiSummaryValue)
      
        quanti.lab<-tklabel(SortieWin,text=.Facto_gettext("Results of the quantitative variables"))
        quanti.check <- tkcheckbutton(SortieWin)
        if(Rquanti) quantiValue <- tclVar("1")
        else quantiValue <- tclVar("0")
        tkconfigure(quanti.check,variable=quantiValue)

        quantisup.lab<-tklabel(SortieWin,text=.Facto_gettext("Results of the supplementary quantitative variables"))
        quantisup.check <- tkcheckbutton(SortieWin)
        if(Rquantisup) quantisupValue <- tclVar("1")
        else quantisupValue <- tclVar("0")
        tkconfigure(quantisup.check,variable=quantisupValue)
      
#        qualiSummary.lab<-tklabel(SortieWin,text=.Facto_gettext("Summary of the qualitative variables"))
#        qualiSummary.check <- tkcheckbutton(SortieWin)
#        if(Rqualisummary) qualiSummaryValue <- tclVar("1")
#        else qualiSummaryValue <- tclVar("0")
#        tkconfigure(qualiSummary.check,variable=qualiSummaryValue)
      
        quali.lab<-tklabel(SortieWin,text=.Facto_gettext("Results of the qualitative variables"))
        quali.check <- tkcheckbutton(SortieWin)
        if(Rquali) qualiValue <- tclVar("1")
        else qualiValue <- tclVar("0")
        tkconfigure(quali.check,variable=qualiValue)

        qualisup.lab<-tklabel(SortieWin,text=.Facto_gettext("Results of the supplementary qualitative variables"))
        qualisup.check <- tkcheckbutton(SortieWin)
        if(Rqualisup) qualisupValue <- tclVar("1")
        else qualisupValue <- tclVar("0")
        tkconfigure(qualisup.check,variable=qualisupValue)
      
      descdim.lab<-tklabel(SortieWin, text=.Facto_gettext("Description of the dimensions"))
      descdim.check<-tkcheckbutton(SortieWin)
      if(Rdescdim) descdimValue<-tclVar("1")
      else descdimValue<-tclVar("0")
      tkconfigure(descdim.check,variable=descdimValue)

      RFichierFrame<-tkframe(SortieWin,borderwidth=2)
      if (is.null(RFichier)) Fichier <- tclVar("")
      else Fichier<-tclVar(RFichier)
      Fichier.entry <-tkentry(RFichierFrame,width="40",textvariable=Fichier)
      tkgrid(tklabel(RFichierFrame,text=.Facto_gettext("Print results on a 'csv' file")),Fichier.entry)

      SortieOK.but<-tkbutton(SortieWin,text="OK",width=16,command=onOK.sortie)

      tkgrid(tklabel(SortieWin, text = .Facto_gettext("Select output options"), fg ="blue"),  columnspan = 2, sticky = "w")
      tkgrid(tklabel(SortieWin, text = " "))
      tkgrid(eig.lab,eig.check,sticky="w")
      tkgrid(groupe.lab,groupe.check,sticky="w")
      tkgrid(ind.lab,ind.check,sticky="w")
      if (!is.null(individuillu)) tkgrid(ind.sup.lab,ind.sup.check,sticky="w")
#      tkgrid(quantiSummary.lab,quantiSummary.check,sticky="w")
      tkgrid(quanti.lab,quanti.check,sticky="w")      
      if (!is.null(variableillu)) tkgrid(quantisup.lab,quantisup.check,sticky="w")
#      tkgrid(qualiSummary.lab,qualiSummary.check,sticky="w")
      tkgrid(quali.lab,quali.check,sticky="w")        
      if (!is.null(variablefact)) tkgrid(qualisup.lab,qualisup.check,sticky="w")
      tkgrid(descdim.lab,descdim.check,sticky="w")
      tkgrid(tklabel(SortieWin, text = " "))
      tkgrid(RFichierFrame)
      tkgrid(SortieOK.but)
      tkgrid(tklabel(SortieWin, text = " "))
   }
    
    SortieFrame<-tkframe(IlluFrame)
    Sortie.but<-tkbutton(SortieFrame, textvariable=.SortieLabel, command=OnSortie, borderwidth=3)
    tkgrid(Sortie.but, sticky="ew")
  })



  #! fonction pour la gestion des options graphiques 
  PLOT.FAMD<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.graph<-0
    .PlotLabel<-tclVar(paste(firstLabel, "", sep=" "))
    #déclaration des variables
    Gchoix<-TRUE
    GTitle<-NULL
    Gcol.var<-Gcol.var.tmp<-"red"
    Gcol.quanti.sup<-Gcol.quanti.sup.tmp<-"darkred"
    Gcol.quali<-Gcol.quali.tmp<-"green"
    Gcol.quali.sup<-Gcol.quali.sup.tmp<-"darkgreen"
    GAxeGrpe<-c(1,2)
    Glabel<-TRUE
        
    Rchoix<-TRUE
    RTitle<-NULL
    Rlabel.indMoy<-TRUE
    Rlabel.quali<-TRUE    
    Rhabillage<-"none"
    Rinvisible<-NULL
    RXlimInd<-NULL
    RYlimInd<-NULL
    
    Wchoix=TRUE
    WTitle<-NULL
    WAxeVar<-c(1,2)
    Wlabel.var<-TRUE
    Wcol.quanti.sup<-Wcol.quanti.sup.tmp<-"blue"
    Wcol.var<-Wcol.var.tmp<-"black"
    Winvisible<-NULL
    Wlim.cos<-0.
    
        
    OnPlot<-function()
    {
      PlotWin<-tktoplevel()
      tkwm.title(PlotWin,.Facto_gettext("Graphical options"))
      tkwm.geometry(PlotWin, "-100+50")
      PlotWin2<-tkframe(PlotWin)

      #création de la fonction onOKsub
      onOKsub<-function()
      {
        assign("compteur.graph", compteur.graph+1, envir=env)
        if(compteur.graph>0) tclvalue(.PlotLabel)<-paste(label, .Facto_gettext(""), sep=" ")
        tkconfigure(Plot.but, fg="blue")

        # gestion des entrées de la partie graphique des Groupes
        if(tclvalue(grpe.check.value)==1) assign("Gchoix", TRUE, envir=env)
        else assign("Gchoix", FALSE, envir=env)

        if(Gchoix) {
          if (tclvalue(GTitre)==" ") assign("GTitle", NULL, envir=env)
          assign("GTitle", tclvalue(GTitre), envir=env)

          assign("Gcol.var", Gcol.var.tmp, envir=env)
          assign("Gcol.quanti.sup", Gcol.quanti.sup.tmp, envir=env)
          assign("Gcol.quali", Gcol.quali.tmp, envir=env)
          assign("Gcol.quali.sup", Gcol.quali.sup.tmp, envir=env)
          
          label.tmp.grpe<-tclvalue(label.grpe.checkValue)
          if(label.tmp.grpe==1) assign("Glabel", TRUE, envir=env)
          else assign("Glabel", FALSE, envir=env)
        }
                    
        # gestion des entrées de la partie graphique des variables
        if(tclvalue(var.check.value)==1) assign("Wchoix", TRUE, envir=env)
        else assign("Wchoix", FALSE, envir=env)

        if(Wchoix) {
          if (tclvalue(WTitre)==" ") assign("WTitle", NULL, envir=env)
          assign("WTitle", tclvalue(WTitre), envir=env)

          assign("Wlim.cos", tclvalue(WlimCosValue), envir=env)

          label.tmp.var<-tclvalue(label.var.checkValue)
          if(label.tmp.var==1) assign("Wlabel.var", TRUE, envir=env)
          else assign("Wlabel.var", FALSE, envir=env)

          assign("Wcol.var", Wcol.var.tmp, envir=env)
          assign("Wcol.quanti.sup", Wcol.quanti.sup.tmp, envir=env)

          if(tclvalue(inv.Value)=="aucun")  assign("Winvisible", NULL, envir=env)
          else assign("Winvisible", tclvalue(inv.Value), envir=env)
        }

        # gestion des entrées de la partie graphique des individus
        if(tclvalue(ind.check.value)==1) assign("Rchoix", TRUE, envir=env)
        else assign("Rchoix", FALSE, envir=env)

        if(Rchoix) {
          if (tclvalue(Titre)==" ") assign("RTitle", NULL, envir=env)
          assign("RTitle", tclvalue(Titre), envir=env)

          label.tmp.indMoy<-tclvalue(label.indMoy.checkValue)
          label.tmp.quali<-tclvalue(label.quali.checkValue)
          if(label.tmp.indMoy==1) assign("Rlabel.indMoy", TRUE, envir=env)
          else assign("Rlabel.indMoy", FALSE, envir=env)
          if(label.tmp.quali==1) assign("Rlabel.quali", TRUE, envir=env)
          else assign("Rlabel.quali", FALSE, envir=env)
          
          habillage.tmp<-listgraph.nom[as.numeric(tkcurselection(listgraph))+1]
          if(length(habillage.tmp)==0) assign("Rhabillage","none", envir=env)
          else assign("Rhabillage", habillage.tmp, envir=env)

          if(tclvalue(XlimIndMin)=="" | tclvalue(XlimIndMax)=="") assign("RXlimInd", NULL, envir=env)
          else assign("RXlimInd", c(as.numeric(tclvalue(XlimIndMin)), as.numeric(tclvalue(XlimIndMax))), envir=env)
          if(tclvalue(YlimIndMin)=="" | tclvalue(YlimIndMax)=="") assign("RYlimInd", NULL, envir=env)
          else assign("RYlimInd", c(as.numeric(tclvalue(YlimIndMin)), as.numeric(tclvalue(YlimIndMax))), envir=env)
          
          inv.ind.tmp<-tclvalue(inv.ind.checkValue)
          inv.ind.sup.tmp<-tclvalue(inv.ind.sup.checkValue)
          inv.quali.tmp<-tclvalue(inv.quali.checkValue)
          assign("Rinvisible", NULL, envir=env)
          if(inv.ind.tmp=="1") assign("Rinvisible", c(Rinvisible, "ind"), envir=env)
          if(inv.ind.sup.tmp=="1") assign("Rinvisible", c(Rinvisible, "ind.sup"), envir=env)          
          if(inv.quali.tmp=="1") assign("Rinvisible", c(Rinvisible, "quali"), envir=env)          
        }
        tkdestroy(PlotWin)
      }
    
      # création l'interface "options graphiques"
    
      ##########################
      # construction de la partie graphique des Groupes
      PlotGrpeFrame<-tkframe(PlotWin2, borderwidth=5, relief="groove")
  
      GchoixFrame<-tkframe(PlotGrpeFrame,borderwidth=2)
      grpe.check<-tkcheckbutton(GchoixFrame)
      if(Gchoix) grpe.check.value<-tclVar("1")
      else grpe.check.value<-tclVar("0")
      tkconfigure(grpe.check, variable=grpe.check.value)
      tkgrid(tklabel(GchoixFrame, text=.Facto_gettext("Graph of all the variables"), font=font2),grpe.check)
      tkgrid(tklabel(GchoixFrame, text="  "))
  
      GTitleFrame<-tkframe(PlotGrpeFrame,borderwidth=2)
      if (is.null(GTitle)) GTitre <- tclVar(" ")
      else GTitre<-tclVar(GTitle)
      GTitre.entry <-tkentry(GTitleFrame,width="40",textvariable=GTitre)
      tkgrid(tklabel(GTitleFrame,text=.Facto_gettext("Title of the graph")),GTitre.entry)
    
    GcolFrame<-tkframe(PlotGrpeFrame,borderwidth=2)
    Gcol.var.value <- Gcol.var
    canvas.var <- tkcanvas(GcolFrame,width="80",height="25",bg=Gcol.var.value)
    ChangeColor.var <- function()
    {
      Gcol.var.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Gcol.var.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Gcol.var.value)>0) {
        tkconfigure(canvas.var,bg=Gcol.var.value)
        assign("Gcol.var.tmp", Gcol.var.value, envir=env)
      }
    }
    ChangeColor.var.button <- tkbutton(GcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.var)
    tkgrid(tklabel(GcolFrame, text=.Facto_gettext("Color of the quantitative variables")),canvas.var,ChangeColor.var.button)
    
    Gcol.quanti.sup.value<-Gcol.quanti.sup
    canvas.quanti.sup <- tkcanvas(GcolFrame,width="80",height="25",bg=Gcol.quanti.sup.value)
    ChangeColor.quanti.sup <- function()
    {
      Gcol.quanti.sup.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Gcol.quanti.sup.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Gcol.quanti.sup.value)>0) {
        tkconfigure(canvas.quanti.sup,bg=Gcol.quanti.sup.value)
        assign("Gcol.quanti.sup.tmp", Gcol.quanti.sup.value, envir=env)
      }
    }
    ChangeColor.quanti.sup.button <- tkbutton(GcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.quanti.sup)
    if(!is.null(variableillu)) tkgrid(tklabel(GcolFrame, text=.Facto_gettext("color for supplementary quantitative variables")),canvas.quanti.sup,ChangeColor.quanti.sup.button)

    Gcol.quali.value<-Gcol.quali
    canvas.quali <- tkcanvas(GcolFrame,width="80",height="25",bg=Gcol.quali.value)
    ChangeColor.quali <- function()
    {
      Gcol.quali.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Gcol.quali.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Gcol.quali.value)>0) {
        tkconfigure(canvas.quali,bg=Gcol.quali.value)
        assign("Gcol.quali.tmp", Gcol.quali.value, envir=env)
      }
    }
    ChangeColor.quali.button <- tkbutton(GcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.quali)
    tkgrid(tklabel(GcolFrame, text=.Facto_gettext("color for qualitative variables")),canvas.quali,ChangeColor.quali.button)

    Gcol.quali.sup.value<-Gcol.quali.sup
    canvas.quali.sup <- tkcanvas(GcolFrame,width="80",height="25",bg=Gcol.quali.sup.value)
    ChangeColor.quali.sup <- function()
    {
      Gcol.quali.sup.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Gcol.quali.sup.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Gcol.quali.sup.value)>0) {
        tkconfigure(canvas.quali.sup,bg=Gcol.quali.sup.value)
        assign("Gcol.quali.sup.tmp", Gcol.quali.sup.value, envir=env)
      }
    }
    ChangeColor.quali.sup.button <- tkbutton(GcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.quali.sup)
    if(!is.null(variablefact)) tkgrid(tklabel(GcolFrame, text=.Facto_gettext("color for supplementary qualitative variables")),canvas.quali.sup,ChangeColor.quali.sup.button)

      GlabelFrame<-tkframe(PlotGrpeFrame,borderwidth=2)
      label.grpe.check<-tkcheckbutton(GlabelFrame)
      if (Glabel) label.grpe.checkValue<-tclVar("1")
      else label.grpe.checkValue<-tclVar("0")
      tkconfigure(label.grpe.check, variable=label.grpe.checkValue)
      tkgrid(tklabel(GlabelFrame, text=.Facto_gettext("Labels for the variables")),label.grpe.check)
       
      #mise en page des différents frames de PlotGrpeFrame
      tkgrid(GchoixFrame)
      tkgrid(GTitleFrame)
      tkgrid(GcolFrame)
      tkgrid(GlabelFrame)
      tkgrid(tklabel(PlotGrpeFrame, text=" "))
      
      
      ########################
      # construction de la partie graphique des variables
      PlotVarFrame<-tkframe(PlotWin2, borderwidth=5, relief="groove")
  
      WchoixFrame<-tkframe(PlotVarFrame,borderwidth=2)
      var.check<-tkcheckbutton(WchoixFrame)
      if(Wchoix) var.check.value<-tclVar("1")
      else var.check.value<-tclVar("0")
      tkconfigure(var.check, variable=var.check.value)
      tkgrid(tklabel(WchoixFrame, text=.Facto_gettext("Graph of the quantitative variables"), font=font2),var.check)
      tkgrid(tklabel(WchoixFrame, text="  "))
  
      WTitleFrame<-tkframe(PlotVarFrame,borderwidth=2)
      if (is.null(WTitle)) WTitre <- tclVar(" ")
      else WTitre<-tclVar(WTitle)
      WTitre.entry <-tkentry(WTitleFrame,width="40",textvariable=WTitre)
      tkgrid(tklabel(WTitleFrame,text=.Facto_gettext("Title of the graph")),WTitre.entry)
  
      WcosFrame<-tkframe(PlotVarFrame,borderwidth=2)
      WlimCosValue<-tclVar(paste(Wlim.cos))
      WlimCos.entry<-tkentry(WcosFrame, width=5, textvariable=WlimCosValue)
      tkgrid(tklabel(WcosFrame,text=.Facto_gettext("Draw variables with a cos2 >:")),WlimCos.entry)
  
      WlabelFrame<-tkframe(PlotVarFrame,borderwidth=2)
      label.var.check<-tkcheckbutton(WlabelFrame)
      if (Wlabel.var) label.var.checkValue<-tclVar("1")
      else label.var.checkValue<-tclVar("0")
      tkconfigure(label.var.check, variable=label.var.checkValue)
      tkgrid(tklabel(WlabelFrame, text=.Facto_gettext("Labels for the variables")),label.var.check)
      
    WcolFrame<-tkframe(PlotVarFrame,borderwidth=2)
    Wcol.var.value <- Wcol.var
    Wcanvas.var <- tkcanvas(WcolFrame,width="80",height="25",bg=Wcol.var.value)
    WChangeColor.var <- function()
    {
      Wcol.var.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Wcol.var.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Wcol.var.value)>0) {
        tkconfigure(Wcanvas.var,bg=Wcol.var.value)
        assign("Wcol.var.tmp", Wcol.var.value, envir=env)
      }
    }
    WChangeColor.var.button <- tkbutton(WcolFrame,text=.Facto_gettext("Change Color"),command=WChangeColor.var)
    tkgrid(tklabel(WcolFrame, text=.Facto_gettext("Color of the active variables")),Wcanvas.var,WChangeColor.var.button)
    
    Wcol.quanti.sup.value<-Wcol.quanti.sup
    Wcanvas.quanti.sup <- tkcanvas(WcolFrame,width="80",height="25",bg=Wcol.quanti.sup.value)
    WChangeColor.quanti.sup <- function()
    {
      Wcol.quanti.sup.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Wcol.quanti.sup.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Wcol.quanti.sup.value)>0) {
        tkconfigure(Wcanvas.quanti.sup,bg=Wcol.quanti.sup.value)
        assign("Wcol.quanti.sup.tmp", Wcol.quanti.sup.value, envir=env)
      }
    }
    WChangeColor.quanti.sup.button <- tkbutton(WcolFrame,text=.Facto_gettext("Change Color"),command=WChangeColor.quanti.sup)
    if(!is.null(variableillu)) tkgrid(tklabel(WcolFrame, text=.Facto_gettext("color for supplementary variables")),Wcanvas.quanti.sup,WChangeColor.quanti.sup.button)
         
      WinvisibleFrame<-tkframe(PlotVarFrame,borderwidth=2)
      inv.aucun.check<-tkradiobutton(WinvisibleFrame)
      inv.act.check<-tkradiobutton(WinvisibleFrame)
      inv.sup.check<-tkradiobutton(WinvisibleFrame)
      if(is.null(Winvisible)) inv.Value<-tclVar("aucun")
      else inv.Value<-tclVar(Winvisible) 
      tkconfigure(inv.aucun.check,variable=inv.Value,value="aucun")
      tkconfigure(inv.act.check,variable=inv.Value, value="actif")
      tkconfigure(inv.sup.check,variable=inv.Value, value="sup")
      tkgrid(tklabel(WinvisibleFrame, text=.Facto_gettext("Hide some elements:")), columnspan=6, sticky="w")
      tkgrid(tklabel(WinvisibleFrame, text="None"),inv.aucun.check, tklabel(WinvisibleFrame, text=.Facto_gettext("active variables")),inv.act.check, tklabel(WinvisibleFrame, text=.Facto_gettext("supplementary variables")),inv.sup.check, sticky="w")
        
      #mise en page des différents frames de PlotVarFrame
      tkgrid(WchoixFrame)
      tkgrid(WTitleFrame)
      tkgrid(WcolFrame)
      tkgrid(WcosFrame)
      tkgrid(WlabelFrame)
      tkgrid(WinvisibleFrame)            
      tkgrid(tklabel(PlotVarFrame, text=" "))
        
      ##########################  
      # construction de la partie graphique des individus
      PlotIndFrame<-tkframe(PlotWin, borderwidth=5, relief="groove")
      
      RchoixFrame<-tkframe(PlotIndFrame,borderwidth=2)
      ind.check<-tkcheckbutton(RchoixFrame)
      if(Rchoix) ind.check.value<-tclVar("1")
      else ind.check.value<-tclVar("0")
      tkconfigure(ind.check, variable=ind.check.value)
      tkgrid(tklabel(RchoixFrame, text=.Facto_gettext("Graph of the individuals"), font=font2),ind.check)
      tkgrid(tklabel(RchoixFrame, text=" "))
      
      RTitleFrame<-tkframe(PlotIndFrame,borderwidth=2)
      if (is.null(RTitle)) Titre <- tclVar(" ")
      else Titre<-tclVar(RTitle)
      Titre.entry <-tkentry(RTitleFrame,width="40",textvariable=Titre)
      tkgrid(tklabel(RTitleFrame,text=.Facto_gettext("Title of the graph")),Titre.entry)
      
      RlabelFrame<-tkframe(PlotIndFrame,borderwidth=2)
      label.indMoy.check<-tkcheckbutton(RlabelFrame)
      if (Rlabel.indMoy) label.indMoy.checkValue<-tclVar("1")
      else label.indMoy.checkValue<-tclVar("0") 
      tkconfigure(label.indMoy.check, variable=label.indMoy.checkValue)
      tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Labels for the mean individuals")),label.indMoy.check)
      label.quali.check<-tkcheckbutton(RlabelFrame)
      if (Rlabel.quali) label.quali.checkValue<-tclVar("1")
      else label.quali.checkValue<-tclVar("0")
      tkconfigure(label.quali.check, variable=label.quali.checkValue)
      tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Labels for the factors")), label.quali.check)
              
      RhabillageFrame<-tkframe(PlotIndFrame,borderwidth=2)
      listgraph<-tklistbox(RhabillageFrame,height=4, selectmode="single",exportselection="FALSE",yscrollcommand=function(...) tkset(scrgraph,...))
      scrgraph <-tkscrollbar(RhabillageFrame,repeatinterval=5,command=function(...)tkyview(listgraph,...))
      listgraph.nom<-c("ind")
      tkinsert(listgraph,"end","by.individual")
      if(Rhabillage=="ind") tkselection.set(listgraph,0)
      indice<-1      
      nbauxli<-c(tclvalue(tkcurselection(listfact)))
      nbaux<-unlist(strsplit(nbauxli,"\\ "))
      varaux = vars.fact[as.numeric(tkcurselection(listfact))+1]
      
      if (!is.null(variablefact)|(length(nbaux)>0)){
        for (j in 1:ncol(donnee)){
         if(vars[j] %in% c(variablefact,varaux)){
          tkinsert(listgraph,"end",vars[j])
          listgraph.nom<-c(listgraph.nom,vars[j])
          if(Rhabillage==vars[j]) tkselection.set(listgraph, indice)
          indice<-indice+1
        }}
      }
      tkgrid(tklabel(RhabillageFrame, text=.Facto_gettext("Select drawing for the individuals")))
      tkgrid(listgraph, scrgraph, sticky = "nw")
      tkgrid.configure(scrgraph, sticky = "wns")
      tkgrid.configure(listgraph, sticky = "ew")
      
      
      RinvisibleFrame<-tkframe(PlotIndFrame,borderwidth=2)
      inv.ind.check<-tkcheckbutton(RinvisibleFrame)
      if ("ind" %in% Rinvisible) inv.ind.checkValue<-tclVar("1")
      else inv.ind.checkValue<-tclVar("0") 
      inv.ind.sup.check<-tkcheckbutton(RinvisibleFrame)
      if ("ind.sup" %in% Rinvisible) inv.ind.sup.checkValue<-tclVar("1")
      else inv.ind.sup.checkValue<-tclVar("0") 
      inv.quali.check<-tkcheckbutton(RinvisibleFrame)      
      if ("quali" %in% Rinvisible) inv.quali.checkValue<-tclVar("1")
      else inv.quali.checkValue<-tclVar("0") 
      tkconfigure(inv.ind.check, variable=inv.ind.checkValue)
      tkconfigure(inv.ind.sup.check, variable=inv.ind.sup.checkValue)
      tkconfigure(inv.quali.check, variable=inv.quali.checkValue)            
      tkgrid(tklabel(RinvisibleFrame, text=.Facto_gettext("Hide some elements:")), columnspan=6, sticky="w")
      tkgrid(tklabel(RinvisibleFrame, text="ind"),inv.ind.check, tklabel(RinvisibleFrame, text="ind sup"),inv.ind.sup.check, tklabel(RinvisibleFrame, text="quali"),inv.quali.check, sticky="w")
      
            
      RlimFrame<-tkframe(PlotIndFrame,borderwidth=2)
      if(is.null(RXlimInd)) XlimIndMin<-tclVar("")
      else XlimIndMin<-tclVar(paste(RXlimInd[1]))
      XlimIndMin.entry <-tkentry(RlimFrame,width="5",textvariable=XlimIndMin)
      if (is.null(RXlimInd)) XlimIndMax<- tclVar("")
      else XlimIndMax<-tclVar(paste(RXlimInd[1]))
      XlimIndMax.entry <-tkentry(RlimFrame,width="5",textvariable=XlimIndMax)
      tkgrid(tklabel(RlimFrame,text=.Facto_gettext("x limits of the graph:")),XlimIndMin.entry, XlimIndMax.entry)
      if(is.null(RYlimInd)) YlimIndMin<- tclVar("")
      else YlimIndMin<-tclVar(paste(RYlimInd[1]))
      YlimIndMin.entry <-tkentry(RlimFrame,width="5",textvariable=YlimIndMin)
      if (is.null(RYlimInd)) YlimIndMax<- tclVar("")
      else YlimIndMax<-tclVar(paste(RYlimInd[2]))
      YlimIndMax.entry <-tkentry(RlimFrame,width="5",textvariable=YlimIndMax)
      tkgrid(tklabel(RlimFrame,text=.Facto_gettext("y limits of the graph:")),YlimIndMin.entry,YlimIndMax.entry)
  
  
      #mise en page des différents frames de PlotIndFrame
      tkgrid(RchoixFrame)
      tkgrid(RTitleFrame)
      tkgrid(RlabelFrame)
      tkgrid(RinvisibleFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))
      tkgrid(RhabillageFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))      
      tkgrid(RlimFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))
      
              
      #mise en page de plotWin
      subOKCancelHelp(PlotWin, "plot.FAMD")
      tkgrid(PlotGrpeFrame)
      tkgrid(PlotVarFrame)
      tkgrid(PlotIndFrame, PlotWin2, sticky="ns")
      tkgrid(subButtonsFrame, sticky="ew", columnspan=2)
    }

    PlotFrame<-tkframe(IlluFrame)
    Plot.but<-tkbutton(PlotFrame, textvariable=.PlotLabel, command=OnPlot, borderwidth=3)
    tkgrid(Plot.but, sticky="ew")
  }) 



    #! fonction associée au bouton Appliquer, execute sans détruire l'interface graphique
  OnAppliquer<-function()
  {
      #liste de l'ensemble des variables créées
      #sur la fenêtre top
#      listQuantiAct
#      listQuantiIllu
#      listQualiAct
#      listQualiIllu
#      resu.val
#      ncp.val
      #pour les individus illustratifs
#      individuillu
      #pour l'affichage
#      Rpropre
#          Rgroupe
#          Rindividu
#        Rindsup
#      Rquantisummary
#      Rquanti
#      Rquantisup  
#          Rqualisummary
#          Rquali
#          Rqualisup
#          Rdescdim
      # pour les graphiques
#      Gchoix
#      GTitle
#      GAxeGrpe
#      Glabel
#        
#      Rchoix
#      RTitle
#      Rlabel.indMoy
#      Rlabel.quali  
#      Rhabillage
#      Rinvisible
#      Rpartial
#      RpartialSouris
#      Rchrono
#      RXlimInd
#      RYlimInd
#    
#      Wchoix
#      WTitle
#      WAxeVar
#      Wlabel.var
#      Winvisible
#      Wlim.cos
#
#      Axe 





    # if(length(as.numeric(tkcurselection(listdesc)))<2) varActives<-listdesc.nom
    # else varActives<-listdesc.nom[as.numeric(tkcurselection(listdesc))+1]
    # varActives <- varActives[!(varActives%in%variableillu)]                       #varActives is the list of selected active variables








  
    # récupération des paramètres de la fenêtre principale
    nom.res<-tclvalue(resu.val)
    if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) justDoIt(paste('remove (',nom.res,')'))       #if object res already exists, it's removed
    ncp<-as.numeric(tclvalue(ncp.val))
    Axe<-c(as.numeric(tclvalue(Axe1)), as.numeric(tclvalue(Axe2)))
    # nbitemlist<-c(tclvalue(tkcurselection(listdesc)))
    # nbitem<-unlist(strsplit(nbitemlist,"\\ "))
    # nbitemlist.q<-c(tclvalue(tkcurselection(listfact)))
    # nbitem.q<-unlist(strsplit(nbitemlist.q,"\\ "))

 
    # gestion du tableau de données pour l'FAMD
    variables <- variables.q <- NULL
    # if (length(nbitem)>0) variables <- vars.desc[as.numeric(tkcurselection(listdesc))+1]
    # if (length(nbitem.q)>0) variables.q <- vars.fact[as.numeric(tkcurselection(listfact))+1]

if(length(as.numeric(tkcurselection(listdesc)))<1) variables <- vars.desc
else variables<-vars.desc[as.numeric(tkcurselection(listdesc))+1]
variables <- variables[!(variables%in%variableillu)]                       #varActives is the list of selected active variables

if(length(as.numeric(tkcurselection(listfact)))<1) variables.q <- vars.fact
else variables.q<-vars.fact[as.numeric(tkcurselection(listfact))+1]
variables.q  <- variables.q[!(variables.q%in%variablefact)]    

    allvariables = c(variables,variables.q,variableillu,variablefact)
    num.group.sup<-NULL
    if (length(variableillu)+length(variablefact)>0) num.group.sup <- ((length(variables)+length(variables.q)+1):length(allvariables))
    #construction du tableau de données.FAMD
      if(!is.null(individuillu)) {
        ind.actif<-rows[-which(rows %in% individuillu)]
        commande.data<-paste(activeDataSet(),'.FAMD', '<-', activeDataSet(),'[c("', paste(ind.actif, collapse='", "'), '", "', paste(individuillu, collapse='", "'), '"),', sep='')
      }
      else commande.data<-paste(activeDataSet(),'.FAMD', '<-', activeDataSet(),'[,', sep='')
      commande.data<-paste(commande.data,' c("',paste(allvariables, collapse='", "'), '")]',sep='')
      
      justDoIt(commande.data)
      logger(commande.data)
      donnee.depart<-activeDataSet()
      activeDataSet(paste(activeDataSet(),'.FAMD', sep=""))

      # gestion de la commande réalisant l'AFM     
      commande.FAMD<-paste(nom.res, '<-FAMD(', activeDataSet(),sep='')
      if(!is.null(individuillu)) commande.FAMD<-paste(commande.FAMD, ', ind.sup=', nrow(get(getRcmdr(".activeDataSet")))-length(individuillu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))),sep='')
      commande.FAMD<-paste(commande.FAMD, ', ncp=', ncp,sep='')
      if (!is.null(num.group.sup)) commande.FAMD<-paste(commande.FAMD, ', sup.var=',num.group.sup[1],':',num.group.sup[length(num.group.sup)],sep='')
      commande.FAMD<-paste(commande.FAMD, ', graph=FALSE)',sep='')
      justDoIt(commande.FAMD)
      logger(commande.FAMD)
	  justDoIt(paste(nom.res,'$call$call <-',deparse(commande.FAMD),sep=""))

	  	      #Commande de la fonction HCPC

    if(Rclassif==TRUE){
      commande.hcpc<-paste(nom.res,'.hcpc', '<-HCPC(', nom.res, ' ,nb.clust=', Rmeth, ',consol=', Rconsolid,',min=', Rminhcpc,',max=',Rmaxhcpc,',graph=', Rgraphhcpc, ')', sep="")
    justDoIt(commande.hcpc)
    logger(commande.hcpc)      
      if(Rreshcpc==TRUE){
        doItAndPrint(paste(nom.res,'.hcpc$data.clust[,ncol(res.hcpc$data.clust),drop=F]', sep=""))
        doItAndPrint(paste(nom.res,'.hcpc$desc.var', sep=""))
        doItAndPrint(paste(nom.res,'.hcpc$desc.axes', sep=""))
        doItAndPrint(paste(nom.res,'.hcpc$desc.ind', sep=""))
      }        
    }

      #gestion des graphiques   
      if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) {if (get(nom.res)$eig[1,2]==100) doItAndPrint(paste('"No graph can be plot: data are unidimensional"'))}
      if((Gchoix)&(length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0)){
      if (get(nom.res)$eig[1,2]!=100) {
        commande.plotG<-paste('plot.FAMD(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="var", new.plot=TRUE, lab.var=', Glabel, sep="")
        commande.plotG <- paste(commande.plotG, ',col.hab = c(',sep='')
        auxi = 0
        if (length(variables)>0){
          commande.plotG <- paste(commande.plotG, 'rep("',Gcol.var,'",',length(variables),')',sep='')
          auxi = 1
        }
        if (length(variables.q)>0){
          if (auxi==1) commande.plotG <- paste(commande.plotG, ',',sep='')
          commande.plotG <- paste(commande.plotG, 'rep("',Gcol.quali,'",',length(variables.q),')',sep='')
          auxi=1
        }
        if (length(variableillu)>0){
          if (auxi==1) commande.plotG <- paste(commande.plotG, ',',sep='')
          commande.plotG <- paste(commande.plotG, 'rep("',Gcol.quanti.sup,'",',length(variableillu),')',sep='')
        }
        if (length(variablefact)>0){
          if (auxi==1) commande.plotG <- paste(commande.plotG, ',',sep='')
          commande.plotG <- paste(commande.plotG, 'rep("',Gcol.quali.sup,'",',length(variablefact),')',sep='')
        }
        commande.plotG <- paste(commande.plotG, ')',sep='')
        if (is.null(GTitle)) commande.plotG <- paste(commande.plotG,')', sep="")
        else {
          if (GTitle ==" ") commande.plotG <- paste(commande.plotG,')', sep="")
          else commande.plotG <- paste(commande.plotG,', title="', GTitle,'")', sep="")
        }
        justDoIt(commande.plotG)
        logger(commande.plotG)
      }}
            
      if((Wchoix)&(length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0)&(length(variables)>0)){
      if (get(nom.res)$eig[1,2]!=100) {
        commande.plotW<-paste('plot.FAMD(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="quanti", new.plot=TRUE, lab.var=', Wlabel.var, ', lim.cos2.var=', Wlim.cos, sep="")

        if (!is.null(Winvisible)) {
          commande.plotW<-paste(commande.plotW, ', invisible=c("', paste(Winvisible, collapse='", "'),'")', sep='')
          if(Winvisible=="actif") commande.plotW<-paste(commande.plotW, ', col.hab=c(rep("', Wcol.quanti.sup,'",length(rownames(', nom.res, '$quanti.var.sup[[1]]))),rep("', Wcol.var, '",length(rownames(', nom.res,'$quanti.var[[1]]))))', sep='')
          else commande.plotW<-paste(commande.plotW, ', col.hab=c(rep("', Wcol.var, '",length(rownames(', nom.res,'$quanti.var[[1]]))),rep("', Wcol.quanti.sup,'",length(rownames(', nom.res, '$quanti.var.sup[[1]]))))', sep='')
        }
        if(is.null(Winvisible)) commande.plotW<-paste(commande.plotW, ', col.hab=c(rep("', Wcol.var, '",length(rownames(', nom.res,'$quanti.var[[1]]))),rep("', Wcol.quanti.sup,'",length(rownames(', nom.res, '$quanti.var.sup[[1]]))))', sep='')

        if (is.null(WTitle)) commande.plotW <- paste(commande.plotW,')', sep="")
        else {
          if (WTitle ==" ") commande.plotW <- paste(commande.plotW,')', sep="")
          else commande.plotW <- paste(commande.plotW,', title="', WTitle,'")', sep="")
        }
        justDoIt(commande.plotW)
        logger(commande.plotW)
      }}
      
      if((Rchoix)&(length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0)){
      if (get(nom.res)$eig[1,2]!=100) {
        if ((Rhabillage!="none") & (Rhabillage!="ind")) {
          Rhabillage<-which(colnames(get(getRcmdr(".activeDataSet")))==Rhabillage)
          if(length(Rhabillage)==0) Rhabillage<-"none"
        }
        if (Rhabillage=="none") Rhabillage<-paste('"', Rhabillage, '"', sep="")
        if (Rhabillage=="ind") Rhabillage<-paste('"', Rhabillage, '"', sep="")
        
          commande.plotI<-paste('plot.FAMD(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="ind", new.plot=TRUE, lab.ind=', Rlabel.indMoy, ', lab.var=', Rlabel.quali, ', habillage=', Rhabillage, sep="") 
          if (!is.null(RXlimInd)) commande.plotI<-paste(commande.plotI, ', xlim=c(', paste(RXlimInd, collapse=", "), ')', sep='')
          if (!is.null(RYlimInd)) commande.plotI<-paste(commande.plotI, ', ylim=c(', paste(RYlimInd, collapse=", "), ')', sep='')
          if (!is.null(Rinvisible)) commande.plotI<-paste(commande.plotI, ', invisible=c("', paste(Rinvisible, collapse='", "'),'")', sep='')
          if (is.null(RTitle)) commande.plotI <- paste(commande.plotI,')', sep="")
          else {
            if (RTitle ==" ") commande.plotI <- paste(commande.plotI,')', sep="")
            else commande.plotI <- paste(commande.plotI,', title="', RTitle,'")', sep="")
          }
        }
        justDoIt(commande.plotI)
        logger(commande.plotI)
      }
      
      # gestion de l'édition de certains resultats
    doItAndPrint(paste('summary(',nom.res,', nb.dec = 3, nbelements=10, nbind = 10, ncp = 3, file="")', sep=""))
    if (RFichier==""){
      if(Rpropre) doItAndPrint(paste( nom.res, '$eig', sep=""))
      if(Rgroupe) doItAndPrint(paste( nom.res, '$group', sep=""))
      if(Rindividu) doItAndPrint(paste( nom.res, '$ind', sep=""))
      if(Rindsup) doItAndPrint(paste( nom.res, '$ind.sup', sep=""))
#      if(Rquantisummary) doItAndPrint(paste( nom.res, '$summary.quanti', sep=""))
      if(Rquanti) doItAndPrint(paste( nom.res, '$quanti.var', sep=""))      
      if(Rquantisup) doItAndPrint(paste( nom.res, '$quanti.var.sup', sep=""))
#      if(Rqualisummary) doItAndPrint(paste( nom.res, '$summary.quali', sep=""))
      if(Rquali) doItAndPrint(paste( nom.res, '$quali.var', sep=""))      
      if(Rqualisup) doItAndPrint(paste( nom.res, '$quali.var.sup', sep=""))
      if(Rdescdim) doItAndPrint(paste('dimdesc(', nom.res, ', axes=1:',ncp,')', sep=""))
    }
    else {
      Fich = RFichier
      if (substr(Fich,1,1)!='"') Fich = paste('"',Fich,sep='')
      if (substr(Fich,nchar(Fich),nchar(Fich))!='"') Fich = paste(Fich,'"',sep='')
      append = FALSE
      if(Rpropre){
        doItAndPrint(paste('write.infile(', nom.res, '$eig, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rgroupe){
        doItAndPrint(paste('write.infile(', nom.res, '$group, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rindividu){
        doItAndPrint(paste('write.infile(', nom.res, '$ind, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rindsup){
        doItAndPrint(paste('write.infile(', nom.res, '$ind.sup, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
#      if(Rquantisummary){
#        doItAndPrint(paste('write.infile(', nom.res, '$summary.quanti, file =',Fich,',append=',append,')', sep=""))
#        append = TRUE
#      }
      if(Rquanti){
        doItAndPrint(paste('write.infile(', nom.res, '$quanti.var, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rquantisup){
        doItAndPrint(paste('write.infile(', nom.res, '$quanti.var.sup, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
#      if(Rqualisummary){
#        doItAndPrint(paste('write.infile(', nom.res, '$summary.quali, file =',Fich,',append=',append,')', sep=""))
#        append = TRUE
#      }
      if(Rquali){
        doItAndPrint(paste('write.infile(', nom.res, '$quali.var, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rqualisup){
        doItAndPrint(paste('write.infile(', nom.res, '$quali.var.sup, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rdescdim) doItAndPrint(paste('write.infile(dimdesc(', nom.res, ', axes=1:',ncp,'), file =',Fich,',append=',append,')', sep=""))
    }

      # Re-chargement du tableau de départ et supression du tableau temporaire
      activeDataSet(donnee.depart)
      justDoIt(paste('remove(',activeDataSet(),'.FAMD)',sep=""))
      logger(paste('remove(',activeDataSet(),'.FAMD)',sep=""))   
  }

  #! fonction HCPC
  
  Hcpc.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    .HcpcLabel<-tclVar(paste(firstLabel, "", sep=" "))    
    compteur.hcpc<-0
    Rclassif<-0
    Rmeth <- -1
    Rconsolid<-0 
    Rgraphhcpc<-1  
    Rreshcpc<-0
    Rminhcpc<-3
    Rmaxhcpc<-10

    OnHCPC <- function()
    {

      HcpcWin<-tktoplevel()
      tkwm.title(HcpcWin, .Facto_gettext("HCPC options"))

      onOKHcpc <- function()
      {
        assign("compteur.hcpc", compteur.hcpc+1, envir=env) 
        if(compteur.hcpc>0) tclvalue(.HcpcLabel)<-paste(label, "", sep=" ")
        tkconfigure(Hcpc.but, fg="blue")
      
        if(tclvalue(methValue)=="0") assign("Rmeth", 0, envir=env)
        else assign("Rmeth", -1, envir=env)
        
        if(tclvalue(consolidValue)=="1") assign("Rconsolid",TRUE, envir=env)
        else assign("Rconsolid",FALSE,envir=env)
        
        if(tclvalue(graphhcpcValue)=="1") assign("Rgraphhcpc",TRUE,envir=env)
        else assign("Rgraphhcpc",FALSE,envir=env)
        
        if(tclvalue(reshcpcValue)=="1") assign("Rreshcpc",TRUE,envir=env)
        else assign("Rreshcpc",FALSE,envir=env)        
        
        assign("Rminhcpc",as.numeric(tclvalue(minhcpc)),envir=env)
        assign("Rmaxhcpc",as.numeric(tclvalue(maxhcpc)),envir=env)
        
        assign("Rclassif",TRUE,envir=env)
        tkdestroy(HcpcWin)
      }
      OKHcpc.but<-tkbutton(HcpcWin, text="OK", width=8,command=onOKHcpc)

      onCancelHcpc <- function()
      {
        assign("Rclassif",FALSE,envir=env)
        tkdestroy(HcpcWin)
      }
      CancelHcpc.but<-tkbutton(HcpcWin, text="Cancel", width=8,command=onCancelHcpc)
 
      tkgrid(tklabel(HcpcWin, text=""))
      tkgrid(tklabel(HcpcWin, text = .Facto_gettext("Hierarchical Clustering on Principal Components"), fg = "darkred"), column=1, columnspan = 8, sticky = "ew")      

      meth1 <- tkradiobutton (HcpcWin)
      meth1.lab <- tklabel(HcpcWin,text=.Facto_gettext("interactive"))
      meth2 <- tkradiobutton (HcpcWin)
      meth2.lab <- tklabel(HcpcWin,text=.Facto_gettext("automatic"))
      methValue <- tclVar(paste(Rmeth))
      meth.lab <- tklabel(HcpcWin,text=.Facto_gettext("Choice of the number of clusters:"))
      tkconfigure(meth1,variable=methValue,value="0")
      tkconfigure(meth2,variable=methValue,value="-1")

      minmaxhcpc.label<-tklabel(HcpcWin,text=.Facto_gettext("The optimal number of clusters is chosen between:"))

      minhcpc<-tclVar(paste(Rminhcpc))
      maxhcpc<-tclVar(paste(Rmaxhcpc))
      minhcpc.entry <-tkentry(HcpcWin,width="3",textvariable=minhcpc)
      maxhcpc.entry <-tkentry(HcpcWin,width="3",textvariable=maxhcpc)

      consolid.lab <- tklabel(HcpcWin,text=.Facto_gettext("Consolidate clusters"))
      consolid.check <- tkcheckbutton(HcpcWin)
      if(Rconsolid) consolidValue<-tclVar("1")
      else consolidValue<-tclVar("0")      
      tkconfigure(consolid.check,variable=consolidValue)
      
      graphhcpc.lab <- tklabel(HcpcWin,text=.Facto_gettext("Print graphs"))
      graphhcpc.check <- tkcheckbutton(HcpcWin)
      if(Rgraphhcpc) graphhcpcValue <- tclVar("1")
      else graphhcpcValue <- tclVar("0")
      tkconfigure(graphhcpc.check,variable=graphhcpcValue)   

      reshcpc.lab <- tklabel(HcpcWin,text=.Facto_gettext("Print results for clusters"))
      reshcpc.check <- tkcheckbutton(HcpcWin)
      if(Rreshcpc) reshcpcValue<-tclVar("1")
      else reshcpcValue <- tclVar("0")
      tkconfigure(reshcpc.check,variable=reshcpcValue)          
    
      tkgrid(tklabel(HcpcWin,text=.Facto_gettext("Options for the clustering"), fg = "blue"), column=1, columnspan=8, sticky="we")
      tkgrid(tklabel(HcpcWin,text=""))
      tkgrid(tklabel(HcpcWin,text=sprintf(.Facto_gettext("Clustering is performed on the first %s dimensions of FAMD"),tclvalue(ncp.val))),column=1,columnspan=4,sticky="w")   #text which takes the nb of dimensions chosen in the main window
      tkgrid(tklabel(HcpcWin,text=.Facto_gettext("(Modify in the main options to change this number)")),column=1,columnspan=4,sticky="w")
      tkgrid(tklabel(HcpcWin,text=""))  
      tkgrid(meth.lab,meth1.lab,meth1)
      tkgrid(meth2.lab,meth2)
      tkgrid(tklabel(HcpcWin,text=""))
      tkgrid(minmaxhcpc.label,minhcpc.entry , maxhcpc.entry)
      tkgrid(tklabel(HcpcWin,text=""))      
      tkgrid(consolid.lab,consolid.check)
      tkgrid(graphhcpc.lab,graphhcpc.check)
      tkgrid(reshcpc.lab,reshcpc.check) 
      tkgrid(tklabel(HcpcWin,text=""))     
      tkgrid(OKHcpc.but, CancelHcpc.but)
      tkgrid(tklabel(HcpcWin, text=""))
      
      tkgrid.configure(minmaxhcpc.label,meth.lab,consolid.lab,graphhcpc.lab,reshcpc.lab,column=1,columnspan=4,sticky="w")
      tkgrid.configure(minhcpc.entry,column=7,columnspan=1,sticky="e")
      tkgrid.configure(maxhcpc.entry,column=8,columnspan=1,sticky="w")
      tkgrid.configure(meth1,meth2,consolid.check,graphhcpc.check,reshcpc.check,column=8,sticky="e")
      tkgrid.configure(meth1.lab,column=6,columnspan=2,sticky="w")
      tkgrid.configure(meth2.lab,column=6,columnspan=2,sticky="w") 
      tkgrid.configure(OKHcpc.but,column=2,columnspan=1,sticky="w")
      tkgrid.configure(CancelHcpc.but,column=6,columnspan=1,sticky="e")

      tkgrid.columnconfigure(HcpcWin,0, minsize=3)
      tkgrid.columnconfigure(HcpcWin,5, minsize=5)
      tkgrid.columnconfigure(HcpcWin,8, minsize=3) 

}      
    Hcpc2Frame<-tkframe(HcpcFrame)
    Hcpc.but<-tkbutton(Hcpc2Frame, textvariable=.HcpcLabel, command=OnHCPC, borderwidth=3)
    tkgrid(Hcpc.but, sticky="ew")
}) 

    #! fonction associée au bouton OK, execute et détruit l'interface graphique
  onOK<-function()
  {
    OnAppliquer()
    tkdestroy(top)     
  }



#                   Création de la fenêtre top                                 #
################################################################################
  

   # création de tous les boutons d'options dans IlluFrame
  IlluFrame<- tkframe(top, borderwidth=2)


    # mise en page de IlluFrame
  Fillu.funct(label=.Facto_gettext("Supplementary factors"), firstLabel=.Facto_gettext("Supplementary factors"))
  Dillu.funct(label=.Facto_gettext("Supplementary quantitative variables"), firstLabel=.Facto_gettext("Supplementary quantitative variables"))
  Iillu.funct(label=.Facto_gettext("Supplementary individuals"), firstLabel=.Facto_gettext("Supplementary individuals"))    
  PLOT.FAMD(label=.Facto_gettext("Graphical options"), firstLabel=.Facto_gettext("Graphical options"))
  Sortie.funct(label=.Facto_gettext("Outputs"), firstLabel=.Facto_gettext("Outputs"))
  tkgrid(DilluFrame, FilluFrame, IilluFrame, columnspan=7)
  tkgrid(tklabel(IlluFrame, text=""))
  tkgrid(PlotFrame, SortieFrame, columnspan=7)

  tkgrid.configure(DilluFrame, column=1, columnspan=1)
  tkgrid.configure(PlotFrame, column=2, columnspan=2,sticky="we")
  tkgrid.configure(SortieFrame, column=4, columnspan=2,sticky="ew")
  tkgrid.configure(FilluFrame, column=3, columnspan=1)
  tkgrid.configure(IilluFrame, column=5, columnspan=1)
  tkgrid.columnconfigure(IlluFrame,0, minsize=25)
  tkgrid.columnconfigure(IlluFrame,2, minsize=40)  
  tkgrid.columnconfigure(IlluFrame,4, minsize=25)  


    # création des options dans OptionFrame  
  OptionFrame<-tkframe(top, borderwidth=2, relief="groove")
  resu.lab<-tklabel(OptionFrame,text=.Facto_gettext("Name of the result object: "))
  resu.val<-tclVar("res")
  resu<-tkentry(OptionFrame,width=10,textvariable=resu.val)
#  reduit.lab<-tklabel(OptionFrame,text=.Facto_gettext("Scale the quantative variables: "))
#  reduit.check <- tkcheckbutton(OptionFrame)
#  reduitValue <- tclVar("1")
#  tkconfigure(reduit.check,variable=reduitValue)
  ncp.lab<-tklabel(OptionFrame,text=.Facto_gettext("Number of dimensions: "))
  ncp.val<-tclVar("5") 
  ncp<-tkentry(OptionFrame,width=5,textvariable=ncp.val)
  Axe.label<-tklabel(OptionFrame,text=.Facto_gettext("Select the dimensions for the graphs:"))
  Axe1<-tclVar("1")
  Axe2<-tclVar("2")
  Axe1.entry <-tkentry(OptionFrame,width="5",textvariable=Axe1)
  Axe2.entry <-tkentry(OptionFrame,width="5",textvariable=Axe2)
  
    # mise en page de OptionFrame
  tkgrid(tklabel(OptionFrame,text=.Facto_gettext("Main options"), fg = "darkred"), columnspan=8, sticky="we") 
  tkgrid(tklabel(OptionFrame,text="")) 
#  tkgrid(reduit.lab,reduit.check)
  tkgrid(ncp.lab, ncp)
  tkgrid(Axe.label,Axe1.entry , Axe2.entry, sticky="w")
  tkgrid(resu.lab, resu)
# tkgrid.configure(ncp.lab, reduit.lab, resu.lab, Axe.label, column=1, columnspan=4, sticky="w")
# tkgrid.configure(ncp, resu, reduit.check, column=6, columnspan=2, sticky="e")
  tkgrid.configure(ncp.lab, resu.lab, Axe.label, column=1, columnspan=4, sticky="w")
  tkgrid.configure(ncp, resu, column=6, columnspan=2, sticky="e")
  tkgrid.configure(Axe1.entry, column=6, columnspan=1, sticky="w")
  tkgrid.configure(Axe2.entry, column=7, columnspan=1, sticky="e")
  tkgrid.columnconfigure(OptionFrame,0, minsize=25)
  tkgrid.columnconfigure(OptionFrame,5, minsize=40)
  tkgrid.columnconfigure(OptionFrame,8, minsize=25)

    #Frame pour HCPC
  HcpcFrame<-tkframe(top, borderwidth=2)
  Hcpc.funct(label=.Facto_gettext("Perform Clustering after FAMD"), firstLabel=.Facto_gettext("Perform Clustering after FAMD")) 
  tkgrid(Hcpc2Frame, columnspan=7)
  tkgrid.configure(Hcpc2Frame,column=4, columnspan=1)
  
  appliquer.but<-tkbutton(top, text=.Facto_gettext("Apply"),width=12,command=OnAppliquer, borderwidth=3, fg="#690f96")
  OKCancelHelp(helpSubject="FAMD",reset="Reinitializ.funct")

  # Mise en page de top
  tkgrid(tklabel(top, text=.Facto_gettext("Factor Analysis of Mixed Data (FAMD)"),font=fontheading), columnspan=3)
  tkgrid(tklabel(top,text=""))
  tkgrid(listFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text=""))
  tkgrid(IlluFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text=""))    
  tkgrid(OptionFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text=""))   
  tkgrid(HcpcFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text="")) # Ligne de blanc
  # tkgrid(appliquer.but, column=1, columnspan=1)
  # tkgrid(tklabel(top,text=""))
  # tkgrid(buttonsFrame, column=1, columnspan=1, sticky="ew" )  
  tkgrid(buttonsFrame, appliquer.but)
  tkgrid.configure(buttonsFrame, column=1,sticky="e")
  tkgrid.configure(appliquer.but, column=2,sticky="w")
}
