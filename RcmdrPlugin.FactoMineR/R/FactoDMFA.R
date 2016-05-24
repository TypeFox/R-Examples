FactoDMFA <-
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
  tkwm.title(top,.Facto_gettext("DMFA"))
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
    lab2 = tklabel(listFrame,text=.Facto_gettext("Select the group variable"),fg="blue")
    lab3 = tklabel(listFrame,text="      ")
    tkgrid(lab1,lab3,lab2)
    tkgrid.configure(lab1,column=1, columnspan=2, sticky = "nw")
    tkgrid.configure(lab2,column=4, columnspan=2, sticky = "ne")
    tkgrid.configure(lab3,column=3, columnspan=1, sticky = "n")

    listdesc<-tklistbox(listFrame,selectmode="extended",exportselection=FALSE,yscrollcommand=function(...) tkset(scr,...))
    scr <- tkscrollbar(listFrame,repeatinterval=5,command=function(...)tkyview(listdesc,...))
    tkselection.set(listdesc,0)
    listfact<-tklistbox(listFrame,selectmode="single",exportselection=FALSE,yscrollcommand=function(...) tkset(scrfact,...))
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
      tkwm.title(IilluWin,.Facto_gettext("Select supplementary individuals"))
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
      tkgrid(tklabel(IilluWin, text = .Facto_gettext("Select supplementary individuals"), fg = "blue"), column=1, columnspan = 1, sticky = "ew")
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
    FactoDMFA()
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
    RXc<-FALSE
    Rvar.partiel<-FALSE
    Rquanti<-FALSE
    Rquantisup<-FALSE    
    RCov<-FALSE
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
        
        if(tclvalue(Xc.Value)=="1") assign("RXc", TRUE, envir=env)
        else assign("RXc", FALSE, envir=env)
        
        if(tclvalue(var.partielValue)=="1") assign("Rvar.partiel", TRUE, envir=env)
        else assign("Rvar.partiel", FALSE, envir=env)
        
        if(tclvalue(quantiValue)=="1") assign("Rquanti", TRUE, envir=env)
        else assign("Rquanti", FALSE, envir=env)
        
        if(tclvalue(quantisupValue)=="1") assign("Rquantisup", TRUE, envir=env)
        else assign("Rquantisup", FALSE, envir=env)
        
        if(tclvalue(CovValue)=="1") assign("RCov", TRUE, envir=env)
        else assign("RCov", FALSE, envir=env)
        
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

      Xc.lab<-tklabel(SortieWin,text=.Facto_gettext("Group of individuals"))
        Xc.check <- tkcheckbutton(SortieWin)
        if(RXc) Xc.Value <- tclVar("1")
        else Xc.Value <- tclVar("0")
        tkconfigure(Xc.check,variable=Xc.Value)
        
      var.partiel.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for partial variables"))
        var.partiel.check <- tkcheckbutton(SortieWin)
        if(Rvar.partiel) var.partielValue <- tclVar("1")
        else var.partielValue <- tclVar("0")
        tkconfigure(var.partiel.check,variable=var.partielValue)
      
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
      
      Cov.lab<-tklabel(SortieWin,text=.Facto_gettext("Covariance matrices by group"))
        Cov.check <- tkcheckbutton(SortieWin)
        if(RCov) CovValue <- tclVar("1")
        else CovValue <- tclVar("0")
        tkconfigure(Cov.check,variable=CovValue)
      
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
      tkgrid(quanti.lab,quanti.check,sticky="w")      
      tkgrid(var.partiel.lab,var.partiel.check,sticky="w")
      if (!is.null(variableillu)) tkgrid(quantisup.lab,quantisup.check,sticky="w")
      if (!is.null(variablefact)) tkgrid(qualisup.lab,qualisup.check,sticky="w")
      tkgrid(Cov.lab,Cov.check,sticky="w")
      tkgrid(Xc.lab,Xc.check,sticky="w")
      tkgrid(descdim.lab,descdim.check,sticky="w")
      tkgrid(tklabel(SortieWin, text = " "))
      tkgrid(RFichierFrame)
      tkgrid(tklabel(SortieWin, text = " "))
      tkgrid(SortieOK.but)
   }
    
    SortieFrame<-tkframe(IlluFrame)
    Sortie.but<-tkbutton(SortieFrame, textvariable=.SortieLabel, command=OnSortie, borderwidth=3)
    tkgrid(Sortie.but, sticky="ew")
  })



  #! fonction pour la gestion des options graphiques 
  PLOT.DMFA<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.graph<-0
    .PlotLabel<-tclVar(paste(firstLabel, "", sep=" "))
    #déclaration des variables
    Gchoix<-TRUE
    GTitle<-NULL
    GAxeGrpe<-c(1,2)
    Glabel<-"all"
        
    Rchoix<-TRUE
    RTitle<-NULL
    Rlabel.indMoy<-"ind"
    Rlabel.quali<-NULL    
    Rhabillage<-"none"
    Rinvisible<-NULL
    RXlimInd<-NULL
    RYlimInd<-NULL
    
    Wchoix=TRUE
    WTitle<-NULL
    WAxeVar<-c(1,2)
    Wlabel.var<-"all"
    #Wcol.quanti.sup<-Wcol.quanti.sup.tmp<-"blue"
    #Wcol.var<-Wcol.var.tmp<-"black"
    #Winvisible<-NULL
    #Wlim.cos<-0.
    
        
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

          label.tmp.grpe<-tclvalue(label.grpe.checkValue)
          if(label.tmp.grpe==1) assign("Glabel", "all", envir=env)
          else assign("Glabel", "none", envir=env)
        }
                    
        # gestion des entrées de la partie graphique des variables
        if(tclvalue(var.check.value)==1) assign("Wchoix", TRUE, envir=env)
        else assign("Wchoix", FALSE, envir=env)

        if(Wchoix) {
          if (tclvalue(WTitre)==" ") assign("WTitle", NULL, envir=env)
          assign("WTitle", tclvalue(WTitre), envir=env)

          #assign("Wlim.cos", tclvalue(WlimCosValue), envir=env)

          label.tmp.var<-tclvalue(label.var.checkValue)
          if(label.tmp.var==1) assign("Wlabel.var", "all", envir=env)
          else assign("Wlabel.var", "none", envir=env)

          #assign("Wcol.var", Wcol.var.tmp, envir=env)
          #assign("Wcol.quanti.sup", Wcol.quanti.sup.tmp, envir=env)

          #if(tclvalue(inv.Value)=="aucun")  assign("Winvisible", NULL, envir=env)
          #else assign("Winvisible", tclvalue(inv.Value), envir=env)
        }

        # gestion des entrées de la partie graphique des individus
        if(tclvalue(ind.check.value)==1) assign("Rchoix", TRUE, envir=env)
        else assign("Rchoix", FALSE, envir=env)

        if(Rchoix) {
          if (tclvalue(Titre)==" ") assign("RTitle", NULL, envir=env)
          assign("RTitle", tclvalue(Titre), envir=env)

          label.tmp.indMoy<-tclvalue(label.indMoy.checkValue)
          label.tmp.quali<-tclvalue(label.quali.checkValue)
          if(label.tmp.indMoy==1) assign("Rlabel.indMoy", "ind", envir=env)
          else assign("Rlabel.indMoy", NULL, envir=env)
          if(label.tmp.quali==1) assign("Rlabel.quali", "quali", envir=env)
          else assign("Rlabel.quali", NULL, envir=env)
          
##          habillage.tmp<-listgraph.nom[as.numeric(tkcurselection(listgraph))+1]
##          if(length(habillage.tmp)==0) assign("Rhabillage","none", envir=env)
##          else assign("Rhabillage", habillage.tmp, envir=env)

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
      tkgrid(tklabel(GchoixFrame, text=.Facto_gettext("Graph the groups"), font=font2),grpe.check)
      tkgrid(tklabel(GchoixFrame, text="  "))
  
      GTitleFrame<-tkframe(PlotGrpeFrame,borderwidth=2)
      if (is.null(GTitle)) GTitre <- tclVar(" ")
      else GTitre<-tclVar(GTitle)
      GTitre.entry <-tkentry(GTitleFrame,width="40",textvariable=GTitre)
      tkgrid(tklabel(GTitleFrame,text=.Facto_gettext("Title of the graph")),GTitre.entry)
    
      GlabelFrame<-tkframe(PlotGrpeFrame,borderwidth=2)
      label.grpe.check<-tkcheckbutton(GlabelFrame)
      if (Glabel=="all") label.grpe.checkValue<-tclVar("1")
      else label.grpe.checkValue<-tclVar("0")
      tkconfigure(label.grpe.check, variable=label.grpe.checkValue)
      tkgrid(tklabel(GlabelFrame, text=.Facto_gettext("Labels for the groups")),label.grpe.check)
       
      #mise en page des différents frames de PlotGrpeFrame
      tkgrid(GchoixFrame)
      tkgrid(GTitleFrame)
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
  
      #WcosFrame<-tkframe(PlotVarFrame,borderwidth=2)
      #WlimCosValue<-tclVar(paste(Wlim.cos))
      #WlimCos.entry<-tkentry(WcosFrame, width=5, textvariable=WlimCosValue)
      #tkgrid(tklabel(WcosFrame,text=.Facto_gettext("Draw variables with a cos2 >:")),WlimCos.entry)
  
      WlabelFrame<-tkframe(PlotVarFrame,borderwidth=2)
      label.var.check<-tkcheckbutton(WlabelFrame)
      if (Wlabel.var=="all") label.var.checkValue<-tclVar("1")
      else label.var.checkValue<-tclVar("0")
      tkconfigure(label.var.check, variable=label.var.checkValue)
      tkgrid(tklabel(WlabelFrame, text=.Facto_gettext("Labels for the variables")),label.var.check)
      
#    WcolFrame<-tkframe(PlotVarFrame,borderwidth=2)
#    Wcol.var.value <- Wcol.var
#    Wcanvas.var <- tkcanvas(WcolFrame,width="80",height="25",bg=Wcol.var.value)
#    WChangeColor.var <- function()
#    {
#      Wcol.var.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Wcol.var.value,title=.Facto_gettext("Choose a color")))
#      if (nchar(Wcol.var.value)>0) {
#        tkconfigure(Wcanvas.var,bg=Wcol.var.value)
#        assign("Wcol.var.tmp", Wcol.var.value, envir=env)
#      }
#    }
#    WChangeColor.var.button <- tkbutton(WcolFrame,text=.Facto_gettext("Change Color"),command=WChangeColor.var)
#    tkgrid(tklabel(WcolFrame, text=.Facto_gettext("Color of the active variables")),Wcanvas.var,WChangeColor.var.button)
    
#    Wcol.quanti.sup.value<-Wcol.quanti.sup
#    Wcanvas.quanti.sup <- tkcanvas(WcolFrame,width="80",height="25",bg=Wcol.quanti.sup.value)
#    WChangeColor.quanti.sup <- function()
#    {
#      Wcol.quanti.sup.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Wcol.quanti.sup.value,title=.Facto_gettext("Choose a color")))
#      if (nchar(Wcol.quanti.sup.value)>0) {
#        tkconfigure(Wcanvas.quanti.sup,bg=Wcol.quanti.sup.value)
#        assign("Wcol.quanti.sup.tmp", Wcol.quanti.sup.value, envir=env)
#      }
#    }
#    WChangeColor.quanti.sup.button <- tkbutton(WcolFrame,text=.Facto_gettext("Change Color"),command=WChangeColor.quanti.sup)
#    if(!is.null(variableillu)) tkgrid(tklabel(WcolFrame, text=.Facto_gettext("color for supplementary variables")),Wcanvas.quanti.sup,WChangeColor.quanti.sup.button)
         
#      WinvisibleFrame<-tkframe(PlotVarFrame,borderwidth=2)
#      inv.aucun.check<-tkradiobutton(WinvisibleFrame)
#      inv.act.check<-tkradiobutton(WinvisibleFrame)
#      inv.sup.check<-tkradiobutton(WinvisibleFrame)
#      if(is.null(Winvisible)) inv.Value<-tclVar("aucun")
#      else inv.Value<-tclVar(Winvisible) 
#      tkconfigure(inv.aucun.check,variable=inv.Value,value="aucun")
#      tkconfigure(inv.act.check,variable=inv.Value, value="actif")
#      tkconfigure(inv.sup.check,variable=inv.Value, value="sup")
#      tkgrid(tklabel(WinvisibleFrame, text=.Facto_gettext("Hide some elements:")), columnspan=6, sticky="w")
#      tkgrid(tklabel(WinvisibleFrame, text="None"),inv.aucun.check, tklabel(WinvisibleFrame, text=.Facto_gettext("active variables")),inv.act.check, tklabel(WinvisibleFrame, text=.Facto_gettext("supplementary variables")),inv.sup.check, sticky="w")
        
      #mise en page des différents frames de PlotVarFrame
      tkgrid(WchoixFrame)
      tkgrid(WTitleFrame)
      #tkgrid(WcolFrame)
      #tkgrid(WcosFrame)
      tkgrid(WlabelFrame)
      #tkgrid(WinvisibleFrame)            
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
      if (!is.null(Rlabel.indMoy)) label.indMoy.checkValue<-tclVar("1")
      else label.indMoy.checkValue<-tclVar("0") 
      tkconfigure(label.indMoy.check, variable=label.indMoy.checkValue)
      tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Labels for the mean individuals")),label.indMoy.check)
      label.quali.check<-tkcheckbutton(RlabelFrame)
      if (!is.null(Rlabel.quali)) label.quali.checkValue<-tclVar("1")
      else label.quali.checkValue<-tclVar("0")
      tkconfigure(label.quali.check, variable=label.quali.checkValue)
      if (!is.null(variablefact)) tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Labels for the factors")), label.quali.check)
              
##      RhabillageFrame<-tkframe(PlotIndFrame,borderwidth=2)
##      listgraph<-tklistbox(RhabillageFrame,height=4, selectmode="single",exportselection="FALSE",yscrollcommand=function(...) tkset(scrgraph,...))
##      scrgraph <-tkscrollbar(RhabillageFrame,repeatinterval=5,command=function(...)tkyview(listgraph,...))
##      listgraph.nom<-c("ind")
##      tkinsert(listgraph,"end","by.individual")
##      if(Rhabillage=="ind") tkselection.set(listgraph,0)
##      indice<-1      
##      nbauxli<-c(tclvalue(tkcurselection(listfact)))
##      nbaux<-unlist(strsplit(nbauxli,"\\ "))
##      varaux = vars.fact[as.numeric(tkcurselection(listfact))+1]
##      
##      if (!is.null(variablefact)|(length(nbaux)>0)){
##        for (j in 1:ncol(donnee)){
##         if(vars[j] %in% c(variablefact,varaux)){
##          tkinsert(listgraph,"end",vars[j])
##          listgraph.nom<-c(listgraph.nom,vars[j])
##          if(Rhabillage==vars[j]) tkselection.set(listgraph, indice)
##          indice<-indice+1
##        }}
##      }
##      tkgrid(tklabel(RhabillageFrame, text=.Facto_gettext("Select drawing for the individuals")))
##      tkgrid(listgraph, scrgraph, sticky = "nw")
##      tkgrid.configure(scrgraph, sticky = "wns")
##      tkgrid.configure(listgraph, sticky = "ew")
      
      
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
      if (!is.null(variablefact)) {
        tkgrid(tklabel(RinvisibleFrame, text=.Facto_gettext("Hide some elements:")), columnspan=6, sticky="w")
#        tkgrid(tklabel(RinvisibleFrame, text="ind"),inv.ind.check, tklabel(RinvisibleFrame, text="ind sup"),inv.ind.sup.check, tklabel(RinvisibleFrame, text="quali"),inv.quali.check, sticky="w")
        tkgrid(tklabel(RinvisibleFrame, text="ind"),inv.ind.check, tklabel(RinvisibleFrame, text="quali"),inv.quali.check, sticky="w")
      }
      
            
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
##      tkgrid(RhabillageFrame)
##      tkgrid(tklabel(PlotIndFrame, text=" "))      
      tkgrid(RlimFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))
      
              
      #mise en page de plotWin
      subOKCancelHelp(PlotWin, "plot.DMFA")
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
#        RXc
#      Rvar.partiel
#      Rquanti
#      Rquantisup  
#          RCov
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
  
    # récupération des paramètres de la fenêtre principale
    nom.res<-tclvalue(resu.val)
    if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) justDoIt(paste('remove (',nom.res,')'))
    ncp<-as.numeric(tclvalue(ncp.val))
    reduction<-TRUE
    if(tclvalue(reduitValue)=="0") reduction<-FALSE
    Axe<-c(as.numeric(tclvalue(Axe1)), as.numeric(tclvalue(Axe2)))

    nbitemlist<-c(tclvalue(tkcurselection(listdesc)))
    nbitem<-unlist(strsplit(nbitemlist,"\\ "))
    nbitemlist.q<-c(tclvalue(tkcurselection(listfact)))
    nbitem.q<-unlist(strsplit(nbitemlist.q,"\\ "))

 
    # gestion du tableau de données pour DMFA
    variables <- variables.q <- NULL
    variables <- vars.desc[as.numeric(tkcurselection(listdesc))+1]
    variables.q <- vars.fact[as.numeric(tkcurselection(listfact))+1]
    allvariables = c(variables,variables.q,variableillu,variablefact)
    num.group.sup<-NULL
    if (length(variableillu)+length(variablefact)>0) num.group.sup <- ((length(variables)+length(variables.q)+1):length(allvariables))
    #construction du tableau de données.DMFA
      if(!is.null(individuillu)) {
        ind.actif<-rows[-which(rows %in% individuillu)]
        commande.data<-paste(activeDataSet(),'.DMFA', '<-', activeDataSet(),'[c("', paste(ind.actif, collapse='", "'), '", "', paste(individuillu, collapse='", "'), '"),', sep='')
      }
      else commande.data<-paste(activeDataSet(),'.DMFA', '<-', activeDataSet(),'[,', sep='')
      commande.data<-paste(commande.data,' c("',paste(allvariables, collapse='", "'), '")]',sep='')
      
      justDoIt(commande.data)
      logger(commande.data)
      donnee.depart<-activeDataSet()
      activeDataSet(paste(activeDataSet(),'.DMFA', sep=""))

      # gestion de la commande réalisant l'AFMD     
      commande.DMFA<-paste(nom.res, '<-DMFA(', activeDataSet(), ', num.fact=',length(variables)+1,', ncp=', ncp,', scale.unit = ', reduction,sep='')
      if(!is.null(individuillu)) commande.DMFA<-paste(commande.DMFA, ', ind.sup=', nrow(get(getRcmdr(".activeDataSet")))-length(individuillu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))),sep='')
      if (!is.null(variableillu)) commande.DMFA<-paste(commande.DMFA, ', quanti.sup =',length(variables)+2,':',length(variables)+1+length(variableillu),sep='')
      if (!is.null(variablefact)) commande.DMFA<-paste(commande.DMFA, ', quali.sup =',length(allvariables)-length(variablefact)+1,':',length(allvariables),sep='')
      commande.DMFA<-paste(commande.DMFA, ', graph=FALSE)',sep='')
      justDoIt(commande.DMFA)
      logger(commande.DMFA)

      #gestion des graphiques   
      if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) {if (get(nom.res)$eig[1,2]==100) doItAndPrint(paste('"No graph can be plot: data are unidimensional"'))}
      if((Gchoix)&(length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0)){
      if (get(nom.res)$eig[1,2]!=100) {
        commande.plotG<-paste('plot.DMFA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="group", new.plot=TRUE', ',label="', Glabel,'"', sep='')
#        if (!Glabel) commande.plotG <- paste(commande.plotG,', label="none"', sep="")
#        else commande.plotG <- paste(commande.plotG,', label="', Glabel,'"', sep="")
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
#        commande.plotW<-paste('plot.DMFA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="var", col.var="', Wcol.var, '", col.quanti.sup="', Wcol.quanti.sup, '", label=c("', paste(Wlabel.var, collapse='", ")'), '"), lim.cos2.var=', Wlim.cos, sep="")
#        commande.plotW<-paste('plot.DMFA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="var", col.var="', Wcol.var, '", col.quanti.sup="', Wcol.quanti.sup, '", label=c("', paste(Wlabel.var, collapse='", ")'), '")', sep="")
        commande.plotW<-paste('plot.DMFA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="var", new.plot=TRUE, label=c("', paste(Wlabel.var, collapse='", ")'), '")', sep="")        

#        if (!is.null(Winvisible)) commande.plotW<-paste(commande.plotW, ', invisible=c("', paste(Winvisible, collapse='", "'),'")', sep='')        
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
        

          commande.plotI<-paste('plot.DMFA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="ind", new.plot=TRUE', sep="") 
          if (!is.null(RXlimInd)) commande.plotI<-paste(commande.plotI, ', xlim=c(', paste(RXlimInd, collapse=", "), ')', sep='')
          if (!is.null(RYlimInd)) commande.plotI<-paste(commande.plotI, ', ylim=c(', paste(RYlimInd, collapse=", "), ')', sep='')
          if (!is.null(Rinvisible)) commande.plotI<-paste(commande.plotI, ', invisible=c("', paste(Rinvisible, collapse='", "'),'")', sep='')
          if (!is.null(Rlabel.indMoy) & !is.null(Rlabel.quali)) commande.plotI<-paste(commande.plotI, ', label = c("',Rlabel.indMoy, '","', Rlabel.quali,'")',sep='')
          if (!is.null(Rlabel.indMoy) & is.null(Rlabel.quali)) commande.plotI<-paste(commande.plotI, ', label = c("',Rlabel.indMoy, '")',sep='')  
          if (is.null(Rlabel.indMoy) & !is.null(Rlabel.quali)) commande.plotI<-paste(commande.plotI, ', label = c("',Rlabel.quali,'")',sep='')   
          if (is.null(Rlabel.indMoy) & is.null(Rlabel.quali)) commande.plotI<-paste(commande.plotI, ', label = "none"',sep='')                          
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
    if (RFichier==""){
      if(Rpropre) doItAndPrint(paste(nom.res, '$eig', sep=""))
      if(Rgroupe) doItAndPrint(paste(nom.res, '$group', sep=""))
      if(Rindividu) doItAndPrint(paste(nom.res, '$ind', sep=""))
      if(Rquanti) doItAndPrint(paste(nom.res, '$var', sep=""))      
      if(Rvar.partiel) doItAndPrint(paste(nom.res, '$var.partiel', sep=""))
      if(Rquantisup) doItAndPrint(paste(nom.res, '$quanti.sup', sep=""))
      if(Rqualisup) doItAndPrint(paste(nom.res, '$quali.sup', sep=""))
      if(RCov) doItAndPrint(paste( nom.res, '$Cov', sep=""))
      if(RXc) doItAndPrint(paste( nom.res, '$Xc', sep=""))
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
      if(Rquanti){
        doItAndPrint(paste('write.infile(', nom.res, '$var, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rvar.partiel){
        doItAndPrint(paste('write.infile(', nom.res, '$var.partiel, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rquantisup){
        doItAndPrint(paste('write.infile(', nom.res, '$quanti.sup, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rqualisup){
        doItAndPrint(paste('write.infile(', nom.res, '$quali.sup, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(RCov){
        doItAndPrint(paste('write.infile(', nom.res, '$Cov, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(RXc){
        doItAndPrint(paste('write.infile(', nom.res, '$Xc, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rdescdim) doItAndPrint(paste('write.infile(dimdesc(', nom.res, ', axes=1:',ncp,'), file =',Fich,',append=',append,')', sep=""))
    }

      # Re-chargement du tableau de départ et supression du tableau temporaire
      activeDataSet(donnee.depart)
      justDoIt(paste('remove(',activeDataSet(),'.DMFA)',sep=""))
      logger(paste('remove(',activeDataSet(),'.DMFA)',sep=""))   
  }


    #! fonction associée au bouton OK, execute et détruit l'interface graphique
  onOK<-function()
  {
    OnAppliquer()
    tkdestroy(top)     
  }



#                   Création de la fenêtre top                                 #
################################################################################
##  top<-tktoplevel(borderwidth=10)
##  tkwm.title(top,.Facto_gettext("DMFA"))
##  tkwm.geometry(top, "-50+50")
##  
##  # définition des polices
##  font2<-tkfont.create(family="times",size=12,weight="bold")
##  fontheading<-tkfont.create(family="times",size=11,weight="bold")
##
##  # récupération du jeu de données actif
##  donnee<-get(getRcmdr(".activeDataSet"))
##  vars<-colnames(donnee)
##  rows<-rownames(donnee)
  

   # création de tous les boutons d'options dans IlluFrame
  IlluFrame<- tkframe(top, borderwidth=2)


    # mise en page de IlluFrame
  Fillu.funct(label=.Facto_gettext("Supplementary factors"), firstLabel=.Facto_gettext("Supplementary factors"))
  Dillu.funct(label=.Facto_gettext("Supplementary quantitative variables"), firstLabel=.Facto_gettext("Supplementary quantitative variables"))
  Iillu.funct(label=.Facto_gettext("Supplementary individuals"), firstLabel=.Facto_gettext("Supplementary individuals"))    
  PLOT.DMFA(label=.Facto_gettext("Graphical options"), firstLabel=.Facto_gettext("Graphical options"))
  Sortie.funct(label=.Facto_gettext("Outputs"), firstLabel=.Facto_gettext("Outputs"))
  tkgrid(DilluFrame, FilluFrame, columnspan=7)
  tkgrid.configure(DilluFrame, column=1, columnspan=1)
  tkgrid.configure(FilluFrame, column=3, columnspan=1)
  tkgrid.columnconfigure(IlluFrame,0, minsize=25)
  tkgrid.columnconfigure(IlluFrame,2, minsize=40)  
  tkgrid.columnconfigure(IlluFrame,4, minsize=25)  
  tkgrid(tklabel(IlluFrame, text=""))

  tkgrid(PlotFrame, SortieFrame, columnspan=7)
  tkgrid.configure(PlotFrame, column=1, columnspan=1)
  tkgrid.configure(SortieFrame, column=3, columnspan=1)


    # création des options dans OptionFrame  
  OptionFrame<-tkframe(top, borderwidth=2, relief="groove")
  resu.lab<-tklabel(OptionFrame,text=.Facto_gettext("Name of the result object:"))
  resu.val<-tclVar("res")
  resu<-tkentry(OptionFrame,width=10,textvariable=resu.val)
  reduit.lab<-tklabel(OptionFrame,text=.Facto_gettext("Scale the quantitative variables:"))
  reduit.check <- tkcheckbutton(OptionFrame)
  reduitValue <- tclVar("1")
  tkconfigure(reduit.check,variable=reduitValue)
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
  tkgrid(reduit.lab,reduit.check)
  tkgrid(ncp.lab, ncp)
  tkgrid(Axe.label,Axe1.entry , Axe2.entry, sticky="w")
  tkgrid(resu.lab, resu)
  tkgrid.configure(ncp.lab, reduit.lab, resu.lab, Axe.label, column=1, columnspan=4, sticky="w")
  tkgrid.configure(ncp, resu, reduit.check, column=6, columnspan=2, sticky="e")
  tkgrid.configure(Axe1.entry, column=6, columnspan=1, sticky="w")
  tkgrid.configure(Axe2.entry, column=7, columnspan=1, sticky="e")
  tkgrid.columnconfigure(OptionFrame,0, minsize=25)
  tkgrid.columnconfigure(OptionFrame,5, minsize=40)
  tkgrid.columnconfigure(OptionFrame,8, minsize=25)

  appliquer.but<-tkbutton(top, text=.Facto_gettext("Apply"),width=12,command=OnAppliquer, borderwidth=3, fg="#690f96")
  OKCancelHelp(helpSubject="DMFA",reset="Reinitializ.funct")

  # Mise en page de top
  tkgrid(tklabel(top, text=.Facto_gettext("Dual Multiple Factor Analysis (DMFA)"),font=fontheading), columnspan=3)
  tkgrid(tklabel(top,text=""))
  tkgrid(listFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text=""))
  tkgrid(IlluFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text=""))    
  tkgrid(OptionFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text=""))   
  # tkgrid(appliquer.but, column=1, columnspan=1)
  # tkgrid(tklabel(top,text=""))
  # tkgrid(buttonsFrame, column=1, columnspan=1, sticky="ew" )  
  tkgrid(buttonsFrame, appliquer.but)
  tkgrid.configure(buttonsFrame, column=1,sticky="e")
  tkgrid.configure(appliquer.but, column=2,sticky="w")
}
