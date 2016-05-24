FactoGPA <-
function()
{

# Fonction pour la gestion des noms ############################################

  nom.correct<-function(text, liste=NULL)
  {
    text<-chartr("^\ ", "...", text)
    if(!is.null(liste)) {
      while(text %in% liste) {
        text<-paste(text, ".bis", sep="")
      }
    }
    return(text)  
  }
################################################################################

#    Création des fonctions pour les options via nouvelle fenêtre graphique   


  #! suppression de groupes quantitatif
  supprimeQuanti.funct<-defmacro(label, expr=
  {
    env<-environment()
    OnSGQ<-function()
    {
      grpeActSupprime<-listQuantiAct.nom[as.numeric(tkcurselection(listQuantiAct))+1]
      if(length(grpeActSupprime)>=1)
      { 
        listQuantiAct.nom.tmp<-listQuantiAct.nom[-which(listQuantiAct.nom %in% grpeActSupprime)]
        assign("listQuantiAct.nom",listQuantiAct.nom.tmp, envir=env)
        tkdelete(listQuantiAct,"0","end")
        if(length(listQuantiAct.nom)>=1) {
          for (grpe in listQuantiAct.nom) tkinsert(listQuantiAct, "end", grpe)
        }  
      }
      nb.grpe<-length(listQuantiAct.nom) 
      if (nb.grpe>1) {
        tclvalue(label.quantiFrame.var)<-paste(nb.grpe, .Facto_gettext("groups"), sep=" ")
        tkconfigure(label.quantiFrame)
      }  
      else
      {
        tclvalue(label.quantiFrame.var)<-paste("0", .Facto_gettext("group"), sep=" ")
        tkconfigure(label.quantiFrame)
      }  
      
    }
    
    SupGpeQuantiFrame<-tkframe(ListeQuantiFrame)
    SupGpeQuanti.but<-tkbutton(SupGpeQuantiFrame, textvariable=tclVar(label), command=OnSGQ, borderwidth=3)
   tkgrid(SupGpeQuanti.but, sticky="ew")  
  })


  #! Ajout d'un groupe quantitatif
  ajoutQuanti.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.GQ<-1
    .AjoutQuantiLabel<-tclVar(paste(firstLabel, "", sep=" "))
    OnAGQ<-function()
    {
      AjoutGpeQuantiWin<-tktoplevel()
      tkwm.title(AjoutGpeQuantiWin,.Facto_gettext("Definition of a group"))
      
      #création de la fonction AGA.OK
      AGQ.OK<-function()
      {
        assign("compteur.GQ", compteur.GQ+1, envir=env)
        nom.groupe<-nom.correct(tclvalue(nomGrpeQuanti.val), liste=c(listQuantiAct.nom))
        if (nom.groupe=="") tkmessageBox(message=.Facto_gettext("Name for the group"), icon="warning", type="ok") 
        else {
          varGroupe<-listVarQuanti.nom[as.numeric(tkcurselection(listVarQuanti))+1]
          if (length(varGroupe)>=1) {
            assign(paste(nom.groupe,".var", sep=""), c(varGroupe), envir=env)
            tkinsert(listQuantiAct,"end",nom.groupe)
            assign("listQuantiAct.nom", c(listQuantiAct.nom, nom.groupe),envir=env)
            if (length(listQuantiAct.nom)==1) tclvalue(label.quantiFrame.var)<-paste(.Facto_gettext("1 group"), sep=" ")
            else tclvalue(label.quantiFrame.var)<-paste(length(listQuantiAct.nom) , .Facto_gettext("groups"), sep=" ")
            tkconfigure(label.quantiFrame)  
            tkdestroy(AjoutGpeQuantiWin)
          }
        }  
      }
      
      
      # choix du nom du groupe
      nomGrpeQuanti.lab<-tklabel(AjoutGpeQuantiWin,text=.Facto_gettext("Name of the group: "))
      nomGrpeQuanti.val<-tclVar(paste("Gc", compteur.GQ, sep=""))
      nomGrpeQuanti<-tkentry(AjoutGpeQuantiWin,width=15,textvariable=nomGrpeQuanti.val)
                  
      # création de la liste pour le choix des variables acives
      listVarQuanti<-tklistbox(AjoutGpeQuantiWin, selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrVarQuanti,...))
      scrVarQuanti<-tkscrollbar(AjoutGpeQuantiWin,repeatinterval=5,command=function(...)tkyview(listVarQuanti,...))
      listVarQuanti.nom<-NULL
      for (i in (1:ncol(donnee))) {
        if (is.numeric(donnee[,i])) {
            tkinsert(listVarQuanti,"end",vars[i])
            listVarQuanti.nom<-c(listVarQuanti.nom, vars[i])
        }
      }
      
      AGQ.but<-tkbutton(AjoutGpeQuantiWin, text="OK", width=16, command=AGQ.OK)
      
      tkgrid(nomGrpeQuanti.lab, nomGrpeQuanti)
      tkgrid.configure(nomGrpeQuanti.lab, column=0, columnspan=2, sticky="w")
      tkgrid.configure(nomGrpeQuanti, column=2, columnspan=3)      
      tkgrid(tklabel(AjoutGpeQuantiWin, text=""))      
      tkgrid(tklabel(AjoutGpeQuantiWin, text = .Facto_gettext("Select the variables for the group"), fg = "blue"), column=0, columnspan = 5, sticky = "w")
      tkgrid(listVarQuanti, scrVarQuanti, sticky = "nw")
      tkgrid.configure(scrVarQuanti, sticky = "wns", column=4,columnspan=1)
      tkgrid.configure(listVarQuanti, sticky = "ew", column=0, columnspan=4)
      tkgrid(tklabel(AjoutGpeQuantiWin, text=""))
      tkgrid(AGQ.but, column=2,columnspan=1, sticky="ew")
      tkgrid.columnconfigure(AjoutGpeQuantiWin,0, minsize=55)
      tkgrid.columnconfigure(AjoutGpeQuantiWin,1, minsize=55)
      tkgrid.columnconfigure(AjoutGpeQuantiWin,2, minsize=55)
      tkgrid.columnconfigure(AjoutGpeQuantiWin,3, minsize=55)
      tkgrid.columnconfigure(AjoutGpeQuantiWin,4, minsize=55)             
   }
   GpeQuantiFrame<-tkframe(ListeQuantiFrame)
   AjoutGpeQuanti.but<-tkbutton(GpeQuantiFrame, textvariable=.AjoutQuantiLabel, command=OnAGQ, borderwidth=3)
   tkgrid(AjoutGpeQuanti.but, sticky="ew")
  })    
  


  #! Modification d'un groupe quantitatif
  modifQuanti.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    .ModifQuantiLabel<-tclVar(paste(firstLabel, "", sep=" "))
    OnMGQ<-function() {
      ModifGpeQuantiWin<-tktoplevel()
      tkwm.title(ModifGpeQuantiWin,.Facto_gettext("Modification of a group"))
      
      #création de la fonction AGA.OK
      MGQ.OK<-function() {
        nom.groupe<-nom.correct(tclvalue(nomModifGrpeQuanti.val), liste=c(listQuantiAct.nom))
        if (nom.groupe=="") tkmessageBox(message=.Facto_gettext("Name for the group"), icon="warning", type="ok") 
        else {
          listQuantiAct.nom.tmp<-listQuantiAct.nom[-which(listQuantiAct.nom== grpeAModifier)]
          assign("listQuantiAct.nom",listQuantiAct.nom.tmp, envir=env)
          tkdelete(listQuantiAct,"0","end")
          for (grpe in listQuantiAct.nom) tkinsert(listQuantiAct, "end", grpe)
          
          varGroupe<-listModifVarQuanti.nom[as.numeric(tkcurselection(listModifVarQuanti))+1]
          if (length(varGroupe)>=1) {
            assign(paste(nom.groupe,".var", sep=""), c(varGroupe), envir=env)
            tkinsert(listQuantiAct,"end",nom.groupe)
            assign("listQuantiAct.nom", c(listQuantiAct.nom, nom.groupe),envir=env)
            tkdestroy(ModifGpeQuantiWin)
          }
        }  
      }
      
      
      if(length(as.numeric(tkcurselection(listQuantiAct)))>=1)  grpeAModifier<-listQuantiAct.nom[as.numeric(tkcurselection(listQuantiAct))+1][1]
      else {
        tkdestroy(ModifGpeQuantiWin)
        return()
      }  
        
      eval(parse(text=paste("grpeAModifier.var<-",paste(grpeAModifier,".var", sep=""),sep="")))
      # choix du nom du groupe
      nomModifGrpeQuanti.lab<-tklabel(ModifGpeQuantiWin,text=.Facto_gettext("Name of the group: "))
      nomModifGrpeQuanti.val<-tclVar(grpeAModifier)
      nomModifGrpeQuanti<-tkentry(ModifGpeQuantiWin,width=15,textvariable=nomModifGrpeQuanti.val)
            
      # création de la liste pour le choix des variables acives
      listModifVarQuanti<-tklistbox(ModifGpeQuantiWin, selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrModifVarQuanti,...))
      scrModifVarQuanti<-tkscrollbar(ModifGpeQuantiWin,repeatinterval=5,command=function(...)tkyview(listModifVarQuanti,...))
      listModifVarQuanti.nom<-NULL
      indice.num<-0
      for (i in (1:ncol(donnee))) {
        if (is.numeric(donnee[,i])) {
            tkinsert(listModifVarQuanti,"end",vars[i])
            listModifVarQuanti.nom<-c(listModifVarQuanti.nom, vars[i])
            if(vars[i] %in% grpeAModifier.var) tkselection.set(listModifVarQuanti, indice.num)
            indice.num<-indice.num+1
        }
      }
      
      MGQ.but<-tkbutton(ModifGpeQuantiWin, text="OK", width=16, command=MGQ.OK)
      
      tkgrid(nomModifGrpeQuanti.lab, nomModifGrpeQuanti)
      tkgrid.configure(nomModifGrpeQuanti.lab, column=0, columnspan=2, sticky="w")
      tkgrid.configure(nomModifGrpeQuanti, column=2, columnspan=3)      
      tkgrid(tklabel(ModifGpeQuantiWin, text=""))      
      tkgrid(tklabel(ModifGpeQuantiWin, text = .Facto_gettext("Select the variables of the group"), fg = "blue"), column=0, columnspan = 5, sticky = "w")
      tkgrid(listModifVarQuanti, scrModifVarQuanti, sticky = "nw")
      tkgrid.configure(scrModifVarQuanti, sticky = "wns", column=4,columnspan=1)
      tkgrid.configure(listModifVarQuanti, sticky = "ew", column=0, columnspan=4)
      tkgrid(tklabel(ModifGpeQuantiWin, text=""))
      tkgrid(MGQ.but, column=2,columnspan=1, sticky="ew")
      tkgrid.columnconfigure(ModifGpeQuantiWin,0, minsize=55)
      tkgrid.columnconfigure(ModifGpeQuantiWin,1, minsize=55)
      tkgrid.columnconfigure(ModifGpeQuantiWin,2, minsize=55)
      tkgrid.columnconfigure(ModifGpeQuantiWin,3, minsize=55)
      tkgrid.columnconfigure(ModifGpeQuantiWin,4, minsize=55)             
   }
   ModifGpeQuantiFrame<-tkframe(ListeQuantiFrame)
   ModifGpeQuanti.but<-tkbutton(ModifGpeQuantiFrame, textvariable=.ModifQuantiLabel, command=OnMGQ, borderwidth=3)
   tkgrid(ModifGpeQuanti.but, sticky="ew")
  })    
  
  
    #! fonction pour la réinitialisation des paramètre
  Reinitializ.funct<-function()
  {
    tkdestroy(top)
    FactoGPA()
  }


  #! fonction pour le choix des éléments de sortie
  Sortie.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.sortie<-0
    #déclaration des variables
    RFichier <- ""
    Rdep<-FALSE
    RRVs<-FALSE
    Rsimi<-FALSE
    Rscaling<-FALSE
    Rconsensus<-TRUE    
    RPANOVA<-FALSE
    RXfin<-FALSE
    Rcorrelations<-FALSE          
    RRV<-TRUE

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
        
        if(tclvalue(depValue)=="1") assign("Rdep", TRUE, envir=env)
        else assign("Rdep", FALSE, envir=env)
        
        if(tclvalue(RVsValue)=="1") assign("RRVs", TRUE, envir=env)
        else assign("RRVs", FALSE, envir=env)
        
        if(tclvalue(simi.Value)=="1") assign("Rsimi", TRUE, envir=env)
        else assign("Rsimi", FALSE, envir=env)
        
        if(tclvalue(scalingValue)=="1") assign("Rscaling", TRUE, envir=env)
        else assign("Rscaling", FALSE, envir=env)
        
        if(tclvalue(consensusValue)=="1") assign("Rconsensus", TRUE, envir=env)
        else assign("Rconsensus", FALSE, envir=env)
        
        if(tclvalue(PANOVAValue)=="1") assign("RPANOVA", TRUE, envir=env)
        else assign("RPANOVA", FALSE, envir=env)
        
        if(tclvalue(XfinValue)=="1") assign("RXfin", TRUE, envir=env)
        else assign("RXfin", FALSE, envir=env)
        
        if(tclvalue(correlationsValue)=="1") assign("Rcorrelations", TRUE, envir=env)
        else assign("Rcorrelations", FALSE, envir=env)
        
        if(tclvalue(RVValue)=="1") assign("RRV", TRUE, envir=env)
        else assign("RRV", FALSE, envir=env)
        
        if (tclvalue(Fichier)=="") assign("RFichier", NULL, envir=env)
        assign("RFichier", tclvalue(Fichier), envir=env)

        tkdestroy(SortieWin)
      
      }
      
        RV.lab<-tklabel(SortieWin, text=.Facto_gettext("RV coefficients between partial configurations"))
      RV.check<-tkcheckbutton(SortieWin)
      if(RRV) RVValue<-tclVar("1")
      else RVValue<-tclVar("0")
      tkconfigure(RV.check,variable=RVValue)

        RVs.lab <-tklabel(SortieWin, text=.Facto_gettext("Standardized RV coefficients between partial configurations"))
        RVs.check <- tkcheckbutton(SortieWin)
        if(RRVs) RVsValue <- tclVar("1")
        else RVsValue <- tclVar("0")
        tkconfigure(RVs.check,variable=RVsValue)

      simi.lab<-tklabel(SortieWin,text=.Facto_gettext("Procrustes similarity indexes between partial configurations"))
        simi.check <- tkcheckbutton(SortieWin)
        if(Rsimi) simi.Value <- tclVar("1")
        else simi.Value <- tclVar("0")
        tkconfigure(simi.check,variable=simi.Value)
        
      scaling.lab<-tklabel(SortieWin,text=.Facto_gettext("Isotropic scaling factors"))
        scaling.check <- tkcheckbutton(SortieWin)
        if(Rscaling) scalingValue <- tclVar("1")
        else scalingValue <- tclVar("0")
        tkconfigure(scaling.check,variable=scalingValue)
      
      dep.lab <-tklabel(SortieWin, text=.Facto_gettext("Initial partial configurations"))
        dep.check <- tkcheckbutton(SortieWin)
        if(Rdep) depValue <- tclVar("1")
        else depValue <- tclVar("0")
        tkconfigure(dep.check,variable=depValue)
        
      consensus.lab<-tklabel(SortieWin,text=.Facto_gettext("Consensus configuration"))
        consensus.check <- tkcheckbutton(SortieWin)
        if(Rconsensus) consensusValue <- tclVar("1")
        else consensusValue <- tclVar("0")
        tkconfigure(consensus.check,variable=consensusValue)
      
      Xfin.lab<-tklabel(SortieWin,text=.Facto_gettext("Partial configurations after transformations"))
        Xfin.check <- tkcheckbutton(SortieWin)
        if(RXfin) XfinValue <- tclVar("1")
        else XfinValue <- tclVar("0")
        tkconfigure(Xfin.check,variable=XfinValue)
      
      PANOVA.lab<-tklabel(SortieWin,text=.Facto_gettext("Procrustes Analysis of Variance tables, per configuration, per objet, per dimension"))
        PANOVA.check <- tkcheckbutton(SortieWin)
        if(RPANOVA) PANOVAValue <- tclVar("1")
        else PANOVAValue <- tclVar("0")
        tkconfigure(PANOVA.check,variable=PANOVAValue)
      
      correlations.lab<-tklabel(SortieWin,text=.Facto_gettext("Correlations between initial partial configurations and consensus dimensions"))
        correlations.check <- tkcheckbutton(SortieWin)
      if(Rcorrelations) correlationsValue <- tclVar("1")
        else correlationsValue <- tclVar("0")
        tkconfigure(correlations.check,variable=correlationsValue)
        
      RFichierFrame<-tkframe(SortieWin,borderwidth=2)
      if (is.null(RFichier)) Fichier <- tclVar("")
      else Fichier<-tclVar(RFichier)
      Fichier.entry <-tkentry(RFichierFrame,width="40",textvariable=Fichier)

      SortieOK.but<-tkbutton(SortieWin,text="OK",width=16,command=onOK.sortie)

      tkgrid(tklabel(SortieWin, text = .Facto_gettext("Select output options"), fg ="blue"),  columnspan = 2, sticky = "w")
      tkgrid(tklabel(SortieWin, text = " "))
      tkgrid(consensus.lab,consensus.check,sticky="w")
      tkgrid(RV.lab,RV.check,sticky="w")
      tkgrid(RVs.lab,RVs.check,sticky="w")
      tkgrid(simi.lab,simi.check,sticky="w")
      tkgrid(scaling.lab,scaling.check,sticky="w")
      tkgrid(dep.lab,dep.check,sticky="w")
      tkgrid(Xfin.lab,Xfin.check,sticky="w")      
      tkgrid(correlations.lab,correlations.check,sticky="w")      
      tkgrid(PANOVA.lab,PANOVA.check,sticky="w")      
      tkgrid(tklabel(SortieWin, text = " "))
      tkgrid(tklabel(RFichierFrame,text=.Facto_gettext("Print results on a 'csv' file")),Fichier.entry)
      tkgrid(RFichierFrame)
      tkgrid(SortieOK.but)
   }
    
    SortieFrame<-tkframe(IlluFrame)
    Sortie.but<-tkbutton(SortieFrame, textvariable=.SortieLabel, command=OnSortie, borderwidth=3)
    tkgrid(Sortie.but, sticky="ew")
  })



  #! fonction pour la gestion des options graphiques 
  PLOT.GPA<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.graph<-0
    .PlotLabel<-tclVar(paste(firstLabel, "", sep=" "))
    #déclaration des variables
        
    Rchoix<-TRUE
    RTitle<-NULL
    Rlabel.indMoy<-TRUE
    Rhabillage<-"group"
    Rpartial<-"all"
    RpartialSouris<-FALSE
    Rchrono<-FALSE
    RXlimInd<-NULL
    RYlimInd<-NULL
    
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
        #récupération des dimensions à représenter
            
        # gestion des entrées de la partie graphique des individus
        if(tclvalue(ind.check.value)==1) assign("Rchoix", TRUE, envir=env)
        else assign("Rchoix", FALSE, envir=env)

        if(Rchoix)
        {
          if (tclvalue(Titre)==" ") assign("RTitle", NULL, envir=env)
          assign("RTitle", tclvalue(Titre), envir=env)

          label.tmp.indMoy<-tclvalue(label.indMoy.checkValue)
          if(label.tmp.indMoy==1) assign("Rlabel.indMoy", TRUE, envir=env)
          else assign("Rlabel.indMoy", FALSE, envir=env)
          
          habillage.tmp<-listgraph.nom[as.numeric(tkcurselection(listgraph))+1]
          if(length(habillage.tmp)==0) assign("Rhabillage","none", envir=env)
          else assign("Rhabillage", habillage.tmp, envir=env)

          if(tclvalue(XlimIndMin)=="" | tclvalue(XlimIndMax)=="") assign("RXlimInd", NULL, envir=env)
          else assign("RXlimInd", c(as.numeric(tclvalue(XlimIndMin)), as.numeric(tclvalue(XlimIndMax))), envir=env)
          if(tclvalue(YlimIndMin)=="" | tclvalue(YlimIndMax)=="") assign("RYlimInd", NULL, envir=env)
          else assign("RYlimInd", c(as.numeric(tclvalue(YlimIndMin)), as.numeric(tclvalue(YlimIndMax))), envir=env)
          
          partial.tmp<-listpartial.nom[as.numeric(tkcurselection(listpartial))+1]
          if(length(partial.tmp)==0) assign("Rpartial",NULL, envir=env)
          else assign("Rpartial", partial.tmp, envir=env)
          
          chrono.tmp<-tclvalue(partial.chrono.checkValue)
          if(chrono.tmp=="1") assign("Rchrono", TRUE, envir=env)
          else assign("Rchrono", FALSE, envir=env)
          
          souris.tmp<-tclvalue(partial.souris.checkValue)
          if(souris.tmp=="1") assign("RpartialSouris", TRUE, envir=env)
          else assign("RpartialSouris", FALSE, envir=env)
          
        }
        tkdestroy(PlotWin)
      }
    
      # création l'interface "options graphiques"
    
              
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
        
      RhabillageFrame<-tkframe(PlotIndFrame,borderwidth=2)
      listgraph<-tklistbox(RhabillageFrame,height=3, selectmode="single",exportselection="FALSE",yscrollcommand=function(...) tkset(scrgraph,...))
      scrgraph <-tkscrollbar(RhabillageFrame,repeatinterval=5,command=function(...)tkyview(listgraph,...))
      listgraph.nom<-c("group","ind")
      tkinsert(listgraph,"end",.Facto_gettext("by.group"))
      tkinsert(listgraph,"end",.Facto_gettext("by.individual"))
      if(Rhabillage=="group") tkselection.set(listgraph,0)
      if(Rhabillage=="ind") tkselection.set(listgraph,1)
##      indice<-2      
##      for (i in 1:ncol(donnee))
##      {
##          if(is.factor(donnee[,i]))
##          {
##          tkinsert(listgraph,"end",vars[i])
##          listgraph.nom<-c(listgraph.nom,vars[i])
##          if(Rhabillage==vars[i]) tkselection.set(listgraph, indice)
##          indice<-indice+1
##        }
##      }
      tkgrid(tklabel(RhabillageFrame, text=.Facto_gettext("Select drawing for the individuals")))
      tkgrid(listgraph, scrgraph, sticky = "nw")
      tkgrid.configure(scrgraph, sticky = "wns")
      tkgrid.configure(listgraph, sticky = "ew")
      
                  
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
  
      RpartialFrame<-tkframe(PlotIndFrame,borderwidth=2)
      listpartial<-tklistbox(RpartialFrame,height=7, selectmode="extended",exportselection="FALSE",yscrollcommand=function(...) tkset(scrpartial,...))
      scrpartial<-tkscrollbar(RpartialFrame,repeatinterval=5,command=function(...)tkyview(listpartial,...))
      listpartial.nom<-NULL
      indice<-0
      for (i in 1:nrow(donnee)) {
          tkinsert(listpartial,"end",rows[i])
        listpartial.nom<-c(listpartial.nom,rows[i])
        if(rows[i] %in% Rpartial) tkselection.set(listpartial, indice)
        indice<-indice+1 
      }
      partial.souris.check<-tkcheckbutton(RpartialFrame)
      if (RpartialSouris) partial.souris.checkValue<-tclVar("1")
      else partial.souris.checkValue<-tclVar("0") 
      partial.chrono.check<-tkcheckbutton(RpartialFrame)
      if (Rchrono) partial.chrono.checkValue<-tclVar("1")
      else partial.chrono.checkValue<-tclVar("0")
      tkconfigure(partial.souris.check, variable=partial.souris.checkValue)
      tkconfigure(partial.chrono.check, variable=partial.chrono.checkValue)      
      tkgrid(tklabel(RpartialFrame, text=.Facto_gettext("Select the individuals for which partial points are drawn")))
      tkgrid(listpartial, scrpartial, sticky = "nw")
      tkgrid.configure(scrpartial, sticky = "wns")
      tkgrid.configure(listpartial, sticky = "ew")
      tkgrid(tklabel(RpartialFrame, text=.Facto_gettext("Interactive selection of the individuals")), partial.souris.check)
      tkgrid(tklabel(RpartialFrame, text=.Facto_gettext("Chronologic representation of the partial points")), partial.chrono.check)
  
      #mise en page des différents frames de PlotIndFrame
      tkgrid(RchoixFrame)
      tkgrid(RTitleFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))
      tkgrid(RlabelFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))
      tkgrid(RhabillageFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))
      tkgrid(RpartialFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))      
      tkgrid(RlimFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))
              
      #mise en page de plotWin
      subOKCancelHelp(PlotWin, "plot.GPA")
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
    # récupération des paramètres de la fenêtre principale
    nom.res<-tclvalue(resu.val)
    if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) justDoIt(paste('remove (',nom.res,')'))       #if object res already exists, it's removed
    nbiter<-as.numeric(tclvalue(nbiter.val))
    scaling<-as.logical(as.numeric(tclvalue(scale.bool)))
    tolerance<-as.numeric(tclvalue(tolerance.val))
    Axe<-c(as.numeric(tclvalue(Axe1)), as.numeric(tclvalue(Axe2)))
    
    # gestion du tableau de données pour la GPA
    group<-NULL
    type<-NULL
    name.group<-NULL
    num.group.sup<-NULL
    variables<-NULL
    indice.grpe<-1
      #récupération des groupes quanti actif
    nb.GQA<-length(listQuantiAct.nom)
    if(nb.GQA>=1) {
      name.group<-c(name.group, listQuantiAct.nom)
      for(i in 1:nb.GQA) {
        eval(parse(text=paste("liste.var.GQA<-", listQuantiAct.nom[i], ".var", sep="")))
        type<-c(type,liste.var.GQA[1])
## modif 7 juin
##        variables<-c(variables, liste.var.GQA[-1])
##        group<-c(group, length(liste.var.GQA)-1)
        group<-c(group, length(liste.var.GQA))
        variables<-c(variables, liste.var.GQA)
        indice.grpe<-indice.grpe+1
      }
    }
    
      #construction du tableau de données.GPA
      commande.data<-paste(activeDataSet(),'.GPA', '<-', activeDataSet(),'[ , c("',paste(variables, collapse='", "'), '")]', sep="")
      
      justDoIt(commande.data)
      logger(commande.data)
      donnee.depart<-activeDataSet()
      activeDataSet(paste(activeDataSet(),'.GPA', sep=""))

      # gestion de la commande réalisant la GPA     
      commande.GPA<-paste(nom.res, '<-GPA(', activeDataSet(), ', group=c(',paste(group, collapse=", "), '), name.group=c("',paste(name.group, collapse='", "'), '"), scale=',scaling,', graph=FALSE',sep="")
      if (nbiter!=200) commande.GPA<-paste(commande.GPA, ', nbiteration =',nbiter,sep='')
      if (tolerance!=1e-10) commande.GPA<-paste(commande.GPA, ', tolerance =',tolerance,sep='')
      commande.GPA<-paste(commande.GPA, ')',sep='')
      justDoIt(commande.GPA)
      logger(commande.GPA)
      
      if((Rchoix)&(length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0)){
        if ((Rhabillage!="none") & (Rhabillage!="ind") & (Rhabillage!="group")) {
          Rhabillage<-which(colnames(get(getRcmdr(".activeDataSet")))==Rhabillage)
          if(length(Rhabillage)==0) Rhabillage<-"none"
        }
        if (Rhabillage=="none") Rhabillage<-paste('"', Rhabillage, '"', sep="")
        if (Rhabillage=="ind") Rhabillage<-paste('"', Rhabillage, '"', sep="")
        if (Rhabillage=="group") Rhabillage<-paste('"', Rhabillage, '"', sep="")
        
        if(RpartialSouris){
          commande.plotI<-paste('plotGPApartial(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), lab.ind.moy=', Rlabel.indMoy, ', habillage=', Rhabillage, sep="")
          if (!is.null(RXlimInd)) commande.plotI<-paste(commande.plotI, ', xlim=c(', paste(RXlimInd, collapse=", "), ')')
          if (!is.null(RYlimInd)) commande.plotI<-paste(commande.plotI, ', ylim=c(', paste(RYlimInd, collapse=", "), ')')
          if (Rchrono) commande.plotI<-paste(commande.plotI, ', chrono=', Rchrono, sep='')
          if (is.null(RTitle)) commande.plotI <- paste(commande.plotI,')', sep="")
          else {
            if (RTitle ==" ") commande.plotI <- paste(commande.plotI,')', sep="")
            else commande.plotI <- paste(commande.plotI,', title="', RTitle,'")', sep="")
          }
        }
        else {
          commande.plotI<-paste('plot.GPA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), lab.ind.moy=', Rlabel.indMoy, ', habillage=', Rhabillage, sep="") 
          if (!is.null(RXlimInd)) commande.plotI<-paste(commande.plotI, ', xlim=c(', paste(RXlimInd, collapse=", "), ')')
          if (!is.null(RYlimInd)) commande.plotI<-paste(commande.plotI, ', ylim=c(', paste(RYlimInd, collapse=", "), ')')
          if (!is.null(Rpartial)) commande.plotI<-paste(commande.plotI, ', partial=c("', paste(Rpartial, collapse='", "'),'")', sep='')
          if (Rchrono) commande.plotI<-paste(commande.plotI, ', chrono=', Rchrono, sep='')
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
      if(RRV) doItAndPrint(paste(nom.res, '$RV', sep=""))
      if(RRVs) doItAndPrint(paste( nom.res, '$RVs', sep=""))
      if(Rsimi) doItAndPrint(paste( nom.res, '$simi', sep=""))
      if(Rscaling) doItAndPrint(paste( nom.res, '$scaling', sep=""))
      if(Rdep) doItAndPrint(paste( nom.res, '$dep', sep=""))
      if(Rconsensus) doItAndPrint(paste( nom.res, '$consensus', sep=""))
      if(RXfin) doItAndPrint(paste( nom.res, '$Xfin', sep=""))
      if(Rcorrelations) doItAndPrint(paste( nom.res, '$correlations', sep=""))
      if(RPANOVA) doItAndPrint(paste( nom.res, '$PANOVA', sep=""))
    }
    else {
      Fich = RFichier
      if (substr(Fich,1,1)!='"') Fich = paste('"',Fich,sep='')
      if (substr(Fich,nchar(Fich),nchar(Fich))!='"') Fich = paste(Fich,'"',sep='')
      append = FALSE
      if(RRV){
        doItAndPrint(paste('write.infile(', nom.res, '$RV, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(RRVs){
        doItAndPrint(paste('write.infile(', nom.res, '$RVs, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rsimi){
        doItAndPrint(paste('write.infile(', nom.res, '$simi, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rscaling){
        doItAndPrint(paste('write.infile(', nom.res, '$scaling, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rdep){
        doItAndPrint(paste('write.infile(', nom.res, '$dep, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rconsensus){
        doItAndPrint(paste('write.infile(', nom.res, '$consensus, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(RXfin){
        doItAndPrint(paste('write.infile(', nom.res, '$Xfin, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rcorrelations){
        doItAndPrint(paste('write.infile(', nom.res, '$correlations, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(RPANOVA) doItAndPrint(paste('write.infile(', nom.res, '$PANOVA, file =',Fich,',append=',append,')', sep=""))
    }

      # Re-chargement du tableau de départ et supression du tableau temporaire
      activeDataSet(donnee.depart)
      justDoIt(paste('remove(',activeDataSet(),'.GPA)',sep=""))
      logger(paste('remove(',activeDataSet(),'.GPA)',sep=""))   
  }


    #! fonction associée au bouton OK, execute et détruit l'interface graphique
  onOK<-function()
  {
    OnAppliquer()
    tkdestroy(top)     
  }



#                   Création de la fenêtre top                                 #
################################################################################
  top<-tktoplevel(borderwidth=10)
  tkwm.title(top,.Facto_gettext("GPA"))
  tkwm.geometry(top, "-50+50")
  
  # définition des polices
  font2<-tkfont.create(family="times",size=12,weight="bold")
  fontheading<-tkfont.create(family="times",size=11,weight="bold")

  # récupération du jeu de données actif
  donnee<-get(getRcmdr(".activeDataSet"))
  vars<-colnames(donnee)
  rows<-rownames(donnee)
  
  # création du frame contenant les listes groupes quanti
  ListeQuantiFrame<- tkframe(top, borderwidth=2, relief="groove")
  label.quantiFrame.var<-tclVar(.Facto_gettext("group"))
  label.quantiFrame<-tklabel(ListeQuantiFrame, textvariable=label.quantiFrame.var,fg = "darkred", font=fontheading)
  # liste des groupes de variables quanti Actives
  listQuantiAct<-tklistbox(ListeQuantiFrame,selectmode="extended",exportselection="FALSE", height=4, yscrollcommand=function(...)tkset(scrQuantiAct,...))
  scrQuantiAct<-tkscrollbar(ListeQuantiFrame,repeatinterval=5,command=function(...)tkyview(listQuantiAct,...))
  listQuantiAct.nom<-NULL
  
    # boutons d'action groupes quantitative
  supprimeQuanti.funct(label=.Facto_gettext("Delete"))
  ajoutQuanti.funct(label=.Facto_gettext("Add 1 group"), firstLabel=.Facto_gettext("Add 1 group"))
  modifQuanti.funct(label=.Facto_gettext("Modify 1 group"), firstLabel=.Facto_gettext("Modify 1 group"))
    # mise en forme de ListeQuantiFrame
  tkgrid(label.quantiFrame, columnspan=11, sticky = "ew")
  tkgrid(listQuantiAct, scrQuantiAct)
  tkgrid.configure(scrQuantiAct, column=3, sticky="wns")
  tkgrid.configure(listQuantiAct, sticky = "ew", column=4, columnspan=2)
  tkgrid.configure(tklabel(ListeQuantiFrame, text=" "))
  tkgrid.configure(GpeQuantiFrame,ModifGpeQuantiFrame, SupGpeQuantiFrame)
  tkgrid.configure(GpeQuantiFrame, sticky = "ew", column=1, columnspan=2)
  tkgrid.configure(ModifGpeQuantiFrame, sticky = "ew", column=4, columnspan=2)
  tkgrid.configure(SupGpeQuantiFrame, sticky = "ew", column=7, columnspan=2)
  tkgrid.columnconfigure(ListeQuantiFrame,0, minsize=25)
  tkgrid.columnconfigure(ListeQuantiFrame,10, minsize=25)
  tkgrid.columnconfigure(ListeQuantiFrame,3, minsize=25)
  tkgrid.columnconfigure(ListeQuantiFrame,9, minsize=25)
  tkgrid.columnconfigure(ListeQuantiFrame,4, minsize=35)
  tkgrid.columnconfigure(ListeQuantiFrame,5, minsize=35)
 
    
   # création de tous les boutons d'options dans IlluFrame
  IlluFrame<- tkframe(top, borderwidth=2)
    # mise en page de IlluFrame
  PLOT.GPA(label=.Facto_gettext("Graphical options"), firstLabel=.Facto_gettext("Graphical options"))
  Sortie.funct(label=.Facto_gettext("Outputs"), firstLabel=.Facto_gettext("Outputs"))
  tkgrid(PlotFrame, SortieFrame, columnspan=7)
  tkgrid.configure(PlotFrame, column=1, columnspan=1)
  tkgrid.configure(SortieFrame, column=3, columnspan=1)
  tkgrid.columnconfigure(IlluFrame,0, minsize=25)
  tkgrid.columnconfigure(IlluFrame,2, minsize=40)  
  tkgrid.columnconfigure(IlluFrame,4, minsize=25)  


    # création des options dans OptionFrame  
  OptionFrame<-tkframe(top, borderwidth=2, relief="groove")
  resu.lab<-tklabel(OptionFrame,text=.Facto_gettext("Name of the result object: "))
  resu.val<-tclVar("res")
  resu<-tkentry(OptionFrame,width=10,textvariable=resu.val)
    scale.check <- tkcheckbutton(top)
    scale.bool <- tclVar("1")
    tkconfigure(scale.check,variable=scale.bool)
    scale.lab<-tklabel(OptionFrame,text=.Facto_gettext("Scale the variables:"))
  nbiter.lab<-tklabel(OptionFrame,text=.Facto_gettext("Maximum number of iteration for the algorithm:"))
  nbiter.val<-tclVar("200") 
  nbiter<-tkentry(OptionFrame,width=5,textvariable=nbiter.val)
  tolerance.lab<-tklabel(OptionFrame,text=.Facto_gettext("Stopping threshold for the algorithm:"))
  tolerance.val<-tclVar("1e-10") 
  tolerance<-tkentry(OptionFrame,width=5,textvariable=tolerance.val)
  Axe.label<-tklabel(OptionFrame,text=.Facto_gettext("Select the dimensions for the graphs:"))
  Axe1<-tclVar("1")
  Axe2<-tclVar("2")
  Axe1.entry <-tkentry(OptionFrame,width="5",textvariable=Axe1)
  Axe2.entry <-tkentry(OptionFrame,width="5",textvariable=Axe2)
  
    # mise en page de OptionFrame
  tkgrid(tklabel(OptionFrame,text=.Facto_gettext("Main options"), fg = "darkred"), columnspan=8, sticky="we") 
  tkgrid(tklabel(OptionFrame,text="")) 
  tkgrid(scale.lab,scale.check,sticky="w")
  tkgrid(Axe.label,Axe1.entry , Axe2.entry, sticky="w")
  tkgrid(nbiter.lab, nbiter)
  tkgrid(tolerance.lab, tolerance)
  tkgrid(resu.lab, resu)
  tkgrid.configure(scale.lab, nbiter.lab, tolerance.lab, resu.lab, Axe.label, column=1, columnspan=4, sticky="w")
  tkgrid.configure(scale.check, tolerance, nbiter, resu, column=6, columnspan=2, sticky="e")
  tkgrid.configure(Axe1.entry, column=6, columnspan=1, sticky="w")
  tkgrid.configure(Axe2.entry, column=7, columnspan=1, sticky="e")
  tkgrid.columnconfigure(OptionFrame,0, minsize=25)
  tkgrid.columnconfigure(OptionFrame,5, minsize=40)
  tkgrid.columnconfigure(OptionFrame,8, minsize=25)

  appliquer.but<-tkbutton(top, text=.Facto_gettext("Apply"),width=12,command=OnAppliquer, borderwidth=3, fg="#690f96")
  OKCancelHelp(helpSubject="GPA",reset="Reinitializ.funct", apply ="OnAppliquer")

  # Mise en page de top
  tkgrid(tklabel(top, text=.Facto_gettext("General Procrustes Analysis (GPA)"),font=fontheading), columnspan=3)
  tkgrid(tklabel(top,text=""))
  tkgrid(ListeQuantiFrame, column=1, columnspan=1, sticky="ew")
  tkgrid(tklabel(top,text=""))    
  tkgrid(IlluFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text=""))    
  tkgrid(OptionFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text="")) # Ligne de blanc       
  # tkgrid(appliquer.but, column=1, columnspan=1)
  # tkgrid(tklabel(top,text="")) # Ligne de blanc  
  # tkgrid(buttonsFrame, column=1, columnspan=1, sticky="ew" )
  tkgrid(buttonsFrame, appliquer.but)
  tkgrid.configure(buttonsFrame, column=1,sticky="e")
  tkgrid.configure(appliquer.but, column=2,sticky="w")



 
}
