FactoMFA <-
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
      grpeIlluSupprime<-listQuantiIllu.nom[as.numeric(tkcurselection(listQuantiIllu))+1]
      if(length(grpeIlluSupprime)>=1)
      { 
        listQuantiIllu.nom.tmp<-listQuantiIllu.nom[-which(listQuantiIllu.nom %in% grpeIlluSupprime)]
        assign("listQuantiIllu.nom",listQuantiIllu.nom.tmp, envir=env)
        tkdelete(listQuantiIllu,"0","end")
        if(length(listQuantiIllu.nom)>=1) {
          for (grpe in listQuantiIllu.nom) tkinsert(listQuantiIllu, "end", grpe)
        }  
      }
      nb.grpe<-length(listQuantiAct.nom) + length(listQuantiIllu.nom)
      if (nb.grpe>=1) {
        #tclvalue(.AjoutQuantiLabel)<-paste(nb.grpe, "Add quanti group(s)", sep=" ")
        #tkconfigure(AjoutGpeQuanti.but, fg="blue")
        tclvalue(label.quantiFrame.var)<-paste(nb.grpe, .Facto_gettext("quantitative groups"), sep=" ")
        tkconfigure(label.quantiFrame)
      }  
      else {
        #tclvalue(.AjoutQuantiLabel)<-paste("Add quanti group", "!", sep=" ")
        #tkconfigure(AjoutGpeQuanti.but, fg="black")
        tclvalue(label.quantiFrame.var)<-paste("0", .Facto_gettext("quantitative group"), sep=" ")
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
    .AjoutQuantiLabel<-tclVar(paste(firstLabel, "!", sep=" "))
    OnAGQ<-function()
    {
      AjoutGpeQuantiWin<-tktoplevel()
      tkwm.title(AjoutGpeQuantiWin,.Facto_gettext("Definition of a quantitatif group"))
      
      #création de la fonction AGA.OK
      AGQ.OK<-function()
      {
        assign("compteur.GQ", compteur.GQ+1, envir=env)
        nom.groupe<-nom.correct(tclvalue(nomGrpeQuanti.val), liste=c(listQuantiAct.nom,listQuantiIllu.nom))
        if (nom.groupe=="") tkmessageBox(message=.Facto_gettext("Name for the group"), icon="warning", type="ok") 
        else {
          varGroupe<-listVarQuanti.nom[as.numeric(tkcurselection(listVarQuanti))+1]
          if (length(varGroupe)>=1) {
            if(tclvalue(norm.Value)=="ok") assign(paste(nom.groupe,".var", sep=""), c("s", varGroupe), envir=env)
            if(tclvalue(norm.Value)=="nok") assign(paste(nom.groupe,".var", sep=""),c("c", varGroupe), envir=env)
            if(tclvalue(etat.Value)=="actif") {  
              tkinsert(listQuantiAct,"end",nom.groupe)
              assign("listQuantiAct.nom", c(listQuantiAct.nom, nom.groupe),envir=env)
            }  
            if(tclvalue(etat.Value)=="illu") {
              tkinsert(listQuantiIllu,"end",nom.groupe)
              assign("listQuantiIllu.nom",c(listQuantiIllu.nom, nom.groupe),envir=env)
              
            }
            #tclvalue(.AjoutQuantiLabel)<-paste(length(listQuantiAct.nom) + length(listQuantiIllu.nom), label, sep=" ")
            #tkconfigure(AjoutGpeQuanti.but, fg="blue")
            tclvalue(label.quantiFrame.var)<-paste(length(listQuantiAct.nom) + length(listQuantiIllu.nom), .Facto_gettext("quantitative groups"), sep=" ")
            tkconfigure(label.quantiFrame)  
            tkdestroy(AjoutGpeQuantiWin)
          }
        }  
      }
      
      
      # choix du nom du groupe
      nomGrpeQuanti.lab<-tklabel(AjoutGpeQuantiWin,text=.Facto_gettext("Name of the group:"))
      nomGrpeQuanti.val<-tclVar(paste("Gc", compteur.GQ, sep=""))
      nomGrpeQuanti<-tkentry(AjoutGpeQuantiWin,width=15,textvariable=nomGrpeQuanti.val)
      
      # choix de l'état actif ou illustratif
      etat.actif.check<-tkradiobutton(AjoutGpeQuantiWin)
      etat.illu.check<-tkradiobutton(AjoutGpeQuantiWin)
      etat.Value<-tclVar("actif")
      tkconfigure(etat.actif.check,variable=etat.Value,value="actif")
      tkconfigure(etat.illu.check,variable=etat.Value, value="illu")
      
      # choix de la normalisation ou non
      norm.ok.check<-tkradiobutton(AjoutGpeQuantiWin)
      norm.nok.check<-tkradiobutton(AjoutGpeQuantiWin)
      norm.Value<-tclVar("ok")
      tkconfigure(norm.ok.check,variable=norm.Value,value="ok")
      tkconfigure(norm.nok.check,variable=norm.Value, value="nok")
            
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
      tkgrid(tklabel(AjoutGpeQuantiWin, text=.Facto_gettext("Status of the group:")), tklabel(AjoutGpeQuantiWin, text="Active"),etat.actif.check, tklabel(AjoutGpeQuantiWin, text=.Facto_gettext("Supplementary")),etat.illu.check, sticky="w")
      tkgrid(tklabel(AjoutGpeQuantiWin, text=""))      
      tkgrid(tklabel(AjoutGpeQuantiWin, text=.Facto_gettext("Scale the variable of the group:")), tklabel(AjoutGpeQuantiWin, text=.Facto_gettext("Yes")),norm.ok.check, tklabel(AjoutGpeQuantiWin, text=.Facto_gettext("No")),norm.nok.check, sticky="w")
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
    OnMGQ<-function()
    {
      ModifGpeQuantiWin<-tktoplevel()
      tkwm.title(ModifGpeQuantiWin,.Facto_gettext("Modification of a quantitative group"))
      
      #création de la fonction AGA.OK
      MGQ.OK<-function()
      {
        nom.groupe<-nom.correct(tclvalue(nomModifGrpeQuanti.val), liste=c(listQuantiAct.nom,listQuantiIllu.nom))
        if (nom.groupe=="") tkmessageBox(message=.Facto_gettext("Name for the group"), icon="warning", type="ok") 
        else {
          if(etat=="actif") {
            listQuantiAct.nom.tmp<-listQuantiAct.nom[-which(listQuantiAct.nom== grpeAModifier)]
            assign("listQuantiAct.nom",listQuantiAct.nom.tmp, envir=env)
            tkdelete(listQuantiAct,"0","end")
            for (grpe in listQuantiAct.nom) tkinsert(listQuantiAct, "end", grpe)
          }
          
          if(etat=="illu")  {
            listQuantiIllu.nom.tmp<-listQuantiIllu.nom[-which(listQuantiIllu.nom== grpeAModifier)]
            assign("listQuantiIllu.nom",listQuantiIllu.nom.tmp, envir=env)
            tkdelete(listQuantiIllu,"0","end")
            for (grpe in listQuantiIllu.nom) tkinsert(listQuantiIllu, "end", grpe)
          }
          
          varGroupe<-listModifVarQuanti.nom[as.numeric(tkcurselection(listModifVarQuanti))+1]
          if (length(varGroupe)>=1) {
            if(tclvalue(normModif.Value)=="ok") assign(paste(nom.groupe,".var", sep=""), c("s", varGroupe), envir=env)
            if(tclvalue(normModif.Value)=="nok") assign(paste(nom.groupe,".var", sep=""),c("c", varGroupe), envir=env)
            if(tclvalue(etatModif.Value)=="actif") {  
              tkinsert(listQuantiAct,"end",nom.groupe)
              assign("listQuantiAct.nom", c(listQuantiAct.nom, nom.groupe),envir=env)
            }  
            if(tclvalue(etatModif.Value)=="illu") {
              tkinsert(listQuantiIllu,"end",nom.groupe)
              assign("listQuantiIllu.nom",c(listQuantiIllu.nom, nom.groupe),envir=env)              
            }
              
            tkdestroy(ModifGpeQuantiWin)
          }
        }  
      }
            
      if(length(as.numeric(tkcurselection(listQuantiAct)))>=1) {
        grpeAModifier<-listQuantiAct.nom[as.numeric(tkcurselection(listQuantiAct))+1][1]
        etat<-"actif"
      }
      else if (length(as.numeric(tkcurselection(listQuantiIllu)))>=1) {
        grpeAModifier<-listQuantiIllu.nom[as.numeric(tkcurselection(listQuantiIllu))+1][1]
        etat<-"illu"
      }
      else  {
        tkdestroy(ModifGpeQuantiWin)
        return()
      }  
        
      eval(parse(text=paste("grpeAModifier.var<-",paste(grpeAModifier,".var", sep=""),sep="")))
      # choix du nom du groupe
      nomModifGrpeQuanti.lab<-tklabel(ModifGpeQuantiWin,text=.Facto_gettext("Name of the group:"))
      nomModifGrpeQuanti.val<-tclVar(grpeAModifier)
      nomModifGrpeQuanti<-tkentry(ModifGpeQuantiWin,width=15,textvariable=nomModifGrpeQuanti.val)
      
      # choix de l'état actif ou illustratif
      etatModif.actif.check<-tkradiobutton(ModifGpeQuantiWin)
      etatModif.illu.check<-tkradiobutton(ModifGpeQuantiWin)
      etatModif.Value<-tclVar(etat)
      tkconfigure(etatModif.actif.check,variable=etatModif.Value,value="actif")
        tkconfigure(etatModif.illu.check,variable=etatModif.Value, value="illu")
      
      # choix de la normalisation ou non
      normModif.ok.check<-tkradiobutton(ModifGpeQuantiWin)
      normModif.nok.check<-tkradiobutton(ModifGpeQuantiWin)
      if(grpeAModifier.var[1]=="s") normModif.Value<-tclVar("ok")
      else normModif.Value<-tclVar("nok")
      tkconfigure(normModif.ok.check,variable=normModif.Value,value="ok")
        tkconfigure(normModif.nok.check,variable=normModif.Value, value="nok")
            
      # création de la liste pour le choix des variables acives
      listModifVarQuanti<-tklistbox(ModifGpeQuantiWin, selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrModifVarQuanti,...))
      scrModifVarQuanti<-tkscrollbar(ModifGpeQuantiWin,repeatinterval=5,command=function(...)tkyview(listModifVarQuanti,...))
      listModifVarQuanti.nom<-NULL
      indice.num<-0
      for (i in (1:ncol(donnee))) {
        if (is.numeric(donnee[,i])) {
            tkinsert(listModifVarQuanti,"end",vars[i])
            listModifVarQuanti.nom<-c(listModifVarQuanti.nom, vars[i])
            if(vars[i] %in% grpeAModifier.var[-1]) tkselection.set(listModifVarQuanti, indice.num)
            indice.num<-indice.num+1
        }
      }
      
      MGQ.but<-tkbutton(ModifGpeQuantiWin, text="OK", width=16, command=MGQ.OK)
      
      tkgrid(nomModifGrpeQuanti.lab, nomModifGrpeQuanti)
      tkgrid.configure(nomModifGrpeQuanti.lab, column=0, columnspan=2, sticky="w")
      tkgrid.configure(nomModifGrpeQuanti, column=2, columnspan=3)      
      tkgrid(tklabel(ModifGpeQuantiWin, text=""))      
      tkgrid(tklabel(ModifGpeQuantiWin, text=.Facto_gettext("Status of the group:")), tklabel(ModifGpeQuantiWin, text=.Facto_gettext("Active")),etatModif.actif.check, tklabel(ModifGpeQuantiWin, text=.Facto_gettext("Supplementary")),etatModif.illu.check, sticky="w")
      tkgrid(tklabel(ModifGpeQuantiWin, text=""))      
      tkgrid(tklabel(ModifGpeQuantiWin, text=.Facto_gettext("Scale the variables of the group:")), tklabel(ModifGpeQuantiWin, text="Yes"),normModif.ok.check, tklabel(ModifGpeQuantiWin, text=.Facto_gettext("No")),normModif.nok.check, sticky="w")
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
  



  #! Suppression de groupes qualitatifs
  supprimeQuali.funct<-defmacro(label, expr=
  {
    env<-environment()
    OnSGQl<-function()
    {
      grpeActSupprime<-listQualiAct.nom[as.numeric(tkcurselection(listQualiAct))+1]
      if(length(grpeActSupprime)>=1) { 
        listQualiAct.nom.tmp<-listQualiAct.nom[-which(listQualiAct.nom %in% grpeActSupprime)]
        assign("listQualiAct.nom",listQualiAct.nom.tmp, envir=env)
        tkdelete(listQualiAct,"0","end")
        if(length(listQualiAct.nom)>=1) {
          for (grpe in listQualiAct.nom) tkinsert(listQualiAct, "end", grpe)
        }  
      }
      grpeIlluSupprime<-listQualiIllu.nom[as.numeric(tkcurselection(listQualiIllu))+1]
      if(length(grpeIlluSupprime)>=1) { 
        listQualiIllu.nom.tmp<-listQualiIllu.nom[-which(listQualiIllu.nom %in% grpeIlluSupprime)]
        assign("listQualiIllu.nom",listQualiIllu.nom.tmp, envir=env)
        tkdelete(listQualiIllu,"0","end")
        if(length(listQualiIllu.nom)>=1) {
          for (grpe in listQualiIllu.nom) tkinsert(listQualiIllu, "end", grpe)
        }  
      }
      nb.grpe<-length(listQualiAct.nom) + length(listQualiIllu.nom)
      if (nb.grpe>=1) {
        tclvalue(label.qualiFrame.var)<-paste(nb.grpe, .Facto_gettext("qualitative groups"), sep=" ")
        tkconfigure(label.qualiFrame)
      }  
      else {
        tclvalue(label.qualiFrame.var)<-paste("0", .Facto_gettext("qualitative group"), sep=" ")
        tkconfigure(label.qualiFrame)
      }  
    }
    
    SupGpeQualiFrame<-tkframe(ListeQualiFrame)
    SupGpeQuali.but<-tkbutton(SupGpeQualiFrame, textvariable=tclVar(label), command=OnSGQl, borderwidth=3)
   tkgrid(SupGpeQuali.but, sticky="ew")  
  })

  
  
  
  
  #! Ajout d'un groupe qualitatif 
  ajoutQuali.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.GQl<-1
    nbGrpeQualiAct<-0
    nbGrpeQualiIllu<-0
    .AjoutQualiLabel<-tclVar(paste(firstLabel, "!", sep=" "))
    OnAGQl<-function()
    {
      AjoutGpeQualiWin<-tktoplevel()
      tkwm.title(AjoutGpeQualiWin,.Facto_gettext("Construction of a qualitative group"))
      
      #création de la fonction AGA.OK
      AGQl.OK<-function()
      {
        assign("compteur.GQl", compteur.GQl+1, envir=env)
        nom.groupe<-nom.correct(tclvalue(nomGrpeQuali.val), liste=c(listQualiAct.nom,listQualiIllu.nom))
        if (nom.groupe=="") tkmessageBox(message=.Facto_gettext("Name for the group"), icon="warning", type="ok") 
        else
        {
          varGroupe<-listVarQuali.nom[as.numeric(tkcurselection(listVarQuali))+1]
          if (length(varGroupe)>=1) {
            assign(paste(nom.groupe,".var", sep=""), c("n",varGroupe), envir=env)
            if(tclvalue(etat.Value)=="actif") {  
              tkinsert(listQualiAct,"end",nom.groupe)
              assign("listQualiAct.nom", c(listQualiAct.nom, nom.groupe),envir=env)
            }  
            if(tclvalue(etat.Value)=="illu") {
              tkinsert(listQualiIllu,"end",nom.groupe)
              assign("listQualiIllu.nom",c(listQualiIllu.nom, nom.groupe),envir=env)
            }  
            tclvalue(label.qualiFrame.var)<-paste(length(listQualiAct.nom) + length(listQualiIllu.nom), "qualitative groups", sep=" ")
            tkconfigure(label.qualiFrame) 
            tkdestroy(AjoutGpeQualiWin)
          }
        }  
      }
      
      
      # choix du nom du groupe
      nomGrpeQuali.lab<-tklabel(AjoutGpeQualiWin,text=.Facto_gettext("Name of the group:"))
      nomGrpeQuali.val<-tclVar(paste("Gq", compteur.GQl, sep=""))
      nomGrpeQuali<-tkentry(AjoutGpeQualiWin,width=15,textvariable=nomGrpeQuali.val)
      
      # choix de l'état actif ou illustratif
      etat.actif.check<-tkradiobutton(AjoutGpeQualiWin)
      etat.illu.check<-tkradiobutton(AjoutGpeQualiWin)
      etat.Value<-tclVar("actif")
      tkconfigure(etat.actif.check,variable=etat.Value,value="actif")
      tkconfigure(etat.illu.check,variable=etat.Value, value="illu")
      
      # création de la liste pour le choix des variables acives
      listVarQuali<-tklistbox(AjoutGpeQualiWin, selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrVarQuali,...))
      scrVarQuali<-tkscrollbar(AjoutGpeQualiWin,repeatinterval=5,command=function(...)tkyview(listVarQuali,...))
      listVarQuali.nom<-NULL
      for (i in (1:ncol(donnee))) {
        if (is.factor(donnee[,i])) {
          tkinsert(listVarQuali,"end",vars[i])
          listVarQuali.nom<-c(listVarQuali.nom, vars[i])
        }
      }
      
      AGQl.but<-tkbutton(AjoutGpeQualiWin, text="OK", width=16, command=AGQl.OK)
      
      tkgrid(nomGrpeQuali.lab, nomGrpeQuali)
      tkgrid.configure(nomGrpeQuali.lab, column=0, columnspan=2, sticky="w")
      tkgrid.configure(nomGrpeQuali, column=2, columnspan=3)      
      tkgrid(tklabel(AjoutGpeQualiWin, text=""))      
      tkgrid(tklabel(AjoutGpeQualiWin, text=.Facto_gettext("Status of the group:")), tklabel(AjoutGpeQualiWin, text=.Facto_gettext("Active")),etat.actif.check, tklabel(AjoutGpeQualiWin, text=.Facto_gettext("Supplementary")),etat.illu.check, sticky="w")
      tkgrid(tklabel(AjoutGpeQualiWin, text=""))
      tkgrid(tklabel(AjoutGpeQualiWin, text = .Facto_gettext("Select the variables of the group"), fg = "blue"), column=0, columnspan = 5, sticky = "w")
      tkgrid(listVarQuali, scrVarQuali, sticky = "nw")
      tkgrid.configure(scrVarQuali, sticky = "wns", column=4,columnspan=1)
      tkgrid.configure(listVarQuali, sticky = "ew", column=0, columnspan=4)
      tkgrid(tklabel(AjoutGpeQualiWin, text=""))
      tkgrid(AGQl.but, column=2,columnspan=1, sticky="ew")
      tkgrid.columnconfigure(AjoutGpeQualiWin,0, minsize=55)
      tkgrid.columnconfigure(AjoutGpeQualiWin,1, minsize=55)
      tkgrid.columnconfigure(AjoutGpeQualiWin,2, minsize=55)
      tkgrid.columnconfigure(AjoutGpeQualiWin,3, minsize=55)
      tkgrid.columnconfigure(AjoutGpeQualiWin,4, minsize=55)              
   }
   GpeQualiFrame<-tkframe(ListeQualiFrame)
   AjoutGpeQuali.but<-tkbutton(GpeQualiFrame, textvariable=.AjoutQualiLabel, command=OnAGQl, borderwidth=3)
   tkgrid(AjoutGpeQuali.but, sticky="ew")
  })    


    #! Modification d'un groupe qualitatif
  modifQuali.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    .ModifQualiLabel<-tclVar(paste(firstLabel, "", sep=" "))
    OnMGQl<-function()
    {
      ModifGpeQualiWin<-tktoplevel()
      tkwm.title(ModifGpeQualiWin,.Facto_gettext("Modification of a qualitative group"))
      
      #création de la fonction AGA.OK
      MGQl.OK<-function()
      {
        nom.groupe<-nom.correct(tclvalue(nomModifGrpeQuali.val), liste=c(listQualiAct.nom,listQualiIllu.nom))
        if (nom.groupe=="") tkmessageBox(message=.Facto_gettext("Give a name for the group"), icon="warning", type="ok") 
        else {
          if(etat=="actif") {
            listQualiAct.nom.tmp<-listQualiAct.nom[-which(listQualiAct.nom== grpeAModifier)]
            assign("listQualiAct.nom",listQualiAct.nom.tmp, envir=env)
            tkdelete(listQualiAct,"0","end")
            for (grpe in listQualiAct.nom) tkinsert(listQualiAct, "end", grpe)
          }          
          if(etat=="illu") {
            listQualiIllu.nom.tmp<-listQualiIllu.nom[-which(listQualiIllu.nom== grpeAModifier)]
            assign("listQualiIllu.nom",listQualiIllu.nom.tmp, envir=env)
            tkdelete(listQualiIllu,"0","end")
            for (grpe in listQualiIllu.nom) tkinsert(listQualiIllu, "end", grpe)
          }          
          varGroupe<-listModifVarQuali.nom[as.numeric(tkcurselection(listModifVarQuali))+1]
          if (length(varGroupe)>=1) {
            assign(paste(nom.groupe,".var", sep=""), c("n", varGroupe), envir=env)
            if(tclvalue(etatModif.Value)=="actif") {  
              tkinsert(listQualiAct,"end",nom.groupe)
              assign("listQualiAct.nom", c(listQualiAct.nom, nom.groupe),envir=env)
            }  
            if(tclvalue(etatModif.Value)=="illu") {
              tkinsert(listQualiIllu,"end",nom.groupe)
              assign("listQualiIllu.nom",c(listQualiIllu.nom, nom.groupe),envir=env)              
            }
              
            tkdestroy(ModifGpeQualiWin)
          }
        }  
      }
      
      
      if(length(as.numeric(tkcurselection(listQualiAct)))>=1) {
        grpeAModifier<-listQualiAct.nom[as.numeric(tkcurselection(listQualiAct))+1][1]
        etat<-"actif"
      }
      else if (length(as.numeric(tkcurselection(listQualiIllu)))>=1) {
        grpeAModifier<-listQualiIllu.nom[as.numeric(tkcurselection(listQualiIllu))+1][1]
        etat<-"illu"
      }
      else  {
        tkdestroy(ModifGpeQualiWin)
        return()
      }  
        
       eval(parse(text=paste("grpeAModifier.var<-",paste(grpeAModifier,".var", sep=""),sep="")))
      # choix du nom du groupe
      nomModifGrpeQuali.lab<-tklabel(ModifGpeQualiWin,text=.Facto_gettext("Name of the group:"))
      nomModifGrpeQuali.val<-tclVar(grpeAModifier)
      nomModifGrpeQuali<-tkentry(ModifGpeQualiWin,width=15,textvariable=nomModifGrpeQuali.val)
      
      # choix de l'état actif ou illustratif
      etatModif.actif.check<-tkradiobutton(ModifGpeQualiWin)
      etatModif.illu.check<-tkradiobutton(ModifGpeQualiWin)
      etatModif.Value<-tclVar(etat)
      tkconfigure(etatModif.actif.check,variable=etatModif.Value,value="actif")
      tkconfigure(etatModif.illu.check,variable=etatModif.Value, value="illu")
      
      # création de la liste pour le choix des variables actives
      listModifVarQuali<-tklistbox(ModifGpeQualiWin, selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrModifVarQuali,...))
      scrModifVarQuali<-tkscrollbar(ModifGpeQualiWin,repeatinterval=5,command=function(...)tkyview(listModifVarQuali,...))
      listModifVarQuali.nom<-NULL
      indice.num<-0
      for (i in (1:ncol(donnee))) {
        if (is.factor(donnee[,i])) {
          tkinsert(listModifVarQuali,"end",vars[i])
          listModifVarQuali.nom<-c(listModifVarQuali.nom, vars[i])
          if(vars[i] %in% grpeAModifier.var[-1]) tkselection.set(listModifVarQuali, indice.num)
          indice.num<-indice.num+1
        }
      }
      
      MGQl.but<-tkbutton(ModifGpeQualiWin, text="OK", width=16, command=MGQl.OK)
      
      tkgrid(nomModifGrpeQuali.lab, nomModifGrpeQuali)
      tkgrid.configure(nomModifGrpeQuali.lab, column=0, columnspan=2, sticky="w")
      tkgrid.configure(nomModifGrpeQuali, column=2, columnspan=3)      
      tkgrid(tklabel(ModifGpeQualiWin, text=""))      
      tkgrid(tklabel(ModifGpeQualiWin, text=.Facto_gettext("Status of the group:")), tklabel(ModifGpeQualiWin, text=.Facto_gettext("Active")),etatModif.actif.check, tklabel(ModifGpeQualiWin, text=.Facto_gettext("Supplementary")),etatModif.illu.check, sticky="w")
      tkgrid(tklabel(ModifGpeQualiWin, text=""))      
      tkgrid(tklabel(ModifGpeQualiWin, text = .Facto_gettext("Select the variables of the group"), fg = "blue"), column=0, columnspan = 5, sticky = "w")
      tkgrid(listModifVarQuali, scrModifVarQuali, sticky = "nw")
      tkgrid.configure(scrModifVarQuali, sticky = "wns", column=4,columnspan=1)
      tkgrid.configure(listModifVarQuali, sticky = "ew", column=0, columnspan=4)
      tkgrid(tklabel(ModifGpeQualiWin, text=""))
      tkgrid(MGQl.but, column=2,columnspan=1, sticky="ew")
      tkgrid.columnconfigure(ModifGpeQualiWin,0, minsize=55)
      tkgrid.columnconfigure(ModifGpeQualiWin,1, minsize=55)
      tkgrid.columnconfigure(ModifGpeQualiWin,2, minsize=55)
      tkgrid.columnconfigure(ModifGpeQualiWin,3, minsize=55)
      tkgrid.columnconfigure(ModifGpeQualiWin,4, minsize=55)              
   }
   ModifGpeQualiFrame<-tkframe(ListeQualiFrame)
   ModifGpeQuali.but<-tkbutton(ModifGpeQualiFrame, textvariable=.ModifQualiLabel, command=OnMGQl, borderwidth=3)
   tkgrid(ModifGpeQuali.but, sticky="ew")
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
    FactoMFA()
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
    Rinertie<-FALSE
    Rindividu<-FALSE
    Rindsup<-FALSE
    Rquantisummary<-FALSE
    Rquanti<-FALSE
    Rquantisup<-FALSE    
    Rqualisummary<-FALSE
    Rquali<-FALSE
    Rqualisup<-FALSE
    Raxepartiel<-FALSE          
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
#        if(compteur.sortie>0) tclvalue(.SortieLabel)<-paste(label, ": Seen", sep=" ")
        if(compteur.sortie>0) tclvalue(.SortieLabel)<-paste(label, "", sep=" ")
        tkconfigure(Sortie.but, fg="blue")
        
        if(tclvalue(eigValue)=="1") assign("Rpropre", TRUE, envir=env)
        else assign("Rpropre", FALSE, envir=env)
        
        if(tclvalue(groupeValue)=="1") assign("Rgroupe", TRUE, envir=env)
        else assign("Rgroupe", FALSE, envir=env)
        
        if(tclvalue(inertieValue)=="1") assign("Rinertie", TRUE, envir=env)
        else assign("Rinertie", FALSE, envir=env)
        
        if(tclvalue(indValue)=="1") assign("Rindividu", TRUE, envir=env)
        else assign("Rindividu", FALSE, envir=env)
        
        if(tclvalue(ind.sup.Value)=="1") assign("Rindsup", TRUE, envir=env)
        else assign("Rindsup", FALSE, envir=env)
        
        if(tclvalue(quantiSummaryValue)=="1") assign("Rquantisummary", TRUE, envir=env)
        else assign("Rquantisummary", FALSE, envir=env)
        
        if(tclvalue(quantiValue)=="1") assign("Rquanti", TRUE, envir=env)
        else assign("Rquanti", FALSE, envir=env)
        
        if(tclvalue(quantisupValue)=="1") assign("Rquantisup", TRUE, envir=env)
        else assign("Rquantisup", FALSE, envir=env)
        
        if(tclvalue(qualiSummaryValue)=="1") assign("Rqualisummary", TRUE, envir=env)
        else assign("Rqualisummary", FALSE, envir=env)
        
        if(tclvalue(qualiValue)=="1") assign("Rquali", TRUE, envir=env)
        else assign("Rquali", FALSE, envir=env)
        
        if(tclvalue(qualisupValue)=="1") assign("Rqualisup", TRUE, envir=env)
        else assign("Rqualisup", FALSE, envir=env)
        
        if(tclvalue(axepartielValue)=="1") assign("Raxepartiel", TRUE, envir=env)
        else assign("Raxepartiel", FALSE, envir=env)
        
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
        
        groupe.lab <-tklabel(SortieWin, text=.Facto_gettext("Results for the groups"))
        groupe.check <- tkcheckbutton(SortieWin)
        if(Rgroupe) groupeValue <- tclVar("1")
        else groupeValue <- tclVar("0")
        tkconfigure(groupe.check,variable=groupeValue)
        
        inertie.lab <-tklabel(SortieWin, text=.Facto_gettext("Inertia"))
        inertie.check <- tkcheckbutton(SortieWin)
        if(Rinertie) inertieValue <- tclVar("1")
        else inertieValue <- tclVar("0")
        tkconfigure(inertie.check,variable=inertieValue)

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
        
      quantiSummary.lab<-tklabel(SortieWin,text=.Facto_gettext("Summary of the quantitative variables"))
        quantiSummary.check <- tkcheckbutton(SortieWin)
        if(Rquantisummary) quantiSummaryValue <- tclVar("1")
        else quantiSummaryValue <- tclVar("0")
        tkconfigure(quantiSummary.check,variable=quantiSummaryValue)
      
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
      
      qualiSummary.lab<-tklabel(SortieWin,text=.Facto_gettext("Summary of the qualitative variables"))
        qualiSummary.check <- tkcheckbutton(SortieWin)
        if(Rqualisummary) qualiSummaryValue <- tclVar("1")
        else qualiSummaryValue <- tclVar("0")
        tkconfigure(qualiSummary.check,variable=qualiSummaryValue)
      
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
      
      axepartiel.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for the partial axes"))
        axepartiel.check <- tkcheckbutton(SortieWin)
        if(Raxepartiel) axepartielValue <- tclVar("1")
        else axepartielValue <- tclVar("0")
        tkconfigure(axepartiel.check,variable=axepartielValue)
        
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
      tkgrid(inertie.lab,inertie.check,sticky="w")
      tkgrid(ind.lab,ind.check,sticky="w")
    nb.GQA<-length(listQuantiAct.nom)
    nb.GQI<-length(listQuantiIllu.nom)
    nb.GQlA<-length(listQualiAct.nom)
    nb.GQlI<-length(listQualiIllu.nom)
      if (!is.null(individuillu))tkgrid(ind.sup.lab,ind.sup.check,sticky="w")
      if (nb.GQA+ nb.GQI>0) tkgrid(quantiSummary.lab,quantiSummary.check,sticky="w")
      if (nb.GQA>0) tkgrid(quanti.lab,quanti.check,sticky="w")      
      if (nb.GQI>0) tkgrid(quantisup.lab,quantisup.check,sticky="w")
      if (nb.GQlA+ nb.GQlI>0) tkgrid(qualiSummary.lab,qualiSummary.check,sticky="w")
      if (nb.GQlA>0) tkgrid(quali.lab,quali.check,sticky="w")        
      if (nb.GQlI>0) tkgrid(qualisup.lab,qualisup.check,sticky="w")
      tkgrid(axepartiel.lab,axepartiel.check,sticky="w")      
      tkgrid(descdim.lab,descdim.check,sticky="w")
      tkgrid(tklabel(SortieWin, text = " "))
      tkgrid(RFichierFrame)
      tkgrid(SortieOK.but)
   }
    
    SortieFrame<-tkframe(IlluFrame)
    Sortie.but<-tkbutton(SortieFrame, textvariable=.SortieLabel, command=OnSortie, borderwidth=3)
    tkgrid(Sortie.but, sticky="ew")
  })



  #! fonction pour la gestion des options graphiques 
  PLOT.MFA<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.graph<-0
    .PlotLabel<-tclVar(paste(firstLabel, "", sep=" "))
    #déclaration des variables
    Gchoix<-TRUE
    GTitle<-NULL
    GAxeGrpe<-c(1,2)
    Glabel<-TRUE
        
    Rchoix<-TRUE
    RTitle<-NULL
    Rlabel.indMoy<-TRUE
    Rlabel.indPar<-TRUE
    Rlabel.quali<-TRUE    
    Rhabillage<-"group"
    Rinvisible<-NULL
    Rpartial<-NULL
    RpartialSouris<-FALSE
    Rchrono<-FALSE
    RXlimInd<-NULL
    RYlimInd<-NULL
    
    Wchoix<-TRUE
    WTitle<-NULL
    WAxeVar<-c(1,2)
    Wlabel.var<-TRUE
    Whabillage<-"group"
    Winvisible<-NULL
    Wlim.cos<-0.
    
    Achoix<-TRUE
    ATitle<-NULL
    AAxeAxe<-c(1,2)
    Ahabillage<-"group"
        
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

          #assign("GAxeGrpe", c(as.numeric(tclvalue(AxeGrpe1)), as.numeric(tclvalue(AxeGrpe2))), envir=env)
          
          label.tmp.grpe<-tclvalue(label.grpe.checkValue)
          if(label.tmp.grpe==1) assign("Glabel", TRUE, envir=env)
          else assign("Glabel", FALSE, envir=env)
        }
        
        # gestion des entrées de la partie graphique des Groupes
        if(tclvalue(axe.check.value)==1) assign("Achoix", TRUE, envir=env)
        else assign("Achoix", FALSE, envir=env)

        if(Achoix) {
          if (tclvalue(ATitre)==" ") assign("ATitle", NULL, envir=env)
          assign("ATitle", tclvalue(ATitre), envir=env)

          habillage.tmp.axe<-tclvalue(Ahabillage.checkValue)
          if(habillage.tmp.axe==1) assign("Ahabillage", "group", envir=env)
          else assign("Ahabillage", "none", envir=env)
        }
            
        # gestion des entrées de la partie graphique des variables
        if(tclvalue(var.check.value)==1) assign("Wchoix", TRUE, envir=env)
        else assign("Wchoix", FALSE, envir=env)

        if(Wchoix) {
          if (tclvalue(WTitre)==" ") assign("WTitle", NULL, envir=env)
          else assign("WTitle", tclvalue(WTitre), envir=env)

          assign("Wlim.cos", tclvalue(WlimCosValue), envir=env)

          if(tclvalue(label.var.checkValue)==1) assign("Wlabel.var", TRUE, envir=env)
          else assign("Wlabel.var", FALSE, envir=env)

          if(tclvalue(Whabillage.checkValue)==1) assign("Whabillage", "group", envir=env)
          else assign("Whabillage", "none", envir=env)
          
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
          label.tmp.indPar<-tclvalue(label.indPar.checkValue)
          label.tmp.quali<-tclvalue(label.quali.checkValue)
          if(label.tmp.indMoy==1) assign("Rlabel.indMoy", TRUE, envir=env)
          else assign("Rlabel.indMoy", FALSE, envir=env)
          if(label.tmp.indPar==1) assign("Rlabel.indPar", TRUE, envir=env)
          else assign("Rlabel.indPar", FALSE, envir=env)
          if(label.tmp.quali==1) assign("Rlabel.quali", TRUE, envir=env)
          else assign("Rlabel.quali", FALSE, envir=env)
          
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
      tkgrid(tklabel(GchoixFrame, text=.Facto_gettext("Graph of the groups"), font=font2),grpe.check)
      tkgrid(tklabel(GchoixFrame, text="  "))
  
      GTitleFrame<-tkframe(PlotGrpeFrame,borderwidth=2)
      if (is.null(GTitle)) GTitre <- tclVar(" ")
      else GTitre<-tclVar(GTitle)
      GTitre.entry <-tkentry(GTitleFrame,width="40",textvariable=GTitre)
      tkgrid(tklabel(GTitleFrame,text=.Facto_gettext("Title of the graph")),GTitre.entry)
    
      GlabelFrame<-tkframe(PlotGrpeFrame,borderwidth=2)
      label.grpe.check<-tkcheckbutton(GlabelFrame)
      if (Glabel) label.grpe.checkValue<-tclVar("1")
      else label.grpe.checkValue<-tclVar("0")
      tkconfigure(label.grpe.check, variable=label.grpe.checkValue)
      tkgrid(tklabel(GlabelFrame, text=.Facto_gettext("Labels for the groups")),label.grpe.check)
  
     
      #mise en page des différents frames de PlotGrpeFrame
      tkgrid(GchoixFrame)
      tkgrid(GTitleFrame)
      tkgrid(GlabelFrame)
      tkgrid(tklabel(PlotGrpeFrame, text=" "))
      
      
      ##########################
      # construction de la partie graphique des Axes partiels
      PlotAxeFrame<-tkframe(PlotWin2, borderwidth=5, relief="groove")
  
      AchoixFrame<-tkframe(PlotAxeFrame,borderwidth=2)
      axe.check<-tkcheckbutton(AchoixFrame)
      if(Achoix) axe.check.value<-tclVar("1")
      else axe.check.value<-tclVar("0")
      tkconfigure(axe.check, variable=axe.check.value)
      tkgrid(tklabel(AchoixFrame, text=.Facto_gettext("Graph of the partial axes"), font=font2),axe.check)
      tkgrid(tklabel(AchoixFrame, text="  "))
  
      ATitleFrame<-tkframe(PlotAxeFrame,borderwidth=2)
      if (is.null(ATitle)) ATitre <- tclVar(" ")
      else ATitre<-tclVar(ATitle)
      ATitre.entry <-tkentry(ATitleFrame,width="40",textvariable=ATitre)
      tkgrid(tklabel(ATitleFrame,text=.Facto_gettext("Title of the graph")),ATitre.entry)
  
      AhabillageFrame<-tkframe(PlotAxeFrame,borderwidth=2)
      Ahabillage.check<-tkcheckbutton(AhabillageFrame)
      if (Ahabillage=="group") Ahabillage.checkValue<-tclVar("1")
      else Ahabillage.checkValue<-tclVar("0")
      tkconfigure(Ahabillage.check, variable=Ahabillage.checkValue)
      tkgrid(tklabel(AhabillageFrame, text=.Facto_gettext("Color the partial axes by group")),Ahabillage.check)
           
      #mise en page des différents frames de PlotGrpeFrame
      tkgrid(AchoixFrame)
      tkgrid(ATitleFrame)
      tkgrid(AhabillageFrame)
      tkgrid(tklabel(PlotAxeFrame, text=" "))
      
      ########################

      # construction de la partie graphique des variables
      PlotVarFrame<-tkframe(PlotWin2, borderwidth=5, relief="groove")
  
      WchoixFrame<-tkframe(PlotVarFrame,borderwidth=2)
      var.check<-tkcheckbutton(WchoixFrame)
      if(Wchoix) var.check.value<-tclVar("1")
      else var.check.value<-tclVar("0")
      tkconfigure(var.check, variable=var.check.value)
      tkgrid(tklabel(WchoixFrame, text=.Facto_gettext("Graph of the variables"), font=font2),var.check)
      tkgrid(tklabel(WchoixFrame, text="  "))
  
      WTitleFrame<-tkframe(PlotVarFrame,borderwidth=2)
      if(is.null(WTitle)) WTitre <- tclVar(" ")
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
      
      WhabillageFrame<-tkframe(PlotVarFrame,borderwidth=2)
      Whabillage.check<-tkcheckbutton(WhabillageFrame)
      if (Whabillage=="group") Whabillage.checkValue<-tclVar("1")
      else Whabillage.checkValue<-tclVar("0")
      tkconfigure(Whabillage.check, variable=Whabillage.checkValue)
      tkgrid(tklabel(WhabillageFrame, text=.Facto_gettext("Color the variables by group")),Whabillage.check)
      
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
      tkgrid(tklabel(WinvisibleFrame, text=.Facto_gettext("None")),inv.aucun.check, tklabel(WinvisibleFrame, text=.Facto_gettext("active variables")),inv.act.check, tklabel(WinvisibleFrame, text=.Facto_gettext("supplementary variables")),inv.sup.check, sticky="w")   
        
      #mise en page des différents frames de PlotVarFrame
      tkgrid(WchoixFrame)
      tkgrid(WTitleFrame)
      tkgrid(WcosFrame)
      tkgrid(WlabelFrame)
      tkgrid(WhabillageFrame)
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
      label.indPar.check<-tkcheckbutton(RlabelFrame)
      if (Rlabel.indPar) label.indPar.checkValue<-tclVar("1")
      else label.indPar.checkValue<-tclVar("0")
      tkconfigure(label.indPar.check, variable=label.indPar.checkValue)
      tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Labels for the partial individuals")),label.indPar.check)
      label.quali.check<-tkcheckbutton(RlabelFrame)
      if (Rlabel.quali) label.quali.checkValue<-tclVar("1")
      else label.quali.checkValue<-tclVar("0")
      tkconfigure(label.quali.check, variable=label.quali.checkValue)
      tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Labels for the factors")), label.quali.check)
  
        
      RhabillageFrame<-tkframe(PlotIndFrame,borderwidth=2)
      listgraph<-tklistbox(RhabillageFrame,height=4, selectmode="single",exportselection="FALSE",yscrollcommand=function(...) tkset(scrgraph,...))
      scrgraph <-tkscrollbar(RhabillageFrame,repeatinterval=5,command=function(...)tkyview(listgraph,...))
      listgraph.nom<-c("group","ind")
      tkinsert(listgraph,"end",.Facto_gettext("by.group"))
      tkinsert(listgraph,"end",.Facto_gettext("by.individual"))
      if(Rhabillage=="group") tkselection.set(listgraph,0)
      if(Rhabillage=="ind") tkselection.set(listgraph,1)
      indice<-2      
##      for (i in 1:ncol(donnee)) {
##          if(is.factor(donnee[,i])) {
##          tkinsert(listgraph,"end",vars[i])
##          listgraph.nom<-c(listgraph.nom,vars[i])
##          if(Rhabillage==vars[i]) tkselection.set(listgraph, indice)
##          indice<-indice+1
##        }
##      }
    if(length(listQualiAct.nom)+length(listQualiIllu.nom)>=1) {
      var.aux = NULL
      if(length(listQualiAct.nom)>=1) {
        for(i in 1:length(listQualiAct.nom)) {
          eval(parse(text=paste("liste.aux.GQlA<-", listQualiAct.nom[i], ".var", sep="")))
          var.aux<-c(var.aux, liste.aux.GQlA[-1])
        }
      }
      if(length(listQualiIllu.nom)>=1) {
        for(i in 1:length(listQualiIllu.nom)) {
          eval(parse(text=paste("liste.aux.GQlI<-", listQualiIllu.nom[i], ".var", sep="")))
          var.aux<-c(var.aux, liste.aux.GQlI[-1])
        }
      }
      for (j in 1:ncol(donnee)){
        if(vars[j] %in% var.aux){
          tkinsert(listgraph,"end",vars[j])
          listgraph.nom<-c(listgraph.nom,vars[j])
          if(Rhabillage==vars[j]) tkselection.set(listgraph, indice)
          indice<-indice+1
        }
      }
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
      tkgrid(RlabelFrame)
      tkgrid(RinvisibleFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))
      tkgrid(RhabillageFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))
      tkgrid(RpartialFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))      
      tkgrid(RlimFrame)
      tkgrid(tklabel(PlotIndFrame, text=" "))
      
              
      #mise en page de plotWin
      subOKCancelHelp(PlotWin, "plot.MFA")
      tkgrid(PlotGrpeFrame)
      tkgrid(PlotAxeFrame)
      if(length(listQuantiAct.nom)+length(listQuantiIllu.nom)>=1) {
        tkgrid(PlotVarFrame)
        Wchoix = TRUE
      }
      tkgrid(PlotIndFrame, PlotWin2, sticky="ns")
      tkgrid(subButtonsFrame, sticky="ew", columnspan=2)
    }

    PlotFrame<-tkframe(IlluFrame)
    Plot.but<-tkbutton(PlotFrame, textvariable=.PlotLabel, command=OnPlot, borderwidth=3)
    tkgrid(Plot.but, sticky="ew")
  }) 


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
#      tkgrid(tklabel(HcpcWin,text=.Facto_gettext(paste('Clustering is performed on the first ', tclvalue(ncp.val), ' dimensions of the MFA',sep=""))),column=1,columnspan=4,sticky="w")
      tkgrid(tklabel(HcpcWin,text=sprintf(.Facto_gettext("Clustering is performed on the first %s dimensions of MFA"),tclvalue(ncp.val))),column=1,columnspan=4,sticky="w")   #text which takes the nb of dimensions chosen in the main window
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
#          Rinertie
#          Rindividu
#        Rindsup
#      Rquantisummary
#      Rquanti
#      Rquantisup  
#          Rqualisummary
#          Rquali
#          Rqualisup
#      Raxepartiel      
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
#      Rlabel.indPar
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
#      Whabillage
#      Winvisible
#      Wlim.cos
#    
#      Achoix
#      ATitle
#      AAxeAxe
#      Ahabillage
#
#      Axe 
  
    # récupération des paramètres de la fenêtre principale
    nom.res<-tclvalue(resu.val)
    if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) justDoIt(paste('remove (',nom.res,')'))       #if object res already exists, it's removed
    ncp<-as.numeric(tclvalue(ncp.val))
    Axe<-c(as.numeric(tclvalue(Axe1)), as.numeric(tclvalue(Axe2)))
    
    # gestion du tableau de données pour l'AFM
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
        group<-c(group, length(liste.var.GQA)-1)
        type<-c(type,liste.var.GQA[1])
        variables<-c(variables, liste.var.GQA[-1])
        indice.grpe<-indice.grpe+1
      }
    }
      #récupération des groupes quanti illustratif
    nb.GQI<-length(listQuantiIllu.nom)
    if(nb.GQI>=1) {
      name.group<-c(name.group, listQuantiIllu.nom)
      for(i in 1:nb.GQI) {
         eval(parse(text=paste("liste.var.GQI<-", listQuantiIllu.nom[i], ".var", sep="")))
        group<-c(group, length(liste.var.GQI)-1)
        type<-c(type,liste.var.GQI[1])
        variables<-c(variables, liste.var.GQI[-1])
        num.group.sup<-c(num.group.sup,indice.grpe)
        indice.grpe<-indice.grpe+1
      }
    }
      #récupération des groupes quali actif
    nb.GQlA<-length(listQualiAct.nom)
    if(nb.GQlA>=1) {
      name.group<-c(name.group, listQualiAct.nom)
      for(i in 1:nb.GQlA) {
        eval(parse(text=paste("liste.var.GQlA<-", listQualiAct.nom[i], ".var", sep="")))
        group<-c(group, length(liste.var.GQlA)-1)
        type<-c(type,liste.var.GQlA[1])
        variables<-c(variables, liste.var.GQlA[-1])
        indice.grpe<-indice.grpe+1
      }
    }
      #récupération des groupes quali illustratif
    nb.GQlI<-length(listQualiIllu.nom)
    if(nb.GQlI>=1) {
      name.group<-c(name.group, listQualiIllu.nom)
      for(i in 1:nb.GQlI) {
        eval(parse(text=paste("liste.var.GQlI<-", listQualiIllu.nom[i], ".var", sep="")))
        group<-c(group, length(liste.var.GQlI)-1)
        type<-c(type,liste.var.GQlI[1])
        variables<-c(variables, liste.var.GQlI[-1])
        num.group.sup<-c(num.group.sup,indice.grpe)
        indice.grpe<-indice.grpe+1
      }
    }
    
      #construction du tableau de données.MFA
      if(!is.null(individuillu)) {
        ind.actif<-rows[-which(rows %in% individuillu)]
        commande.data<-paste(activeDataSet(),'.MFA', '<-', activeDataSet(),'[c("', paste(ind.actif, collapse='", "'), '", "', paste(individuillu, collapse='", "'), '"), c("',paste(variables, collapse='", "'), '")]', sep="")
      }
      else commande.data<-paste(activeDataSet(),'.MFA', '<-', activeDataSet(),'[, c("',paste(variables, collapse='", "'), '")]', sep="")
      
      justDoIt(commande.data)
      logger(commande.data)
      donnee.depart<-activeDataSet()
      activeDataSet(paste(activeDataSet(),'.MFA', sep=""))

      # gestion de la commande réalisant l'AFM     
      if(!is.null(individuillu)) {
        ind.actif<-rows[-which(rows %in% individuillu)]
        commande.MFA<-paste(nom.res, '<-MFA(', activeDataSet(), ', group=c(',paste(group, collapse=", "), '), type=c("', paste(type, collapse='", "'),'"), ind.sup=', nrow(get(getRcmdr(".activeDataSet")))-length(individuillu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))), ', ncp=', ncp, ', name.group=c("',paste(name.group, collapse='", "'), '"), num.group.sup=c(',paste(num.group.sup, collapse=", "), '), graph=FALSE)',sep="")
      }
      else commande.MFA<-paste(nom.res, '<-MFA(', activeDataSet(), ', group=c(',paste(group, collapse=", "), '), type=c("', paste(type, collapse='", "'),'"), ncp=', ncp, ', name.group=c("',paste(name.group, collapse='", "'), '"), num.group.sup=c(',paste(num.group.sup, collapse=", "), '), graph=FALSE)',sep="")
      justDoIt(commande.MFA)
      logger(commande.MFA)
	  justDoIt(paste(nom.res,'$call$call <-',deparse(commande.MFA),sep=""))


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
        commande.plotG<-paste('plot.MFA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="group", new.plot=TRUE, lab.grpe=', Glabel, sep="")
        if (is.null(GTitle)) commande.plotG <- paste(commande.plotG,')', sep="")
        else {
          if (GTitle ==" ") commande.plotG <- paste(commande.plotG,')', sep="")
          else commande.plotG <- paste(commande.plotG,', title="', GTitle,'")', sep="")
        }
        justDoIt(commande.plotG)
        logger(commande.plotG)
      }}
      
      if((Achoix)&(length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0)){
      if (get(nom.res)$eig[1,2]!=100) {
        commande.plotA<-paste('plot.MFA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="axes", new.plot=TRUE, habillage="', Ahabillage, '"', sep="")
        if (is.null(ATitle)) commande.plotA <- paste(commande.plotA,')', sep="")
        else {
          if (ATitle ==" ") commande.plotA <- paste(commande.plotA,')', sep="")
          else commande.plotA <- paste(commande.plotA,', title="', ATitle,'")', sep="")
        }
        justDoIt(commande.plotA)
        logger(commande.plotA) 
      }}
      
      if((Wchoix)&(length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0)&(nb.GQI+nb.GQA>0)){
      if (get(nom.res)$eig[1,2]!=100) {
        commande.plotW<-paste('plot.MFA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="var", new.plot=TRUE, lab.var=', Wlabel.var, ', habillage="', Whabillage, '", lim.cos2.var=', Wlim.cos, sep="")
        if (!is.null(Winvisible)) commande.plotW<-paste(commande.plotW, ', invisible=c("', paste(Winvisible, collapse='", "'),'")', sep='')
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
        if ((Rhabillage!="none") & (Rhabillage!="ind") & (Rhabillage!="group")) {
          Rhabillage<-which(colnames(get(getRcmdr(".activeDataSet")))==Rhabillage)
          if(length(Rhabillage)==0) Rhabillage<-"none"
        }
        if (Rhabillage=="none") Rhabillage<-paste('"', Rhabillage, '"', sep="")
        if (Rhabillage=="ind") Rhabillage<-paste('"', Rhabillage, '"', sep="")
        if (Rhabillage=="group") Rhabillage<-paste('"', Rhabillage, '"', sep="")
        
        if(RpartialSouris){
          commande.plotI<-paste('plotMFApartial(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), lab.ind=', Rlabel.indMoy, ', lab.par=', Rlabel.indPar, ', lab.var=', Rlabel.quali,', habillage=', Rhabillage, sep="")
          if (!is.null(RXlimInd)) commande.plotI<-paste(commande.plotI, ', xlim=c(', paste(RXlimInd, collapse=", "), ')', sep='')
          if (!is.null(RYlimInd)) commande.plotI<-paste(commande.plotI, ', ylim=c(', paste(RYlimInd, collapse=", "), ')', sep='')
          if (!is.null(Rinvisible)) commande.plotI<-paste(commande.plotI, ', invisible=c("', paste(Rinvisible, collapse='", "'),'")', sep='')
          if (Rchrono) commande.plotI<-paste(commande.plotI, ', chrono=', Rchrono, sep='')
          if (is.null(RTitle)) commande.plotI <- paste(commande.plotI,')', sep='')
          else {
            if (RTitle ==" ") commande.plotI <- paste(commande.plotI,')', sep="")
            else commande.plotI <- paste(commande.plotI,', title="', RTitle,'")', sep="")
          }
        }
        else {
          commande.plotI<-paste('plot.MFA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="ind", new.plot=TRUE, lab.ind=', Rlabel.indMoy, ', lab.par=', Rlabel.indPar, ', lab.var=', Rlabel.quali,', habillage=', Rhabillage, sep="") 
          if (!is.null(RXlimInd)) commande.plotI<-paste(commande.plotI, ', xlim=c(', paste(RXlimInd, collapse=", "), ')', sep='')
          if (!is.null(RYlimInd)) commande.plotI<-paste(commande.plotI, ', ylim=c(', paste(RYlimInd, collapse=", "), ')', sep='')
          if (!is.null(Rinvisible)) commande.plotI<-paste(commande.plotI, ', invisible=c("', paste(Rinvisible, collapse='", "'),'")', sep='')
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
      }}
      
      # gestion de l'édition de certains resultats
    doItAndPrint(paste('summary(',nom.res,', nb.dec = 3, nbelements=10, nbind = 10, ncp = 3, file="")', sep=""))
    if (RFichier==""){
      if(Rpropre) doItAndPrint(paste( nom.res, '$eig', sep=""))
      if(Rgroupe) doItAndPrint(paste( nom.res, '$group', sep=""))
      if(Rinertie) doItAndPrint(paste( nom.res, '$inertia.ratio', sep=""))
      if(Rindividu) doItAndPrint(paste( nom.res, '$ind', sep=""))
      if(Rindsup) doItAndPrint(paste( nom.res, '$ind.sup', sep=""))
      if(Rquantisummary) doItAndPrint(paste( nom.res, '$summary.quanti', sep=""))
      if(Rquanti) doItAndPrint(paste( nom.res, '$quanti.var', sep=""))      
      if(Rquantisup) doItAndPrint(paste( nom.res, '$quanti.var.sup', sep=""))
      if(Rqualisummary) doItAndPrint(paste( nom.res, '$summary.quali', sep=""))
      if(Rquali) doItAndPrint(paste( nom.res, '$quali.var', sep=""))      
      if(Rqualisup) doItAndPrint(paste( nom.res, '$quali.var.sup', sep=""))
      if(Raxepartiel) doItAndPrint(paste( nom.res, '$partial.axes', sep=""))
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
      if(Rinertie){
        doItAndPrint(paste('write.infile(', nom.res, '$inertia.ratio, file =',Fich,',append=',append,')', sep=""))
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
      if(Rquantisummary){
        doItAndPrint(paste('write.infile(', nom.res, '$summary.quanti, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rquanti){
        doItAndPrint(paste('write.infile(', nom.res, '$quanti.var, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rquantisup){
        doItAndPrint(paste('write.infile(', nom.res, '$quanti.var.sup, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rqualisummary){
        doItAndPrint(paste('write.infile(', nom.res, '$summary.quali, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rquali){
        doItAndPrint(paste('write.infile(', nom.res, '$quali.var, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rqualisup){
        doItAndPrint(paste('write.infile(', nom.res, '$quali.var.sup, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Raxepartiel){
        doItAndPrint(paste('write.infile(', nom.res, '$partial.axes, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rdescdim) doItAndPrint(paste('write.infile(dimdesc(', nom.res, ', axes=1:',ncp,'), file =',Fich,',append=',append,')', sep=""))
    }

      # Re-chargement du tableau de départ et supression du tableau temporaire
      activeDataSet(donnee.depart)
      justDoIt(paste('remove(',activeDataSet(),'.MFA)',sep=""))
      logger(paste('remove(',activeDataSet(),'.MFA)',sep=""))   
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
  tkwm.title(top,.Facto_gettext("MFA"))
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
  label.quantiFrame.var<-tclVar(.Facto_gettext("Quantitative groups"))
  label.quantiFrame<-tklabel(ListeQuantiFrame, textvariable=label.quantiFrame.var,fg = "darkred", font=fontheading)
    # liste des groupes de variables quanti Actives
  listQuantiAct<-tklistbox(ListeQuantiFrame,selectmode="extended",exportselection="FALSE", height=4, yscrollcommand=function(...)tkset(scrQuantiAct,...))
  scrQuantiAct<-tkscrollbar(ListeQuantiFrame,repeatinterval=5,command=function(...)tkyview(listQuantiAct,...))
  listQuantiAct.nom<-NULL
  
    # liste des groupes de variables quanti Actives
  listQuantiIllu<-tklistbox(ListeQuantiFrame,selectmode="extended",exportselection="FALSE", height=4, yscrollcommand=function(...)tkset(scrQuantiIllu,...))
  scrQuantiIllu<-tkscrollbar(ListeQuantiFrame,repeatinterval=5,command=function(...)tkyview(listQuantiIllu,...))
  listQuantiIllu.nom<-NULL
  
    # boutons d'action groupes quantitative
  supprimeQuanti.funct(label=.Facto_gettext("Delete"))  
  ajoutQuanti.funct(label=.Facto_gettext("Add quanti. group"), firstLabel=.Facto_gettext("Add quanti group"))
  modifQuanti.funct(label=.Facto_gettext("Modify 1 group"), firstLabel=.Facto_gettext("Modify 1 group"))
    # mise en forme de ListeQuantiFrame
  #tkgrid(tklabel(ListeQuantiFrame, text = "Quantitative groups",fg = "darkred"), columnspan=11, sticky = "ew")
  tkgrid(label.quantiFrame, columnspan=11, sticky = "ew")
  tkgrid(listQuantiAct, scrQuantiAct, listQuantiIllu, scrQuantiIllu)
  tkgrid.configure(scrQuantiAct, column=3, sticky="wns")
  tkgrid.configure(scrQuantiIllu, column=9, sticky="wns")  
  tkgrid.configure(listQuantiAct, sticky = "ew", column=1, columnspan=2)
  tkgrid.configure(listQuantiIllu, sticky = "ew", column=7, columnspan=2)
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
  
 # création du frame contenant les listes groupes quali
  ListeQualiFrame<- tkframe(top, borderwidth=2, relief="groove")
  label.qualiFrame.var<-tclVar(.Facto_gettext("Qualitative groups"))
  label.qualiFrame<-tklabel(ListeQualiFrame, textvariable=label.qualiFrame.var,fg = "darkred", font=fontheading)
    # liste des groupes de variables quali Actives
  listQualiAct<-tklistbox(ListeQualiFrame,selectmode="extended",exportselection="TRUE", height=4, yscrollcommand=function(...)tkset(scrQualiAct,...))
  scrQualiAct<-tkscrollbar(ListeQualiFrame,repeatinterval=5,command=function(...)tkyview(listQualiAct,...))
  listQualiAct.nom<-NULL
  
  # liste des groupes de variables quali Actives
  listQualiIllu<-tklistbox(ListeQualiFrame,selectmode="extended",exportselection="TRUE", height=4, yscrollcommand=function(...)tkset(scrQualiIllu,...))
  scrQualiIllu<-tkscrollbar(ListeQualiFrame,repeatinterval=5,command=function(...)tkyview(listQualiIllu,...))
  listQualiIllu.nom<-NULL

      # boutons d'action groupes qualitatif
  supprimeQuali.funct(label=.Facto_gettext("Delete"))
  ajoutQuali.funct(label=.Facto_gettext("Add quali. group"), firstLabel=.Facto_gettext("Add quali group")) 
  modifQuali.funct(label=.Facto_gettext("Modify 1 group"), firstLabel=.Facto_gettext("Modify 1 group"))
  
    # mise en forme de ListeQualiFrame
  #tkgrid(tklabel(ListeQualiFrame, text = "Qualitative groups",fg = "darkred"), columnspan=11, sticky = "ew")
  tkgrid(label.qualiFrame, columnspan=11, sticky = "ew")
  tkgrid(listQualiAct, scrQualiAct, listQualiIllu, scrQualiIllu)
  tkgrid.configure(scrQualiAct, column=3, columnspan=1, sticky="wns")
  tkgrid.configure(scrQualiIllu, column=9, columnspan=1, sticky="wns")  
  tkgrid.configure(listQualiAct, sticky = "ew", column=1, columnspan=2)
  tkgrid.configure(listQualiIllu, sticky = "ew", column=7, columnspan=2)
  tkgrid.configure(tklabel(ListeQualiFrame, text=" "))
  tkgrid.configure(GpeQualiFrame, ModifGpeQualiFrame, SupGpeQualiFrame)
  tkgrid.configure(GpeQualiFrame, sticky = "ew", column=1, columnspan=2)
  tkgrid.configure(ModifGpeQualiFrame, sticky = "ew", column=4, columnspan=2)  
  tkgrid.configure(SupGpeQualiFrame, sticky = "ew", column=7, columnspan=2)
  tkgrid.columnconfigure(ListeQualiFrame,0, minsize=25)
  tkgrid.columnconfigure(ListeQualiFrame,10, minsize=25)
  tkgrid.columnconfigure(ListeQualiFrame,3, minsize=25)
  tkgrid.columnconfigure(ListeQualiFrame,9, minsize=25)
  tkgrid.columnconfigure(ListeQualiFrame,4, minsize=35)
  tkgrid.columnconfigure(ListeQualiFrame,5, minsize=35)


   # création de tous les boutons d'options dans IlluFrame
  IlluFrame<- tkframe(top, borderwidth=2)
    # mise en page de IlluFrame
  Iillu.funct(label=.Facto_gettext("Supplementary individuals"), firstLabel=.Facto_gettext("Supplementary individuals"))    
  PLOT.MFA(label=.Facto_gettext("Graphical options"), firstLabel=.Facto_gettext("Graphical options"))
  Sortie.funct(label=.Facto_gettext("Outputs"), firstLabel=.Facto_gettext("Outputs"))
  tkgrid(IilluFrame, PlotFrame, SortieFrame, columnspan=7)
  tkgrid.configure(PlotFrame, column=3, columnspan=1)
  tkgrid.configure(SortieFrame, column=5, columnspan=1)
  tkgrid.configure(IilluFrame, column=1, columnspan=1)
  tkgrid.columnconfigure(IlluFrame,0, minsize=25)
  tkgrid.columnconfigure(IlluFrame,2, minsize=40)  
  tkgrid.columnconfigure(IlluFrame,4, minsize=25)  


    # création des options dans OptionFrame  
  OptionFrame<-tkframe(top, borderwidth=2, relief="groove")
  resu.lab<-tklabel(OptionFrame,text=.Facto_gettext("Name of the result object:"))
  resu.val<-tclVar("res")
  resu<-tkentry(OptionFrame,width=10,textvariable=resu.val)
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
  tkgrid(resu.lab, resu)
  tkgrid(ncp.lab, ncp)
  tkgrid(Axe.label,Axe1.entry , Axe2.entry, sticky="w")
  tkgrid.configure(ncp.lab, resu.lab, Axe.label, column=1, columnspan=4, sticky="w")
  tkgrid.configure(ncp, resu, column=6, columnspan=2, sticky="e")
  tkgrid.configure(Axe1.entry, column=6, columnspan=1, sticky="w")
  tkgrid.configure(Axe2.entry, column=7, columnspan=1, sticky="e")
  tkgrid.columnconfigure(OptionFrame,0, minsize=25)
  tkgrid.columnconfigure(OptionFrame,5, minsize=40)
  tkgrid.columnconfigure(OptionFrame,8, minsize=25)

  #Frame pour HCPC
  HcpcFrame<-tkframe(top, borderwidth=2)
  Hcpc.funct(label=.Facto_gettext("Perform Clustering after MFA"), firstLabel=.Facto_gettext("Perform Clustering after MFA")) 
  tkgrid(Hcpc2Frame, columnspan=7)
  tkgrid.configure(Hcpc2Frame,column=4, columnspan=1)


  appliquer.but<-tkbutton(top, text=.Facto_gettext("Apply"),width=12,command=OnAppliquer, borderwidth=3, fg="#690f96")
  OKCancelHelp(helpSubject="MFA",reset="Reinitializ.funct")


  # Mise en page de top
  tkgrid(tklabel(top, text=.Facto_gettext("Multiple Factor Analysis (MFA)"),font=fontheading), columnspan=3)
  tkgrid(tklabel(top,text=""))
  if (length(listNumeric())>0){
    tkgrid(ListeQuantiFrame, column=1, columnspan=1, sticky="ew")
    tkgrid(tklabel(top,text=""))    
  }
  if (length(listFactors())>0) {
    tkgrid(ListeQualiFrame, column=1, columnspan=1, sticky="ew")
    tkgrid(tklabel(top,text=""))    
  }
  tkgrid(IlluFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text=""))    
  tkgrid(OptionFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text="")) # Ligne de blanc  
  tkgrid(HcpcFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text="")) # Ligne de blanc        
  # tkgrid(appliquer.but, column=1, columnspan=1)
  # tkgrid(tklabel(top,text="")) # Ligne de blanc  
  # tkgrid(buttonsFrame, column=1, columnspan=1, sticky="ew" )
  tkgrid(buttonsFrame, appliquer.but)
  tkgrid.configure(buttonsFrame, column=1,sticky="e")
  tkgrid.configure(appliquer.but, column=2,sticky="w")
  
}
