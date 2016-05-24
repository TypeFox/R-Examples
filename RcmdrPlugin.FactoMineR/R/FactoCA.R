FactoCA <- function(){

################################################################################
#    Création des fonctions pour les options via nouvelle fenêtre graphique    #
################################################################################

 #! fonction pour le choix des variables colonnes supplémentaires 
  Cillu.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    variableColIllu<-NULL
    .CilluLabel<-tclVar(paste(firstLabel, "", sep=" "))
    OnCillu<-function()
    { 
      CilluWin<-tktoplevel()
      tkwm.title(CilluWin,.Facto_gettext("Select supplementary column(s)"))
      #création de la fonction COK.funct
      COK.funct<-function()
      {
        vsup.select<-listvar.nom[as.numeric(tkcurselection(listvar))+1]
        if(length(vsup.select)==0)
        {
          assign("variableColIllu", NULL, envir=env)
          tclvalue(.CilluLabel)<-paste(firstLabel, "", sep=" ")
          tkconfigure(Cillu.but, fg="black")
          tkdestroy(CilluWin)
          return()
        }
        assign("variableColIllu", vsup.select, envir=env)
#        tclvalue(.CilluLabel)<-paste(label, ": OK", sep=" ")
        tclvalue(.CilluLabel)<-paste(label, "", sep=" ")
        tkconfigure(Cillu.but, fg="blue")
        tkdestroy(CilluWin)
      }
      
      # création et mise en page de la fenetre Cillu
      listvar<-tklistbox(CilluWin,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrvar,...)) # Liste vide
      scrvar <-tkscrollbar(CilluWin,repeatinterval=5,command=function(...)tkyview(listvar,...)) 
      listvar.nom<-NULL
      indice<-0
      for (i in (1:ncol(donnee))){
        if (is.numeric(donnee[,i])) {
          tkinsert(listvar,"end",vars[i]) # On renseigne la liste
          listvar.nom<-c(listvar.nom,vars[i])
          if(vars[i] %in% variableColIllu) tkselection.set(listvar, indice)
          indice<-indice+1
        }
      }
  
      COK.but<-tkbutton(CilluWin, text="OK", width=16,command=COK.funct)

      tkgrid(tklabel(CilluWin, text=""))
        tkgrid(tklabel(CilluWin, text = .Facto_gettext("Select supplementary column(s)"), fg = "blue"), column=1, columnspan = 1, sticky = "ew")
        tkgrid(listvar, scrvar, sticky = "nw")
        tkgrid.configure(scrvar, sticky = "ens", columnspan=1)
        tkgrid.configure(listvar, sticky = "ew", column=1, columnspan=1)
        tkgrid(tklabel(CilluWin, text=""))
        tkgrid(COK.but, column=1,columnspan=1, sticky="ew")
        tkgrid(tklabel(CilluWin, text=""))
        tkgrid.columnconfigure(CilluWin,0, minsize=25)
      tkgrid.columnconfigure(CilluWin,2, minsize=25)
  }  

   CilluFrame<-tkframe(IlluFrame)
   Cillu.but<-tkbutton(CilluFrame, textvariable=.CilluLabel, command=OnCillu, borderwidth=3, width=35)
   tkgrid(Cillu.but, sticky="ew")
  })
  

  #! fonction pour le choix des variables lignes supplémentaires 
  Lillu.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    variableLigneIllu<-NULL
    .LilluLabel<-tclVar(paste(firstLabel, "", sep=" "))
    OnLillu<-function()
    {   
      LilluWin<-tktoplevel()
      tkwm.title(LilluWin,.Facto_gettext("Select supplementary rows"))
      #création de la fonction LOK.funct
      LOK.funct<-function()
      {
        Ligne.select<-rows[as.numeric(tkcurselection(listLigne))+1]
        if(length(Ligne.select)==0)
        {
          assign("individuillu", NULL, envir=env)
          tclvalue(.LilluLabel)<-paste(firstLabel, "", sep=" ")
          tkconfigure(Lillu.but, fg="black")
          tkdestroy(LilluWin)
          return()
        }
        assign("variableLigneIllu", Ligne.select, envir=env)
#        tclvalue(.LilluLabel)<-paste(label, ": OK", sep=" ")
        tclvalue(.LilluLabel)<-paste(label, "", sep=" ")
        tkconfigure(Lillu.but, fg="blue")
        tkdestroy(LilluWin)
      }
      
      # création et mise en page de la fenetre Lillu
      listLigne<-tklistbox(LilluWin,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrLigne,...)) # Liste vide
      scrLigne <-tkscrollbar(LilluWin,repeatinterval=5,command=function(...)tkyview(listLigne,...)) 
      indice<-0
      for (i in (1:nrow(donnee)))
      {
          tkinsert(listLigne,"end",rows[i]) # On renseigne la liste
          if(rows[i] %in% variableLigneIllu) tkselection.set(listLigne, indice)
          indice<-indice+1
        }
  
      LOK.but<-tkbutton(LilluWin, text="OK", width=16,command=LOK.funct)

      tkgrid(tklabel(LilluWin, text=""))
        tkgrid(tklabel(LilluWin, text = .Facto_gettext("Select supplementary row(s)"), fg = "blue"), column=1, columnspan = 1, sticky = "ew")
        tkgrid(listLigne, scrLigne, sticky = "nw")
        tkgrid.configure(scrLigne, sticky = "ens", columnspan=1)
        tkgrid.configure(listLigne, sticky = "ew", column=1, columnspan=1)
        tkgrid(tklabel(LilluWin, text=""))
        tkgrid(LOK.but, column=1,columnspan=1, sticky="ew")
        tkgrid(tklabel(LilluWin, text=""))
        tkgrid.columnconfigure(LilluWin,0, minsize=25)
      tkgrid.columnconfigure(LilluWin,2, minsize=25)
  }  

   LilluFrame<-tkframe(IlluFrame)
   Lillu.but<-tkbutton(LilluFrame, textvariable=.LilluLabel, command=OnLillu, borderwidth=3, width=35)
   tkgrid(Lillu.but, sticky="ew")
  })



  #! fonction pour la gestion des options graphiques 
  PLOT.CA<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.graph<-0
    #déclaration des variables
    Rchoix<-TRUE
    RTitle<-NULL
    Rlabel<-c("col", "col.sup", "row", "row.sup")
    Rcol.col<-Rcol.col.tmp<-"blue"
    Rcol.col.sup<-Rcol.col.sup.tmp<-"darkblue"
    Rcol.row<-Rcol.row.tmp<-"red"
    Rcol.row.sup<-Rcol.row.sup.tmp<-"darkred"
    Rinvis<-c("")
    RXlimInd<-NULL
    RYlimInd<-NULL
	Rell <- c("")
    
    .PlotLabel<-tclVar(paste(firstLabel, "", sep=" "))
    
    OnPlot<-function()
    {
      PlotWin<-tktoplevel()
      tkwm.title(PlotWin,.Facto_gettext("Outputs"))
      tkwm.geometry(PlotWin, "-100+50")
      
      #création de la fonction onOKsub
      onOKsub<-function()
      {
        assign("compteur.graph", compteur.graph+1, envir=env)
#        if(compteur.graph>0) tclvalue(.PlotLabel)<-paste(label, ": Seen", sep=" ")
        if(compteur.graph>0) tclvalue(.PlotLabel)<-paste(label, "", sep=" ")
        tkconfigure(Plot.but, fg="blue")

        # gestion des entrées de la partie graphique
        if(tclvalue(ind.check.value)==1) assign("Rchoix", TRUE, envir=env)
        else assign("Rchoix", FALSE, envir=env)

        if(Rchoix)
        {
          if (tclvalue(Titre)==" ") assign("RTitle", NULL, envir=env)
          assign("RTitle", tclvalue(Titre), envir=env)

          label.tmp.col<-tclvalue(label.col.checkValue)
          label.tmp.col.sup<-tclvalue(label.col.sup.checkValue)
          label.tmp.row<-tclvalue(label.row.checkValue)
          label.tmp.row.sup<-tclvalue(label.row.sup.checkValue)
          ell.tmp.row<-tclvalue(ell.row.checkValue)
          ell.tmp.col<-tclvalue(ell.col.checkValue)
          assign("Rlabel", NULL, envir=env)
          if(label.tmp.col==1) assign("Rlabel", c(Rlabel, "col"), envir=env)
          if(label.tmp.col.sup==1) assign("Rlabel", c(Rlabel, "col.sup"), envir=env)
          if(label.tmp.row==1) assign("Rlabel", c(Rlabel, "row"), envir=env)
          if(label.tmp.row.sup==1) assign("Rlabel", c(Rlabel, "row.sup"), envir=env)
          assign("Rell", NULL, envir=env)
          if(ell.tmp.row==1) assign("Rell", c(Rell, "row"), envir=env)
          if(ell.tmp.col==1) assign("Rell", c(Rell, "col"), envir=env)
          
          invis.tmp.row<-tclvalue(invis.row.checkValue)
          invis.tmp.row.sup<-tclvalue(invis.row.sup.checkValue)
          invis.tmp.col.sup<-tclvalue(invis.col.sup.checkValue)
          invis.tmp.col<-tclvalue(invis.col.checkValue)
          assign("Rinvis", NULL, envir=env)
          if(invis.tmp.row==0) assign("Rinvis", c(Rinvis, "row"), envir=env)
          if(invis.tmp.row.sup==0) assign("Rinvis", c(Rinvis, "row.sup"), envir=env)
          if(invis.tmp.col.sup==0) assign("Rinvis", c(Rinvis, "col.sup"), envir=env)
          if(invis.tmp.col==0) assign("Rinvis", c(Rinvis, "col"), envir=env)

          assign("Rcol.col", Rcol.col.tmp, envir=env)
          assign("Rcol.col.sup", Rcol.col.sup.tmp, envir=env)
          assign("Rcol.row", Rcol.row.tmp, envir=env)
          assign("Rcol.row.sup", Rcol.row.sup.tmp, envir=env)

          if(tclvalue(XlimIndMin)=="" | tclvalue(XlimIndMax)=="") assign("RXlimInd", NULL, envir=env)
          else assign("RXlimInd", c(as.numeric(tclvalue(XlimIndMin)), as.numeric(tclvalue(XlimIndMax))), envir=env)
          if(tclvalue(YlimIndMin)=="" | tclvalue(YlimIndMax)=="") assign("RYlimInd", NULL, envir=env)
          else assign("RYlimInd", c(as.numeric(tclvalue(YlimIndMin)), as.numeric(tclvalue(YlimIndMax))), envir=env)
        }
        
         tkdestroy(PlotWin)
      }
        
      RchoixFrame<-tkframe(PlotWin,borderwidth=2)
      ind.check<-tkcheckbutton(RchoixFrame)
      if(Rchoix) ind.check.value<-tclVar("1")
      else ind.check.value<-tclVar("0")
      tkconfigure(ind.check, variable=ind.check.value)
      tkgrid(tklabel(RchoixFrame, text=.Facto_gettext("Graphical output"), font=font2),ind.check)
      tkgrid(tklabel(RchoixFrame, text="  "))
      
      RTitleFrame<-tkframe(PlotWin,borderwidth=2)
      if (is.null(RTitle)) Titre <- tclVar(" ")
      else Titre<-tclVar(RTitle)
      Titre.entry <-tkentry(RTitleFrame,width="40",textvariable=Titre)
      tkgrid(tklabel(RTitleFrame,text=.Facto_gettext("Title of the graph")),Titre.entry)
        
      RlabelFrame<-tkframe(PlotWin,borderwidth=2)
      tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("  ")),tklabel(RlabelFrame,text=.Facto_gettext("Plot")),tklabel(RlabelFrame, text=.Facto_gettext("Label")))
      label.row.check<-tkcheckbutton(RlabelFrame)
      if ("row" %in% Rlabel) label.row.checkValue<-tclVar("1")
      else label.row.checkValue<-tclVar("0") 
      tkconfigure(label.row.check, variable=label.row.checkValue)
      invis.row.check<-tkcheckbutton(RlabelFrame)
      if ("row" %in% Rinvis) invis.row.checkValue<-tclVar("0")
      else invis.row.checkValue<-tclVar("1") 
      tkconfigure(invis.row.check, variable=invis.row.checkValue)
      tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Active rows")),invis.row.check,label.row.check)
      label.row.sup.check<-tkcheckbutton(RlabelFrame)
      if ("row.sup" %in% Rlabel) label.row.sup.checkValue<-tclVar("1")
      else label.row.sup.checkValue<-tclVar("0")
      tkconfigure(label.row.sup.check, variable=label.row.sup.checkValue)
      invis.row.sup.check<-tkcheckbutton(RlabelFrame)
      if ("row.sup" %in% Rinvis) invis.row.sup.checkValue<-tclVar("0")
      else invis.row.sup.checkValue<-tclVar("1")
      tkconfigure(invis.row.sup.check, variable=invis.row.sup.checkValue)
      if(!is.null(variableLigneIllu)) tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Supplementary rows")), invis.row.sup.check, label.row.sup.check)
      label.col.check<-tkcheckbutton(RlabelFrame)
      if ("col" %in% Rlabel) label.col.checkValue<-tclVar("1")
      else label.col.checkValue<-tclVar("0") 
      tkconfigure(label.col.check, variable=label.col.checkValue)
      invis.col.check<-tkcheckbutton(RlabelFrame)
      if ("col" %in% Rinvis) invis.col.checkValue<-tclVar("0")
      else invis.col.checkValue<-tclVar("1") 
      tkconfigure(invis.col.check, variable=invis.col.checkValue)
      tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Active columns")), invis.col.check,label.col.check)
      label.col.sup.check<-tkcheckbutton(RlabelFrame)
      if ("col.sup" %in% Rlabel) label.col.sup.checkValue<-tclVar("1")
      else label.col.sup.checkValue<-tclVar("0")
      tkconfigure(label.col.sup.check, variable=label.col.sup.checkValue)
      invis.col.sup.check<-tkcheckbutton(RlabelFrame)
      if ("col.sup" %in% Rinvis) invis.col.sup.checkValue<-tclVar("0")
      else invis.col.sup.checkValue<-tclVar("1")
      tkconfigure(invis.col.sup.check, variable=invis.col.sup.checkValue)
      if(!is.null(variableColIllu)) tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Supplementary columns")),invis.col.sup.check,label.col.sup.check)
  
      RellipseFrame<-tkframe(PlotWin, borderwidth=2)
      ell.row.check<-tkcheckbutton(RellipseFrame)
      if ("row" %in% Rell) ell.row.checkValue<-tclVar("1")
      else ell.row.checkValue<-tclVar("0") 
      tkconfigure(ell.row.check, variable=ell.row.checkValue)
      tkgrid(tklabel(RellipseFrame, text=.Facto_gettext("Confidence ellipses around rows")),ell.row.check)
      ell.col.check<-tkcheckbutton(RellipseFrame)
      if ("col" %in% Rell) ell.col.checkValue<-tclVar("1")
      else ell.col.checkValue<-tclVar("0") 
      tkconfigure(ell.col.check, variable=ell.col.checkValue)
      tkgrid(tklabel(RellipseFrame, text=.Facto_gettext("Confidence ellipses around columns")),ell.col.check)

      RcolFrame<-tkframe(PlotWin, borderwidth=2)
      Rcol.row.value <- Rcol.row
      canvas.row <- tkcanvas(RcolFrame,width="80",height="25",bg=Rcol.row.value)
      ChangeColor.row <- function()
      {
        Rcol.row.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Rcol.row.value,title="Select color"))
        if (nchar(Rcol.row.value)>0)
        {
          tkconfigure(canvas.row,bg=Rcol.row.value)
          assign("Rcol.row.tmp", Rcol.row.value, envir=env)
        }
      }
      ChangeColor.row.button <- tkbutton(RcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.row)
      tkgrid(tklabel(RcolFrame, text=.Facto_gettext("Color for active rows")),canvas.row,ChangeColor.row.button, sticky="w")
      
      Rcol.row.sup.value<-Rcol.row.sup
      canvas.row.sup <- tkcanvas(RcolFrame,width="80",height="25",bg=Rcol.row.sup.value)
      ChangeColor.row.sup <- function()
      {
        Rcol.row.sup.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Rcol.row.sup.value,title="Select color"))
        if (nchar(Rcol.row.sup.value)>0)
        {
          tkconfigure(canvas.row.sup,bg=Rcol.row.sup.value)
          assign("Rcol.row.sup.tmp", Rcol.row.sup.value, envir=env)
        }
      }
      ChangeColor.row.sup.button <- tkbutton(RcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.row.sup)
      if(!is.null(variableLigneIllu)) tkgrid(tklabel(RcolFrame, text=.Facto_gettext("Color for supplementary rows")),canvas.row.sup,ChangeColor.row.sup.button, sticky="w")
      Rcol.col.value <- Rcol.col
      canvas.col <- tkcanvas(RcolFrame,width="80",height="25",bg=Rcol.col.value)
      ChangeColor.col <- function()
      {
        Rcol.col.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Rcol.col.value,title="Select color"))
        if (nchar(Rcol.col.value)>0)
        {
          tkconfigure(canvas.col,bg=Rcol.col.value)
          assign("Rcol.col.tmp", Rcol.col.value, envir=env)
        }
      }
      ChangeColor.col.button <- tkbutton(RcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.col)
      tkgrid(tklabel(RcolFrame, text=.Facto_gettext("Color for active columns")),canvas.col,ChangeColor.col.button, sticky="w")
      
      Rcol.col.sup.value<-Rcol.col.sup
      canvas.col.sup <- tkcanvas(RcolFrame,width="80",height="25",bg=Rcol.col.sup.value)
      ChangeColor.col.sup <- function()
      {
        Rcol.col.sup.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Rcol.col.sup.value,title="Select color"))
        if (nchar(Rcol.col.sup.value)>0)
        {
          tkconfigure(canvas.col.sup,bg=Rcol.col.sup.value)
          assign("Rcol.col.sup.tmp", Rcol.col.sup.value, envir=env)
        }
      }
      ChangeColor.col.sup.button <- tkbutton(RcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.col.sup)
      if(!is.null(variableColIllu)) tkgrid(tklabel(RcolFrame, text=.Facto_gettext("Color for supplementary columns")),canvas.col.sup,ChangeColor.col.sup.button, sticky="w")
                 
      RlimFrame<-tkframe(PlotWin, borderwidth=2)
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
      tkgrid(tklabel(PlotWin, text=" "))
      tkgrid(RTitleFrame)
      tkgrid(tklabel(PlotWin, text=" "))
      tkgrid(RlabelFrame)
      tkgrid(tklabel(PlotWin, text=" "))
      tkgrid(RellipseFrame)
      tkgrid(tklabel(PlotWin, text=" "))
      tkgrid(RcolFrame)
      tkgrid(tklabel(PlotWin, text=" "))
      tkgrid(RlimFrame)
      tkgrid(tklabel(PlotWin, text=" "))
      
      subOKCancelHelp(PlotWin, "plot.CA")
      tkgrid(subButtonsFrame, sticky="ew")
    }
    PlotFrame<-tkframe(IlluFrame)
    Plot.but<-tkbutton(PlotFrame, textvariable=.PlotLabel, command=OnPlot, borderwidth=3, width=35)
    tkgrid(Plot.but, sticky="ew")
  })

    
     #! fonction pour le choix des éléments de sortie
  Sortie.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.sortie<-0
    #déclaration des variables
    Rpropre<-FALSE
    RFichier <- ""
    Rcol<-FALSE
    Rcolsup<-FALSE
    Rrow<-FALSE
    Rrowsup<-FALSE
    Rdescdim<-FALSE

    .SortieLabel<-tclVar(paste(firstLabel, "", sep=" "))

    OnSortie<-function()
    {
      SortieWin<-tktoplevel()
      tkwm.title(SortieWin,.Facto_gettext("Outputs options"))

      #création de la fonction onOKsub
      onOK.sortie<-function()
      {
        assign("compteur.sortie", compteur.sortie+1, envir=env)
#        if(compteur.sortie>0) tclvalue(.SortieLabel)<-paste(label, ": Seen", sep=" ")
        if(compteur.sortie>0) tclvalue(.SortieLabel)<-paste(label, "", sep=" ")
        tkconfigure(Sortie.but, fg="blue")
        
        if(tclvalue(eigValue)=="1") assign("Rpropre", TRUE, envir=env)
        else assign("Rpropre", FALSE, envir=env)
        
        if(tclvalue(colValue)=="1") assign("Rcol", TRUE, envir=env)
        else assign("Rcol", FALSE, envir=env)
        
        if(tclvalue(colsupValue)=="1") assign("Rcolsup", TRUE, envir=env)
        else assign("Rcolsup", FALSE, envir=env)
        
        if(tclvalue(rowValue)=="1") assign("Rrow", TRUE, envir=env)
        else assign("Rrow", FALSE, envir=env)
        
        if(tclvalue(rowsupValue)=="1") assign("Rrowsup", TRUE, envir=env)
        else assign("Rrowsup", FALSE, envir=env)
        
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

      collab<-tklabel(SortieWin,text=.Facto_gettext("Results for column variables"))
        col.check <- tkcheckbutton(SortieWin)
        if(Rcol) colValue <- tclVar("1")
        else colValue <- tclVar("0")
        tkconfigure(col.check,variable=colValue)

      colsup.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for column supplementary variables"))
        colsup.check <- tkcheckbutton(SortieWin)

        if(Rcolsup) colsupValue <- tclVar("1")
        else colsupValue <- tclVar("0")
        tkconfigure(colsup.check,variable=colsupValue)

      row.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for row variables"))
        row.check <- tkcheckbutton(SortieWin)
        if(Rrow) rowValue <- tclVar("1")
        else rowValue <- tclVar("0")
        tkconfigure(row.check,variable=rowValue)


      rowsup.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for row supplementary variables"))
        rowsup.check <- tkcheckbutton(SortieWin)
        if(Rrowsup) rowsupValue <- tclVar("1")
        else rowsupValue <- tclVar("0")
        tkconfigure(rowsup.check,variable=rowsupValue)
        
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

        tkgrid(tklabel(SortieWin, text = .Facto_gettext("Select the outputs options"), fg ="blue"),  columnspan = 2, sticky = "w")
        tkgrid(tklabel(SortieWin, text = " "))
        tkgrid(eig.lab,eig.check,sticky="w")
        tkgrid(collab,col.check,sticky="w")
        if (!is.null(variableColIllu)) tkgrid(colsup.lab,colsup.check,sticky="w")
        tkgrid(row.lab,row.check,sticky="w")
        if (!is.null(variableLigneIllu)) tkgrid(rowsup.lab,rowsup.check,sticky="w")
        tkgrid(tklabel(SortieWin, text = " "))
        tkgrid(descdim.lab,descdim.check,sticky="w")
        tkgrid(RFichierFrame)
        tkgrid(SortieOK.but)
   }
    
    SortieFrame<-tkframe(IlluFrame)
    Sortie.but<-tkbutton(SortieFrame, textvariable=.SortieLabel, command=OnSortie, borderwidth=3, width=35)
    tkgrid(Sortie.but, sticky="ew")
  })


  #! fonction HCPC
  
  Hcpc.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    .HcpcLabel<-tclVar(paste(firstLabel, "", sep=" "))    
    compteur.hcpc<-0
    Rclassif<-0
    RclassifCA<-"rows"
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
      
        if(tclvalue(classifCAValue)=="rows") assign("RclassifCA", "rows", envir=env)
        else assign("RclassifCA", "columns", envir=env)
      
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

      lignes <- tkradiobutton (HcpcWin)
      lignes.lab <- tklabel(HcpcWin,text=.Facto_gettext("rows"))
      colonnes <- tkradiobutton (HcpcWin)
      colonnes.lab <- tklabel(HcpcWin,text=.Facto_gettext("columns"))
      classifCAValue <- tclVar(paste(RclassifCA))
      classifCA.lab <- tklabel(HcpcWin,text=.Facto_gettext("Perform clustering on: "))
      tkconfigure(lignes,variable=classifCAValue,value="rows")
      tkconfigure(colonnes,variable=classifCAValue,value="columns")

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
#      tkgrid(tklabel(HcpcWin,text=.Facto_gettext(paste('Clustering is performed on the first ', tclvalue(ncp.val), ' dimensions of the CA',sep=""))),column=1,columnspan=4,sticky="w")
      tkgrid(tklabel(HcpcWin,text=sprintf(.Facto_gettext("Clustering is performed on the first %s dimensions of CA"),tclvalue(ncp.val))),column=1,columnspan=4,sticky="w")   #text which takes the nb of dimensions chosen in the main window
      tkgrid(tklabel(HcpcWin,text=.Facto_gettext("(Modify in the main options to change this number)")),column=1,columnspan=4,sticky="w")
      tkgrid(tklabel(HcpcWin,text=""))
      tkgrid(classifCA.lab,lignes.lab,lignes)  
      tkgrid(colonnes.lab,colonnes)
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
      
      tkgrid.configure(minmaxhcpc.label,classifCA.lab,meth.lab,consolid.lab,graphhcpc.lab,reshcpc.lab,column=1,columnspan=4,sticky="w")
      tkgrid.configure(minhcpc.entry,column=7,columnspan=1,sticky="e")
      tkgrid.configure(maxhcpc.entry,column=8,columnspan=1,sticky="w")
      tkgrid.configure(lignes,colonnes,meth1,meth2,consolid.check,graphhcpc.check,reshcpc.check,column=8,sticky="e")
      tkgrid.configure(meth1.lab,lignes.lab,column=6,columnspan=2,sticky="w")
      tkgrid.configure(meth2.lab,colonnes.lab,column=6,columnspan=2,sticky="w") 
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

  
  #! fonction pour la réinitialisation des paramètre
  Reinitializ.funct<-function()
  {
    tkdestroy(top)
    FactoCA()
  }

  
#! fonction associée au bouton Appliquer, execute sans détruire l'interface graphique
  OnAppliquer<-function()
  {

        # liste de toutes les variables interne créées      (** mise en forme incomplète)
      # sur la fenetre principale
#         listColonne      **
#         listLigne        **
#         resu.val         **
#         ncp.val          **
#         Axe1
#         Axe2
      # dans les boutons des fenêtres illustratives
#         variableColIllu    **
#         variableLigneIllu  **
      # dans le bouton Plot PCA
#         Rchoix
#         RTitle
#         Rlabel
#         Rcol.ind
#         Rcol.ind.sup
#         Rcol.quali
#         Rinvisible       
#         RXlimInd
#         RYlimInd
      # dans le bouton affichage sortie
#         Rpropre
#             Rcol
#           Rcolsup
#           Rrow
#           Rrowsup


    # récupération des paramètres de la fenêtre principale
    nom.res<-tclvalue(resu.val)
    if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) justDoIt(paste('remove (',nom.res,')'))
    if(length(as.numeric(tkcurselection(listColonne)))<2) colActives<-listColonne.nom
    else colActives<-listColonne.nom[as.numeric(tkcurselection(listColonne))+1]
    colActives <- colActives[!(colActives%in%variableColIllu)]
    
    if(length(as.numeric(tkcurselection(listLigne)))<2) rowActives<-listLigne.nom
    else rowActives<-listLigne.nom[as.numeric(tkcurselection(listLigne))+1]
    rowActives <- rowActives[!(rowActives%in%variableLigneIllu)]
    
    Axe<-c(as.numeric(tclvalue(Axe1)), as.numeric(tclvalue(Axe2)))
    
    # gestion du tableau de données pour l'AFC
    if(!is.null(variableColIllu)) {
      if(!is.null(variableLigneIllu)) commande.data<-paste(activeDataSet(),'.CA', '<-', activeDataSet(), '[c("', paste(rowActives, collapse='", "'), '", "', paste(variableLigneIllu, collapse='", "'), '") ,c("', paste(colActives, collapse='", "'), '", "', paste(variableColIllu, collapse='", "'), '")]', sep="")
      else  commande.data<-paste(activeDataSet(),'.CA', '<-', activeDataSet(), '[c("', paste(rowActives, collapse='", "'), '") ,c("', paste(colActives, collapse='", "'), '", "', paste(variableColIllu, collapse='", "'), '")]', sep="")
    }
    else {
      if(!is.null(variableLigneIllu)) commande.data<-paste(activeDataSet(),'.CA', '<-', activeDataSet(), '[c("', paste(rowActives, collapse='", "'), '", "', paste(variableLigneIllu, collapse='", "'), '") ,c("', paste(colActives, collapse='", "'), '")]', sep="")
      else commande.data<-paste(activeDataSet(),'.CA', '<-', activeDataSet(), '[c("', paste(rowActives, collapse='", "'), '") ,c("', paste(colActives, collapse='", "'), '")]', sep="")
    }
    justDoIt(commande.data)
    logger(commande.data)
    donnee.depart<-activeDataSet()
    activeDataSet(paste(activeDataSet(),'.', 'CA', sep=""))
    
    # gestion de la commande réalisant l'AFC
    ncp<-as.numeric(tclvalue(ncp.val))
   
    if(!is.null(variableColIllu)) {
      if(!is.null(variableLigneIllu)) commande.ca<-paste(nom.res, '<-CA(', activeDataSet(), ', ncp=', ncp, ', row.sup=c(', nrow(get(getRcmdr(".activeDataSet")))-length(variableLigneIllu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))), '), col.sup=c(', ncol(get(getRcmdr(".activeDataSet")))-length(variableColIllu)+1, ': ', ncol(get(getRcmdr(".activeDataSet"))), '), graph = FALSE)', sep="")
      else commande.ca<-paste(nom.res, '<-CA(', activeDataSet(), ', ncp=', ncp, ', row.sup=NULL, col.sup=c(', ncol(get(getRcmdr(".activeDataSet")))-length(variableColIllu)+1, ': ', ncol(get(getRcmdr(".activeDataSet"))), '), graph = FALSE)', sep="")
    }
    else {
      if(!is.null(variableLigneIllu)) commande.ca<-paste(nom.res, '<-CA(', activeDataSet(), ', ncp=', ncp, ', row.sup=c(', nrow(get(getRcmdr(".activeDataSet")))-length(variableLigneIllu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))), '), col.sup=NULL, graph = FALSE)', sep="")
      else commande.ca<-paste(nom.res, '<-CA(', activeDataSet(), ', ncp=', ncp, ', row.sup=NULL, col.sup=NULL, graph = FALSE)', sep="")
    }
    justDoIt(commande.ca)
    logger(commande.ca)
	justDoIt(paste(nom.res,'$call$call <-',deparse(commande.ca),sep=""))

    # gestion des graphiques
    if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) {if (get(nom.res)$eig[1,2]==100) doItAndPrint(paste('"No graph can be plot: data are unidimensional"'))}
    if((Rchoix)&(length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0)){
    if (get(nom.res)$eig[1,2]!=100) {
      if (is.null(Rell)) commande.plot<-paste('plot.CA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), col.row="', Rcol.row, '", col.col="', Rcol.col, '", label=c("', paste(Rlabel, collapse='", "'), '")', sep="")
	  else commande.plot<-paste('ellipseCA(', nom.res,', ellipse=c("', paste(Rell, collapse='", "'), '"), axes=c(', paste(Axe, collapse=", "), '), col.row="', Rcol.row, '", col.col="', Rcol.col, '", label=c("', paste(Rlabel, collapse='", "'), '")', sep="")
      if (!is.null(RXlimInd)) commande.plot<-paste(commande.plot, ', xlim=c(', paste(RXlimInd, collapse=", "), ')')
      if (!is.null(RYlimInd)) commande.plot<-paste(commande.plot, ', ylim=c(', paste(RYlimInd, collapse=", "), ')')
      if (!is.null(Rinvis)&&(Rinvis!="")) commande.plot<-paste(commande.plot, ', invisible=c("', paste(Rinvis, collapse='", "'), '")',sep='')
      if(!is.null(variableLigneIllu)) commande.plot<-paste(commande.plot, ', col.row.sup="', Rcol.row.sup, '"',sep='')
      if(!is.null(variableColIllu))  commande.plot<-paste(commande.plot, ', col.col.sup="', Rcol.col.sup, '"',sep='')
      if (is.null(RTitle)) commande.plot <- paste(commande.plot,')', sep="")
      else {
        if (RTitle ==" ") commande.plot <- paste(commande.plot,')', sep="")
        else commande.plot <- paste(commande.plot,', title="', RTitle,'")', sep="")
      }
      justDoIt(commande.plot)
      logger(commande.plot)
    }}

      #Commande de la fonction HCPC

    if(Rclassif==TRUE){
      commande.hcpc<-paste(nom.res,'.hcpc', '<-HCPC(', nom.res, ' ,nb.clust=', Rmeth, ',consol=', Rconsolid,',min=', Rminhcpc,',max=',Rmaxhcpc,',cluster.CA="',RclassifCA,'",graph=', Rgraphhcpc, ')', sep="")
    justDoIt(commande.hcpc)
    logger(commande.hcpc)      
      if(Rreshcpc==TRUE){
        doItAndPrint(paste(nom.res,'.hcpc$data.clust[,ncol(res.hcpc$data.clust),drop=F]', sep=""))
        doItAndPrint(paste(nom.res,'.hcpc$desc.var', sep=""))
        doItAndPrint(paste(nom.res,'.hcpc$desc.axes', sep=""))
        doItAndPrint(paste(nom.res,'.hcpc$desc.ind', sep=""))
      }        
    }    
    
    # gestion de l'édition de certains resultats
    doItAndPrint(paste('summary(',nom.res,', nb.dec = 3, nbelements=10, nbind = 10, ncp = 3, file="")', sep=""))
    if (RFichier==""){
      if(Rpropre) doItAndPrint(paste(nom.res, '$eig', sep=""))
      if(Rcol) doItAndPrint(paste(nom.res, '$col', sep=""))
      if(Rcolsup & !is.null(variableColIllu))doItAndPrint(paste(nom.res, '$col.sup', sep=""))
      if(Rrow) doItAndPrint(paste(nom.res, '$row', sep=""))
      if(Rrowsup & !is.null(variableLigneIllu)) doItAndPrint(paste(nom.res, '$row.sup', sep=""))
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
      if(Rcol){
        doItAndPrint(paste('write.infile(', nom.res, '$col, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rcolsup){
        doItAndPrint(paste('write.infile(', nom.res, '$col.sup, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rrow){
        doItAndPrint(paste('write.infile(', nom.res, '$row, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rrowsup){
        doItAndPrint(paste('write.infile(', nom.res, '$row.sup, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rdescdim) doItAndPrint(paste('write.infile(dimdesc(', nom.res, ', axes=1:',ncp,'), file =',Fich,',append=',append,')', sep=""))
    }

    # Re-chargement du tableau de départ
    activeDataSet(donnee.depart)
    justDoIt(paste('remove(',activeDataSet(),'.CA)',sep=""))
    logger(paste('remove(',activeDataSet(),'.CA)',sep=""))    
  }
  
  #! fonction associée au bouton OK, execute et détruit l'interface graphique
  onOK<-function()
  {
    OnAppliquer()   
    # destuction de la fenêtre Top
    tkdestroy(top)
  
  }


################################################################################
#                   Création de la fenêtre top                                 #
################################################################################
  top<-tktoplevel(borderwidth=10)
  tkwm.title(top,.Facto_gettext("CA"))
  tkwm.geometry(top, "-100+50")
        
  # définition des polices
  font2<-tkfont.create(family="times",size=12,weight="bold")
  fontheading<-tkfont.create(family="times",size=18,weight="bold")

  # récupération du jeu de données actif
  donnee<-get(getRcmdr(".activeDataSet"))
  vars<-colnames(donnee)
  rows<-rownames(donnee)

  # création du frame contenant les listes colonnes et lignes
  ListeFrame<- tkframe(top, borderwidth=2)
    # liste des variables colonnes
  listColonne<-tklistbox(ListeFrame,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrColonne,...))
  scrColonne <-tkscrollbar(ListeFrame,repeatinterval=5,command=function(...)tkyview(listColonne,...))
  listColonne.nom<-NULL
  for (i in (1:ncol(donnee))) {
      if (is.numeric(donnee[,i])) {
          tkinsert(listColonne,"end",vars[i])
          listColonne.nom<-c(listColonne.nom, vars[i])
      }
  }
    # liste des variables lignes
  listLigne<-tklistbox(ListeFrame,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrLigne,...))
  scrLigne <-tkscrollbar(ListeFrame,repeatinterval=5,command=function(...)tkyview(listLigne,...))
  listLigne.nom<-NULL
  for (i in (1:nrow(donnee))) {
      tkinsert(listLigne,"end",rows[i])
        listLigne.nom<-c(listLigne.nom, rows[i])
  }
    # mise en forme de ListeFrame
  tkgrid(tklabel(ListeFrame, text = .Facto_gettext("Select the active rows and the active columns. \nBy default all rows and all columns are active"),fg = "darkred"), columnspan=5, sticky = "ew")
  tkgrid(listLigne, scrLigne, listColonne, scrColonne)
  tkgrid.configure(scrColonne, column=4, sticky="wns")
  tkgrid.configure(scrLigne, column=1, sticky="wns")
  tkgrid.configure(listLigne, sticky = "ew", column=0, columnspan=1)
  tkgrid.configure(listColonne, sticky = "ew", column=3, columnspan=1)
  tkgrid.columnconfigure(ListeFrame,2, minsize=95)
  
  
  # création de tous les boutons d'options dans IlluFrame
  IlluFrame<- tkframe(top, borderwidth=2)
  Cillu.funct(label=.Facto_gettext("Supplementary columns"), firstLabel=.Facto_gettext("Supplementary columns"))
  Lillu.funct(label=.Facto_gettext("Supplementary rows"), firstLabel=.Facto_gettext("Supplementary rows"))    
  PLOT.CA(label=.Facto_gettext("Graphical options"), firstLabel=.Facto_gettext("Graphical options"))
  Sortie.funct(label=.Facto_gettext("Outputs"), firstLabel=.Facto_gettext("Outputs"))
  tkgrid(LilluFrame, CilluFrame)
  tkgrid(tklabel(IlluFrame, text=""))
  tkgrid(PlotFrame, SortieFrame)
  tkgrid(tklabel(IlluFrame, text=""))
  tkgrid.configure(LilluFrame, PlotFrame, column=1, columnspan=2, sticky="ew")
  tkgrid.configure(CilluFrame, SortieFrame, column=4, columnspan=2, sticky="ew")
  tkgrid.columnconfigure(IlluFrame,0, minsize=25)
  tkgrid.columnconfigure(IlluFrame,6, minsize=25)
  tkgrid.columnconfigure(IlluFrame,3, minsize=35)  

  #Frame pour HCPC
  HcpcFrame<-tkframe(top, borderwidth=2)
  Hcpc.funct(label=.Facto_gettext("Perform Clustering after CA"), firstLabel=.Facto_gettext("Perform Clustering after CA")) 
  tkgrid(Hcpc2Frame, columnspan=7)
  tkgrid.configure(Hcpc2Frame,column=4, columnspan=1)

  # création des options dans OptionFrame  
  OptionFrame<-tkframe(top, borderwidth=2, relief="groove")
  resu.lab<-tklabel(OptionFrame,text=.Facto_gettext("Name of the result object: "))
  resu.val<-tclVar("res")
  resu<-tkentry(OptionFrame,width=10,textvariable=resu.val)
  ncp.lab<-tklabel(OptionFrame,text=.Facto_gettext("Number of dimensions: "))
  ncp.val<-tclVar("5") 
  ncp<-tkentry(OptionFrame,width=5,textvariable=ncp.val)
  Axe.label<-tklabel(OptionFrame,text=.Facto_gettext("Graphical output: select the dimensions"))
  Axe1<-tclVar("1")
  Axe2<-tclVar("2")
  Axe1.entry <-tkentry(OptionFrame,width="5",textvariable=Axe1)
  Axe2.entry <-tkentry(OptionFrame,width="5",textvariable=Axe2)
  
  
    # mise en page de OptionFrame
  tkgrid(tklabel(OptionFrame,text=.Facto_gettext("Main options"), fg = "darkred"), columnspan=8, sticky="we") 
  tkgrid(tklabel(OptionFrame,text="")) # Ligne de blanc
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

  appliquer.but<-tkbutton(top, text=.Facto_gettext("Apply"),width=12,command=OnAppliquer, borderwidth=3, fg="#690f96")
  OKCancelHelp(helpSubject="CA",reset="Reinitializ.funct")

  # Mise en page de top
  tkgrid(tklabel(top, text=.Facto_gettext("Correspondence Analysis (CA)"),font=fontheading), columnspan=3)
  tkgrid(tklabel(top,text=""))
  tkgrid(ListeFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text=""))
  tkgrid(IlluFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text=""))
  tkgrid(OptionFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text="")) # Ligne de blanc
  tkgrid(HcpcFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text=""))  
  # tkgrid(appliquer.but, column=1, columnspan=1)
  # tkgrid(tklabel(top,text="")) # Ligne de blanc
  # tkgrid(buttonsFrame, column=1, columnspan=1, sticky="ew" )
  # tkgrid(tklabel(top,text="")) # Ligne de blanc
  tkgrid(buttonsFrame, appliquer.but)
  tkgrid.configure(buttonsFrame, column=1,sticky="e")
  tkgrid.configure(appliquer.but, column=2,sticky="w")
  
}
