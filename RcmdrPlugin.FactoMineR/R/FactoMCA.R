FactoMCA <-
function()
{

#    Création des fonctions pour les options via nouvelle fenêtre graphique


  #! fonction pour le choix des variables qualitatives supplémentaires
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
        if(length(fact.select)==0)
        {
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
   Fillu.but<-tkbutton(FilluFrame, textvariable=.FilluLabel, command=OnFillu, borderwidth=3)
   tkgrid(Fillu.but, sticky="ew")
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
          if(rows[i] %in% individuillu) tkselection.set(listind,indice)
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


  #! fonction pour la gestion des options graphiques
  PLOT.MCA<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.graph<-0
    #déclaration des variables
    Rchoix<-TRUE
    RTitle<-NULL
    Rlabel<-c("ind", "ind.sup", "quali.sup","var")
    Rinvis<- ""
    Rcol.ind<-Rcol.ind.tmp<-"black"
    Rcol.ind.sup<-Rcol.ind.sup.tmp<-"blue"
    Rcol.quali<-Rcol.quali.tmp<-"darkred"
    Rcol.qualisup<-Rcol.qualisup.tmp<-"darkgreen"
    RXlimInd<-NULL
    RYlimInd<-NULL

    Wchoix=TRUE
    WTitle<-NULL
    Wlabel<-c("quanti.sup")
    #Wlim.cos<-0.
    Wcol.quanti.sup<-Wcol.quanti.sup.tmp<-"blue"
    Wcol.var<-Wcol.var.tmp<-"black"
    WXlimVar<-NULL
    WYlimVar<-NULL

    Vchoix=TRUE
    VTitle<-NULL
    Vlabel<-c("var","quali.sup")
    Vinvis<- ""
    Vcol.quali.sup<-Vcol.quali.sup.tmp<-"darkgreen"
    Vcol.var<-Vcol.var.tmp<-"darkred"
    VXlimVar<-NULL
    VYlimVar<-NULL

    .PlotLabel<-tclVar(paste(firstLabel, "", sep=" "))

    OnPlot<-function()
    {
      PlotWin<-tktoplevel()
      tkwm.title(PlotWin, .Facto_gettext("Select graphical options"))
      tkwm.geometry(PlotWin, "-100+50")

      #création de la fonction onOKsub
      onOKsub<-function()
      {
        assign("compteur.graph", compteur.graph+1, envir=env)
        if(compteur.graph>0) tclvalue(.PlotLabel)<-paste(label, "", sep=" ")
        tkconfigure(Plot.but, fg="blue")

        # gestion des entrées de la partie graphique des individus
        if(tclvalue(ind.check.value)==1) assign("Rchoix", TRUE, envir=env)
        else assign("Rchoix", FALSE, envir=env)

        if(Rchoix)
        {
          if (tclvalue(Titre)==" ") assign("RTitle", NULL, envir=env)
          assign("RTitle", tclvalue(Titre), envir=env)

          label.tmp.ind<-tclvalue(label.ind.checkValue)
          label.tmp.ind.sup<-tclvalue(label.ind.sup.checkValue)
          label.tmp.quali.sup<-tclvalue(label.quali.sup.checkValue)
          label.tmp.var<-tclvalue(label.var.checkValue)
          assign("Rlabel", NULL, envir=env)
          if(label.tmp.ind==1) assign("Rlabel", c(Rlabel, "ind"), envir=env)
          if(label.tmp.ind.sup==1) assign("Rlabel", c(Rlabel, "ind.sup"), envir=env)
          if(label.tmp.quali.sup==1) assign("Rlabel", c(Rlabel, "quali.sup"), envir=env)
          if(label.tmp.var==1) assign("Rlabel", c(Rlabel, "var"), envir=env)

          invis.tmp.ind<-tclvalue(invis.ind.checkValue)
          invis.tmp.ind.sup<-tclvalue(invis.ind.sup.checkValue)
          invis.tmp.quali.sup<-tclvalue(invis.quali.sup.checkValue)
          invis.tmp.var<-tclvalue(invis.var.checkValue)
          assign("Rinvis", NULL, envir=env)
          if(invis.tmp.ind==0) assign("Rinvis", c(Rinvis, "ind"), envir=env)
          if(invis.tmp.ind.sup==0) assign("Rinvis", c(Rinvis, "ind.sup"), envir=env)
          if(invis.tmp.quali.sup==0) assign("Rinvis", c(Rinvis, "quali.sup"), envir=env)
          if(invis.tmp.var==0) assign("Rinvis", c(Rinvis, "var"), envir=env)

          assign("Rcol.ind", Rcol.ind.tmp, envir=env)
          assign("Rcol.ind.sup", Rcol.ind.sup.tmp, envir=env)
          assign("Rcol.quali", Rcol.quali.tmp, envir=env)
          assign("Rcol.qualisup", Rcol.qualisup.tmp, envir=env)
          assign("Wcol.var", Wcol.var.tmp, envir=env)

          if(tclvalue(XlimIndMin)=="" | tclvalue(XlimIndMax)=="") assign("RXlimInd", NULL, envir=env)
          else assign("RXlimInd", c(as.numeric(tclvalue(XlimIndMin)), as.numeric(tclvalue(XlimIndMax))), envir=env)
          if(tclvalue(YlimIndMin)=="" | tclvalue(YlimIndMax)=="") assign("RYlimInd", NULL, envir=env)
          else assign("RYlimInd", c(as.numeric(tclvalue(YlimIndMin)), as.numeric(tclvalue(YlimIndMax))), envir=env)
        }

        if(tclvalue(quanti.var.check.value)==1) assign("Wchoix", TRUE, envir=env)
        else assign("Wchoix", FALSE, envir=env)

        if(Wchoix)
        {
          label.tmp.quanti.sup<-tclvalue(label.quanti.sup.checkValue)
          if(label.tmp.quanti.sup==1) assign("Wlabel", c(Wlabel, "quanti.sup"), envir=env)

          if (tclvalue(WTitre)==" ") assign("WTitle", NULL, envir=env)
          assign("WTitle", tclvalue(WTitre), envir=env)
          #assign("Wlim.cos", tclvalue(WlimCosValue), envir=env)
          label.tmp.quanti.sup<-tclvalue(label.quanti.sup.checkValue)
          assign("Wlabel", NULL, envir=env)
          if(label.tmp.quanti.sup==1) assign("Wlabel", c(Wlabel, "quanti.sup"), envir=env)
          assign("Wcol.quanti.sup", Wcol.quanti.sup.tmp, envir=env)
               
        }

        if(tclvalue(var.check.value)==1) assign("Vchoix", TRUE, envir=env)
        else assign("Vchoix", FALSE, envir=env)

        if(Vchoix)
        {
          if (tclvalue(VTitre)==" ") assign("VTitle", NULL, envir=env)
          assign("VTitle", tclvalue(VTitre), envir=env)

          labelV.tmp.quali.sup<-tclvalue(labelV.quali.sup.checkValue)
          labelV.tmp.var<-tclvalue(labelV.var.checkValue)
          assign("Vlabel", NULL, envir=env)
          if(labelV.tmp.quali.sup==1) assign("Vlabel", c(Vlabel, "quali.sup"), envir=env)
          if(labelV.tmp.var==1) assign("Vlabel", c(Vlabel, "var"), envir=env)

          invisV.tmp.quali.sup<-tclvalue(invisV.quali.sup.checkValue)
          invisV.tmp.var<-tclvalue(invisV.var.checkValue)
          assign("Vinvis", NULL, envir=env)
          if(invisV.tmp.quali.sup==0) assign("Vinvis", c(Vinvis, "quali.sup"), envir=env)
          if(invisV.tmp.var==0) assign("Vinvis", c(Vinvis, "var"), envir=env)

          assign("Vcol.quali.sup", Vcol.quali.sup.tmp, envir=env)
          assign("Vcol.var", Vcol.var.tmp, envir=env)
               
        }
        tkdestroy(PlotWin)
      }

    # construction de la partie graphique des individus
    PlotIndFrame<-tkframe(PlotWin, borderwidth=5, relief="groove")

    RchoixFrame<-tkframe(PlotIndFrame,borderwidth=2)
    ind.check<-tkcheckbutton(RchoixFrame)
    if(Rchoix) ind.check.value<-tclVar("1")
    else ind.check.value<-tclVar("0")
    tkconfigure(ind.check, variable=ind.check.value)
    tkgrid(tklabel(RchoixFrame, text=.Facto_gettext("Graphical output"), font=font2),ind.check)
    tkgrid(tklabel(RchoixFrame, text="  "))

    RTitleFrame<-tkframe(PlotIndFrame,borderwidth=2)
    if (is.null(RTitle)) Titre <- tclVar(" ")
    else Titre<-tclVar(RTitle)
    Titre.entry <-tkentry(RTitleFrame,width="40",textvariable=Titre)
    tkgrid(tklabel(RTitleFrame,text=.Facto_gettext("Title of the graph")),Titre.entry)

    RlabelFrame<-tkframe(PlotIndFrame,borderwidth=2)
    tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("  ")),tklabel(RlabelFrame,text=.Facto_gettext("Plot")),tklabel(RlabelFrame, text=.Facto_gettext("Label")))
    label.ind.check<-tkcheckbutton(RlabelFrame)
    if ("ind" %in% Rlabel) label.ind.checkValue<-tclVar("1")
    else label.ind.checkValue<-tclVar("0")
    invis.ind.check<-tkcheckbutton(RlabelFrame)
    if ("ind" %in% Rinvis) invis.ind.checkValue<-tclVar("0")
    else invis.ind.checkValue<-tclVar("1")
    tkconfigure(label.ind.check, variable=label.ind.checkValue)
    tkconfigure(invis.ind.check, variable=invis.ind.checkValue)
    tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Active individuals")), invis.ind.check,label.ind.check)
    label.ind.sup.check<-tkcheckbutton(RlabelFrame)
    if ("ind.sup" %in% Rlabel) label.ind.sup.checkValue<-tclVar("1")
    else label.ind.sup.checkValue<-tclVar("0")
    invis.ind.sup.check<-tkcheckbutton(RlabelFrame)
    if ("ind.sup" %in% Rinvis) invis.ind.sup.checkValue<-tclVar("0")
    else invis.ind.sup.checkValue<-tclVar("1")
    tkconfigure(label.ind.sup.check, variable=label.ind.sup.checkValue)
    tkconfigure(invis.ind.sup.check, variable=invis.ind.sup.checkValue)
    if(!is.null(individuillu)) tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Supplementary individuals")),invis.ind.sup.check,label.ind.sup.check)
    label.var.check<-tkcheckbutton(RlabelFrame)
    if ("var" %in% Rlabel) label.var.checkValue<-tclVar("1")
    else label.var.checkValue<-tclVar("0")
    invis.var.check<-tkcheckbutton(RlabelFrame)
    if ("var" %in% Rinvis) invis.var.checkValue<-tclVar("0")
    else invis.var.checkValue<-tclVar("1")
    tkconfigure(label.var.check, variable=label.var.checkValue)
    tkconfigure(invis.var.check, variable=invis.var.checkValue)
    tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Factors")),invis.var.check, label.var.check)
    label.quali.sup.check<-tkcheckbutton(RlabelFrame)
    if ("quali.sup" %in% Rlabel) label.quali.sup.checkValue<-tclVar("1")
    else label.quali.sup.checkValue<-tclVar("0")
    invis.quali.sup.check<-tkcheckbutton(RlabelFrame)
    if ("quali.sup" %in% Rinvis) invis.quali.sup.checkValue<-tclVar("0")
    else invis.quali.sup.checkValue<-tclVar("1")
    tkconfigure(label.quali.sup.check, variable=label.quali.sup.checkValue)
    tkconfigure(invis.quali.sup.check, variable=invis.quali.sup.checkValue)
    if(!is.null(variablefact)) tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Supplementary factors")),invis.quali.sup.check, label.quali.sup.check)

    RcolFrame<-tkframe(PlotIndFrame,borderwidth=2)
    Rcol.ind.value <- Rcol.ind
    canvas.ind <- tkcanvas(RcolFrame,width="80",height="25",bg=Rcol.ind.value)
    ChangeColor.ind <- function()
    {
      Rcol.ind.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Rcol.ind.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Rcol.ind.value)>0)
      {
        tkconfigure(canvas.ind,bg=Rcol.ind.value)
        assign("Rcol.ind.tmp", Rcol.ind.value, envir=env)
      }
    }
    ChangeColor.ind.button <- tkbutton(RcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.ind)
    tkgrid(tklabel(RcolFrame, text=.Facto_gettext("Color for active individuals")),canvas.ind,ChangeColor.ind.button)

    Rcol.ind.sup.value<-Rcol.ind.sup
    canvas.ind.sup <- tkcanvas(RcolFrame,width="80",height="25",bg=Rcol.ind.sup.value)
    ChangeColor.ind.sup <- function()
    {
      Rcol.ind.sup.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Rcol.ind.sup.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Rcol.ind.sup.value)>0)
      {
        tkconfigure(canvas.ind.sup,bg=Rcol.ind.sup.value)
        assign("Rcol.ind.sup.tmp", Rcol.ind.sup.value, envir=env)
      }
    }
    ChangeColor.ind.sup.button <- tkbutton(RcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.ind.sup)
    if(!is.null(individuillu)) tkgrid(tklabel(RcolFrame, text=.Facto_gettext("color for supplementary individuals")),canvas.ind.sup,ChangeColor.ind.sup.button)

    Rcol.quali.value<- Rcol.quali
    canvas.quali <- tkcanvas(RcolFrame,width="80",height="25",bg=Rcol.quali.value)
    ChangeColor.quali <- function()
    {
      Rcol.quali.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Rcol.quali.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Rcol.quali.value)>0)
      {
        tkconfigure(canvas.quali,bg=Rcol.quali.value)
        assign("Rcol.quali.tmp", Rcol.quali.value, envir=env)
      }
    }
    ChangeColor.quali.button <- tkbutton(RcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.quali)
    tkgrid(tklabel(RcolFrame, text=.Facto_gettext("Color for factors")),canvas.quali,ChangeColor.quali.button)

    Rcol.qualisup.value<- Rcol.qualisup
    canvas.qualisup <- tkcanvas(RcolFrame,width="80",height="25",bg=Rcol.qualisup.value)
    ChangeColor.qualisup <- function()
    {
      Rcol.qualisup.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Rcol.qualisup.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Rcol.qualisup.value)>0)
      {
        tkconfigure(canvas.qualisup,bg=Rcol.qualisup.value)
        assign("Rcol.qualisup.tmp", Rcol.qualisup.value, envir=env)
      }
    }
    ChangeColor.qualisup.button <- tkbutton(RcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.qualisup)
    if(!is.null(variablefact)) tkgrid(tklabel(RcolFrame, text=.Facto_gettext("Color for supplementary factors")),canvas.qualisup,ChangeColor.qualisup.button)

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
    #tkgrid(tklabel(PlotIndFrame, text=.Facto_gettext("Individuals graph"), font=font2))
    #tkgrid(tklabel(PlotIndFrame, text=" "))
    tkgrid(RchoixFrame)
    tkgrid(RTitleFrame)
    tkgrid(tklabel(PlotIndFrame, text=" "))
    tkgrid(RlabelFrame)
    tkgrid(tklabel(PlotIndFrame, text=" "))
    tkgrid(RcolFrame)
    tkgrid(tklabel(PlotIndFrame, text=" "))
    tkgrid(RlimFrame)
    tkgrid(tklabel(PlotIndFrame, text=" "))

    PlotVarFrame<-tkframe(PlotWin, borderwidth=5, relief="groove")
    WchoixFrame<-tkframe(PlotVarFrame,borderwidth=2)
    quanti.var.check<-tkcheckbutton(WchoixFrame)
    if(Wchoix) quanti.var.check.value<-tclVar("1")
    else quanti.var.check.value<-tclVar("0")
    tkconfigure(quanti.var.check, variable=quanti.var.check.value)
    tkgrid(tklabel(WchoixFrame, text=.Facto_gettext("Plot supplementary variables graph"), font=font2),quanti.var.check)
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
    label.quanti.sup.check<-tkcheckbutton(WlabelFrame)
    if ("quanti.sup" %in% Wlabel) label.quanti.sup.checkValue<-tclVar("1")
    else label.quanti.sup.checkValue<-tclVar("0")
    tkconfigure(label.quanti.sup.check, variable=label.quanti.sup.checkValue)
    tkgrid(tklabel(WlabelFrame, text=.Facto_gettext("Labels for the supplementary quantitative variables")),label.quanti.sup.check)

    WcolFrame<-tkframe(PlotVarFrame,borderwidth=2)
    Wcol.quanti.sup.value<-Wcol.quanti.sup
    canvas.quanti.sup <- tkcanvas(WcolFrame,width="80",height="25",bg=Wcol.quanti.sup.value)
    ChangeColor.quanti.sup <- function()
    {
      Wcol.quanti.sup.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Wcol.quanti.sup.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Wcol.quanti.sup.value)>0) {
        tkconfigure(canvas.quanti.sup,bg=Wcol.quanti.sup.value)
        assign("Wcol.quanti.sup.tmp", Wcol.quanti.sup.value, envir=env)
      }
    }
    ChangeColor.quanti.sup.button <- tkbutton(WcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.quanti.sup)
    tkgrid(tklabel(WcolFrame, text=.Facto_gettext("Color for supplementary variables")),canvas.quanti.sup,ChangeColor.quanti.sup.button)

    #mise en page des différents frames de PlotVarFrame
    tkgrid(WchoixFrame)
    tkgrid(WTitleFrame)
    #tkgrid(WcosFrame)
    tkgrid(WlabelFrame)
    tkgrid(tklabel(PlotVarFrame, text=" "))
    tkgrid(WcolFrame)
    tkgrid(tklabel(PlotVarFrame, text=" "))
    
################Début à changer
    PlotVVarFrame<-tkframe(PlotWin, borderwidth=5, relief="groove")
    VchoixFrame<-tkframe(PlotVVarFrame,borderwidth=2)
    var.check<-tkcheckbutton(VchoixFrame)
    if(Vchoix) var.check.value<-tclVar("1")
    else var.check.value<-tclVar("0")
    tkconfigure(var.check, variable=var.check.value)
    tkgrid(tklabel(VchoixFrame, text=.Facto_gettext("Plot variables graph"), font=font2),var.check)
    tkgrid(tklabel(VchoixFrame, text="  "))

    VTitleFrame<-tkframe(PlotVVarFrame,borderwidth=2)
    if (is.null(VTitle)) VTitre <- tclVar(" ")
    else VTitre<-tclVar(VTitle)
    VTitre.entry <-tkentry(VTitleFrame,width="40",textvariable=VTitre)
    tkgrid(tklabel(VTitleFrame,text=.Facto_gettext("Title of the graph")),VTitre.entry)

    VlabelFrame<-tkframe(PlotVVarFrame,borderwidth=2)
    tkgrid(tklabel(VlabelFrame, text=.Facto_gettext("  ")),tklabel(VlabelFrame,text=.Facto_gettext("Plot")),tklabel(VlabelFrame, text=.Facto_gettext("Label")))
    labelV.var.check<-tkcheckbutton(VlabelFrame)
##    label.quali.sup.check<-tkcheckbutton(VlabelFrame)
##    if ("quali.sup" %in% Vlabel) label.quali.sup.checkValue<-tclVar("1")
##    else label.quali.sup.checkValue<-tclVar("0")
##    tkconfigure(label.quali.sup.check, variable=label.quali.sup.checkValue)
##    tkgrid(tklabel(VlabelFrame, text=.Facto_gettext("Labels for the supplementary variables")),label.quali.sup.check)
    if ("var" %in% Vlabel) labelV.var.checkValue<-tclVar("1")
    else labelV.var.checkValue<-tclVar("0")
    invisV.var.check<-tkcheckbutton(VlabelFrame)
    if ("var" %in% Vinvis) invisV.var.checkValue<-tclVar("0")
    else invisV.var.checkValue<-tclVar("1")
    tkconfigure(labelV.var.check, variable=labelV.var.checkValue)
    tkconfigure(invisV.var.check, variable=invisV.var.checkValue)
    tkgrid(tklabel(VlabelFrame, text=.Facto_gettext("Active variables")),invisV.var.check, labelV.var.check)
    labelV.quali.sup.check<-tkcheckbutton(VlabelFrame)
    if ("quali.sup" %in% Vlabel) labelV.quali.sup.checkValue<-tclVar("1")
    else labelV.quali.sup.checkValue<-tclVar("0")
    invisV.quali.sup.check<-tkcheckbutton(VlabelFrame)
    if ("quali.sup" %in% Vinvis) invisV.quali.sup.checkValue<-tclVar("0")
    else invisV.quali.sup.checkValue<-tclVar("1")
    tkconfigure(labelV.quali.sup.check, variable=labelV.quali.sup.checkValue)
    tkconfigure(invisV.quali.sup.check, variable=invisV.quali.sup.checkValue)
    if(!is.null(variablefact)) tkgrid(tklabel(VlabelFrame, text=.Facto_gettext("Supplementary factors")),invisV.quali.sup.check, labelV.quali.sup.check)

    VcolFrame<-tkframe(PlotVVarFrame,borderwidth=2)
    Vcol.var.value<-Vcol.var
    canvas.varV <- tkcanvas(VcolFrame,width="80",height="25",bg=Vcol.var.value)
    ChangeColor.varV <- function()
    {
      Vcol.var.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Vcol.var.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Vcol.var.value)>0) {
        tkconfigure(canvas.varV,bg=Vcol.var.value)
        assign("Vcol.var.tmp", Vcol.var.value, envir=env)
      }
    }
    ChangeColor.varV.button <- tkbutton(VcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.varV)
    tkgrid(tklabel(VcolFrame, text=.Facto_gettext("Color for active variables")),canvas.varV,ChangeColor.varV.button)
###
    Vcol.quali.sup.value<-Vcol.quali.sup
    canvas.quali.sup <- tkcanvas(VcolFrame,width="80",height="25",bg=Vcol.quali.sup.value)
    ChangeColor.quali.sup <- function()
    {
      Vcol.quali.sup.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Vcol.quali.sup.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Vcol.quali.sup.value)>0) {
        tkconfigure(canvas.quali.sup,bg=Vcol.quali.sup.value)
        assign("Vcol.quali.sup.tmp", Vcol.quali.sup.value, envir=env)
      }
    }
    ChangeColor.quali.sup.button <- tkbutton(VcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.quali.sup)
    tkgrid(tklabel(VcolFrame, text=.Facto_gettext("Color for supplementary variables")),canvas.quali.sup,ChangeColor.quali.sup.button)

    #mise en page des différents frames de PlotVarFrame
    tkgrid(VchoixFrame)
    tkgrid(VTitleFrame)
    tkgrid(VlabelFrame)
    tkgrid(tklabel(PlotVVarFrame, text=" "))
    tkgrid(VcolFrame)
    tkgrid(tklabel(PlotVVarFrame, text=" "))
    
################Fin à changer

    # construction de la partie graphique des variables

    subOKCancelHelp(PlotWin, "plot.MCA")
    tkgrid(PlotIndFrame,PlotVVarFrame)
    tkgrid.configure(PlotVVarFrame,sticky="ne")
    if (length(variableillu)>0) {
      tkgrid(PlotVarFrame)
      tkgrid.configure(PlotVarFrame, sticky="se")
    }
    tkgrid(subButtonsFrame, sticky="ew", columnspan=2)

    }
    PlotFrame<-tkframe(IlluFrame)
    Plot.but<-tkbutton(PlotFrame, textvariable=.PlotLabel, command=OnPlot, borderwidth=3)
    tkgrid(Plot.but, sticky="ew")
  })


  #! fonction pour la réinitialisation des paramètre
  Reinitializ.funct<-function()
  {
    tkdestroy(top)
        FactoMCA()
  }


  #! fonction pour le choix des éléments de sortie
  Sortie.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.sortie<-0
    #déclaration des variables
    Rpropre<-FALSE
    RFichier <- ""
        Rvariable<-FALSE
        Rindividu<-FALSE
        Rindsup<-FALSE
        Rquantisup<-FALSE
        Rqualisup<-FALSE
        Rdescdim<-FALSE

    .SortieLabel<-tclVar(paste(firstLabel, "", sep=" "))

    OnSortie<-function()
    {
      SortieWin<-tktoplevel()
      tkwm.title(SortieWin,"Displayed outputs")

      #création de la fonction onOKsub
      onOK.sortie<-function()
      {
        assign("compteur.sortie", compteur.sortie+1, envir=env)
        if(compteur.sortie>0) tclvalue(.SortieLabel)<-paste(label, "", sep=" ")
        tkconfigure(Sortie.but, fg="blue")

        if(tclvalue(eigValue)=="1") assign("Rpropre", TRUE, envir=env)
        else assign("Rpropre", FALSE, envir=env)

        if(tclvalue(varValue)=="1") assign("Rvariable", TRUE, envir=env)
        else assign("Rvariable", FALSE, envir=env)

        if(tclvalue(indValue)=="1") assign("Rindividu", TRUE, envir=env)
        else assign("Rindividu", FALSE, envir=env)

        if(tclvalue(ind.sup.Value)=="1") assign("Rindsup", TRUE, envir=env)
        else assign("Rindsup", FALSE, envir=env)

        if(tclvalue(quanti.sup.Value)=="1") assign("Rquantisup", TRUE, envir=env)
        else assign("Rquantisup", FALSE, envir=env)

        if(tclvalue(quali.sup.Value)=="1") assign("Rqualisup", TRUE, envir=env)
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

      var.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for active variables"))
        var.check <- tkcheckbutton(SortieWin)
        if(Rvariable) varValue <- tclVar("1")
        else varValue <- tclVar("0")
        tkconfigure(var.check,variable=varValue)

      ind.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for active individuals"))
        ind.check <- tkcheckbutton(SortieWin)
        if(Rindividu) indValue <- tclVar("1")
        else indValue <- tclVar("0")
        tkconfigure(ind.check,variable=indValue)

      ind.sup.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for supplementary individuals"))
        ind.sup.check <- tkcheckbutton(SortieWin)
        if(Rindsup) ind.sup.Value <- tclVar("1")
        else ind.sup.Value <- tclVar("0")
        tkconfigure(ind.sup.check,variable=ind.sup.Value)

      quanti.sup.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for supplementary variables"))
        quanti.sup.check <- tkcheckbutton(SortieWin)
        if(Rquantisup) quanti.sup.Value <- tclVar("1")
        else quanti.sup.Value <- tclVar("0")
        tkconfigure(quanti.sup.check,variable=quanti.sup.Value)

      quali.sup.lab<-tklabel(SortieWin,text=.Facto_gettext("Results for supplementary factors"))
        quali.sup.check <- tkcheckbutton(SortieWin)
      if(Rqualisup) quali.sup.Value <- tclVar("1")
        else quali.sup.Value <- tclVar("0")
        tkconfigure(quali.sup.check,variable=quali.sup.Value)

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
        tkgrid(var.lab,var.check,sticky="w")
        tkgrid(ind.lab,ind.check,sticky="w")
        if (!is.null(individuillu)) tkgrid(ind.sup.lab,ind.sup.check,sticky="w")
        if (!is.null(variableillu)) tkgrid(quanti.sup.lab,quanti.sup.check,sticky="w")
        if (!is.null(variablefact)) tkgrid(quali.sup.lab,quali.sup.check,sticky="w")
        tkgrid(descdim.lab,descdim.check,sticky="w")
        tkgrid(tklabel(SortieWin, text = " "))
        tkgrid(RFichierFrame)
        tkgrid(SortieOK.but)
   }

    SortieFrame<-tkframe(IlluFrame)
    Sortie.but<-tkbutton(SortieFrame, textvariable=.SortieLabel, command=OnSortie, borderwidth=3)
    tkgrid(Sortie.but, sticky="ew")
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
#      tkgrid(tklabel(HcpcWin,text=.Facto_gettext(paste('Clustering is performed on the first ', tclvalue(ncp.val), ' dimensions of the MCA',sep=""))),column=1,columnspan=4,sticky="w")
      tkgrid(tklabel(HcpcWin,text=sprintf(.Facto_gettext("Clustering is performed on the first %s dimensions of MCA"),tclvalue(ncp.val))),column=1,columnspan=4,sticky="w")   #text which takes the nb of dimensions chosen in the main window
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
onOK <- function(){
  done = OnAppliquer()
  tkdestroy(top)
}

  #! fonction associer au bouton Appliquer, execute sans détruire la fenêtre top
  OnAppliquer<-function()
  {
    # liste de toutes les variables interne créées      (** mise en forme incomplète)
      # sur la fenetre principale
#         listdesc         **
#         resu.val         **
#         ncp.val          **
      # dans les boutons des fenêtres illustratives
#         variablefact     **
#         variableillu     **
#         individuillu     **
      # dans le bouton Plot MCA
#         Rchoix
#         RTitle
#         Rlabel
#         Rcol.ind
#         Rcol.ind.sup
#         Rcol.quali
#         Rhabillage       **
#         RXlimInd
#         RYlimInd
#         Wcol.var
      # dans le bouton affichage sortie
#         Rpropre
#             Rvariable
#           Rindividu
#           Rindsup
#           Rquantisup
#           Rqualisup


    # récupération des paramètres de la fenêtre principale
    nom.res<-tclvalue(resu.val)
    if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) justDoIt(paste('remove (',nom.res,')'))
    if(length(as.numeric(tkcurselection(listdesc)))<2) varActives<-listdesc.nom
    else varActives<-listdesc.nom[as.numeric(tkcurselection(listdesc))+1]
    varActives <- varActives[!(varActives%in%variablefact)]
    ncp<-as.numeric(tclvalue(ncp.val))
    Axe<-c(as.numeric(tclvalue(Axe1)), as.numeric(tclvalue(Axe2)))

    # gestion du tableau de données pour l'ACM
    if(!is.null(individuillu)) {
      ind.actif<-rows[-which(rows %in% individuillu)]
      if(!is.null(variableillu)) {
        if(!is.null(variablefact)) commande.data<-paste(activeDataSet(),'.', 'MCA', '<-', activeDataSet(), '[c("', paste(ind.actif, collapse='", "'), '", "', paste(individuillu, collapse='", "'), '") ,c("', paste(varActives, collapse='", "'), '", "', paste(variableillu, collapse='", "'), '", "', paste(variablefact, collapse='", "'), '")]', sep="")
        else commande.data<-paste(activeDataSet(),'.', 'MCA', '<-', activeDataSet(), '[c("', paste(ind.actif, collapse='", "'), '", "', paste(individuillu, collapse='", "'), '") ,c("', paste(varActives, collapse='", "'), '", "', paste(variableillu, collapse='", "'), '")]', sep="")
      }
      else {
        if(!is.null(variablefact)) commande.data<-paste(activeDataSet(),'.', 'MCA', '<-', activeDataSet(), '[c("', paste(ind.actif, collapse='", "'), '", "', paste(individuillu, collapse='", "'), '") ,c("', paste(varActives, collapse='", "'), '", "', paste(variablefact, collapse='", "'), '")]', sep="")
        else commande.data<-paste(activeDataSet(),'.', 'MCA', '<-', activeDataSet(), '[c("', paste(ind.actif, collapse='", "'), '", "', paste(individuillu, collapse='", "'), '") ,c("', paste(varActives, collapse='", "'), '")]', sep="")
      }
    }
    else {
      if(!is.null(variableillu)) {
        if(!is.null(variablefact)) commande.data<-paste(activeDataSet(),'.', 'MCA', '<-', activeDataSet(), '[, c("', paste(varActives, collapse='", "'), '", "', paste(variableillu, collapse='", "'), '", "', paste(variablefact, collapse='", "'), '")]', sep="")
        else commande.data<-paste(activeDataSet(),'.', 'MCA', '<-', activeDataSet(), '[, c("', paste(varActives, collapse='", "'), '", "', paste(variableillu, collapse='", "'), '")]', sep="")
      }
      else {
        if(!is.null(variablefact)) commande.data<-paste(activeDataSet(),'.', 'MCA', '<-', activeDataSet(), '[, c("', paste(varActives, collapse='", "'), '", "', paste(variablefact, collapse='", "'), '")]', sep="")
        else commande.data<-paste(activeDataSet(),'.', 'MCA', '<-', activeDataSet(), '[, c("', paste(varActives, collapse='", "'), '")]', sep="")
      }
    }
    justDoIt(commande.data)
    logger(commande.data)
    donnee.depart<-activeDataSet()
    activeDataSet(paste(activeDataSet(),'.', 'MCA', sep=""))

    # gestion de la commande réalisant l'ACM
    if(!is.null(individuillu)) {
      ind.actif<-rows[-which(rows %in% individuillu)]
      if(!is.null(variableillu)) {
        if(!is.null(variablefact))    commande.acm<-paste(nom.res, '<-MCA(', activeDataSet(), ', ncp=', ncp, ', ind.sup=', nrow(get(getRcmdr(".activeDataSet")))-length(individuillu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))), ', quanti.sup=', length(varActives)+1, ': ', length(varActives)+ length(variableillu), ', quali.sup=', length(varActives)+length(variableillu)+1, ': ', length(varActives)+length(variableillu)+length(variablefact), ', graph = FALSE)', sep="")
        else commande.acm<-paste(nom.res, '<-MCA(', activeDataSet(), ', ncp=', ncp, ', ind.sup=', nrow(get(getRcmdr(".activeDataSet")))-length(individuillu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))), ', quanti.sup=', length(varActives)+1, ': ', length(varActives)+ length(variableillu), ', graph = FALSE)', sep="")
      }
      else {
        if(!is.null(variablefact))  commande.acm<-paste(nom.res, '<-MCA(', activeDataSet(), ', ncp=', ncp, ', ind.sup=', nrow(get(getRcmdr(".activeDataSet")))-length(individuillu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))), ', quali.sup=', length(varActives)+1, ': ', length(varActives)+length(variablefact), ', graph = FALSE)', sep="")
        else commande.acm<-paste(nom.res, '<-MCA(', activeDataSet(), ', ncp=', ncp, ', ind.sup=c(', nrow(get(getRcmdr(".activeDataSet")))-length(individuillu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))), '), graph = FALSE)', sep="")
      }
    }
    else
    {
      if(!is.null(variableillu)) {
        if(!is.null(variablefact))  commande.acm<-paste(nom.res, '<-MCA(', activeDataSet(), ', ncp=', ncp, ', ind.sup=NULL, quanti.sup=', length(varActives)+1, ': ', length(varActives)+ length(variableillu), ', quali.sup=', length(varActives)+length(variableillu)+1, ': ', length(varActives)+length(variableillu)+length(variablefact), ', graph = FALSE)', sep="")
        else commande.acm<-paste(nom.res, '<-MCA(', activeDataSet(), ', ncp=', ncp, ', ind.sup=NULL, quanti.sup=', length(varActives)+1, ': ', length(varActives)+ length(variableillu), ', graph = FALSE)', sep="")
      }
      else
      {
        if(!is.null(variablefact)) commande.acm<-paste(nom.res, '<-MCA(', activeDataSet(), ', ncp=', ncp, ', quali.sup=', length(varActives)+1, ': ', length(varActives)+length(variablefact), ', graph = FALSE)', sep="")
        else commande.acm<-paste(nom.res, '<-MCA(', activeDataSet(), ', ncp=', ncp, ', graph = FALSE)', sep="")
      }
    }
    justDoIt(commande.acm)
    logger(commande.acm)
	justDoIt(paste(nom.res,'$call$call <-',deparse(commande.acm),sep=""))

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

    # gestion des graphiques

    if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) {if (get(nom.res)$eig[1,2]==100) doItAndPrint(paste('"No graph can be plot: data are unidimensional"'))}
    if((Rchoix)&(length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0)){
    if (get(nom.res)$eig[1,2]!=100) {
      if (is.null(Rinvis)) Rinvis <- ""
	  if (Rinvis!="") commande.plotInd<-paste('plot.MCA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), new.plot=TRUE, col.ind="', Rcol.ind, '", col.ind.sup="', Rcol.ind.sup, '", col.var="', Rcol.quali, '", col.quali.sup="', Rcol.qualisup, '", label=c("', paste(Rlabel, collapse='", "'), '"), invisible=c("', paste(Rinvis, collapse='", "'), '")', sep="")
      else  commande.plotInd<-paste('plot.MCA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), new.plot=TRUE, col.ind="', Rcol.ind, '", col.ind.sup="', Rcol.ind.sup, '", col.var="', Rcol.quali, '", col.quali.sup="', Rcol.qualisup, '", label=c("', paste(Rlabel, collapse='", "'), '")', sep="")
      if (!is.null(RTitle)) {
        if (RTitle !=" ") commande.plotInd <- paste(commande.plotInd,', title="', RTitle,'"', sep="")
      }
      if (!is.null(RXlimInd)) commande.plotInd <- paste(commande.plotInd,', xlim=c(', paste(RXlimInd, collapse=", "), ')', sep="")
      if (!is.null(RYlimInd)) commande.plotInd <- paste(commande.plotInd,', ylim=c(', paste(RYlimInd, collapse=", "), ')', sep="")
      commande.plotInd <- paste(commande.plotInd,')', sep="")
      justDoIt(commande.plotInd)
      logger(commande.plotInd)
    }}

    if((Vchoix)&(length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0)){
    if (get(nom.res)$eig[1,2]!=100) {
      if (is.null(Vinvis)) Vinvis <- ""
	  if (Vinvis!="") commande.plotInd<-paste('plot.MCA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), new.plot=TRUE, choix="var", col.var="',Vcol.var,'", col.quali.sup="',Vcol.quali.sup,'", label=c("', paste(Vlabel, collapse='", "'), '"), invisible=c("', paste(Vinvis, collapse='", "'), '")', sep="")
      else commande.plotInd<-paste('plot.MCA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), new.plot=TRUE, choix="var", col.var="',Vcol.var,'", col.quali.sup="',Vcol.quali.sup,'", label=c("', paste(Vlabel, collapse='", "'),'")', sep="")
	  if (!is.null(VTitle)) {
        if (VTitle !=" ") commande.plotInd <- paste(commande.plotInd,', title="', VTitle,'"', sep="")
      }
      commande.plotInd <- paste(commande.plotInd,')', sep="")
      justDoIt(commande.plotInd) 
      logger(commande.plotInd)
    }}

    if((Wchoix)&(length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0)){
    if (get(nom.res)$eig[1,2]!=100) {
      commande.plotInd<-paste('plot.MCA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), new.plot=TRUE, choix="quanti.sup", col.quanti.sup="',Wcol.quanti.sup,'"',', label=c("', paste(Wlabel, collapse='", "'),'")', sep="")
      if (!is.null(WTitle)) {
        if (WTitle !=" ") commande.plotInd <- paste(commande.plotInd,', title="', WTitle,'"', sep="")
      }
#      if ("quanti.sup"%in%Wlabel) commande.plotInd <- paste(commande.plotInd, ',label=c("quanti.sup")',sep='') 
      commande.plotInd <- paste(commande.plotInd,')', sep="")
      justDoIt(commande.plotInd) 
      logger(commande.plotInd)
    }}

    # gestion de l'édition de certains resultats
    doItAndPrint(paste('summary(',nom.res,', nb.dec = 3, nbelements=10, nbind = 10, ncp = 3, file="")', sep=""))
    if (RFichier==""){
      if(Rpropre) doItAndPrint(paste( nom.res, '$eig', sep=""))
      if(Rvariable) doItAndPrint(paste( nom.res, '$var', sep=""))
      if(Rindividu) doItAndPrint(paste( nom.res, '$ind', sep=""))
      if(Rindsup & !is.null(individuillu)) doItAndPrint(paste( nom.res, '$ind.sup', sep=""))
      if(Rquantisup & !is.null(variableillu)) doItAndPrint(paste( nom.res, '$quanti.sup', sep=""))
      if(Rqualisup & !is.null(variablefact)) doItAndPrint(paste( nom.res, '$quali.sup', sep=""))
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
      if(Rvariable){
        doItAndPrint(paste('write.infile(', nom.res, '$var, file =',Fich,',append=',append,')', sep=""))
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
      if(Rquantisup){
        doItAndPrint(paste('write.infile(', nom.res, '$quanti.sup, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rqualisup){
        doItAndPrint(paste('write.infile(', nom.res, '$quali.sup, file =',Fich,',append=',append,')', sep=""))
        append = TRUE
      }
      if(Rdescdim) doItAndPrint(paste('write.infile(dimdesc(', nom.res, ', axes=1:',ncp,'), file =',Fich,',append=',append,')', sep=""))
    }

    # Re-chargement du tableau de départ et supression du tableau temporaire
    activeDataSet(donnee.depart)
    justDoIt(paste('remove(',activeDataSet(),'.MCA)',sep=""))
    logger(paste('remove(',activeDataSet(),'.MCA)',sep=""))
  }


################################################################################
#                   Création de la fenêtre top                                 #
################################################################################
  top<-tktoplevel(borderwidth=10)
  tkwm.title(top, .Facto_gettext("MCA"))
  tkwm.geometry(top, "-100+50")

  # définition des polices
  font2<-tkfont.create(family="times",size=12,weight="bold")
  fontheading<-tkfont.create(family="times",size=18,weight="bold")

  # récupération du jeu de données actif
  donnee<-get(getRcmdr(".activeDataSet"))
  vars<-colnames(donnee)
  rows<-rownames(donnee)

  # création de la liste pour le choix des variables acives
  listdesc<-tklistbox(top,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scr,...))
  scr <-tkscrollbar(top,repeatinterval=5,command=function(...)tkyview(listdesc,...))
  listdesc.nom<-NULL
  for (i in (1:ncol(donnee))) {
      if (is.factor(donnee[,i])) {
          tkinsert(listdesc,"end",vars[i])
          listdesc.nom<-c(listdesc.nom, vars[i])
      }
  }

  # création de tous les boutons d'options dans IlluFrame
  IlluFrame<- tkframe(top, borderwidth=2)
       # mise en page de IlluFrame

  Fillu.funct(label=.Facto_gettext("Supplementary factors"), firstLabel=.Facto_gettext("Supplementary factors"))
  Dillu.funct(label=.Facto_gettext("Supplementary quantitative variables"), firstLabel=.Facto_gettext("Supplementary quantitative variables"))
  Iillu.funct(label=.Facto_gettext("Supplementary individuals"), firstLabel=.Facto_gettext("Supplementary individuals"))
  PLOT.MCA(label=.Facto_gettext("Graphical options"), firstLabel=.Facto_gettext("Graphical options"))
  Sortie.funct(label=.Facto_gettext("Outputs"), firstLabel=.Facto_gettext("Outputs"))
  tkgrid(FilluFrame, DilluFrame, IilluFrame, columnspan=7)
  tkgrid(tklabel(IlluFrame, text=""))
  tkgrid(PlotFrame, SortieFrame, columnspan=7)
  tkgrid.configure(FilluFrame, column=1, columnspan=1)
  tkgrid.configure(PlotFrame, column=2, columnspan=2,sticky="we")
  tkgrid.configure(SortieFrame, column=4, columnspan=2,sticky="ew")
  tkgrid.configure(DilluFrame, column=3, columnspan=1)
  tkgrid.configure(IilluFrame, column=5, columnspan=1)
  tkgrid.columnconfigure(IlluFrame,0, minsize=25)
  tkgrid.columnconfigure(IlluFrame,7, minsize=25)
  tkgrid.columnconfigure(IlluFrame,2, minsize=40)
  tkgrid.columnconfigure(IlluFrame,4, minsize=40)

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
  tkgrid(ncp.lab, ncp)
  tkgrid(Axe.label,Axe1.entry , Axe2.entry, sticky="w")
  tkgrid(resu.lab, resu)
  tkgrid(tklabel(OptionFrame,text="")) # Ligne de blanc
  tkgrid.configure(ncp.lab, resu.lab, Axe.label, column=1, columnspan=4, sticky="w")
  tkgrid.configure(ncp, resu, column=6, columnspan=2, sticky="e")
  tkgrid.configure(Axe1.entry, column=6, columnspan=1, sticky="w")
  tkgrid.configure(Axe2.entry, column=7, columnspan=1, sticky="e")
  tkgrid.columnconfigure(OptionFrame,0, minsize=25)
  tkgrid.columnconfigure(OptionFrame,5, minsize=40)
  tkgrid.columnconfigure(OptionFrame,8, minsize=25)

  #Frame pour HCPC
  HcpcFrame<-tkframe(top, borderwidth=2)
  Hcpc.funct(label=.Facto_gettext("Perform Clustering after MCA"), firstLabel=.Facto_gettext("Perform Clustering after MCA")) 
  tkgrid(Hcpc2Frame, columnspan=7)
  tkgrid.configure(Hcpc2Frame,column=4, columnspan=1)
  
  appliquer.but<-tkbutton(top, text=.Facto_gettext("Apply"),width=12,command=OnAppliquer, borderwidth=3, fg="#690f96")
  OKCancelHelp(helpSubject="MCA",reset="Reinitializ.funct")

  #TOP
  tkgrid(tklabel(top, text=.Facto_gettext("Multiple Correspondence Analysis (MCA)"),font=fontheading),columnspan=3)
  tkgrid(tklabel(top,text=""))
  tkgrid(tklabel(top, text = .Facto_gettext("Select active variables (by default all the variables are active)"),fg = "darkred"), column=1,columnspan=2, sticky = "w")
  tkgrid(listdesc, scr, sticky = "nw")
  tkgrid.configure(scr, sticky = "ens",column=2)
  tkgrid.configure(listdesc, sticky = "ew", column=1, columnspan=2)
  tkgrid(tklabel(top,text="")) # Ligne de blanc
  tkgrid(IlluFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text="")) # Ligne de blanc
  tkgrid(OptionFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text="")) # Ligne de blanc
  tkgrid(HcpcFrame, column=1, columnspan=1)
  tkgrid(tklabel(top,text="")) # Ligne de blanc 
  # tkgrid(appliquer.but, column=1, columnspan=1)
  # tkgrid(tklabel(top,text="")) # Ligne de blanc
  # tkgrid(buttonsFrame, column=1, columnspan=1, sticky="ew" )
  # tkgrid(tklabel(top,text="")) # Ligne de blanc
  tkgrid(buttonsFrame, appliquer.but)
  tkgrid.configure(buttonsFrame, column=1,sticky="e")
  tkgrid.configure(appliquer.but, column=2,sticky="w")

}
