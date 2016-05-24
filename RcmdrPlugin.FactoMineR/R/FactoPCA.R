FactoPCA <-
function()                                                                                        #FactoPCA function which will ba called by OK button of the PCA window
{

#    Création des fonctions pour les options via nouvelle fenêtre graphique


  #! fonction pour le choix des variables qualitatives supplémentaires
  Fillu.funct<-defmacro(label, firstLabel, expr=
  {                                                                                                         #function to choose supplementary categorical variables
    env<-environment()
    variablefact<-NULL                                                                                      #creation of an empty object: variablefact
    .FilluLabel<-tclVar(paste(firstLabel, "", sep=" "))
    .factors<-Factors()
    OnFillu<-function()                                                                                     #function launched by clicking on button "select supplementary factors"
    {
      if(length(.factors)==0) errorCondition(recall=NULL, message=.Facto_gettext("No Factor available"))      #message if no factor in the dataset

      FilluWin<-tktoplevel()                                                                                #creation of a new tcltk window
      tkwm.title(FilluWin,.Facto_gettext("Choice of supplementary factors"))                                  #title of the window
      #création de la fonction FOK.funct
      FOK.funct<-function()                                                                                 #function launched by clicking on OK after having chosen sup fact
      {
        fact.select<-listfact.nom[as.numeric(tkcurselection(listfact))+1]                                   #creation of a list box with factors
        if(length(fact.select)==0) {                                                                        #if nothing selected, no factor affected to variablefact and the window is closed
          assign("variablefact", NULL, envir=env)
          tclvalue(.FilluLabel)<-paste(firstLabel, "", sep=" ")
          tkconfigure(Fillu.but, fg="black")
          tkdestroy(FilluWin)
          return()
        }
        assign("variablefact", fact.select, envir=env)                                                      #else, selected factors are affected to variablefact then the window is closed
        tclvalue(.FilluLabel)<-paste(label, "", sep=" ")
        tkconfigure(Fillu.but, fg="blue")                                                                   #button "select sup fact" becomes blue
        tkdestroy(FilluWin)
      }

      # création et mise en page de la fenetre Fillu
      listfact<-tklistbox(FilluWin,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrfact,...)) #empty list box
      scrfact <-tkscrollbar(FilluWin,repeatinterval=5,command=function(...)tkyview(listfact,...))                                #scroll bar
      listfact.nom<-NULL
      indice<-0
      for (i in (1:ncol(donnee))) {
        if (is.factor(donnee[,i])) {                                                                                             #we're interested only in categorical variables
          tkinsert(listfact,"end",vars[i])                                                                                       #list box is completed with all available factors
          listfact.nom<-c(listfact.nom,vars[i])
          if(vars[i] %in% variablefact) tkselection.set(listfact, indice)                                                        #factors' selection
          indice<-indice+1
        }
      }

      FOK.but<-tkbutton(FilluWin, text="OK", width=16,command=FOK.funct)                                                         #OK button to launch the function above
                                                                                                                                                      #page setting
      tkgrid(tklabel(FilluWin, text=""))                                                                                                              #blank line
        tkgrid(tklabel(FilluWin, text = .Facto_gettext("Select supplementary factor(s)"), fg = "blue"), column=1, columnspan = 1, sticky = "ew")        #text
        tkgrid(listfact, scrfact, sticky = "nw")                                                                                                      #list box and scroll bar side by side
        tkgrid.configure(scrfact, sticky = "ens", columnspan=1)                                                                                       #scroll bar's configuration
        tkgrid.configure(listfact, sticky = "ew", column=1, columnspan=1)                                                                             #list box's config
        tkgrid(tklabel(FilluWin, text=""))                                                                                                            #blank line
        tkgrid(FOK.but, column=1,columnspan=1, sticky="ew")                                                                                           #OK button
        tkgrid(tklabel(FilluWin, text=""))                                                                                                            #blank line
        tkgrid.columnconfigure(FilluWin,0, minsize=25)                                                                                                #window configuration: column nb 0 with size 25
      tkgrid.columnconfigure(FilluWin,2, minsize=25)                                                                                                  #window configuration: column nb 2 with size 25
  }

   FilluFrame<-tkframe(IlluFrame)                                                                                                                     #new frame of main window to put "select sup fact" button
   if(length(.factors)==0){                                                                                                                           #if no factors in the dataset, grey button with message
     Fillu.but<-tkbutton(FilluFrame, text=.Facto_gettext("No factors available"), borderwidth=3)
     tkconfigure(Fillu.but, fg="grey")
   }
   else Fillu.but<-tkbutton(FilluFrame, textvariable=.FilluLabel, command=OnFillu, borderwidth=3)                                                     #else, function OnFillu is launched
   tkgrid(Fillu.but, sticky="ew")                                                                                                                     #setting of this button in the new frame

##   Fillu.but<-tkbutton(FilluFrame, textvariable=.FilluLabel, command=OnFillu, borderwidth=3)
##   tkgrid(Fillu.but, sticky="ew")
  })

  #! fonction pour le choix des variables quantitatives supplémentaires
  Dillu.funct<-defmacro(label, firstLabel, expr=                                                                                         #function to choose supplementary continuous variables
  {
    env<-environment()
    variableillu<-NULL                                                                                                                   #creation of an empty object: variableillu
    .DilluLabel<-tclVar(paste(firstLabel, "", sep=" "))
    OnDillu<-function()                                                                                                                  #function launched by clicking on "select sup variables"
    {
      DilluWin<-tktoplevel()                                                                                                             #creation of a new tcltk window
      tkwm.title(DilluWin,.Facto_gettext("Select supplementary variables"))                                                                #title of the window
      #création de la fonction DOK.funct
      DOK.funct<-function()                                                                                                              #function launched by clicking on OK after habing chosen sup cont variables
      {
        vsup.select<-listvar.nom[as.numeric(tkcurselection(listvar))+1]                                                                  #list box
        if(length(vsup.select)==0)                                                                                                       #if nothing selected, no variable affected to variableillu and the window is closed
        {
          assign("variableillu", NULL, envir=env)
          tclvalue(.DilluLabel)<-paste(firstLabel, "", sep=" ")
          tkconfigure(Dillu.but, fg="black")
          tkdestroy(DilluWin)
          return()
        }
        assign("variableillu", vsup.select, envir=env)                                                                                   #else, chosen variables are affected to variableillu and then the window is closed
        tclvalue(.DilluLabel)<-paste(label, "", sep=" ")
        tkconfigure(Dillu.but, fg="blue")                                                                                                #"select sup var" button  becomes blue
        tkdestroy(DilluWin)
      }

      # création et mise en page de la fenetre Dillu
      listvar<-tklistbox(DilluWin,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrvar,...))            #empty list box
      scrvar <-tkscrollbar(DilluWin,repeatinterval=5,command=function(...)tkyview(listvar,...))                                           #scroll bar
      listvar.nom<-NULL
      indice<-0
      for (i in (1:ncol(donnee))) {
          if (is.numeric(donnee[,i])) {                                                                                                   #we're interested only in numerical variables
            tkinsert(listvar,"end",vars[i])                                                                                               #continuous variables are put in list box
            listvar.nom<-c(listvar.nom,vars[i])
            if(vars[i] %in% variableillu) tkselection.set(listvar, indice)                                                                #selected cont var
            indice<-indice+1
          }
      }

      DOK.but<-tkbutton(DilluWin, text="OK", width=16,command=DOK.funct)                                                                  #OK button to launch the function above

        tkgrid(tklabel(DilluWin, text=""))                                                                                                           #page setting: blank line
        tkgrid(tklabel(DilluWin, text = .Facto_gettext("Select supplementary variable(s)"), fg = "blue"), column=1, columnspan = 1, sticky = "ew")     #text
        tkgrid(listvar, scrvar, sticky = "nw")                                                                                                       #list box and scroll bar side by side
        tkgrid.configure(scrvar, sticky = "ens", columnspan=1)                                                                                       #scroll bar's configuration
        tkgrid.configure(listvar, sticky = "ew", column=1, columnspan=1)                                                                             #list box's configuration
        tkgrid(tklabel(DilluWin, text=""))                                                                                                           #blank line
        tkgrid(DOK.but, column=1,columnspan=1, sticky="ew")                                                                                          #OK button
        tkgrid(tklabel(DilluWin, text=""))                                                                                                           #blank line
        tkgrid.columnconfigure(DilluWin,0, minsize=25)                                                                                               #window configuration: column 0 will be of size 25
        tkgrid.columnconfigure(DilluWin,2, minsize=25)                                                                                               #window configuration: column 2 will be of size 25
  }

   DilluFrame<-tkframe(IlluFrame)                                                                                                                    #new frame on the main window to put "select sup var" button
   Dillu.but<-tkbutton(DilluFrame, textvariable=.DilluLabel, command=OnDillu, borderwidth=3)                                                         #"select sup var" button
   tkgrid(Dillu.but, sticky="ew")                                                                                                                    #"select sup var" button put in its frame
  })


  #! fonction pour le choix des individus supplémentaires
  Iillu.funct<-defmacro(label, firstLabel, expr=                                                                                  #function to choose sup individuals
  {
    env<-environment()
    individuillu<-NULL                                                                                                            #creation of an empty object: individuillu
    .IilluLabel<-tclVar(paste(firstLabel, "", sep=" "))
    OnIillu<-function()                                                                                                           #function launched by clicking on "select sup ind"
    {
      IilluWin<-tktoplevel()                                                                                                      #creation of a new tcltk window
      tkwm.title(IilluWin,.Facto_gettext("Select supplementary individuals"))                                                       #title of the window
      #création de la fonction IOK.funct
      IOK.funct<-function()                                                                                                       #function to select sup ind
      {
        ind.select<-rows[as.numeric(tkcurselection(listind))+1]                                                                   #list box
        if(length(ind.select)==0)                                                                                                 #if no ind selected, nothing put in individuillu and the window is closed
        {
          assign("individuillu", NULL, envir=env)
          tclvalue(.IilluLabel)<-paste(firstLabel, "", sep=" ")
          tkconfigure(Iillu.but, fg="black")
          tkdestroy(IilluWin)
          return()
        }
        assign("individuillu", ind.select, envir=env)                                                                             #else selected ind are put in individuillu then the window is closed
#        tclvalue(.IilluLabel)<-paste(label, ": OK", sep=" ")
        tclvalue(.IilluLabel)<-paste(label, "", sep=" ")
        tkconfigure(Iillu.but, fg="blue")                                                                                         #"select sup ind" button becomes blue
        tkdestroy(IilluWin)
      }

      # création et mise en page de la fenetre Fillu
      listind<-tklistbox(IilluWin,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scrind,...))    #empty list box
      scrind <-tkscrollbar(IilluWin,repeatinterval=5,command=function(...)tkyview(listind,...))                                   #scroll bar
      indice<-0
      for (i in (1:nrow(donnee))) {
          tkinsert(listind,"end",rows[i]) # On renseigne la liste                                                                 #list box is completed with individuals
          if(rows[i] %in% individuillu) tkselection.set(listind,indice)                                                           #selected individuals
          indice<-indice+1
      }

      IOK.but<-tkbutton(IilluWin, text="OK", width=16,command=IOK.funct)                                                          #OK button to launch the function above

      tkgrid(tklabel(IilluWin, text=""))                                                                                                             #page setting: blank line
      tkgrid(tklabel(IilluWin, text = .Facto_gettext("Select supplementary individual(s)"), fg = "blue"), column=1, columnspan = 1, sticky = "ew")     #text
      tkgrid(listind, scrind, sticky = "nw")                                                                                                         #list box and scroll bar side by side
      tkgrid.configure(scrind, sticky = "ens", columnspan=1)                                                                                         #scroll bar's config
      tkgrid.configure(listind, sticky = "ew", column=1, columnspan=1)                                                                               #list box's config
      tkgrid(tklabel(IilluWin, text=""))                                                                                                             #blank line
      tkgrid(IOK.but, column=1,columnspan=1, sticky="ew")                                                                                            #OK button
      tkgrid(tklabel(IilluWin, text=""))                                                                                                             #blank line
      tkgrid.columnconfigure(IilluWin,0, minsize=25)                                                                                                 #configuration of column nb 0: size 25
      tkgrid.columnconfigure(IilluWin,2, minsize=25)                                                                                                 #configuration of column nb 2: size 25
  }

   IilluFrame<-tkframe(IlluFrame)                                                                                                                    #new frame in the main window
   Iillu.but<-tkbutton(IilluFrame, textvariable=.IilluLabel, command=OnIillu, borderwidth=3)                                                         #creation of "select sup ind" button
   tkgrid(Iillu.but, sticky="ew")                                                                                                                    #new button put in new frame
  })


  #! fonction pour la gestion des options graphiques
  PLOT.PCA<-defmacro(label, firstLabel, expr=                                      #function to choose graphical options
  {
    env<-environment()
    compteur.graph<-0                                                              #variables initialization
    #déclaration des variables
    Rchoix<-TRUE                                                                   #plotting graph of individuals
    RTitle<-NULL                                                                   #title for graph of ind
    Rinvisible<-NULL                                                               #invisible elements for graph of ind
    Rlabel<-c("ind", "ind.sup", "quali")                                           #different lables for graph of ind
    Rcol.ind<-Rcol.ind.tmp<-"black"                                                #color for ind
    Rcol.ind.sup<-Rcol.ind.sup.tmp<-"blue"                                         #color for sup ind
    Rcol.quali<-Rcol.quali.tmp<-"magenta"                                          #color for sup fact
    Rhabillage<-"none"                                                             #color according...
    RXlimInd<-NULL                                                                 #limit for x axes
    RYlimInd<-NULL                                                                 #limit for y axes

    Wchoix=TRUE                                                                    #plotting graph of variables
    WTitle<-NULL                                                                   #title of graph of var
    Wlabel<-c("var", "quanti.sup")                                                 #labels
    Wlim.cos<-0.                                                                   #selection on cos2
    Wcol.quanti.sup<-Wcol.quanti.sup.tmp<-"blue"                                   #color for sup var
    Wcol.var<-Wcol.var.tmp<-"black"                                                #color for var
    WXlimVar<-NULL                                                                 #limit for x axes
    WYlimVar<-NULL                                                                 #limit for y axes

    .PlotLabel<-tclVar(paste(firstLabel, "", sep=" "))

    OnPlot<-function()                                                                        #function to be launched by clicking on "graphical options" button
    {
      PlotWin<-tktoplevel()                                                                   #new tcltk window
      tkwm.title(PlotWin, .Facto_gettext("Graphical options"))                                  #title
      tkwm.geometry(PlotWin, "-100+50")                                                       #size

      #création de la fonction onOKsub
      onOKsub<-function()                                                                     #function launched by clicking on OK
      {
        assign("compteur.graph", compteur.graph+1, envir=env)                                 #after clicking on OK, "graph opt" button becomes blue
#        if(compteur.graph>0) tclvalue(.PlotLabel)<-paste(label, ": Seen", sep=" ")
        if(compteur.graph>0) tclvalue(.PlotLabel)<-paste(label, "", sep=" ")
        tkconfigure(Plot.but, fg="blue")

        # gestion des entrées de la partie graphique des individus
        if(tclvalue(ind.check.value)==1) assign("Rchoix", TRUE, envir=env)                    #if the user has chosen to plot the graph of individuals, TRUE is affected to Rchoix, else FALSE is affected
        else assign("Rchoix", FALSE, envir=env)

        if(Rchoix)                                                                            #if Rchoix is TRUE
        {
          if (tclvalue(Titre)==" ") assign("RTitle", NULL, envir=env)                         #if the user did not write anything for the title, then RTitle becomes NULL, else it takes what the user wrote
          assign("RTitle", tclvalue(Titre), envir=env)

          label.tmp.ind<-tclvalue(label.ind.checkValue)
          label.tmp.ind.sup<-tclvalue(label.ind.sup.checkValue)
          label.tmp.quali.sup<-tclvalue(label.quali.sup.checkValue)
          assign("Rlabel", NULL, envir=env)                                                   #at first Rlabel is NULL, then it is completed with "ind", "ind.sup" and "quali" if the user has chosen to plot those labels
          if(label.tmp.ind==1) assign("Rlabel", c(Rlabel, "ind"), envir=env)
          if(label.tmp.ind.sup==1) assign("Rlabel", c(Rlabel, "ind.sup"), envir=env)
          if(label.tmp.quali.sup==1) assign("Rlabel", c(Rlabel, "quali"), envir=env)

          assign("Rcol.ind", Rcol.ind.tmp, envir=env)                                         #Rcol.ind takes the color the user has chosen for individuals
          assign("Rcol.ind.sup", Rcol.ind.sup.tmp, envir=env)                                 #Rcol.ind.sup takes the color the user has chosen for sup individuals
          assign("Rcol.quali", Rcol.quali.tmp, envir=env)                                     #Rcol.quali takes the color the user has chosen for sup factors

          inv.ind.tmp<-tclvalue(inv.ind.checkValue)
          inv.ind.sup.tmp<-tclvalue(inv.ind.sup.checkValue)
          inv.quali.tmp<-tclvalue(inv.quali.checkValue)
          assign("Rinvisible", NULL, envir=env)                                               #at first Rinvisible is NULL, then it is completed with "ind", "ind.sup" and "quali" if the user has chosen to hide those elements
          if(inv.ind.tmp=="1") assign("Rinvisible", c(Rinvisible, "ind"), envir=env)
          if(inv.ind.sup.tmp=="1") assign("Rinvisible", c(Rinvisible, "ind.sup"), envir=env)
          if(inv.quali.tmp=="1") assign("Rinvisible", c(Rinvisible, "quali"), envir=env)

          habillage.tmp<-listgraph.nom[as.numeric(tkcurselection(listgraph))+1]               #individuals are colored according to what the user chose
          if(length(habillage.tmp)==0) assign("Rhabillage","none", envir=env)                 #if nothing chosen, hbillage ="none"
          else assign("Rhabillage", habillage.tmp, envir=env)

          if(tclvalue(XlimIndMin)=="" | tclvalue(XlimIndMax)=="") assign("RXlimInd", NULL, envir=env)                   #when no choice for x limits, RXlimInd=NULL
          else assign("RXlimInd", c(as.numeric(tclvalue(XlimIndMin)), as.numeric(tclvalue(XlimIndMax))), envir=env)     #else RXlimInd is what the user chose
          if(tclvalue(YlimIndMin)=="" | tclvalue(YlimIndMax)=="") assign("RYlimInd", NULL, envir=env)                   #when no choice for y limits, RYlimInd=NULL
          else assign("RYlimInd", c(as.numeric(tclvalue(YlimIndMin)), as.numeric(tclvalue(YlimIndMax))), envir=env)     #else RYlimInd is what the user chose
        }


        # gestion des entrées de la partie graphique des variables
        if(tclvalue(var.check.value)==1) assign("Wchoix", TRUE, envir=env)                     #if the user has chosen to plot the graph of variables, TRUE is affected to Wchoix, else FALSE is affected
        else assign("Wchoix", FALSE, envir=env)

        if(Wchoix)                                                                             #if Wchoix=T
        {
          if (tclvalue(WTitre)==" ") assign("WTitle", NULL, envir=env)                         #if the user did not write anything for the title, then WTitle becomes NULL, else it takes what the user wrote
          assign("WTitle", tclvalue(WTitre), envir=env)

          assign("Wlim.cos", tclvalue(WlimCosValue), envir=env)                                #Wlim.cos takes the value the user chose

          label.tmp.var<-tclvalue(label.var.checkValue)
          label.tmp.quanti.sup<-tclvalue(label.quanti.sup.checkValue)
          assign("Wlabel", NULL, envir=env)                                                    #at first Wlabel is NULL, then it is completed with "var" and "quanti.sup" if the user has chosen to plot those labels
          if(label.tmp.var==1) assign("Wlabel", c(Wlabel, "var"), envir=env)
          if(label.tmp.quanti.sup==1) assign("Wlabel", c(Wlabel, "quanti.sup"), envir=env)

          assign("Wcol.var", Wcol.var.tmp, envir=env)                                          #Wcol.var takes the color the user has chosen for active variables
          assign("Wcol.quanti.sup", Wcol.quanti.sup.tmp, envir=env)                            #Wcol.quanti.sup takes the color the user has chosen for supplementary variables

        }

        tkdestroy(PlotWin)                                                                     #window is closed
      }

    # construction de la partie graphique des individus
    PlotIndFrame<-tkframe(PlotWin, borderwidth=5, relief="groove")                             #new frame for all the choices concerning graph of individuals

    RchoixFrame<-tkframe(PlotIndFrame,borderwidth=2)                                           #frame for the button to choose whether or not to plot the graph of individuals
    ind.check<-tkcheckbutton(RchoixFrame)                                                      #creation of a check button
    if(Rchoix) ind.check.value<-tclVar("1")                                                    #ind.check.value=1 if Rchoix=T
    else ind.check.value<-tclVar("0")                                                          #else ind.check.value=0
    tkconfigure(ind.check, variable=ind.check.value)                                           #assigning variable to check button: will appear as clicked  or unclicked according to ind.check.value (1 or 0)
    tkgrid(tklabel(RchoixFrame, text=.Facto_gettext("Plot individuals graph"), font=font2),ind.check)       #positionning text and check button side by side
    tkgrid(tklabel(RchoixFrame, text="  "))                                                               #blank line

    RTitleFrame<-tkframe(PlotIndFrame,borderwidth=2)                                            #frame for the title
    if (is.null(RTitle)) Titre <- tclVar(" ")                                                   #default value is nothing
    else Titre<-tclVar(RTitle)                                                                  #if the user writes sth, it becomes the value of RTitle
    Titre.entry <-tkentry(RTitleFrame,width="40",textvariable=Titre)                            #case to enter the title
    tkgrid(tklabel(RTitleFrame,text=.Facto_gettext("Title of the graph")),Titre.entry)            #positionning text and case side by side

    RinvisibleFrame<-tkframe(PlotIndFrame,borderwidth=2)                                        #frame for invisible elements
    inv.ind.check<-tkcheckbutton(RinvisibleFrame)                                               #check button for individuals
    if ("ind" %in% Rinvisible) inv.ind.checkValue<-tclVar("1")                                  #inv.ind.checkValue=1 if "ind" is in Rinvisible
    else inv.ind.checkValue<-tclVar("0")
    inv.ind.sup.check<-tkcheckbutton(RinvisibleFrame)                                           #check button for sup individuals
    if ("ind.sup" %in% Rinvisible) inv.ind.sup.checkValue<-tclVar("1")                          #inv.ind.sup.checkValue=1 if "ind.sup" is in Rinvisible
    else inv.ind.sup.checkValue<-tclVar("0")
    inv.quali.check<-tkcheckbutton(RinvisibleFrame)                                             #check button for sup individuals
    if ("quali" %in% Rinvisible) inv.quali.checkValue<-tclVar("1")                              #inv.quali.checkValue<=1 if "quali" is in Rinvisible
    else inv.quali.checkValue<-tclVar("0")
    tkconfigure(inv.ind.check, variable=inv.ind.checkValue)                                     #assigning inv.ind.checkValue to check button, button will appear clicked or unclicked according to it
    tkconfigure(inv.ind.sup.check, variable=inv.ind.sup.checkValue)                             #assigning inv.ind.sup.checkValue to check button
    tkconfigure(inv.quali.check, variable=inv.quali.checkValue)                                 #assigning inv.ind.sup.checkValue to check button
    if (!is.null(variablefact)|!is.null(individuillu)){                                         #if there are supplementary factors or individuals
      tkgrid(tklabel(RinvisibleFrame, text=.Facto_gettext("Hide some elements:")), columnspan=6, sticky="w")           #positionning text
      tkgrid(tklabel(RinvisibleFrame, text="ind"),inv.ind.check, tklabel(RinvisibleFrame, text="ind sup"),inv.ind.sup.check, tklabel(RinvisibleFrame, text="quali"),inv.quali.check, sticky="w")       #and check buttons
    }

    RlabelFrame<-tkframe(PlotIndFrame,borderwidth=2)                                                 #frame for labels
    label.ind.check<-tkcheckbutton(RlabelFrame)                                                      #check button for labels of individuals
    if ("ind" %in% Rlabel) label.ind.checkValue<-tclVar("1")                                         #label.ind.checkValue=1 if "ind" is in Rlabel
    else label.ind.checkValue<-tclVar("0")
    tkconfigure(label.ind.check, variable=label.ind.checkValue)                                      #assigning label.ind.checkValue to the button
    tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Label for the active individuals")),label.ind.check)    #positionning text and button
    label.ind.sup.check<-tkcheckbutton(RlabelFrame)                                                        #check button for labels of sup individuals
    if ("ind.sup" %in% Rlabel) label.ind.sup.checkValue<-tclVar("1")                                       #label.ind.sup.checkValue=1 if "ind.sup" is in Rlabel
    else label.ind.sup.checkValue<-tclVar("0")
    tkconfigure(label.ind.sup.check, variable=label.ind.sup.checkValue)                                    #assigning label.ind.sup.heckValue to the button
    if(!is.null(individuillu)) tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Label for the supplementary individuals")),label.ind.sup.check) #positionning text and button
    label.quali.sup.check<-tkcheckbutton(RlabelFrame)                                                         #check button for labels of sup fact
    if ("quali" %in% Rlabel) label.quali.sup.checkValue<-tclVar("1")                                          #label.quali.sup.checkValue=1 if "ind.sup" is in Rlabel
    else label.quali.sup.checkValue<-tclVar("1")
    tkconfigure(label.quali.sup.check, variable=label.quali.sup.checkValue)                                   #assigning label.quali.sup.checkValue to the button
    if(!is.null(variablefact)) tkgrid(tklabel(RlabelFrame, text=.Facto_gettext("Label for the supplementary factor")), label.quali.sup.check) #positionning text and button

    RcolFrame<-tkframe(PlotIndFrame,borderwidth=2)                                                             #frame for choice of colors
    Rcol.ind.value <- Rcol.ind
    canvas.ind <- tkcanvas(RcolFrame,width="80",height="25",bg=Rcol.ind.value)                                 #button to change color of individuals
    ChangeColor.ind <- function()                                                                              #function to change color of individuals
    {
      Rcol.ind.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Rcol.ind.value,title=.Facto_gettext("Choose a color")))   #assigning chosen color to Rcol.ind.value
      if (nchar(Rcol.ind.value)>0)
      {
        tkconfigure(canvas.ind,bg=Rcol.ind.value)                                                                        #canva's color becomes the chosen one
        assign("Rcol.ind.tmp", Rcol.ind.value, envir=env)                                                                #assigning Rcol.ind.value to Rcol.ind.tmp
      }
    }
    ChangeColor.ind.button <- tkbutton(RcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.ind)              #button to launch function above
    tkgrid(tklabel(RcolFrame, text=.Facto_gettext("Color of the active individuals")),canvas.ind,ChangeColor.ind.button)   #settings: text, canva and button

    Rcol.ind.sup.value<-Rcol.ind.sup                                                                                     #same for supplementary individuals
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

    Rcol.quali.value<- Rcol.quali                                                                                          #same for sup factors
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
    if(!is.null(variablefact)) tkgrid(tklabel(RcolFrame, text=.Facto_gettext("Color for factors")),canvas.quali,ChangeColor.quali.button)

    RhabillageFrame<-tkframe(PlotIndFrame,borderwidth=2)                                                                                               #new frame for coloring individuals according to sth
    listgraph<-tklistbox(RhabillageFrame,height=4, selectmode="single",exportselection="FALSE",yscrollcommand=function(...) tkset(scrgraph,...))       #empty list box
    scrgraph <-tkscrollbar(RhabillageFrame,repeatinterval=5,command=function(...)tkyview(listgraph,...))            #scroll bar
    listgraph.nom<-"ind"
    tkinsert(listgraph,"end",.Facto_gettext("by.individual"))                                                         #putting "by.individual" in the list box
    if(Rhabillage=="ind") tkselection.set(listgraph, 0)                                                             #color by individual if the user chose it
    indice<-1
##    for (i in 1:ncol(donnee)) {
##      if(is.factor(donnee[,i])) {
##        tkinsert(listgraph,"end",vars[i])
##        listgraph.nom<-c(listgraph.nom,vars[i])
##        if(vars[i]==Rhabillage) tkselection.set(listgraph, indice)
##        indice<-indice+1
##      }
##    }
    if (!is.null(variablefact)){                                                                                    #if there are supplementary variables, they are added to the list box
      for (j in 1:ncol(donnee)){
       if(vars[j] %in% variablefact){
        tkinsert(listgraph,"end",vars[j])
        listgraph.nom<-c(listgraph.nom,vars[j])
        if(Rhabillage==vars[j]) tkselection.set(listgraph, indice)                                                  #selection of the chosen factor
        indice<-indice+1
      }}
    }
    tkgrid(tklabel(RhabillageFrame, text=.Facto_gettext("Coloring for individuals")))                                #settings: text
    tkgrid(listgraph, scrgraph, sticky = "nw")                                                                     #list box and scroll bar
    tkgrid.configure(scrgraph, sticky = "wns")                                                                     #configuration for scroll bar
    tkgrid.configure(listgraph, sticky = "ew")                                                                     #config for list box

    RlimFrame<-tkframe(PlotIndFrame,borderwidth=2)                                                                 #new frame for x and y limits
    if(is.null(RXlimInd)) XlimIndMin<-tclVar("")                                                                   #if nothing is written, XlimIndMin=NULL
    else XlimIndMin<-tclVar(paste(RXlimInd[1]))                                                                    #else it equal to RXlimInd
    XlimIndMin.entry <-tkentry(RlimFrame,width="5",textvariable=XlimIndMin)                                        #entry for x minimum limit value
    if (is.null(RXlimInd)) XlimIndMax<- tclVar("")                                                                 #same for x maximum limit
    else XlimIndMax<-tclVar(paste(RXlimInd[1]))
    XlimIndMax.entry <-tkentry(RlimFrame,width="5",textvariable=XlimIndMax)
    tkgrid(tklabel(RlimFrame,text=.Facto_gettext("x limits of the graph:")),XlimIndMin.entry, XlimIndMax.entry)      #positionning of text, entry for min and entry for max limits for axis x
    if(is.null(RYlimInd)) YlimIndMin<- tclVar("")                                                                  #same for y axis
    else YlimIndMin<-tclVar(paste(RYlimInd[1]))
    YlimIndMin.entry <-tkentry(RlimFrame,width="5",textvariable=YlimIndMin)
    if (is.null(RYlimInd)) YlimIndMax<- tclVar("")
    else YlimIndMax<-tclVar(paste(RYlimInd[2]))
    YlimIndMax.entry <-tkentry(RlimFrame,width="5",textvariable=YlimIndMax)
    tkgrid(tklabel(RlimFrame,text=.Facto_gettext("y limits of the graph:")),YlimIndMin.entry,YlimIndMax.entry)

    #mise en page des différents frames de PlotIndFrame
    tkgrid(RchoixFrame)                                                                                            #positionning of the different frames:
    tkgrid(RTitleFrame)                                                                                            #title frame
    tkgrid(RinvisibleFrame)                                                                                        #invisible elements frame
    tkgrid(RlabelFrame)                                                                                            #labels frame
    tkgrid(tklabel(PlotIndFrame, text=" "))                                                                        #blank line
    tkgrid(RcolFrame)                                                                                              #colors frame
    tkgrid(RhabillageFrame)                                                                                        #coloring frame
    tkgrid(tklabel(PlotIndFrame, text=" "))                                                                        #blank line
    tkgrid(RlimFrame)                                                                                              #axes limits frame
    tkgrid(tklabel(PlotIndFrame, text=" "))                                                                        #blank line

    # construction de la partie graphique des variables
    PlotVarFrame<-tkframe(PlotWin, borderwidth=5, relief="groove")                                                 #frame for all the choices for the variables graph

    WchoixFrame<-tkframe(PlotVarFrame,borderwidth=2)
    var.check<-tkcheckbutton(WchoixFrame)                                                                          #frame for choosing whether or not to plot the variables graph
    if(Wchoix) var.check.value<-tclVar("1")                                                                        #1 if yes
    else var.check.value<-tclVar("0")                                                                              #else 0
    tkconfigure(var.check, variable=var.check.value)                                                               #assigning value to check button (appears clicked or unclicked
    tkgrid(tklabel(WchoixFrame, text=.Facto_gettext("Plot variables graph"), font=font2),var.check)                  #positionning
    tkgrid(tklabel(WchoixFrame, text="  "))                                                                        #blank line

    WTitleFrame<-tkframe(PlotVarFrame,borderwidth=2)                                                               #title frame
    if (is.null(WTitle)) WTitre <- tclVar(" ")
    else WTitre<-tclVar(WTitle)
    WTitre.entry <-tkentry(WTitleFrame,width="40",textvariable=WTitre)
    tkgrid(tklabel(WTitleFrame,text=.Facto_gettext("Title of the graph")),WTitre.entry)

    WcosFrame<-tkframe(PlotVarFrame,borderwidth=2)                                                                 #cos2 frame
    WlimCosValue<-tclVar(paste(Wlim.cos))                                                                          #WlimCosValue=default value at first
    WlimCos.entry<-tkentry(WcosFrame, width=5, textvariable=WlimCosValue)                                          #creation of the entry for the user to type a value
    tkgrid(tklabel(WcosFrame,text=.Facto_gettext("Draw variables with a cos2 >:")),WlimCos.entry)                    #positionning text and entry in the frame

    WlabelFrame<-tkframe(PlotVarFrame,borderwidth=2)                                                               #labels frame
    label.var.check<-tkcheckbutton(WlabelFrame)
    if ("var" %in% Wlabel) label.var.checkValue<-tclVar("1")
    else label.var.checkValue<-tclVar("0")
    tkconfigure(label.var.check, variable=label.var.checkValue)
    tkgrid(tklabel(WlabelFrame, text=.Facto_gettext("Labels for the active variables")),label.var.check)
    label.quanti.sup.check<-tkcheckbutton(WlabelFrame)
    if ("quanti.sup" %in% Wlabel) label.quanti.sup.checkValue<-tclVar("1")
    else label.quanti.sup.checkValue<-tclVar("0")
    tkconfigure(label.quanti.sup.check, variable=label.quanti.sup.checkValue)
    if(!is.null(variableillu)) tkgrid(tklabel(WlabelFrame, text=.Facto_gettext("Labels for the supplementary variables")),label.quanti.sup.check)

    WcolFrame<-tkframe(PlotVarFrame,borderwidth=2)                                                                  #color frame
    Wcol.var.value <- Wcol.var
    canvas.var <- tkcanvas(WcolFrame,width="80",height="25",bg=Wcol.var.value)
    ChangeColor.var <- function()
    {
      Wcol.var.value<-tclvalue(tcl("tk_chooseColor",initialcolor=Wcol.var.value,title=.Facto_gettext("Choose a color")))
      if (nchar(Wcol.var.value)>0) {
        tkconfigure(canvas.var,bg=Wcol.var.value)
        assign("Wcol.var.tmp", Wcol.var.value, envir=env)
      }
    }
    ChangeColor.var.button <- tkbutton(WcolFrame,text=.Facto_gettext("Change Color"),command=ChangeColor.var)
    tkgrid(tklabel(WcolFrame, text=.Facto_gettext("Color for active variables")),canvas.var,ChangeColor.var.button)

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
    if(!is.null(variableillu)) tkgrid(tklabel(WcolFrame, text=.Facto_gettext("Color for supplementary variables")),canvas.quanti.sup,ChangeColor.quanti.sup.button)

    #mise en page des différents frames de PlotVarFrame                                                          #positionning all frames
    tkgrid(WchoixFrame)                                                                                          #plotting or not frame
    tkgrid(WTitleFrame)                                                                                          #title frame
    tkgrid(WcosFrame)                                                                                            #cos2 frame
    tkgrid(WlabelFrame)                                                                                          #labels frame
    tkgrid(tklabel(PlotVarFrame, text=" "))                                                                      #blank line
    tkgrid(WcolFrame)                                                                                            #color frame
    tkgrid(tklabel(PlotVarFrame, text=" "))                                                                      #blank line

    subOKCancelHelp(PlotWin, "plot.PCA")                                                                         #creating OK, Cancel and Help buttons. Help being associated to ?plot.PCA
    tkgrid(PlotIndFrame,PlotVarFrame)                                                                            #positionning of individuals and variables frames side by side
    tkgrid.configure(PlotVarFrame, sticky="n")
    tkgrid(subButtonsFrame, sticky="ew", columnspan=2)                                                           #positionning buttons OK, Cancel and Help
    }

    PlotFrame<-tkframe(IlluFrame)                                                                                #new frame in the main window
    Plot.but<-tkbutton(PlotFrame, textvariable=.PlotLabel, command=OnPlot, borderwidth=3)                        #creation of a new button in that new frame. Clicking on this button will open the window with graphical options
    tkgrid(Plot.but, sticky="ew")                                                                                #positionning the new button
  })


  #! fonction pour la réinitialisation des paramètres
  Reinitializ.funct<-function()                                                                                  #function to re-initialize everything
  {
    tkdestroy(top)                                                                                               #window is closed
    FactoPCA()                                                                                                   #and re-opened
  }


  #! fonction pour le choix des éléments de sortie                                                               #choice of output options
  Sortie.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    compteur.sortie<-0
    #déclaration des variables
    RFichier <- ""                                                                                               #write results in a file
    Rpropre<-FALSE                                                                                               #print eigenvalues
    Rvariable<-FALSE                                                                                             #print results for active variables
    Rindividu<-FALSE                                                                                             #print results for active individuals
    Rindsup<-FALSE                                                                                               #print results for supplementary individuals
    Rquantisup<-FALSE                                                                                            #print results for supplementary variables
    Rqualisup<-FALSE                                                                                             #print results for supplementary factors
    Rdescdim<-FALSE                                                                                              #print results of the description of dimensions

    .SortieLabel<-tclVar(paste(firstLabel, "", sep=" "))

    OnSortie<-function()                                                                                         #function launched by clicking on the "outputs" button
    {
      SortieWin<-tktoplevel()                                                                                    #new tcltk window
      tkwm.title(SortieWin,.Facto_gettext("Outputs"))                                                                            #title

      #création de la fonction onOKsub
      onOK.sortie<-function()                                                                                    #function launched by clicking on OK
      {
        assign("compteur.sortie", compteur.sortie+1, envir=env)                                                  #when clicking on OK, "compteur.sortie" goes from 0 to 1 and the "output" button of the main window becomes blue
        if(compteur.sortie>0) tclvalue(.SortieLabel)<-paste(label, "", sep=" ")
        tkconfigure(Sortie.but, fg="blue")

        if(tclvalue(eigValue)=="1") assign("Rpropre", TRUE, envir=env)                                           #if "eigenvalues" check button is checked, Rpropre is TRUE
        else assign("Rpropre", FALSE, envir=env)                                                                 #else it is FALSE

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
                                                                                                                
      eig.lab <-tklabel(SortieWin, text=.Facto_gettext("Eigenvalues"))                                            #label for the "eigenvalues" check button
        eig.check <- tkcheckbutton(SortieWin)                                                                   #creation of the check button
        if(Rpropre) eigValue <- tclVar("1")                                                                     #variable eigValue is equal to 1 or 0 according to the value or Rpropre
        else eigValue <- tclVar("0")
        tkconfigure(eig.check,variable=eigValue)                                                                #the button appears checked or unchecked according to eigValue

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

      SortieOK.but<-tkbutton(SortieWin,text="OK",width=16,command=onOK.sortie)                                                #OK button associated to the function above

        tkgrid(tklabel(SortieWin, text = .Facto_gettext("Select output options"), fg ="blue"),  columnspan = 2, sticky = "w")   #settings. Text
        tkgrid(tklabel(SortieWin, text = " "))                                                                                #blank line
        tkgrid(eig.lab,eig.check,sticky="w")                                                                                  #"eigenvalues" check button
        tkgrid(var.lab,var.check,sticky="w")                                                                                  #"results for var" check button
        tkgrid(ind.lab,ind.check,sticky="w")                                                                                  #"results for ind" check button
        if (!is.null(individuillu)) tkgrid(ind.sup.lab,ind.sup.check,sticky="w")                                              #"results for sup ind" check button if sup ind have been selected
        if (!is.null(variableillu)) tkgrid(quanti.sup.lab,quanti.sup.check,sticky="w")                                        #"results for sup var" check button if sup var have been selected
        if (!is.null(variablefact)) tkgrid(quali.sup.lab,quali.sup.check,sticky="w")                                          #"results for sup fact" check button if sup fact have been selected
        tkgrid(descdim.lab,descdim.check,sticky="w")                                                                          #"results for description of the dimension" check button
        tkgrid(tklabel(SortieWin, text = " "))                                                                                #blank line
        tkgrid(RFichierFrame)                                                                                                 #frame for the option to print or not results in a file
        tkgrid(SortieOK.but)                                                                                                  #OK button
   }

    SortieFrame<-tkframe(IlluFrame)                                                                                           #new frame in the main window
    Sortie.but<-tkbutton(SortieFrame, textvariable=.SortieLabel, command=OnSortie, borderwidth=3)                             #"output" button
    tkgrid(Sortie.but, sticky="ew")                                                                                           #new button positionned in the new frame
  })

  #! fonction HCPC                                                                                           

  Hcpc.funct<-defmacro(label, firstLabel, expr=
  {
    env<-environment()
    .HcpcLabel<-tclVar(paste(firstLabel, "", sep=" "))
    compteur.hcpc<-0
    Rclassif<-0                                                                                      #variable for the choice: perform HCPC or not
    Rmeth <- -1                                                                                      #variable for the choice of the method
    Rconsolid<-0                                                                                     #variable for the choice of consolidation
    Rgraphhcpc<-1                                                                                    #variable for graphs
    Rreshcpc<-0                                                                                      #variable for printing results
    Rminhcpc<-3                                                                                      #variable for the minimum nb of clusters
    Rmaxhcpc<-10                                                                                     #variable for the maximum nb of clusters

    OnHCPC <- function()                                                                             #function launched by clicking on "perform HCPC after PCA"
    {

      HcpcWin<-tktoplevel()                                                                          #creation of a tcltk window
      tkwm.title(HcpcWin, .Facto_gettext("HCPC options"))                                              #title of the window

      onOKHcpc <- function()                                                                         #function launched by clicking on OK
      {
        assign("compteur.hcpc", compteur.hcpc+1, envir=env)                                          #"perform HCPC after PCA" button becomes blue after licking on OK
        if(compteur.hcpc>0) tclvalue(.HcpcLabel)<-paste(label, "", sep=" ")
        tkconfigure(Hcpc.but, fg="blue")

        if(tclvalue(methValue)=="0") assign("Rmeth", 0, envir=env)                         #assigning value to Rmeth: 0 or -1
        else assign("Rmeth", -1, envir=env)

        if(tclvalue(consolidValue)=="1") assign("Rconsolid",TRUE, envir=env)                         #assigning value to Rconsolid: TRUE or FALSE
        else assign("Rconsolid",FALSE,envir=env)

        if(tclvalue(graphhcpcValue)=="1") assign("Rgraphhcpc",TRUE,envir=env)                        #assigning value to Rgraphhcpc: TRUE or FALSE
        else assign("Rgraphhcpc",FALSE,envir=env)

        if(tclvalue(reshcpcValue)=="1") assign("Rreshcpc",TRUE,envir=env)                            #assigning value to Rreshcpc: TRUE or FALSE
        else assign("Rreshcpc",FALSE,envir=env)

        assign("Rminhcpc",as.numeric(tclvalue(minhcpc)),envir=env)                                   #assigning value to Rminhcpc: numeric value
        assign("Rmaxhcpc",as.numeric(tclvalue(maxhcpc)),envir=env)                                   #assigning value to Rmaxhcpc: numeric value

        assign("Rclassif",TRUE,envir=env)                                                            #assigning value to Rclassif: TRUE
        tkdestroy(HcpcWin)                                                                           #window is closed
      }
      OKHcpc.but<-tkbutton(HcpcWin, text="OK", width=8,command=onOKHcpc)  #OK button

      onCancelHcpc <- function()                                          #function launched by clicking on Cancel button
      {
        assign("Rclassif",FALSE,envir=env)                                #Rclassif=FALSE
        tkdestroy(HcpcWin)                                                #the window is closed
      }
      CancelHcpc.but<-tkbutton(HcpcWin, text="Cancel", width=8,command=onCancelHcpc)   #Cancel button

      tkgrid(tklabel(HcpcWin, text=""))                                   #settings: blank line
      tkgrid(tklabel(HcpcWin, text = .Facto_gettext("Hierarchical Clustering on Principal Components"), fg = "darkred"), column=1, columnspan = 8, sticky = "ew")                   #settings: text

      meth1 <- tkradiobutton (HcpcWin)                                    #creation of a radio button
      meth1.lab <- tklabel(HcpcWin,text=.Facto_gettext("interactive"))      #label of the button
      meth2 <- tkradiobutton (HcpcWin)                                    #other radio button
      meth2.lab <- tklabel(HcpcWin,text=.Facto_gettext("automatic"))        #label
      methValue <- tclVar(paste(Rmeth))                                  #variable associated with the radio buttons
      meth.lab <- tklabel(HcpcWin,text=.Facto_gettext("Choice of the number of clusters:"))  #main label
      tkconfigure(meth1,variable=methValue,value="0")                             #assigning value to the first radio button
      tkconfigure(meth2,variable=methValue,value="-1")                               #assigning value to the second radio button
      minmaxhcpc.label<-tklabel(HcpcWin,text=.Facto_gettext("The optimal number of clusters is chosen between:"))                                                              #text for the choice of the nb of clusters

      minhcpc<-tclVar(paste(Rminhcpc))
      maxhcpc<-tclVar(paste(Rmaxhcpc))
      minhcpc.entry <-tkentry(HcpcWin,width="3",textvariable=minhcpc)     #text zone for the min nb
      maxhcpc.entry <-tkentry(HcpcWin,width="3",textvariable=maxhcpc)     #text zone for the max nb

      consolid.lab <- tklabel(HcpcWin,text=.Facto_gettext("Consolidate clusters"))   #label for consolidation
      consolid.check <- tkcheckbutton(HcpcWin)                                      #check button
      if(Rconsolid) consolidValue<-tclVar("1")                                      #assigning value according to Rconsolid's value
      else consolidValue<-tclVar("0")                                               
      tkconfigure(consolid.check,variable=consolidValue)                            #assigning consolidValue to the check button

      graphhcpc.lab <- tklabel(HcpcWin,text=.Facto_gettext("Print graphs"))          #label for graphs
      graphhcpc.check <- tkcheckbutton(HcpcWin)                                     #check button
      if(Rgraphhcpc) graphhcpcValue <- tclVar("1")                                  #value depending on Rgraphhcpc
      else graphhcpcValue <- tclVar("0")
      tkconfigure(graphhcpc.check,variable=graphhcpcValue)                          #assigning graphhcpcValue to the check button

      reshcpc.lab <- tklabel(HcpcWin,text=.Facto_gettext("Print results for clusters"))  #label
      reshcpc.check <- tkcheckbutton(HcpcWin)                                           #check button
      if(Rreshcpc) reshcpcValue<-tclVar("1")                                            #value depending on Rreshcpc
      else reshcpcValue <- tclVar("0")
      tkconfigure(reshcpc.check,variable=reshcpcValue)                                  #assigning reshcpcValue to the check button

      tkgrid(tklabel(HcpcWin,text=.Facto_gettext("Options for the clustering"), fg = "blue"), column=1, columnspan=8, sticky="we")                                                       #settings: text
      tkgrid(tklabel(HcpcWin,text=""),column=1,columnspan=4,sticky="w")                                            #blank line
      tkgrid(tklabel(HcpcWin,text=sprintf(.Facto_gettext("Clustering is performed on the first %s dimensions of PCA"),tclvalue(ncp.val))),column=1,columnspan=4,sticky="w")   #text which takes the nb of dimensions chosen in the main window
#      tkgrid(tklabel(HcpcWin,text=paste(.Facto_gettext("Clustering is performed on the first "),tclvalue(ncp.val), .Facto_gettext(" dimensions of PCA"),sep="")),column=1,columnspan=4,sticky="w")   #text which takes the nb of dimensions chosen in the main window
      tkgrid(tklabel(HcpcWin,text=.Facto_gettext("(Modify in the main options to change this number)")),column=1,columnspan=4,sticky="w")                                     #text
      tkgrid(tklabel(HcpcWin,text=""))                                            #blank line
      tkgrid(meth.lab,meth1.lab,meth1)                                            #positionning text concerning the choice of the method and the first radio button side by side)
      tkgrid(meth2.lab,meth2)                                                     #positionning radio button nb 2 and its label
      tkgrid(tklabel(HcpcWin,text=""))                                            #blank line
      tkgrid(minmaxhcpc.label,minhcpc.entry , maxhcpc.entry)                      #positionning everything that concerns the nb of clusters
      tkgrid(tklabel(HcpcWin,text=""))                                            #blank line
      tkgrid(consolid.lab,consolid.check)                                         #positionning text and button for consolidation
      tkgrid(graphhcpc.lab,graphhcpc.check)                                       #positionning label and button for graphs
      tkgrid(reshcpc.lab,reshcpc.check)                                           #positionning label and button for results
      tkgrid(tklabel(HcpcWin,text=""))                                            #blank line
      tkgrid(OKHcpc.but, CancelHcpc.but)                                          #OK and Cancel buttons
      tkgrid(tklabel(HcpcWin, text=""))                                           #blank line

      tkgrid.configure(minmaxhcpc.label,meth.lab,consolid.lab,graphhcpc.lab,reshcpc.lab,column=1,columnspan=4,sticky="w")                                                         #configuration of everything with the number of the column it will appear in (8 columns in a window) and the position in the column (w: west, e: east)
      tkgrid.configure(minhcpc.entry,column=7,columnspan=1,sticky="e")
      tkgrid.configure(maxhcpc.entry,column=8,columnspan=1,sticky="w")
      tkgrid.configure(meth1,meth2,consolid.check,graphhcpc.check,reshcpc.check,column=8,sticky="e")
      tkgrid.configure(meth1.lab,column=6,columnspan=2,sticky="w")
      tkgrid.configure(meth2.lab,column=6,columnspan=2,sticky="w")
      tkgrid.configure(OKHcpc.but,column=2,columnspan=1,sticky="w")
      tkgrid.configure(CancelHcpc.but,column=6,columnspan=1,sticky="e")

      tkgrid.columnconfigure(HcpcWin,0, minsize=3)                                #configuration of the size of some columns
      tkgrid.columnconfigure(HcpcWin,5, minsize=5)
      tkgrid.columnconfigure(HcpcWin,8, minsize=3)

}
    Hcpc2Frame<-tkframe(HcpcFrame)                                                #new frame in the HCPC frame
    Hcpc.but<-tkbutton(Hcpc2Frame, textvariable=.HcpcLabel, command=OnHCPC, borderwidth=3)  #button in this new frame
    tkgrid(Hcpc.but, sticky="ew")                                                 #positionning the button
})


  #! fonction associée au bouton Appliquer, execute sans détruire la fenêtre top
  OnAppliquer<-function()                                                         #funtion launched by clicking on "Apply"
  {

    # liste de toutes les variables interne créées      (** mise en forme incomplète)
      # sur la fenetre principale
#         listdesc         **
#         resu.val         **
#         ncp.val          **
#         reduitValue      **
#         Axe1
#         Axe2
      # dans les boutons des fenêtres illustratives
#         variablefact     **
#         variableillu     **
#         individuillu     **
      # dans le bouton Plot PCA
#         Rchoix
#         RTitle
#         Rlabel
#         Rcol.ind
#         Rcol.ind.sup
#         Rcol.quali
#         Rhabillage       **
#         RXlimInd
#         RYlimInd
#         Wchoix
#         WTitle
#         Wlabel
#         Wlim.cos
#         Wcol.quanti.sup
#         Wcol.var
#         WXlimVar
#         WYlimVar
      # dans le bouton affichage sortie
#         Rpropre
#             Rvariable
#           Rindividu
#           Rindsup
#           Rquantisup
#           Rqualisup


    # récupération des paramètres de la fenêtre principale
    nom.res<-tclvalue(resu.val)                                                   #putting resu.val in nom.res
    if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) justDoIt(paste('remove (',nom.res,')'))       #if object res already exists, it's removed
    if(length(as.numeric(tkcurselection(listdesc)))<2) varActives<-listdesc.nom
    else varActives<-listdesc.nom[as.numeric(tkcurselection(listdesc))+1]
    varActives <- varActives[!(varActives%in%variableillu)]                       #varActives is the list of selected active variables
    reduction<-TRUE                                                               #scaling=T by default
    if(tclvalue(reduitValue)=="0") reduction<-FALSE                               #if "sclaing" check button is unchecked, scaling=F
    ncp<-as.numeric(tclvalue(ncp.val))                                            #putting chosen nb of dim in ncp
    Axe<-c(as.numeric(tclvalue(Axe1)), as.numeric(tclvalue(Axe2)))                #putting chosen axes in Axe

    # gestion du tableau de données pour l'ACP

    if(!is.null(individuillu)) {                                                  #if sup ind have been chosen, active ind are the remaining ones
      ind.actif<-rows[-which(rows %in% individuillu)]
      if(!is.null(variableillu)) {                                                #a new dataset is created by the function commande.data which uses everything concerning active and supplementary variables and individuals. Different cases are separated according to presence or absence of those elements
        if(!is.null(variablefact)) commande.data<-paste(activeDataSet(),'.', 'PCA', '<-', activeDataSet(), '[c("', paste(ind.actif, collapse='", "'), '", "', paste(individuillu, collapse='", "'), '") ,c("', paste(varActives, collapse='", "'), '", "', paste(variableillu, collapse='", "'), '", "', paste(variablefact, collapse='", "'), '")]', sep="")
        else commande.data<-paste(activeDataSet(),'.', 'PCA', '<-', activeDataSet(), '[c("', paste(ind.actif, collapse='", "'), '", "', paste(individuillu, collapse='", "'), '") ,c("', paste(varActives, collapse='", "'), '", "', paste(variableillu, collapse='", "'), '")]', sep="")
      }
      else {
        if(!is.null(variablefact)) commande.data<-paste(activeDataSet(),'.', 'PCA', '<-', activeDataSet(), '[c("', paste(ind.actif, collapse='", "'), '", "', paste(individuillu, collapse='", "'), '") ,c("', paste(varActives, collapse='", "'), '", "', paste(variablefact, collapse='", "'), '")]', sep="")
        else commande.data<-paste(activeDataSet(),'.', 'PCA', '<-', activeDataSet(), '[c("', paste(ind.actif, collapse='", "'), '", "', paste(individuillu, collapse='", "'), '") ,c("', paste(varActives, collapse='", "'), '")]', sep="")
      }
    }
    else
    {
      if(!is.null(variableillu)) {
        if(!is.null(variablefact)) commande.data<-paste(activeDataSet(),'.', 'PCA', '<-', activeDataSet(), '[, c("', paste(varActives, collapse='", "'), '", "', paste(variableillu, collapse='", "'), '", "', paste(variablefact, collapse='", "'), '")]', sep="")
        else commande.data<-paste(activeDataSet(),'.', 'PCA', '<-', activeDataSet(), '[, c("', paste(varActives, collapse='", "'), '", "', paste(variableillu, collapse='", "'), '")]', sep="")
      }
      else {
        if(!is.null(variablefact)) commande.data<-paste(activeDataSet(),'.', 'PCA', '<-', activeDataSet(), '[, c("', paste(varActives, collapse='", "'), '", "', paste(variablefact, collapse='", "'), '")]', sep="")
        else commande.data<-paste(activeDataSet(),'.', 'PCA', '<-', activeDataSet(), '[, c("', paste(varActives, collapse='", "'), '")]', sep="")
      }
    }
    justDoIt(commande.data)                                       #commande.data is launched
    logger(commande.data)
    donnee.depart<-activeDataSet()
    activeDataSet(paste(activeDataSet(),'.', 'PCA', sep=""))      #the new dataset (where columns and rows have been re-ordered) is called dataset_name.PCA

    # gestion de la commande réalisant l'ACP
    if(!is.null(individuillu)) {
      ind.actif<-rows[-which(rows %in% individuillu)]
      if(!is.null(variableillu)) {                               #PCA is performed by commande.acp which uses the name of the active dataset (created right above), scale.unit, ncp, ind.sup, quanti.sup and quali.sup. graph is FALSE by default
        if(!is.null(variablefact)) commande.acp<-paste(nom.res, '<-PCA(', activeDataSet(), ', scale.unit=', reduction, ', ncp=', ncp, ', ind.sup=c(', nrow(get(getRcmdr(".activeDataSet")))-length(individuillu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))), '), quanti.sup=c(', length(varActives)+1, ': ', length(varActives)+ length(variableillu), '), quali.sup=c(', length(varActives)+length(variableillu)+1, ': ', length(varActives)+length(variableillu)+length(variablefact), '), graph = FALSE)', sep="")
        else commande.acp<-paste(nom.res, '<-PCA(', activeDataSet(), ', scale.unit=', reduction, ', ncp=', ncp, ', ind.sup=c(', nrow(get(getRcmdr(".activeDataSet")))-length(individuillu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))), '), quanti.sup=c(', length(varActives)+1, ': ', length(varActives)+ length(variableillu), '), graph = FALSE)', sep="")
      }
      else {
        if(!is.null(variablefact)) commande.acp<-paste(nom.res, '<-PCA(', activeDataSet(), ', scale.unit=', reduction, ', ncp=', ncp, ', ind.sup=c(', nrow(get(getRcmdr(".activeDataSet")))-length(individuillu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))), '), quali.sup=c(', length(varActives)+1, ': ', length(varActives)+length(variablefact), '), graph = FALSE)', sep="")
        else commande.acp<-paste(nom.res, '<-PCA(', activeDataSet(), ', scale.unit=', reduction, ', ncp=', ncp, ', ind.sup=c(', nrow(get(getRcmdr(".activeDataSet")))-length(individuillu)+1, ': ', nrow(get(getRcmdr(".activeDataSet"))), '), graph = FALSE)', sep="")
      }
    }
    else {
      if(!is.null(variableillu)) {
        if(!is.null(variablefact)) commande.acp<-paste(nom.res, '<-PCA(', activeDataSet(), ' , scale.unit=', reduction, ', ncp=', ncp, ', quanti.sup=c(', length(varActives)+1, ': ', length(varActives)+ length(variableillu), '), quali.sup=c(', length(varActives)+length(variableillu)+1, ': ', length(varActives)+length(variableillu)+length(variablefact), '), graph = FALSE)', sep="")
        else commande.acp<-paste(nom.res, '<-PCA(', activeDataSet(), ' , scale.unit=', reduction, ', ncp=', ncp, ', quanti.sup=c(', length(varActives)+1, ': ', length(varActives)+ length(variableillu), '), graph = FALSE)', sep="")
      }
      else{
        if(!is.null(variablefact)) commande.acp<-paste(nom.res, '<-PCA(', activeDataSet(), ' , scale.unit=', reduction, ', ncp=', ncp, ', quali.sup=c(', length(varActives)+1, ': ', length(varActives)+length(variablefact), '), graph = FALSE)', sep="")
        else commande.acp<-paste(nom.res, '<-PCA(', activeDataSet(), ' , scale.unit=', reduction, ', ncp=', ncp, ', graph = FALSE)', sep="")
      }
    }
    justDoIt(commande.acp)                         #commande.acp is performed
    logger(commande.acp)
	justDoIt(paste(nom.res,'$call$call <-',deparse(commande.acp),sep=""))


    #Commande de la fonction HCPC

    if(Rclassif==TRUE){                            #if Rclassif is TRUE then HCPC is performed with options chosen by the user
      commande.hcpc<-paste(nom.res,'.hcpc', '<-HCPC(', nom.res, ' ,nb.clust=', Rmeth, ',consol=', Rconsolid,',min=', Rminhcpc,',max=',Rmaxhcpc,',graph=', Rgraphhcpc, ')', sep="")
    justDoIt(commande.hcpc)
    logger(commande.hcpc)
      if(Rreshcpc==TRUE){                          #if Rreshcpc is TRUE then results are printed
        doItAndPrint(paste(nom.res,'.hcpc$data.clust[,ncol(res.hcpc$data.clust),drop=F]', sep=""))
        doItAndPrint(paste(nom.res,'.hcpc$desc.var', sep=""))
        doItAndPrint(paste(nom.res,'.hcpc$desc.axes', sep=""))
        doItAndPrint(paste(nom.res,'.hcpc$desc.ind', sep=""))
      }
    }

    # gestion des graphiques

    if (length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0) {if (get(nom.res)$eig[1,2]==100) doItAndPrint(paste('"No graph can be plot: data are unidimensional"'))}     #error message in case data are unidimensional
    if((Rchoix)&length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0){  #command for the graph of individuals
    if (get(nom.res)$eig[1,2]!=100) {
      if ((Rhabillage!="none") & (Rhabillage!="ind")) {
        Rhabillage<-which(colnames(get(getRcmdr(".activeDataSet")))==Rhabillage)
        if(length(Rhabillage)==0) Rhabillage<-"none"
      }
      if (Rhabillage=="none") Rhabillage<-paste('"', Rhabillage, '"', sep="")
      if (Rhabillage=="ind") Rhabillage<-paste('"', Rhabillage, '"', sep="")

      commande.plotInd<-paste('plot.PCA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="ind", habillage=', Rhabillage, ', col.ind="', Rcol.ind, '", col.ind.sup="', Rcol.ind.sup, '", col.quali="', Rcol.quali, '", label=c("', paste(Rlabel, collapse='", "'), '"),new.plot=TRUE', sep="")
      if (!is.null(RXlimInd)) commande.plotInd<-paste(commande.plotInd, ', xlim=c(', paste(RXlimInd, collapse=", "), ')')
      if (!is.null(RYlimInd)) commande.plotInd<-paste(commande.plotInd, ', ylim=c(', paste(RYlimInd, collapse=", "), ')')
      if (!is.null(Rinvisible)) commande.plotInd<-paste(commande.plotInd, ', invisible=c("', paste(Rinvisible, collapse='", "'),'")', sep='')
      if (is.null(RTitle)) commande.plotInd <- paste(commande.plotInd,')', sep="")
      else {
        if (RTitle ==" ") commande.plotInd <- paste(commande.plotInd,')', sep="")
        else commande.plotInd <- paste(commande.plotInd,', title="', RTitle,'")', sep="")
      }
      justDoIt(commande.plotInd)
      logger(commande.plotInd)
    }}


    if((Wchoix)&length(which(ls(envir = .GlobalEnv, all.names = TRUE)==nom.res))>0){  #command for the graph of variables
    if (get(nom.res)$eig[1,2]!=100) {
      commande.plotVar<-paste('plot.PCA(', nom.res, ', axes=c(', paste(Axe, collapse=", "), '), choix="var", new.plot=TRUE, col.var="', Wcol.var, '", col.quanti.sup="', Wcol.quanti.sup, '", label=c("', paste(Wlabel, collapse='", "'), '"), lim.cos2.var=', Wlim.cos, sep="")
      if (is.null(WTitle)) commande.plotVar <- paste(commande.plotVar,')', sep="")
      else {
        if (WTitle ==" ") commande.plotVar <- paste(commande.plotVar,')', sep="")
        else commande.plotVar <- paste(commande.plotVar,', title="', WTitle,'")', sep="")
      }
      justDoIt(commande.plotVar)
      logger(commande.plotVar)
    }}


    # gestion de l'édition de certains resultats                             
    doItAndPrint(paste('summary(',nom.res,', nb.dec = 3, nbelements=10, nbind = 10, ncp = 3, file="")', sep=""))
    if (RFichier==""){                                                       #printing the results
      if(Rpropre) doItAndPrint(paste( nom.res, '$eig', sep=""))
      if(Rvariable) doItAndPrint(paste( nom.res, '$var', sep=""))
      if(Rindividu) doItAndPrint(paste( nom.res, '$ind', sep=""))
      if(Rindsup) doItAndPrint(paste( nom.res, '$ind.sup', sep=""))
      if(Rquantisup) doItAndPrint(paste( nom.res, '$quanti.sup', sep=""))
      if(Rqualisup) doItAndPrint(paste( nom.res, '$quali.sup', sep=""))
      if(Rdescdim) doItAndPrint( paste('dimdesc(', nom.res, ', axes=1:',ncp,')', sep=""))
    }
    else {                                                                  #or saving them in a file
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
    activeDataSet(donnee.depart)                                                    #reloading of the "real" dataset
    justDoIt(paste('remove(',activeDataSet(),'.PCA)',sep=""))                       #the new one is removed
    logger(paste('remove(',activeDataSet(),'.PCA)',sep=""))
  }


    #! fonction associée au bouton OK, execute et détruit l'interface graphique    #function associated with the OK button
  onOK<-function()
  {
    OnAppliquer()                                                                 #OnAppliquer function is performed
    closeDialog(top)                                                              #then window is closed
  }


################################################################################
#                   Création de la fenêtre top                                 #  #creation of the main window
################################################################################
  top<-tktoplevel(borderwidth=10)                                                 #tcltk window
  tkwm.title(top, .Facto_gettext("PCA"))                                            #title
  tkwm.geometry(top, "-100+50")                                                   #size

  # définition des polices
  font2<-tkfont.create(family="times",size=12,weight="bold")                      #font for subtitles
  fontheading<-tkfont.create(family="times",size=18,weight="bold")                #font for titles

  # récupération du jeu de données actif
  donnee<-get(getRcmdr(".activeDataSet"))                                                     #getting active dataset
  vars<-colnames(donnee)                                                          #getting variables names
  rows<-rownames(donnee)                                                          #getting individuals names

  # création de la liste pour le choix des variables acives
  listdesc<-tklistbox(top,selectmode="extended",exportselection="FALSE",yscrollcommand=function(...)tkset(scr,...))                                                                  #empty list box for active variables
  scr <-tkscrollbar(top,repeatinterval=5,command=function(...)tkyview(listdesc,...))  #scroll bar
  listdesc.nom<-NULL
  for (i in (1:ncol(donnee))) {                                                   #putting variables in the list box
      if (is.numeric(donnee[,i])) {
          tkinsert(listdesc,"end",vars[i])
          listdesc.nom<-c(listdesc.nom, vars[i])
      }
  }

  # création de tous les boutons d'options dans IlluFrame
  IlluFrame<- tkframe(top, borderwidth=2)                                         #new frame in the top window
       # mise en page de IlluFrame

  Fillu.funct(label=.Facto_gettext("Supplementary factors"), firstLabel=.Facto_gettext("Supplementary factors"))                                                         #labels for the "select sup fact" button: firstlabel is the default label, label is the label which appears after the selection of sup fact(after clicking on OK)
  Dillu.funct(label=.Facto_gettext("Supplementary quantitative variables"), firstLabel=.Facto_gettext("Supplementary quantitative variables"))                                                       #same for sup var
  Iillu.funct(label=.Facto_gettext("Supplementary individuals"), firstLabel=.Facto_gettext("Supplementary individuals"))                                                     #same for sup ind
  PLOT.PCA(label=.Facto_gettext("Graphical options"), firstLabel=.Facto_gettext("Graphical options"))  
                                                                                  #same for graphical options
  Sortie.funct(label=.Facto_gettext("Outputs"), firstLabel=.Facto_gettext("Outputs")) #same for outputs
  tkgrid(FilluFrame, DilluFrame, IilluFrame, columnspan=7)                        #positionning buttons for supplementary elements
  tkgrid(tklabel(IlluFrame, text=""))                                             #blank line
  tkgrid(PlotFrame, SortieFrame, columnspan=7)                   #buttons for graphical and output options and "reinitialize" button
  tkgrid.configure(FilluFrame, column=1, columnspan=1)                            #configuration of each button
  tkgrid.configure(PlotFrame, column=2, columnspan=2,sticky="we")
  tkgrid.configure(SortieFrame, column=4, columnspan=2,sticky="ew")
  tkgrid.configure(DilluFrame, column=3, columnspan=1)
  tkgrid.configure(IilluFrame, column=5, columnspan=1)
  tkgrid.columnconfigure(IlluFrame,0, minsize=25)                                 #configuration of columns
  tkgrid.columnconfigure(IlluFrame,7, minsize=25)
  tkgrid.columnconfigure(IlluFrame,2, minsize=40)
  tkgrid.columnconfigure(IlluFrame,4, minsize=40)

  # création des options dans OptionFrame                                         #options frame
  OptionFrame<-tkframe(top, borderwidth=2, relief="groove")                       #creation of a new frame
  resu.lab<-tklabel(OptionFrame,text=.Facto_gettext("Name of the result object:")) #label for result's name
  resu.val<-tclVar("res")                                                         #default value for result's name
  resu<-tkentry(OptionFrame,width=10,textvariable=resu.val)                       #entry for result's name
  ncp.lab<-tklabel(OptionFrame,text=.Facto_gettext("Number of dimensions:"))       #label for the choice of the nb of dimensions
  ncp.val<-tclVar("5")                                                            #default value
  ncp<-tkentry(OptionFrame,width=5,textvariable=ncp.val)                          #entry
  reduit.lab<-tklabel(OptionFrame,text=.Facto_gettext("Scale the variables:"))     #label for scaling
  reduit.check <- tkcheckbutton(OptionFrame)                                      #checkbutton
  reduitValue <- tclVar("1")                                                      #default value
  tkconfigure(reduit.check,variable=reduitValue)                                  #configuration
  Axe.label<-tklabel(OptionFrame,text=.Facto_gettext("Graphical output: select the dimensions")) #label for axes
  Axe1<-tclVar("1")                                                               #default value for first axes
  Axe2<-tclVar("2")                                                               #default value for second axes
  Axe1.entry <-tkentry(OptionFrame,width="5",textvariable=Axe1)                   #entry for first axes
  Axe2.entry <-tkentry(OptionFrame,width="5",textvariable=Axe2)                   #entry for second axes

    # mise en page de OptionFrame                                                 #settings
  tkgrid(tklabel(OptionFrame,text=.Facto_gettext("Main options"), fg = "darkred"), columnspan=8, sticky="we")                                                                             #text
  tkgrid(tklabel(OptionFrame,text=""))                                            #blank line
  tkgrid(resu.lab, resu)
  tkgrid(ncp.lab, ncp)
  tkgrid(reduit.lab,reduit.check)
  tkgrid(Axe.label,Axe1.entry , Axe2.entry, sticky="w")

  tkgrid.configure(ncp.lab, reduit.lab, resu.lab, Axe.label, column=1, columnspan=4, sticky="w")  #configuration
  tkgrid.configure(ncp, reduit.check, resu, column=6, columnspan=2, sticky="e")
  tkgrid.configure(Axe1.entry, column=6, columnspan=1, sticky="w")
  tkgrid.configure(Axe2.entry, column=7, columnspan=1, sticky="e")

  tkgrid.columnconfigure(OptionFrame,0, minsize=25)                              #columns' configuration
  tkgrid.columnconfigure(OptionFrame,5, minsize=40)
  tkgrid.columnconfigure(OptionFrame,8, minsize=25)

  #Frame pour HCPC                                                               #HCPC frame
  HcpcFrame<-tkframe(top, borderwidth=2)
  Hcpc.funct(label=.Facto_gettext("Perform Clustering after PCA"), firstLabel=.Facto_gettext("Perform Clustering after PCA"))
  tkgrid(Hcpc2Frame, columnspan=7)
  tkgrid.configure(Hcpc2Frame,column=4, columnspan=1)

  appliquer.but<-tkbutton(top, text=.Facto_gettext("Apply"),width=12,command=OnAppliquer, borderwidth=3, fg="#690f96")                                                                   #"Apply" button
  OKCancelHelp(helpSubject="PCA",reset="Reinitializ.funct")                                                #OK, Cancel and Help button. Help leads to ?PCA

  #TOP                                                                           #settings of the main window
  tkgrid(tklabel(top, text=.Facto_gettext("Principal Component Analysis (PCA)"),font=fontheading),columnspan=3)                                                                   #title
  tkgrid(tklabel(top,text=""))                                                   #blank line
  tkgrid(tklabel(top, text = .Facto_gettext("Select active variables (by default all the variables are active)"),fg = "darkred"), column=1,columnspan=2, sticky = "w")                 #text
  tkgrid(listdesc, scr, sticky = "nw")                                           #list box with available continuous variables and scroll bar
  tkgrid.configure(scr, sticky = "ens",column=2)
  tkgrid.configure(listdesc, sticky = "ew", column=1, columnspan=2)
  tkgrid(tklabel(top,text=""))                                                   #blank line
  tkgrid(IlluFrame, column=1, columnspan=1)                                      #buttons for sup elements, graphical and output options and reinitialization
  tkgrid(tklabel(top,text=""))                                                   #blank line
  tkgrid(OptionFrame, column=1, columnspan=1)                                    #options frame
  tkgrid(tklabel(top,text=""))                                                   #blank line
  tkgrid(HcpcFrame, column=1, columnspan=1)                                      #HCPC frame
  tkgrid(tklabel(top,text=""))                                                   #blank line
  # tkgrid(appliquer.but, column=1, columnspan=1)                                  #"Apply" button
  # tkgrid(tklabel(top,text=""))                                                   #blank line
  # tkgrid(buttonsFrame, column=1, columnspan=1, sticky="ew" )                     #OK, Cancel and Help buttons
  # tkgrid(tklabel(top,text=""))                                                   #blank line
  tkgrid(buttonsFrame, appliquer.but)
  tkgrid.configure(buttonsFrame, column=1,sticky="e")
  tkgrid.configure(appliquer.but, column=2,sticky="w")

}
