selectdirectory <-
function ( )  {

 Envir$selectsave.WD <- tclvalue(tkchooseDirectory(mustexist=TRUE,title="Choose a save directory, then click OK."))
        if(Envir$selectsave.WD == "") {  } 
        else { Envir$save.WD <- Envir$selectsave.WD }
    
##### copier ici a partir de 'Save.directory <- function...' jusqu'a 'open <- tkgrid(HELP4.but...'  #####   

     Save.directory <- function() { selectdirectory() }                                                                                                  # appel la fonction selectdirectory.R
     save.select <- tk2button(Envir$rawdata, text=" Select your save directory ", command=Save.directory)                          # bouton select save path
     tkgrid(save.select, column=1, row=3, sticky="w")

     fixdata1 <- function() { fixdata( ) }                                                                                                                                                        # appel la fonction fixdata.R
     imgFixdata <- tclVar()
     tcl("image","create","photo",imgFixdata,file=file.path(path.package("TTAinterfaceTrendAnalysis"),"aide","imgFixdata.gif",fsep=.Platform$file.sep))
     fixdata.but <- tk2button(Envir$rawdata, image=imgFixdata, text="Fix Data ", compound="right", command=fixdata1, width=8)         # bouton fixdata
     tkgrid(fixdata.but, column=1, row=7, sticky="w")

     showdata <- function()  { showData(Envir$Data) }                                                                                                                                                   # appel la fonction showData
     imgShowdata <- tclVar()
     tcl("image","create","photo",imgShowdata,file=file.path(path.package("TTAinterfaceTrendAnalysis"),"aide","imgShowdata.gif",fsep=.Platform$file.sep))
     showdata <- tk2button(Envir$rawdata, image=imgShowdata, text= "Show Data ", compound="right", command=showdata, width=10)      # bouton showdata
     tkgrid(showdata, column=1, row=5, sticky="w")

#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Frame avec texte variable
     invisibleframe <- tkwidget(Envir$rawdata, "labelframe", text= "")
     tkconfigure(invisibleframe, borderwidth=0)
     tkgrid(invisibleframe, column=1, row=2, columnspan=2, sticky="w")
     Envir$Text <- tk2label(invisibleframe,text= paste("Current active file: "))
     tkgrid(Envir$Text, column=1, row=1, sticky="w")
     Envir$activebox <- tklistbox(invisibleframe, selectmode="single", activestyle="dotbox", height=1, width=70, background="white")
     tkconfigure(Envir$activebox, background="grey95")
     tkgrid(Envir$activebox, column=2, row=1, sticky="w")

     tkdelete(Envir$activebox, 0)
     tkinsert(Envir$activebox, "end", Envir$Name.split)

     tktitle(Envir$tt) <- paste("Temporal Trend Analysis interface (v", Envir$pversion, ") : ", Envir$Name.split, sep="")                   # nom du fichier dans le top panel

     invisibleframe2 <- tkwidget(Envir$rawdata, "labelframe", text= "")
     tkconfigure(invisibleframe2, borderwidth=0)
     tkgrid(invisibleframe2, column=1, row=4, columnspan=2, sticky="w")
     Envir$Text2 <- tk2label(invisibleframe2, text= paste("Current save directory: "))       # texte montrant chemin de sauvegarde
     tkgrid(Envir$Text2, column=1, row=1, sticky="w")
     Envir$savebox <- tklistbox(invisibleframe2, selectmode="single", activestyle="dotbox", height=1, width=66, background="white")
     tkconfigure(Envir$savebox, background="grey95")
     tkgrid(Envir$savebox, column=2, row=1, sticky="w")

     tkdelete(Envir$savebox, 0)
     tkinsert(Envir$savebox, "end", Envir$save.WD)
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ Fin de frame

     STAT1 <- function () { FULLoption( param, depth, sal, site, rawdata="YES")  }                                                              # bouton appelant l'argument rawdata
     STAT1.but <- tk2button(Envir$rawdata, text="Summary ", width=8, command=STAT1)
     tkgrid(STAT1.but, column=1, row=6, sticky="w")

     tkgrid(tklabel(Envir$rawdata, text="  "), column=1, row=8)                                                                # lignes vides (espaces)
     tkgrid(tklabel(Envir$rawdata, text="  "), column=1, row=9)

     AdviceFrame <- tkwidget(Envir$rawdata,"labelframe",text="Important : how to well import your data in 9 steps",padx=30,pady=8, relief = "groove")      # cadre de texte
     tkgrid(AdviceFrame, column=1, row=10, sticky="w")                                                                                  # du premier onglet
     tkgrid(tklabel(AdviceFrame, text="Your data must be in a .txt file (with tab as column separator)"), sticky="w")
     tkgrid(tklabel(AdviceFrame, text="Decimal separtor must be '.'"), sticky="w")
     tkgrid(tklabel(AdviceFrame, text="Missing value must be empty case"), sticky="w")
     tkgrid(tklabel(AdviceFrame, text="Dates must be in  'yyyy-mm-dd' (ISO 8601 time format)"), sticky="w")
     tkgrid(tklabel(AdviceFrame, text=" Dates column -> Dates"), sticky="w")
     tkgrid(tklabel(AdviceFrame, text=" Categorical factors column (Taxa, Stations...) -> Category"), sticky="w")
     tkgrid(tklabel(AdviceFrame, text=" Salinity column -> Salinity"), sticky="w")
     tkgrid(tklabel(AdviceFrame, text=" Depth column -> Depth"), sticky="w")
     tkgrid(tklabel(AdviceFrame, text="If your parameters don't appear, select them as 'numeric' with 'Fix Data'"), sticky="w")
     Adv1 <- tklabel(AdviceFrame, text="You can use the 'Fix Data' button to change column labels and data category")
     tkconfigure(Adv1, foreground="red")
     tkgrid(Adv1, sticky="w")

     tkgrid(tklabel(Envir$rawdata, text="      "), column=1, row=12)

     HELP1.but <- tk2button(Envir$HelpFrame, image=Envir$imgHelp, text=" Help ", compound="right", command=function() { Aide1() })                       # fichier txt a aller chercher
     tkgrid(HELP1.but, column=2, row=1, sticky="nw")

#_______________________________________________________________________________________________________________________Onglet Selection des parametres

#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _Liste de category                                                                                                                    # cree un environnement pour les listes
     st <- levels(as.factor(Envir$Data$Category))                                                                            # label des SERIES
     Env$variables <- c("-All-", st)                                                                                         # on ajoute '-All-' a la liste
     Env$variables.selectionnees.temp <- integer()

     liste1 <- tklistbox(Envir$Select,selectmode="extended", activestyle="dotbox", height=10, width=27, background="white")       # cree la premiere liste
                for (i in 1:length(Env$variables)) { tkinsert(liste1,"end",Env$variables[i]) }                                                                         # rempli la premiere liste avec le nom des SERIES
     liste2 <- tklistbox(Envir$Select,selectmode="extended", activestyle="dotbox", height=10, width=27, background="white")       # cree la deuxieme liste (de selection)

     imgArrowright <- tclVar()
     tcl("image","create","photo",imgArrowright,file=file.path(path.package("TTAinterfaceTrendAnalysis"),"aide","imgArrowright.gif",fsep=.Platform$file.sep))
     imgArrowleft <- tclVar()
     tcl("image","create","photo",imgArrowleft,file=file.path(path.package("TTAinterfaceTrendAnalysis"),"aide","imgArrowleft.gif",fsep=.Platform$file.sep))
     imgArrowright2 <- tclVar()
     tcl("image","create","photo",imgArrowright2,file=file.path(path.package("TTAinterfaceTrendAnalysis"),"aide","imgArrowright2.gif",fsep=.Platform$file.sep))

     bouton1 <- tk2button(Envir$Select,image=imgArrowright, command=function() {                                               # bouton de selection
                if (tclvalue(tkcurselection(liste1))!="") {
                      val <- unlist(strsplit(tclvalue(tkcurselection(liste1)), "\\ "))                                         # les SERIES selectionnees sont stockees dans 'val'
                      selection <- as.numeric(val)+1                                                                           # puis sont transformees (+1 car la premiere valeur de la liste = 0)
                      for (i in min(selection):max(selection)) { tkinsert(liste2,"end",Env$variables[i]) }                     # insertion des valeurs selectionnees dans la liste 2
                      Env$variables.selectionnees.temp <- c(Env$variables.selectionnees.temp, selection) }                     # les valeurs finales selectionnees sont stockees
                else { tkmessageBox(message="No categorical factor selected !",type="ok",icon="info", title="!Warning!") }})
     bouton2 <- tk2button(Envir$Select,image=imgArrowleft, command=function() {                                                # bouton de deselection
                if (tclvalue(tkcurselection(liste2))!="") {
                      val <- unlist(strsplit(tclvalue(tkcurselection(liste2)), "\\ "))
                      selection <- as.numeric(val)+1
                      for (i in 1:length(selection)) { tkdelete(liste2, min(selection)-1) }
                      for (i in 1:length(selection)) { Env$variables.selectionnees.temp <- Env$variables.selectionnees.temp[-min(selection)]  }}
                else { tkmessageBox(message="No categorical factor selected !",type="ok",icon="info", title="!Warning!") } })

     tkgrid(tklabel(Envir$Select, text="Select your categorical factor(s)"), row=0, column=0,rowspan=2)                        # titres des listes
     tkgrid(tklabel(Envir$Select, text="Selected categorical factor(s)"), row=0, column=2,rowspan=2)                           # mise en page...
     tkgrid(liste1,row=2,column=0,rowspan=2)
     tkgrid(bouton1,row=2,column=1)
     tkgrid(bouton2,row=3,column=1)
     tkgrid(liste2,row=2,column=2,rowspan=2)
     tkgrid(tklabel(Envir$Select,text=" "), column=0, row=4)

#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _Liste des parametres
     ListParam <- NULL                                                                                                                    # cree une liste vide a remplir avec les parametres present dans le csv
     for (i in 1:ncol(Envir$Data)) {                                                                                                      # i entre 1 et le nbre de colonne
     if (class(Envir$Data[ ,i]) == "numeric") {                                                                                           # pour les colonnes dans lesquelles les valeurs sont 'numeric'
     d <- names(Envir$Data[i])                                                                                                            # on garde le nom de ces colonnes
     ListParam <- as.character(c(ListParam, d)) } }                                                                                       #    et on en fait la liste des parametres
                                                                                                                                          # on cree les listes comme pour les SERIES
     Env2$variables <- ListParam
     Env2$variables.selectionnees.temp <- integer()

     liste3 <- tklistbox(Envir$Select,selectmode="single", activestyle="dotbox", height=6, width=27, background="white")
     tkdelete(liste3, "0", "end")
               for (i in 1:length(Env2$variables)) { tkinsert(liste3,"end",Env2$variables[i]) }
     liste4 <- tklistbox(Envir$Select,selectmode="single", activestyle="dotbox", height=6, width=27, background="white")

     bouton3 <- tk2button(Envir$Select, image=imgArrowright2, command=function() {
                if (tclvalue(tkcurselection(liste3))!="") {
                   val2 <- unlist(strsplit(tclvalue(tkcurselection(liste3)), "\\ "))
                   selection2 <- as.numeric(val2)+1
                   for (i in min(selection2):max(selection2)) {
                   tkdelete(liste4, 0)                                                                                # supprime la valeur precedemment presente dans la liste 2
                   tkinsert(liste4,"end",Env2$variables[i]) }                                   #     et la remplace par la nouvelle valeur
                   Env2$variables.selectionnees.temp <- c(selection2) }
                else { tkmessageBox(message="No variable selected !",type="ok",icon="info", title="!Warning!") }})

     tkgrid(tklabel(Envir$Select, text="Select a parameter"), row=11, column=0,rowspan=2)            # titre des listes
     tkgrid(tklabel(Envir$Select, text="Selected parameter"), row=11, column=2,rowspan=2)            #       "
     tkgrid(liste3, row=13, column=0)
     tkgrid(bouton3, row=13, column=1)
     tkgrid(liste4, row=13, column=2)

     tkgrid(tklabel(Envir$Select, text="      "), column=0, row=15)
     tkgrid(tklabel(Envir$Select, text="      "), column=2, row=15)

     if (class(Env2$variables[1]) != "character")   {
         editPopupMenu <- tkmenu(liste3, tearoff = FALSE)
         tkadd(editPopupMenu, "command", label = "A parameter doesn't appear? -> Click and change the parameter(s) type to numeric", command= fixdata1)
	       tkconfigure( editPopupMenu, foreground="steelblue4")
	       LeftClick <- function(x, y) {
         rootx <- as.integer(tkwinfo("rootx", liste3))
         rooty <- as.integer(tkwinfo("rooty", liste3))
         xTxt <- as.integer(x) + rootx
         yTxt <- as.integer(y) + rooty
         .Tcl(paste("tk_popup", .Tcl.args(editPopupMenu, xTxt, yTxt))) }
         tkbind(liste3, "<Button-3>", LeftClick)	}
     else { if (length(Env2$variables) == ncol(Envir$Data)-2  ) {} else {
	          editPopupMenu <- tkmenu(liste3, tearoff = FALSE)
            tkadd(editPopupMenu, "command", label = "A parameter doesn't appear? -> Click and change the parameter(s) type to numeric", command= fixdata1)
	          tkconfigure( editPopupMenu, foreground="steelblue4")
	          LeftClick <- function(x, y) {
            rootx <- as.integer(tkwinfo("rootx", liste3))
            rooty <- as.integer(tkwinfo("rooty", liste3))
            xTxt <- as.integer(x) + rootx
            yTxt <- as.integer(y) + rooty
            .Tcl(paste("tk_popup", .Tcl.args(editPopupMenu, xTxt, yTxt))) }
            tkbind(liste3, "<Button-3>", LeftClick)	}}

#_______________________________________________________________________________________________________________________Sliders salinite et profondeur
#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _Sliders pour les salinites
     if (any(colnames(Envir$Data) == "Salinity") & any(!is.na(Envir$Data$Salinity))==TRUE) {                                                     # si une colonne S est presente
        if (class(Envir$Data$Salinity) == "numeric") {                                                                                                                                  #      et si c'est bien une valeur numerique  :
          sal1 <- tclVar(min(Envir$Data$Salinity, na.rm=TRUE))                                                                                                                    # prend la valeur mini de la salinite
          slider1 <- tk2scale(Envir$Select, from=min(Envir$Data$Salinity, na.rm=TRUE), to=max(Envir$Data$Salinity, na.rm=TRUE),        # slider coulissant entre...
                     variable=sal1,  orientation="horizontal", command=function(...) {                                                                                              #    ...les valeurs min et max de la colonne S
                     tkconfigure(Label1,text=paste(" Salinity min :", format(floor(as.numeric(tclvalue(sal1))*2)/2, nsmall=1), "psu ")) })                      # titre des sliders
          Label1 <- tklabel(Envir$Select,text=paste(" Salinity min :", format(floor(as.numeric(tclvalue(sal1))*2)/2, nsmall=1), "psu "))
          tkgrid(Label1, column=0, row=16)
          tkgrid(slider1, column=0, row=17)

          sal2 <- tclVar(max(Envir$Data$Salinity, na.rm=TRUE))                                                                                                                    # idem pour la valeur de salinite max
          slider2 <- tk2scale(Envir$Select, from=min(Envir$Data$Salinity, na.rm=TRUE), to=max(Envir$Data$Salinity, na.rm=TRUE),
                     variable=sal2,
                     orientation="horizontal", command=function(...) {
                     tkconfigure(Label2,text=paste(" Salinity max :", format(ceiling(as.numeric(tclvalue(sal2))*2)/2, nsmall=1), "psu ")) })
          Label2 <- tklabel(Envir$Select,text=paste(" Salinity max :", format(ceiling(as.numeric(tclvalue(sal2))*2)/2, nsmall=1), "psu "))
          tkgrid(Label2, column=0, row=18)
          tkgrid(slider2, column=0, row=19) }

     else{ tkmessageBox(message=paste("A salinity column is present but don't contain any 'numeric' values", "\n\n"   # message d'erreur si probleme de formattage...
                         , "You can fix it using the fix data button (use Help for more info)", sep="")                                                #   ...des donnees dans la colonne S
                         , icon = "warning", type = "ok", title="!Warning!") } }

     else {sal1 <- tclVar(0)                                                                                          # s'il n'y a pas de colonne Salinity
          slider1 <- tk2scale(Envir$Select, from=0, to=0,
                     variable=sal1,  orientation="horizontal", command=function(...) {
                     tkconfigure(Label1,text=paste("   No Salinity column !   ")) } )
          Label1 <- tklabel(Envir$Select,text=paste("   No Salinity column !   "))
          tkgrid(Label1, column=0, row=16)
          tkgrid(slider1, column=0, row=17)

          sal2 <- tclVar(0)
          slider2 <- tk2scale(Envir$Select, from=0, to=0,
                     variable=sal2,
                     orientation="horizontal", command=function(...) {
                     tkconfigure(Label2,text=paste("   No Salinity column !   ")) })
          Label2 <- tklabel(Envir$Select,text=paste("   No Salinity column !   "))
          tkgrid(Label2, column=0, row=18)
          tkgrid(slider2, column=0, row=19) }

     tkgrid(tklabel(Envir$Select, text="                                  "), column=0, row=20)

#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _Sliders pour les profondeurs
     if (any(colnames(Envir$Data) == "Depth") & any(!is.na(Envir$Data$Depth))==TRUE ) {
       if (class(Envir$Data$Depth) == "numeric" ) {
          depth1 <- tclVar(min(Envir$Data$Depth, na.rm=TRUE))
          slider3 <- tk2scale(Envir$Select, from=min(Envir$Data$Depth, na.rm=TRUE), to=max(Envir$Data$Depth, na.rm=TRUE),
                     variable=depth1, orientation="horizontal", command=function(...) {
                     tkconfigure(Label3,text=paste("Depth min :", format(floor(as.numeric(tclvalue(depth1))*2)/2, nsmall=1), "m")) })
          Label3 <- tklabel(Envir$Select,text=paste("Depth min :", format(floor(as.numeric(tclvalue(depth1))*2)/2, nsmall=1), "m"))
          tkgrid(Label3, column=0, row=21)
          tkgrid(slider3, column=0, row=22)

          depth2 <- tclVar(max(Envir$Data$Depth, na.rm=TRUE))
          slider4 <- tk2scale(Envir$Select, from=min(Envir$Data$Depth, na.rm=TRUE), to=max(Envir$Data$Depth, na.rm=TRUE),
                     variable=depth2,  orientation="horizontal", command=function(...) {
                     tkconfigure(Label4,text=paste("Depth max :", format(ceiling(as.numeric(tclvalue(depth2))*2)/2, nsmall=1), "m")) })
          Label4 <- tklabel(Envir$Select,text=paste("Depth max :", format(ceiling(as.numeric(tclvalue(depth2))*2)/2, nsmall=1), "m"))
          tkgrid(Label4, column=0, row=23)
          tkgrid(slider4, column=0, row=24)  }

       else{ tkmessageBox(message=paste("A depth column is present but don't contain any 'numeric' values", "\n\n"
                         , "You can fix it using the fix data button (use Help for more info)", sep="")
                         , icon = "warning", type = "ok", title="!Warning!") } }

       else { depth1 <- tclVar(0)
          slider3 <- tk2scale(Envir$Select, from=0, to=0,
                     variable=depth1, orientation="horizontal", command=function(...) {
                     tkconfigure(Label3,text=paste("   No Depth column !   ")) })
          Label3 <- tklabel(Envir$Select,text=paste("   No Depth column !   "))
          tkgrid(Label3, column=0, row=21)
          tkgrid(slider3, column=0, row=22)

          depth2 <- tclVar(0)
          slider4 <- tk2scale(Envir$Select, from=0, to=0,
                     variable=depth2,  orientation="horizontal", command=function(...) {
                     tkconfigure(Label4,text=paste("   No Depth column !   ")) })
          Label4 <- tklabel(Envir$Select,text=paste("   No Depth column !   "))
          tkgrid(Label4, column=0, row=23)
          tkgrid(slider4, column=0, row=24)  }

#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _Choix des dates (annees et mois)
     if (any(colnames(Envir$Data) == "Dates")) {                                                                                                         # idem pour les annees

            if (is.numeric(Envir$Data$Dates) == TRUE ) { years <- Envir$Data$Dates } else {
            d <- as.Date(Envir$Data$Dates, format="%Y-%m-%d")                                                                                      # formattage de la date
            years <- as.numeric(format(d, format = "%Y")) }                                                                                             # on recupere les annees

            if ( is.na(d) == FALSE ) {} else { return(tkmessageBox(message="Wrong date format: must be 'yyy-mm-dd'"
                                                                   ,icon = "warning", type = "ok", title="!Warning!")) }

            year1 <- tclVar(min(years, na.rm=TRUE))                                                                                                          # boite de texte avec fleches
            SpinBox1 <- tk2spinbox(Envir$Select
                                 ,from=min(years, na.rm=TRUE)
                                 ,to=max(years, na.rm=TRUE)
				                         ,relief="groove"                                                                                                                                # aspect de la boite "groove" "flat" "solid"
                                 ,readonlybackground="white"
                                 ,width=16
                                 ,textvariable=year1
				                         ,background="white")
            tkgrid(tklabel(Envir$Select, text="       Start year at :       "), column=2, row=16)
            tkgrid(SpinBox1, column=2, row=17)

            year2 <- tclVar(max(years, na.rm=TRUE))
            SpinBox2 <- tk2spinbox(Envir$Select
                                 ,from=min(years, na.rm=TRUE)
                                 ,to=max(years, na.rm=TRUE)
				                         ,relief="groove"
                                 ,readonlybackground="white"
                                 ,width=16
                                 ,textvariable=year2
				                         ,background="white")
           tkgrid(tklabel(Envir$Select, text="        End year at :        "), column=2, row=18)
           tkgrid(SpinBox2, column=2, row=19)

           emptybox20 <- tklabel(Envir$Select, text="                                          ")
           tkgrid(emptybox20, column=2, row=20)

           if (is.numeric(Envir$Data$Dates) == TRUE) { mois <- NULL } else {

           m <- as.factor(format(d, format = "%m"))
           mois <- tclVar(as.numeric(levels(m)))                                                                                     # uniquement les mois presents dans la base de donnees (separes par un espace)
           entry.mois <-tk2entry(Envir$Select,width="21",textvariable=mois , background="white")
           tkgrid(tklabel(Envir$Select,text="     Select month(s)     "), column=2, row=21)
           tkgrid(entry.mois, column=2, row=22) } }

     else {tkgrid(tklabel(Envir$Select, text="                      "), column=2, row=16)
	         emptybox10 <- tklabel(Envir$Select, text="!No date in your data!")                       # s'il n'y a pas de colonne Dates
           emptybox30 <- tklabel(Envir$Select, text="  !Or wrong format!  ")
	         tkconfigure(emptybox10, foreground="red")
	         tkconfigure(emptybox30, foreground="red")
           tkgrid(emptybox10, column=2, row=17)
	         tkgrid(emptybox30, column=2, row=18)
           tkgrid(tklabel(Envir$Select, text="                                    "), column=2, row=19)
           tkgrid(tklabel(Envir$Select, text="   Check the column   "), column=2, row=20)
           tkgrid(tklabel(Envir$Select, text="  label and the Date   "), column=2, row=21)
           tkgrid(tklabel(Envir$Select, text="                  format                  "), column=2, row=22)
           tkgrid(tklabel(Envir$Select, text="                                          "), column=2, row=23) }

#_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _champs separateur
	separator1 <- tk2label(Envir$Select, text=" | ", foreground="steelblue3")
   tkgrid(separator1, column=1, row=17)
	separator2 <- tk2label(Envir$Select, text=" | ", foreground="steelblue3")
   tkgrid(separator2, column=1, row=18)
  separator3 <- tk2label(Envir$Select, text=" | ", foreground="steelblue3")
   tkgrid(separator3, column=1, row=19)
	separator4 <- tk2label(Envir$Select, text=" | ", foreground="steelblue3")
   tkgrid(separator4, column=1, row=20)
  separator5 <- tk2label(Envir$Select, text=" | ", foreground="steelblue3")
   tkgrid(separator5, column=1, row=21)
  separator6 <- tk2label(Envir$Select, text=" | ", foreground="steelblue3")
   tkgrid(separator6, column=1, row=22)
  separator7 <- tk2label(Envir$Select, text=" | ", foreground="steelblue3")
   tkgrid(separator7, column=1, row=23)
  separator8 <- tk2label(Envir$Select, text=" | ", foreground="steelblue3")
   tkgrid(separator8, column=1, row=24)

#_______________________________________________________________________________bouton d'appel de la fonction 'select' dans FULLoption
     STAT2 <- function () {                                                                                                                # on rentre tout les arguments precedemment selectionne dans la fonction STAT2
            param <- Env2$variables[unique(Env2$variables.selectionnees.temp)]               # les parametres
            site <- Env$variables[unique(Env$variables.selectionnees.temp)]                     # les category (anciennement les SITES)
            if (any(site == "-All-")) { site <- st  }                                                                        # si ALL est selectionnee
            if (length(param)==0 | length(site)==0 )
               { return(tkmessageBox(message="Please select a parameter and a category before to proceed"      # erreur si aucun argument n'est selectionnede
                                     , icon = "warning", type = "ok", title="!Warning!")) }
            else{}
           mix <- "YES"
           if (any(colnames(Envir$Data) == "Salinity") & any(!is.na(Envir$Data$Salinity))==TRUE) {                     # valeurs de salinites
            Sal1 <- floor(as.numeric(tclvalue(sal1))*2)/2
            Sal2 <- ceiling(as.numeric(tclvalue(sal2))*2)/2
            sal <- c(Sal1, Sal2) }
           else{ sal <- "NULL" }
           if (any(colnames(Envir$Data) == "Depth") & any(!is.na(Envir$Data$Depth)) == TRUE) {                               # valeurs de profondeurs
            Depth1 <- floor(as.numeric(tclvalue(depth1))*2)/2
            Depth2 <- ceiling(as.numeric(tclvalue(depth2))*2)/2
            depth <- c(Depth1, Depth2) }
           else{ depth <- "NULL" }
           if (any(colnames(Envir$Data) == "Dates")) {
            start <- as.numeric(tclvalue(year1))                                                                                                                           # annee de debut
            end <- as.numeric(tclvalue(year2))                                                                                                                               # annee de fin
            if (is.numeric(Envir$Data$Dates) ==TRUE) { } else {
            months <- as.numeric(unlist(strsplit((tclvalue(mois)),"\\ "))) } }                                                                  # mois
           else { tkmessageBox(message="No date selected!", icon = "warning", type = "ok", title="!Warning!")}        # erreur si pas de dates
           FULLoption(param, depth, sal, site, rawdata="NO", select="YES", resume.reg="NO", test.normality="NO",       # applique la fonction FULLoption...
                     plotB="NO", selectBox, log.trans="NO", plotZ="NO", datashow="NO",                                                               #    ...avec les arguments
                     help.timestep="NO", auto.timestep="NO", time.step="NULL", help.aggreg="NO", auto.aggreg="NO", aggreg="NULL",
                     mix, outliers.re="NO", na.replace="NO", start, end, months)   }
     STAT2.but <- tk2button(Envir$Select, text="Summary",command=STAT2)                                                                           # bouton d'appel de la fonction STAT2
     tkgrid(STAT2.but, column=2, row=24)
#_______________________________________________________________________________fin du bouton

     tkgrid(tklabel(Envir$Select, text="      "), column=0, row=26)
     tkgrid(tklabel(Envir$Select, text="      "), column=0, row=27)

     HELP2.but <- tk2button(Envir$Select, image=Envir$imgHelp, text=" Help ", compound="right",command=function() { Aide2() } )     # bouton d'aide 2
     tkgrid(HELP2.but, column=0, row=28, sticky="w")

#_________________________________________________________________________________________________________________________________Onglet regularisation

     LabeledFrame1 <- tkwidget(Envir$datam,"labelframe",text="Data interactions",padx=25,pady=10)                      # cadre
     tkgrid(LabeledFrame1, column=0, row=0, sticky="w")

         cb7 <- tk2checkbutton(LabeledFrame1)                                                                                                                                # check button pour la transformation log
         cb7Value <- tclVar("0")
         tkconfigure(cb7,variable=cb7Value)
         tkgrid(tklabel(LabeledFrame1,text="Log10(x+1) transform"), column=0, row=0, sticky="w")
         tkgrid(cb7, column=1, row=0)

         cb2 <- tk2checkbutton(LabeledFrame1)                                                                                                                               # check button pour missing values
         cb2Value <- tclVar("0")
         tkconfigure(cb2,variable=cb2Value)
         tkgrid(tklabel(LabeledFrame1,text="Replace missing values ?"), column=0, row=1, sticky="w")
         tkgrid(cb2, column=1, row=1)

         cb3 <- tk2checkbutton(LabeledFrame1)                                                                                                                               # check button pour outliers
         cb3Value <- tclVar("0")
         tkconfigure(cb3,variable=cb3Value)
         tkgrid(tklabel(LabeledFrame1,text="Remove outliers ?"), column=0, row=2, sticky="w")
         tkgrid(cb3, column=1, row=2)

         selectBox<-NULL                                                                                                                                                                    # radio button pour les boxplots
         rb21 <- tk2radiobutton(LabeledFrame1)
         rb22 <- tk2radiobutton(LabeledFrame1)
         rb20Value <- tclVar("ByYears")
         tkconfigure(rb21,variable=rb20Value,value="ByYears")
         tkconfigure(rb22,variable=rb20Value,value="ByMonths")
         tkgrid(tklabel(LabeledFrame1,text="   |->By years "),rb21, row = 4, sticky="w")
         tkgrid(tklabel(LabeledFrame1,text="   |->By months "),rb22, row= 5, sticky="w")

#____________________________________________________________________________________________bouton d'appel de la boxplot (argument plotB)
     BoxPlot <- function()  {
            param <- Env2$variables[unique(Env2$variables.selectionnees.temp)]
            site <- Env$variables[unique(Env$variables.selectionnees.temp)]
            if (any(site == "-All-")) { site <- st  }
            if (length(param)==0 | length(site)==0 )
               { return(tkmessageBox(message="Please select a parameter and a category before to proceed"
                                    , icon = "warning", type = "ok", title="!Warning!")) }
            else{}
           mix <- "YES"
           if (any(colnames(Envir$Data) == "Salinity") & any(!is.na(Envir$Data$Salinity))==TRUE) {                     # valeurs de salinites
            Sal1 <- floor(as.numeric(tclvalue(sal1))*2)/2
            Sal2 <- ceiling(as.numeric(tclvalue(sal2))*2)/2
            sal <- c(Sal1, Sal2) }
           else{ sal <- "NULL" }
           if (any(colnames(Envir$Data) == "Depth") & any(!is.na(Envir$Data$Depth)) == TRUE) {                               # valeurs de profondeurs
            Depth1 <- floor(as.numeric(tclvalue(depth1))*2)/2
            Depth2 <- ceiling(as.numeric(tclvalue(depth2))*2)/2
            depth <- c(Depth1, Depth2) }
           else{ depth <- "NULL" }
           if (any(colnames(Envir$Data) == "Dates")) {
            start <- as.numeric(tclvalue(year1))
            end <- as.numeric(tclvalue(year2))
            if (is.numeric(Envir$Data$Dates) ==TRUE) { } else {
            months <- as.numeric(unlist(strsplit((tclvalue(mois)),"\\ "))) }  }
           else { tkmessageBox(message="no date selected", icon = "warning", type = "ok", title="!Warning!")}
           rb20Value <- as.character(tclvalue(rb20Value))
           if (rb20Value=="ByYears"){ selectBox <- "ByYears" }
           if (rb20Value=="ByMonths"){ selectBox <- "ByMonths" }
           FULLoption(param, depth, sal, site, rawdata="NO", select="NO", resume.reg="NO",test.normality="NO",
                     plotB="YES", selectBox, log.trans="NO", plotZ="NO", datashow="NO",
                     help.timestep="NO", auto.timestep="NO", time.step="NULL", help.aggreg="NO", auto.aggreg="NO", aggreg="NULL",
                     mix, outliers.re="NO", na.replace="NO", start, end, months)  }
     BoxPlot.but <- tk2button(LabeledFrame1, text=" Show boxplot ",command=BoxPlot)
     tkgrid(BoxPlot.but, column=0, row=3, sticky="w")
#_________________________________________________________________________________________________fin du bouton Boxplot

     tkgrid(tklabel(Envir$datam, text="              "), column=0, row=1)

     LabeledFrame2 <- tkwidget(Envir$datam,"labelframe",text="Select the data frequency in your final time series",padx=25,pady=10, relief = "groove")
     tkgrid(LabeledFrame2, column=0, row=2, sticky="w")

         rb1 <- tk2radiobutton(LabeledFrame2)                                                                     # radio buttons pour le time step
         rb2 <- tk2radiobutton(LabeledFrame2)
         rb3 <- tk2radiobutton(LabeledFrame2)
         rb4 <- tk2radiobutton(LabeledFrame2)
         rb5 <- tk2radiobutton(LabeledFrame2)
         rb6 <- tk2radiobutton(LabeledFrame2)
         rb7 <- tk2radiobutton(LabeledFrame2)
         rb8 <- tk2radiobutton(LabeledFrame2)
         if (is.numeric(Envir$Data$Dates) ==TRUE) { rb1Value <- tclVar("Annual") } else {
         rb1Value <- tclVar("auto") }
         tkconfigure(rb1,variable=rb1Value,value="Annual")                                         # valeur de l'argument pour chaque bouton
         tkconfigure(rb2,variable=rb1Value,value="Monthly")
         tkconfigure(rb3,variable=rb1Value,value="Fortnight")
         tkconfigure(rb4,variable=rb1Value,value="Semi-fortnight")
         tkconfigure(rb5,variable=rb1Value,value="Mono-mensual")
         tkconfigure(rb6,variable=rb1Value,value="Daily")
         tkconfigure(rb7,variable=rb1Value,value="help")
         tkconfigure(rb8,variable=rb1Value,value="auto")
         tkgrid(tklabel(LabeledFrame2,text="Daily "),rb6, sticky="w")                     # texte afficher a cote de chaque bouton
         tkgrid(tklabel(LabeledFrame2,text="Semi-fortnightly "),rb4, sticky="w")
         tkgrid(tklabel(LabeledFrame2,text="Fortnightly "),rb3, sticky="w")
         tkgrid(tklabel(LabeledFrame2,text="Monthly "),rb2, sticky="w")
         tkgrid(tklabel(LabeledFrame2,text="Yearly "),rb1, sticky="w")
         tkgrid(tklabel(LabeledFrame2,text="Monthly - Climato (see help)"),rb5, sticky="w")
         tkgrid(tklabel(LabeledFrame2,text="Guidance to choose the frequency "),rb7, sticky="w")
         tkgrid(tklabel(LabeledFrame2,text="Auto "),rb8, sticky="w")

     tkgrid(tklabel(Envir$datam, text="              "), column=0, row=3)

     LabeledFrame3 <- tkwidget(Envir$datam,"labelframe",text="Select the method to aggregate your data",padx=25,pady=10, relief = "groove")
     tkgrid(LabeledFrame3, column=0, row=4, sticky="w")

         rb1 <- tk2radiobutton(LabeledFrame3)                                                                      # radio buttons pour la methode mathematique
         rb2 <- tk2radiobutton(LabeledFrame3)
         rb3 <- tk2radiobutton(LabeledFrame3)
         rb4 <- tk2radiobutton(LabeledFrame3)
         rb5 <- tk2radiobutton(LabeledFrame3)
         rb6 <- tk2radiobutton(LabeledFrame3)
         if (is.numeric(Envir$Data$Dates) ==TRUE) { rb2Value <- tclVar("Mean") } else {
         rb2Value <- tclVar("auto") }
         tkconfigure(rb1,variable=rb2Value,value="Mean")
         tkconfigure(rb2,variable=rb2Value,value="Median")
         tkconfigure(rb3,variable=rb2Value,value="Quantile")
         tkconfigure(rb4,variable=rb2Value,value="Max")
         tkconfigure(rb5,variable=rb2Value,value="help")
         tkconfigure(rb6,variable=rb2Value,value="auto")
         tkgrid(tklabel(LabeledFrame3,text="Mean "),rb1, sticky="w")
         tkgrid(tklabel(LabeledFrame3,text="Median "),rb2, sticky="w")
         tkgrid(tklabel(LabeledFrame3,text="Quantile 0.9 "),rb3, sticky="w")
         tkgrid(tklabel(LabeledFrame3,text="Maximum "),rb4, sticky="w")
         tkgrid(tklabel(LabeledFrame3,text="Guidance to choose the method      "),rb5, sticky="w")
         tkgrid(tklabel(LabeledFrame3,text="Auto "),rb6, sticky="w")

         tkgrid(tklabel(Envir$datam, text="           "), column=2, row=2)

     LabeledFrame7 <- tkwidget(Envir$datam,"labelframe",text="Show regularised time series", padx=10, pady=0, relief = "groove")
     tkgrid(LabeledFrame7, column=3, row=2, sticky="n")
     tkgrid(tklabel(LabeledFrame7, text="                                                                                                                                                                            "
                                 , font=tkfont.create(size=1)), column=0, row=0, columnspan=3)
     tkgrid(tklabel(LabeledFrame7, text=" ", font=tkfont.create(size=1)), column=0, row=2)

#_______________________________________________________________________________bouton d'appel de la figure de la serie regularisee (plotZ)
     OnOK <- function()  {
            param <- Env2$variables[unique(Env2$variables.selectionnees.temp)]
            site <- Env$variables[unique(Env$variables.selectionnees.temp)]
            if (any(site == "-All-")) { site <- st  }
            if (length(param)==0 | length(site)==0 )
               { return(tkmessageBox(message="Please select a parameter and a category before to proceed"
                                    , icon = "warning", type = "ok", title="!Warning!")) }
            else{}
           mix <- "YES"
           if (any(colnames(Envir$Data) == "Salinity") & any(!is.na(Envir$Data$Salinity))==TRUE) {                     # valeurs de salinites
            Sal1 <- floor(as.numeric(tclvalue(sal1))*2)/2
            Sal2 <- ceiling(as.numeric(tclvalue(sal2))*2)/2
            sal <- c(Sal1, Sal2) }
           else{ sal <- "NULL" }
           if (any(colnames(Envir$Data) == "Depth") & any(!is.na(Envir$Data$Depth)) == TRUE) {                               # valeurs de profondeurs
            Depth1 <- floor(as.numeric(tclvalue(depth1))*2)/2
            Depth2 <- ceiling(as.numeric(tclvalue(depth2))*2)/2
            depth <- c(Depth1, Depth2) }
           else{ depth <- "NULL" }
           if (any(colnames(Envir$Data) == "Dates")) {
            start <- as.numeric(tclvalue(year1))
            end <- as.numeric(tclvalue(year2))
            if (is.numeric(Envir$Data$Dates) ==TRUE) { } else {
            months <- as.numeric(unlist(strsplit((tclvalue(mois)),"\\ "))) } }
           else { tkmessageBox(message="no date selected", icon = "warning", type = "ok", title="!Warning!")}

          cb7Value <- as.character(tclvalue(cb7Value))                           # valeur de l'argument log transform
          if (cb7Value=="1"){ log.trans <- "YES" }
           else { log.trans <- "NO" }
          cb2Value <- as.character(tclvalue(cb2Value))                           # valeur de l'argument na.replace
          if (cb2Value=="1"){ na.replace <- "YES" }
           else { na.replace <- "NO" }
          cb3Value <- as.character(tclvalue(cb3Value))                           # valeur de l'argument outliers.re
          if (cb3Value=="1"){ outliers.re <- "YES" }
           else { outliers.re <- "NO" }

          rb1Value <- as.character(tclvalue(rb1Value))
          if (rb1Value=="Annual"){ time.step <- "Annual" }                       # valeur de l'argument time.step pour chaque bouton
          if (rb1Value=="Monthly"){ time.step <- "Monthly" }
          if (rb1Value=="Fortnight"){ time.step <- "Fortnight" }
          if (rb1Value=="Semi-fortnight"){ time.step <- "Semi-fortnight" }
          if (rb1Value=="Mono-mensual"){ time.step <- "Mono-mensual" }
          if (rb1Value=="Daily"){ time.step <- "Daily" }
          if (rb1Value=="help"){ help.timestep <- "YES"
                                 time.step <- "NULL" }
          else{ help.timestep <- "N0" }
          if (rb1Value=="auto"){ auto.timestep <- "YES"
                                 time.step <- "NULL" }
          else{ auto.timestep <- "NO" }

          rb2Value <- as.character(tclvalue(rb2Value))
          if (rb2Value=="Mean"){ aggreg <- "Mean" }                              # valeur de l'argument aggreg pour chauqe bouton
          if (rb2Value=="Median"){ aggreg <- "Median" }
          if (rb2Value=="Quantile"){ aggreg <- "Quantile" }
          if (rb2Value=="Max"){ aggreg <- "Max" }
          if (rb2Value=="help"){ help.aggreg <- "YES"
                                 aggreg <- "NULL" }
          else{ help.aggreg <- "N0" }
          if (rb2Value=="auto"){ auto.aggreg <- "YES"
                                 aggreg <- "NULL" }
          else{ auto.aggreg <- "N0" }
          FULLoption(param, depth, sal, site, rawdata="NO", select="NO", resume.reg="NO", test.normality="NO",
                 plotB="NO", selectBox, log.trans, plotZ="YES", datashow="NO",
                 help.timestep, auto.timestep, time.step, help.aggreg, auto.aggreg, aggreg,
                 mix, outliers.re, na.replace, start, end, months, norm="NO", npsu,
                 autocorr = "NO", spectrum="NO", anomaly="NO", a.barplot="NO", zsmooth="NO", local.trend="NO", test="NO") }
     OK.but <- tk2button(LabeledFrame7, text="  Plot  ",command=OnOK, width=9)
     tkgrid(OK.but, column=0, row=1, sticky="w")
#______________________________________________________________________________________________bouton d'appel du tableau des donnees regularisees
     OnOK1 <- function()  {
            param <- Env2$variables[unique(Env2$variables.selectionnees.temp)]
            site <- Env$variables[unique(Env$variables.selectionnees.temp)]
            if (any(site == "-All-")) { site <- st  }
            if (length(param)==0 | length(site)==0 )
               { return(tkmessageBox(message="Please select a parameter and a category before to proceed"
                                     , icon = "warning", type = "ok", title="!Warning!")) }
            else{}
           mix <- "YES"
           if (any(colnames(Envir$Data) == "Salinity") & any(!is.na(Envir$Data$Salinity))==TRUE) {                     # valeurs de salinites
            Sal1 <- floor(as.numeric(tclvalue(sal1))*2)/2
            Sal2 <- ceiling(as.numeric(tclvalue(sal2))*2)/2
            sal <- c(Sal1, Sal2) }
           else{ sal <- "NULL" }
           if (any(colnames(Envir$Data) == "Depth") & any(!is.na(Envir$Data$Depth)) == TRUE) {                               # valeurs de profondeurs
            Depth1 <- floor(as.numeric(tclvalue(depth1))*2)/2
            Depth2 <- ceiling(as.numeric(tclvalue(depth2))*2)/2
            depth <- c(Depth1, Depth2) }
           else{ depth <- "NULL" }
           if (any(colnames(Envir$Data) == "Dates")) {
            start <- as.numeric(tclvalue(year1))
            end <- as.numeric(tclvalue(year2))
            if (is.numeric(Envir$Data$Dates) ==TRUE) { } else {
            months <- as.numeric(unlist(strsplit((tclvalue(mois)),"\\ "))) } }
           else { tkmessageBox(message="no date selected", icon = "warning", type = "ok", title="!Warning!")}

          cb7Value <- as.character(tclvalue(cb7Value))
          if (cb7Value=="1"){ log.trans <- "YES" }
           else { log.trans <- "NO" }
          cb2Value <- as.character(tclvalue(cb2Value))
          if (cb2Value=="1"){ na.replace <- "YES" }
           else { na.replace <- "NO" }
          cb3Value <- as.character(tclvalue(cb3Value))
          if (cb3Value=="1"){ outliers.re <- "YES" }
           else { outliers.re <- "NO" }

          rb1Value <- as.character(tclvalue(rb1Value))
          if (rb1Value=="Annual"){ time.step <- "Annual" }
          if (rb1Value=="Monthly"){ time.step <- "Monthly" }
          if (rb1Value=="Fortnight"){ time.step <- "Fortnight" }
          if (rb1Value=="Semi-fortnight"){ time.step <- "Semi-fortnight" }
          if (rb1Value=="Mono-mensual"){ time.step <- "Mono-mensual" }
          if (rb1Value=="Daily"){ time.step <- "Daily" }
          if (rb1Value=="help"){ help.timestep <- "YES"
                                 time.step <- "NULL" }
          else{ help.timestep <- "N0" }
          if (rb1Value=="auto"){ auto.timestep <- "YES"
                                 time.step <- "NULL" }
          else{ auto.timestep <- "NO" }

          rb2Value <- as.character(tclvalue(rb2Value))
          if (rb2Value=="Mean"){ aggreg <- "Mean" }
          if (rb2Value=="Median"){ aggreg <- "Median" }
          if (rb2Value=="Quantile"){ aggreg <- "Quantile" }
          if (rb2Value=="Max"){ aggreg <- "Max" }
          if (rb2Value=="help"){ help.aggreg <- "YES"
                                 aggreg <- "NULL" }
          else{ help.aggreg <- "N0" }
          if (rb2Value=="auto"){ auto.aggreg <- "YES"
                                 aggreg <- "NULL" }
          else{ auto.aggreg <- "N0" }
          FULLoption(param, depth, sal, site, rawdata="NO", select="NO", resume.reg="NO", test.normality="NO",
                 plotB="NO", selectBox, log.trans, plotZ="NO", datashow="YES",
                 help.timestep, auto.timestep, time.step, help.aggreg, auto.aggreg, aggreg,
                 mix, outliers.re, na.replace, start, end, months, norm="NO", npsu,
                 autocorr = "NO", spectrum="NO", anomaly="NO", a.barplot="NO", zsmooth="NO", local.trend="NO", test="NO")  }
     OK1.but <- tk2button(LabeledFrame7, text=" Table ",command=OnOK1, width=9)
     tkgrid(OK1.but, column=1, row=1, sticky="w")

#______________________________________________________________________________________bouton d'appel des stats desciptives sur donnees regularisees
     OnOK2 <- function()  {
            param <- Env2$variables[unique(Env2$variables.selectionnees.temp)]
            site <- Env$variables[unique(Env$variables.selectionnees.temp)]
            if (any(site == "-All-")) { site <- st  }
            if (length(param)==0 | length(site)==0 )
               { return(tkmessageBox(message="Please select a parameter and a category before to proceed"
                                    , icon = "warning", type = "ok", title="!Warning!")) }
            else{}
           mix <- "YES"
           if (any(colnames(Envir$Data) == "Salinity") & any(!is.na(Envir$Data$Salinity))==TRUE) {                     # valeurs de salinites
            Sal1 <- floor(as.numeric(tclvalue(sal1))*2)/2
            Sal2 <- ceiling(as.numeric(tclvalue(sal2))*2)/2
            sal <- c(Sal1, Sal2) }
           else{ sal <- "NULL" }
           if (any(colnames(Envir$Data) == "Depth") & any(!is.na(Envir$Data$Depth)) == TRUE) {                               # valeurs de profondeurs
            Depth1 <- floor(as.numeric(tclvalue(depth1))*2)/2
            Depth2 <- ceiling(as.numeric(tclvalue(depth2))*2)/2
            depth <- c(Depth1, Depth2) }
           else{ depth <- "NULL" }
           if (any(colnames(Envir$Data) == "Dates")) {
            start <- as.numeric(tclvalue(year1))
            end <- as.numeric(tclvalue(year2))
            if (is.numeric(Envir$Data$Dates) ==TRUE) { } else {
            months <- as.numeric(unlist(strsplit((tclvalue(mois)),"\\ "))) } }
           else { tkmessageBox(message="no date selected", icon = "warning", type = "ok", title="!Warning!")}

          cb7Value <- as.character(tclvalue(cb7Value))
          if (cb7Value=="1"){ log.trans <- "YES" }
           else { log.trans <- "NO" }
          cb2Value <- as.character(tclvalue(cb2Value))
          if (cb2Value=="1"){ na.replace <- "YES" }
           else { na.replace <- "NO" }
          cb3Value <- as.character(tclvalue(cb3Value))
          if (cb3Value=="1"){ outliers.re <- "YES" }
           else { outliers.re <- "NO" }

          rb1Value <- as.character(tclvalue(rb1Value))
          if (rb1Value=="Annual"){ time.step <- "Annual" }
          if (rb1Value=="Monthly"){ time.step <- "Monthly" }
          if (rb1Value=="Fortnight"){ time.step <- "Fortnight" }
          if (rb1Value=="Semi-fortnight"){ time.step <- "Semi-fortnight" }
          if (rb1Value=="Mono-mensual"){ time.step <- "Mono-mensual" }
          if (rb1Value=="Daily"){ time.step <- "Daily" }
          if (rb1Value=="help"){ help.timestep <- "YES"
                                 time.step <- "NULL" }
          else{ help.timestep <- "N0" }
          if (rb1Value=="auto"){ auto.timestep <- "YES"
                                 time.step <- "NULL" }
          else{ auto.timestep <- "NO" }

          rb2Value <- as.character(tclvalue(rb2Value))
          if (rb2Value=="Mean"){ aggreg <- "Mean" }
          if (rb2Value=="Median"){ aggreg <- "Median" }
          if (rb2Value=="Quantile"){ aggreg <- "Quantile" }
          if (rb2Value=="Max"){ aggreg <- "Max" }
          if (rb2Value=="help"){ help.aggreg <- "YES"
                                 aggreg <- "NULL" }
          else{ help.aggreg <- "N0" }
          if (rb2Value=="auto"){ auto.aggreg <- "YES"
                                 aggreg <- "NULL" }
          else{ auto.aggreg <- "N0" }
          FULLoption(param, depth, sal, site, rawdata="NO", select="NO", resume.reg="YES", test.normality="NO",
                 plotB="NO", selectBox, log.trans, plotZ="NO", datashow="NO",
                 help.timestep, auto.timestep, time.step, help.aggreg, auto.aggreg, aggreg,
                 mix, outliers.re, na.replace, start, end, months, norm="NO", npsu,
                 autocorr = "NO", spectrum="NO", anomaly="NO", a.barplot="NO", zsmooth="NO", local.trend="NO", test="NO")   }
     OK2.but <- tk2button(LabeledFrame7, text=" Summary ",command=OnOK2, width=9)
     tkgrid(OK2.but, column=2, row=1, sticky="w")
#_____________________________________________________________________________________fin du bouton des stats descriptives

     tkgrid(tklabel(Envir$datam, text="      "), column=0, row=5)

     HELP3.but <- tk2button(Envir$datam, image=Envir$imgHelp, text=" Help ", compound="right", command=function() { Aide3() })       # bouton d'aide 3
     tkgrid(HELP3.but, column=0, row=6, sticky="w")

#_______________________________________________________________________________________________________________________________________Onglet analyses

     LabeledFrame4 <- tkwidget(Envir$trend,"labelframe",text="Diagnostics (optional)",padx=25,pady=10, relief = "groove")
     tkgrid(LabeledFrame4, column=0, row=0, sticky="w")

     diagFrame <- tkwidget(LabeledFrame4,"labelframe",padx=0,pady=0)
     tkconfigure(diagFrame, borderwidth=0)
     tkpack(diagFrame, side="left")

     rb17 <- tk2radiobutton(diagFrame)                                                                              # radio bouton du choix du diagnostic
     rb18 <- tk2radiobutton(diagFrame)
     rb19 <- tk2radiobutton(diagFrame)
     rb20 <- tk2radiobutton(diagFrame)
     rb21 <- tk2radiobutton(diagFrame)
     rb22 <- tk2radiobutton(diagFrame)
     rb5Value <- tclVar("2")
     tkconfigure(rb17,variable=rb5Value,value="1")                                                    # valeur donnee a chaque bouton
     tkconfigure(rb18,variable=rb5Value,value="2")
     tkconfigure(rb19,variable=rb5Value,value="3")
     tkconfigure(rb20,variable=rb5Value,value="4")
     tkconfigure(rb22,variable=rb5Value,value="6")
     tkconfigure(rb21,variable=rb5Value,value="5")
     tkgrid(tklabel(diagFrame,text="Spectrum analysis*    "),rb17, sticky="w")                    # texte afficher avec chaque bouton
     tkgrid(tklabel(diagFrame,text="Autocorrelation    "),rb18, sticky="w")
     tkgrid(tklabel(diagFrame,text="Shapiro normality test    "),rb19, sticky="w")
     tkgrid(tklabel(diagFrame,text="Anomaly (color.plot)    "),rb20, sticky="w")
     tkgrid(tklabel(diagFrame,text="Anomaly (barplot)    "),rb22, sticky="w")
     tkgrid(tklabel(diagFrame,text="Seasonal decomposition*  "),rb21, sticky="w")

     imgProcess <- tclVar()
     tcl("image","create","photo",imgProcess,file=file.path(path.package("TTAinterfaceTrendAnalysis"),"aide","imgProcess.gif",fsep=.Platform$file.sep))

#_______________________________________________________________________________bouton de diagnostic
     OnOK3 <- function()  {
            param <- Env2$variables[unique(Env2$variables.selectionnees.temp)]
            site <- Env$variables[unique(Env$variables.selectionnees.temp)]
            if (any(site == "-All-")) { site <- st  }
            if (length(param)==0 | length(site)==0 )
               { return(tkmessageBox(message="Please select a parameter and a category before to proceed"
                                    , icon = "warning", type = "ok", title="!Warning!")) }
            else{}
           mix <- "YES"
           if (any(colnames(Envir$Data) == "Salinity") & any(!is.na(Envir$Data$Salinity))==TRUE) {                     # valeurs de salinites
            Sal1 <- floor(as.numeric(tclvalue(sal1))*2)/2
            Sal2 <- ceiling(as.numeric(tclvalue(sal2))*2)/2
            sal <- c(Sal1, Sal2) }
           else{ sal <- "NULL" }
           if (any(colnames(Envir$Data) == "Depth") & any(!is.na(Envir$Data$Depth)) == TRUE) {                               # valeurs de profondeurs
            Depth1 <- floor(as.numeric(tclvalue(depth1))*2)/2
            Depth2 <- ceiling(as.numeric(tclvalue(depth2))*2)/2
            depth <- c(Depth1, Depth2) }
           else{ depth <- "NULL" }
           if (any(colnames(Envir$Data) == "Dates")) {
            start <- as.numeric(tclvalue(year1))
            end <- as.numeric(tclvalue(year2))
            if (is.numeric(Envir$Data$Dates) ==TRUE) { } else {
            months <- as.numeric(unlist(strsplit((tclvalue(mois)),"\\ "))) } }
           else { tkmessageBox(message="no date selected", icon = "warning", type = "ok", title="!Warning!")}

          cb7Value <- as.character(tclvalue(cb7Value))
          if (cb7Value=="1"){ log.trans <- "YES" }
           else { log.trans <- "NO" }
          cb2Value <- as.character(tclvalue(cb2Value))
          if (cb2Value=="1"){ na.replace <- "YES" }
           else { na.replace <- "NO" }
          cb3Value <- as.character(tclvalue(cb3Value))
          if (cb3Value=="1"){ outliers.re <- "YES" }
           else { outliers.re <- "NO" }

          rb1Value <- as.character(tclvalue(rb1Value))
          if (rb1Value=="Annual"){ time.step <- "Annual" }
          if (rb1Value=="Monthly"){ time.step <- "Monthly" }
          if (rb1Value=="Fortnight"){ time.step <- "Fortnight" }
          if (rb1Value=="Semi-fortnight"){ time.step <- "Semi-fortnight" }
          if (rb1Value=="Mono-mensual"){ time.step <- "Mono-mensual" }
          if (rb1Value=="Daily"){ time.step <- "Daily" }
          if (rb1Value=="help"){ help.timestep <- "YES"
                                 time.step <- "NULL" }
          else{ help.timestep <- "N0" }
          if (rb1Value=="auto"){ auto.timestep <- "YES"
                                 time.step <- "NULL" }
          else{ auto.timestep <- "NO" }

          rb2Value <- as.character(tclvalue(rb2Value))
          if (rb2Value=="Mean"){ aggreg <- "Mean" }
          if (rb2Value=="Median"){ aggreg <- "Median" }
          if (rb2Value=="Quantile"){ aggreg <- "Quantile" }
          if (rb2Value=="Max"){ aggreg <- "Max" }
          if (rb2Value=="help"){ help.aggreg <- "YES"
                                 aggreg <- "NULL" }
          else{ help.aggreg <- "N0" }
          if (rb2Value=="auto"){ auto.aggreg <- "YES"
                                 aggreg <- "NULL" }
          else{ auto.aggreg <- "N0" }

          rb5Val <- as.character(tclvalue(rb5Value))

          if (rb5Val=="1"){ spectrum <- "YES" }
           else { spectrum <- "NO" }
          if (rb5Val=="2"){ autocorr <- "YES" }
           else { autocorr <- "NO" }
          if (rb5Val=="3"){ test.normality <- "YES" }
          else { test.normality <- "NO" }
          if (rb5Val=="4"){ anomaly <- "YES" }
          else { anomaly <- "NO" }
          if (rb5Val=="6"){ a.barplot <- "YES" }
          else { a.barplot <- "NO" }
          if (rb5Val=="5"){ zsmooth <- "YES" }
          else { zsmooth <- "NO" }

          FULLoption(param, depth, sal, site, rawdata="NO", select="NO", resume.reg="NO",test.normality,
                 plotB="NO", selectBox, log.trans, plotZ="NO", datashow="NO",
                 help.timestep, auto.timestep, time.step, help.aggreg, auto.aggreg, aggreg,
                 mix, outliers.re, na.replace, start, end, months, norm="NO", npsu,
                 autocorr, spectrum, anomaly, a.barplot, zsmooth, local.trend="NO", test="NO")  }
     OK3.but <- tk2button(LabeledFrame4, image=imgProcess, text="Run ", compound="right", command=OnOK3, width=6)
     tkpack(OK3.but, side="right")
#_______________________________________________________________________________fin du bouton de diagnostic

     tkgrid(tklabel(Envir$trend, text="* Cannot be perform with missing values", font=tkfont.create(size=7)), column=0, row=2, sticky="w")
     tkgrid(tklabel(Envir$trend, text="      "), column=0, row=3)

     LabeledFrame5 <- tkwidget(Envir$trend,"labelframe",text="Trend Analyses",padx=25,pady=10, relief = "groove")
     tkgrid(LabeledFrame5, column=0, row=4, sticky="w")

     testFrame <- tkwidget(LabeledFrame5,"labelframe",padx=0,pady=0)
     tkconfigure(testFrame, borderwidth=0)
     tkpack(testFrame, side="left")

     cb10 <- tk2checkbutton(testFrame)                                                        # check button pour cusum
     cb10Value <- tclVar("0")
     tkconfigure(cb10,variable=cb10Value)
     tkgrid(tklabel(testFrame,text="Cumulative sum*"), column=0, row=1, sticky="w")
     tkgrid(cb10, column=1, row=1)

     rb12 <- tk2radiobutton(testFrame)                                                        # radio button du choix de l'analyse
     rb13 <- tk2radiobutton(testFrame)
     rb14 <- tk2radiobutton(testFrame)
     rb15 <- tk2radiobutton(testFrame)
     rb16 <- tk2radiobutton(testFrame)
     if (is.numeric(Envir$Data$Dates) ==TRUE) { rb4Value <- tclVar("MannKen") } else {
     rb4Value <- tclVar("seasonMann") }
     tkconfigure(rb12,variable=rb4Value,value="seasonMann")                                                # valeur de chaque bouton
     tkconfigure(rb13,variable=rb4Value,value="MannKen")
     tkconfigure(rb14,variable=rb4Value,value="MixingDiagram")
     tkconfigure(rb15,variable=rb4Value,value="Extended")
     tkconfigure(rb16,variable=rb4Value,value="Lowess")
     tkgrid(tklabel(testFrame,text="Seasonal Trend "),rb12, row=2, sticky="w")               # texte affiche
     tkgrid(tklabel(testFrame,text="Global Trend "),rb13, row=3, sticky="w")
     tkgrid(tklabel(testFrame,text="Using Mixing Diagram"),rb14, row=6, sticky="w")
     tkgrid(tklabel(testFrame,text="Trend based on LOESS      "),rb16, row=4, sticky="w")

     Npsu <- tclVar(c(30))                                                                                                                   # choix de la valeur de la salinite de standardisation (30 par defaut)
     entry.npsu <-tk2entry(testFrame, width="4", textvariable=Npsu, background="white")
     tkgrid(tklabel(testFrame, text="      --> select psu"), row= 7, column=0, sticky="w")
     tkgrid(entry.npsu, row=7, column =0, sticky="e")

#_______________________________________________________________________________bouton d'appel de l'analyse temporelle
     OnOK4 <- function()  {
            param <- Env2$variables[unique(Env2$variables.selectionnees.temp)]
            site <- Env$variables[unique(Env$variables.selectionnees.temp)]
            if (any(site == "-All-")) { site <- st  }
            if (length(param)==0 | length(site)==0 )
               { return(tkmessageBox(message="Please select a parameter and a category before to proceed"
                                     , icon = "warning", type = "ok", title="!Warning!")) }
            else{}
           mix <- "YES"
           if (any(colnames(Envir$Data) == "Salinity") & any(!is.na(Envir$Data$Salinity))==TRUE) {                     # valeurs de salinites
            Sal1 <- floor(as.numeric(tclvalue(sal1))*2)/2
            Sal2 <- ceiling(as.numeric(tclvalue(sal2))*2)/2
            sal <- c(Sal1, Sal2) }
           else{ sal <- "NULL" }
           if (any(colnames(Envir$Data) == "Depth") & any(!is.na(Envir$Data$Depth)) == TRUE) {                               # valeurs de profondeurs
            Depth1 <- floor(as.numeric(tclvalue(depth1))*2)/2
            Depth2 <- ceiling(as.numeric(tclvalue(depth2))*2)/2
            depth <- c(Depth1, Depth2) }
           else{ depth <- "NULL" }
           if (any(colnames(Envir$Data) == "Dates")) {
            start <- as.numeric(tclvalue(year1))
            end <- as.numeric(tclvalue(year2))
            if (is.numeric(Envir$Data$Dates) ==TRUE) { } else {
            months <- as.numeric(unlist(strsplit((tclvalue(mois)),"\\ "))) } }
           else { tkmessageBox(message="no date selected", icon = "warning", type = "ok", title="!Warning!")}

          cb7Value <- as.character(tclvalue(cb7Value))
          if (cb7Value=="1"){ log.trans <- "YES" }
           else { log.trans <- "NO" }
          cb2Value <- as.character(tclvalue(cb2Value))
          if (cb2Value=="1"){ na.replace <- "YES" }
           else { na.replace <- "NO" }
          cb3Value <- as.character(tclvalue(cb3Value))
          if (cb3Value=="1"){ outliers.re <- "YES" }
           else { outliers.re <- "NO" }

          rb1Value <- as.character(tclvalue(rb1Value))
          if (rb1Value=="Annual"){ time.step <- "Annual" }
          if (rb1Value=="Monthly"){ time.step <- "Monthly" }
          if (rb1Value=="Fortnight"){ time.step <- "Fortnight" }
          if (rb1Value=="Semi-fortnight"){ time.step <- "Semi-fortnight" }
          if (rb1Value=="Mono-mensual"){ time.step <- "Mono-mensual" }
          if (rb1Value=="Daily"){ time.step <- "Daily" }
          if (rb1Value=="help"){ help.timestep <- "YES"
                                 time.step <- "NULL" }
          else{ help.timestep <- "N0" }
          if (rb1Value=="auto"){ auto.timestep <- "YES"
                                 time.step <- "NULL" }
          else{ auto.timestep <- "NO" }

          rb2Value <- as.character(tclvalue(rb2Value))
          if (rb2Value=="Mean"){ aggreg <- "Mean" }
          if (rb2Value=="Median"){ aggreg <- "Median" }
          if (rb2Value=="Quantile"){ aggreg <- "Quantile" }
          if (rb2Value=="Max"){ aggreg <- "Max" }
          if (rb2Value=="help"){ help.aggreg <- "YES"
                                 aggreg <- "NULL" }
          else{ help.aggreg <- "N0" }
          if (rb2Value=="auto"){ auto.aggreg <- "YES"
                                 aggreg <- "NULL" }
          else{ auto.aggreg <- "N0" }

          cb10Value <- as.character(tclvalue(cb10Value))
          if (cb10Value=="1"){ local.trend <- "YES" }
           else { local.trend <- "NO" }

          rb4Value <- as.character(tclvalue(rb4Value))                           # valeur de l'argument pris en fonction du bouton
          if (rb4Value=="seasonMann"){ test <- "SMK" }
          if (rb4Value=="MannKen"){ test <- "MK" }
          if (rb4Value=="Extended"){ test <- "ELM" }
          if (rb4Value=="Lowess"){ test <- "LOWESS" }
          if (rb4Value=="MixingDiagram"){ norm <- "YES"
                                          test <- "NO" }
          else {  norm <- "NULL" }
          npsu <- as.numeric(tclvalue(Npsu))
          FULLoption(param, depth, sal, site, rawdata="NO", select="NO", resume.reg="NO",test.normality="NO",
                 plotB="NO", selectBox, log.trans, plotZ="NO", datashow="NO",
                 help.timestep, auto.timestep, time.step, help.aggreg, auto.aggreg, aggreg,
                 mix, outliers.re, na.replace, start, end, months, norm, npsu,
                 autocorr = "NO", spectrum="NO",anomaly="NO", a.barplot="NO", zsmooth="NO", local.trend, test)   }

     OK4.but <- tk2button(LabeledFrame5, image=imgProcess, text="Run ", compound="right", command=OnOK4, width=6)
     tkpack(OK4.but, side="right")
#__________________________________________________________________________________________fin du bouton d'analyse temporelle

     tkgrid(tklabel(Envir$trend, text="* Selected periods should be longer than 1 year", font=tkfont.create(size=7)), column=0, row=5, sticky="w")
     tkgrid(tklabel(Envir$trend, text="      "), column=0, row=6)

     HELP4.but <- tk2button(Envir$trend, image=Envir$imgHelp, text=" Help ", compound="right", command=function() {  Aide4() })      # bouton d'aide 4
     tkgrid(HELP4.but, column=0, row=7, sticky="w")
}