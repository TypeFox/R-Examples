#' @title First tab of the uHMM interface
#' @description A user-friendly interface to detect usual or extreme events in a dataset and to characterize their dynamic, by building an unsupervised Hidden Markov Model
#' @param mainWindow the mainWindow of the uHMMinterface.
#' @param tm a one row dataframe containing text to display in the interface.
#' @param language character, either "en" (for English-speaking user) of "fr" (for French-speaking user).
#' @param uHMMenv an environment in which data and results will be stored.
#' @param leftMargin left magin size of interface tabs.
#' @param hscaleGraphicFrame the hscale parameter value of the tkplot function used to create the graphic frame.
#' @param vscaleGraphicFrame the vscale parameter value of the tkplot function used to create the graphic frame.
#' @return Results are saved in the directory chosen by the user.
#' @source Rousseeuw, Kevin, et al. "Hybrid hidden Markov model for marine environment monitoring." Selected Topics in Applied Earth Observations and Remote Sensing, IEEE Journal of 8.1 (2015): 204-213.
#' @import tcltk tcltk2 
#' @importFrom tkrplot tkrplot


.firstTab<-function(language,uHMMenv,leftMargin=30,hscaleGraphicFrame=1.2,vscaleGraphicFrame=1.2){  
  tm<-.languageManagement(language)
  if (language=="fr"){fontOverviewIntro<-tkfont.create(family = "Arial", size = 12, slant = "italic")
      tcl("image","create","photo", "HMMscheme", file =file.path(path.package("uHMM"),"figures", "schema_hmm_FR.gif", fsep = .Platform$file.sep))
    }else{  fontOverviewIntro<-tkfont.create(family = "Arial", size = 16, slant = "italic")
      tcl("image","create","photo", "HMMscheme", file =file.path(path.package("uHMM"),"figures", "schema_hmm_EN.gif", fsep = .Platform$file.sep))
  }
  
  mainWindow<-tktoplevel()
  tktitle(mainWindow)<-tm$mainWindowTitle  
  
  tkwm.protocol(mainWindow, "WM_DELETE_WINDOW", function() {
    response <- tkmessageBox(title = tm$confirmLabel, icon = "question", 
                             message = tm$confirmCloseMsg, type = "yesno", 
                             parent = mainWindow)
    if (as.character(response) == "yes"){
      tkdestroy(mainWindow)
      rm(list=ls())
    }
  })
  
# Affichage logos
  .logoFrame(mainWindow)
  
# Msg Frame  
  msgFrame<-tkwidget(mainWindow,"labelframe",borderwidth = 0)
  tkgrid(msgFrame, column=2,row=1)
  #tkwm.geometry(msgFrame, "600x180")
  
  # msg frame title
  titleMsgFrameFrame<-tkwidget(msgFrame,"labelframe",borderwidth = 0) #,background = "#ffffff"
  tkgrid(titleMsgFrameFrame,sticky="w")
  tkgrid(tk2label(titleMsgFrameFrame,text="Messages :"))
  
  scrx <- tk2scrollbar(msgFrame, orientation = "horizontal", command = function(...) tkxview(console, ...))
  scry <- tk2scrollbar(msgFrame, orientation = "vertical", command = function(...) tkyview(console, ...))
  console <- tk2text(msgFrame, width = 70, height = 12, wrap = "none", xscrollcommand = function(...) tkset(scrx, ...), yscrollcommand = function(...) tkset(scry, ...))
  tkgrid(console, scry, sticky = "nsew",pady=c(0,0))
  tkgrid.rowconfigure(msgFrame, console, weight = 1)
  tkgrid.columnconfigure(msgFrame, console, weight = 1)
  tkgrid(scrx, sticky = "ew")
  
# Graphic Frame
  graphicFrame<-tkwidget(mainWindow,"labelframe",borderwidth = 0)
  tkgrid(graphicFrame, column=2,row=2,rowspan=2,sticky="w")
  GF<-tkrplot(graphicFrame,hscale=hscaleGraphicFrame,vscale=vscaleGraphicFrame,function(){par(bg = "white");plot(1,col="white",axes=FALSE,xlab=NA,ylab=NA)})
  tkgrid(GF,row=0,column=1,sticky="w")  

# Tabs
  win1<-tkwidget(mainWindow,"labelframe",borderwidth = 0)
  tkgrid(win1, column=1,row=1,rowspan=2)
  
  win1$env$nb <- tk2notebook(win1, tabs = c(tm$overviewTabLabel,tm$importTabLabel,tm$variableTabLabel,tm$classificationTabLabel,tm$modelingTabLabel,tm$predictionTabLabel))
  tkpack(win1$env$nb, fill = "both", expand = TRUE)
  
  win1$env$overview <- tk2notetab(win1$env$nb, tm$overviewTabLabel)
  win1$env$import <- tk2notetab(win1$env$nb, tm$importTabLabel)
  win1$env$variables <- tk2notetab(win1$env$nb, tm$variableTabLabel)
  win1$env$classification <- tk2notetab(win1$env$nb, tm$classificationTabLabel)
  win1$env$modelisation <- tk2notetab(win1$env$nb, tm$modelingTabLabel)
  win1$env$prediction <- tk2notetab(win1$env$nb, tm$predictionTabLabel)
  
  
  ### Overview tab
  tkgrid(tklabel(win1$env$overview, text=tm$overviewIntro, font = fontOverviewIntro),column=0,row=1, padx=c(20,20),pady = c(50,50),columnspan=10)
  
  # Scheme  
  schemeFont<-tkfont.create(slant = "italic",size=10)
  schemeFrame<-tkwidget(win1$env$overview,"labelframe",text=tm$titleScheme)
  tkgrid(schemeFrame, columnspan=1,rowspan=3,column=0, row=4,padx=leftMargin, sticky="w")
  
  # import
  cadreVisualisation<-tkwidget(schemeFrame,"labelframe")
  tkgrid(cadreVisualisation, columnspan=1, row=4,padx=c(10,10),pady=c(10,0))
  tkgrid(tklabel(cadreVisualisation, text=tm$importScheme,font = schemeFont),pady=3,column = 0, row = 4)
  tkgrid(ttklabel(schemeFrame, image="imageID", compound="image"), columnspan=1, row=5)
  
  # variable selection
  cadreVariables<-tkwidget(schemeFrame,"labelframe")
  tkgrid(cadreVariables, columnspan=1, row=6)
  tkgrid(tk2label(cadreVariables, text=tm$variableScheme,font = schemeFont),pady=3,column = 0, row = 4)
  tkgrid(ttklabel(schemeFrame, image="imageID", compound="image"), columnspan=1, row=7)
  
  # classification
  cadreClassification<-tkwidget(schemeFrame,"labelframe")
  tkgrid(cadreClassification, columnspan=1, row=8)
  tkgrid(tklabel(cadreClassification, text=tm$classificationScheme,font = schemeFont),pady=3,column = 0, row = 4)
  tkgrid(ttklabel(schemeFrame, image="imageID", compound="image"), columnspan=1, row=9)
  
  # modelisation
  cadreModeling<-tkwidget(schemeFrame,"labelframe")
  tkgrid(cadreModeling, columnspan=1, row=10)
  tkgrid(tklabel(cadreModeling, text=tm$modelingScheme,font = schemeFont),pady=3,column = 0, row = 4)
  tkgrid(ttklabel(schemeFrame, image="imageID", compound="image"), columnspan=1, row=11)
  
  # prediction
  cadrePrediction<-tkwidget(schemeFrame,"labelframe")
  tkgrid(cadrePrediction, columnspan=1, row=12,pady=c(0,10))
  tkgrid(tklabel(cadrePrediction, text=tm$predictionScheme,font = schemeFont),pady=3,column = 0, row = 4)
  
  #Recommended equipment Frame
  recoFrame<-tkwidget(win1$env$overview,"labelframe",text=tm$titleEquipment)
  tkgrid(recoFrame, columnspan=1,column=1, row=4,sticky="w")
  tkgrid(tk2label(recoFrame,text=tm$textEquipment), columnspan=1, row=4)  
  
  #Language selection block
  
  languageBlock <- tkwidget(win1$env$overview,"labelframe",text=tm$titleLanguageLabel)
  tkgrid(languageBlock,column=1,row=5,sticky="w") 
  
  languageVar2 <- tclVar(language)

  englishButton <- tkradiobutton(languageBlock) 
  frenchButton <- tkradiobutton(languageBlock)
  
  # config des boutons radio. Une seule variable tcl pour 2 boutons
  tkconfigure(englishButton,variable=languageVar2,value="en", text=tm$englishLabel)
  tkconfigure(frenchButton,variable=languageVar2,value="fr", text=tm$frenchLabel)
  tkgrid(englishButton, row=2,padx=20, sticky="w")
  tkgrid(frenchButton, row=3,padx=20, sticky="w")
  
  #User type Frame
  
  userFrame <- tkwidget(win1$env$overview,"labelframe",text=tm$titleUserLabel)
  tkgrid(userFrame,column=1,row=6,sticky="w") 
  
  user <- tclVar("default")
  
  defaultButton <- tkradiobutton(userFrame) 
  expertButton <- tkradiobutton(userFrame)
  
  # config des boutons radio. Une seule variable tcl pour 2 boutons
  tkconfigure(defaultButton,variable=user,value="default", text=tm$standardLabel)
  tkconfigure(expertButton,variable=user,value="expert", text=tm$expertLabel)
  tkgrid(defaultButton, row=2,padx=20, sticky="w")
  tkgrid(expertButton, row=3,padx=20, sticky="w")
  
  # Start button  
  StartButton<-tk2button(win1$env$overview,text=tm$startLabel,image="run",compound = "left",
                         command=function(){
                           if (language!=tclvalue(languageVar2)){
                             tm<-.languageManagement(tclvalue(languageVar2))
                           }
                           tk2notetab.select(win1$env$nb, tm$importTabLabel)
                           userType<-tclvalue(user)
                           .importTab(leftMargin=leftMargin,userType=userType,tm=tm,
                                      console=console,graphicFrame=graphicFrame,win1=win1,uHMMenv=uHMMenv
                           )
                         })
  tkgrid(StartButton,column = 1, row = 14,padx = leftMargin, pady = c(50,120),sticky="w")
  
}  
  