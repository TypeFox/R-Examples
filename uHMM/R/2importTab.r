#' @title Import tab of the uHMM interface
#' @description This function generates the import tab of the \code{\link{uHMMinterface}}, in which the user can import a data file and select a save directory.
#' @param userType character : either "default" or "expert"
#' @param leftMargin left magin size of interface tabs.
#' @param tm a one row dataframe containing text to display in the interface.
#' @param console frame of the uHMM interface in which messages should be displayed. 
#' @param graphicFrame frame of the uHMM interface in which graphics should be dispayed.
#' @param win1 frame of the uHMM interface containing main tabs.
#' @param uHMMenv environment in which data and intermediate results are stored.
#' @import tcltk tcltk2

.importTab<-function(userType,tm,leftMargin=30,
                    console,graphicFrame,win1,uHMMenv){
  
  tkfocus(win1$env$import)
  
  #Mise en place du bouton pour rechercher la base de donnees au format TXT
  browseButton<-tk2button(win1$env$import,text=tm$importTxtLabel,image="data",compound = "left",command=function(){
    imp<-.importFile(win1$env$import,dispFileRow=0,fixSumRow=5,tm=tm,
               console=console,win1=win1)

    if(exists("imp")){
      assign("rawData",imp$impData,envir=uHMMenv)
      assign("rawMoments",.dateProcessing(uHMMenv$rawData),envir=uHMMenv)
      assign("firstDataFile",imp$fileName,envir=uHMMenv)
    
      tclvalue(checkMarel)<-0
      tkconfigure(importMarel, variable = checkMarel)
    }
  })
  tkgrid(browseButton,row=0,column=0,sticky="w",padx=c(leftMargin,0),pady=c(20,0))
  
  
  ## Import MarelCarnot data
  checkMarel <- tclVar("0")
  importMarel <- tk2checkbutton(win1$env$import,variable=checkMarel , text = tm$IstepImportExample)

  set.entry.state<-function(){
    on.off<-tclvalue(checkMarel)
    print(on.off)
    if(on.off==0){
      assign("firstDataFile","MarelCarnot",uHMMenv)
      assign("rawData",MarelCarnot,envir=uHMMenv)
      assign("rawMoments",.dateProcessing(uHMMenv$rawData),envir=uHMMenv)
      
      tkgrid(tklabel(win1$env$import,text=uHMMenv$firstDataFile),row=0,column=1,sticky="w")
    }
  }  
  tkgrid(importMarel,column=0,row=1,padx=c(leftMargin*2,0),pady=c(0,10),sticky="w")
  tkbind(importMarel,"<1>",set.entry.state)
  
  
  # Directory selection
  directoryButton <- tk2button(win1$env$import,text=tm$directoryButtonLabel,image="save",compound = "left", command=function(){
    
    .saveDirectory(win1$env$import, output="saveDirectory",tm=tm,uHMMenv=uHMMenv)#,fileName=firstDataFile,fixSumRow=5)
    
    # Fix button
    fixButton<-tk2button(win1$env$import,text=tm$fixDataLabel,image="fix",compound = "left",command=function(){
      if(!exists("saveDirectory",where=uHMMenv)){
        tkmessageBox(message=tm$noDirectoryMsg,type="ok",icon="info", title=tm$warningLabel)
        
      }else if(!exists("rawData",where=uHMMenv)){
        tkmessageBox(message=tm$noFileMsg,type="ok",icon="info", title=tm$warningLabel)
        
      }else if(uHMMenv$saveDirectory==""){
        tkmessageBox(message=tm$noFileMsg,type="ok",icon="info", title=tm$warningLabel)
        
      }else{
        fixedData<-.unFixdata(impData=uHMMenv$rawData,dispTab=win1$env$import,fileName=uHMMenv$firstDataFile,tm=tm,
                               output=paste(uHMMenv$saveDirectory,tm$rawDataRepertory,tm$tablesRepertory,sep=""),console=console,win1=win1,rowNum=0)
        assign("rawData",fixedData$impData,envir=uHMMenv)
        assign("firstDataFile",fixedData$newFileName,envir=uHMMenv)
      }
    })
    tkgrid(fixButton,row=5,column=1)
        
      
  }) 
  tkgrid(directoryButton, column=0, row=2,sticky="w",padx=c(leftMargin,0))
  
  #Creation du bouton permettant de faire un summary
  summaryButton<-tk2button(win1$env$import,text=tm$summaryLabel,image="loupe",compound = "left",command=function(){
    
    if(!exists("rawData",where=uHMMenv)){
      tkmessageBox(message=tm$noFileMsg,type="ok",icon="info", title=tm$warningLabel)
      
    }else if(exists("saveDirectory",where=uHMMenv)){
        .unSummary(uHMMenv$rawData,directory=paste(uHMMenv$saveDirectory,tm$rawDataRepertory,tm$tablesRepertory,sep=""),summaryLabel=tm$summaryLabel)
      
    }else{
        .unSummary(uHMMenv$rawData,summaryLabel=tm$summaryLabel)
        
    }
  })
  tkgrid(summaryButton,row=5,column=0)
  
  
  #Creation du bloc precisant ce que l'on veut comme fichier TXT
  
  AdviceFrame <- tkwidget(win1$env$import,"labelframe",text=tm$titleImportFrame,padx=30,pady=8, relief = "groove") # cadre de texte
  tkgrid(AdviceFrame, columnspan=3, row=3, sticky="w",padx=c(leftMargin,0),pady=c(10,10))
  tkgrid(tk2label(AdviceFrame, text =tm$textImportFrame)) 
  

  #Creation du cadre de validation du fichier demande
  tkgrid(tklabel(win1$env$import, text="      "), column=1, row=6)
  
  nextTabFrame <- tkwidget(win1$env$import,"labelframe",text=tm$titleNextFrame,borderwidth = 0)
  tkgrid(nextTabFrame , column=0,row=7,padx=c(leftMargin,0),pady=c(50,0),sticky="w")
  
  
  classifButton<-tk2button(nextTabFrame,text=tm$classificationTabLabel,image="run",compound = "left",command=function(){
    
    if(!exists("rawData",where=uHMMenv)){
      tkmessageBox(message=tm$noFileMsg,type="ok",icon="info", title=tm$warningLabel)
      
    }else if(!exists("saveDirectory",where=uHMMenv)){
      tkmessageBox(message=tm$noDirectoryMsg,type="ok",icon="info", title=tm$warningLabel)
      
    }else if(uHMMenv$saveDirectory==""){
      tkmessageBox(message=tm$noDirectoryMsg,type="ok",icon="info", title=tm$warningLabel)
        
    }else{
      assign("dataCopy",uHMMenv$rawData,envir=uHMMenv)
      .variableTab(leftMargin=leftMargin,userType=userType,tm=tm,
                  console=console,graphicFrame=graphicFrame,win1=win1,
                  uHMMenv=uHMMenv)
      tk2notetab.select(win1$env$nb, tm$variableTabLabel)
      tkinsert(console,"1.0","\n---------------------------------------\n")
    }
  })
  tkgrid(classifButton,row=7,column=1)
  
  
  modelingButton<-tk2button(nextTabFrame,text=tm$modelingTabLabel,image="run",compound = "left",command=function(){
    
    if(!exists("rawData",where=uHMMenv)){
      tkmessageBox(message=tm$noFileMsg,type="ok",icon="info", title=tm$warningLabel)
      
    }else if(!exists("saveDirectory",where=uHMMenv)){
      tkmessageBox(message=tm$noDirectoryMsg,type="ok",icon="info", title=tm$warningLabel)
      
    }else if(uHMMenv$saveDirectory==""){
        tkmessageBox(message=tm$noDirectoryMsg,type="ok",icon="info", title=tm$warningLabel)
        
    }else{
      assign("dataCopy",uHMMenv$rawData,envir=uHMMenv)
      .modelingTab(leftMargin=leftMargin,tm=tm,
                  console=console,graphicFrame=graphicFrame,win1=win1,
                  uHMMenv=uHMMenv)
      tk2notetab.select(win1$env$nb, tm$modelingTabLabel)
      tkinsert(console,"1.0","\n---------------------------------------\n")
     }
  })
  tkgrid(modelingButton,row=8,column=1)
  
  
  predictButton<-tk2button(nextTabFrame,text=tm$predictionTabLabel,image="run",compound = "left",command=function(){
    
    if(!exists("rawData",where=uHMMenv)){
      tkmessageBox(message=tm$noFileMsg,type="ok",icon="info", title=tm$warningLabel)
      
    }else if(!exists("saveDirectory",where=uHMMenv)){
      tkmessageBox(message=tm$noDirectoryMsg,type="ok",icon="info", title=tm$warningLabel)
      
    }else if(uHMMenv$saveDirectory==""){
      tkmessageBox(message=tm$noDirectoryMsg,type="ok",icon="info", title=tm$warningLabel)
        
    }else{
      assign("dataCopy",uHMMenv$rawData,envir=uHMMenv)
      assign("rawValidData",uHMMenv$rawData,envir=uHMMenv)
      assign("rawValidMoments",uHMMenv$rawMoments,envir=uHMMenv)
      .predictionTab(leftMargin=leftMargin,tm=tm,
                    console=console,graphicFrame=graphicFrame,win1=win1,
                    uHMMenv=uHMMenv)
      .periodSelectionFrame(data=uHMMenv$rawValidData,tm=tm,
                           leftMargin=leftMargin,uHMMenv=uHMMenv,win1=win1)
      tk2notetab.select(win1$env$nb, tm$predictionTabLabel)

      fileLab<-tklabel(win1$env$prediction,text=uHMMenv$firstDataFile)
      tkgrid(fileLab,row=3,column=1,sticky="w")
      
      tkinsert(console,"1.0","\n---------------------------------------\n")
    }
  })
  tkgrid(predictButton,row=9,column=1)
  
  
  tkgrid(tklabel(win1$env$import, text="      "), column=1, row=8)
  

  
}
  