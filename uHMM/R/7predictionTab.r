#' @title Prediction tab
#' @description This function generates the prediction tab of the \code{\link{uHMMinterface}}, in which the user can launch the prediction of new data states.
#' @param leftMargin left magin size of interface tabs.
#' @param tm a one row dataframe containing text to display in the interface.
#' @param console frame of the uHMM interface in which messages should be displayed. 
#' @param graphicFrame frame of the uHMM interface in which graphics should be dispayed.
#' @param win1 frame of the uHMM interface containing main tabs.
#' @param uHMMenv environment in which data and intermediate results are stored.
#' @import tcltk tcltk2
#' @importFrom chron chron

.predictionTab<-function(tm,leftMargin=30,
                        console,graphicFrame,win1,uHMMenv){
  
  #Mise en place du bouton pour rechercher resultats modelisation
  importModelingButton<-tk2button(win1$env$prediction,text=tm$modelingResultsButtonLabel,image="import",compound = "left",command=function(){
    .importResults(noFileMsg=tm$noFileMsg,win1=win1,tab=win1$env$prediction,envir=uHMMenv)
    })
  tkgrid(importModelingButton,row=1,column=0,sticky="w",padx=c(leftMargin,0),pady=c(20,0))
  
  #Creation du bloc precisant ce que l'on veut comme fichier de resultats de modelisation
  AdviceFrame <- tkwidget(win1$env$prediction,"labelframe",text=tm$titleModelingImportFrame,padx=15,pady=8, relief = "groove") # cadre de texte
  tkgrid(AdviceFrame, row=2, columnspan=10,sticky="w",padx=c(leftMargin,10),pady=c(10,10))
  tkgrid(tk2label(AdviceFrame, text=tm$textModelingImportFrame), sticky="w")
  
  #Mise en place du bouton pour importer nouvelles donnees
  importDataPredButton<-tk2button(win1$env$prediction,text=tm$importPredictionTxtLabel,image="data",compound = "left",command=function(){
    
      imp<-.importFile(win1$env$prediction,dispFileRow=3,fixSumRow=5,tm=tm,console=console,win1=win1)  # import donneesBrutesAPredire
      if(exists("imp")){
        assign("rawValidData",imp$impData,envir=uHMMenv)
        assign("rawValidMoments",.dateProcessing(uHMMenv$rawValidData),envir=uHMMenv)
        assign("ValidDataFile",imp$fileName,envir=uHMMenv)
        
        .periodSelectionFrame(data=uHMMenv$rawValidData,tm=tm,
                              leftMargin=leftMargin,uHMMenv=uHMMenv,win1=win1)
        
        tclvalue(checkReimport)<-0
        tkconfigure(importFirstDatasetB, variable = checkReimport)
      }
    })
  tkgrid(importDataPredButton,row=3,column=0,sticky="w",padx=c(leftMargin,0))  
  
  
  ## Reutiliser premier fichier importe
  checkReimport <- tclVar("0")
  importFirstDatasetB  <- tkcheckbutton(win1$env$prediction,variable=checkReimport,text=tm$importFirstDataset)
  
  set.entry.state<-function(){
    on.off<-tclvalue(checkReimport)
    print(on.off)
    if(on.off==0){
      assign("rawValidData",uHMMenv$dataCopy,envir=uHMMenv)
      assign("rawValidMoments",.dateProcessing(uHMMenv$rawValidData),envir=uHMMenv)
      assign("ValidDataFile",uHMMenv$firstDataFile,envir=uHMMenv)
      .periodSelectionFrame(data=uHMMenv$rawValidData,tm=tm,leftMargin=leftMargin,uHMMenv=uHMMenv,win1=win1)
      
      tkgrid(tklabel(win1$env$prediction,text=uHMMenv$ValidDataFile),row=3,column=1,sticky="w")
      
    }
  }  
  tkgrid(importFirstDatasetB,column=0,row=4,padx=c(leftMargin*2,0),pady=c(0,10),sticky="w")
  tkbind(importFirstDatasetB,"<1>",set.entry.state)

  
  # Fix button
  fixButton<-tk2button(win1$env$prediction,text=tm$fixDataLabel,image="fix",compound = "left",command=function(){
    
    if(!exists("rawValidData",where=uHMMenv)){
      
      tkmessageBox(message=tm$noFileMsg,type="ok",icon="info", title=tm$warningLabel)
      
    }else{
      fixedData<-.unFixdata(impData=uHMMenv$rawValidData,dispTab=win1$env$prediction,fileName=uHMMenv$ValidDataFile,tm=tm,
                           output=uHMMenv$saveDirectory,console=console,win1=win1,rowNum=1)
      assign("rawValidData",fixedData$impData,envir=uHMMenv)
      assign("ValidDataFile",fixedData$newFileName,envir=uHMMenv)
      
    }
  })
  tkgrid(fixButton,row=6,column=1)

  # Summary button
  summaryButton<-tk2button(win1$env$prediction,text=tm$summaryLabel,image="loupe",compound = "left",command=function(){
  
   if(!exists("rawValidData",where=uHMMenv)){
    
     tkmessageBox(message=tm$noFileMsg,type="ok",icon="info", title=tm$warningLabel)
    
    }else{
       .unSummary(uHMMenv$rawValidData,paste(uHMMenv$saveDirectory,tm$predictionRepertory,tm$tablesRepertory,sep=""),summaryLabel=tm$summaryLabel)
   }
  })
  tkgrid(summaryButton,row=6,column=0)
  

  #Creation du bloc precisant ce que l'on veut comme fichier TXT
  AdviceFrameTXT <- tkwidget(win1$env$prediction,"labelframe",text=tm$titleImportPredictionFrame,padx=15,pady=8, relief = "groove") # cadre de texte
  tkgrid(AdviceFrameTXT, columnspan=10, row=5, sticky="w",padx=c(leftMargin,10),pady=c(10,10)) # du premier onglet
  tkgrid(tk2label(AdviceFrameTXT, text=tm$textImportPredictionFrame,font=tkfont.create(weight = "bold",size=10)), sticky="w")
  tkgrid(tk2label(AdviceFrameTXT, text=tm$textImportFrame),sticky="w")
         

  #Creation du bouton de validation du fichier demande
  runButton<-tk2button(win1$env$prediction,text=tm$RunPredictionButtonLabel,image="run",compound = "left",command=function(){
    if(!exists("rawValidData",envir = uHMMenv)){
      tkmessageBox(message=tm$noFileMsg,type="ok",icon="info", title=tm$warningLabel)
      
    }else if(!exists("estimatedHMM",envir = uHMMenv)){
      tkmessageBox(message=tm$PstepModelingWarning,
                   type="ok",icon="info", title=tm$warningLabel)
    }else{
      PredictionBeginningTime<-Sys.time()
      
      #Recuperation des dates minimum et maximum selectionnees par l'utilisation
      predDmin<-paste(substring(tclvalue(tkget(uHMMenv$PredFirstDate)),9,10),"/",substring(tclvalue(tkget(uHMMenv$PredFirstDate)),6,7),"/",substring(tclvalue(tkget(uHMMenv$PredFirstDate)),1,4),sep="")
      predDmax<-paste(substring(tclvalue(tkget(uHMMenv$PredLastDate)),9,10),"/",substring(tclvalue(tkget(uHMMenv$PredLastDate)),6,7),"/",substring(tclvalue(tkget(uHMMenv$PredLastDate)),1,4),sep="")
      
      predMinTime<-as.numeric(chron(predDmin,tclvalue(tkget(uHMMenv$PredFirstTime)),format=c("d/m/y","h:m:s")))
      predMaxTime<-as.numeric(chron(predDmax,tclvalue(tkget(uHMMenv$PredLastTime)),format=c("d/m/y","h:m:s")))
      selectedTimePred=(uHMMenv$rawValidMoments>=predMinTime & uHMMenv$rawValidMoments<=predMaxTime);
      
      #Prise en comptes de ces dates pour selectionner la portion de donnees a retenir
      assign("validationRows",selectedTimePred,envir=uHMMenv)
      assign("validationPeriod",uHMMenv$rawValidMoments[selectedTimePred],envir=uHMMenv)
      assign("validationSet",uHMMenv$rawValidData[selectedTimePred,uHMMenv$selectedNames],envir=uHMMenv) 
      
      PstateSeq<-.uHMMprediction(uHMMenv$validationSet,uHMMenv$estimatedHMM,uHMMenv$selectedNames,symbols=uHMMenv$symbCentersNorm,
                                normParams=uHMMenv$normParams,tm=tm,console=console)

      # Graphics creation
      # display in console
      tkinsert(console,"1.0",tm$PstepGraphics)
      tcl("update","idletasks")
      
      .graphicsByState(data=uHMMenv$rawValidData[uHMMenv$validationRows,],period=uHMMenv$validationPeriod,stateSeq=PstateSeq,step="prediction",
                       directory=paste(uHMMenv$saveDirectory,tm$predictionRepertory,sep=""),uHMMinterface=TRUE,tm=tm,graphicFrame=graphicFrame)
      
      # display in console
      tkinsert(console,"1.0",tm$PstepDone)

      # result saving
      assign("PstateSeq",PstateSeq,envir = uHMMenv)
      PredictionSavingTime<-format(Sys.time(), "%d_%m_%Y_%Hh%Mmin")
      lastPredictionName<-paste("MarkovPredictEstimation_",PredictionSavingTime,sep="")
      
      save(PstateSeq,selectedNames,validationPeriod,envir=uHMMenv,
           file=paste(uHMMenv$saveDirectory,tm$predictionRepertory,tm$rfilesRepertory,lastPredictionName,sep=""))

      #Display the name of the file in which results are saved
      displayLastModel<-tklabel(win1$env$prediction,text=lastPredictionName)
      tkgrid(displayLastModel,row=1,column=1,sticky="w")
      tkinsert(console,"1.0",paste(tm$PstepSavedIn,lastPredictionName,"\n",sep=""))

      PredictionEndTime<-Sys.time()
      
      PredictionDuration<-ceiling(as.numeric(PredictionEndTime)-as.numeric(PredictionBeginningTime))
      
      tkinsert(console,"1.0",paste("\n---------------------------------------\n",
                                   tm$PstepDuration,PredictionDuration,tm$secondsLabel,"\n",
                                   tm$CstepVariables,length(uHMMenv$selectedNames),"\n",
                                   tm$CstepObservations,nrow(uHMMenv$validationSet)-sum(apply(is.na(uHMMenv$validationSet),MARGIN=1, FUN=any)) ,"\n",
                                   "\n\n---------------------------------------\n",sep=""))
      
      # export en csv
      #.exportXLS(uHMMenv$rawValidData[uHMMenv$validationRows,],PstateSeq,paste(uHMMenv$saveDirectory,tm$predictionRepertory,tm$tablesRepertory,sep=""),tm=tm)
      varNums<-which(tolower(colnames(uHMMenv$rawValidData)) %in% c("dates","hours","latitude","longitude","latitudes","longitudes"))
      seq<-uHMMenv$PstateSeq-1
      seq[seq==0]<-NA
      toExport<-cbind(uHMMenv$rawValidData[uHMMenv$validationRows,varNums],seq)
      write.table(toExport,file=paste(uHMMenv$saveDirectory,tm$predictionRepertory,tm$tablesRepertory,gsub("/","",tm$predictionRepertory),"_",gsub(".jpg", paste(PredictionSavingTime,".xls",sep=""), tm$seqClassifFileTitle),sep=""),row.names =FALSE)
      
      #tkmessageBox(message=tm$textPredictionResultsWindow,type="yesno",icon="info", title=tm$titlePredictionResultsWindow)

    }
  })
  tkgrid(runButton,row=9,column=0,sticky="w",padx=c(leftMargin,0),pady=c(10,10))
}