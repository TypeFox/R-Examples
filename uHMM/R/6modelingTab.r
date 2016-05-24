#' logo window
#' @title Modeling tab
#' @description This function generates the modeling tab of the \code{\link{uHMMinterface}}, in which the user can launch the estimation of model parameters.
#' @param leftMargin left magin size of interface tabs.
#' @param tm a one row dataframe containing text to display in the interface.
#' @param console frame of the uHMM interface in which messages should be displayed. 
#' @param graphicFrame frame of the uHMM interface in which graphics should be dispayed.
#' @param win1 frame of the uHMM interface containing main tabs.
#' @param uHMMenv environment in which data and intermediate results are stored.
#' @import tcltk2

.modelingTab<-function(tm,leftMargin=30,
                      console,graphicFrame,win1,uHMMenv){
  
  #Mise en place du bouton pour rechercher resultats classification
  importClassifButton<-tk2button(win1$env$modelisation,text=tm$classificationResultsButtonLabel,image="import",compound = "left",command=function(){
    .importResults(noFileMsg=tm$noFileMsg,win1=win1,tab=win1$env$modelisation,envir=uHMMenv)
    trainingSet<-uHMMenv$rawData[uHMMenv$trainingRows,uHMMenv$selectedNames]
    assign("trainingPeriod",uHMMenv$rawMoments[uHMMenv$trainingRows],envir=uHMMenv)
    assign("trainingSet",trainingSet,envir=uHMMenv)
  })
  tkgrid(importClassifButton,row=1,column=0,sticky="w",padx=c(leftMargin,0),pady=c(20,0))
  
  
  
  #Creation du bloc precisant ce que l'on veut comme fichier TXT
  
  AdviceFrame <- tkwidget(win1$env$modelisation,"labelframe",text=tm$titleClassifImportFrame,padx=30,pady=8, relief = "groove") # cadre de texte
  tkgrid(AdviceFrame, columnspan=2, row=3, sticky="w",padx=c(leftMargin,0),pady=c(10,10))
  tkgrid(tk2label(AdviceFrame, text=tm$textClassifImportFrame), sticky="w")
  
  
  
  #Creation du bouton de validation du fichier demande
  
  runButton<-tk2button(win1$env$modelisation,text=tm$runLabel,image="run",compound = "left",command=function(){

    if(!exists("nbSymbols",where=uHMMenv)){
      
      tkmessageBox(message=tm$MstepClassificationWarning,type="ok",icon="info", title=tm$warningLabel)
      
    }else{
      
      ModelingBeginningTime<-Sys.time()
      
      
      # Estimation des parametres (matrices de transitions, emissions, et vecteur de probas initiales)
      #Display in console
      tkinsert(console,"1.0",paste(tm$MstepParameterEstimation,"\n",sep=""))
      tcl("update","idletasks")
      
      params<-HMMparams(stateSeq=uHMMenv$CstateSeq,symbolSeq=uHMMenv$CsymbolSeq)


      #Display in console
      tkinsert(console,"1.0",paste(tm$MstepAssignSymbols,"\n",sep=""))
      tcl("update","idletasks")
      
      # pour le moment les symboles sont affectes en fonction des distances dans l'espace des donnees normalisees
      MsymbolSeq<-.kppvObsCentre(data=uHMMenv$trainingSet,gravityCenters=uHMMenv$symbCentersNorm,normParams=uHMMenv$normParams)
      #MsymbolSeq<-kppvObsCentre(data=uHMMenv$trainingSet,gravityCenters=uHMMenv$symbCenters,normalization=FALSE)
      
      #Estimation des etats des nouvelles donnees
      #Display in console
      tkinsert(console,"1.0",paste(tm$MstepViterbi,"\n",sep=""))
      tcl("update","idletasks")
      print("OK")
      MVE<-.MarkovViterbiEstimation(symbolSeq=MsymbolSeq,nbStates=uHMMenv$gap, nbSymbols=uHMMenv$nbSymbols, transProbs=params$trans, startProb=params$startProb, emissionProbs=params$emis)
      print("OK")
      assign("estimatedHMM",MVE$estimatedHMM,envir=uHMMenv)
      assign("MstateSeq",MVE$stateSeq,envir=uHMMenv)
      
      #Creation des figures pour l'interpretation des donnees
      #Display in console
      tkinsert(console,"1.0",tm$MstepDone)
      tkinsert(console,"1.0",tm$CstepGraphics)
      tcl("update","idletasks")  
      
      print(length(uHMMenv$trainingPeriod))
      .graphicsByState(data=uHMMenv$rawData[uHMMenv$trainingRows,],period=uHMMenv$trainingPeriod,stateSeq=uHMMenv$MstateSeq,step="modeling",
                       directory=paste(uHMMenv$saveDirectory,tm$modelingRepertory,sep=""),uHMMinterface=TRUE,tm=tm,graphicFrame=graphicFrame)

      
      # result saving
      ModelingSavingTime<-format(Sys.time(), "%d_%m_%Y_%Hh%Mmin")
      lastModelName<-paste("MarkovModelEstimation_",ModelingSavingTime,sep="")
      
      save(MstateSeq,selectedNames,estimatedHMM,symbCentersNorm,symbCenters,gap,normParams,envir=uHMMenv,
           file=paste(uHMMenv$saveDirectory,tm$modelingRepertory,tm$rfilesRepertory,lastModelName,sep=""))
      
      
      #Display the name of the file in which results are saved
      displayLastModel<-tklabel(win1$env$prediction,text=lastModelName)
      tkgrid(displayLastModel,row=1,column=1,sticky="w")
      tkinsert(console,"1.0",paste(tm$MstepEstimateModel,lastModelName,"\n",sep=""))
      
      ModelingEndTime<-Sys.time()
      
      ModelingDuration<-ceiling(as.numeric(ModelingEndTime)-as.numeric(ModelingBeginningTime))
      
      tkinsert(console,"1.0",paste("\n---------------------------------------\n",
                                   tm$MstepDuration,ModelingDuration,tm$secondsLabel,"\n",
                                   tm$CstepVariables,length(uHMMenv$selectedNames),"\n",
                                   tm$CstepObservations,nrow(uHMMenv$trainingSet)-sum(apply(is.na(uHMMenv$trainingSet),MARGIN=1, FUN=any)) ,"\n",
                                   tm$MstepStates,length(unique(uHMMenv$CstateSeq))-1,"\n",
                                   tm$MstepSymbols,nrow(uHMMenv$symbCenters)-1,"\n",
                                   "\n\n---------------------------------------\n",sep=""))
      
      # export en csv
      #.exportXLS(uHMMenv$rawData[uHMMenv$trainingRows,],uHMMenv$MstateSeq,paste(uHMMenv$saveDirectory,tm$modelingRepertory,tm$tablesRepertory,sep=""),tm=tm)
      varNums<-which(tolower(colnames(uHMMenv$rawData)) %in% c("dates","hours","latitude","longitude","latitudes","longitudes"))
      seq<-uHMMenv$MstateSeq-1
      seq[seq==0]<-NA
      toExport<-cbind(uHMMenv$rawData[uHMMenv$trainingRows,varNums],seq)
      write.table(toExport,file=paste(uHMMenv$saveDirectory,tm$modelingRepertory,tm$tablesRepertory,gsub("/","",tm$modelingRepertory),"_",gsub(".jpg", paste(ModelingSavingTime,".xls",sep=""), tm$seqClassifFileTitle),sep=""),row.names =FALSE)
      
      #tkmessageBox(message=tm$textModelingResultsWindow,type="yesno",icon="info", title=tm$titleModelingResultsWindow)
      
      .predictionTab(leftMargin=leftMargin,tm=tm,
                    console=console,graphicFrame=graphicFrame,win1=win1,
                    uHMMenv=uHMMenv
                    )
      tk2notetab.select(win1$env$nb, tm$predictionTabLabel)
      
      
    }
    
  })
  tkgrid(runButton,row=7,column=0,padx=c(leftMargin,0),pady=c(10,10),sticky="w")
  
  
  HMM_scheme_frame<-tkwidget(win1$env$modelisation ,"labelframe",borderwidth = 0)
  tkgrid(HMM_scheme_frame,row=8,columnspan=10)
  HMM_scheme<- ttklabel(HMM_scheme_frame, image="HMMscheme", compound="image")
  
  tkgrid(HMM_scheme,row=0,column=0)
  
}