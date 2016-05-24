#' @title Classification tab
#' @description This function generates the classification tab of the \code{\link{uHMMinterface}}, in which the user can launch the classification (only for standard users).
#' @param leftMargin left magin size of interface tabs.
#' @param tm a one row dataframe containing text to display in the interface.
#' @param console frame of the uHMM interface in which messages should be displayed. 
#' @param graphicFrame frame of the uHMM interface in which graphics should be dispayed.
#' @param win1 frame of the uHMM interface containing main tabs.
#' @param uHMMenv environment in which data and intermediate results are stored.
#' @import tcltk tcltk2
#' @importFrom stats sd
#' 

.classificationTab_standard<-function(tm,leftMargin=30,
                                     console,graphicFrame,win1,uHMMenv){
  
  # Frame informing about options by default
  infoFrame <- tkwidget(win1$env$classification,"labelframe",borderwidth=0)
  tkgrid(infoFrame , row=1,padx=c(leftMargin,0),pady=c(50,50),sticky="w")
  tkgrid(tk2label(infoFrame,text=tm$textClassificationDefault,font=tkfont.create(size=10)),padx=20,pady=10,row=1,sticky="w")
  
  lancerCalcul<-tclVar() 
  
  
  boutonRunClassifDefault<-tk2button(win1$env$classification,text=tm$runLabel,image="run",compound = "left",command=function(){  

      # message box asking whether the classification must be done      
      lancerCalcul<-tkmessageBox(title=tm$titleClassifWarning,message=tm$textClassifWarning,icon="question",type="yesno",default="yes")      
        
        if(tclvalue(lancerCalcul)=="yes"){
          
          ClassifBeginningTime<-Sys.time()
          
          # display hourglass and change the cursor in the classification frame
          classificationHourglass <- ttklabel(win1$env$classification, image="hourglass", compound="image")
          tkgrid(classificationHourglass,row=3)
          tkconfigure(win1$env$classification, cursor = "watch")
          tcl("update","idletasks")        
          
          assign("normParams",list(mean=apply(uHMMenv$trainingSet,2,mean,na.rm=TRUE),sd=apply(uHMMenv$trainingSet,2,sd,na.rm=TRUE)),envir=uHMMenv)
        
        # debut calculs classif
          if(nrow(uHMMenv$trainingSet)>20000){
            SC<-.spectralClassifGT20000(data=uHMMenv$trainingSet,normalization=TRUE,varEx=0.95,stateNb=0,uHMMinterface=TRUE,console=console,tm=tm);
          } else{
            SC<-.spectralClassifLT20000(data=uHMMenv$trainingSet,normalization=TRUE,varEx=0.95,stateNb=0,uHMMinterface=TRUE,console=console,tm=tm);
          }
          
          assign("CstateSeq",SC$CstateSeq,envir=uHMMenv)
          assign("CsymbolSeq",SC$CsymbolSeq,envir=uHMMenv)
          assign("symbCentersNorm",SC$symbCentersNorm,envir=uHMMenv)
          assign("symbCenters",SC$symbCenters,envir=uHMMenv)
          assign("nbSymbols",SC$nbSymbols,envir=uHMMenv)
          assign("gap",SC$gap,envir=uHMMenv)
          
          
          #Creation graphiques
          tkinsert(console,"1.0",tm$CstepDone)
          tkinsert(console,"1.0","\n",tm$CstepGraphics)
          tcl("update","idletasks")  
          
          print(paste("gap",uHMMenv$gap))
          .graphicsByState(data=uHMMenv$rawData[uHMMenv$trainingRows,],period=uHMMenv$trainingPeriod,stateSeq=uHMMenv$CstateSeq,step="classification",
                           directory=paste(uHMMenv$saveDirectory,tm$classificationRepertory,sep=""),uHMMinterface=TRUE,tm=tm,graphicFrame=graphicFrame)
          
          # result saving
          ClassifSavingTime<-format(Sys.time(), "%d_%m_%Y_%Hh%Mmin")
          lastClassifName<-paste("Classification_",ClassifSavingTime,sep="")
          
          save(CstateSeq,CsymbolSeq,symbCentersNorm,symbCenters,nbSymbols,gap,selectedNames,trainingRows,normParams,envir=uHMMenv,
               file=paste(uHMMenv$saveDirectory,tm$classificationRepertory,tm$rfilesRepertory,lastClassifName,sep=""))
          
          # fin calculs classif
          tkconfigure(win1$env$classification, cursor = "left_ptr")
          tkdestroy(classificationHourglass)  
          
          ClassifEndTime<-Sys.time()
          ClassifDuration<-ceiling(as.numeric(ClassifEndTime)-as.numeric(ClassifBeginningTime))
          
          #Demande a l'utilisateur si il faut afficher le rapport
          #report<-tclVar() 
          #report<-tkmessageBox(title=tm$titleModelingResultsWindow,message=tm$textClassifResultsWindow,icon="question",type="yesno",default="yes") 
          
          #if(tclvalue(report)=="yes"){
            #### A completer
          #}else{
            #### A completer
          #}
          
          #Demande a l'utilisateur si la modelisation Markovienne pour l'estimation des etats de nouvelles donnees doit etre realisee
          Markov<-tclVar() 
          Markov<-tkmessageBox(title=tm$titleEstimateMMCWindow,message=tm$textEstimateMMCWindow,icon="question",type="yesno",default="yes") 
          tk2notetab.select(win1$env$nb, tm$modelingTabLabel)  
          
          if(tclvalue(Markov)=="yes"){
            
            # opening next tab
            tk2notetab.select(win1$env$nb, tm$modelingTabLabel)
            .modelingTab(tm=tm,leftMargin=leftMargin,
                        console=console,graphicFrame=graphicFrame,win1=win1,
                        uHMMenv=uHMMenv
                       )
            
            #Display the name of the file in which results are saved
            displayLastClassif<-tklabel(win1$env$modelisation,text=lastClassifName)
            tkgrid(displayLastClassif,row=1,column=1,sticky="w")
            tkinsert(console,"1.0",paste(tm$CstepClassifResults,lastClassifName,sep=""))
            tcl("update","idletasks") 
            
          }else{
            
            tk2notetab.select(win1$env$nb, tm$modelingTabLabel)
            tkmessageBox(title=tm$titleProgramIsOverWindow,message=tm$textProgramIsOverWindow, icon="info",type="ok") 
            
          }
          
          tkinsert(console,"1.0",paste("\n---------------------------------------\n",
                                       tm$CstepDuration,ClassifDuration,tm$secondsLabel,"\n",
                                       tm$CstepVariables,length(uHMMenv$selectedNames),"\n",
                                       tm$CstepObservations,nrow(uHMMenv$trainingSet)-sum(apply(is.na(uHMMenv$trainingSet),MARGIN=1, FUN=any)) ," (",sum(apply(is.na(uHMMenv$trainingSet),MARGIN=1, FUN=any)),tm$CstepNAnumber,"\n",
                                       tm$CstepDetectedStates,length(unique(uHMMenv$CstateSeq))-1,"\n",
                                       "\n---------------------------------------\n"),sep="")
          tcl("update","idletasks")   
          
          # export en csv
          #.exportXLS(uHMMenv$rawData[uHMMenv$trainingRows,],uHMMenv$CstateSeq,paste(uHMMenv$saveDirectory,tm$classificationRepertory,tm$tablesRepertory,sep=""),tm=tm)
          varNums<-which(tolower(colnames(uHMMenv$rawData)) %in% c("dates","hours","latitude","longitude","latitudes","longitudes"))
          seq<-uHMMenv$CstateSeq-1
          seq[seq==0]<-NA
          toExport<-cbind(uHMMenv$rawData[uHMMenv$trainingRows,varNums],seq)
          write.table(toExport,file=paste(uHMMenv$saveDirectory,tm$classificationRepertory,tm$tablesRepertory,gsub("/","",tm$classificationRepertory),"_",gsub(".jpg", paste(ClassifSavingTime,".xls",sep=""), tm$seqClassifFileTitle),sep=""),row.names =FALSE)
          
        } else{ # tclvalue(lancerCalcul)=="no"
          tkmessageBox(title=tm$titleProgramIsOverWindow,message=tm$textProgramIsOverWindow, icon="info",type="ok") 
        }
      
  })
  
  tkgrid(boutonRunClassifDefault,row=2)
  
}