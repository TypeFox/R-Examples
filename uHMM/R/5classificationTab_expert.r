#' @title Classification tab
#' @description This function generates the classification tab of the \code{\link{uHMMinterface}}, in which the user can select the classification method (only for expert users).
#' @param leftMargin left magin size of interface tabs.
#' @param tm a one row dataframe containing text to display in the interface.
#' @param console frame of the uHMM interface in which messages should be displayed. 
#' @param graphicFrame frame of the uHMM interface in which graphics should be dispayed.
#' @param win1 frame of the uHMM interface containing main tabs.
#' @param uHMMenv environment in which data and intermediate results are stored.
#' @import tcltk tcltk2
#' 

.classificationTab_expert<-function(tm,leftMargin=30,
                                   console,graphicFrame,win1,uHMMenv){
  
  #Normalization
  nomalizationFrame <- tkwidget(win1$env$classification,"labelframe",text=tm$normalizationFrameTitle)
  tkgrid(nomalizationFrame , row=2,padx=20,pady=10,sticky="w") 
  
  tkgrid(tklabel(nomalizationFrame,text=tm$normalizationFrameText),row=1,column=0)
  
  nomalizationTCL <- tclVar(TRUE)
  
  normButton1 <- tkradiobutton(nomalizationFrame) 
  normButton2 <- tkradiobutton(nomalizationFrame)
  
  # config des boutons radio. Une seule variable tcl pour 2 boutons
  tkconfigure(normButton1,variable=nomalizationTCL,value=TRUE, text=tm$yesLabel)
  tkconfigure(normButton2,variable=nomalizationTCL,value=FALSE, text=tm$noLabel)
  tkgrid(normButton1, row=2,padx=20, sticky="w")
  tkgrid(normButton2, row=3,padx=20, sticky="w")
  
  tkgrid(tklabel(win1$env$classification,text=" "),row=3, sticky="w")


  ########## Method selection frame
  
  methodFrame <- tkwidget(win1$env$classification,"labelframe",text=tm$selectMethodTitle)
  tkgrid(methodFrame , row=3,padx=20,pady=10) 
  
  method <-tclVar(NA)
  tkgrid(tk2label(methodFrame, text = tm$SpectralAdviceLabel,font=tkfont.create(slant = "italic",size=8)), row=3)
  
  displayParamsKM<-function(){
    tkgrid(runButton,row=11)
    
    if (exists("stopCrit",envir = uHMMenv)) tkdestroy(uHMMenv$stopCrit)
    if (exists("stateNbFrame",envir = uHMMenv)) tkdestroy(uHMMenv$stateNbFrame)
    if (exists("varExFrame",envir = uHMMenv)) tkdestroy(uHMMenv$varExFrame)
    
    #Quantification Vectorielle
    assign("varExFrame",tkwidget(methodFrame,"labelframe"),envir = uHMMenv)
    tkgrid(uHMMenv$varExFrame , row=6) 
    
    textQuantification <-tklabel(uHMMenv$varExFrame ,text=tm$vectorialQuantization)
    tkgrid(textQuantification ,row=1, columnspan=2)
    
    textVarianceExpliquee <-tklabel(uHMMenv$varExFrame ,text=tm$explainedVariance)
    tkgrid(textVarianceExpliquee ,row=2,column=0)
    assign("varExTCL", tkentry(uHMMenv$varExFrame, width=3, textvariable=tclVar("95"),background = "#ffffff"), envir=uHMMenv)
    tkgrid(uHMMenv$varExTCL, row=2, column=1, sticky="w")
    
    #Choix du nombre de groupes
    assign("stateNbFrame", tkwidget(methodFrame,"labelframe"), envir=uHMMenv)
    tkgrid(uHMMenv$stateNbFrame , row=8) 
    
    StateNbText <-tklabel(uHMMenv$stateNbFrame ,text=tm$stateNumberLabel)
    tkgrid(StateNbText ,row=1,column=0)
    assign("stateNbTCL", tkentry(uHMMenv$stateNbFrame , width=3, textvariable=tclVar("0"),background = "#ffffff"), envir=uHMMenv)
    tkgrid(uHMMenv$stateNbTCL, row=1, column=1, sticky="w")
    autoStateNbText <-tklabel(uHMMenv$stateNbFrame ,text=tm$autoStateNumberLabel)
    tkgrid(autoStateNbText, row=2, columnspan=2)
  }    

  displayParamsSP<-function(){
    tkgrid(runButton,row=11)
    
    if (exists("stopCrit",envir = uHMMenv)) tkdestroy(uHMMenv$stopCrit)
    if (exists("stateNbFrame",envir = uHMMenv)) tkdestroy(uHMMenv$stateNbFrame)
    if (exists("varExFrame",envir = uHMMenv)) tkdestroy(uHMMenv$varExFrame)
    
    # Si methode spectrale, choix du critere pour selectionner nombre de clusters (et, pour le moment, K et varExpl aussi)
    
    assign("stopCrit",tkwidget(methodFrame,"labelframe",text=tm$stopCriteriaLabel),envir = uHMMenv)
    tkgrid(uHMMenv$stopCrit , row=10)
    
    eigen <- tkradiobutton(uHMMenv$stopCrit)
    gap <- tkradiobutton(uHMMenv$stopCrit)
    
    crit <-tclVar(NA)
    
    tkconfigure(eigen, variable = crit, value = "eigen")
    tkconfigure(gap, variable = crit, value = "gap")
    
    tkgrid(tk2label(uHMMenv$stopCrit, text = tm$principalEigenValueLabel), eigen,row=1)
    tkgrid(tk2label(uHMMenv$stopCrit, text = tm$gapLabel), gap,row=2)
    
    
    #Choix du nombre de groupes
    assign("stateNbFrame",tkwidget(methodFrame,"labelframe"),envir = uHMMenv)
    tkgrid(uHMMenv$stateNbFrame , row=8) 
    
    StateNbText <-tklabel(uHMMenv$stateNbFrame ,text=tm$stateNumberLabel)
    tkgrid(StateNbText ,row=1,column=0)
    assign("stateNbTCL", tkentry(uHMMenv$stateNbFrame , width=3, textvariable=tclVar("0"),background = "#ffffff"), envir=uHMMenv)
    tkgrid(uHMMenv$stateNbTCL, row=1, column=1, sticky="w")
    autoStateNbText <-tklabel(uHMMenv$stateNbFrame ,text=tm$autoStateNumberLabel)
    tkgrid(autoStateNbText, row=2, columnspan=2)
    
    
    #Quantification Vectorielle (utilise pour spectral ? pour symboles ?)
    assign("varExFrame",tkwidget(methodFrame,"labelframe"),envir=uHMMenv)
    tkgrid(uHMMenv$varExFrame , row=9) 
    
    textQuantification <-tklabel(uHMMenv$varExFrame ,text=tm$vectorialQuantization)
    tkgrid(textQuantification ,row=1, columnspan=2)
    
    textVarianceExpliquee <-tklabel(uHMMenv$varExFrame ,text=tm$explainedVariance)
    tkgrid(textVarianceExpliquee ,row=2,column=0)
    assign("varExTCL", tkentry(uHMMenv$varExFrame, width=3, textvariable=tclVar("95"),background = "#ffffff"), envir=uHMMenv)
    tkgrid(uHMMenv$varExTCL, row=2, column=1, sticky="w")
  }
  
  
  
  # bouton Kmeans
  rb_kmeans <- tkradiobutton(methodFrame)
  tkconfigure(rb_kmeans, variable = method, value = "Kmeans", text = tm$KmeansLabel)
  tkgrid(tk2label(methodFrame),rb_kmeans, row=4,column=0,padx=20,sticky="w")
  tkbind(rb_kmeans,"<1>",displayParamsKM)
  
  # bouton spetral
  rb_spectral <- tkradiobutton(methodFrame)
  tkconfigure(rb_spectral, variable = method, value = "Spectral",text = tm$SpectralLabel)
  tkgrid(tk2label(methodFrame),rb_spectral,row=5,column=0,padx=20,sticky="w")
  tkbind(rb_spectral,"<1>",displayParamsSP)
  

  

  
  
  runButton<-tk2button(win1$env$classification,text=tm$runLabel,image="run",compound = "left",command=function(){
    
    varEx<-as.numeric(tclvalue(tkget(uHMMenv$varExTCL)));
    stateNb<-as.numeric(tclvalue(tkget(uHMMenv$stateNbTCL)))
    normalization<-as.numeric(tclvalue(nomalizationTCL))
    
      #Verification que le pourcentage n'est pas superieur a 100
      if(varEx<=100){
        if(varEx>1){

          # message box asking whether the classification must be done
          launch<-tkmessageBox(title=tm$titleClassifWarning,message=tm$textClassifWarning,icon="question",type="yesno",default="yes")
          
          if(tclvalue(launch)=="yes"){
            
            # display hourglass and change the cursor in the classification frame
            classificationHourglass <- ttklabel(win1$env$classification, image="hourglass", compound="image")
            tkgrid(classificationHourglass,row=3,rowspan=4, column=2)
            tkconfigure(win1$env$classification, cursor = "watch")
            
            ClassifBeginningTime<-Sys.time()
            
            # je sauvegarde moyenne et ecart-type pour projeter l'echantillon de validation dans le bon espace
            if(normalization){
              assign("normParams",list(mean=apply(uHMMenv$trainingSet,2,mean,na.rm=TRUE),sd=apply(uHMMenv$trainingSet,2,sd,na.rm=TRUE)),envir=uHMMenv)
            }else{
              assign("normParams",list(mean=rep(0,ncol(uHMMenv$trainingSet)),sd=rep(1,ncol(uHMMenv$trainingSet))),envir=uHMMenv)
            }
         
            # calculs classif
            if (tclvalue(method)=="Kmeans"){   
              
                SC<-.classifKmeans(data=uHMMenv$trainingSet,
                                  varEx=varEx*0.01,stateNb=stateNb,normalization=normalization,uHMMinterface=TRUE,console=console,tm=tm)
                
              }else{ # si classification spectrale
                if(nrow(uHMMenv$trainingSet)>20000){
                  SC<-.spectralClassifGT20000(data=uHMMenv$trainingSet,
                                             varEx=varEx*0.01,stateNb=stateNb,normalization=normalization,uHMMinterface=TRUE,console=console,tm=tm)
                } else{
                  SC<-.spectralClassifLT20000(data=uHMMenv$trainingSet,
                                             varEx=varEx*0.01,stateNb=stateNb,normalization=normalization,uHMMinterface=TRUE,console=console,tm=tm)
                }
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
            
            tkconfigure(win1$env$classification, cursor = "left_ptr")
            tkdestroy(classificationHourglass)  
            
            ClassifEndTime<-Sys.time()
            
            ClassifDuration<-ceiling(as.numeric(ClassifEndTime)-as.numeric(ClassifBeginningTime))

          #Demande a l'utilisateur si il faut afficher le rapport
          #  report<-tclVar() 
          #  report<-tkmessageBox(title=tm$titleModelingResultsWindow,message=tm$textClassifResultsWindow,icon="question",type="yesno",default="yes") 
            
          #  if(tclvalue(report)=="yes"){
              ### A COMPLETER
          #  }else{
              ### A COMPLETER
          #  }
            
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
            
            
            
          } else{
            tkmessageBox(title=tm$titleProgramIsOverWindow,message=tm$textProgramIsOverWindow, icon="info",type="ok") 
          }
          
          
          
        } else{
          tkmessageBox(title=tm$titleVarianceErrorWindow,message=tm$textVarianceErrorWindow, icon="info",type="ok")        
        }
      } else{      
        tkmessageBox(title=tm$titleVarianceErrorWindow,message=tm$textVarianceErrorWindow, icon="info",type="ok")      
      }
      
    
  })
  
  
  
  
}