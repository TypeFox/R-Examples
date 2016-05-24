#' @title representationPlot function
#' @description This function is used by the MMCNS interface to display plots of selected variables
#' @param tm a one row dataframe containing text to display in the interface.
#' @param data dataframe containing data that should be plotted.
#' @param varNames caracter vector of selected variable names.
#' @param period numeric vector of selected observation times.
#' @import tcltk tcltk2
#' @importFrom graphics plot
#' @importFrom tkrplot tkrplot
#' @importFrom chron chron


.representationPlot<-function(data,period,tm,varNames){
	
	plotW<-tktoplevel() 
	tktitle(plotW)<-tm$titlePlotWindow
	
	selectionFrame<-tkwidget(plotW,"labelframe",borderwidth = 0)
	tkgrid(selectionFrame,column=0, row=0)
	
	varList<-tk2listbox(selectionFrame, height=20, selectmode="single", background="white")
	
	
	for(i in 1:length(varNames)){

		tkinsert(varList,"end",varNames[i])
	
	}

	tkselection.set(varList, 0) 
	tkgrid(varList, row=0,column=0)
	
	displayButton<-tk2button(selectionFrame,text=tm$displayButtonLabel,command=function(){

		if (tclvalue(tkcurselection(varList))!="") {
                    
			#Recherche du parametre selectionne
			selection<-as.numeric(tkcurselection(varList))+1
					
			tempSelect<-varNames[selection]
			
			myPlot<-tkrplot(plotW,hscale=2.5,vscale=2,function()plot(data[,tempSelect]~chron(period),ylab="Level",xlab="Date",main=tempSelect))
			tkgrid(myPlot,row=0,column=1,sticky="w")
			tkconfigure(myPlot, bg="white")
			
					
        } else { tkmessageBox(message=tm$noParamMsg,type="ok",icon="info", title=tm$warningLabel) } #On avertit l'utilisateur qu'aucun parametre n'est selectionne

	})
	tkgrid(displayButton,row=1,column=0)
	
}







#' @title savePlot function
#' @keywords internal 
#' @description This function is used by the MMCNS interface to save plots of selected variables
#' @param data dataframe containing data that should be plotted.
#' @param period numeric vector of selected observation times.
#' @param output directory path in which graphics should be saved.
#' @param tm a one row dataframe containing text to display in the interface.
#' @import tcltk tcltk2
#' @importFrom grDevices jpeg dev.off
#' @importFrom graphics plot


.savePlot<-function(data,period,tm,output){
  
  for(i in 1:ncol(data)){
    
    figure<-paste(output,colnames(data)[i],".jpg",sep="")
    jpeg(figure,750,600)
    plot(data[,i]~chron(period),xlab="Dates",ylab="Level",main=colnames(data)[i])
    dev.off()
    
    figure<-paste(output,"boxplot_",colnames(data)[i],".jpg",sep="")
    jpeg(figure,750,600)
    boxplot(data[,i],main=colnames(data)[i])
    dev.off()
  }
  
}



# boxplotWindow function
#' @title Display a boxplot from MMCNS interface
#' @description Function used by the MMCNS interface to create a window in which the user can display boxplots of selected variables
#' @param tm a one row dataframe containing text to display in the interface.
#' @param data dataframe containing data that should be plotted.
#' @param varNames caracter vector of selected variable names.
#' @import tcltk tcltk2
#' @importFrom tkrplot tkrplot
#' @importFrom graphics boxplot

.boxplotWindow<-function(data,varNames,tm){
  
  boxplotWin<-tktoplevel() 
  tktitle(boxplotWin)<-tm$titleBoxplotWindow
  
  selectionFrame<-tkwidget(boxplotWin,"labelframe",borderwidth = 0)
  tkgrid(selectionFrame,column=0, row=0)
  
  varList<-tk2listbox(selectionFrame, height=20, selectmode="single", background="white")
  
  
  for(i in 1:length(varNames)){tkinsert(varList,"end",varNames[i])}
  
  tkselection.set(varList, 0) 
  tkgrid(varList, row=0,column=0)
  
  displayButton<-tk2button(selectionFrame,text=tm$displayButtonLabel,command=function(){
    
    if (tclvalue(tkcurselection(varList))!="") {
      
      #Recherche du parametre selectionne
      selection<-as.numeric(tkcurselection(varList))+1
      
      tempVarNames<-varNames[selection]
      
      myBoxplot<-tkrplot(boxplotWin,hscale=2.5,vscale=2,function(){boxplot(data[,tempVarNames],main=tempVarNames)})
      tkgrid(myBoxplot,row=0,column=1)
      tkconfigure(myBoxplot, bg="white")
      
      
    } else { tkmessageBox(message=tm$noParamMsg,type="ok",icon="info", title=tm$warningLabel) } #On avertit l'utilisateur qu'aucun parametre n'est selectionne
    
  })
  tkgrid(displayButton,row=1,column=0)
  
}








#' @title PCAwindow function
#' @description This function is used by the MMCNS interface to display the correlation circles from the PCA performed on imported data.
#' @param tm a one row dataframe containing text to display in the interface.
#' @param quanti.sup a vector indicating the indexes of the quantitative supplementary variables.
#' @param data dataframe containing data that should be plotted.
#' @import tcltk tcltk2
#' @importFrom tkrplot tkrplot
#' @importFrom FactoMineR PCA
#' @importFrom graphics plot

.PCAwindow<-function(data,quanti.sup,tm){
  
  fenetrePlotACP<-tktoplevel() 
  tktitle(fenetrePlotACP)<-tm$titlePcaCircle
  
  ToRemove<-apply(is.na(data),MARGIN=1, FUN=any); #margin=1 travail sur les lignes any si au moins 1 TRUE
  data <- data[!ToRemove,]
  
  if (is.null(quanti.sup)){
    PCAdudi<-PCA(X=data,scale.unit = TRUE,ncp=ncol(data),graph=FALSE)
  }else{
    PCAdudi<-PCA(X=data,quanti.sup=quanti.sup,scale.unit = TRUE,ncp=ncol(data),graph=FALSE)
  }
  
  monACP<-tkrplot(fenetrePlotACP,hscale=2.5,vscale=2,function()plot(PCAdudi, axes = c(1, 2), choix = "var",cex.axis=1.5,cex.lab=1.5,cex=1.5))
  tkgrid(monACP,row=0,sticky="w")
  tkconfigure(monACP, bg="white")
  
  
  #Cadre pour la selection des dimensions
  textblanc<-tklabel(fenetrePlotACP,text=" ")
  tkgrid(textblanc,row=1, sticky="w")
  
  dimensionFrame <- tkwidget(fenetrePlotACP,"labelframe",borderwidth = 0)
  tkgrid(dimensionFrame , row=2) 
  
  textDim1 <-tklabel(dimensionFrame ,text="x : Dim ")
  tkgrid(textDim1 ,row=0, column=0)
  
  laDim1<-tkentry(dimensionFrame, width=3, textvariable=tclVar("1")) 
  tkgrid(laDim1, row=0, column=1, sticky="w")
  
  textDim1 <-tklabel(dimensionFrame ,text="      ")
  tkgrid(textDim1 ,row=0, column=2)
  
  textDim2 <-tklabel(dimensionFrame ,text="y : Dim ")
  tkgrid(textDim2 ,row=0, column=3)
  
  laDim2<-tkentry(dimensionFrame, width=3, textvariable=tclVar("2")) 
  tkgrid(laDim2, row=0, column=4, sticky="w")
  
  textDim1 <-tklabel(dimensionFrame ,text="      ")
  tkgrid(textDim1 ,row=0, column=5)
  
  
  
  boutonVisualiser<-tk2button(dimensionFrame,text=tm$displayButtonLabel,command=function(){
    
    #Recherche du parametre selectionne
    Dim1Ecrit<-as.numeric(tclvalue(tkget(laDim1)));
    Dim2Ecrit<-as.numeric(tclvalue(tkget(laDim2)));
    
    if(Dim1Ecrit>0){
      if(Dim2Ecrit>0){
        if(Dim1Ecrit<(ncol(data)+1)){
          if(Dim2Ecrit<(ncol(data)+1)){
            
            #Affichage de l'ACP avec les dimensions selectionnees
            monACP<-tkrplot(fenetrePlotACP,hscale=2.5,vscale=2,function()plot(PCAdudi, axes = c(Dim1Ecrit, Dim2Ecrit), choix = "var",cex.axis=1.5,cex.lab=1.5,cex=1.5))
            tkgrid(monACP,row=0,sticky="w")
            tkconfigure(monACP, bg="white")
            
          } else{
            tkmessageBox(message=paste(tm$dimensionWarning,ncol(data),"!",sep=""),type="ok",icon="info", title=tm$warningLabel) #On avertit l'utilisateur que les dimensions sont incorrects
          }
        } else{
          tkmessageBox(message=paste(tm$dimensionWarning,ncol(data),"!",sep=""),type="ok",icon="info", title=tm$warningLabel) #On avertit l'utilisateur que les dimensions sont incorrects
        }
      } else{
        tkmessageBox(message=paste(tm$dimensionWarning,ncol(data),"!",sep=""),type="ok",icon="info", title=tm$warningLabel)#On avertit l'utilisateur que les dimensions sont incorrects
      }
    } else{
      tkmessageBox(message=paste(tm$dimensionWarning,ncol(data),"!",sep=""),type="ok",icon="info", title=tm$warningLabel) #On avertit l'utilisateur que les dimensions sont incorrects
    }
    
  })
  tkgrid(boutonVisualiser,row=0,column=6)
  
  return(PCAdudi)
}

