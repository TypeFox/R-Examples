### Choix de la periode d'echantillonnage
#' @title Select validation sample
#' @description This function is used by the MMCNS interface to let the user select the validation sample.
#' @param data dataframe from which first and last available dates should be extracted.
#' @param leftMargin the padx value of tk widgets in the first column of the tab.
#' @param tm a one row dataframe containing text to display in the interface.
#' @param win1 frame of the uHMM interface containing main tabs.
#' @param uHMMenv environment in which data and intermediate results are stored.

.periodSelectionFrame<-function(data,tm,leftMargin=30,win1,uHMMenv){

dateFrameP <- tkwidget(win1$env$prediction,"labelframe",text=tm$titlePeriodFrame2)
tkgrid(dateFrameP , row=7,padx=c(leftMargin,0),pady=c(20,20),sticky="w")


#Recuperation des dates minimum et maximum
firstObsDate<-min(.dateProcessing(data))
lastObsDate<-max(.dateProcessing(data))

minYear=min(as.numeric(substring(data[,"Dates"],1,4)))
maxYear=max(as.numeric(substring(data[,"Dates"],1,4)))

minMonth=min(as.numeric(format(chron(firstObsDate),"%m")))
maxMonth=max(as.numeric(format(chron(lastObsDate),"%m")))

minDay=min(as.numeric(format(chron(firstObsDate),"%d")))
maxDay=max(as.numeric(format(chron(lastObsDate),"%d")))

#Adaptation de l'ecriture des jours et mois s'ils sont inferieurs a 10
#ex: "1" devient "01"; "2" devient "02", etc...

if(minMonth<10){
  minMonth<-paste("0",minMonth,sep="")
}

if(maxMonth<10){
  maxMonth<-paste("0",maxMonth,sep="")
}

if(minDay<10){
  minDay<-paste("0",minDay,sep="")
}

if(maxDay<10){
  maxDay<-paste("0",maxDay,sep="")
}

#Recuperation de l'heure minimum et maximum
Hmin<-substring(chron(firstObsDate),11,18)
Hmax<-substring(chron(lastObsDate),11,18)

#Affichage de la date minimum
tkgrid(tklabel(dateFrameP,text=tm$fromLabel),row=10,column=0,padx=c(20,20))

PredFirstDate<-tkentry(dateFrameP, width=10, textvariable=tclVar(paste(minYear,"-",minMonth,"-",minDay,sep="")),background = "#ffffff");
tkgrid(PredFirstDate,row=10,column=1)

PredFirstTime<-tkentry(dateFrameP, width=7, textvariable=tclVar(Hmin),background = "#ffffff");
tkgrid(PredFirstTime,row=10,column=2)

#Affichage de la date maximum
tkgrid(tklabel(dateFrameP,text=tm$toLabel),row=11,column=0)

PredLastDate<-tkentry(dateFrameP, width=10, textvariable=tclVar(paste(maxYear,"-",maxMonth,"-",maxDay,sep="")),background = "#ffffff");
tkgrid(PredLastDate,row=11,column=1) 

PredLastTime<-tkentry(dateFrameP, width=7, textvariable=tclVar(Hmax),background = "#ffffff");
tkgrid(PredLastTime,row=11,column=2)

assign("PredFirstDate",PredFirstDate,uHMMenv)
assign("PredFirstTime",PredFirstTime,uHMMenv)
assign("PredLastDate",PredLastDate,uHMMenv)
assign("PredLastTime",PredLastTime,uHMMenv)
}

