#' @title Graphical Interface to Build an uHMM
#' @description A user-friendly interface to detect usual or extreme events in a dataset and to characterize their dynamic,
#'  by building an unsupervised Hidden Markov Model.
#' @param uHMMenv an environment in which data and results will be stored. If NULL, a local environment will be created.
#' @return Results are saved in the directory chosen by the user.
#' @references Rousseeuw, Kevin, et al. "Hybrid hidden Markov model for marine environment monitoring." Selected Topics in Applied Earth Observations and Remote Sensing, IEEE Journal of 8.1 (2015): 204-213.
#' @import tcltk tcltk2 
#' @importFrom tkrplot tkrplot
#' @export

uHMMinterface<-function(uHMMenv=NULL){
  
  if (is.null(uHMMenv)){
    uHMMenv <- new.env(hash = TRUE, size = NA)
  }
  
  # language selection window
  languageWindow<-tktoplevel()
  tktitle(languageWindow)<-"uHMM interface"
  .logoFrame(languageWindow)
  
  # Language selection frame
  languageBlock <- tkwidget(languageWindow,"labelframe",text=paste("Please select your language/\nVeuillez s",intToUtf8(0233),"lectionner votre langue",sep=""),borderwidth = 0)
  tkgrid(languageBlock,padx=c(100,100),pady=c(50,50),column=1,row=1) 
  
  englishButton<-tkbutton(languageBlock,text="English/Anglais",compound = "left", command=function(){
    tkdestroy(languageWindow)  
    .firstTab("en",uHMMenv,leftMargin=30,hscaleGraphicFrame=1.2,vscaleGraphicFrame=1.2)
  })
  tkgrid(englishButton,padx=20)
  
  
  frenchButton<-tkbutton(languageBlock,text=paste("French/Fran",intToUtf8(0231),"ais",sep=""),compound = "left",command=function(){
    tkdestroy(languageWindow) 
    .firstTab("fr",uHMMenv,leftMargin=30,hscaleGraphicFrame=1.2,vscaleGraphicFrame=1.2)
  })
  tkgrid(frenchButton,padx=20)

  
}