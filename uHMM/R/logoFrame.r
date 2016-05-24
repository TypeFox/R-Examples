#' @title Display logos in uHMM interface
#' @param window the window or the frame in which logos should be displayed.

.logoFrame<-function(window){
  
  tcl("image","create","photo", "imageID" , file =file.path(path.package("uHMM"),"figures", "flecheBas.gif", fsep = .Platform$file.sep))
  tcl("image","create","photo", "run" , file =file.path(path.package("uHMM"),"figures", "run.gif", fsep = .Platform$file.sep))
  tcl("image","create","photo", "data" , file =file.path(path.package("uHMM"),"figures", "data.gif", fsep = .Platform$file.sep))
  tcl("image","create","photo", "save" , file =file.path(path.package("uHMM"),"figures", "save.gif", fsep = .Platform$file.sep))
  tcl("image","create","photo", "loupe" , file =file.path(path.package("uHMM"),"figures", "loupe.gif", fsep = .Platform$file.sep))
  tcl("image","create","photo", "fix" , file =file.path(path.package("uHMM"),"figures", "fix.gif", fsep = .Platform$file.sep))     
  tcl("image","create","photo", "import" , file =file.path(path.package("uHMM"),"figures", "import.gif", fsep = .Platform$file.sep))  
  tcl("image", "create", "photo", "UlcoLogo" , file =file.path(path.package("uHMM"),"figures", "UlcoLogo.gif", fsep = .Platform$file.sep))
  tcl("image","create","photo", "LisicLogo", file= file.path(path.package("uHMM"),"figures", "LisicLogo.gif", fsep = .Platform$file.sep))
  tcl("image","create","photo", "IfremerLogo" , file= file.path(path.package("uHMM"),"figures", "IfremerLogo.gif", fsep = .Platform$file.sep))
  tcl("image","create","photo", "AeapLogo" , file= file.path(path.package("uHMM"),"figures", "AeapLogo.gif", fsep = .Platform$file.sep))
  tcl("image","create","photo", "hourglass" , file =file.path(path.package("uHMM"),"figures", "hourglass.gif", fsep = .Platform$file.sep))
  
  # Logo Frame  
  logoFrame<-tkwidget(window,"labelframe",borderwidth = 0) #,background = "#ffffff"
  tkgrid(logoFrame, column=1,row=3,pady=c(0,10)) #,sticky = "s"

  padxLogo<-10
  tkgrid(ttklabel(logoFrame, image="UlcoLogo", compound="image"),row=0,column=1,padx=padxLogo)
  tkgrid(ttklabel(logoFrame, image="LisicLogo", compound="image"),row=0,column=2,padx=padxLogo)
  tkgrid(ttklabel(logoFrame, image="IfremerLogo", compound="image"),row=0,column=3,padx=padxLogo)
  tkgrid(ttklabel(logoFrame, image="AeapLogo", compound="image"),row=0,column=4,padx=padxLogo)

}