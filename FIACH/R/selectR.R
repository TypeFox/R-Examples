selectR<-function(update=TRUE){
  tt<-tktoplevel()
 
  tktitle(tt)<-"File Selector"
  options(show.error.messages = FALSE) 
  #########################################
  ####### FUNCTIONS AND VARIABLES #########
  #########################################
  changeDir<-function(){
    CurserInd<-as.integer(tkcurselection(tl))+1 ## get indices of the cursor 
    selection<-init[CurserInd]                     ## new direcory is the index as specified by cursor 
    
        
    if(selection==".."){newDir<-dirname(current)}else{newDir<-selection}    ## if going up a directory get previous directory name       
    
    tclvalue(tlv)<<-""                             ## delete everything in directory frame
    if(wind){
      if(selection == ".." & any(dirname(current)==drives)){init<<-drives}else{  ## if on windows and want to enter the "drive" folder list drives
        newFolds<-list.files(newDir,full.names = TRUE)                                  ## list the files in the new directory
        log<-dir.exists(newFolds)
        init<<-c("..",newFolds[log])                       ## show directories only
      }}else{
        newFolds<-list.files(newDir,full.names = TRUE)                                  ## list the files in the new directory
        log<-dir.exists(newFolds)
        init<<-c("..",newFolds[log])  
      }
    
    for (i in 1:length(init)){                                                        ## loop over folders and add them to listbox1
      tkinsert(tl,"end",basename(init[i]))
    }
    current<<-newDir
    tclvalue(tl2v)<-""
    if(exists("newFolds")){
      init2<<-newFolds[!log]}else{
        init2<<-""
      }
    
    if(length(init2)==0){tkinsert(tl2,"end","")}else{
      for (i in 1:length(init2)){
        tkinsert(tl2,"end",basename(init2[i]))
      }
    }
  }
  fileFilter<-function(){
    l2<-length(init2)-1
    for(i in 0:l2){
      tkdelete(tl2,0)  
    }
    
    init2<<-list.files(current,full.names = TRUE,include.dirs = FALSE,pattern=tclvalue(expr))
    if(length(init2)==0){tkinsert(tl2,"end","")}else{
      for (i in 1:length(init2)){
        tkinsert(tl2,"end",basename(init2[i]))
      }
    }
  }
  passFiles<-function(){
    ind2<-as.integer(tkcurselection(tl2))+1
    ind2d<-as.integer(tkcurselection(tl2))
    
    if(length(ind2)>0){
      for (i in 1:length(ind2)){
        tkinsert(tl3,"end",init2[ind2[i]])
      }
    }
    
    tclvalue(num)<<-paste(length(tclFile2R(tl3v)),"Files Selected")
  }
  deleteFiles<-function(){
    l<-length(as.integer(tkcurselection(tl3)))
    while(l>0){
      ind2<-as.integer(tkcurselection(tl3))
      tkdelete(tl3,ind2[1])
      l<-length(as.integer(tkcurselection(tl3)))
    }
    
    tclvalue(num)<<-paste(length(tclFile2R(tl3v)),"Files Selected")
  }
  doneFunc<-function(){
    tclvalue(ret)<<-tclvalue(tl3v)
    
    
  }
  tclFile2R<-function(tclvar){
    test<-tclvalue(tclvar)
    
    if(grepl("([{]|[}])",test)){
      spFiles<-strsplit(test,"([{]|[}])")[[1]]
      nSpaInd<-unlist(lapply(spFiles,function(x){grepl("^ ",x)|grepl(" $",x)}))
      nspa<-spFiles[nSpaInd]    ## no spaces
      spa<-spFiles[!nSpaInd]    ## spaces
      if(sum(nSpaInd)==0){spl<-""}else{
        spl<-strsplit(nspa," ")[[1]]
      }
      tot<-c(spa[nchar(spa)>0],spl[nchar(spl)>0])
    }else{
      tot<-strsplit(test," ")[[1]]
    }
    return(tot)
  }
  helpBox<-function(){
    txt<-paste("Top Left Window = left click to navigate directories",
               "Top Right Window = highlight and right click to select files",
               "Bottom Window = highlight and right click to remove files",sep="\n")
    tkmessageBox(title = "Help",message = txt, icon = "info", type = "ok")
  }
  
  
  expr<-tclVar("^.*nii$")
  tlv<-tclVar()
  tl2v<-tclVar()
  tl3v<-tclVar()
  ret<-tclVar(0)
  num<-tclVar("0 Files Selected")
  
  ##########################################
  ######### SCROLLBARS #####################
  ##########################################
  scry <- tkscrollbar(tt, repeatinterval=5, command=function(...)tkyview(tl,...))
  scrx <- tkscrollbar(tt, repeatinterval=5,orient="horizontal", command=function(...)tkxview(tl,...))
  
  scry2 <- tkscrollbar(tt, repeatinterval=5, command=function(...)tkyview(tl2,...))
  scrx2 <- tkscrollbar(tt, repeatinterval=5,orient="horizontal", command=function(...)tkxview(tl2,...))
  
  scry3 <- tkscrollbar(tt, repeatinterval=5, command=function(...)tkyview(tl3,...))
  scrx3 <- tkscrollbar(tt, repeatinterval=5,orient="horizontal", command=function(...)tkxview(tl3,...))
  
  ###########################################
  ######### DEFINE WIDGETS ##################
  ###########################################
  tl<-tklistbox(tt,height=10,width=40,selectmode="single",background="white",yscrollcommand=function(...)tkset(scry,...),xscrollcommand = function(...)tkset(scrx,...),listvariable=tlv)
  tl2<-tklistbox(tt,height=10,width=40,selectmode="extended",background="white",yscrollcommand=function(...)tkset(scry2,...),xscrollcommand = function(...)tkset(scrx2,...),listvariable=tl2v)
  tl3<-tklistbox(tt,height=10,width=80,selectmode="extended",background="white",listvariable=tl3v,yscrollcommand=function(...)tkset(scry3,...),xscrollcommand = function(...)tkset(scrx3,...))
  reg<-tkentry(parent = tt,textvariable = expr)
  done<-tkbutton(parent = tt,text="Done",command = doneFunc)
  help<-tkbutton(parent = tt,text="Help",command = helpBox)
  filtLab<-ttklabel(tt,text="Filter")
  numFiles<-ttklabel(tt,textvariable = num)
  ###########################################
  ######## LAYOUT WIDGETS ###################
  ###########################################
  tkgrid(tl,scry,tl2,scry2)
  tkgrid(scrx,tklabel(tt,text="    "),scrx2)
  tkgrid(done,filtLab,reg)
  
  tkgrid.configure(filtLab,sticky="e")
  tkgrid.configure(reg,sticky="w")
  tkgrid(tl3,scry3,columnspan=3)
  tkgrid(scrx3,columnspan=3)
  tkgrid(numFiles,help)
  tkgrid.configure(scry,sticky="wns")
  tkgrid.configure(scrx,sticky="we")
  
  tkgrid.configure(scry2,sticky="ns")
  tkgrid.configure(scrx2,sticky="ew")
  
  tkgrid.configure(scry3,sticky="wns")
  tkgrid.configure(scrx3,sticky="ew")
  
  ###########################################
  ########## INITIALIZE WIDGETS #############
  ###########################################
  current<-getwd()
  ## get directories
  all<-list.files(current,full.names = TRUE)
  init<-c("..",all[file.info(all)$isdir])
  wind<-Sys.info()['sysname']=="Windows"
  if(wind){
    possibleDrives<-paste(LETTERS,":/",sep="")
    driveInd<-unlist(lapply(possibleDrives,function(x){length(list.files(x))>0}))
    drives<-possibleDrives[driveInd]
    if(grepl("^\\\\",current)){init<-drives}
  }
  for (i in 1:length(init)){
    tkinsert(tl,"end",basename(init[i]))
  }
  init2<-all[!file.info(all)$isdir]
  
  for (i in 1:length(init2)){
    tkinsert(tl2,"end",basename(init2[i]))
  }
  ###########################################
  ########## KEY BINDINGS ###################
  ###########################################
  tkbind(reg, "<Return>", fileFilter)
  tkbind(tl, "<ButtonRelease-1>", changeDir)
  tkbind(tl2,"<ButtonRelease-3>",passFiles)
  tkbind(tl3, "<ButtonRelease-3>", deleteFiles)
  tkbind(tt, "<Destroy>", doneFunc)
  ###########################################
  ########## END ############################
  ###########################################
  tkwm.resizable(tt,FALSE,FALSE)
  tkfocus(tt)
  tkwait.variable(ret)
  tkdestroy(tt)
  options(show.error.messages = TRUE) 
  if(update){setwd(current)}
   return(tclFile2R(ret))
  
}