ic50<-function(){
  cat("NOTE: Preliminary change of the workspace directory to the folder\ncontaining the data will remarkably reduce the number of mouse clicks.\n")
  InhibPerc<-NULL
  defaultfiles.write()
  require(tcltk)
  files<-character(0)
  main<-tktoplevel()
  tkwm.title(main,"Automatic screen evaluation v1.4.1")
  
  #Add option boxes for plate format
  op96<-tkradiobutton(main,text="96 wells  ")
  op384<-tkradiobutton(main,text="384 wells ")
  plate<-tclVar("384")
  tkconfigure(op96,variable=plate,value="96")
  tkconfigure(op384,variable=plate,value="384")
  
  #Add textbox for number of plates 
  plates<-tclVar("2")
  PlatesInput<-tkentry(main,width=4,textvariable=plates)
                      
  #Add option boxes for normalization method
  opNrmMean <- tkradiobutton(main,text="Mean      ")
  opNrmSingle <- tkradiobutton(main,text="Single    ")
  normalize<-tclVar("single")
  tkconfigure(opNrmMean,variable=normalize,value="mean")
  tkconfigure(opNrmSingle,variable=normalize,value="single")
  
  #Add option boxes for graphics output
  opGrMean<-tkradiobutton(main,text="Mean      ")
  opGrFitted<-tkradiobutton(main,text="Fitted    ")
  opGrSingle<-tkradiobutton(main,text="Single    ")
  graphics<-tclVar("mean")
  tkconfigure(opGrMean,variable=graphics,value="mean")
  tkconfigure(opGrFitted,variable=graphics,value="fitted")
  tkconfigure(opGrSingle,variable=graphics,value="single")
  
  #Entry widgets for configuration files
  indir<-tclVar(".")
  measure<-tclVar("")
  control<-tclVar("")
  dilution<-tclVar("")
  outdir<-tclVar("results")
  IndirInput<-tkentry(main,width="45",textvariable=indir)
  MeasureInput<-tkentry(main,width="45",textvariable=measure)
  ControlInput<-tkentry(main,width="45",textvariable=control)
  DilutionInput<-tkentry(main,width="45",textvariable=dilution)
  OutdirInput<-tkentry(main,width="45",textvariable=outdir)
  
  #Create buttons for file specification
  onIndirBrowse<-function() tclvalue(indir)<-tclvalue(tkchooseDirectory(title="Input directory"))
  IndirBrowse.but<-tkbutton(main,text="Browse ...",width=10,command=onIndirBrowse)
  onMeasureBrowse<-function() tclvalue(measure)<-tclvalue(tkgetOpenFile(title="Measure file"))
  MeasureBrowse.but<-tkbutton(main,text="Browse ...",width=10,command=onMeasureBrowse)
  onControlBrowse<-function() tclvalue(control)<-tclvalue(tkgetOpenFile(title="Control file"))
  ControlBrowse.but<-tkbutton(main,text="Browse ...",width=10,command=onControlBrowse)
  onDilutionBrowse<-function() tclvalue(dilution)<-tclvalue(tkgetOpenFile(title="Dilution file"))
  DilutionBrowse.but<-tkbutton(main,text="Browse ...",width=10,command=onDilutionBrowse)
  onOutdirBrowse<-function() tclvalue(outdir)<-tclvalue(tkchooseDirectory(title="Output directory"))
  OutdirBrowse.but<-tkbutton(main,text="Browse ...",width=10,command=onOutdirBrowse)
  
  #Create Compute and Cancel buttons
  onCompute<-function(){
    if(!is.null(InhibPerc)) InhibPerc<-as.numeric(InhibPerc[,2])
    
    if(as.character(tclvalue(indir))=="") stop("No input data available.")
    if(as.character(tclvalue(outdir))==""){
      warning("No output directory chosen, writing to current directory.")
      outdir<-tclVar(".")
    }

    if(as.character(tclvalue(plate))=="384"){
      if(as.character(tclvalue(measure))=="") measure<-".last384_measure.txt"
      else measure<-as.character(tclvalue(measure))
      if(as.character(tclvalue(control))=="") control<-".last384_control.txt"
      else control<-as.character(tclvalue(control))
      if(as.character(tclvalue(dilution))=="") dilution<-".last384_dilution.txt"
      else dilution<-as.character(tclvalue(dilution))

      hts.384(indir=as.character(tclvalue(indir)),plates=as.numeric(as.character(tclvalue(plates))),
               measure=measure,control=control,dilution=dilution,inhib=InhibPerc,
               normalize=as.character(tclvalue(normalize)),
               graphics=as.character(tclvalue(graphics)),
               outdir=as.character(tclvalue(outdir))
              )
    }
    if(as.character(tclvalue(plate))=="96"){
      if(as.character(tclvalue(measure))=="") measure<-".last96_measure.txt"
      else measure<-as.character(tclvalue(measure))
      if(as.character(tclvalue(control))=="") control<-".last96_control.txt"
      else control<-as.character(tclvalue(control))
      if(as.character(tclvalue(dilution))=="") dilution<-".last96_dilution.txt"
      else dilution<-as.character(tclvalue(dilution))

      hts.96(indir=as.character(tclvalue(indir)),plates=as.numeric(as.character(tclvalue(plates))),
               measure=measure,control=control,dilution=dilution,inhib=InhibPerc,
               normalize=as.character(tclvalue(normalize)),
               graphics=as.character(tclvalue(graphics)),
               outdir=as.character(tclvalue(outdir))
              )
    }
    tkdestroy(main)
  }
  Compute.but<-tkbutton(main,text="Compute",width=15,command=onCompute)
  Cancel.but<-tkbutton(main,text="Cancel",width=15,command=function() tkdestroy(main))

  #Create buttons for design editing
  onMeasureEdit<-function() measure.edit(as.character(tclvalue(measure)),as.character(tclvalue(plate)))
  MeasureEdit.but<-tkbutton(main,text="Edit measure file",width=15,command=onMeasureEdit)
  onControlEdit<-function() control.edit(as.character(tclvalue(control)),as.character(tclvalue(plate)))
  ControlEdit.but<-tkbutton(main,text="Edit control file",width=15,command=onControlEdit)
  onDilutionEdit<-function() dilution.edit(as.character(tclvalue(dilution)),as.character(tclvalue(plate)))
  DilutionEdit.but<-tkbutton(main,text="Edit dilution file",width=15,command=onDilutionEdit)

  #Create button for change of inhibitory percentage
  onInhibPerc<-function(){
    if(as.character(tclvalue(dilution))==""){
      last<-paste(".last",as.character(tclvalue(plate)),"_dilution.txt",sep="")
      compounds<-levels(as.factor(read.delim(last,header=FALSE)[,1]))
    }
    else{
      compounds<-levels(as.factor(read.delim(as.character(tclvalue(dilution)),header=FALSE)[,1]))
    }
    InhibPerc<-data.frame(compounds,rep(0.5,length(compounds)))
    InhibPerc<<-edit(InhibPerc)
  }
  InhibPerc.but<-tkbutton(main,text="IC percentage ...",width=15,command=onInhibPerc)
  
  #Design GUI window
  tkgrid(tklabel(main,text=" "),row=1,column=1)
  tkgrid(tklabel(main,text="Input directory"),row=2,column=1)
  tkgrid(IndirInput,row=2,column=2,columnspan=3)
  tkgrid(IndirBrowse.but,row=2,column=5)
  
  tkgrid(tklabel(main,text=" "),row=3,column=1)
  tkgrid(tklabel(main,text="Output directory"),row=4,column=1)
  tkgrid(OutdirInput,row=4,column=2,columnspan=3)
  tkgrid(OutdirBrowse.but,row=4,column=5)
  
  tkgrid(tklabel(main,text=" "),row=5,column=1)
  tkgrid(tklabel(main,text="Plate format"),row=6,column=1)
  tkgrid(op96,row=6,column=2)
  tkgrid(op384,row=6,column=3)

  tkgrid(tklabel(main,text=" "),row=7,column=1)
  tkgrid(tklabel(main,text="Plates/measurement"),row=8,column=1)
  tkgrid(PlatesInput,row=8,column=2)
  tkgrid(InhibPerc.but,row=8,column=3)

  tkgrid(tklabel(main,text=" "),row=9,column=1)
  tkgrid(tklabel(main,text="Normalize by"),row=10,column=1)
  tkgrid(opNrmMean,row=10,column=2)
  tkgrid(opNrmSingle,row=10,column=3)

  tkgrid(tklabel(main,text=" "),row=11,column=1)
  tkgrid(tklabel(main,text="Graphical output"),row=12,column=1)
  tkgrid(opGrMean,row=12,column=2)
  tkgrid(opGrFitted,row=12,column=3)
  tkgrid(opGrSingle,row=12,column=4)
  
  tkgrid(tklabel(main,text=" "),row=13,column=1)
  tkgrid(tklabel(main,text="Measure file"),row=14,column=1)
  tkgrid(MeasureInput,row=14,column=2,columnspan=3)
  tkgrid(MeasureBrowse.but,row=14,column=5)
  tkgrid(MeasureEdit.but,row=14,column=6)
                              
  tkgrid(tklabel(main,text=" "),row=15,column=1)
  tkgrid(tklabel(main,text="Control file"),row=16,column=1)
  tkgrid(ControlInput,row=16,column=2,columnspan=3)
  tkgrid(ControlBrowse.but,row=16,column=5)
  tkgrid(ControlEdit.but,row=16,column=6)

  tkgrid(tklabel(main,text=" "),row=17,column=1)
  tkgrid(tklabel(main,text="Dilution file"),row=18,column=1)
  tkgrid(DilutionInput,row=18,column=2,columnspan=3)
  tkgrid(DilutionBrowse.but,row=18,column=5)
  tkgrid(DilutionEdit.but,row=18,column=6)
                     
  tkgrid(tklabel(main,text=" "),row=19,column=1)
  tkgrid(Compute.but,row=20,column=2)
  tkgrid(Cancel.but,row=20,column=3)
  tkgrid(tklabel(main,text=" "),row=21,column=1)
  
  tkfocus(main)
  tkwait.window(main)
}
