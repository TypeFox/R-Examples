GUI<-function(){
  disp<-function(){
    viewR()
  }
  fiachGUI<-function(){
    base<-tktoplevel()
    
    tktitle(base) <- "FIACH" 
    master<-tkframe(base,width=400,height=400*1.61,relief="sunken")
    top<-tkframe(master,borderwidth=0,width=380,height=190)
    middle<-tkframe(master,borderwidth=0,width=380,height=190)
    bottom<-tkframe(master,borderwidth=0,width=380,height=190)
    ################################
    ### VARIABLES ##################
    ################################
    field<-tclVar(init = "")
    tr<-tclVar(init="")
    te<-tclVar(init="")
    maxgap<-tclVar(init="1")
    freq<-tclVar(init="128")
    nmads<-tclVar(init="1.96")
    save<-tclVar(init=TRUE)
    #################################
    ####### FUNCTIONS ###############
    #################################
    conditionSelectR<-function(){ 
      
      currentId<-as.character(tkselect(treeWidget))
      if(!is.null(currentId) && grepl("(rp|ff)",currentId) &&length(currentId)==1){
        
        
        allSubTags <-as.character(tcl(treeWidget ,"children" ,""))
        numSub<-strtoi(gsub("s","",allSubTags))
        allffTags<-paste("ff",numSub,sep="")
        allrpTags<-paste("rp",numSub,sep="")
        allContTags<-c(allffTags,allrpTags)
        
        currentIds<-unlist(lapply(allContTags,function(x){as.character(tcl(treeWidget,"children",x))}))
        if(length(currentIds)>0){
          currentTags<-lapply(currentIds,function(x){as.character(tcl(treeWidget,"item",x))[10]})
          biggestTag<-max(strtoi(gsub("[a-z][a-z][0-9].*file","",currentTags)))
          
          
          files<-selectR()
          if(length(files)>0){
            tag<-paste(currentId,"file",(1:length(files))+biggestTag,sep="")
            for(i in 1:length(files)){
              tkinsert(treeWidget,currentId,"end",values= as.tclObj(files[i]),tag=tag[i])
              tktag.configure(treeWidget, tag[i],background = "white",foreground="black")
            }
          }
        }else{
          files<-selectR()
          if(length(files)>0){
            tag<-paste(currentId,"file",1:length(files),sep="")
            for(i in 1:length(files)){
              tkinsert(treeWidget,currentId,"end",values= as.tclObj(files[i]),tag=tag[i])
              tktag.configure(treeWidget, tag[i],background = "white")
            }
          }
          
        }
      }
    }
    deleteFiles<-function(){
      currentId<-as.character(tkselect(treeWidget))
      
      if(length(currentId)>0 && !grepl("(rp|ff)",currentId)){tkdelete(treeWidget,currentId)}
    }
    newSubject<-function(){
      children <- as.character(tcl(treeWidget , "children" , "" ))
      if(length(children)>0){
        lastChild<-children[length(children)]
        numSub<-strtoi(gsub("s",replacement = "",x = lastChild))+1
        id<-paste(c("s","ff","rp"),numSub,sep="")
        tex<-paste("Subject",numSub)
        
        tkinsert(treeWidget,"","end",id=id[1],text=tex,tag=id[1])
        tktag.configure(treeWidget, id[1],background = "white",foreground="black")
        tkinsert(treeWidget,id[1],"end",id=id[2],text="Functional Files",tag=id[2])
        tktag.configure(treeWidget, id[2],background = "white",foreground="black")
        tkinsert(treeWidget,id[1],"end",id=id[3],text="Realignment Parameters",tag=id[3])
        tktag.configure(treeWidget, id[3],background = "white",foreground="black")
      }else{
        tkinsert(treeWidget,"","end",id="s1",text="Subject 1",tag="s1")
        tktag.configure(treeWidget, "s1",background = "white",foreground="black")
        tkinsert(treeWidget,"s1","end",id="ff1",text="Functional Files",tag="ff1")
        tktag.configure(treeWidget, "ff1",background = "white",foreground="black")
        tkinsert(treeWidget,"s1","end",id="rp1",text="Realignment Parameters",tag="rp1")
        tktag.configure(treeWidget, "rp1",background = "white",foreground="black")
      }
    }
    
    highlight<-function(){
      
      allSubTags <-as.character(tcl(treeWidget ,"children" ,""))
      numSub<-strtoi(gsub("s","",allSubTags))
      allffTags<-paste("ff",numSub,sep="")
      allrpTags<-paste("rp",numSub,sep="")
      allContTags<-c(allffTags,allrpTags)
      nFiles<-unlist(lapply(allContTags,function(x){length(as.character(tcl(treeWidget,"children",x)))}))
      
      currentIds<-unlist(lapply(allContTags,function(x){as.character(tcl(treeWidget,"children",x))}))
      
      
      if(sum(nFiles)>0){
        fileTags<-unlist(lapply(currentIds,function(x){as.character(tcl(treeWidget,"item",x))[10]}))
        biggestTag<-max(strtoi(gsub("[a-z][a-z][0-9].*file","",fileTags)))
        allTags<-c(allSubTags,allContTags,fileTags)
      }else{allTags<-c(allSubTags,allContTags)
      }
      for (i in 1:length(allTags)) {
        tktag.configure(
          widget = treeWidget,allTags[i],background = "white",foreground = "black"
        )
      }
      sel<- as.character(tkselect(treeWidget))
      if(length(sel)>0){
        for(i in 1:length(sel)){
          cs<-as.character(tcl(treeWidget, "item" ,as.character(sel)[i]))
          tag<-cs[length(cs)]
          conf<-paste(as.character(tktag.configure(treeWidget,tag)),collapse=" ")
          if(grepl("background blue",conf)){
            tktag.configure(widget = treeWidget,tag,background="white",foreground="black")
          }else{
            tktag.configure(widget = treeWidget,tag,background="cornflowerblue",foreground="white")
          }
        }
      }
    }
    openScript<-function(){
      file<-selectR()
      if(length(file)==1){
        children <- as.character(tcl(treeWidget , "children" , "" ))
        isfiachFile<-grepl("fiachBatch",file)
        
        if(isfiachFile){
          tkdelete(treeWidget,children)
          test<-readLines(file,n = 1)
        }else{
            stop("This does not appear to be a script generated by FIACH or it has been edited.")
          }
          isFIACH<-grepl("FIACH script generated on:",test)
        if(isFIACH){
          txt<-readLines(file)
        }else{
          print("This does not appear to be a script generated by FIACH or it has been edited.")
        }
        
        startData<-grep("input<-list",txt)
        startRp<-grep("rp<-list",txt)
        ends<-grep("^\\)$",txt)
        
        data<-eval(parse(text = txt[startData:ends[1]]))
        rps<-eval(parse(text = txt[startRp:ends[2]]))
        
        
        ffNames<-c()
        rpNames<-c()
        for(i in 1:length(data)){
          newSubject()
          ffNames[i]<-paste("ff",i,sep="")
          rpNames[i]<-paste("rp",i,sep="")
        }
        for(i in 1:length(data)){
          lapply(data[[i]],function(x){tkinsert(treeWidget,ffNames[i],"end",values=x)})
          lapply(rps[[i]],function(x){tkinsert(treeWidget,rpNames[i],"end",values=x)})
        }
      }else{print("More than one file was selected. Please only select one file.")}
    }
    genSave<-function(){
      subjects<- as.character(tcl(treeWidget , "children" , "" ))
      if(length(subjects)==0){
        print("No Subjects have been defined")
      }else{
        numbers<-gsub("s","",subjects)
        funcNames<-paste("ff",numbers,sep="")
        funcIds<-lapply(funcNames,function(x){as.character(tcl(treeWidget,"children",x))})
        if(any(lapply(funcIds,length)==0)){
          print("Not all Subjects have functional files")
        }else{
          rpNames<-paste("rp",numbers,sep="")
          rpIds<-lapply(rpNames,function(x){as.character(tcl(treeWidget,"children",x))})
          
          B0<-as.numeric(tclvalue(field))
          TR<-as.numeric(tclvalue(tr))
          TE<-as.numeric(tclvalue(te))
          nMads<-as.numeric(tclvalue(nmads))
          FREQ<-as.numeric(tclvalue(freq))
          MAXGAP<-as.numeric(tclvalue(maxgap))
          
          if(is.na(B0)||is.na(TR)||is.na(TE)||is.na(nMads)||is.na(FREQ)||is.na(MAXGAP)){
            print("Not all study specific variables have been defined or were defined incorrectly.")
          }else{
            
            dat<-list()
            rp<-list()
            
            for(i in 1:length(funcIds)){
              if(length(funcIds[[i]])>0){
                dat[[i]]<-unlist(lapply(funcIds[[i]],function(x){as.character(tcl(treeWidget,"item",x,"-values"))}))
              }else{
                dat[[i]]<-NA
              }
            }
            
            
            for(i in 1:length(rpIds)){
              if(length(rpIds[[i]])>0){
                rp[[i]]<-unlist(lapply(rpIds[[i]],function(x){as.character(tcl(treeWidget,"item",x,"-values"))}))
              }else{
                rp[[i]]<-NA
              }
            }
            if(any(lapply(rp,length)>1)){
              print("More than one realignment parameter has been specfiied")
            }else{
              dateChunk<-paste("#FIACH script generated on:",Sys.time())
              B0Chunk<-paste("B0<-",B0,sep="")
              TRChunk<-paste("tr<-",TR,sep="")
              TEChunk<-paste("te<-",TE,sep="")
              TChunk<-paste("t<-boldContrast(B0,te,plot=FALSE)",sep="")
              madsChunk<-paste("nMads<-",nMads)
              freqChunk<-paste("freq<-",FREQ)
              gapChunk<-paste("maxgap<-",MAXGAP)
              processChunk<-paste(
                "for(i in 1:length(rp)){",
                "print(paste('Beginning subject',i))",
                "if(rp[[i]]=='NA'){fiach(input=input[[i]],tr=tr,t=t,freq=freq,maxgap=maxgap,nMads=nMads)}else{",
                "fiach(input=input[[i]],tr=tr,t=t,rp=rp[[i]],freq=freq,maxgap=maxgap,nMads=nMads)}",
                "print(paste('Completed subject',i))",
                "}",sep="\n")
              
              a<-"c("
              
              allFunc<-list()
              for(i in 1:length(dat)){
                meat<-paste(paste("'",dat[[i]],"'",sep=""),collapse=",\n")
                allFunc[[i]]<-paste(a,meat,")\n",collapse="\n") 
              }
              
              func<-paste(unlist(allFunc),collapse=", \n")
              funcChunk<-paste("input<-list(\n",func,"\n)")
              
              allRps<-list()
              for(i in 1:length(rp)){
                meat<-paste(paste("'",rp[[i]],"'",sep=""),collapse=",\n")
                allRps[[i]]<-paste(a,meat,")\n",collapse="\n") 
              }
              datarp<-paste(unlist(allRps),collapse=", \n")
              rpChunk<-paste("\n rp<-list(\n",datarp,"\n)")
              
              
              combine<-paste(dateChunk,B0Chunk,TEChunk,TRChunk,TChunk,madsChunk,freqChunk,gapChunk,funcChunk,rpChunk,processChunk,sep="\n")
              direc<-tkchooseDirectory()
              if(length(direc)>0){
              direc<-paste(as.character(direc),collapse=" ")
                output<-paste(direc,"/",gsub("(:|-| )","_",Sys.time()),"_","fiachBatch.r",sep="")
                writeLines(combine,output)
              }else{
                print("Save Cancelled")
              }
            }
          }
        }
      }
    }
    genRun<-function(){
      subjects<- as.character(tcl(treeWidget , "children" , "" ))
      if(length(subjects)==0){
        print("No Subjects have been defined")
      }else{
        numbers<-gsub("s","",subjects)
        funcNames<-paste("ff",numbers,sep="")
        funcIds<-lapply(funcNames,function(x){as.character(tcl(treeWidget,"children",x))})
        if(any(lapply(funcIds,length)==0)){
          print("Not all Subjects have functional files")
        }else{
          rpNames<-paste("rp",numbers,sep="")
          rpIds<-lapply(rpNames,function(x){as.character(tcl(treeWidget,"children",x))})
          
          B0<-as.numeric(tclvalue(field))
          TR<-as.numeric(tclvalue(tr))
          TE<-as.numeric(tclvalue(te))
          nMads<-as.numeric(tclvalue(nmads))
          FREQ<-as.numeric(tclvalue(freq))
          MAXGAP<-as.numeric(tclvalue(maxgap))
          
          if(is.na(B0)||is.na(TR)||is.na(TE)||is.na(nMads)||is.na(FREQ)||is.na(MAXGAP)){
            print("Not all study specific variables have been defined or were defined incorrectly.")
          }else{
            
            dat<-list()
            rp<-list()
            
            for(i in 1:length(funcIds)){
              if(length(funcIds[[i]])>0){
                dat[[i]]<-unlist(lapply(funcIds[[i]],function(x){as.character(tcl(treeWidget,"item",x,"-values"))}))
              }else{
                dat[[i]]<-NA
              }
            }
            
            
            for(i in 1:length(rpIds)){
              if(length(rpIds[[i]])>0){
                rp[[i]]<-unlist(lapply(rpIds[[i]],function(x){as.character(tcl(treeWidget,"item",x,"-values"))}))
              }else{
                rp[[i]]<-NA
              }
            }
            if(any(lapply(rp,length)>1)){
              print("More than one realignment parameter has been specfiied")
            }else{
              dateChunk<-paste("#FIACH script generated on:",Sys.time())
              B0Chunk<-paste("B0<-",B0,sep="")
              TRChunk<-paste("tr<-",TR,sep="")
              TEChunk<-paste("te<-",TE,sep="")
              TChunk<-paste("t<-boldContrast(B0,te,plot=FALSE)",sep="")
              madsChunk<-paste("nMads<-",nMads)
              freqChunk<-paste("freq<-",FREQ)
              gapChunk<-paste("maxgap<-",MAXGAP)
              processChunk<-paste(
                "for(i in 1:length(rp)){",
                "print(paste('Beginning subject',i))",
                "if(rp[[i]]=='NA'){fiach(input=input[[i]],tr=tr,t=t,freq=freq,maxgap=maxgap,nMads=nMads)}else{",
                "fiach(input=input[[i]],tr=tr,t=t,rp=rp[[i]],freq=freq,maxgap=maxgap,nMads=nMads)}",
                "print(paste('Completed subject',i))",
                "}",sep="\n")
              
              a<-"c("
              
              allFunc<-list()
              for(i in 1:length(dat)){
                meat<-paste(paste("'",dat[[i]],"'",sep=""),collapse=",\n")
                allFunc[[i]]<-paste(a,meat,")\n",collapse="\n") 
              }
              
              func<-paste(unlist(allFunc),collapse=", \n")
              funcChunk<-paste("input<-list(\n",func,"\n)")
              
              allRps<-list()
              for(i in 1:length(rp)){
                meat<-paste(paste("'",rp[[i]],"'",sep=""),collapse=",\n")
                allRps[[i]]<-paste(a,meat,")\n",collapse="\n") 
              }
              datarp<-paste(unlist(allRps),collapse=", \n")
              rpChunk<-paste("\n rp<-list(\n",datarp,"\n)")
              
              
              combine<-paste(dateChunk,B0Chunk,TEChunk,TRChunk,TChunk,madsChunk,freqChunk,gapChunk,funcChunk,rpChunk,processChunk,sep="\n")
              
              
              
              output<-file.path(normalizePath(tempdir(),winslash = "/"),paste(gsub("(:|-| )","_",Sys.time()),"_","fiachBatch",".r",sep=""))
              writeLines(combine,output)
              source(output)
              unlink(output)
            }
          }
        }
      }
    }
    output<-function(){
      base<-tktoplevel()
      tktitle(base)<-"Output Files"
      textF<-paste("Nothing is returned to R but the filtered files are written to the directory they came from with the prefix filt_. A directory is also created to store the various diagnostic images and plots produced by this method(median, rTSNR, mask and rTSNR histogram). The regressors to be used in further analysis are in the noise_basis6 file. The global median signal is also outputted in gs.txt file. The framewise displacement is also outputted appended to the noise regressors in fd_noise.txt if movement regressors are supplied. The peak rTSNR and the % data changed by FIACH is in the metrics.txt file." ,sep="")
      mean<-tktext(parent = base,width = 80,height=20,wrap="word")
      tkgrid(mean)
      tkinsert(mean,"end",textF)
      tkconfigure(mean,state="disabled")
      tkwm.resizable(base,FALSE,FALSE)
    }
    leave<-function(){
      tkconfigure(HELP,state="normal")
      tkdelete(HELP,"1.0","end")
      tkinsert(HELP,"1.0","Hover over a section to see the help file.")
      tkconfigure(HELP,state="disabled")
    }
    helpB0<-function(){
      tkconfigure(HELP,state="normal")
      tkdelete(HELP,"1.0","end")
      tkinsert(HELP,"1.0","B0 is the field strength of the MRI. This value is specified in Tesla.")
      tkconfigure(HELP,state="disabled")
    }
    helpTR<-function(){
      tkconfigure(HELP,state="normal")
      tkdelete(HELP,"1.0","end")
      tkinsert(HELP,"1.0","TR is the time between volumes. This value is specified in seconds.")
      tkconfigure(HELP,state="disabled")
    }
    helpTE<-function(){
      tkconfigure(HELP,state="normal")
      tkdelete(HELP,"1.0","end")
      tkinsert(HELP,"1.0","TE is the echo time. This value is specified in milliseconds.")
      tkconfigure(HELP,state="disabled")
    }
    helpMADS<-function(){
      tkconfigure(HELP,state="normal")
      tkdelete(HELP,"1.0","end")
      tkinsert(HELP,"1.0","No. MADS is the number of median absolute deviations added to the filtering threshold. Lower values produce more aggressive filtering.")
      tkconfigure(HELP,state="disabled")
    }
    helpHP<-function(){
      tkconfigure(HELP,state="normal")
      tkdelete(HELP,"1.0","end")
      tkinsert(HELP,"1.0","High-Pass cutoff is the period(seconds) beyond which low frequency noise is removed. This default value is the same as in SPM and should rarely be changed.")
      tkconfigure(HELP,state="disabled")
    }
    helpMG<-function(){
      tkconfigure(HELP,state="normal")
      tkdelete(HELP,"1.0","end")
      tkinsert(HELP,"1.0","Max Gap is the number of consecutive volumes over which interpolation is allowed. For designs with TR<2 seconds this value should be approximately ceiling(2/TR).")
      tkconfigure(HELP,state="disabled")
    }
    helpRS<-function(){
      tkconfigure(HELP,state="normal")
      tkdelete(HELP,"1.0","end")
      tkinsert(HELP,"1.0","The script is run using the information entered in the GUI.")
      tkconfigure(HELP,state="disabled")
    }
    helpSS<-function(){
      tkconfigure(HELP,state="normal")
      tkdelete(HELP,"1.0","end")
      tkinsert(HELP,"1.0","Save script saves the information entered in the GUI into a directory of choice. The script is not run.")
      tkconfigure(HELP,state="disabled")
    }
    helpSelect<-function(){
      tkconfigure(HELP,state="normal")
      tkdelete(HELP,"1.0","end")
      help<-paste("1. Highlight Functional Files or Realignment Parameters and double right click to add files.",
                  "2. Realignment files are optional. If not specified framewise displacement will not be calculated.",
                  "3. Functional Files are not optional.",
                  "4. Functional Files must be realigned.",
                  "5. Ctrl+C = Remove files or subjects",
                  "6. Ctrl+V = Add new Subject",sep="\n"
      )
      tkinsert(HELP,"1.0",help)
      tkconfigure(HELP,state="disabled")
    }
    helpOPEN<-function(){
      tkconfigure(HELP,state="normal")
      tkdelete(HELP,"1.0","end")
      tkinsert(HELP,"1.0","Open script opens a FIACH script. WARNING!! pressing this button removes all currently specified files.")
      tkconfigure(HELP,state="disabled")
    }
    
    #################################
    ###### WIDGETS ##################
    #################################
    addScrollbars<-function(parent,widget){
      xscr <- ttkscrollbar ( parent , orient  = "horizontal",command = function(...) tkxview(widget, ...))
      yscr <- ttkscrollbar ( parent , orient = "vertical" ,command = function(...) tkyview(widget,...))
      tkconfigure(widget, xscrollcommand = function(...)tkset(xscr,...), yscrollcommand = function(...)tkset(yscr,...))
      tkgrid(widget, row = 0 , column = 0 , sticky = "news" )
      tkgrid (yscr, row = 0, column = 1 , sticky = "ns" )
      tkgrid (xscr, row = 1, column = 0 , sticky = "ew" )
      tkgrid.columnconfigure(parent, 0 , weight = 1 )
      tkgrid.rowconfigure(parent, 0, weight = 1 )
      
    }                                                               
    
    B0<-tkentry(parent = top,textvariable=field,width=5)
    B0LAB<-tklabel(parent=top,text="B0")
    TR<-tkentry(parent = top,textvariable=tr,width=5)
    TRLAB<-tklabel(parent=top,text="TR")
    TE<-tkentry(parent = top,textvariable=te,width=5)
    TELAB<-tklabel(parent=top,text="TE")
    MAXGAP<-tkentry(parent = top,textvariable=maxgap,width=5)
    MAXGAPLAB<-tklabel(parent=top,text="Max Gap")
    FREQ<-tkentry(parent = top,textvariable=freq,width=5)
    FREQLAB<-tklabel(parent=top,text="High-Pass Cutoff")
    nMADS<-tkentry(parent = top,textvariable=nmads,width=5)
    nMADSLAB<-tklabel(parent=top,text="No. MADS")
    SAVE<-tkbutton(parent = top,text="Save Script",command=genSave)
    OUTPUT<-tkbutton(parent=top,text="What is the output?",command=output)
    RUN<-tkbutton(parent = top,text="Run Script",command=genRun)
    OPEN<-tkbutton(parent = top,text="Open Script",command=openScript)
    HELP<-tktext(parent = bottom,width=45,height= 11,wrap="word")
    HELPLAB<-tklabel(parent=master,text="Help Section")
    STUDYLAB<-tklabel(parent=master,text="Study Specific Variables")
    SUBLAB<-tklabel(parent=master,text="Subject Specific Variables")
    tkinsert(HELP,"1.0","Hover over a section to see the help file")
    tkconfigure(HELP,state="disabled")
    
    tkgrid(master)
    tkplace(top,x=10,y=25)
    tkplace(middle,x=10,y=235,width=380,height=200)
    tkplace(bottom,x=10,y=455)
    
    
    
    treeWidget <-ttktreeview(middle,columns = 2,height = 8,selectmode="extended")
    xscr <- ttkscrollbar (middle , orient  = "horizontal",command = function(...) tkxview(treeWidget, ...))
    yscr <- ttkscrollbar (middle, orient = "vertical" ,command = function(...) tkyview(treeWidget,...))
    tkconfigure(treeWidget, xscrollcommand = function(...)tkset(xscr,...), yscrollcommand = function(...)tkset(yscr,...))
    tkplace(treeWidget,x=0,y=0,width=360,height=180)
    
    tkplace(yscr, x=360, y=0,width=20,height=200)
    tkplace(xscr,x=0,y=180,width=360,height=20)
    tcl(treeWidget , "column" ,2,width = 1000, stretch = FALSE , anchor = "w" )
    tcl(treeWidget , "heading" , 2 , text = "Files" , anchor = "w" )
    
    tkinsert(treeWidget,"","end",id="s1",text="Subject 1",tag="s1")
    tktag.configure(treeWidget, "s1",background = "white",foreground="black")
    tkinsert(treeWidget,"s1","end",id="ff1",text="Functional Files",tag="ff1")
    tktag.configure(treeWidget, "ff1",background = "white",foreground="black")
    tkinsert(treeWidget,"s1","end",id="rp1",text="Realignment Parameters",tag="rp1")
    tktag.configure(treeWidget, "rp1",background = "white",foreground="black")
    
    tkbind(treeWidget, "<Double-Button-3>", conditionSelectR)
    tkbind(treeWidget, "<Control-c>", deleteFiles)
    tkbind(treeWidget, "<Control-v>", newSubject)
    tkbind(B0, "<Enter>", helpB0)
    tkbind(B0, "<Leave>", leave)
    tkbind(TR, "<Enter>", helpTR)
    tkbind(TR, "<Leave>", leave)
    tkbind(TE, "<Enter>", helpTE)
    tkbind(TE, "<Leave>", leave)
    tkbind(nMADS, "<Enter>", helpMADS)
    tkbind(nMADS, "<Leave>", leave)
    tkbind(MAXGAP, "<Enter>", helpMG)
    tkbind(MAXGAP, "<Leave>", leave)
    tkbind(FREQ, "<Enter>", helpHP)
    tkbind(FREQ, "<Leave>", leave)
    tkbind(RUN, "<Enter>", helpRS)
    tkbind(RUN, "<Leave>", leave)
    tkbind(SAVE, "<Enter>", helpSS)
    tkbind(SAVE, "<Leave>", leave)
    tkbind(OPEN, "<Enter>", helpOPEN)
    tkbind(OPEN, "<Leave>", leave)
    tkbind(treeWidget, "<Enter>", helpSelect)
    tkbind(treeWidget,"<ButtonRelease-1>",highlight)
    
    
    ####################################
    ######### PLACEMENT ###################
    #####################################
    tkplace(SUBLAB,x=10,y=211)
    tkplace(STUDYLAB,x=10,y=5)
    tkplace(HELPLAB,x=10,y=432)
    tkplace(B0,x=30,y=40)
    tkplace(B0LAB,x=10,y=40)
    tkplace(TR,x=30,y=90)
    tkplace(TRLAB,x=10,y=90)
    tkplace(TE,x=30,y=140)
    tkplace(TELAB,x=10,y=140)
    
    tkplace(nMADS,x=200,y=40)
    tkplace(nMADSLAB,x=135,y=40)
    tkplace(FREQ,x=200,y=90)
    tkplace(FREQLAB,x=95,y=90)
    tkplace(MAXGAP,x=200,y=140)
    tkplace(MAXGAPLAB,x=138,y=140)
    
    tkplace(OUTPUT,x=260,y=40)
    tkplace(OPEN,x=260,y=75)
    tkplace(RUN,x=260,y=105)
    tkplace(SAVE,x=260,y=140)
    
    tkgrid(HELP)
    tkwm.resizable(base,FALSE,FALSE)
  }
  
  myFont<-tkfont.create(family="calibri",size=43,weight="bold")
  base <- tktoplevel()
  tkfocus()
  tktitle(base) <- "FIACH"
  master<-tkframe(base,width=474,height=600,relief="sunken",background="lightblue")
  top<-tkframe(base,borderwidth=0,width=474,height=282)
  file<-system.file("extdata","fiach.gif",package = "FIACH")
  
  meaningHelp<-function(){
    base<-tktoplevel()
    tktitle(base)<-"The Meaning of FIACH"
    textF<-paste("Officially FIACH (approximately pronounced: FEE-UCK; IPA: \u02C8f\u02B2i\u0259x) is an acronym for Functional Image Artefact Correction Heuristic. However the word fiach in  the Irish language has a number of meanings, some of which are listed below.",
                "\nFiach = To Hunt",
                "Fiach = Raven",
                "Gadhar Fiaigh = Hunting Dog",
                "\nThere are three reasons that a black labrador is the emblem of this sofware.", 
                "1. The dog is a raven coloured dog.",
                "2. He was descended from a  long line of hunting dogs.", 
                "3. His name is also Fiach.",
                sep="\n")
    mean<-tktext(parent = base,width = 80,height=20,wrap="word")
    tkgrid(mean)
    tkinsert(mean,"end",textF)
    tkconfigure(mean,state="disabled")
    tkwm.resizable(base,FALSE,FALSE)
  }
  
  wtext<-function(){
    base<-tktoplevel()
    tktitle(base)<-"Warranty"
    textF<-paste(" Copyright (C) <2015>  <Tim Tierney>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA")
    mean<-tktext(parent = base,width = 80,height=20,wrap="word")
    tkgrid(mean)
    tkinsert(mean,"end",textF)
    tkconfigure(mean,state="disabled")
    tkwm.resizable(base,FALSE,FALSE)
  }
  
  
  labelText <- tclVar("Developer: Tim Tierney")
  by <- tklabel(master,text=tclvalue(labelText),background="lightblue")
  tkconfigure(by,textvariable=labelText)
  
  tcl("image","create","photo", "imageID", file=file)
  l <- tklabel(top, image="imageID", compound="center")
  viewBut<-tkbutton(parent = master,text="Display",width=13,height=3,state="normal",command=disp)
  fiachBut<-tkbutton(parent = master,text="FIACH",width=13,height=3,state="normal",command=fiachGUI)
  meaning<-tkbutton(parent = top,text="The Meaning of FIACH",width=20,height=1,state="normal",command=meaningHelp)
  warranty<-tkbutton(parent = master,text="Warranty",width=7,height=1,state="normal",command=wtext)
  tkgrid(master,column=0,row=0)
  tkplace(top,x=0,y=0)
  tkplace(warranty,rely=.95,relx=.02)
  tkplace(by,relx=.65,y=570)
  tkgrid(l)
  tkplace(viewBut,relx=.2,rely=.7)
  tkplace(fiachBut,relx=.6,rely=.7)
  tkplace(meaning,relx=0,rely=.9)
  tkwm.resizable(base,FALSE,FALSE)
  tkfocus(base)
   }