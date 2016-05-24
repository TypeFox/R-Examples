ConstructDemo <- function(file=""){
  # require(tcltk) # 140306
  if(file==""){ cat("sorry, please call ConstructDemo with name of R-, rev- or a Snw-file\n");            
                return() } 
  if(  ( 0 == length(grep(".R$",file))) || !file.exists(file) ){ 
      # .R must be generated&   
      file <- sub("\\.R$",".rev",file)  # rev name of not existing file if .R file
      if(!file.exists(file)){
        file <- sub("\\.rev$","",file)  # base name of file if .rev file
        if(!file.exists(file)) file <- paste(file,".rev",sep="")  
        if(!file.exists(file)) file <- sub("\\.rev$",".Snw",file)  
        if(!file.exists(file)) file <- sub("\\.Snw$",".snw",file) 
        if(!file.exists(file)) file <- sub("\\.snw$",".nw",file)   
        if(!file.exists(file)) file <- sub("\\.nw$",".SNW",file)  
        if(!file.exists(file)) file <- sub("\\.SNW$","",file)   
        if(!file.exists(file)){
          cat(paste("sorry, error: File",file,"not found!!!"));return()
        }
      } 
      cat("function will work on file",file,"\n") 
      # construct a .R file
      try(tangleR(file))
      # file <- sub(".rev$","",file); file <- sub(".[sS]{0,1}[nN][wW]$","",file)
      file <- sub("\\.[a-zA-Z]+$","",file)
      file <- paste(file,".R",sep=""); cat(file,"generated\n")
  }


  workname <- file; chunks <- where <- ""
  txt <- scan(workname,"",sep="\n")
  fh<-function(file,no=1,start=TRUE){
    # initial code of demo function
    # require("tcltk"); where<-environment() # 140306
  }
; fh <- deparse(fh)
  fh<-c(fh[-length(fh)], "'require(\"tcltk\")'; where<-environment()")
  ft<-function(){
    # end of code of demo function
    ## activate start chunk
    if(start==TRUE){
      no.0<-"0"
      no.start.0<-grep(paste("^#",no.0,":$",sep=""),chunks)
      no.end.0<-grep(paste("^#:",no.0,"$",sep=""),chunks)
      code.0<-chunks[no.start.0:no.end.0]
      eval(parse(text=code.0),envir=where)
    }
    ## activate chunk no
    # eval(parse(text=code),envir=where)
    secno<-tclVar("1") # 0
    show.next.number<-function(...){
     no<-as.character(as.numeric(tclvalue(secno))+1)
     no.start<-grep(paste("^#",no,":$",sep=""),chunks)
     no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
     if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
          is.nan(no.end)||is.nan(no.start)){
          cat("# sorry, chunk number '",no,"' wrong!\n"); return()
     }
     ### cat("# aktueller chunk:",no,"\n")
     code<-paste(chunks[no.start:no.end],collapse="\n")
     if(0<length(code)) {
      tkdelete(ttext,"0.0","end")  
      tkinsert(ttext,"0.0",code)
      tclvalue(secno)<-as.character(no)
      }
    }
    show.back.number<-function(...){
     no<-as.character(as.numeric(tclvalue(secno))-1)
     no.start<-grep(paste("^#",no,":$",sep=""),chunks)
     no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
     if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
          is.nan(no.end)||is.nan(no.start)){
          cat("# sorry, chunk number '",no,"' wrong!\n"); return()
     }
     ### cat("# aktueller chunk:",no,"\n")
     code<-paste(chunks[no.start:no.end],collapse="\n")
     if(0<length(code)) {
      tkdelete(ttext,"0.0","end")  
      tkinsert(ttext,"0.0",code)
      tclvalue(secno)<-as.character(no)
     }
    }
    show.number<-function(...){
     no<-as.character(as.numeric(tclvalue(secno)))
     no.start<-grep(paste("^#",no,":$",sep=""),chunks)
     no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
     if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
          is.nan(no.end)||is.nan(no.start)){
          cat("# sorry, chunk number '",no,"' wrong!\n"); return()
     }
     ### cat("# aktueller chunk:",no,"\n")
     code<-paste(chunks[no.start:no.end],collapse="\n")
     if(0<length(code)) {
      tkdelete(ttext,"0.0","end")  
      tkinsert(ttext,"0.0",code)
      tclvalue(secno)<-as.character(no)
     }
    }
    eval.code<-function(...){
      code<-tclvalue(tkget(ttext,"0.0","end"))
      code.orig<-code<-unlist(strsplit(code,"\n"))
      code<-code[!substring(code,1,1)=="#"]
      ##?## code<-unlist(strsplit(code,";"))  ##110429 Fehler bei sep=";"
      if(length(code)==0){ cat("ok\n"); return() }
      result<-try(eval(parse(text=code),envir=where))
      code.orig<-sub("#([0-9]+):","##wnt-Code-Chunk:\\1-begin#",code.orig)
      code.orig<-sub("#:([0-9]+)","##wnt-Code-Chunk:\\1-end#",code.orig)
      h<-get("allcodechunks",envir=where)
      h<-c(h,paste("<","<*>",">=",sep=""),code.orig,"\n@\n")
      assign("allcodechunks",h,envir=where)
      code<-sub("^ *","",code)
      code<-code[nchar(code)>0]
      lexpr<-rev(code)[1]; lexpr<-substring(lexpr,1,4)
      if(length(code)==0||is.null(lexpr)||is.na(lexpr)) return()
      plot.res<-c("plot","boxp","par(","abli","pie(","hist","axis","show",
             "lsfi","pair","ylab","help",
             "qqli","qqno","qqpl","rug(","lege","segm","text","xlab", 
             "poin","line","titl","eda(","imag","vgl.","curv")
      if(any(plot.res==lexpr)){
         cat("Plot erstellt\n"); return()
      }
      if(is.null(result)||is.na(result)||lexpr=="prin"||lexpr=="cat("){
        cat("ok\n"); return() }
      if(is.list(result)&& length(names(result))> 0 && 
                                 names(result)[1]=="ID") return()
      ## if(is.list(result)&& TRUE) return()
      no<-as.character(as.numeric(tclvalue(secno)))
      cat("Result of code chunk",no,":\n") 
     if(class(result)=="try-error"){
       class(result)<-"character"
       cat(result,"\n")
     }else{
      print(result)
     }
     cat("ok\n")
     }
    exit.function<-function(){
       tkdestroy(top)
       filename<-tkgetSaveFile(filetypes="{{Paper Files} {.rev}}",
                                             title="Do you want to save the activated R statements?")
       if(!is.character(filename)) filename<-tclvalue(filename)
       if(filename==""){
         cat("Demo function stopped without saving\n")
         return()
       }
       if(0==length(grep("rev$",filename))) filename<-paste(filename,".rev",sep="")
       h<-get("allcodechunks",envir=where)
       try(cat(h,sep="\n",file=filename))
       cat(paste("Remark: activated statements saved in\n   ",filename,"\n"))
       return()
    }
    allcodechunks<-paste(
            "@\nReport of activated R-chunks from: ",date(),
            "\n(demo function constructed by relax (c) Peter Wolf 2007)\n\n  ", sep="")
    no<-0
    no.start<-grep(paste("^#",no,":$",sep=""),chunks)
    no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
    if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
         is.nan(no.end)||is.nan(no.start)){
         cat("# sorry, chunk number '",no,"' wrong!\n"); return()
    }
    ### cat("# aktueller chunk:",no,"\n")
    code<-paste(chunks[no.start:no.end],collapse="\n")
    h<-paste(rep("#",60),collapse="")
    code<-sub("#0:",h,code); code<-sub("#:0",h,code)
    allcodechunks<-c(allcodechunks,"\n@\n<<start>>=",code,"\n@")
    assign("allcodechunks",allcodechunks,envir=where)
    top<-tktoplevel()
    ttext<-tktext(top,height=19,background="#f7fffF",
                  font="-Adobe-courier-Medium-R-Normal--18-180-*")
    tf<-tkframe(top) 
    tkwm.title(top, "demo of file WoRkNaMe, constructed by relax (c) Peter Wolf 2011")
    tkpack(tf,side="bottom"); tkpack(tf,ttext,side="bottom",fill="both",expand="y") # 091026
    tkevent.add("<<Paste>>",   "<Control_L><v>")
  ## ok:  tkbind(ttext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
  if(substring(version$os,1,6)=="darwin"  ){
    mac.paste<-function(...){  ## Ctrl-V  # mac-PASTE
       try({.Tcl("clipboard append hello"); .Tcl("clipboard clear")
            news<-base::scan(file=pipe("pbpaste","r"),what="",sep="\n",blank.lines.skip=FALSE)
            tkinsert(ttext,"insert",paste(news,collapse="\n"))})
            tksee(ttext,"insert - 7 lines");  tksee(ttext,"insert + 7 lines") #090706 
    } 
    tkbind(ttext,"<Control_L><v>",mac.paste)
    mac.copy<-function(...){  ## Ctrl-C  # mac-COPY
       news<-""
       try(news<-tclvalue(.Tcl("if {[catch {clipboard get}]} {set aa empty} {set aa full}")))
       if(news=="empty") return()
       try({news<-tclvalue(.Tcl("set aaa [selection get -selection CLIPBOARD]"))
            tmp.file.name <- tempfile("rt-tmp")
            base::cat(news,file=tmp.file.name); system(paste("pbcopy < ",tmp.file.name))
            .Tcl("clipboard append hello"); .Tcl("clipboard clear")})
    }
    tkbind(ttext,"<Control_L><c>",mac.copy)
    tkevent.add("<<extract>>",  "<Control_L><c><KeyRelease>")   # mac-extract
    tkbind(ttext,"<<extract>>",mac.copy)
  }else{
    tkevent.add("<<Paste>>",   "<Control_L><v>")
    tkbind(ttext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
  }


    bexit<-tkbutton(tf,text="QUIT",width=9)
    beval<-tkbutton(tf,text="EVALUATE",width=9)
    bnext<-tkbutton(tf,text="NEXT",width=9)
    bback<-tkbutton(tf,text="BACK",width=9)
    tkbind(top, "<<EvalRCode>>", eval.code)
    if( ! substring(version$os,1,6)=="darwin"  ) {
      tkbind(top, "<<Next>>", show.next.number)
      tkbind(top, "<<Back>>", show.back.number)
    }
    lno  <-tkentry(tf,textvariable=secno,width=9)
    linfo<-tklabel(tf,text="chunk number:")
    tkpack(linfo,lno,beval,bnext,bback,bexit,side="left")
    tkconfigure(bexit,command=exit.function)
    tkconfigure(bnext,command=show.next.number)
    tkconfigure(bback,command=show.back.number)
    tkconfigure(beval,command=eval.code)
    # tkbind(lno,"<Return>",show.number)
    tkbind(lno,"<KeyRelease>",show.number)
    tclvalue(secno)<-as.character(no)
    show.number()
    ### tkwait.window(top)
  }

  ft <- deparse(ft)[-(1:2)]; ft <- sub("WoRkNaMe",workname,ft)
  fns.name <- gsub("[^A-Za-z]","",sub(".R$",".demo",workname))
  txt <- gsub("\\\\","\\\\\\\\",txt); txt <- gsub("\"","\\\\\"",txt)
  txt <- paste(',"',txt,'"',sep=""); txt[1]<-substring(txt[1],2)
  txt <- c(paste(fns.name,"<-"),fh,"chunks<-c(",txt,")",ft
           ,paste("cat(\"Demo will be started by > ",fns.name,"()\\n\")",sep="")
           ,paste(fns.name,"()\n",sep="")
      ) 
  workname<-sub(".R$",".demo.R",workname)
  cat(txt,file=workname,sep="\n")
  cat("Demo saved in file ",workname,", load demo by source(\"",workname,"\") \n",sep="")
}

CodeChunkPlayer<-function(file=""){
  # require(tcltk) # 140306
  if(file==""){ cat("sorry, please call player with name of R-, rev- or a Snw-file\n");            
                return() } 
  if(  ( 0 == length(grep(".R$",file))) || !file.exists(file) ){ 
      # .R must be generated&   
      file <- sub("\\.R$",".rev",file)  # rev name of not existing file if .R file
      if(!file.exists(file)){
        file <- sub("\\.rev$","",file)  # base name of file if .rev file
        if(!file.exists(file)) file <- paste(file,".rev",sep="")  
        if(!file.exists(file)) file <- sub("\\.rev$",".Snw",file)  
        if(!file.exists(file)) file <- sub("\\.Snw$",".snw",file) 
        if(!file.exists(file)) file <- sub("\\.snw$",".nw",file)   
        if(!file.exists(file)) file <- sub("\\.nw$",".SNW",file)  
        if(!file.exists(file)) file <- sub("\\.SNW$","",file)   
        if(!file.exists(file)){
          cat(paste("sorry, error: File",file,"not found!!!"));return()
        }
      } 
      cat("function will work on file",file,"\n") 
      # construct a .R file
      try(tangleR(file))
      # file <- sub(".rev$","",file); file <- sub(".[sS]{0,1}[nN][wW]$","",file)
      file <- sub("\\.[a-zA-Z]+$","",file)
      file <- paste(file,".R",sep=""); cat(file,"generated\n")
  }

  if( !file.exists(file) ){
    cat("sorry, .R-file not found or unable to construct .R-file\n");return() 
  } 
  workname <- file; chunks<-""; where<-""
  txt<-scan(workname,"",sep="\n")
  # construct structure of function
  fh<-function(file,no=1,start=TRUE){
    # initial code of demo function
    # require("tcltk"); where<-environment() # 140306
  }

  fh<-deparse(fh)
  fh<-c(fh[-length(fh)], "'require(\"tcltk\")'; where<-environment()")
  ft<-function(){
    # end of code of demo function
    ## activate start chunk
    if(start==TRUE){
      no.0<-"0"
      no.start.0<-grep(paste("^#",no.0,":$",sep=""),chunks)
      no.end.0<-grep(paste("^#:",no.0,"$",sep=""),chunks)
      code.0<-chunks[no.start.0:no.end.0]
      eval(parse(text=code.0),envir=where)
    }
    ## activate chunk no
    # eval(parse(text=code),envir=where)
    secno<-tclVar("1") # 0
    show.next.number<-function(...){
     no<-as.character(as.numeric(tclvalue(secno))+1)
     no.start<-grep(paste("^#",no,":$",sep=""),chunks)
     no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
     if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
          is.nan(no.end)||is.nan(no.start)){
          cat("# sorry, chunk number '",no,"' wrong!\n"); return()
     }
     ### cat("# aktueller chunk:",no,"\n")
     code<-paste(chunks[no.start:no.end],collapse="\n")
     if(0<length(code)) {
      tkdelete(ttext,"0.0","end")  
      tkinsert(ttext,"0.0",code)
      tclvalue(secno)<-as.character(no)
      }
    }
    show.back.number<-function(...){
     no<-as.character(as.numeric(tclvalue(secno))-1)
     no.start<-grep(paste("^#",no,":$",sep=""),chunks)
     no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
     if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
          is.nan(no.end)||is.nan(no.start)){
          cat("# sorry, chunk number '",no,"' wrong!\n"); return()
     }
     ### cat("# aktueller chunk:",no,"\n")
     code<-paste(chunks[no.start:no.end],collapse="\n")
     if(0<length(code)) {
      tkdelete(ttext,"0.0","end")  
      tkinsert(ttext,"0.0",code)
      tclvalue(secno)<-as.character(no)
     }
    }
    show.number<-function(...){
     no<-as.character(as.numeric(tclvalue(secno)))
     no.start<-grep(paste("^#",no,":$",sep=""),chunks)
     no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
     if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
          is.nan(no.end)||is.nan(no.start)){
          cat("# sorry, chunk number '",no,"' wrong!\n"); return()
     }
     ### cat("# aktueller chunk:",no,"\n")
     code<-paste(chunks[no.start:no.end],collapse="\n")
     if(0<length(code)) {
      tkdelete(ttext,"0.0","end")  
      tkinsert(ttext,"0.0",code)
      tclvalue(secno)<-as.character(no)
     }
    }
    eval.code<-function(...){
      code<-tclvalue(tkget(ttext,"0.0","end"))
      code.orig<-code<-unlist(strsplit(code,"\n"))
      code<-code[!substring(code,1,1)=="#"]
      ##?## code<-unlist(strsplit(code,";"))  ##110429 Fehler bei sep=";"
      if(length(code)==0){ cat("ok\n"); return() }
      result<-try(eval(parse(text=code),envir=where))
      code.orig<-sub("#([0-9]+):","##wnt-Code-Chunk:\\1-begin#",code.orig)
      code.orig<-sub("#:([0-9]+)","##wnt-Code-Chunk:\\1-end#",code.orig)
      h<-get("allcodechunks",envir=where)
      h<-c(h,paste("<","<*>",">=",sep=""),code.orig,"\n@\n")
      assign("allcodechunks",h,envir=where)
      code<-sub("^ *","",code)
      code<-code[nchar(code)>0]
      lexpr<-rev(code)[1]; lexpr<-substring(lexpr,1,4)
      if(length(code)==0||is.null(lexpr)||is.na(lexpr)) return()
      plot.res<-c("plot","boxp","par(","abli","pie(","hist","axis","show",
             "lsfi","pair","ylab","help",
             "qqli","qqno","qqpl","rug(","lege","segm","text","xlab", 
             "poin","line","titl","eda(","imag","vgl.","curv")
      if(any(plot.res==lexpr)){
         cat("Plot erstellt\n"); return()
      }
      if(is.null(result)||is.na(result)||lexpr=="prin"||lexpr=="cat("){
        cat("ok\n"); return() }
      if(is.list(result)&& length(names(result))> 0 && 
                                 names(result)[1]=="ID") return()
      ## if(is.list(result)&& TRUE) return()
      no<-as.character(as.numeric(tclvalue(secno)))
      cat("Result of code chunk",no,":\n") 
     if(class(result)=="try-error"){
       class(result)<-"character"
       cat(result,"\n")
     }else{
      print(result)
     }
     cat("ok\n")
     }
    exit.function<-function(){
       tkdestroy(top)
       filename<-tkgetSaveFile(filetypes="{{Paper Files} {.rev}}",
                                             title="Do you want to save the activated R statements?")
       if(!is.character(filename)) filename<-tclvalue(filename)
       if(filename==""){
         cat("Demo function stopped without saving\n")
         return()
       }
       if(0==length(grep("rev$",filename))) filename<-paste(filename,".rev",sep="")
       h<-get("allcodechunks",envir=where)
       try(cat(h,sep="\n",file=filename))
       cat(paste("Remark: activated statements saved in\n   ",filename,"\n"))
       return()
    }
    allcodechunks<-paste(
            "@\nReport of activated R-chunks from: ",date(),
            "\n(demo function constructed by relax (c) Peter Wolf 2007)\n\n  ", sep="")
    no<-0
    no.start<-grep(paste("^#",no,":$",sep=""),chunks)
    no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
    if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
         is.nan(no.end)||is.nan(no.start)){
         cat("# sorry, chunk number '",no,"' wrong!\n"); return()
    }
    ### cat("# aktueller chunk:",no,"\n")
    code<-paste(chunks[no.start:no.end],collapse="\n")
    h<-paste(rep("#",60),collapse="")
    code<-sub("#0:",h,code); code<-sub("#:0",h,code)
    allcodechunks<-c(allcodechunks,"\n@\n<<start>>=",code,"\n@")
    assign("allcodechunks",allcodechunks,envir=where)
    top<-tktoplevel()
    ttext<-tktext(top,height=19,background="#f7fffF",
                  font="-Adobe-courier-Medium-R-Normal--18-180-*")
    tf<-tkframe(top) 
    tkwm.title(top, "demo of file WoRkNaMe, constructed by relax (c) Peter Wolf 2011")
    tkpack(tf,side="bottom"); tkpack(tf,ttext,side="bottom",fill="both",expand="y") # 091026
    tkevent.add("<<Paste>>",   "<Control_L><v>")
  ## ok:  tkbind(ttext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
  if(substring(version$os,1,6)=="darwin"  ){
    mac.paste<-function(...){  ## Ctrl-V  # mac-PASTE
       try({.Tcl("clipboard append hello"); .Tcl("clipboard clear")
            news<-base::scan(file=pipe("pbpaste","r"),what="",sep="\n",blank.lines.skip=FALSE)
            tkinsert(ttext,"insert",paste(news,collapse="\n"))})
            tksee(ttext,"insert - 7 lines");  tksee(ttext,"insert + 7 lines") #090706 
    } 
    tkbind(ttext,"<Control_L><v>",mac.paste)
    mac.copy<-function(...){  ## Ctrl-C  # mac-COPY
       news<-""
       try(news<-tclvalue(.Tcl("if {[catch {clipboard get}]} {set aa empty} {set aa full}")))
       if(news=="empty") return()
       try({news<-tclvalue(.Tcl("set aaa [selection get -selection CLIPBOARD]"))
            tmp.file.name <- tempfile("rt-tmp")
            base::cat(news,file=tmp.file.name); system(paste("pbcopy < ",tmp.file.name))
            .Tcl("clipboard append hello"); .Tcl("clipboard clear")})
    }
    tkbind(ttext,"<Control_L><c>",mac.copy)
    tkevent.add("<<extract>>",  "<Control_L><c><KeyRelease>")   # mac-extract
    tkbind(ttext,"<<extract>>",mac.copy)
  }else{
    tkevent.add("<<Paste>>",   "<Control_L><v>")
    tkbind(ttext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
  }


    bexit<-tkbutton(tf,text="QUIT",width=9)
    beval<-tkbutton(tf,text="EVALUATE",width=9)
    bnext<-tkbutton(tf,text="NEXT",width=9)
    bback<-tkbutton(tf,text="BACK",width=9)
    tkbind(top, "<<EvalRCode>>", eval.code)
    if( ! substring(version$os,1,6)=="darwin"  ) {
      tkbind(top, "<<Next>>", show.next.number)
      tkbind(top, "<<Back>>", show.back.number)
    }
    lno  <-tkentry(tf,textvariable=secno,width=9)
    linfo<-tklabel(tf,text="chunk number:")
    tkpack(linfo,lno,beval,bnext,bback,bexit,side="left")
    tkconfigure(bexit,command=exit.function)
    tkconfigure(bnext,command=show.next.number)
    tkconfigure(bback,command=show.back.number)
    tkconfigure(beval,command=eval.code)
    # tkbind(lno,"<Return>",show.number)
    tkbind(lno,"<KeyRelease>",show.number)
    tclvalue(secno)<-as.character(no)
    show.number()
    ### tkwait.window(top)
  }

  ft<-deparse(ft)[-(1:2)]; ft<-sub("WoRkNaMe",workname,ft)
  # modify function
  txt<-c("CCPlayer<-",fh,"    chunks<-readLines(file)",ft,"CCPlayer()")
  # set file argument
  txt[2] <- sub("file",paste("file = '",workname,"'",sep=""),txt[2]) 
  # set no argument as start chunk number
  idx <- grep("^    no .. 0",txt); txt[idx] <- sub("0","no #set start chunk",txt[idx])
  # modify title
  idx <- grep("tkwm.title",txt); txt[idx] <- sub("Peter Wolf.","",txt[idx])
  txt[idx] <- sub("demo of file","Code-Chunk-Player, file:",txt[idx])
  # modify exit message
  idx <- grep("Demo function stopped",txt); 
  txt[idx] <- sub("Demo function stopped without","Code-Chunk-Player stopped without",txt[idx])
  # print viewer function for debugging:    print(txt)
  # start of Code-Chunk-Player
  try(eval(parse(text=txt)))
  invisible(NULL)
}

playground<-function(playground.envir=NULL,code=NULL){
  # (tcltk) # 140306
  pg<-tktoplevel(); tkwm.geometry(pg,"+100+100")
  tkwm.title(pg, "playground for testing R code (error messages appear in the R Console)")
  pgtext<-tktext(pg,height=19,background="#f7fffF",
                 font="-Adobe-courier-Medium-R-Normal--18-180-*")
  ## ok: tkbind(pgtext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
if(substring(version$os,1,6)=="darwin"  ){
  mac.paste<-function(...){  ## Ctrl-V  # mac-PASTE
     try({.Tcl("clipboard append hello"); .Tcl("clipboard clear")
          news<-base::scan(file=pipe("pbpaste","r"),what="",sep="\n",blank.lines.skip=FALSE)
          tkinsert(pgtext,"insert",paste(news,collapse="\n"))})
          tksee(pgtext,"insert - 7 lines");  tksee(pgtext,"insert + 7 lines") #090706 
  } 
  tkbind(pgtext,"<Control_L><v>",mac.paste)
  mac.copy<-function(...){  ## Ctrl-C  # mac-COPY
     news<-""
     try(news<-tclvalue(.Tcl("if {[catch {clipboard get}]} {set aa empty} {set aa full}")))
     if(news=="empty") return()
     try({news<-tclvalue(.Tcl("set aaa [selection get -selection CLIPBOARD]"))
          tmp.file.name <- tempfile("rt-tmp")
          base::cat(news,file=tmp.file.name); system(paste("pbcopy < ",tmp.file.name))
          .Tcl("clipboard append hello"); .Tcl("clipboard clear")})
  }
  tkbind(pgtext,"<Control_L><c>",mac.copy)
  tkevent.add("<<extract>>",  "<Control_L><c><KeyRelease>")   # mac-extract
  tkbind(pgtext,"<<extract>>",mac.copy)
}else{
  tkevent.add("<<Paste>>",   "<Control_L><v>")
  tkbind(pgtext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
}

  pgbutfr<-tkframe(pg)
  tkpack(pgtext, side="top",fill="both",expand="y")
  tkpack(pgbutfr,side="top",fill="both")
  
  bexit<-tkbutton(pgbutfr,text="QUIT",width=9)
  beval<-tkbutton(pgbutfr,text="EVALUATE",width=9)
  tkpack(beval,side="left")
  tkpack(bexit,side="right")

  if(!is.null(code)) try(tkinsert(pgtext,"0.0",paste(code,collapse="\n")))

  eval.code<-function(){
    code<-tclvalue(tkget(pgtext,"0.0","end"))
    code.orig<-code<-unlist(strsplit(code,"\n"))
    if(length(code)==0){ cat("warning: no code found!\n"); return() }    
    if(!is.null(playground.envir)){  ## exists("revive.env") ##
      if(0 < grep("revive.sys", ls(envir=playground.envir))){
        revive.sys <- get("revive.sys",playground.envir)
        if(0<length(code)){
          tld <- paste("<","<*>",">=",sep="")
          rh <- c( get("relax.history",envir=revive.sys), list(c("@",tld,code)) )
          assign("relax.history",rh,envir=revive.sys)
        }
   ## stores code in revive.sys
      }
      result<-try(eval(parse(text=code),envir=playground.envir)) 
    }else{ 
      result<-try(eval(parse(text=code),envir=pos.to.env(1)))
    }
    if(class(result)=="try-error"){
      class(result)<-"character"; cat(result,"\n"); return()
    } else { 
      #idx <- which("relax.fns"==search()) # 121214
      #if( 0 < length(idx) ) print <- get("print",pos=idx)
      print(result) 
    }  
    NULL
  }
  tkbind(pg, "<<EvalRCode>>", eval.code)
  exit.function<-function(){
     tkdestroy(pg)
     return()
  }
  tkconfigure(bexit,command=exit.function)
  tkconfigure(beval,command=eval.code)
  invisible(NULL)
}

readline<-function(prompt=""){
  # for str, cat, ... it is necessary that revive.sys and revive.env exists
  revive.sys <- "noenv"
  if (!exists("revive.env")) revive.env <- "noenv"
  if (is.environment(revive.env))
      revive.sys <- get("revive.sys", envir=revive.env)
  # dummy function used in chunk to append info to output widget
  melde <- function(...) {"relax"} 
  # cat("revive.env"); print(revive.env); cat("revive.sys"); print(revive.sys)


  if( !((is.environment(revive.env) && 
         "toutwin" %in% ls( pattern="toutwin", envir=r.sys<-get("revive.sys",revive.env)) && "1"==tclvalue(tkwinfo("exists",get("toutwin",r.sys)))) ) ){  ### 111103 # 130408
    readline<-base::readline
    return(readline(prompt=prompt))
  }
  set.tclvalue <- get("set.tclvalue",envir=revive.sys) # 121214
  linfo <- get("linfo",revive.sys)
  linfo.tmp <- get("linfo.tmp",revive.sys)
  einfo.tmp <- get("einfo.tmp",revive.sys)
  TopW <- get("TopW",revive.sys)
  if(!exists("tworkwin"))
    tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

  frage <- if(prompt=="") "Input?" else prompt
  set.tclvalue("tvinfo","")
  tkconfigure(linfo.tmp,text=frage)
  # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
  tkpack("forget",linfo); Sys.sleep(0.01)
  tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
  tkfocus(einfo.tmp)
  tkselection.range(einfo.tmp,"0","end") ## 051219
  tkbind(TopW,"<Escape>",function(){
      tkbind(TopW,"<Return>","")
      tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
      # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
      tkpack(linfo,side="left",fill="x",expand="yes")


    }
  )


  tkbind(TopW,"<Return>", function(){
      tkbind(TopW,"<Return>","")
      tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
      # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
      tkpack(linfo,side="left",fill="x",expand="yes")


      input<-tclvalue("tvinfo")
      set.tclvalue("tvreadlinedone",1)
      mess<-paste("relax"); set.tclvalue("tvmess",mess)
    } # end of function
  )
  tkbind(TopW,"<Escape>",function(){ # Esc should not have an effect
      "relax" 
    }
  )
  set.tclvalue("tvreadlinedone",0); tkwait.variable("tvreadlinedone")
  out<-tclvalue("tvinfo")
  if(0<nchar(out)){
    news<-if(prompt!="") prompt else "readline Input:"
    news<-paste("\n",news,"\n", out,sep="")
    if(!exists("toutwin"))
      toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
    pos.to.insert<-"end"
    ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
    news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
    try(tkinsert(toutwin,pos.to.insert,news))
    tksee(toutwin,"end - 0 lines")
    melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

    tcl("update","idletasks") # 090710
    NULL
  } 
  tkbind(TopW,"<Return>","")
  tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
  # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
  tkpack(linfo,side="left",fill="x",expand="yes")


  return(out)
}
## assign("readline",readline,pos=#<stelle Nummer von [[relax]] im Suchpfad fest>#)
## assign("readline",readline, pos=pos.of.relax.fns)
## melde("readline saved",3)

menu<-function(choices, graphics=FALSE, title=""){
  # for str, cat, ... it is necessary that revive.sys and revive.env exists
  revive.sys <- "noenv"
  if (!exists("revive.env")) revive.env <- "noenv"
  if (is.environment(revive.env))
      revive.sys <- get("revive.sys", envir=revive.env)
  # dummy function used in chunk to append info to output widget
  melde <- function(...) {"relax"} 
  # cat("revive.env"); print(revive.env); cat("revive.sys"); print(revive.sys)


  if( !((is.environment(revive.env) && 
         "toutwin" %in% ls( pattern="toutwin", envir=r.sys<-get("revive.sys",revive.env)) && "1"==tclvalue(tkwinfo("exists",get("toutwin",r.sys)))) ) ){ # 130408
    menu<-utils::menu
    return(menu(choices=choices,graphics=graphics,title=title))
  }else{
  set.tclvalue <- get("set.tclvalue",envir=revive.sys) # 121214
  TopN<-tktoplevel()
  if(title=="") title<-"select item, press Return (or exit by Esc)"
  tkwm.geometry(TopN,"+0+15"); tkwm.title(TopN, title)
  if(missing(choices)||length(choices)==0) choices<-"relax"
  choices<-paste(seq(choices),": ",choices, sep="")
  nc<-length(choices<- c(choices, "EXIT"))
  scr <- tkscrollbar(TopN, command=function(...)tkyview(tl,...))
  tl<-tklistbox(TopN,height=min(30,length(choices)),width=60,
                selectmode="single",yscrollcommand=
                  function(...)tkset(scr,...),background="white")
  for(ch in choices) tkinsert(tl,"end",ch)
  tkpack(tl,side="left",expand="yes",fill="y")
  tkpack(scr,side="left",expand="yes",fill="y")
  set.tclvalue("choice","0"); tkbind(TopN,"<Escape>",function(){tkdestroy(TopN)})
  tkbind(TopN,"<Return>",function(){
     choice<-as.numeric(tkcurselection(tl))+1
     set.tclvalue("choice",choice); tkdestroy(TopN)
  })
  tkwait.window(TopN)
  choice<-tclvalue("choice")
  choice<- if(choice==length(choices)) 0 else as.numeric(choice)
  news<-paste("",paste(choices,collapse="\n"),"\nSelection: ",choice, "",sep="")
  if(!exists("toutwin"))
    toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
  pos.to.insert<-"end"
  ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
  news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
  try(tkinsert(toutwin,pos.to.insert,news))
  tksee(toutwin,"end - 0 lines")
  melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

  tcl("update","idletasks") # 090710
  NULL
  tkfocus(get("tworkwin",envir=get("revive.sys",envir=revive.env))) # wichtig!!
  return(choice)
 }
}
## assign("menu",menu,pos=#<stelle Nummer von [[relax]] im Suchpfad fest>#)
## assign("menu",menu, pos=pos.of.relax.fns)
## melde("menu saved",3)

scan<-function(file = "", what = double(), nmax = -1L, n = -1L, sep = "",
             quote = if (identical(sep, "\n")) "" else "'\"",
             dec = ".", skip = 0L, nlines = 0L, na.strings = "NA",
             flush = FALSE, fill = FALSE, strip.white = FALSE,
             quiet = FALSE, blank.lines.skip = TRUE,
             multi.line = TRUE, comment.char = "", allowEscapes = FALSE, 
             fileEncoding = "", encoding = "unknown", text = NA){
    # for str, cat, ... it is necessary that revive.sys and revive.env exists
    revive.sys <- "noenv"
    if (!exists("revive.env")) revive.env <- "noenv"
    if (is.environment(revive.env))
        revive.sys <- get("revive.sys", envir=revive.env)
    # dummy function used in chunk to append info to output widget
    melde <- function(...) {"relax"} 
    # cat("revive.env"); print(revive.env); cat("revive.sys"); print(revive.sys)


    scan<-get("scan",pos="package:base")

    if(file!="" || !is.na(text) || !((is.environment(revive.env) && 
                                      "toutwin" %in% ls( pattern="toutwin", envir=r.sys<-get("revive.sys",revive.env)) && "1"==tclvalue(tkwinfo("exists",get("toutwin",r.sys)))) ) ){ # 130408
      scan.args <- list(file=file,what=what)
      if(!missing(nmax)) scan.args <- c(scan.args, list(nmax=nmax))
      if(!missing(n))    scan.args <- c(scan.args, list(n=n))
      scan.args <- c(scan.args, 
                     list(sep=sep,quote=quote,dec=dec,skip=skip,                  
                        nlines=nlines,na.strings=na.strings,
                        flush=flush,fill=fill,strip.white=strip.white,quiet=quiet,
                        blank.lines.skip=blank.lines.skip,
                        multi.line = multi.line, comment.char = comment.char,
                        allowEscapes = allowEscapes))
      if(paste(version$major,version$minor) > "212"){
        scan.args <- c(scan.args, list(fileEncoding = fileEncoding))
      }
      scan.args <- c(scan.args, list(encoding = encoding))
      if(paste(version$major,version$minor) > "214"){
        scan.args <- c(scan.args, list(text=text))
      }
      worktext <- do.call(scan, scan.args)

    } else {
      set.tclvalue <- get("set.tclvalue",envir=revive.sys) # 121214
      outfont.sys <- get("outfont.sys",envir=revive.sys) # 121214
      typ<-if(is.numeric(what)) "Zahlen-Eingabe" else "Text-Eingabe"
      .newl<-tktoplevel(); tkwm.title(.newl, typ); tkwm.geometry(.newl,"+0+15")
      revive.sys<-get("revive.sys",envir=revive.env)
      assign(".newl",.newl,envir=revive.sys)
      if(TRUE){ ## exists("running.function") && running.function==#<relax>#){ ### 111103
        tkpack(minfo<-tkmessage(.newl,width="1000",justify="left",relief="raised"))
        if(!exists("toutwin"))
          toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
        news<-strsplit(tclvalue(tkget(toutwin,"0.0","end")),"\n")[[1]]
        if(length(news)>10) news<-rev(rev(news)[1:10])
        # kaum relevant
        if(length(h<-grep("^output", news))>0) news<-news[-h]
        if(length(h<-grep("^@", news))>0)      news<-news[-h]
        if(length(h<-grep("verbatim",news))>0) news<-news[-h]
        if(length(news)>0){
          news<-rev(rev(news)[1:min(20,length(news))])
          news<-paste(news,collapse="\n")
          tkconfigure(minfo,text=news)
        }

      }else{
        tkpack(minfo<-tkmessage(.newl,width="1000",justify="left",relief="raised"),
               fill="both",expand="true")
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-worktext[worktext!=""]
        line<-rev(grep("^<<(.*)>>=",worktext))[1]
        if(!is.na(line)&length(line)>0) worktext<-worktext[-(1:line)]
        line<-grep("^@",worktext)[1]
        if(!is.na(line)&length(line)>0) worktext<-worktext[-(1:line)]
        if(length(h<-grep("^output", worktext))>0) worktext<-worktext[-h]
        if(length(h<-grep("^@", worktext))>0)      worktext<-worktext[-h]
        if(length(h<-grep("verbatim",worktext))>0) worktext<-worktext[-h]
        if(length(worktext)>0){
          worktext<-rev(rev(worktext)[1:min(20,length(worktext))])
          worktext<-paste(worktext,collapse="\n")
          tkconfigure(minfo,text=worktext)
        }

      }

      if(!missing(n)&&n==1){
        sp<-"10"; zei<-"1"
        tkbind(.newl,"<Return>",function() set.tclvalue("tvscandone",2))
      }else{sp<-"50"; zei<-"3"}
      tkbind(.newl,"<Return><Return>",function()set.tclvalue("tvscandone",2))
      tscan<-tktext(.newl,width=sp,height=zei)
      twin<-tscan; f.sonderzeichen<-function(zeichen){
                     function(){
                       tkinsert(twin,"insert",zeichen)
                       if((.Platform$OS.type=="windows"))tkdelete(twin,"insert-1chars")
                       if(substring(version$os,1,6)=="darwin" )tkdelete(twin,"insert-1chars")
                     }
                   }
                   f.umlaut<-function(zeichen){
                     function(){
                       return()
                       #char337<-eval(parse(text='"\\337"'))
                       #if(zeichen==char337 & tclvalue(tkget(twin,"insert-1chars","insert"))=="\\") return()
                       #tkinsert(twin,"insert",zeichen); tkdelete(twin,"insert-2chars")
                     }
                   }
                   tkbind(twin,"<<LKeckig>>", f.sonderzeichen("["))
                   tkbind(twin,"<<RKeckig>>", f.sonderzeichen("]"))
                   tkbind(twin,"<<Tilde>>",   f.sonderzeichen("~"))
                   tkbind(twin,"<<LKgeschw>>",f.sonderzeichen("{"))
                   tkbind(twin,"<<RKgeschw>>",f.sonderzeichen("}"))
                   tkbind(twin,"<<Klammera>>",f.sonderzeichen("@"))
                   tkbind(twin,"<<Pipe>>",    f.sonderzeichen("|"))
                   tkbind(twin,"<<Backsl>>",  f.sonderzeichen("\\"))
                   renewhighlighting<-function(){
                     tworkwin<-get("tworkwin",envir=revive.sys)
                     melde("ak texthervor",1)
                     tcl("markclear",tworkwin)
                     tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
                     tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
                     tcl("marklinetypes",tworkwin)
                     melde("ak texthervor",2)

                   }
                   # tkbind(twin,"<<Klammeraffe>>",renewhighlighting)
                   tkbind(twin,"<Return>",renewhighlighting)

      bexit<-tkbutton(.newl, text="Exit: end of input")
      tkpack(tscan, bexit); tkfocus(tscan)
      tkconfigure(bexit,command=function()set.tclvalue("tvscandone",2))
      set.tclvalue("tvscandone",0); tkwait.variable("tvscandone")
      worktext<-tclvalue(tkget(tscan,"0.0","end"))
      if(!exists("tworkwin"))
        tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

      tkfocus(tworkwin); tkdestroy(.newl)
      news<-paste("\nscan-Eingabe:\n",worktext,sep="")
      if(!exists("toutwin"))
        toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
      pos.to.insert<-"end"
      ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
      news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
      try(tkinsert(toutwin,pos.to.insert,news))
      tksee(toutwin,"end - 0 lines")
      melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      tcl("update","idletasks") # 090710
      NULL
           worktext<-strsplit(worktext,"\n")[[1]]
      worktext<-strsplit(paste(worktext,collapse=" ")," ")[[1]]
      worktext<-worktext[worktext!=""]
      if(typ=="Zahlen-Eingabe"){
        try.res<-try(as.numeric(worktext))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok) worktext<-try.res else NULL
      }


    }
    worktext
}
## assign("scan",scan,pos=#<stelle Nummer von [[relax]] im Suchpfad fest>#)
## assign("scan",scan, pos=pos.of.relax.fns)
## melde("scan saved",3)

print<-function(x, ...){ # 050614
  # for str, cat, ... it is necessary that revive.sys and revive.env exists
  revive.sys <- "noenv"
  if (!exists("revive.env")) revive.env <- "noenv"
  if (is.environment(revive.env))
      revive.sys <- get("revive.sys", envir=revive.env)
  # dummy function used in chunk to append info to output widget
  melde <- function(...) {"relax"} 
  # cat("revive.env"); print(revive.env); cat("revive.sys"); print(revive.sys)


  if( (is.environment(revive.env) && 
       "toutwin" %in% ls( pattern="toutwin", envir=r.sys<-get("revive.sys",revive.env)) && "1"==tclvalue(tkwinfo("exists",get("toutwin",r.sys))))  && !is.null(x) ){ # 130408
    sink(get("tmp.file.name",envir=revive.sys)
); base::print(x, ...); sink()
    news<-paste(# "", # 111123
                paste(scan(file=get("tmp.file.name",envir=revive.sys)
,what="",sep="\n"),collapse="\n"),
                "",sep="\n" )
    if(!exists("toutwin"))
      toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
    pos.to.insert<-"end"
    ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
    news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
    try(tkinsert(toutwin,pos.to.insert,news))
    tksee(toutwin,"end - 0 lines")
    melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

    tcl("update","idletasks") # 090710
    NULL
    if(all(class(x) != "help_files_with_topic")) base::print(x, ...)
  }  else base::print(x, ...)
  invisible(NULL)
}
## assign("print",print,pos=#<stelle Nummer von [[relax]] im Suchpfad fest>#)
## assign("print",print, pos=pos.of.relax.fns)
## melde("print saved",3)

cat<-function(...,file="",sep=" ",fill=FALSE,labels=NULL,append=FALSE){
  # for str, cat, ... it is necessary that revive.sys and revive.env exists
  revive.sys <- "noenv"
  if (!exists("revive.env")) revive.env <- "noenv"
  if (is.environment(revive.env))
      revive.sys <- get("revive.sys", envir=revive.env)
  # dummy function used in chunk to append info to output widget
  melde <- function(...) {"relax"} 
  # cat("revive.env"); print(revive.env); cat("revive.sys"); print(revive.sys)


  if(file=="" &&  (is.environment(revive.env) && 
                   "toutwin" %in% ls( pattern="toutwin", envir=r.sys<-get("revive.sys",revive.env)) && "1"==tclvalue(tkwinfo("exists",get("toutwin",r.sys))))  ){ # 130408
     base::cat(...,file=get("tmp.file.name",envir=revive.sys)
,sep=sep,fill=fill,labels=labels,append=append)
     news<-paste( # "\n", # 111123
                 paste(scan(file=get("tmp.file.name",envir=revive.sys)
,what="",sep="\n"),collapse="\n"), "", # 111123
                 sep="\n")
     if(!exists("toutwin"))
       toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
     pos.to.insert<-"end"
     ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
     news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
     try(tkinsert(toutwin,pos.to.insert,news))
     tksee(toutwin,"end - 0 lines")
     melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

     tcl("update","idletasks") # 090710
     NULL
     invisible(NULL)
  }
  base::cat(...,file=file,sep=sep,fill=fill,labels=labels,append=append)
}
# assign("cat",cat,pos=#<stelle Nummer von [[relax]] im Suchpfad fest>#)
# assign("cat",cat, pos=pos.of.relax.fns)
# melde("cat saved",3)

str<-function(object,...){
  # for str, cat, ... it is necessary that revive.sys and revive.env exists
  revive.sys <- "noenv"
  if (!exists("revive.env")) revive.env <- "noenv"
  if (is.environment(revive.env))
      revive.sys <- get("revive.sys", envir=revive.env)
  # dummy function used in chunk to append info to output widget
  melde <- function(...) {"relax"} 
  # cat("revive.env"); print(revive.env); cat("revive.sys"); print(revive.sys)


  if( (is.environment(revive.env) && 
       "toutwin" %in% ls( pattern="toutwin", envir=r.sys<-get("revive.sys",revive.env)) && "1"==tclvalue(tkwinfo("exists",get("toutwin",r.sys))))  ){ # 130408
     fname <- get("tmp.file.name",envir=revive.sys)
; fname <- gsub("\\\\","/",fname)
     base::cat(file=fname,"")
     if(is.data.frame(object)) 
        base::cat(file=fname,
                  paste("'data.frame':  ", 
                        paste(dim(object)[1],"obs. of",
                              dim(object)[2],"variables:\n")))
     a <- deparse(getS3method("str", "default"))
     a <- gsub("cat[(]", paste("base::cat(file=\"", fname, 
               "\",append=TRUE,", sep = ""), a) # ))
     a <- sub(" str[(]", " mystr(", a) # ))
     mystr <- eval(parse(text = a))
     mystr(object,...)  # 111123
     news <- scan(file = fname,
                  what = "", sep = "\n")
     news <- sub("chr .data.frame.", "", news)
     if(0<length(ind<-grep("^ -",news))) news<-news[-ind]
     ## in data frames additional "List of ?" appear and should be removed
     if(is.data.frame(object) && 0<length(ind<-grep("^List of",news))) news<-news[-ind]
     news<-paste("\n", # date(),
                 paste(news,collapse="\n"),"",sep="\n")
     if(!exists("toutwin"))
       toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
     pos.to.insert<-"end"
     ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
     news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
     try(tkinsert(toutwin,pos.to.insert,news))
     tksee(toutwin,"end - 0 lines")
     melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

     tcl("update","idletasks") # 090710
     NULL
  } else { utils::str(object,...) }
  invisible(NULL)
}
# assign("str",str,pos=#<stelle Nummer von [[relax]] im Suchpfad fest>#)
# assign("str",str, pos=pos.of.relax.fns)
# melde("str saved",3)

step.dep<-deparse(stats::step)
# 120828
step.dep<-gsub("cat[(]",   "relax::cat(",   step.dep) # ))  
step.dep<-gsub("print[(]", "relax::print(", step.dep) # ))  
## step.dep<-gsub("cat[(]", 'get("cat",  pos=which(search()=="relax.fns"))(',step.dep) # ))  
## step.dep<-gsub("print[(]", 'get("print",pos=which(search()=="relax.fns"))(',step.dep) # ))  

###step.dep<-gsub("cat[(]", "relax::cat(",step.dep) # ))  
###step.dep<-gsub("print[(]", "relax::print(",step.dep) # ))  
mystep <- eval(parse(text = step.dep))
step<- mystep
## assign("step",mystep,pos=#<stelle Nummer von [[relax]] im Suchpfad fest>#)
## assign("step",step, pos=pos.of.relax.fns)
## melde("step saved",3)


#.First.lib<-function(lib, pkg) {
#      options(warn=1)
#      options(warning.expression={cat("WARN: Warning-Info see Console window")})
#      cat(".First.lib erledigt!\n")
#}
relax<-function(file.name,no.plots=FALSE,cmds="",but.Wizardry="all"){
 # Copyright (C) 2005--2008 Hans Peter Wolf
  options(warn=1)
 #      options(warning.expression={cat("WARNING: Warning-Info see Console window")})
  relax.path<-path.package("relax"
)  # 130325 .path.package defunct
  relax.pos<-grep("^package:relax$",search())[1]

  ## running.function<- #<relax>#  assign("running.function",#<relax>#, pos=relax.pos) ##111103
  revive.env<-"1"; rm(revive.env)
  if(!exists("revive.env")){
    revive.env<<-rev.env<-new.env(parent=.GlobalEnv) ### 111103
  } else {
    revive.sys<-get("revive.sys",envir=revive.env)
    rm(list=ls(envir=revive.sys),  envir=revive.sys)
  }
  revive.sys<-environment()
  assign("revive.sys", revive.sys, envir=revive.env)
  # if(0<length(grep("relax.fns",search()))) detach("relax.fns") # 121214 ???
  # attach(what=NULL,name="relax.fns")
  # pos.of.relax.fns <- grep("relax.fns",search())
  # if(!exists("path.package") & exists(".path.package")) path.package <- .path.package # 130516
  # if(!exists("find.package") & exists(".find.package")) find.package <- .find.package
  if(version$major <= 2 && version$minor < 13){
    if(!exists("path.package")) eval(parse(text="path.package <- .path.package"),envir=revive.env)
    if(!exists("find.package")) eval(parse(text="find.package <- .find.package"),envir=revive.env)
  }  


  myhead.menu<-function(item="Test",code=function()cat("Menu-Test"),
                        title="Menue",rm.menu=FALSE,menu.no=1){
    set.tclvalue<-function(name,value)tclvalue(name)<-as.character(value) # 121217
    out.msg<-NULL
    menu.widget.name<-paste("mrtrevive",   menu.no,sep="")
    menu.item.name  <-paste("mmyhead.menu",menu.no,sep="")
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    if( !( exists(menu.widget.name,envir=revive.sys)
 )){
      if( rm.menu==TRUE && menu.no!=0){
        tkmessageBox(title="Warnung", icon="warning",
                     message="Achtung:  kein Kopfmenue vorhanden")
        return("Error")
      }
      fhead<-get("fhead",envir=revive.sys)
      mrtrevive<-tkmenubutton(fhead,text=title,
                              font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"
,
                              relief="flat", width=10
)
      assign(menu.widget.name,mrtrevive,envir=revive.sys)
      tkbind(mrtrevive,"<Enter>",function() set.tclvalue("tvmess","R cmd menue"))
      tkbind(mrtrevive,"<Leave>",function(){ set.tclvalue("tvmess","relax") })
      tkpack(mrtrevive,side="left")
      mmyhead.menu<-tkmenu(mrtrevive)
      assign(menu.item.name,mmyhead.menu,envir=revive.sys)
      tkconfigure(mrtrevive, menu=get(menu.item.name,envir=revive.sys))
      out.msg<-c(out.msg, paste("Menue im Kopf eingerichtet!"))

    }
    if( rm.menu==TRUE ){
      mrtrevive<-get(menu.widget.name,envir=revive.sys)
      tkpack("forget",mrtrevive)
      tkdestroy(mrtrevive)
      remove(list=menu.widget.name,envir=revive.sys)
      out.msg<-c(out.msg, paste("Menue im Kopf entfernt!"))

    }else{
      item<-as.character(item[[1]])
      mmyhead.menu<-get(menu.item.name,envir=revive.sys)
      revive.env <- get("revive.env",envir=revive.sys)
      if(is.function(code))  code<-deparse(body(code))
      item.code<-function()eval(parse(text=code),envir=revive.env)
      tkadd(mmyhead.menu, "command", label=item, command=item.code,
              font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"

      )
      out.msg<-c(out.msg, paste("Item",item,"im Kopfmenue eingerichtet!"))

    }
    return(out.msg)
  }
  # assign("myhead.menu",myhead.menu, pos=relax.pos) ### 111103
  # assign("myhead.menu",myhead.menu, pos=pos.of.relax.fns)

  DEBUG<-"1"; rm(DEBUG) 
  melde<-function(report.msg,typ=0,...){
    cat<-get("cat",pos="package:base")
    # if(exists("last.warning"))get("print",pos="package:base")(last.warning)
    if(typ==0) cat(report.msg,...,"\n")
    if(exists("DEBUG") && TRUE==DEBUG){
      if(typ==1 | typ==2) cat(c(a="Start:",z="Ende:")[typ], report.msg,"\n")
      if(typ==3){cat(report.msg); cat(...); cat("\n")}
    }
    if(typ=="cmd.msg"){
      if((.Platform$OS.type=="windows")) flush.console()
      if(exists("report.msg.to.file.flag")) cat(report.msg, "\n", file="report.msg", append=TRUE)
    }
    invisible()
  }
  # assign("melde",melde, pos=relax.pos) ## 111103
  #assign("melde",melde, pos=pos.of.relax.fns) ## 121214 ???

  #############################################
  # configuration file of relax
  #############################################
  #
  ## general parameters of relax
  # language for messages:
  language<-"english"
  # language<-"german"
  # maximal characters of output to be shown
  maxol.sys<-4000
  # size of fonts # set of values: 1,2, ..., 7
  initial.font.size<-4
  # height of included postscript plot in LaTeX document:
  psheight.sys<-"10cm"
  # height parameter for postscript device
  psdesignheight.sys<-6
  # width parameter for postscript device
  psdesignwidth.sys<-6
  # rotation of PS graphics
  pshorizontal.sys<-FALSE
  # size parameter for jpeg device in inch
  jpgdesignsize.sys<-4
  # resolution for ppm-plot within report widget
  ppmresolution.sys<-15
  # height of relax window ( e.g.: "500" or "1200" ) 
  relaxwindow.height.sys<-"700"
  # width of relax window ( e.g.: "500" or "1200" ) 
  relaxwindow.width.sys<-"700"
  # use of "tab"-key for name completion
  name.complete.sys<-TRUE
  # replace german umlaute
  replace.umlaute.sys<-TRUE # FALSE
  #############################################
  ## special settings for windows
  # how to call LaTeX
  latex.command.windows<-"echo q | latex"
  # view command for viewing .dvi files
  view.command.windows<-"yap"
  # how to call dvipdf
  dvipdf.command.windows<-"dvipdfm"
  # text editor
  text.editor.windows<-"notepad"
  # browser
  browser.windows<-" "
  # browser.windows<-"c:/Programme/\"Mozilla Firefox\"/firefox.exe "
  pdfview.windows<-"acroread"
  # ghostscript-program e.g.: 
  #   ghostscript<-"c:\\gs\\gs7.04\\bin\\gswin32c"
  #   ghostscript<-"C:\\Programme\\gs\\gs8.71\\bin\\gswin32c
  ghostscript<-" "
  ## defunc: Path to Tcl/Tk-img-package: # imgpath.sys <- "C:/Tcl/lib"
  #############################################
  ## special settings for linux
  # how to call LaTeX
  latex.command.linux<-"echo q | latex"
  # view command for viewing .dvi files
  view.command.linux<-"xdvi"
  # how to call dvipdf
  dvipdf.command.linux<-"dvipdf"
  # text editor ?????????
  text.editor.linux<-"kwrite"
  # browser
  browser.linux<-"konqueror"
  # pdf viewer
  pdfview.linux<-"acroread"
  ## defunc: path.tcltk.package.img<-"/usr/local/lib"
  #############################################
  ## special settings for mac
  # how to call LaTeX
  latex.command.mac<-"echo q | pdflatex" 
  # view command for viewing .dvi files
  view.command.mac<-"open " 
  # how to call dvipdf
  dvipdf.command.mac<-"dvipdf" 
  # text editor
  text.editor.mac<-"nano"
  # browser
  browser.mac<-"safari"

  try({
    settings<-scan(file=file.path(relax.path,"config/settings.relax"),what="",sep="\n")
    for(i in seq(settings)) eval(parse(text=settings[i]))
  })
  editor.sys<-text.editor.mac
  if((.Platform$OS.type=="windows")) {
    editor.sys <- text.editor.windows
  }
  if(substring(version$os,1,5)=="linux") {
    editor.sys <- text.editor.linux
  }
  browser.sys<-""
  if((.Platform$OS.type=="windows")) {
    browser.sys<-"start "
    if(exists("browser.windows")&&nchar(browser.windows)>0&&0<length(grep("[a-zA-Z]",browser.windows))) 
      browser.sys <- browser.windows
  }
  if(substring(version$os,1,5)=="linux") {
    if(exists("browser.linux")) browser.sys <- browser.linux
  }
  if(substring(version$os,1,6)=="darwin" ){
    browser.sys<-"open "
    latex.command.linux<-latex.command.mac
    view.command.linux<-view.command.mac
    dvipdf.command.linux<-dvipdf.command.mac
    text.editor.linux<-text.editor.mac
  }
  if(!exists("initial.font.size")) initial.font.size<-4
  initial.font.size<-initial.font.size[1]
  if(is.na(initial.font.size) || all(initial.font.size!=(1:7))) initial.font.size<-4
  if(!exists("relaxwindow.width.sys")) relaxwindow.width.sys<-"700"
  if(!exists("relaxwindow.height.sys")) relaxwindow.height.sys<-"700"
  if(!exists("name.complete.sys")) name.complete.sys<-TRUE
  if((.Platform$OS.type=="windows") && nchar(ghostscript)<=1){
    if(Sys.getenv("R_GSCMD")!="") ghostscript <- Sys.getenv("R_GSCMD")
    else { path.sep <- ";"
      # find pathes where a search will take place and generate all parts of the pathes
      s.path<-unlist(strsplit(Sys.getenv("PATH"),path.sep))
      s.path<-unique(unlist(lapply(s.path,function(x){
        dirs<-unlist(strsplit(x,"\\\\"))
        sapply(1:length(dirs), 
          function(x) paste(dirs[1:x],collapse=.Platform$file.sep))[-1]
      })))
      for(i in 1){
        # append "gs" to the pathes and look for existing pathes
        s.path<-file.path(s.path,"gs"); s.path<-s.path[file.exists(s.path)]
        if(0==length(s.path)) break
        # look for files (dirs) of the form gs7.81 in the generated pathes
        pattern <- "^gs[1-9]" # pattern <- "^gs[bj1-9]" # for testing
        idx <- sapply(s.path, function(x) 0<length(grep(pattern,list.files(x))))
        s.path <- s.path[idx]
        if(0==length(s.path)) break
        # extract the first path that matches the conditions and append gswin32c.exe
        dir <- lapply(s.path, function(x) (grep("^gs[bj1-9]",list.files(x),value=TRUE))[1])
        s.path <- file.path(s.path,dir,"bin","gswin32c.exe")[1]    
        if(!file.exists(s.path)){
          s.path<-NULL; cat("relax warning: gswin32c.exe not found")
        }
      }
      if(0<length(s.path)) {
        ghostscript <- s.path; Sys.setenv("R_GSCMD"=ghostscript)
      }
    }
  }

  if(substring(version$os,1,5)=="linux" && nchar(ghostscript)<=1){
    if(Sys.getenv("R_GSCMD")!="") ghostscript <- Sys.getenv("R_GSCMD")
    else { 
      gsexe <- Sys.getenv("R_GSCMD")
      if (is.null(gsexe) || !nzchar(gsexe)) {
        gsexe <- system("which gs",intern=TRUE)
      }
      if (!is.null(gsexe) && 0<length(gsexe)) {
        ghostscript <- gsexe 
      }
    }
  }

  ##lade das [[Img]]-Paket für [[Tcl/Tk]]>> 121212 
  set.tclvalue<-function(name,value)  tclvalue(name)<-as.character(value)
  # assign("set.tclvalue",set.tclvalue, pos=relax.pos) ### 111103
  # assign("set.tclvalue",set.tclvalue, pos=pos.of.relax.fns) ## 121214  ok ???
    
  .Tcl("set XYZ [encoding system]")
  UTF<-is.UTF<- 0<length(grep("utf",tclvalue("XYZ")))

  # myscan<-get("scan",pos="package:base"); formals(myscan)$comment.char<-""
  myscan<-function(file,what,sep="\n", blank.lines.skip=FALSE){ readLines(file) } #2.1.0

  WinToTcl.read<-function(x){  # 130507 Umlautbehandlung
    try(if(is.UTF && replace.umlaute.sys){  
        #utf8.umlaute <- "[\xa4\xb6\xbc\x84\x96\x9c\x9f]" # utf-Umlaute ohne xc3
        #umlaute <- "[\xe4\xf6\xfc\xc4\xd6\xdc\xdf]" # latin1-Umlaute
        #umlaute <- iconv(umlaute,"latin1","")       # locale-Umlaute
        utf8.umlaute <- iconv("\x84\x94\x81\x8e\x99\x9a\xe1","cp850","utf8")
        utf8.umlaute <- paste(strsplit(utf8.umlaute,"")[[1]],collapse="|")
        uml.from.other <- c(
  #         length(grep(utf8.umlaute,useBytes=TRUE,x)),
  #         length(grep(utf8.umlaute,useBytes=TRUE,iconv(x,"latin1",""))),
  #         length(grep(utf8.umlaute,useBytes=TRUE,iconv(x,"macintosh",""))),
  #         length(grep(utf8.umlaute,useBytes=TRUE,iconv(x,"cp850",""))) 
           length(grep(utf8.umlaute,useBytes=!TRUE,iconv(x,"utf8",""))),
           length(grep(utf8.umlaute,useBytes=!TRUE,iconv(x,"latin1",""))),
           length(grep(utf8.umlaute,useBytes=!TRUE,iconv(x,"macintosh",""))),
           length(grep(utf8.umlaute,useBytes=!TRUE,iconv(x,"cp850",""))), 
           length(grep(utf8.umlaute,useBytes=!TRUE,iconv(x,"HPROMAN8","")))
        )
        melde(paste(utf8.umlaute,"uml.from",c("locale","latin1","mac","cp850",
                                              "HPROMAN8"),uml.from.other),3)
        idx.max <- which.max(uml.from.other)
        from <- c("","LATIN1", "macintosh", "cp850", "HPROMAN8")[idx.max]
        if( 1 < idx.max){
            melde(paste("source encoding maybe:",from),3)
            res<-tkmessageBox(
              message=if(language=="german") 
                        paste("Dokument-Encoding wahrscheinlich:", from, 
                              "Soll es zum lokalen umgewandelt werden?")
                      else 
                        paste("Encoding of document maybe:", from, 
                              "Do you want to change it to the locale one?"),
              title="Encoding -- German Umlaute", icon="warning", 
              type="yesnocancel", default="yes"
            )
            if("externalptr"==mode(res))  res<-tclvalue(res)
            if(res=="yes"){
              xx <- sapply(x, function(x) iconv(x,from,""))
              idx<-is.na(xx); xx[idx]<-x[idx]; x<-xx
              cat("========================================================================\n")
              cat("WARNING:",from,"encoding of document has been changed to locale encoding\n")
              cat("========================================================================\n")
              Sys.sleep(1)
            }
        }

    })
    try(if(!is.UTF&&replace.umlaute.sys){ 
        umlaute <- "[\xe4\xf6\xfc\xc4\xd6\xdc\xdf]" # latin1-Umlaute
        umlaute <- iconv(umlaute,"latin1","")       # locale-Umlaute
        uml.from.other <- c(
           length(grep(umlaute,useBytes=TRUE,x)),
           length(grep(umlaute,useBytes=TRUE,iconv(x,"utf8",""))),
           length(grep(umlaute,useBytes=TRUE,iconv(x,"macintosh",""))), 
           length(grep(umlaute,useBytes=TRUE,iconv(x,"cp850",""))),
           length(grep(umlaute,useBytes=TRUE,iconv(x,"HPROMAN8","")))
        )
        melde(paste(umlaute,"uml.from",c("locale","utf8","mac","cp850","HPROMAN8"),uml.from.other),3)
        idx.max <- which.max(uml.from.other)
        from <- c("","utf8", "macintosh", "cp850", "HPROMAN8")[idx.max]
        if( 1 < idx.max){
            melde(paste("source encoding maybe:",from),3)
            res<-tkmessageBox(
              message=if(language=="german") 
                        paste("Dokument-Encoding wahrscheinlich:", from, 
                              "Soll es zum lokalen umgewandelt werden?")
                      else 
                        paste("Encoding of document maybe:", from, 
                              "Do you want to change it to the locale one?"),
              title="Encoding -- German Umlaute", icon="warning", 
              type="yesnocancel", default="yes"
            )
            if("externalptr"==mode(res))  res<-tclvalue(res)
            if(res=="yes"){
              xx <- sapply(x, function(x) iconv(x,from,""))
              idx<-is.na(xx); xx[idx]<-x[idx]; x<-xx
              cat("========================================================================\n")
              cat("WARNING:",from,"encoding of document has been changed to locale encoding\n")
              cat("========================================================================\n")
              Sys.sleep(1)
            }
        }

    })
    return(x)
  }

  if(is.UTF)  TcltoWin.write<-function(x){ return(x) } else {
    # Latin1-Umwandlung von Umlauten
    pc<-eval(parse(text='"\\303"'))  # UTF-8-pre-char, old: 283
    uml.utf.8 <-eval(parse(text='"\\244\\266\\274\\204\\226\\234\\237"'))
    uml.latin1<-eval(parse(text='"\\344\\366\\374\\304\\326\\334\\337"'))
    TcltoWin.write<-function(x){
      if(replace.umlaute.sys){
        x<-chartr(uml.utf.8,uml.latin1,gsub(pc,"",x))
      }
      return(x)
    }
  }

## attach(relax.fns)
  tmp.file.name <- tempfile("rt-tmp")
  assign(tmp.file.name,"tmp.file.name",envir=revive.sys)

  tmp.sink.name <- tempfile("rt-sink")
  assign(tmp.sink.name,"tmp.sink.name",envir=revive.sys)
  ##definiere [[SaveAsHtml]]##           
  ##definiere Testknopf\-funktion##
  fEvalRCode<-function(){
    melde("fEvalRCode",1)
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    if(any(h<-(substring(worktext,1,1)==">"))){ # 091026
      worktext[h]<-sub("^.","@SpLiTtE?<<*>>=SpLiTtE?",worktext[h])
      worktext<-substring(unlist(strsplit(paste("",worktext),"SpLiTtE")),2)
      if(is.vector(line)) line<-line+2*sum(which(h)<=line)
    }

    code.start<-grep("^<<(.*)>>=",worktext)
    if(0==length(code.start)){cat("Warning: no code found!!!\n");return()}
    code.start<-code.start[code.start<=line]
    if(0<length(code.start)){
      if((code.start<-max(code.start)+1)>length(worktext)) return()
      code.end  <-c(grep("^@",worktext),1+length(worktext))
      code.end  <-min(code.end[code.end>code.start])-1
      code<-worktext[code.start:code.end]
      code<-code[code!=""]
      if(length(weg.ab<-grep("^<<(.*)>>=",code))>0) code<-code[-(weg.ab:length(code))]
      if(length(code)==0 || code[1]=="@") code<-" "
      melde("code:",3,code,"\n")

    }
    if(0<length(code)){
      melde("vor tangle - Code:\n",3,code)
      if(length(grep("<<(.*)>>",code))>0 || length(grep(">#",code))>0){
        code.a    <- grep("^<<(.*)>>=",worktext)
        code.z    <- grep("^@",worktext)
        code.z    <- unlist(sapply(code.a ,function(x,y) y[y>x][1], code.z))
        if(any(h<-is.na(code.z))) code.z<-code.z[!h]
        ################
        # code.n    <- length(worktext)
        # change    <- rep(0,code.n); change[c(code.a ,code.z)]<-1
        # code.ch   <- worktext[1==(cumsum(change)%%2)]
        ## 080311
        ind<-rep(0,length(worktext)); copy<-0
        for(i in seq(ind)){
          if(i %in% code.z) copy<-0
          if(i %in% code.a) copy<-1
          ind[i]<-copy
        }
        code.ch   <- worktext[ind==1]
        ## cat("input-text, Anfaenge, Code.chunks"); print(worktext); print(code.a); print(code.ch)

        code.n    <- length(code.ch)

        code.ch<-gsub("@>>","DoSpCloseKl-esc",gsub("@<<","DoSpOpenKl-esc",code.ch))

        code.ch<-gsub("(.*)<<(.*)>>=(.*)","cOdEdEf\\2",code.ch)
        repeat{
          if(0==length(cand<-grep("<<(.*)>>",code.ch))) break
          code.ch<-unlist(strsplit(gsub("(.*)<<(.*)>>(.*)",
                     "\\1bReAkuSeChUnK\\2bReAk\\3",code.ch),"bReAk"))
        }
        code.ch<-code.ch[code.ch!=""]
        code.n<-length(code.ch)
        melde("code.ch:",3,code.ch,"\n")

        line.typ  <-rep("C",code.n)
        code.a    <-grep("cOdEdEf",code.ch)
        code.ch[code.a]<-substring(code.ch[code.a],8)
        line.typ[code.a]<-"D"
        code.use    <-grep("uSeChUnK",code.ch)
        code.ch[code.use]<-substring(code.ch[code.use],9)
        line.typ[code.use]<-"U"
        code.ext  <-grep("#<file",code.ch)
        line.typ[code.ext]<-"E"
        melde("code.ch:",3,code.ch,"\n")

        code.out<-"##act:##"

        def.names<-code.ch[code.a]
        use.names<- if(length(code.use)>0) code.ch[code.use] else NULL
        code.z<-c(if(length(code.a)>1) code.a[-1]-1, code.n)
        code.ch<-paste(line.typ,code.ch,sep="")
        melde("code.ch:",3,code.ch,"\n")

        melde("vor expand - Code:\n",3,code)
        melde("bearbeite aktuellen Chunk\n",3)
        ###<bestimme Cursorzeile [[line]] von [[tworkwin]]> not a good idea###
        ch.no<-length(grep("^<<(.*)>>=",worktext[1:line]))

        rows      <-c((code.a[ch.no]+1),code.z[ch.no])
        if(all(!is.na(rows))&&rows[1]<=rows[2]){
          rows<-rows[1]:rows[2]
          code.stack<-code.ch[rows]
          max.depth.refinements<-500; i<-1
          repeat{
             if((i<-i+1)>max.depth.refinements){ 
                 cat("ERROR: maximal number of expandations (",max.depth.refinements,
                     ") exceeded\n --- perhaps a unintended recursion ???")
                 return()
             }
             if(0==length(code.stack))break
             typ<-substring(code.stack[1],1,1)
             if("C"==typ||"E"==typ){
               n.lines<-sum(cumprod("C"==substring(code.stack,1,1)))
               code.out<-c(code.out, substring(code.stack[1:n.lines],2))
               code.stack<-code.stack[-(1:n.lines)]
             }
             if(length(code.stack)>0 && "U"==substring(code.stack[1],1,1)){
               if(any(found<-def.names==substring(code.stack[1],2))){
                 found<-seq(along=def.names)[found]; rows<-NULL
                 for(no in found){
                   if((code.a[no]+1)<=code.z[no]) rows<-c(rows,(code.a[no]+1):code.z[no])
                 }
                 code.stack<-c(code.ch[rows],code.stack[-1])
                 melde("found",0,found)
               }else{code.stack<-code.stack[-1]}
             }

          }
        }
        if(length(code.ext)>0){
          code.out<-code.out[code.out!=""]
          code.ext<-rev(grep(">#",code.out))
          found<-TRUE
          repeat{
            if(length(code.ext)==0) break

            if(!found){
              code.out[code.ext[1]]<-paste("# ??",code.out[code.ext[1]])
              cat("ERROR: External Chunk",code.out[code.ext[1]],"not found!!!\n")
              code.ext<-code.ext[-1]
            }

            found<-TRUE
            ext.name <- rev(unlist(strsplit(code.out[code.ext[1]],"#<file:")))[1]
            ext.name <- unlist(strsplit(unlist(strsplit(ext.name,">#"))[1],":"))
            ext.chunk<-ext.name[2]; ext.name <-ext.name[1]
            ext.name.n<-nchar(ext.name)
            if(ext.name.n >4 && ".rev"==substring(ext.name,ext.name.n-3,ext.name.n)){
              ext.name<-substring(ext.name,1,ext.name.n-4)
            }

            if(is.na(as.numeric(ext.chunk))){
              # tld untersuchen
              filename<-paste(ext.name[1],".rev",sep="")
              if(!file.exists(filename)){
                cat("ERROR: file",filename,"for expansion of code chunk not found!!!\n")
                ext.file<-"Error"
              }else{
                ext.file<-try(myscan(file=filename,what="",sep="\n"))
              }
              if("Error"==substring(unlist(ext.file)[1],1,5)){
                found <-FALSE; next
              }
              ext.file <-ext.file[grep("^<<(.*)>>=",ext.file)]
              if(!is.null(ext.file)){ found<-FALSE; next }
              ext.chunk<-grep(ext.chunk,ext.file)
            }

            filename<-paste(ext.name[1],".R",sep="")
            if(!file.exists(filename)){
              cat("Warning: file",filename,"not found!!!\n")
              cat("         file",filename,"is now generated!!!\n")
              try(tangleR(ext.name[1]))
            }
            if(!file.exists(filename)){
              ext.file<-"Error"
            }else{
              ext.file<-try(myscan(file=filename,what="",sep="\n"))
            }
            if("Error"==substring(unlist(ext.file)[1],1,5)){
              found <-FALSE; next
            }

            ext.chunk<-as.numeric(ext.chunk)
            a        <-grep(paste("#", ext.chunk,":",sep=""),ext.file)[1]
            z        <-grep(paste("#:",ext.chunk,    sep=""),ext.file)[1]
            if(is.na(a)){
              found <- FALSE; next
            }
            if(a<=z) ext.file <-ext.file[a:z]

            code.out <-c(code.out[1:(code.ext[1]-1)], ext.file,
                         if(length(code.out)>(code.ext[1]+1))
                           code.out[(code.ext[1]+1):length(code.out)]
                       )

            code.ext<-code.ext[-1]

          }

        }
        code.out<-c(code.out,"##:act##")

        code.out<-gsub("DoSpCloseKl-esc",">>",gsub("DoSpOpenKl-esc","<<",code.out))

        melde("Ende Rtangle-last\n",3)
        code<-code.out[code.out!=""]
        melde("nach expand\n",3,code)
      }

      if(0<length(code)){
        tld <- paste("<","<*>",">=",sep="")
        rh <- c( get("relax.history",envir=revive.sys), list(c("@",tld,code)) )
        assign("relax.history",rh,envir=revive.sys)
      }

      code<-c("options(warn=2)",code)
      try.res <- try(eval(parse(text=code),envir=revive.env))
      options(warn=1)

      if(is.function(try.res)){
        ok <- "OK"
      } else {
        if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
        if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
        if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
        if(!is.character(ok)) { ok <- "OK" }
      }
      if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
        ok<-FALSE
        error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
        cat(error.msg,"\n")
        if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
           cat("A warning message stopped the evaluation!",
                 "If you want to\nevaluate the code anyway",
                 "evaluate code by:\n>WarnEval<")
        # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
      } else { ok<-TRUE }


      if(ok){
        if(!is.null(try.res)&&0<length(try.res)){
          if(!exists("tworkwin"))
            tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

          worktext<-TcltoWin.write(tclvalue(tkget(tworkwin,"0.0","end")))
          get("cat","package:base")(worktext,file="report-UnDo-bak.rev")

          sink(get("tmp.file.name",envir=revive.sys)
);get("print",pos="package:base")(try.res);sink()
          news<-paste(myscan(get("tmp.file.name",envir=revive.sys)
,"",sep="\n"),collapse="\n")
          if(nchar(news)>maxol.sys){
              news<-paste(substring(news,1,maxol.sys),"...",sep="\n")
          }
          news<-paste("", date(), news,"",sep="\n")
          if(!exists("toutwin"))
            toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
          pos.to.insert<-"end"
          ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
          news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
          try(tkinsert(toutwin,pos.to.insert,news))
          tksee(toutwin,"end - 0 lines")
          melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

          ##zeige Ergebnisse in neuem Fenster an>>
        }
      } else { cat("sorry, evaluation not successful!!!\n") }
    } else { cat("no code found!!!\n") }
    melde("event wird generiert",3)
    tkevent.generate(TopW,"<<Acticmds>>")

    melde("ak texthervor",1)
    tcl("markclear",tworkwin)
    tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
    tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
    tcl("marklinetypes",tworkwin)
    melde("ak texthervor",2)

    melde("fEvalRCode",2)
  }
  fFindText<-function(){
    melde("fFindText",1)
    frage<-"search string?"; set.tclvalue("tvinfo",string.sys)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>", function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        set.tclvalue("tvmess","relax")
        such.orig<-such<-string.sys<-tclvalue("tvinfo")
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        repl.pat<-gsub("(.)","\\\\\\1","^!$%&/()=?{}}+*#,.-;:\\_[") ## 070830
        repl.pat<-paste("([",repl.pat,"])",collapse="")
        such<-gsub(repl.pat,"\\\\\\1",such)
        assign("string.sys",string.sys,envir=revive.sys)
        if(nchar(such)==0) {
          tcl("findclear",tworkwin); tkfocus(tworkwin); return()
        }
 
        if(length(found<-grep(such,worktext))>0){
          cline.pos <- as.character(tkindex(tworkwin,"insert"))
          cline <- floor(as.numeric(cline.pos))
          cline.pos<-as.numeric(sub(".*\\.","",cline.pos))

          # if nothing found in the cursor line or later go to the beginning of the document
          line <- if(any(h<-(cline<found))) found[h][1] else found[1]
          # look always for the first occurence, index origin: 1
          h<-worktext[line]; hh<-1:nchar(h)
          line.pos<- -1+which(such.orig==substring(h,hh,hh-1+nchar(such.orig)))[1]
          # if cursor line and match line are equal you have to look at the positions
          if(any(found==cline)){
            h<-worktext[cline]; hh<-1:nchar(h)
            found.pos<- -1+which(such.orig==substring(h,hh,hh-1+nchar(such.orig)))
            if(any(h<-(cline.pos<found.pos))) {line<-cline; line.pos<-found.pos[h][1]}
          }

          tksee(tworkwin,paste(line,".1",sep="")); h<-paste(line,".",line.pos,sep="")
          tkmark.set(tworkwin, "insert", h); tkfocus(tworkwin)

          tktag.configure(tworkwin,"found",background="#000999fff",relief="raised")
          tcl("findmatches",tworkwin,such)
        } else set.tclvalue("tvmess",paste("Warning: search string >",
                                           such.orig,"< not found!!!"))
      } # end of function
    )
    melde("fFindText",2)
  }

  proc<-c(  # 050628
     "proc findmatches {w pattern} {",
     "$w tag remove found    1.0 end",
     "  scan [$w index end] %d numLines",
     "  for {set i 1} {$i < $numLines} {incr i} {",
     "    $w mark set last $i.0",
     "    while {[regexp -indices $pattern \\",
     "       [$w get last \"last lineend\"] indices]} {",
     "     $w mark set first \"last + [lindex $indices 0] chars\"",
     "     $w mark set last \"last  + 1 chars\\",
     "                              + [lindex $indices 1] chars\"",
     "     uplevel [$w tag add found first last]",
     "     }",
     "  }",
     "}",
      "proc findclear w {",
        "$w tag remove found    1.0 end",
      "}"
    )
  .Tcl(paste(proc,collapse="\n"))



  fExamples<-function(){
    melde("fExamples",1)
    if(paste(R.version$major,R.version$minor,sep="")<211){
      cat("sorry, for this version not implemented, have a look at the help ..."); return()
    }
    frage<-"name of R function?"; set.tclvalue("tvinfo",string.sys)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>", function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        topic<-tclvalue("tvinfo")   # topic <- "plot"
        help.text <- do.call("example",list(topic=topic,give.line=TRUE))
        help.text <- as.character(help.text)
        if(0<length(help.text)){
          help.text <- paste(help.text,collapse="\n")
          playground(code=
              paste(paste("# help examples of",topic,"-> last output will be printed:"),
                    help.text,sep="\n"))
        } else help.text <- "\nrelax warning: no examples found"
        if(is.function(get(topic))) {
            args<-deparse(args(get(topic)))
            args<-paste(args[-length(args)],collapse="\n")
            args[1] <- paste("# header of",topic,":\n#",args[1],"\n")
            args <- paste(args,collapse="\n")
        } else args <- NULL
        if(0<length(args) && is.character(args) && 0<nchar(args)) { #121217
          news <- args
          if(!exists("toutwin"))
            toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
          pos.to.insert<-"end"
          ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
          news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
          try(tkinsert(toutwin,pos.to.insert,news))
          tksee(toutwin,"end - 0 lines")
          melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

        }
        return()
    })
    melde("fExamples",2)
  }

  fHelp.R<-function(){
    melde("fHelp.R",1)
    frage<-"name of R function?"; set.tclvalue("tvinfo",string.sys)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>", function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        # help.start(browser=browser.sys); Sys.sleep(0.1);
        fns<-tclvalue("tvinfo")
        mess<-paste("documentation of",fns,"appears in the R window;",
                              "click \"ok\" and change to the R window")
        res<-tkmessageBox(message=mess,title="Help",icon="info",type="ok")
        try.res<-try(eval(parse(text=
  ##        paste("get(\"print\",\"package:base\")(help(\"",fns,"\",htmlhelp=FALSE))",sep="")
  ## wegen Argument-Änderung
          paste("get(\"print\",\"package:base\")(help(\"",fns,"\",help_type=\"text\"))",sep="")
        )))
        if(length(try.res)==0){
          mess<- paste("Warning: no documentation for",fns,"found!")
          res<-tkmessageBox(message=mess,title="Help",icon="info",type="ok")
        }
        mess<-paste("relax"); set.tclvalue("tvmess",mess)
      } # end of function
    )
    melde("fHelp.R",2)
  }

  ##definiere Funktion für Knopf: [[InsertPlot]]##
  InsertTeX<-function(){
    melde("InsertTeX",1)
    h<-gsub(" ","",gsub(":","",date()))
    h<-substring(h,4:nchar(h),4:nchar(h))
    bildname<-paste(c("p",h[6:11],h[4:5],h[1:2],h[14:15]),collapse="")
    if(!exists("tworkwin"))
      tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    start<-grep("begin-tex-code",worktext)
    start<-rev(start[start<line])[1]
    ende<-grep("end-tex-code",worktext)
    ende<-ende[ende>start][1]
    h<-worktext[start+1:ende-1]
    h<-c("\\documentclass{article}\\begin{document}",h,"\\end{document}")
    cat(h,file="tmptmp.tex",sep="\n")
    if( substring(version$os,1,5)=="linux" || substring(version$os,1,6)=="darwin" ) 
       system("echo q|latex tmptmp.tex;dvips tmptmp.dvi")
    ##lege Bild unter [[bildname]] ab##
    ##zeige Bilder im Textfenster an##
    ##show.single.plot(tworkwin,line,jpgname) 
    melde("InsertTeX",2)
  }

  fPlanRCode<-function(){
     melde("fPlanRCode",1)
     if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
     tworkwin<-get("tworkwin",envir=revive.sys)
     worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
     if(nchar(worktext)<10000){
       worktext<-strsplit(worktext,"\n")[[1]]
     }else{
       base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
       worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
     }

     news <- "\n@\n\n<<*>>=\n\n"
     ##hole ggf. [[tworkwin]]>>
     line <-floor(as.numeric(tkindex(tworkwin,"insert")))

     ##lese Arbeitsfenster auf [[worktext]] ein>>
     textstart<-grep("^@",worktext)-1; textstart<-textstart[textstart>=line][1]
     codestart<-grep("^<<(.*)>>=",worktext)-1; codestart<-codestart[codestart>=line][1]
     if(is.na(codestart))codestart<-Inf; if(is.na(textstart))textstart<-Inf
     insertline<-if(codestart==textstart) NA else min(codestart,textstart)
     anzrows<-length(unlist(strsplit(news,"\n")))
     if(is.na(insertline)){
         insertline<-"end"
         try(tkinsert(tworkwin,"end","\n"))
         try(tkinsert(tworkwin,"end",paste(news,collapse="\n")))
         tkmark.set(tworkwin, "insert","end - 2 lines")
         tksee(tworkwin,"end")  # paste(insertline+anzrows,"0",sep="."))
         insertline<-length(worktext)
     }else{
       # in einem Text-Chunks muss ein Kl-Affe eingebaut werden.
         if(length(grep("<<\\*>>=",news[1]))>0 && codestart < textstart) news<-c(news,"@\n")
         try(tkinsert(tworkwin,paste(insertline+1,"0",sep="."),paste(news,collapse="\n")))
         tkmark.set(tworkwin, "insert", paste(insertline+anzrows,"0",sep="."))
         tksee(tworkwin,paste(insertline+anzrows,"0",sep="."))
     }
     ## melde(insertline)
     melde("ak texthervor",1)
     tcl("markclear",tworkwin)
     tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
     tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
     tcl("marklinetypes",tworkwin)
     melde("ak texthervor",2)

     ##zeige Bilder im Textfenster an##
     tkfocus(tworkwin)
     melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

     if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
     tworkwin<-get("tworkwin",envir=revive.sys)
     worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
     if(nchar(worktext)<10000){
       worktext<-strsplit(worktext,"\n")[[1]]
     }else{
       base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
       worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
     }

     line <-floor(as.numeric(tkindex(tworkwin,"insert")))

     code.start<-grep("^<<(.*)>>=",worktext)
     try(if(0<length(code.start)){ 
            worktext[code.start]<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext[code.start])
            worktext[code.start]<-paste(worktext[code.start]," (",1:length(code.start),")",sep="")
     })
     if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
     tkdelete(tworkwin,"0.0","end")
     try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
     tksee(tworkwin,"end")
     melde("ak texthervor",1)
     tcl("markclear",tworkwin)
     tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
     tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
     tcl("marklinetypes",tworkwin)
     melde("ak texthervor",2)

     if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


     tkmark.set(tworkwin, "insert", paste(line,"0",sep="."))
     tksee(tworkwin,paste(line,"0",sep="."))
     tkfocus(tworkwin)

     ##zeige Bilder im Textfenster an##
     melde("fPlanRCode",2)
  }

  fRemoveOut<-function(){
    melde("fRemoveOut",1)
    worktext<-""
    if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
    tkdelete(toutwin,"0.0","end")
    try(tkinsert(toutwin,"0.0",paste(worktext,collapse="\n")))
    melde("fRemoveOut",2)
  }

  fInsert<-function(){
    melde("fInsert",1)
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    news<-tclvalue(tkget(toutwin,"0.0","end"))
    if(1<nchar(news)){
      if(length(grep("egin[{]table[}]",news))>0 &&
         length(grep("generated.*xtable.*package",news))>0){
          news<-sub("(\n%.*begin[{]table[}])","\noutput-end\n\\1",news)       
          news<-sub("(.end[{]table[}])","\\1\noutput-start",news)       
          news<-paste("\n@",paste("output-start",news,sep=""),"output-end\n", sep="\n")
          news<-sub("output-start\n+output-end","",news)
      }else{    
        news<-paste("\n@","output-start",news,"output-end\n",sep="\n")
      }
      news<-gsub("\n+","\n",news)
      ##hole ggf. [[tworkwin]]>>
      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      ##lese Arbeitsfenster auf [[worktext]] ein>>
      textstart<-grep("^@",worktext)-1; textstart<-textstart[textstart>=line][1]
      codestart<-grep("^<<(.*)>>=",worktext)-1; codestart<-codestart[codestart>=line][1]
      if(is.na(codestart))codestart<-Inf; if(is.na(textstart))textstart<-Inf
      insertline<-if(codestart==textstart) NA else min(codestart,textstart)
      anzrows<-length(unlist(strsplit(news,"\n")))
      if(is.na(insertline)){
          insertline<-"end"
          try(tkinsert(tworkwin,"end","\n"))
          try(tkinsert(tworkwin,"end",paste(news,collapse="\n")))
          tkmark.set(tworkwin, "insert","end - 2 lines")
          tksee(tworkwin,"end")  # paste(insertline+anzrows,"0",sep="."))
          insertline<-length(worktext)
      }else{
        # in einem Text-Chunks muss ein Kl-Affe eingebaut werden.
          if(length(grep("<<\\*>>=",news[1]))>0 && codestart < textstart) news<-c(news,"@\n")
          try(tkinsert(tworkwin,paste(insertline+1,"0",sep="."),paste(news,collapse="\n")))
          tkmark.set(tworkwin, "insert", paste(insertline+anzrows,"0",sep="."))
          tksee(tworkwin,paste(insertline+anzrows,"0",sep="."))
      }
      ## melde(insertline)
      melde("ak texthervor",1)
      tcl("markclear",tworkwin)
      tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
      tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
      tcl("marklinetypes",tworkwin)
      melde("ak texthervor",2)

      ##zeige Bilder im Textfenster an##
      tkfocus(tworkwin)
      melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      worktext<-""
      if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
      tkdelete(toutwin,"0.0","end")
      try(tkinsert(toutwin,"0.0",paste(worktext,collapse="\n")))
    }
    melde("fInsert",2)
  }

  fSavePlot<-function(){
    melde("fSavePlot",1)
    h<-gsub(" +"," ",date()); h<-strsplit(gsub(":","",h)," ")[[1]]
    bildname<-paste("p",h[5],"-",h[2],h[3],"-",h[4],".ps",sep="")
    if(!is.null(bildname)&&nchar(bildname)>0){
    # check name of picture
      n<-nchar(bildname<-gsub(" ","",bildname))
      bildname<-sub(".ps$","",bildname)
    # postscript:
      psname <-paste(bildname,".ps", sep="")
      try.res<-try({dev.copy(postscript,psname,horizontal=pshorizontal.sys,
                             width=psdesignwidth.sys,height=psdesignheight.sys);dev.off()})
      if(is.function(try.res)){
        ok <- "OK"
      } else {
        if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
        if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
        if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
        if(!is.character(ok)) { ok <- "OK" }
      }
      if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
        ok<-FALSE
        error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
        cat(error.msg,"\n")
        if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
           cat("A warning message stopped the evaluation!",
                 "If you want to\nevaluate the code anyway",
                 "evaluate code by:\n>WarnEval<")
        # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
      } else { ok<-TRUE }


      if(!ok){ cat("Error: *ps file not generated by dev.copy!!!\n"); return() }
      news<-paste("@\n \\begin{center}","\\includegraphics[",
                  "height=",psheight.sys,"]{",bildname,"}\\end{center}\n",sep="") #081121
    # jpeg:
      jpgname<-paste(bildname,".jpg",sep="")
      if((.Platform$OS.type=="windows")){ # width=width in pixel, 72 dpi
        try.res<-try({dev.copy(jpeg,jpgname,width=jpgdesignsize.sys*72,
                               height=jpgdesignsize.sys*72,quality=100,pointsize=7);dev.off()})
      }else{
        try.res<-try({dev.copy(bitmap,type="jpeg",jpgname,
             width=jpgdesignsize.sys,height=jpgdesignsize.sys);dev.off()})
      }
      if(is.function(try.res)){
        ok <- "OK"
      } else {
        if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
        if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
        if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
        if(!is.character(ok)) { ok <- "OK" }
      }
      if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
        ok<-FALSE
        error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
        cat(error.msg,"\n")
        if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
           cat("A warning message stopped the evaluation!",
                 "If you want to\nevaluate the code anyway",
                 "evaluate code by:\n>WarnEval<")
        # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
      } else { ok<-TRUE }


      if(!ok) cat("Error: *jpg file not generated by dev.copy!!!\n")
      news<-paste(news,'\n% <p><img src="',jpgname,'">\n@\n', sep="" )
    # ppm
      ppmname<-paste(bildname,".ppm",sep="")
      # if( <das OS ist Windows>  && 2<=nchar(ghostscript) && !Img.package.found && !no.plots){
      if((.Platform$OS.type=="windows") && 2<=nchar(ghostscript) && !no.plots){
        try.res<-try({dev.copy(bitmap,type="ppmraw",ppmname,res=ppmresolution.sys);dev.off()})
      }
      # if( <das OS ist Linux> && 2<=nchar(ghostscript) && !Img.package.found && !no.plots){ }
      if(substring(version$os,1,5)=="linux" && 2<=nchar(ghostscript) && !no.plots){ # 121113
        try.res<-try({dev.copy(bitmap,type="ppmraw",ppmname,res=ppmresolution.sys);dev.off()})
      }
    # gif+ppm:
      if(substring(version$os,1,6)=="darwin"  && !no.plots){
        gifname<-paste(bildname,".gif",sep="")
        try({dev2bitmap(type="ppmraw",ppmname,res=ppmresolution.sys); dev.off()}) #121218
        try.res<-try({system(paste("convert",jpgname,gifname))})
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(!ok) cat("Error: gif file not generated by dev.copy!!!\n")
      }
    }

    # include links
    if(!is.null(bildname)&&nchar(bildname)>0){
      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      ##hole ggf. [[tworkwin]]>>
      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      ##lese Arbeitsfenster auf [[worktext]] ein>>
      textstart<-grep("^@",worktext)-1; textstart<-textstart[textstart>=line][1]
      codestart<-grep("^<<(.*)>>=",worktext)-1; codestart<-codestart[codestart>=line][1]
      if(is.na(codestart))codestart<-Inf; if(is.na(textstart))textstart<-Inf
      insertline<-if(codestart==textstart) NA else min(codestart,textstart)
      anzrows<-length(unlist(strsplit(news,"\n")))
      if(is.na(insertline)){
          insertline<-"end"
          try(tkinsert(tworkwin,"end","\n"))
          try(tkinsert(tworkwin,"end",paste(news,collapse="\n")))
          tkmark.set(tworkwin, "insert","end - 2 lines")
          tksee(tworkwin,"end")  # paste(insertline+anzrows,"0",sep="."))
          insertline<-length(worktext)
      }else{
        # in einem Text-Chunks muss ein Kl-Affe eingebaut werden.
          if(length(grep("<<\\*>>=",news[1]))>0 && codestart < textstart) news<-c(news,"@\n")
          try(tkinsert(tworkwin,paste(insertline+1,"0",sep="."),paste(news,collapse="\n")))
          tkmark.set(tworkwin, "insert", paste(insertline+anzrows,"0",sep="."))
          tksee(tworkwin,paste(insertline+anzrows,"0",sep="."))
      }
      ## melde(insertline)
      melde("ak texthervor",1)
      tcl("markclear",tworkwin)
      tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
      tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
      tcl("marklinetypes",tworkwin)
      melde("ak texthervor",2)

      ##zeige Bilder im Textfenster an##
      tkfocus(tworkwin)
      melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      melde(paste("p", psname), "cmd.msg")
    }
    insertline<-insertline+3 # insertline is assigned during writing the pic link in the report
    base::cat(bildname,"is included into the report\n")
    if(!no.plots){
      if(substring(version$os,1,6)=="darwin" ){ picname <- gifname; type <- "gif" } else 
      ##  if(Img.package.found)      { picname <- jpgname; type <- "jpg" } else ## 121114 jpg
          if(2<=nchar(ghostscript)){ picname <- ppmname; type <- "ppm" } else type <- NA
      # there had been situations in which the (windows) file system was to slow !!
      if(!is.na(type) ) Sys.sleep(0.2) 
      if(!is.na(type) && !file.exists(picname) && file.exists(sub("....$",".ps",picname))) Sys.sleep(0.2) 
      if(!is.na(type) && !file.exists(picname) && file.exists(sub("....$",".ps",picname))) Sys.sleep(0.4) 
      if(substring(version$os,1,6)=="darwin" ) createandshow.single.plot(tworkwin,insertline,picname,type=type) else 
      ##  if(Img.package.found) createandshow.single.plot(tworkwin,insertline,picname,type=type) else ## 121113
          if(2<=nchar(ghostscript)) createandshow.single.plot(tworkwin,insertline,picname,type=type)
    }
    melde("fSavePlot",2)
  }

  # a new plot should be shown in a Tcl/Tk widget
  .Tcl( paste(""
    ,"proc createandshowimagegif {w im place} {"
    ,"global imageno"
    ,"set imageno [image create photo -format gif -file $im]" 
    ,"$w image create $place -image $imageno"
    ,"}"
    ,"proc createandshowimagejpg {w im place} {"
    ,"global imageno"
    ,"set imageno [image create photo -format jpeg -file $im]" 
    ,"$w image create $place -image $imageno"
    ,"}"
    ,"proc createandshowimageppm {w im place} {"
    ,"global imageno"
    ,"set imageno [image create photo -format ppm -file $im]" 
    ,"$w image create $place -image $imageno"
    ,"}"
  , sep="\n")) 
  assign("pic.list.sys",NULL,envir=revive.sys)
  # --------------------------------------------------------------------
  # show an image 
  try(.Tcl( paste(
    "proc showsingleimage {w imageno place} {",
       "$w image create $place -image $imageno",
    "}", sep="\n") ) )

  createandshow.single.plot<-function(textwidget,row.of.plot,picname,type="jpg"){
    melde("createandshow.single.plot",1)
    # base::cat(picname,"will be included into report field\n") 
    if(!file.exists(picname)){cat("relax warning:",picname,"not found"); return()}
    # create and show image
    place<-paste(row.of.plot,".0",sep="")
    if(type=="jpg") try.res<-try(tcl("createandshowimagejpg",textwidget,picname,place))
    if(type=="gif") try.res<-try(tcl("createandshowimagegif",textwidget,picname,place))
    if(type=="ppm") try.res<-try(tcl("createandshowimageppm",textwidget,picname,place))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    # update pic list
    if(try.res=="ok"){
      imageno<-tclvalue(.Tcl("set imageno"))
      pic.list<-rbind(get("pic.list.sys",envir=revive.sys),c(picname,imageno))
      assign("pic.list.sys",pic.list,envir=revive.sys)
    }
    melde("createandshow.single.plot",2)
  }

  show.plots.again<-function(tworkwin){
      melde("show.plots.again",1)
      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      rows   <-grep("<img src=",worktext); if(0==length(rows))return()
      pic.names.report  <-unlist(lapply(strsplit(worktext[rows],"<img src=\""), 
                               function(x){ x<-x[2]; strsplit(x,"\"")[[1]][1] }))
      pic.list<-get("pic.list.sys",envir=revive.sys)
      ind<-match(sub("....$","",pic.names.report),sub("....$","",pic.list[,1]))
      image.ok<-!is.na(ind); ind<-ind[image.ok]
      if(0==length(ind))return()
      place<-paste(rows[image.ok]-1,".0",sep=""); imageno<-pic.list[ind,2]
      for(i in seq(imageno)){ 
           try(tcl("showsingleimage",tworkwin,imageno[i],place[i]))
      }
      melde("show.plots.again",2)
  }

  # create and show all plots
  createandshow.all.plots<-function(tworkwin){
     melde("createandshow.all.plots",1)
      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      rows   <-grep("<img src=",worktext); if(0==length(rows))return()
      pic.names.report<-unlist(lapply(strsplit(worktext[rows],"<img src=\""), 
                               function(x){ x<-x[2]; strsplit(x,"\"")[[1]][1] }))
      place<-paste(rows-1,".0",sep="")
      if(0==length(pic.names.report)) return()
      # fix pic name / file 
      if(substring(version$os,1,6)=="darwin" ) picname<-sub(".jpg$",".gif",pic.names.report) else 
      ##  if(Img.package.found) picname<-pic.names.report  else ## 121113 jpg
          if(2<=nchar(ghostscript)) picname<-sub(".jpg$",".ppm",pic.names.report)
      # create missing files by postscript versions of the images
      if(!substring(version$os,1,6)=="darwin" ){ #12114
        for(i in seq(picname)){
           melde(paste("generation of ppm",picname[i]),3)
           # ggf. new generation of plot
           if(file.exists(sub("....$",".ps",picname[i])) && (!file.exists(picname[i]))){
             ##  if(Img.package.found) 
             ##         pstojpg(sub("....$",".ps",picname[i]),picname[i]) else
             res <- try(pstoppm(sub("....$",".ps",picname[i]),picname[i])) #121114
           }
        }
      }
      # 
     image.nos<-rep("xx",length(picname))
     for(i in seq(picname)){
       melde(paste("integration of ppm",picname[i]),3) # resolution
       if(!file.exists(picname[i])){
         next
       }
       if(substring(version$os,1,6)=="darwin" ) 
         try.res<-try(tcl("createandshowimagegif",tworkwin,picname[i],place[i])) 
       else
         ## if(Img.package.found)
         ##   try.res<-try(tcl("createandshowimagejpg",tworkwin,picname[i],place[i]))
         ## else ## 121113 jpg
         if(2<=nchar(ghostscript)) 
           try.res<-try(tcl("createandshowimageppm",tworkwin,picname[i],place[i])) 
       if(is.function(try.res)){
         ok <- "OK"
       } else {
         if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
         if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
         if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
         if(!is.character(ok)) { ok <- "OK" }
       }
       if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
         ok<-FALSE
         error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
         cat(error.msg,"\n")
         if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
            cat("A warning message stopped the evaluation!",
                  "If you want to\nevaluate the code anyway",
                  "evaluate code by:\n>WarnEval<")
         # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
       } else { ok<-TRUE }


       if(try.res=="ok") image.nos[i]<-tclvalue(.Tcl("set imageno"))
     }
     pic.list<-cbind(picname,image.nos)
     pic.list<-pic.list[pic.list[,2]!="xx",,drop=FALSE]
     assign("pic.list.sys",pic.list,envir=revive.sys)
     melde("createandshow.all.plots",2)
  } 

  exclude.plots<-function(tworkwin){
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    rows<-grep("<img src=",worktext); if(0==length(rows))return()
    names  <-unlist(lapply(strsplit(worktext[rows],"<img src=\""), 
                    function(x){ x<-x[2]; strsplit(x,"\"")[[1]][1] }))
    rows <- rows - 1   # + 1
    for(i in seq(along=rows)){
      anf<-paste(rows[i],".0",sep=""); end<-paste(rows[i],".end",sep="")
      a<-tclvalue(tkget(tworkwin,anf,end))
      tkdelete(tworkwin,anf,end)
      tkinsert(tworkwin,anf,paste(a,collapse="\n")) # 040922
  #    tkdelete(tworkwin,paste(rows[i],".0",sep=""),paste(rows[i],".1",sep=""))
    }
  }

  pstojpg<-function(psname,jpgname){
    return() # im Moment nicht im Dienst
    if(substring(version$os,1,5)=="linux"){
      gsexe <- Sys.getenv("R_GSCMD")
      if(is.null(gsexe) || nchar(gsexe) == 0) {
        gsexe <- "gs"; rc <- system(paste(gsexe, "-help > /dev/null"))
        if (rc != 0) return()
      }
    }
    if((.Platform$OS.type=="windows")){
      if(exists("ghostscript")&&2<=nchar(ghostscript)&&0<length(grep("[A-Za-z]",ghostscript))) 
        gsexe<-ghostscript else return()
    }
    type<-"jpeg"; width<-height<-13; res<-30
    cmd <- paste(gsexe, " -dNOPAUSE -dBATCH -q -sDEVICE=", type,
                        " -r", res,
                        " -g",ceiling(res*width),"x",ceiling(res*height),
                        " -sOutputFile=", jpgname, "  ",psname, sep = "")
    try(system(cmd)); invisible()
  }
  pstoppm<-function(psname,ppmname){
    # return() # im Moment nicht im Dienst // ginge auch mit convert
    melde("pstoppm startet",3)
    if(substring(version$os,1,5)=="linux"){
      gsexe <- Sys.getenv("R_GSCMD")
      if(is.null(gsexe) || nchar(gsexe) == 0) {
        gsexe <- "gs"; rc <- system(paste(gsexe, "-help > /dev/null"))
        if (rc != 0) return()
      }
    }
    if((.Platform$OS.type=="windows")){
      if(exists("ghostscript")&&2<=nchar(ghostscript)&&0<length(grep("[A-Za-z]",ghostscript))) 
        gsexe<-ghostscript else return()
    }
    type<-"ppmraw"; res<-ppmresolution.sys # 121114
    cmd <- paste(gsexe, " -dNOPAUSE -dBATCH -q -sDEVICE=", type,
                        " -r", res,
                        " -sOutputFile=", ppmname, "  ",psname, sep = "")
    melde(cmd,3)
    try(system(cmd))
    if(substring(version$os,1,5)=="linux" && 0 < nchar(system("which pnmcrop",intern=TRUE))){ # 121114
      system(paste("mv ",ppmname," tmptmp.ppm; pnmcrop tmptmp.ppm > ",ppmname,"; rm tmptmp.ppm"))    
    }                       
    invisible()
  }

  ##definiere Funktion für Knopf: [[CopyToEnd]]##
  fTrashROutput<-function(){
    melde("fTrashROutput",1)
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    out.end  <-grep("^output-end",worktext)
    out.end  <-out.end[out.end>=(line-1)][1]
    if(is.na(out.end)) return()
    out.start <-grep("^output-start",worktext)
    out.start<-rev(out.start[out.start<out.end])[1]
    if(is.na(out.end)) return()
    code.start   <-grep("^<<(.*)>>=",worktext)
    code.start <- rev(code.start[code.start<out.end])[1]
    if(is.na(code.start)) return()
    if(code.start>out.start) return()
    if("@"==worktext[out.start-1]) out.start<-out.start-1
    if(""==worktext[out.start-1]) out.start<-out.start-1
    if(""==worktext[out.start-1]) out.start<-out.start-1
    trash<-worktext[out.start:out.end]
    worktext<-worktext[-(out.start:out.end)]
    if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
    tkdelete(tworkwin,"0.0","end")
    try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
    tksee(tworkwin,"end")
    melde("ak texthervor",1)
    tcl("markclear",tworkwin)
    tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
    tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
    tcl("marklinetypes",tworkwin)
    melde("ak texthervor",2)

    if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


    tkfocus(tworkwin)
    tkmark.set(tworkwin,"insert",paste(out.start-1,".0",sep=""))
    tksee(tworkwin,"insert")
    melde("fTrashROutput",2)
  }

  fWarnEval<-function(){ # vorher EvalCursorChunk
    melde("fWarnEval",1)
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    if(any(h<-(substring(worktext,1,1)==">"))){ # 091026
      worktext[h]<-sub("^.","@SpLiTtE?<<*>>=SpLiTtE?",worktext[h])
      worktext<-substring(unlist(strsplit(paste("",worktext),"SpLiTtE")),2)
      if(is.vector(line)) line<-line+2*sum(which(h)<=line)
    }

    code.start<-grep("^<<(.*)>>=",worktext)
    if(0==length(code.start)){cat("Warning: no code found!!!\n");return()}
    code.start<-code.start[code.start<=line]
    if(0<length(code.start)){
      if((code.start<-max(code.start)+1)>length(worktext)) return()
      code.end  <-c(grep("^@",worktext),1+length(worktext))
      code.end  <-min(code.end[code.end>code.start])-1
      code<-worktext[code.start:code.end]
      code<-code[code!=""]
      if(length(weg.ab<-grep("^<<(.*)>>=",code))>0) code<-code[-(weg.ab:length(code))]
      if(length(code)==0 || code[1]=="@") code<-" "
      melde("code:",3,code,"\n")

    }
    if(0<length(code)){
      melde("vor tangle - Code:\n",3,code)
      if(length(grep("<<(.*)>>",code))>0 || length(grep(">#",code))>0){
        code.a    <- grep("^<<(.*)>>=",worktext)
        code.z    <- grep("^@",worktext)
        code.z    <- unlist(sapply(code.a ,function(x,y) y[y>x][1], code.z))
        if(any(h<-is.na(code.z))) code.z<-code.z[!h]
        ################
        # code.n    <- length(worktext)
        # change    <- rep(0,code.n); change[c(code.a ,code.z)]<-1
        # code.ch   <- worktext[1==(cumsum(change)%%2)]
        ## 080311
        ind<-rep(0,length(worktext)); copy<-0
        for(i in seq(ind)){
          if(i %in% code.z) copy<-0
          if(i %in% code.a) copy<-1
          ind[i]<-copy
        }
        code.ch   <- worktext[ind==1]
        ## cat("input-text, Anfaenge, Code.chunks"); print(worktext); print(code.a); print(code.ch)

        code.n    <- length(code.ch)

        code.ch<-gsub("@>>","DoSpCloseKl-esc",gsub("@<<","DoSpOpenKl-esc",code.ch))

        code.ch<-gsub("(.*)<<(.*)>>=(.*)","cOdEdEf\\2",code.ch)
        repeat{
          if(0==length(cand<-grep("<<(.*)>>",code.ch))) break
          code.ch<-unlist(strsplit(gsub("(.*)<<(.*)>>(.*)",
                     "\\1bReAkuSeChUnK\\2bReAk\\3",code.ch),"bReAk"))
        }
        code.ch<-code.ch[code.ch!=""]
        code.n<-length(code.ch)
        melde("code.ch:",3,code.ch,"\n")

        line.typ  <-rep("C",code.n)
        code.a    <-grep("cOdEdEf",code.ch)
        code.ch[code.a]<-substring(code.ch[code.a],8)
        line.typ[code.a]<-"D"
        code.use    <-grep("uSeChUnK",code.ch)
        code.ch[code.use]<-substring(code.ch[code.use],9)
        line.typ[code.use]<-"U"
        code.ext  <-grep("#<file",code.ch)
        line.typ[code.ext]<-"E"
        melde("code.ch:",3,code.ch,"\n")

        code.out<-"##act:##"

        def.names<-code.ch[code.a]
        use.names<- if(length(code.use)>0) code.ch[code.use] else NULL
        code.z<-c(if(length(code.a)>1) code.a[-1]-1, code.n)
        code.ch<-paste(line.typ,code.ch,sep="")
        melde("code.ch:",3,code.ch,"\n")

        melde("vor expand - Code:\n",3,code)
        melde("bearbeite aktuellen Chunk\n",3)
        ###<bestimme Cursorzeile [[line]] von [[tworkwin]]> not a good idea###
        ch.no<-length(grep("^<<(.*)>>=",worktext[1:line]))

        rows      <-c((code.a[ch.no]+1),code.z[ch.no])
        if(all(!is.na(rows))&&rows[1]<=rows[2]){
          rows<-rows[1]:rows[2]
          code.stack<-code.ch[rows]
          max.depth.refinements<-500; i<-1
          repeat{
             if((i<-i+1)>max.depth.refinements){ 
                 cat("ERROR: maximal number of expandations (",max.depth.refinements,
                     ") exceeded\n --- perhaps a unintended recursion ???")
                 return()
             }
             if(0==length(code.stack))break
             typ<-substring(code.stack[1],1,1)
             if("C"==typ||"E"==typ){
               n.lines<-sum(cumprod("C"==substring(code.stack,1,1)))
               code.out<-c(code.out, substring(code.stack[1:n.lines],2))
               code.stack<-code.stack[-(1:n.lines)]
             }
             if(length(code.stack)>0 && "U"==substring(code.stack[1],1,1)){
               if(any(found<-def.names==substring(code.stack[1],2))){
                 found<-seq(along=def.names)[found]; rows<-NULL
                 for(no in found){
                   if((code.a[no]+1)<=code.z[no]) rows<-c(rows,(code.a[no]+1):code.z[no])
                 }
                 code.stack<-c(code.ch[rows],code.stack[-1])
                 melde("found",0,found)
               }else{code.stack<-code.stack[-1]}
             }

          }
        }
        if(length(code.ext)>0){
          code.out<-code.out[code.out!=""]
          code.ext<-rev(grep(">#",code.out))
          found<-TRUE
          repeat{
            if(length(code.ext)==0) break

            if(!found){
              code.out[code.ext[1]]<-paste("# ??",code.out[code.ext[1]])
              cat("ERROR: External Chunk",code.out[code.ext[1]],"not found!!!\n")
              code.ext<-code.ext[-1]
            }

            found<-TRUE
            ext.name <- rev(unlist(strsplit(code.out[code.ext[1]],"#<file:")))[1]
            ext.name <- unlist(strsplit(unlist(strsplit(ext.name,">#"))[1],":"))
            ext.chunk<-ext.name[2]; ext.name <-ext.name[1]
            ext.name.n<-nchar(ext.name)
            if(ext.name.n >4 && ".rev"==substring(ext.name,ext.name.n-3,ext.name.n)){
              ext.name<-substring(ext.name,1,ext.name.n-4)
            }

            if(is.na(as.numeric(ext.chunk))){
              # tld untersuchen
              filename<-paste(ext.name[1],".rev",sep="")
              if(!file.exists(filename)){
                cat("ERROR: file",filename,"for expansion of code chunk not found!!!\n")
                ext.file<-"Error"
              }else{
                ext.file<-try(myscan(file=filename,what="",sep="\n"))
              }
              if("Error"==substring(unlist(ext.file)[1],1,5)){
                found <-FALSE; next
              }
              ext.file <-ext.file[grep("^<<(.*)>>=",ext.file)]
              if(!is.null(ext.file)){ found<-FALSE; next }
              ext.chunk<-grep(ext.chunk,ext.file)
            }

            filename<-paste(ext.name[1],".R",sep="")
            if(!file.exists(filename)){
              cat("Warning: file",filename,"not found!!!\n")
              cat("         file",filename,"is now generated!!!\n")
              try(tangleR(ext.name[1]))
            }
            if(!file.exists(filename)){
              ext.file<-"Error"
            }else{
              ext.file<-try(myscan(file=filename,what="",sep="\n"))
            }
            if("Error"==substring(unlist(ext.file)[1],1,5)){
              found <-FALSE; next
            }

            ext.chunk<-as.numeric(ext.chunk)
            a        <-grep(paste("#", ext.chunk,":",sep=""),ext.file)[1]
            z        <-grep(paste("#:",ext.chunk,    sep=""),ext.file)[1]
            if(is.na(a)){
              found <- FALSE; next
            }
            if(a<=z) ext.file <-ext.file[a:z]

            code.out <-c(code.out[1:(code.ext[1]-1)], ext.file,
                         if(length(code.out)>(code.ext[1]+1))
                           code.out[(code.ext[1]+1):length(code.out)]
                       )

            code.ext<-code.ext[-1]

          }

        }
        code.out<-c(code.out,"##:act##")

        code.out<-gsub("DoSpCloseKl-esc",">>",gsub("DoSpOpenKl-esc","<<",code.out))

        melde("Ende Rtangle-last\n",3)
        code<-code.out[code.out!=""]
        melde("nach expand\n",3,code)
      }

      if(0<length(code)){
        tld <- paste("<","<*>",">=",sep="")
        rh <- c( get("relax.history",envir=revive.sys), list(c("@",tld,code)) )
        assign("relax.history",rh,envir=revive.sys)
      }

      try.res <- try(eval(parse(text=code),envir=revive.env))
      #options(warn=0)
      #options(warn.expression=NULL)

      if(is.function(try.res)){
        ok <- "OK"
      } else {
        if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
        if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
        if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
        if(!is.character(ok)) { ok <- "OK" }
      }
      if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
        ok<-FALSE
        error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
        cat(error.msg,"\n")
        if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
           cat("A warning message stopped the evaluation!",
                 "If you want to\nevaluate the code anyway",
                 "evaluate code by:\n>WarnEval<")
        # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
      } else { ok<-TRUE }


      if(ok){
        if(!is.null(try.res)&&0<length(try.res)){
          if(!exists("tworkwin"))
            tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

          worktext<-TcltoWin.write(tclvalue(tkget(tworkwin,"0.0","end")))
          get("cat","package:base")(worktext,file="report-UnDo-bak.rev")

          sink(get("tmp.file.name",envir=revive.sys)
);get("print",pos="package:base")(try.res);sink()
          news<-paste(myscan(get("tmp.file.name",envir=revive.sys)
,"",sep="\n"),collapse="\n")
          if(nchar(news)>maxol.sys){
              news<-paste(substring(news,1,maxol.sys),"...",sep="\n")
          }
          news<-paste("", date(), news,"",sep="\n")
          if(!exists("toutwin"))
            toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
          pos.to.insert<-"end"
          ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
          news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
          try(tkinsert(toutwin,pos.to.insert,news))
          tksee(toutwin,"end - 0 lines")
          melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

          ##zeige Ergebnisse in neuem Fenster an>>
        }
      } else { cat("sorry, evaluation not successful!!!\n") }
    } else { cat("no code found!!!\n") }
    melde("event wird generiert",3)
    tkevent.generate(TopW,"<<Acticmds>>")

    melde("ak texthervor",1)
    tcl("markclear",tworkwin)
    tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
    tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
    tcl("marklinetypes",tworkwin)
    melde("ak texthervor",2)

    ##lese Arbeitsfenster auf [[worktext]] ein##
    ##aktualisiere Code-Chunk-Zähler und schreibe [[worktext]]##
    ##generiere ein Ereignis zum Anzeigen von Warnungen##
    ##zeige ggf. Warnungen an#
    melde("fWarnEval",2)
  }

  fUp<-function(){
    melde("fUp",1)
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    code.start<-grep("^<<(.*)>>=",worktext)
    code.start<-rev(code.start[code.start<line])[1]
    if(!is.na(code.start)) tkmark.set(tworkwin, "insert", paste(code.start,".0",sep=""))
    tksee(tworkwin,"insert - 7 lines")  #;tksee(tworkwin,"insert + 7 lines")
    tksee(tworkwin,"insert + 4 lines")
    #090706 # tksee(tworkwin,"insert")
    melde("fUp",2)
  }

  fDown<-function(){
    melde("fDown",1)
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    code.start<-grep("^<<(.*)>>=",worktext)
    code.start<-code.start[code.start>line][1]
    if(!is.na(code.start)) tkmark.set(tworkwin, "insert", paste(code.start,".0",sep=""))
    tksee(tworkwin,"insert - 7 lines")  #;tksee(tworkwin,"insert + 7 lines") 
    tksee(tworkwin,"insert + 4 lines") 
    #090706 # tksee(tworkwin,"insert")
    melde("fDown",2)
  }

  Exit<-function(){
    melde("Exit",1)
    res<-tkmessageBox(message=
              if(language=="german") "Report-Manager ohne erneute Speicherung beenden?"
              else "Quit RELAX without saving again?",
                      title="Exit",icon="warning",type="yesnocancel",default="no")
    if("externalptr"==mode(res))  res<-tclvalue(res)
    if(res=="cancel") return()
    if(res=="no") SaveReport()
    set.tclvalue("tvexit","fertig"); tkdestroy(TopW)
    ## remove("print",pos=which(path.package("relax")==searchpaths())) ## 111103 # 130325 .path.package defunct
    melde("==================================================\n")
    melde("RELAX --- EXIT                    \n")
    # melde("copy of report saved as: \n")
    # melde(" report-UnDo-bak.rev                           \n")
    melde("restart RELAX by: \n relax()           \n")
    melde("==================================================\n")
    melde("q","cmd.msg")
    melde("Exit",2)
  }

  SetPSDesignWidth<-function(){
    melde("SetPSDesignWidth",1)
    frage<-"ps design width?"; set.tclvalue("tvinfo",psdesignwidth.sys)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>",
      function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        h<-as.numeric(tclvalue("tvinfo"))
        assign("psdesignwidth.sys",h,envir=revive.sys)
        melde(paste('> assign("psdesignwidth.sys","',psdesignwidth.sys,'",envir=revive.sys)',sep="")
                  , "cmd.msg")
      }
    )
    melde("SetPSDesignWidth",2)
  }
  SetPlotHeight<-function(){
    melde("SetPlotHeight",1)
    frage<-"ps height?"; set.tclvalue("tvinfo",psheight.sys)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>",
      function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        psheight.sys<-tclvalue("tvinfo")
        assign("psheight.sys",psheight.sys,envir=revive.sys)
        melde(paste('> assign("psheight.sys","',psheight.sys,'",envir=revive.sys)',sep="")
              ,"cmd.msg")
      }
    )
    melde("SetPlotHeight",2)
  }
  SetPSDesignHeight<-function(){
    melde("SetPSDesignHeight",1)
    frage<-"ps design width?"; set.tclvalue("tvinfo",psdesignheight.sys)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>",
      function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        h<-as.numeric(tclvalue("tvinfo"))
        assign("psdesignheight.sys",h,envir=revive.sys)
        melde(paste('> assign("psdesignheight.sys","',psdesignheight.sys,'",envir=revive.sys)',sep="")
                  , "cmd.msg")
      }
    )
    melde("SetPSDesignHeight",2)
  }

  SetPSRotation<-function(){
    melde("SetPSRotation",1)
    choices<-c("vertical / standard","horizontal") 
    activate<-1+pshorizontal.sys
    menuwidget <- tkmenu(TopW); set.tclvalue("tvchoice","0")
    for(i in choices) tkadd(menuwidget,"radiobutton",label=i,variable="tvchoice")
    tkpost(menuwidget,"0","0"); tkactivate(menuwidget,activate)
    tkbind(menuwidget,"<Escape>",function(){
           tkdestroy(menuwidget);set.tclvalue("tvchoice","0")})
    if( (! (.Platform$OS.type=="windows")) ) tkwait.variable("tvchoice")

    choice <- tclvalue("tvchoice")
    choice <- if(choice!="0") choice<-which(choice==choices) else 0

    if(choice!=0){
        pshorizontal.sys<-2==choice
        assign("pshorizontal.sys",pshorizontal.sys,envir=revive.sys)
        melde(paste(
             '> assign("pshorizontal.sys","',pshorizontal.sys,'",envir=revive.sys)'
             ,sep=""),"cmd.msg")
    }
    melde("SetPSRotation",2)
  }
  SetPPMResolution<-function(){ #121114
    melde("SetPPMResolution",1)
    frage<-"ppm resolution (10..100)?"; set.tclvalue("tvinfo",ppmresolution.sys)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>",
      function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        h<-as.numeric(tclvalue("tvinfo"))
        assign("ppmresolution.sys",h,envir=revive.sys)
        melde(paste(
    '> assign("ppmresolution.sys","',ppmresolution.sys,'",envir=revive.sys)',
           sep=""), "cmd.msg")
      }
    )
    melde("SetPPMResolution",2)
  }
  SetJPGSize<-function(){
    melde("SetJPGSize",1)
    frage<-"jpg design size (3..8)?"; set.tclvalue("tvinfo",jpgdesignsize.sys)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>",
      function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        h<-as.numeric(tclvalue("tvinfo"))
        assign("jpgdesignsize.sys",h,envir=revive.sys)
        melde(paste(
    '> assign("jpgdesignsize.sys","',jpgdesignsize.sys,'",envir=revive.sys)',
           sep=""), "cmd.msg")
      }
    )
    melde("SetJPGSize",2)
  }

  SetFontType<-function(){
    melde("SetFontType",1)
    choices<-c("helvetica","courier","new century schoolbook","times")
    activate<-which(strsplit(tfont.sys,"-")[[1]][3]==choices)
    menuwidget <- tkmenu(TopW); set.tclvalue("tvchoice","0")
    for(i in choices) tkadd(menuwidget,"radiobutton",label=i,variable="tvchoice")
    tkpost(menuwidget,"0","0"); tkactivate(menuwidget,activate)
    tkbind(menuwidget,"<Escape>",function(){
           tkdestroy(menuwidget);set.tclvalue("tvchoice","0")})
    if( (! (.Platform$OS.type=="windows")) ) tkwait.variable("tvchoice")

    choice <- tclvalue("tvchoice")
    choice <- if(choice!="0") choice<-which(choice==choices) else 0

    if(choice!=0){
      font<-choices[choice]
      tfont.sys<-sub("be-.+-Me",paste("be-",font,"-Me",sep=""),tfont.sys)
      tkconfigure(tworkwin,font=tfont.sys)
      assign("tfont.sys",tfont.sys,envir=revive.sys)
      ## if(exists("trevwin")) tkconfigure(trevwin, font=tfont.sys)
    }
    melde(paste("> font type changed"))
    melde("SetFontType",2)
  }
  SetFontSize<-function(){
    melde("SetFontSize",1)
    sizes<-c("8-80","10-100","12-120","14-140","18-180","24-240","*-2000")
    choices<-c("1 tiny","2 very small","3 small","4 normal","5 large",
               "6 very large","7 huge")
    activate<-which(substring(sub("(.*)Normal..","",tfont.sys),1,3)==substring(sizes,1,3))
    menuwidget <- tkmenu(TopW); set.tclvalue("tvchoice","0")
    for(i in choices) tkadd(menuwidget,"radiobutton",label=i,variable="tvchoice")
    tkpost(menuwidget,"0","0"); tkactivate(menuwidget,activate)
    tkbind(menuwidget,"<Escape>",function(){
           tkdestroy(menuwidget);set.tclvalue("tvchoice","0")})
    if( (! (.Platform$OS.type=="windows")) ) tkwait.variable("tvchoice")

    choice <- tclvalue("tvchoice")
    choice <- if(choice!="0") choice<-which(choice==choices) else 0

    if(choice!=0){
      size<-sizes[choice]
      tfont.sys<-sub("al--.+-.+-\\*",paste("al--",size,"-*",sep=""),tfont.sys)
      tkconfigure(tworkwin,font=tfont.sys)
      assign("tfont.sys",tfont.sys,envir=revive.sys)
      outfont.sys<-sub("al--.+-.+-\\*",paste("al--",size,"-*",sep=""),outfont.sys)
      assign("outfont.sys",outfont.sys,envir=revive.sys)
      melde("ak texthervor",1)
      tcl("markclear",tworkwin)
      tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
      tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
      tcl("marklinetypes",tworkwin)
      melde("ak texthervor",2)

      tkconfigure(toutwin,font=outfont.sys)
      melde(paste("> font size changed"))
    }
    melde("SetFontSize",2)
  }

  SetRelaxWinSize<-function(){
    melde("SetRelaxWinSize",1)
    frage<-"Size of relax window (width x height)?"
    set.tclvalue("tvinfo",paste(tclvalue(tkwinfo("width",TopW)),tclvalue(tkwinfo("height",TopW)),sep="x"))
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>",
      function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        h<-tclvalue("tvinfo"); h <- gsub("[^0-9x]*","",h); h <- unlist(strsplit(h,"x"))
        if(2!=length(h) || any(0==nchar(h)) ) return()
        h <- try(as.numeric(h)); if(class(h)=="try-error") return(); h <- paste(pmax(h,250),collapse="x")
        tkwm.geometry(TopW,h); cat("relax window set to size",h)
        melde(paste("> maxol.sys <-",maxol.sys),"cmd.msg")
      }
    )
    melde("SetRelaxWinSize",2)
  }

  SetOutputLength<-function(){
    melde("SetOutputLength",1)
    frage<-"maximum number of output character?"; set.tclvalue("tvinfo",maxol.sys)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>",
      function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        h<-as.numeric(tclvalue("tvinfo"))
        if(!is.na(h)) assign("maxol.sys",h,envir=revive.sys)
        melde(paste("> maxol.sys <-",maxol.sys),"cmd.msg")
      }
    )
    melde("SetOutputLength",2)
  }

  ConfigRelax<-function(){
    melde("ConfigRelax",1)
        .newl<-tktoplevel();tkwm.geometry(.newl,"+0+15");tkpack(tt<-tktext(.newl))
        tkpack(fcmd<-tkframe(.newl))
        Save<-tkbutton(fcmd,width=15,text="save settings",command=function(){
                   newsettings<-paste(tclvalue(tkget(tt,"0.0","end")),"\n")
                   filename<-file.path(path.package("relax"),"config/settings.relax") # 130325 .path.package defunct
                   try(cat(file=filename,newsettings))
                   tkmessageBox(title="", icon="warning",
                            message="to get new settings work, relax has to be restarted")
                   print("new settings saved")
        })
        Exit<-tkbutton(fcmd,width=15,text="Exit",command=function(){
                   tkdestroy(.newl); set.tclvalue("tvscandone",2)
        })
        ReloadOld<-tkbutton(fcmd,width=15,text="reload old",command=function(){
                   settings<-get("settings",envir=revive.sys)
                   tkdelete(tt,"0.0","end")
                   try(tkinsert(tt,"0.0",paste(settings,collapse="\n")))
        })
        LoadOrig<-tkbutton(fcmd,width=15,text="load original",command=function(){
          filename<-file.path(path.package("relax"),"config/settings.init") # 130325 .path.package defunct
          settings<-scan(file=filename,what="",sep="\n")
          try(tkinsert(tt,"0.0",paste(settings,collapse="\n")))
        })
        tkpack(Save,ReloadOld,LoadOrig,Exit,side="left")
        tkwm.title(.newl,"view and change settings of relax: config/settings.relax")
        filename<-file.path(path.package("relax"),"config/settings.relax") # 130325 .path.package defunct
        settings<-scan(file=filename,what="",sep="\n")
        assign("settings",settings,envir=revive.sys)
        try(tkinsert(tt,"0.0",paste(settings,collapse="\n")))
        tkbind(.newl,"<Escape>", function(){
                   tkdestroy(.newl); set.tclvalue("tvscandone",2)
                 })
        tkfocus(.newl)
        set.tclvalue("tvmess","relax")
    melde("ConfigRelax",2)
  }
  FindReportText<-fFindText
  GoToLine<-function(){
    melde("GoToLine",1)
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    frage<-"GOTO line?"; set.tclvalue("tvinfo",line)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>", function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        line<-string.sys<-tclvalue("tvinfo")[1]
        assign("string.sys",string.sys,envir=revive.sys)
        line<-as.numeric(line)
        if(!is.na(line)){
          tksee(tworkwin,h<-paste(line,".1",sep=""))
          tkmark.set(tworkwin, "insert", h); tkfocus(tworkwin)
        } else set.tclvalue("tvmess",paste("Warning: line",tclvalue("tvinfo"),"< not found!!!"))
      } # end of function
    )
    melde("GoToLine",2)
  }
  Replace<-function(){ #100917
    melde("Replace",1)
    tworkwin<-get("tworkwin",envir=revive.sys)
    .Tcl(paste("searchrep",tworkwin$ID))
    melde("Replace",2)
  }

  ### if(but.Wizardry=="all"){ }
   WebReport<-function(){
    melde("WebReport",1)
        if(file.exists(filename<-file.path(getwd(),workname.sys))){
          res<-tkmessageBox(message=
                     if(language=="german") 
                        paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                     else paste(filename,"exists. Do you want to replace it?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
          if("externalptr"==mode(res))  res<-tclvalue(res)
          if(res=="cancel")return()
          if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
        }
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
        worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
        worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
        if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
           worktext<-c(worktext,"@\n\\end{document}")
        worktext<-TcltoWin.write(worktext)
        try.res <- try(cat(worktext,file=filename,sep="\n"))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){ cat(paste("file",filename,"saved\n")) }

        if(ok){
          melde(paste("report file",filename,"saved\n"),0)
          melde(paste("w",workname.sys),"cmd.msg")
        } else {
          cat("ERROR: write operation failed!!!\n"); return()
        }





    try(weaveR(workname.sys,replace.umlaute=replace.umlaute.sys))
    cat(sub(".rev$",".tex",workname.sys),"generated\n")
    try(tangleR(workname.sys))
    cat(sub(".rev$",".R",workname.sys),"generated\n")
    melde("WebReport",2)
   }
   WeaveReport<-function(){
    melde("WeaveReport",1)
        if(file.exists(filename<-file.path(getwd(),workname.sys))){
          res<-tkmessageBox(message=
                     if(language=="german") 
                        paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                     else paste(filename,"exists. Do you want to replace it?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
          if("externalptr"==mode(res))  res<-tclvalue(res)
          if(res=="cancel")return()
          if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
        }
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
        worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
        worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
        if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
           worktext<-c(worktext,"@\n\\end{document}")
        worktext<-TcltoWin.write(worktext)
        try.res <- try(cat(worktext,file=filename,sep="\n"))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){ cat(paste("file",filename,"saved\n")) }

        if(ok){
          melde(paste("report file",filename,"saved\n"),0)
          melde(paste("w",workname.sys),"cmd.msg")
        } else {
          cat("ERROR: write operation failed!!!\n"); return()
        }





    try(weaveR(workname.sys,replace.umlaute=replace.umlaute.sys))
    cat(sub(".rev$",".tex",workname.sys),"generated\n")
    melde("WeaveReport",2)
   }
   ProcessWithSexpr<-function(){
    melde("ProcessWithSexpr",1)
        if(file.exists(filename<-file.path(getwd(),workname.sys))){
          res<-tkmessageBox(message=
                     if(language=="german") 
                        paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                     else paste(filename,"exists. Do you want to replace it?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
          if("externalptr"==mode(res))  res<-tclvalue(res)
          if(res=="cancel")return()
          if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
        }
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
        worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
        worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
        if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
           worktext<-c(worktext,"@\n\\end{document}")
        worktext<-TcltoWin.write(worktext)
        try.res <- try(cat(worktext,file=filename,sep="\n"))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){ cat(paste("file",filename,"saved\n")) }

        if(ok){
          melde(paste("report file",filename,"saved\n"),0)
          melde(paste("w",workname.sys),"cmd.msg")
        } else {
          cat("ERROR: write operation failed!!!\n"); return()
        }





    try(weaveR(workname.sys,replace.umlaute=replace.umlaute.sys,eval_Sexpr = TRUE))
    cat(sub(".rev$",".tex",workname.sys),"generated\n")
    LatexReport()
    melde("ProcessWithSexpr",2)
   }
   WeaveReportNoCode<-function(){
    melde("WeaveReportNoCode",1)
        if(file.exists(filename<-file.path(getwd(),workname.sys))){
          res<-tkmessageBox(message=
                     if(language=="german") 
                        paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                     else paste(filename,"exists. Do you want to replace it?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
          if("externalptr"==mode(res))  res<-tclvalue(res)
          if(res=="cancel")return()
          if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
        }
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
        worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
        worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
        if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
           worktext<-c(worktext,"@\n\\end{document}")
        worktext<-TcltoWin.write(worktext)
        try.res <- try(cat(worktext,file=filename,sep="\n"))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){ cat(paste("file",filename,"saved\n")) }

        if(ok){
          melde(paste("report file",filename,"saved\n"),0)
          melde(paste("w",workname.sys),"cmd.msg")
        } else {
          cat("ERROR: write operation failed!!!\n"); return()
        }





    try(weaveR(workname.sys,show.code=FALSE,replace.umlaute=replace.umlaute.sys))
    cat(sub(".rev$",".tex",workname.sys),"generated\n")
    melde("WeaveReportNoCode",2)
   }
   WeaveReportEchoCode<-function(){
    melde("WeaveReportNoCode",1)
        if(file.exists(filename<-file.path(getwd(),workname.sys))){
          res<-tkmessageBox(message=
                     if(language=="german") 
                        paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                     else paste(filename,"exists. Do you want to replace it?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
          if("externalptr"==mode(res))  res<-tclvalue(res)
          if(res=="cancel")return()
          if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
        }
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
        worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
        worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
        if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
           worktext<-c(worktext,"@\n\\end{document}")
        worktext<-TcltoWin.write(worktext)
        try.res <- try(cat(worktext,file=filename,sep="\n"))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){ cat(paste("file",filename,"saved\n")) }

        if(ok){
          melde(paste("report file",filename,"saved\n"),0)
          melde(paste("w",workname.sys),"cmd.msg")
        } else {
          cat("ERROR: write operation failed!!!\n"); return()
        }





    try(weaveR(workname.sys,show.code="echo",replace.umlaute=replace.umlaute.sys))
    cat(sub(".rev$",".tex",workname.sys),"generated\n")
    melde("WeaveReportNoCode",2)
   }
   WeaveReportNoText<-function(){
    melde("WeaveReportNoText",1)
        if(file.exists(filename<-file.path(getwd(),workname.sys))){
          res<-tkmessageBox(message=
                     if(language=="german") 
                        paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                     else paste(filename,"exists. Do you want to replace it?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
          if("externalptr"==mode(res))  res<-tclvalue(res)
          if(res=="cancel")return()
          if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
        }
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
        worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
        worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
        if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
           worktext<-c(worktext,"@\n\\end{document}")
        worktext<-TcltoWin.write(worktext)
        try.res <- try(cat(worktext,file=filename,sep="\n"))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){ cat(paste("file",filename,"saved\n")) }

        if(ok){
          melde(paste("report file",filename,"saved\n"),0)
          melde(paste("w",workname.sys),"cmd.msg")
        } else {
          cat("ERROR: write operation failed!!!\n"); return()
        }





    try(weaveR(workname.sys,show.text=FALSE,replace.umlaute=replace.umlaute.sys))
    cat(sub(".rev$",".tex",workname.sys),"generated\n")
    melde("WeaveReportNoText",2)
   }
   TangleReport<-function(){
    melde("TangleReport",1)
        if(file.exists(filename<-file.path(getwd(),workname.sys))){
          res<-tkmessageBox(message=
                     if(language=="german") 
                        paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                     else paste(filename,"exists. Do you want to replace it?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
          if("externalptr"==mode(res))  res<-tclvalue(res)
          if(res=="cancel")return()
          if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
        }
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
        worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
        worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
        if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
           worktext<-c(worktext,"@\n\\end{document}")
        worktext<-TcltoWin.write(worktext)
        try.res <- try(cat(worktext,file=filename,sep="\n"))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){ cat(paste("file",filename,"saved\n")) }

        if(ok){
          melde(paste("report file",filename,"saved\n"),0)
          melde(paste("w",workname.sys),"cmd.msg")
        } else {
          cat("ERROR: write operation failed!!!\n"); return()
        }





    try(tangleR(workname.sys))
    cat(sub(".rev$",".R",workname.sys),"generated\n")
    melde("TangleReport",2)
   }
   TangleReportNoComments<-function(){
    melde("TangleReportNoComments",1)
        if(file.exists(filename<-file.path(getwd(),workname.sys))){
          res<-tkmessageBox(message=
                     if(language=="german") 
                        paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                     else paste(filename,"exists. Do you want to replace it?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
          if("externalptr"==mode(res))  res<-tclvalue(res)
          if(res=="cancel")return()
          if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
        }
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
        worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
        worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
        if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
           worktext<-c(worktext,"@\n\\end{document}")
        worktext<-TcltoWin.write(worktext)
        try.res <- try(cat(worktext,file=filename,sep="\n"))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){ cat(paste("file",filename,"saved\n")) }

        if(ok){
          melde(paste("report file",filename,"saved\n"),0)
          melde(paste("w",workname.sys),"cmd.msg")
        } else {
          cat("ERROR: write operation failed!!!\n"); return()
        }





    try(tangleR(workname.sys,insert.comments=FALSE))
    cat(sub(".rev$",".R",workname.sys),"generated\n")
    melde("TangleReportNoComments",2)
   }
   TangleReportChunk<-function(){
    melde("TangleReportChunk",1)
        if(file.exists(filename<-file.path(getwd(),workname.sys))){
          res<-tkmessageBox(message=
                     if(language=="german") 
                        paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                     else paste(filename,"exists. Do you want to replace it?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
          if("externalptr"==mode(res))  res<-tclvalue(res)
          if(res=="cancel")return()
          if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
        }
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
        worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
        worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
        if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
           worktext<-c(worktext,"@\n\\end{document}")
        worktext<-TcltoWin.write(worktext)
        try.res <- try(cat(worktext,file=filename,sep="\n"))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){ cat(paste("file",filename,"saved\n")) }

        if(ok){
          melde(paste("report file",filename,"saved\n"),0)
          melde(paste("w",workname.sys),"cmd.msg")
        } else {
          cat("ERROR: write operation failed!!!\n"); return()
        }





    revfile<-readLines(workname.sys)
    codechunks<-grep(paste("^<","<.*",">",">","=",sep=""),revfile) 
    if(length(codechunks)==0) return()
    codechunks<-revfile[codechunks]; codechunks<-unique(codechunks)
    if(0==(ind<-menu(codechunks))) return()
    TLD<-codechunks[ind]; TLD<-sub(paste("^<","<(.*)",">",">","=",sep=""),"\\1",TLD)
    try(tangleR(workname.sys, expand.roots=TLD,expand.root.start=FALSE))
    cat("Remark: processing without saving of",workname.sys,"\n")
    cat(sub(".rev$",".R",workname.sys),"generated\n")
    melde("TangleReportChunk",2)
   }
   LaTeX.head<-function(){
    melde("LaTeX.head",1)
    # news<-scan(file=paste(path.package("relax"),"lib/LaTeX-head.tex",sep="/"),what="",sep="\n")
    news<-readLines(paste(path.package("relax"),"lib/LaTeX-head.tex",sep="/")) # 2.1.0 # 130325 .path.package defunct
    news<-sub("NOWEBSTYLEFILE",file.path(relax.path,"lib","noweb"),news)
    news<-sub("JOBPATH",getwd(),news); news<-c(news,"")
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    line<-"2.0"
    if(!exists("tworkwin"))
      tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    try(tkinsert(tworkwin,line,paste(news,collapse="\n")))
    melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))


    melde("LaTeX.head",2)
   }
  # LatexReport<-function(){
    LatexReport<-function(...,filename){ ## test???
    melde("LatexReport",1)
    if(!missing(filename)) workname.sys <- filename
    n<-nchar(filename<-workname.sys)
    filename<-sub("rev$","tex",filename)
    if(is.null(n)||5>n||substring(filename,n-3,n)!=".tex"){
      cat("ERROR: file",filename,"not compatible!!!\n");return()
    }
    if(!file.exists(filename)){
      cat("ERROR: file",filename,"not found!!!\n");return()
    }
    melde("latex process starts, messages: see output window\n") 
    # 060116: show logfile
    if((.Platform$OS.type=="windows")){
        latex<-paste(latex.command.windows," ",filename,sep="")
        ok<-try(shell(latex,wait=TRUE))
    }else{ 
        latex<-paste(latex.command.linux," ",filename,sep="")
        ok<-try(system(latex))
    }
    if(ok!=0){ 
          fname<-sub("tex$","log",filename); txt<-readLines(fname)
          ind<-c(grep("(^[!])",txt)[1],grep("^[?].q",txt)[1])
          ind<-ind[!is.na(ind)]
          if(0==length(ind)) ind<-length(txt)-10
          ind<-(-8+min(ind)):(-1+max(ind)); txt<-txt[ind]
          line<-grep("(^l[.][0-9])",txt,value=TRUE)
          line<-line[!is.na(line)]
          if(0<length(line)){
            line<-paste("......\n-> find error line by pressing <Ctrl/Strg> <G>  + line number:",
                        sub("..([0-9]*).*","\\1",line))
          }else{
            line<-"-> Wizardry -> ShowLogFile "
          }
          txt<-c(
                  "!!!ERROR in LaTeX process !!!\n-- for details look at the log file\n",
                  " \nExtract of LaTeX log file:\n......",txt,line
               )
          cat(txt,sep="\n")
    }
    cat("latex process finished:",filename)
    melde("LatexReport",2)
   }
   ShowLogFile<-function(){
    melde("ShowLogFile",1)
    n<-nchar(filename<-workname.sys)
    filename<-sub("rev$","log",filename)
    if(is.null(n)||5>n||substring(filename,n-3,n)!=".log"){
      cat("ERROR: file",filename,"not compatible!!!\n");return()
    }
    if(!file.exists(filename)){ cat("ERROR: file",filename,"not found!!!\n");return() }
    cmd <-paste(editor.sys, filename)
    if((.Platform$OS.type=="windows")){
        try(shell(cmd,wait=FALSE))
    }else{
      if( substring(version$os,1,6)=="darwin"  ){
        cmd<-paste(view.command.mac,filename)
      } 
      try(system(paste(cmd," &")))
    }
    melde("ShowLogFile",2)
   }
   ViewReport<-function(){
    melde("ViewReport",1)
    n<-nchar(filename<-workname.sys)
    filename<-sub("rev$","dvi",filename)
    if(is.null(n)||5>n||substring(filename,n-3,n)!=".dvi"){
      cat("ERROR: file",filename,"not compatible!!!\n");return()
    }
    if(substring(version$os,1,6)=="darwin"  ){
          filename<-sub("dvi$","pdf",filename)   
    }    
    if(!file.exists(filename)){
      cat("ERROR: file",filename,"not found!!!\n");return()
    }
    if((.Platform$OS.type=="windows")){
        view<-paste(view.command.windows," ",filename,sep="")
        try(shell(view,wait=FALSE))
    }else{
        if(substring(version$os,1,5)=="linux") 
          view<-paste(view.command.linux,"  ",filename," &",sep="")
        if(substring(version$os,1,6)=="darwin"  ){
          view<-paste(view.command.mac,"  ",filename," &",sep="")
        }
        try(system(view))
    }
    melde("ViewReport",2)
   }
   DvipdfReport<-function(){  #050607
    melde("DvipdfReport",1)
    n<-nchar(filename<-workname.sys)
    filename<-sub("rev$","dvi",filename)
    if(is.null(n)||5>n||substring(filename,n-3,n)!=".dvi"){
      cat("ERROR: file",filename,"not compatible!!!\n");return()
    }
    if(!file.exists(filename)){
        cat("ERROR: file",filename,"not found!!!\n");return()
    }
    if((.Platform$OS.type=="windows")){
        dvipdf<-paste(dvipdf.command.windows," ",filename,sep="")
        try(shell(dvipdf,wait=FALSE))
    }else{
        dvipdf<-paste(dvipdf.command.linux,"  ",filename," &",sep="")
        try(system(dvipdf))
    }
    cat("\"",dvipdf,"\" has been started!\n")
    melde("DvipdfReport",2)
   }
   ProcessReport<-function(){
          if(file.exists(filename<-file.path(getwd(),workname.sys))){
            res<-tkmessageBox(message=
                       if(language=="german") 
                          paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                       else paste(filename,"exists. Do you want to replace it?"),
                     title="Save File",icon="warning",type="yesnocancel",default="yes")
            if("externalptr"==mode(res))  res<-tclvalue(res)
            if(res=="cancel")return()
            if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
          }
          if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
          tworkwin<-get("tworkwin",envir=revive.sys)
          worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
          if(nchar(worktext)<10000){
            worktext<-strsplit(worktext,"\n")[[1]]
          }else{
            base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
            worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
          }

          worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
          worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
          worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
          if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
             worktext<-c(worktext,"@\n\\end{document}")
          worktext<-TcltoWin.write(worktext)
          try.res <- try(cat(worktext,file=filename,sep="\n"))
          if(is.function(try.res)){
            ok <- "OK"
          } else {
            if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
            if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
            if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
            if(!is.character(ok)) { ok <- "OK" }
          }
          if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
            ok<-FALSE
            error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
            cat(error.msg,"\n")
            if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
               cat("A warning message stopped the evaluation!",
                     "If you want to\nevaluate the code anyway",
                     "evaluate code by:\n>WarnEval<")
            # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
          } else { ok<-TRUE }


          if(ok){ cat(paste("file",filename,"saved\n")) }

          if(ok){
            melde(paste("report file",filename,"saved\n"),0)
            melde(paste("w",workname.sys),"cmd.msg")
          } else {
            cat("ERROR: write operation failed!!!\n"); return()
          }




    
      if(0==length(grep("documentclass", worktext)))
        cat("Warning: there is no \\documentclass-command in the file:",
             "LaTeX needs a preamble for processing tex-files!",sep="\n")
      Sys.sleep(0.25)
      try(weaveR(workname.sys,replace.umlaute=replace.umlaute.sys))
      Sys.sleep(0.5)
      LatexReport()
   }
   ProcessChunk<-function(){
     filename <- "local-chunk.rev"
     if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
     tworkwin<-get("tworkwin",envir=revive.sys)
     worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
     if(nchar(worktext)<10000){
       worktext<-strsplit(worktext,"\n")[[1]]
     }else{
       base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
       worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
     }

     line <-floor(as.numeric(tkindex(tworkwin,"insert")))

     code.start <- grep(paste("^<","<.*>",">=",sep=""),worktext)
     text.start <- grep("^@",worktext)
     doc.start <- grep("^.begin.document.",worktext)
     if(0 == doc.start) { cat("relax warning: no begin document found"); return()}
     if(0 < length(code.start)){  # if 0 no code at all => process complete report
       code.h <- code.start[code.start <= line] 
       if(0 == length(code.h)){ 
         # no code in front of cursor => process from beginning to first code chunk
         code.start <- code.start[1]
         text.h <- c(text.start,length(worktext)+1)
         code.end <- text.h[code.start < text.h][1] - 1
         worktext <- worktext[1:code.end]
       } else {
         cursor.in.code <- TRUE
         text.h <- text.start[text.start <= line] 
         if(0<length(text.h) && code.h[length(code.h)] < text.h[length(text.h)] ){
           # cursor is in text chunk => take text beginning after the end of the last code chunk
           extract.start <- text.h[code.h[length(code.h)] < text.h][1]
           code.h <- code.start[line < code.start] 
           if(0 == length(code.h)){
             # no code after cursor
             extract.end <- length(worktext)
           } else {
             text.h <- text.start[code.h[1] < text.start]
             extract.end <- if(0 == length(text.h)) length(worktext) else text.h[1] - 1
           }
           worktext <- c(worktext[1:doc.start],worktext[extract.start:extract.end])
         } else {
           # cursor is in code chunk => take code end before text chunk up to the end of the code chunk
           if(1 == length(code.h)) extract.start <- doc.start + 1 else {
             code.h <- code.h[length(code.h) - 1]
             extract.start <- text.h[code.h < text.h][1]        
           }
           text.h <- text.start[line < text.start] 
           extract.end <- if(0 == length(text.h)) length(worktext) else text.h[1]-1
           worktext <- c(worktext[1:doc.start],worktext[extract.start:extract.end])
         }
       }
     }
     TcltoWin.write <- get("TcltoWin.write",envir=revive.sys)
     worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
     worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
     worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
     if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
        worktext<-c(worktext,"@\n\\end{document}")
     worktext<-TcltoWin.write(worktext)
     try.res <- try(cat(worktext,file=filename,sep="\n"))
     if(is.function(try.res)){
       ok <- "OK"
     } else {
       if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
       if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
       if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
       if(!is.character(ok)) { ok <- "OK" }
     }
     if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
       ok<-FALSE
       error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
       cat(error.msg,"\n")
       if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
          cat("A warning message stopped the evaluation!",
                "If you want to\nevaluate the code anyway",
                "evaluate code by:\n>WarnEval<")
       # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
     } else { ok<-TRUE }


     if(ok){ cat(paste("file",filename,"saved\n")) }

     Sys.sleep(0.25)
     try(weaveR(filename,replace.umlaute=replace.umlaute.sys))
     Sys.sleep(0.5)
     ## LatexReport()
     LatexReport(filename=filename)
     cat("look at formated file or log file of  'local-chunk'")  
   }

   ProcessBeginEndEnv<-function(){
     filename <- "local-chunk.rev"
     if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
     tworkwin<-get("tworkwin",envir=revive.sys)
     worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
     if(nchar(worktext)<10000){
       worktext<-strsplit(worktext,"\n")[[1]]
     }else{
       base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
       worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
     }

     line <-floor(as.numeric(tkindex(tworkwin,"insert")))

       worktext <- sub("%.*$","",worktext)
       doc.start   <- grep("^.begin.document.",worktext)
       if(0 == length(doc.start)) { cat("relax warning: no begin document found"); return()}
       env.start <- grep("\\\\begin\\{",worktext) ; env.start <- env.start[env.start != doc.start] # } for symmetry
       env.end   <- grep("\\\\end\\{",worktext)  # } for symmetry
       idx.start <- idx.end <- rep(0,length(worktext))
       idx.start[env.start] <- 1; idx.end  [env.end  ] <- 1
       counter <- 0*idx.start; depth <- 0
       for( i in 1:length(counter)){
         counter[i] <- depth <- max(0, depth - idx.end[i] + idx.start[i])
       }
       #  counter <- cumsum(idx.start) - cumsum(idx.end) # fail if there are too much ends
       extract.start <- which(counter == 0 & seq(along=counter) < line)
       extract.start <- if(0 < length(extract.start)) max(extract.start) + 1 else doc.start + 1
       extract.end   <- which(counter == 0 & seq(along=counter) >= line)
       extract.end   <- if(0 < length(extract.start)) min(extract.end)       else length(worktext)  
       worktext <- c(worktext[1:doc.start],worktext[extract.start:extract.end],"\\end{document}")

     worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
     worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
     worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
     if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
        worktext<-c(worktext,"@\n\\end{document}")
     worktext<-TcltoWin.write(worktext)
     try.res <- try(cat(worktext,file=filename,sep="\n"))
     if(is.function(try.res)){
       ok <- "OK"
     } else {
       if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
       if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
       if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
       if(!is.character(ok)) { ok <- "OK" }
     }
     if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
       ok<-FALSE
       error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
       cat(error.msg,"\n")
       if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
          cat("A warning message stopped the evaluation!",
                "If you want to\nevaluate the code anyway",
                "evaluate code by:\n>WarnEval<")
       # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
     } else { ok<-TRUE }


     if(ok){ cat(paste("file",filename,"saved\n")) }

     Sys.sleep(0.25)
     try(weaveR(filename,replace.umlaute=replace.umlaute.sys))
     Sys.sleep(0.5)
     LatexReport(filename=filename)
     cat("look at formated file or log file of  'local-chunk'")  
   }

   SWEAVE<-function(){ ###060112
      workname.Rnw<-sub("rev$","rnw",workname.sys)
      if(file.exists(filename<-file.path(getwd(),workname.sys))){
        res<-tkmessageBox(message=
                   if(language=="german") 
                      paste("Datei",workname.sys,"oder",workname.Rnw,
                                "existiert. Soll(en) sie ersetzt werden?")
                   else paste(filename,"or",workname.Rnw,
                                   "exists. Do you want to replace it/them?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
        if("externalptr"==mode(res))  res<-tclvalue(res)
        if(res=="cancel") return(); if(res=="no") return() 
      }
      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
      worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
      worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
      if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
         worktext<-c(worktext,"@\n\\end{document}")
      worktext<-TcltoWin.write(worktext)
      try.res <- try(cat(worktext,file=filename,sep="\n"))
      if(is.function(try.res)){
        ok <- "OK"
      } else {
        if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
        if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
        if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
        if(!is.character(ok)) { ok <- "OK" }
      }
      if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
        ok<-FALSE
        error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
        cat(error.msg,"\n")
        if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
           cat("A warning message stopped the evaluation!",
                 "If you want to\nevaluate the code anyway",
                 "evaluate code by:\n>WarnEval<")
        # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
      } else { ok<-TRUE }


      if(ok){ cat(paste("file",filename,"saved\n")) }

      Sys.sleep(0.25)
      file.copy(workname.sys,workname.Rnw,overwrite = TRUE)
      frage<-"add additional arguments for Sweave then return:"
      set.tclvalue("tvinfo",sweave.args.sys)
      tkconfigure(linfo.tmp,text=frage)
      # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
      tkpack("forget",linfo); Sys.sleep(0.01)
      tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
      tkfocus(einfo.tmp)
      tkselection.range(einfo.tmp,"0","end") ## 051219
      tkbind(TopW,"<Escape>",function(){
          tkbind(TopW,"<Return>","")
          tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
          # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
          tkpack(linfo,side="left",fill="x",expand="yes")


        }
      )


      tkbind(TopW,"<Return>", function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        set.tclvalue("tvmess","relax")
        sweave.args.sys<-tclvalue("tvinfo")
        assign("sweave.args.sys",sweave.args.sys,envir=revive.sys)
        if(nchar(sweave.args.sys)!=0) {
          sweave.call<-paste("Sweave(\"",workname.Rnw,"\",",sweave.args.sys,")",sep="")
        } else sweave.call<-paste("Sweave(\"",workname.Rnw,"\")",sep="")
        try(eval(parse(text=sweave.call),envir=revive.env))
        Sys.sleep(0.5); LatexReport()
        cat("check R-window for messages!\n")
      } # end of function
    )
   } # end of SWEAVE
  SWEAVEB<-function(){ ###060112
      workname.Rnw<-sub("rev$","rnw",workname.sys)
      if(file.exists(filename<-file.path(getwd(),workname.sys))){
        res<-tkmessageBox(message=
                   if(language=="german") 
                      paste("Datei",workname.sys,"oder",workname.Rnw,
                                "existiert. Soll(en) sie ersetzt werden?")
                   else paste(filename,"or",workname.Rnw,
                                   "exists. Do you want to replace it/them?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
        if("externalptr"==mode(res))  res<-tclvalue(res)
        if(res=="cancel") return(); if(res=="no") return() 
      }
      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
      worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
      worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
      if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
         worktext<-c(worktext,"@\n\\end{document}")
      worktext<-TcltoWin.write(worktext)
      try.res <- try(cat(worktext,file=filename,sep="\n"))
      if(is.function(try.res)){
        ok <- "OK"
      } else {
        if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
        if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
        if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
        if(!is.character(ok)) { ok <- "OK" }
      }
      if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
        ok<-FALSE
        error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
        cat(error.msg,"\n")
        if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
           cat("A warning message stopped the evaluation!",
                 "If you want to\nevaluate the code anyway",
                 "evaluate code by:\n>WarnEval<")
        # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
      } else { ok<-TRUE }


      if(ok){ cat(paste("file",filename,"saved\n")) }

      Sys.sleep(0.25)
      file.copy(workname.sys,workname.Rnw,overwrite = TRUE)
      frage<-"add additional arguments for Sweave then return:"; set.tclvalue("tvinfo",sweave.args.sys)
      tkconfigure(linfo.tmp,text=frage)
      # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
      tkpack("forget",linfo); Sys.sleep(0.01)
      tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
      tkfocus(einfo.tmp)
      tkselection.range(einfo.tmp,"0","end") ## 051219
      tkbind(TopW,"<Escape>",function(){
          tkbind(TopW,"<Return>","")
          tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
          # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
          tkpack(linfo,side="left",fill="x",expand="yes")


        }
      )


      tkbind(TopW,"<Return>", function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        set.tclvalue("tvmess","relax")
        sweave.args.sys<-tclvalue("tvinfo")
        assign("sweave.args.sys",sweave.args.sys,envir=revive.sys)
        if(nchar(sweave.args.sys)!=0) {
          sweave.call<-paste("Sweave(\"",workname.Rnw,"\",",sweave.args.sys,")",sep="")
        } else sweave.call<-paste("Sweave(\"",workname.Rnw,"\")",sep="")
        infile<-tempfile("sweavein"); outfile<-tempfile("sweaveout")
        cat(file=infile,sweave.call,"\n")
     cmd<-paste(R.home(),"/bin/Rcmd BATCH ",infile," ",
                  outfile,sep="")
  if((.Platform$OS.type=="windows")){
     shell(cmd)
  }else{
      system(cmd)
  }
        news<-readLines(outfile)
        if(!exists("toutwin"))
          toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
        pos.to.insert<-"end"
        ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
        news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
        try(tkinsert(toutwin,pos.to.insert,news))
        tksee(toutwin,"end - 0 lines")
        melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

        Sys.sleep(0.5); LatexReport()
        cat("check R-window for messages!\n")
      } # end of function
    )
   } # end of SWEAVEB


  ConstructDemoFunction<-function(){
    melde("ConstructDemoFunction",1)
        if(file.exists(filename<-file.path(getwd(),workname.sys))){
          res<-tkmessageBox(message=
                     if(language=="german") 
                        paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                     else paste(filename,"exists. Do you want to replace it?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
          if("externalptr"==mode(res))  res<-tclvalue(res)
          if(res=="cancel")return()
          if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
        }
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
        worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
        worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
        if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
           worktext<-c(worktext,"@\n\\end{document}")
        worktext<-TcltoWin.write(worktext)
        try.res <- try(cat(worktext,file=filename,sep="\n"))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){ cat(paste("file",filename,"saved\n")) }

        if(ok){
          melde(paste("report file",filename,"saved\n"),0)
          melde(paste("w",workname.sys),"cmd.msg")
        } else {
          cat("ERROR: write operation failed!!!\n"); return()
        }




    
    if(!file.exists(workname.sys)){
      cat(paste("Error: File",workname.sys,"not found!!!"));return()
    }
  #  cat("Remark: processing without saving of",workname.sys,"\n")
    try(tangleR(workname.sys))
    workname<-sub(".rev$",".R",workname.sys)
    cat(workname,"generated\n")
    chunks<-where<-""
    txt<-scan(workname,"",sep="\n")
    fh<-function(file,no=1,start=TRUE){
      # initial code of demo function
      # require("tcltk"); where<-environment() # 140306
    }

    fh<-deparse(fh)
    fh<-c(fh[-length(fh)], "'require(\"tcltk\")'; where<-environment()")
    ft<-function(){
      # end of code of demo function
      ## activate start chunk
      if(start==TRUE){
        no.0<-"0"
        no.start.0<-grep(paste("^#",no.0,":$",sep=""),chunks)
        no.end.0<-grep(paste("^#:",no.0,"$",sep=""),chunks)
        code.0<-chunks[no.start.0:no.end.0]
        eval(parse(text=code.0),envir=where)
      }
      ## activate chunk no
      # eval(parse(text=code),envir=where)
      secno<-tclVar("1") # 0
      show.next.number<-function(...){
       no<-as.character(as.numeric(tclvalue(secno))+1)
       no.start<-grep(paste("^#",no,":$",sep=""),chunks)
       no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
       if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
            is.nan(no.end)||is.nan(no.start)){
            cat("# sorry, chunk number '",no,"' wrong!\n"); return()
       }
       ### cat("# aktueller chunk:",no,"\n")
       code<-paste(chunks[no.start:no.end],collapse="\n")
       if(0<length(code)) {
        tkdelete(ttext,"0.0","end")  
        tkinsert(ttext,"0.0",code)
        tclvalue(secno)<-as.character(no)
        }
      }
      show.back.number<-function(...){
       no<-as.character(as.numeric(tclvalue(secno))-1)
       no.start<-grep(paste("^#",no,":$",sep=""),chunks)
       no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
       if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
            is.nan(no.end)||is.nan(no.start)){
            cat("# sorry, chunk number '",no,"' wrong!\n"); return()
       }
       ### cat("# aktueller chunk:",no,"\n")
       code<-paste(chunks[no.start:no.end],collapse="\n")
       if(0<length(code)) {
        tkdelete(ttext,"0.0","end")  
        tkinsert(ttext,"0.0",code)
        tclvalue(secno)<-as.character(no)
       }
      }
      show.number<-function(...){
       no<-as.character(as.numeric(tclvalue(secno)))
       no.start<-grep(paste("^#",no,":$",sep=""),chunks)
       no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
       if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
            is.nan(no.end)||is.nan(no.start)){
            cat("# sorry, chunk number '",no,"' wrong!\n"); return()
       }
       ### cat("# aktueller chunk:",no,"\n")
       code<-paste(chunks[no.start:no.end],collapse="\n")
       if(0<length(code)) {
        tkdelete(ttext,"0.0","end")  
        tkinsert(ttext,"0.0",code)
        tclvalue(secno)<-as.character(no)
       }
      }
      eval.code<-function(...){
        code<-tclvalue(tkget(ttext,"0.0","end"))
        code.orig<-code<-unlist(strsplit(code,"\n"))
        code<-code[!substring(code,1,1)=="#"]
        ##?## code<-unlist(strsplit(code,";"))  ##110429 Fehler bei sep=";"
        if(length(code)==0){ cat("ok\n"); return() }
        result<-try(eval(parse(text=code),envir=where))
        code.orig<-sub("#([0-9]+):","##wnt-Code-Chunk:\\1-begin#",code.orig)
        code.orig<-sub("#:([0-9]+)","##wnt-Code-Chunk:\\1-end#",code.orig)
        h<-get("allcodechunks",envir=where)
        h<-c(h,paste("<","<*>",">=",sep=""),code.orig,"\n@\n")
        assign("allcodechunks",h,envir=where)
        code<-sub("^ *","",code)
        code<-code[nchar(code)>0]
        lexpr<-rev(code)[1]; lexpr<-substring(lexpr,1,4)
        if(length(code)==0||is.null(lexpr)||is.na(lexpr)) return()
        plot.res<-c("plot","boxp","par(","abli","pie(","hist","axis","show",
               "lsfi","pair","ylab","help",
               "qqli","qqno","qqpl","rug(","lege","segm","text","xlab", 
               "poin","line","titl","eda(","imag","vgl.","curv")
        if(any(plot.res==lexpr)){
           cat("Plot erstellt\n"); return()
        }
        if(is.null(result)||is.na(result)||lexpr=="prin"||lexpr=="cat("){
          cat("ok\n"); return() }
        if(is.list(result)&& length(names(result))> 0 && 
                                   names(result)[1]=="ID") return()
        ## if(is.list(result)&& TRUE) return()
        no<-as.character(as.numeric(tclvalue(secno)))
        cat("Result of code chunk",no,":\n") 
       if(class(result)=="try-error"){
         class(result)<-"character"
         cat(result,"\n")
       }else{
        print(result)
       }
       cat("ok\n")
       }
      exit.function<-function(){
         tkdestroy(top)
         filename<-tkgetSaveFile(filetypes="{{Paper Files} {.rev}}",
                                               title="Do you want to save the activated R statements?")
         if(!is.character(filename)) filename<-tclvalue(filename)
         if(filename==""){
           cat("Demo function stopped without saving\n")
           return()
         }
         if(0==length(grep("rev$",filename))) filename<-paste(filename,".rev",sep="")
         h<-get("allcodechunks",envir=where)
         try(cat(h,sep="\n",file=filename))
         cat(paste("Remark: activated statements saved in\n   ",filename,"\n"))
         return()
      }
      allcodechunks<-paste(
              "@\nReport of activated R-chunks from: ",date(),
              "\n(demo function constructed by relax (c) Peter Wolf 2007)\n\n  ", sep="")
      no<-0
      no.start<-grep(paste("^#",no,":$",sep=""),chunks)
      no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
      if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
           is.nan(no.end)||is.nan(no.start)){
           cat("# sorry, chunk number '",no,"' wrong!\n"); return()
      }
      ### cat("# aktueller chunk:",no,"\n")
      code<-paste(chunks[no.start:no.end],collapse="\n")
      h<-paste(rep("#",60),collapse="")
      code<-sub("#0:",h,code); code<-sub("#:0",h,code)
      allcodechunks<-c(allcodechunks,"\n@\n<<start>>=",code,"\n@")
      assign("allcodechunks",allcodechunks,envir=where)
      top<-tktoplevel()
      ttext<-tktext(top,height=19,background="#f7fffF",
                    font="-Adobe-courier-Medium-R-Normal--18-180-*")
      tf<-tkframe(top) 
      tkwm.title(top, "demo of file WoRkNaMe, constructed by relax (c) Peter Wolf 2011")
      tkpack(tf,side="bottom"); tkpack(tf,ttext,side="bottom",fill="both",expand="y") # 091026
      tkevent.add("<<Paste>>",   "<Control_L><v>")
    ## ok:  tkbind(ttext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
    if(substring(version$os,1,6)=="darwin"  ){
      mac.paste<-function(...){  ## Ctrl-V  # mac-PASTE
         try({.Tcl("clipboard append hello"); .Tcl("clipboard clear")
              news<-base::scan(file=pipe("pbpaste","r"),what="",sep="\n",blank.lines.skip=FALSE)
              tkinsert(ttext,"insert",paste(news,collapse="\n"))})
              tksee(ttext,"insert - 7 lines");  tksee(ttext,"insert + 7 lines") #090706 
      } 
      tkbind(ttext,"<Control_L><v>",mac.paste)
      mac.copy<-function(...){  ## Ctrl-C  # mac-COPY
         news<-""
         try(news<-tclvalue(.Tcl("if {[catch {clipboard get}]} {set aa empty} {set aa full}")))
         if(news=="empty") return()
         try({news<-tclvalue(.Tcl("set aaa [selection get -selection CLIPBOARD]"))
              tmp.file.name <- tempfile("rt-tmp")
              base::cat(news,file=tmp.file.name); system(paste("pbcopy < ",tmp.file.name))
              .Tcl("clipboard append hello"); .Tcl("clipboard clear")})
      }
      tkbind(ttext,"<Control_L><c>",mac.copy)
      tkevent.add("<<extract>>",  "<Control_L><c><KeyRelease>")   # mac-extract
      tkbind(ttext,"<<extract>>",mac.copy)
    }else{
      tkevent.add("<<Paste>>",   "<Control_L><v>")
      tkbind(ttext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
    }


      bexit<-tkbutton(tf,text="QUIT",width=9)
      beval<-tkbutton(tf,text="EVALUATE",width=9)
      bnext<-tkbutton(tf,text="NEXT",width=9)
      bback<-tkbutton(tf,text="BACK",width=9)
      tkbind(top, "<<EvalRCode>>", eval.code)
      if( ! substring(version$os,1,6)=="darwin"  ) {
        tkbind(top, "<<Next>>", show.next.number)
        tkbind(top, "<<Back>>", show.back.number)
      }
      lno  <-tkentry(tf,textvariable=secno,width=9)
      linfo<-tklabel(tf,text="chunk number:")
      tkpack(linfo,lno,beval,bnext,bback,bexit,side="left")
      tkconfigure(bexit,command=exit.function)
      tkconfigure(bnext,command=show.next.number)
      tkconfigure(bback,command=show.back.number)
      tkconfigure(beval,command=eval.code)
      # tkbind(lno,"<Return>",show.number)
      tkbind(lno,"<KeyRelease>",show.number)
      tclvalue(secno)<-as.character(no)
      show.number()
      ### tkwait.window(top)
    }

    ft<-deparse(ft)[-(1:2)]
    ft<-sub("WoRkNaMe",workname,ft)
    fname<-gsub("[^A-Za-z]","",sub(".R$",".demo",workname))
    txt<-gsub("\\\\","\\\\\\\\",txt)
    txt<-gsub("\"","\\\\\"",txt)
    txt<-paste(',"',txt,'"',sep=""); txt[1]<-substring(txt[1],2)
    txt<-c(paste(fname,"<-"),fh,"chunks<-c(",txt,")",ft
           ,paste("cat(\"Demo will be started by > ",fname,"()\\n\")",sep="")
           ,paste(fname,"()\n",sep="")
        )
    workname<-sub(".R$",".demo.R",workname)
    cat(txt,file=workname,sep="\n")
    cat("Demo in file ",workname," saved, load demo by source(\"",workname,"\") \n",sep="")
    melde("ConstructDemoFunction",2)
  }

  StartCodeChunkPlayer<-function(){
    melde("StartCodeChunkPlayer",1)
    cat("Code-Chunk-Player saves actual version of",workname.sys,"\n") 
        if(file.exists(filename<-file.path(getwd(),workname.sys))){
          res<-tkmessageBox(message=
                     if(language=="german") 
                        paste("Datei",filename,"existiert. Soll sie ersetzt werden?")
                     else paste(filename,"exists. Do you want to replace it?"),
                   title="Save File",icon="warning",type="yesnocancel",default="yes")
          if("externalptr"==mode(res))  res<-tclvalue(res)
          if(res=="cancel")return()
          if(res=="no"){msg<-SaveReport(); if(msg=="cancel") return()}  #050607
        }
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
        worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
        worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
        if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
           worktext<-c(worktext,"@\n\\end{document}")
        worktext<-TcltoWin.write(worktext)
        try.res <- try(cat(worktext,file=filename,sep="\n"))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){ cat(paste("file",filename,"saved\n")) }

        if(ok){
          melde(paste("report file",filename,"saved\n"),0)
          melde(paste("w",workname.sys),"cmd.msg")
        } else {
          cat("ERROR: write operation failed!!!\n"); return()
        }





    try(tangleR(workname.sys)); workname<-sub(".rev$",".R",workname.sys); cat(workname,"generated\n")
    chunks<-where<-""; txt<-scan(workname,"",sep="\n")
    # construct structure of function
    fh<-function(file,no=1,start=TRUE){
      # initial code of demo function
      # require("tcltk"); where<-environment() # 140306
    }

    fh<-deparse(fh)
    fh<-c(fh[-length(fh)], "'require(\"tcltk\")'; where<-environment()")
    ft<-function(){
      # end of code of demo function
      ## activate start chunk
      if(start==TRUE){
        no.0<-"0"
        no.start.0<-grep(paste("^#",no.0,":$",sep=""),chunks)
        no.end.0<-grep(paste("^#:",no.0,"$",sep=""),chunks)
        code.0<-chunks[no.start.0:no.end.0]
        eval(parse(text=code.0),envir=where)
      }
      ## activate chunk no
      # eval(parse(text=code),envir=where)
      secno<-tclVar("1") # 0
      show.next.number<-function(...){
       no<-as.character(as.numeric(tclvalue(secno))+1)
       no.start<-grep(paste("^#",no,":$",sep=""),chunks)
       no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
       if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
            is.nan(no.end)||is.nan(no.start)){
            cat("# sorry, chunk number '",no,"' wrong!\n"); return()
       }
       ### cat("# aktueller chunk:",no,"\n")
       code<-paste(chunks[no.start:no.end],collapse="\n")
       if(0<length(code)) {
        tkdelete(ttext,"0.0","end")  
        tkinsert(ttext,"0.0",code)
        tclvalue(secno)<-as.character(no)
        }
      }
      show.back.number<-function(...){
       no<-as.character(as.numeric(tclvalue(secno))-1)
       no.start<-grep(paste("^#",no,":$",sep=""),chunks)
       no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
       if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
            is.nan(no.end)||is.nan(no.start)){
            cat("# sorry, chunk number '",no,"' wrong!\n"); return()
       }
       ### cat("# aktueller chunk:",no,"\n")
       code<-paste(chunks[no.start:no.end],collapse="\n")
       if(0<length(code)) {
        tkdelete(ttext,"0.0","end")  
        tkinsert(ttext,"0.0",code)
        tclvalue(secno)<-as.character(no)
       }
      }
      show.number<-function(...){
       no<-as.character(as.numeric(tclvalue(secno)))
       no.start<-grep(paste("^#",no,":$",sep=""),chunks)
       no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
       if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
            is.nan(no.end)||is.nan(no.start)){
            cat("# sorry, chunk number '",no,"' wrong!\n"); return()
       }
       ### cat("# aktueller chunk:",no,"\n")
       code<-paste(chunks[no.start:no.end],collapse="\n")
       if(0<length(code)) {
        tkdelete(ttext,"0.0","end")  
        tkinsert(ttext,"0.0",code)
        tclvalue(secno)<-as.character(no)
       }
      }
      eval.code<-function(...){
        code<-tclvalue(tkget(ttext,"0.0","end"))
        code.orig<-code<-unlist(strsplit(code,"\n"))
        code<-code[!substring(code,1,1)=="#"]
        ##?## code<-unlist(strsplit(code,";"))  ##110429 Fehler bei sep=";"
        if(length(code)==0){ cat("ok\n"); return() }
        result<-try(eval(parse(text=code),envir=where))
        code.orig<-sub("#([0-9]+):","##wnt-Code-Chunk:\\1-begin#",code.orig)
        code.orig<-sub("#:([0-9]+)","##wnt-Code-Chunk:\\1-end#",code.orig)
        h<-get("allcodechunks",envir=where)
        h<-c(h,paste("<","<*>",">=",sep=""),code.orig,"\n@\n")
        assign("allcodechunks",h,envir=where)
        code<-sub("^ *","",code)
        code<-code[nchar(code)>0]
        lexpr<-rev(code)[1]; lexpr<-substring(lexpr,1,4)
        if(length(code)==0||is.null(lexpr)||is.na(lexpr)) return()
        plot.res<-c("plot","boxp","par(","abli","pie(","hist","axis","show",
               "lsfi","pair","ylab","help",
               "qqli","qqno","qqpl","rug(","lege","segm","text","xlab", 
               "poin","line","titl","eda(","imag","vgl.","curv")
        if(any(plot.res==lexpr)){
           cat("Plot erstellt\n"); return()
        }
        if(is.null(result)||is.na(result)||lexpr=="prin"||lexpr=="cat("){
          cat("ok\n"); return() }
        if(is.list(result)&& length(names(result))> 0 && 
                                   names(result)[1]=="ID") return()
        ## if(is.list(result)&& TRUE) return()
        no<-as.character(as.numeric(tclvalue(secno)))
        cat("Result of code chunk",no,":\n") 
       if(class(result)=="try-error"){
         class(result)<-"character"
         cat(result,"\n")
       }else{
        print(result)
       }
       cat("ok\n")
       }
      exit.function<-function(){
         tkdestroy(top)
         filename<-tkgetSaveFile(filetypes="{{Paper Files} {.rev}}",
                                               title="Do you want to save the activated R statements?")
         if(!is.character(filename)) filename<-tclvalue(filename)
         if(filename==""){
           cat("Demo function stopped without saving\n")
           return()
         }
         if(0==length(grep("rev$",filename))) filename<-paste(filename,".rev",sep="")
         h<-get("allcodechunks",envir=where)
         try(cat(h,sep="\n",file=filename))
         cat(paste("Remark: activated statements saved in\n   ",filename,"\n"))
         return()
      }
      allcodechunks<-paste(
              "@\nReport of activated R-chunks from: ",date(),
              "\n(demo function constructed by relax (c) Peter Wolf 2007)\n\n  ", sep="")
      no<-0
      no.start<-grep(paste("^#",no,":$",sep=""),chunks)
      no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
      if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
           is.nan(no.end)||is.nan(no.start)){
           cat("# sorry, chunk number '",no,"' wrong!\n"); return()
      }
      ### cat("# aktueller chunk:",no,"\n")
      code<-paste(chunks[no.start:no.end],collapse="\n")
      h<-paste(rep("#",60),collapse="")
      code<-sub("#0:",h,code); code<-sub("#:0",h,code)
      allcodechunks<-c(allcodechunks,"\n@\n<<start>>=",code,"\n@")
      assign("allcodechunks",allcodechunks,envir=where)
      top<-tktoplevel()
      ttext<-tktext(top,height=19,background="#f7fffF",
                    font="-Adobe-courier-Medium-R-Normal--18-180-*")
      tf<-tkframe(top) 
      tkwm.title(top, "demo of file WoRkNaMe, constructed by relax (c) Peter Wolf 2011")
      tkpack(tf,side="bottom"); tkpack(tf,ttext,side="bottom",fill="both",expand="y") # 091026
      tkevent.add("<<Paste>>",   "<Control_L><v>")
    ## ok:  tkbind(ttext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
    if(substring(version$os,1,6)=="darwin"  ){
      mac.paste<-function(...){  ## Ctrl-V  # mac-PASTE
         try({.Tcl("clipboard append hello"); .Tcl("clipboard clear")
              news<-base::scan(file=pipe("pbpaste","r"),what="",sep="\n",blank.lines.skip=FALSE)
              tkinsert(ttext,"insert",paste(news,collapse="\n"))})
              tksee(ttext,"insert - 7 lines");  tksee(ttext,"insert + 7 lines") #090706 
      } 
      tkbind(ttext,"<Control_L><v>",mac.paste)
      mac.copy<-function(...){  ## Ctrl-C  # mac-COPY
         news<-""
         try(news<-tclvalue(.Tcl("if {[catch {clipboard get}]} {set aa empty} {set aa full}")))
         if(news=="empty") return()
         try({news<-tclvalue(.Tcl("set aaa [selection get -selection CLIPBOARD]"))
              tmp.file.name <- tempfile("rt-tmp")
              base::cat(news,file=tmp.file.name); system(paste("pbcopy < ",tmp.file.name))
              .Tcl("clipboard append hello"); .Tcl("clipboard clear")})
      }
      tkbind(ttext,"<Control_L><c>",mac.copy)
      tkevent.add("<<extract>>",  "<Control_L><c><KeyRelease>")   # mac-extract
      tkbind(ttext,"<<extract>>",mac.copy)
    }else{
      tkevent.add("<<Paste>>",   "<Control_L><v>")
      tkbind(ttext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
    }


      bexit<-tkbutton(tf,text="QUIT",width=9)
      beval<-tkbutton(tf,text="EVALUATE",width=9)
      bnext<-tkbutton(tf,text="NEXT",width=9)
      bback<-tkbutton(tf,text="BACK",width=9)
      tkbind(top, "<<EvalRCode>>", eval.code)
      if( ! substring(version$os,1,6)=="darwin"  ) {
        tkbind(top, "<<Next>>", show.next.number)
        tkbind(top, "<<Back>>", show.back.number)
      }
      lno  <-tkentry(tf,textvariable=secno,width=9)
      linfo<-tklabel(tf,text="chunk number:")
      tkpack(linfo,lno,beval,bnext,bback,bexit,side="left")
      tkconfigure(bexit,command=exit.function)
      tkconfigure(bnext,command=show.next.number)
      tkconfigure(bback,command=show.back.number)
      tkconfigure(beval,command=eval.code)
      # tkbind(lno,"<Return>",show.number)
      tkbind(lno,"<KeyRelease>",show.number)
      tclvalue(secno)<-as.character(no)
      show.number()
      ### tkwait.window(top)
    }

    ft<-deparse(ft)[-(1:2)]; ft<-sub("WoRkNaMe",workname,ft)
    # modify function for StartCodeChunkPlayer
    txt<-c("CCPlayer<-",fh,"    chunks<-readLines(file)",ft,"CCPlayer()")
    # set file argument
    txt[2] <- sub("file",paste("file = '",workname,"'",sep=""),txt[2]) 
    # use argument no as start chunk number
    idx <- grep("^    no .. 0",txt); txt[idx] <- sub("0","no #set start chunk",txt[idx])
    # modify title
    idx <- grep("tkwm.title",txt); txt[idx] <- sub("Peter Wolf.","",txt[idx])
    txt[idx] <- sub("demo of file","Code-Chunk-Player, file:",txt[idx])
    # modify exit message
    idx <- grep("Demo function stopped",txt); 
    txt[idx] <- sub("Demo function stopped without","Code-Chunk-Player stopped without",txt[idx])
    # set evaluation should take place in revive.env
    idx <- grep("environment",txt); txt[idx] <- sub("environment..","revive.env",txt[idx])
    # print viewer function for debugging
    melde(paste(txt,collapse="\n"),2)
    # start of history viewer
    try(eval(parse(text=txt),envir=revive.env))
    melde("StartCodeChunkPlayer",2)
  }

  assign("relax.history",
         c(if(0<length(ls(pattern="relax.history",envir=revive.sys))) # to extend existing history
           get("relax.history",envir=revive.sys) else NULL,           # to extend existing history
           list(c("@",paste("<","<*>",">=",sep=""),paste("# relax history file --",date())))
         )                                                           # to extend existing history
         ,envir=revive.sys)
  ShowHistory<-function(){
    melde("ShowHistory",1)
      # check history
      if(0==length(grep("relax.history",ls(envir=revive.sys)))){ 
        cat("relax warning: no history found"); return() 
      }
      # save history in file
      relax.history<-get("relax.history",envir=revive.sys); txt<-unlist(relax.history)
      workname<-sub(".rev$",".history.rev",get("workname.sys",envir=revive.sys))
      base::cat(txt,file=workname,sep="\n")
      # tangle history file
      try(tangleR(workname)); workname<-sub(".rev$",".R",workname)
      cat(sub(".R$",".rev",workname),"/",workname,"generated\n")

    # construct history viewer
    chunks<-where<-"" # to reduce warnings during package construction
    fh<-function(file,no=1,start=TRUE){
      # initial code of demo function
      # require("tcltk"); where<-environment() # 140306
    }

    fh<-deparse(fh)
    fh<-c(fh[-length(fh)], "'require(\"tcltk\")'; where<-environment()")
    ft<-function(){
      # end of code of demo function
      ## activate start chunk
      if(start==TRUE){
        no.0<-"0"
        no.start.0<-grep(paste("^#",no.0,":$",sep=""),chunks)
        no.end.0<-grep(paste("^#:",no.0,"$",sep=""),chunks)
        code.0<-chunks[no.start.0:no.end.0]
        eval(parse(text=code.0),envir=where)
      }
      ## activate chunk no
      # eval(parse(text=code),envir=where)
      secno<-tclVar("1") # 0
      show.next.number<-function(...){
       no<-as.character(as.numeric(tclvalue(secno))+1)
       no.start<-grep(paste("^#",no,":$",sep=""),chunks)
       no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
       if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
            is.nan(no.end)||is.nan(no.start)){
            cat("# sorry, chunk number '",no,"' wrong!\n"); return()
       }
       ### cat("# aktueller chunk:",no,"\n")
       code<-paste(chunks[no.start:no.end],collapse="\n")
       if(0<length(code)) {
        tkdelete(ttext,"0.0","end")  
        tkinsert(ttext,"0.0",code)
        tclvalue(secno)<-as.character(no)
        }
      }
      show.back.number<-function(...){
       no<-as.character(as.numeric(tclvalue(secno))-1)
       no.start<-grep(paste("^#",no,":$",sep=""),chunks)
       no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
       if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
            is.nan(no.end)||is.nan(no.start)){
            cat("# sorry, chunk number '",no,"' wrong!\n"); return()
       }
       ### cat("# aktueller chunk:",no,"\n")
       code<-paste(chunks[no.start:no.end],collapse="\n")
       if(0<length(code)) {
        tkdelete(ttext,"0.0","end")  
        tkinsert(ttext,"0.0",code)
        tclvalue(secno)<-as.character(no)
       }
      }
      show.number<-function(...){
       no<-as.character(as.numeric(tclvalue(secno)))
       no.start<-grep(paste("^#",no,":$",sep=""),chunks)
       no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
       if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
            is.nan(no.end)||is.nan(no.start)){
            cat("# sorry, chunk number '",no,"' wrong!\n"); return()
       }
       ### cat("# aktueller chunk:",no,"\n")
       code<-paste(chunks[no.start:no.end],collapse="\n")
       if(0<length(code)) {
        tkdelete(ttext,"0.0","end")  
        tkinsert(ttext,"0.0",code)
        tclvalue(secno)<-as.character(no)
       }
      }
      eval.code<-function(...){
        code<-tclvalue(tkget(ttext,"0.0","end"))
        code.orig<-code<-unlist(strsplit(code,"\n"))
        code<-code[!substring(code,1,1)=="#"]
        ##?## code<-unlist(strsplit(code,";"))  ##110429 Fehler bei sep=";"
        if(length(code)==0){ cat("ok\n"); return() }
        result<-try(eval(parse(text=code),envir=where))
        code.orig<-sub("#([0-9]+):","##wnt-Code-Chunk:\\1-begin#",code.orig)
        code.orig<-sub("#:([0-9]+)","##wnt-Code-Chunk:\\1-end#",code.orig)
        h<-get("allcodechunks",envir=where)
        h<-c(h,paste("<","<*>",">=",sep=""),code.orig,"\n@\n")
        assign("allcodechunks",h,envir=where)
        code<-sub("^ *","",code)
        code<-code[nchar(code)>0]
        lexpr<-rev(code)[1]; lexpr<-substring(lexpr,1,4)
        if(length(code)==0||is.null(lexpr)||is.na(lexpr)) return()
        plot.res<-c("plot","boxp","par(","abli","pie(","hist","axis","show",
               "lsfi","pair","ylab","help",
               "qqli","qqno","qqpl","rug(","lege","segm","text","xlab", 
               "poin","line","titl","eda(","imag","vgl.","curv")
        if(any(plot.res==lexpr)){
           cat("Plot erstellt\n"); return()
        }
        if(is.null(result)||is.na(result)||lexpr=="prin"||lexpr=="cat("){
          cat("ok\n"); return() }
        if(is.list(result)&& length(names(result))> 0 && 
                                   names(result)[1]=="ID") return()
        ## if(is.list(result)&& TRUE) return()
        no<-as.character(as.numeric(tclvalue(secno)))
        cat("Result of code chunk",no,":\n") 
       if(class(result)=="try-error"){
         class(result)<-"character"
         cat(result,"\n")
       }else{
        print(result)
       }
       cat("ok\n")
       }
      exit.function<-function(){
         tkdestroy(top)
         filename<-tkgetSaveFile(filetypes="{{Paper Files} {.rev}}",
                                               title="Do you want to save the activated R statements?")
         if(!is.character(filename)) filename<-tclvalue(filename)
         if(filename==""){
           cat("Demo function stopped without saving\n")
           return()
         }
         if(0==length(grep("rev$",filename))) filename<-paste(filename,".rev",sep="")
         h<-get("allcodechunks",envir=where)
         try(cat(h,sep="\n",file=filename))
         cat(paste("Remark: activated statements saved in\n   ",filename,"\n"))
         return()
      }
      allcodechunks<-paste(
              "@\nReport of activated R-chunks from: ",date(),
              "\n(demo function constructed by relax (c) Peter Wolf 2007)\n\n  ", sep="")
      no<-0
      no.start<-grep(paste("^#",no,":$",sep=""),chunks)
      no.end<-grep(paste("^#:",no,"$",sep=""),chunks)
      if(length(no.end)==0||is.na(no.end) ||is.na(no.start)||
           is.nan(no.end)||is.nan(no.start)){
           cat("# sorry, chunk number '",no,"' wrong!\n"); return()
      }
      ### cat("# aktueller chunk:",no,"\n")
      code<-paste(chunks[no.start:no.end],collapse="\n")
      h<-paste(rep("#",60),collapse="")
      code<-sub("#0:",h,code); code<-sub("#:0",h,code)
      allcodechunks<-c(allcodechunks,"\n@\n<<start>>=",code,"\n@")
      assign("allcodechunks",allcodechunks,envir=where)
      top<-tktoplevel()
      ttext<-tktext(top,height=19,background="#f7fffF",
                    font="-Adobe-courier-Medium-R-Normal--18-180-*")
      tf<-tkframe(top) 
      tkwm.title(top, "demo of file WoRkNaMe, constructed by relax (c) Peter Wolf 2011")
      tkpack(tf,side="bottom"); tkpack(tf,ttext,side="bottom",fill="both",expand="y") # 091026
      tkevent.add("<<Paste>>",   "<Control_L><v>")
    ## ok:  tkbind(ttext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
    if(substring(version$os,1,6)=="darwin"  ){
      mac.paste<-function(...){  ## Ctrl-V  # mac-PASTE
         try({.Tcl("clipboard append hello"); .Tcl("clipboard clear")
              news<-base::scan(file=pipe("pbpaste","r"),what="",sep="\n",blank.lines.skip=FALSE)
              tkinsert(ttext,"insert",paste(news,collapse="\n"))})
              tksee(ttext,"insert - 7 lines");  tksee(ttext,"insert + 7 lines") #090706 
      } 
      tkbind(ttext,"<Control_L><v>",mac.paste)
      mac.copy<-function(...){  ## Ctrl-C  # mac-COPY
         news<-""
         try(news<-tclvalue(.Tcl("if {[catch {clipboard get}]} {set aa empty} {set aa full}")))
         if(news=="empty") return()
         try({news<-tclvalue(.Tcl("set aaa [selection get -selection CLIPBOARD]"))
              tmp.file.name <- tempfile("rt-tmp")
              base::cat(news,file=tmp.file.name); system(paste("pbcopy < ",tmp.file.name))
              .Tcl("clipboard append hello"); .Tcl("clipboard clear")})
      }
      tkbind(ttext,"<Control_L><c>",mac.copy)
      tkevent.add("<<extract>>",  "<Control_L><c><KeyRelease>")   # mac-extract
      tkbind(ttext,"<<extract>>",mac.copy)
    }else{
      tkevent.add("<<Paste>>",   "<Control_L><v>")
      tkbind(ttext,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
    }


      bexit<-tkbutton(tf,text="QUIT",width=9)
      beval<-tkbutton(tf,text="EVALUATE",width=9)
      bnext<-tkbutton(tf,text="NEXT",width=9)
      bback<-tkbutton(tf,text="BACK",width=9)
      tkbind(top, "<<EvalRCode>>", eval.code)
      if( ! substring(version$os,1,6)=="darwin"  ) {
        tkbind(top, "<<Next>>", show.next.number)
        tkbind(top, "<<Back>>", show.back.number)
      }
      lno  <-tkentry(tf,textvariable=secno,width=9)
      linfo<-tklabel(tf,text="chunk number:")
      tkpack(linfo,lno,beval,bnext,bback,bexit,side="left")
      tkconfigure(bexit,command=exit.function)
      tkconfigure(bnext,command=show.next.number)
      tkconfigure(bback,command=show.back.number)
      tkconfigure(beval,command=eval.code)
      # tkbind(lno,"<Return>",show.number)
      tkbind(lno,"<KeyRelease>",show.number)
      tclvalue(secno)<-as.character(no)
      show.number()
      ### tkwait.window(top)
    }

    ft<-deparse(ft)[-(1:2)]; ft<-sub("WoRkNaMe",workname,ft)
    txt<-c("relax.hist<-",fh,"    chunks<-readLines(file)",ft,"relax.hist()")
    # adapt arguments
    txt[2] <- sub("TRUE","FALSE",txt[2])
    txt[2] <- sub("[1]",length(relax.history),txt[2])
    txt[2] <- sub("file",paste("file = '",workname,"'",sep=""),txt[2])
    # starting chunk: newest one
    idx <- grep("^    no .. 0",txt); txt[idx] <- sub("0","no #set start chunk",txt[idx])
    # modify title
    idx <- grep("tkwm.title",txt); txt[idx] <- sub("demo of file","history saved in file",txt[idx])
    # simplify exit.function of viewer
    idx <- grep("exit.function.*function",txt)
    txt[idx] <- paste(txt[idx]," tkdestroy(top);return() }; exit.demo <- function(){")
    # evaluation should take place in revive.env
    idx <- grep("environment",txt); txt[idx] <- sub("environment..","revive.env",txt[idx])
    # print viewer function during debugging
    melde(paste(txt,collapse="\n"),2)
    # start of history viewer
    try(eval(parse(text=txt),envir=revive.env))
    melde("ShowHistory",2)
  }

  SaveHistory<-function(){
    melde("SaveHistory",1)
      # check history
      if(0==length(grep("relax.history",ls(envir=revive.sys)))){ 
        cat("relax warning: no history found"); return() 
      }
      # save history in file
      relax.history<-get("relax.history",envir=revive.sys); txt<-unlist(relax.history)
      workname<-sub(".rev$",".history.rev",get("workname.sys",envir=revive.sys))
      base::cat(txt,file=workname,sep="\n")
      # tangle history file
      try(tangleR(workname)); workname<-sub(".rev$",".R",workname)
      cat(sub(".R$",".rev",workname),"/",workname,"generated\n")

    melde("SaveHistory",2)
  }

  FindReportChunk<-function(){
    melde("FindReportChunk",1)
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    no <- grep("^<<(.*)>>=",worktext); if(0==length(no)) return()

    cchoices<-paste(worktext[no],worktext[no+1],sep=":: ")
    cchoices<-substring(cchoices,1,pmin(nchar(cchoices),60))
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    actual.chunk<-min(max(1,c(which(line<no),length(no)+1)[1]-1),length(no))

    newtop<-tktoplevel();tkwm.title(newtop,"code chunk? -- Esc=Exit, Return=Selection")
    scr <- tkscrollbar(newtop, command=function(...)tkyview(tl,...))
    tl<-tklistbox(newtop,height=min(20,length(no)),width=60,selectmode="single",
                yscrollcommand=function(...)tkset(scr,...),background="white")
    for(ch in cchoices) tkinsert(tl,"end",ch)
    tkselection.set(tl,actual.chunk-1)  # Default
    tksee(tl,max(actual.chunk-10,0))
    tkpack(tl,side="left",expand="yes",fill="y"); tkpack(scr,side="left",expand="yes",fill="y")
    tkbind(newtop,"<Escape>",function()tkdestroy(newtop))
    tkbind(newtop,"<Return>",function(){
       h<-as.numeric(tkcurselection(tl))+1; cchoice<-cchoices[h]; 
       cchoice<-sub("::.*$","",cchoice)
       if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
       tworkwin<-get("tworkwin",envir=revive.sys)
       worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
       if(nchar(worktext)<10000){
         worktext<-strsplit(worktext,"\n")[[1]]
       }else{
         base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
         worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
       }

       line<-which(cchoice==substring(worktext,1,nchar(cchoice)))[1]
       # line<-no[as.numeric(tkcurselection(tl))+1]; tkdestroy(newtop)
       if(!is.na(line)){
         tksee(tworkwin,paste(line,".1",sep=""))
         tkmark.set(tworkwin, "insert", paste(line,".1",sep=""))
         tkfocus(tworkwin)
       }
    })
    melde("FindReportChunk",2)
  }
  InsertLaTeXEnv<-function(){
    melde("InsertLaTeXEnv",1)
    newtop<-tktoplevel();tkwm.title(newtop,"LaTeX-environment? -- Esc=Exit, Return=Selection")
    tkwm.geometry(newtop,"+20+15")
    scr <- tkscrollbar(newtop, command=function(...)tkyview(tl,...))
    choices<-c("center","quote","itemize","enumerate","eqnarray*","verbatim",
                "\\[..\\]", "\\section", "\\subsection", "\\subsubsection", "\\paragraph",
                "\\item","\\includegraphics","\\emph","\\textbf","\\texttt")
    tl<-tklistbox(newtop,height=length(choices),width=20,selectmode="single",
                yscrollcommand=function(...)tkset(scr,...),background="white")
    for(ch in choices) tkinsert(tl,"end",ch)
    tkpack(tl,side="left",expand="yes",fill="y"); tkpack(scr,side="left",expand="yes",fill="y")
    tkbind(newtop,"<Escape>",function()tkdestroy(newtop))
    tkbind(newtop,"<Return>",function(){
       choice<-as.numeric(tkcurselection(tl))+1; # tkdestroy(newtop)
       if(!is.na(choice) && any(choice==1:16)){
         news<-c("\n\\begin{center}\n\n\\end{center}",
                 "\n\\begin{quote}\n\n\\end{quote}",
                 "\n\\begin{itemize}\n\\item \n\\item \n\n\\end{itemize}",
                 "\n\\begin{enumerate}\n\\item \n\\item \n\n\\end{enumerate}",
                 "\n\\begin{eqnarray*}\n\n\\end{eqnarray*}",
                 "\n\\begin{verbatim}\n\n\\end{verbatim}",
                 "\n\\[\n\n\\]",
                 "\n@\n\\section{}",
                 "\n@\n\\subsection{}",
                 "\n@\n\\subsubsection{}",
                 "\n@\n\\paragraph{}"
                        ); n.env <- length(news)
         news<-c(news,  # 140306
                 "\n\\item  ",
                 "\n\\includegraphics{}",
                 "\\emph{}",
                 "\\textbf{}",
                 "\\texttt{}"
               )[choice]
         # <bestimme Cursorzeile [[line]] von [[tworkwin]]> #
         pos <-as.character(tkindex(tworkwin,"insert"))
         line <-floor(as.numeric(pos))
         pos <- as.numeric(sub("^[0-9]*[.]","",pos))

  #       insert.line <- c(1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 0, 0, 0)[choice] + line
         insert.line <-  line
  #       insert.pos  <- c(0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 1, 1, 1)[choice]*pos 
         insert.pos  <- c(rep("end",13), rep(pos,3))[choice]
         cursor.line <- c(2,2,2,2,2,2,2, 2, 2, 2, 2, 1, 1, 0, 0, 0)[choice] + line
         cursor.pos <-  c(0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 1, 1, 1)[choice]*pos +
                        c(0,0,7,7,0,0,0, 9,12,15,11, 6,17, 6, 8, 8)[choice]
         # line <- if(choice <= n.env) paste(line+1,"0",sep=".") else paste(line,pos-1,sep=".")
         line <- paste(insert.line,insert.pos,sep=".")
         if(!exists("tworkwin"))
           tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

         if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
         tworkwin<-get("tworkwin",envir=revive.sys)
         worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
         if(nchar(worktext)<10000){
           worktext<-strsplit(worktext,"\n")[[1]]
         }else{
           base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
           worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
         }

         try(tkinsert(tworkwin,line,paste(news,collapse="\n")))
         melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))


       }
       if(!is.na(line)){
         line <- paste(cursor.line,cursor.pos,sep=".")
         tksee(tworkwin,line)
         tkmark.set(tworkwin, "insert", line)
         tkfocus(tworkwin)
       }
    })
    melde("InsertLaTeXEnv",2)
  }
  FindLaTeXSection<-function(){ # 091023
    melde("FindLaTeXSection",1)
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    no <- grep("\\\\s(.*)ection",worktext); if(0==length(no)) return()

    tchoices<-paste(worktext[no],worktext[no+1],sep=":: ")
    tchoices<-substring(tchoices,1,pmin(nchar(tchoices),60))
    tchoices<-sub("\\\\subsubsection","->->",tchoices)
    tchoices<-sub("\\\\subsection","->",tchoices)
    tchoices<-sub("\\\\section","",tchoices)
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    actual.chunk<-min(max(1,c(which(line<no),length(no)+1)[1]-1),length(no))

    newtop<-tktoplevel();tkwm.title(newtop,"section? -- Esc=Exit, Return=Selection")
    scr <- tkscrollbar(newtop, command=function(...)tkyview(tl,...))
    tl<-tklistbox(newtop,height=min(20,length(no)),width=60,selectmode="single",
                  yscrollcommand=function(...)tkset(scr,...),background="white")
    for(ch in tchoices) tkinsert(tl,"end",ch)
    tkselection.set(tl,actual.chunk-1)  # Default
    tksee(tl,max(actual.chunk-10,0))
    tkpack(tl,side="left",expand="yes",fill="y"); tkpack(scr,side="left",expand="yes",fill="y")
    tkbind(newtop,"<Escape>",function()tkdestroy(newtop))
    tkbind(newtop,"<Return>",function(){
       h<-as.numeric(tkcurselection(tl))+1
       tchoice<-tchoices[h]; 
       tchoice<-paste("section",tchoice,sep="")
       tchoice<-paste("\\",sub("^section->","subsection",
                               sub("^section->->","subsubsection",tchoice)),sep="")
       tchoice<-sub("::.*$","",tchoice)
       if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
       tworkwin<-get("tworkwin",envir=revive.sys)
       worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
       if(nchar(worktext)<10000){
         worktext<-strsplit(worktext,"\n")[[1]]
       }else{
         base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
         worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
       }

       line<-which(tchoice==substring(worktext,1,nchar(tchoice)))[1]
       # line<-no[as.numeric(tkcurselection(tl))+1]; tkdestroy(newtop)
       if(!is.na(line)){
         tksee(tworkwin,paste(line,".1",sep=""))
         tkmark.set(tworkwin, "insert", paste(line,".1",sep=""))
         tkfocus(tworkwin)
       }
    })
    melde("FindLaTeXSection",2)
  }
  ShowAboutRelax<-function(){
    melde("ShowAboutRelax",1)
    doc<-c("=================================================",
           "RELAX -- R Editor for Literate Analysis and lateX",
           "=================================================","",
           paste("version:","relax 1.3.15 - 140310"),"",
                 "relax() is designed to support processes of data analysis, report writing,",
                 "presentation, and programming. On start relax() creates a Tk/Tcl window",
                 "consisting a report and an output field. You can write down your ideas as",
                 "well as some chunks of code into the report field. The R statements will be",
                 "activated by pressing the button <EvalRCode> and the output will occur in",
                 "the output field. By clicking <Insert> the output is appended to the",
                 "report field and you can add your interpretations of the results.","",
                 "This style of work will lead to correct reports and after formatting by",
                 "LaTex you may get a very nice documentation of your analysis.","",
                 "The saved reports can be reloaded again in cause of checken, modification",
                 "and presentation of the results.","", 
                 "Have a look at the help page of relax and visit",
                 "http://www.wiwi.uni-bielefeld.de/com/wolf/software/relax.html","",
                 "Copyright (C) 2005 Hans Peter Wolf.","",
                 "This software is licensed under the GPL (version 2 or newer) terms.",
                 "It is free software and comes with ABSOLUTELY NO WARRANTY.","",
                # "The windows version may use some ingredients",
                # "of the noweb system",
                # " (Norman Ramsey -- http://www.eecs.harvard.edu/~nr/noweb/intro.html),",
                # "of the Img package of Jan Nijtmans (see:",
                # "http://wiki.tcl.tk/1404,",
                # "http://home.kpnplanet.nl/~J.Nijtmans@kpnplanet.nl/img.html",
                # "and http://members.chello.nl/~j.nijtmans/img.html",
                # "the package is found in http://sourceforge.net/projects/tkimg).",
                # "and gawk (http://www.gnu.org/software/gawk/gawk.html).",
                # "If used the license terms concerning the img package are found in the source", 
                # "file of the package, see: relax/src/tkimg1.3.tar.gz.","",
                 "---------------------------------------------------------","",""
                 )
    try(doc<-c(doc,scan(file=paste(path.package("relax"),"lib/gpl.txt",sep="/"),what="",sep="\n"))) # 130325 .path.package defunct
    .newl<-tktoplevel();tkwm.geometry(.newl,"+0+15")
    tkpack(tt<-tktext(.newl,height="33")) ## ((length(doc)))
    tkwm.title(.newl,paste("What's relax? ","Exit by Return or Escape!"))
    try(tkinsert(tt,"0.0",paste(doc,collapse="\n")))
    abbruch<-function(){tkdestroy(.newl); set.tclvalue("tvscandone",2)}
    tkbind(.newl,"<Escape>", function(){
                 tkdestroy(.newl)
    })
    tkbind(.newl,"<Return>", function(){
                 tkdestroy(.newl)
    })
    tkfocus(.newl)
    melde("ShowAboutRelax",2)
  }

  ShowShortCuts<-function(){
    melde("ShowShortCuts",1)
    keys<-c("Shortcuts of relax:","",
            if((.Platform$OS.type=="windows") | substring(version$os,1,5)=="linux")
                c("Alt-D:  move cursor one code chunk DOWN",
                  "Alt-U:  move cursor one code chunk UP",
                  "Alt-P:  plan new code chunk",
                  "Alt-E:  eval code chunk",
                  "Alt-W:  eval code chunk / ignore warnings",
                  "Alt-I:  insert output into report",
                  "Alt-S:  copy plot and generate postscript / JPEG file",
                  "Alt-R:  clear output field",
                  "Alt-T:  delete inserted output from report",
                  "Alt-H:  R help",
                  "Alt-F, Crtl-F:  Find text in report field"),
            if(substring(version$os,1,6)=="darwin" )
                c(#"Alt-D:  move cursor one code chunk DOWN",
                  #"Alt-U:  move cursor one code chunk UP",
                  "Apple-P:  plan new code chunk",
                  "Apple-E:  eval code chunk",
                  # "Apple-W:  eval code chunk / ignore warnings", # will destroy relax window
                  "Apple-I:  insert output into report",
                  "Apple-S:  copy plot and generate postscript / JPEG file",
                  "Apple-R:  clear output field",
                  "Apple-T:  delete output of text field",
                  # "Alt-H:  R help",
                  "Crtl-G:  Goto line number",
                  "<Fn><F3>:  Find text in report field again",
                  "Alt-F, Crtl-F:  Find text in report field") 
             else
                c("Crtl-G:  Goto line number",
                  "<F3>:  Find text in report field again"),  ## end of if
            "Crtl-S:  Save Report",
            "Crtl-P:  Process Report: save, weave, and latex Report",
                  "","Shortcuts of tcltk, may depend on OS / tcl/tk version:","",
                  "Crtl-Z:  UNDO",
                  "Crtl-C:  copy marked text",
                  "Crtl-V:  paste marked text",
                  "Crtl-Y:  paste marked text (=Crtl-V)",
                  "Crtl-X:  delete marked text",
                  "Crtl-W:  copy and delete marked text",
                  "Crtl-Shift-7:  mark all",
                 # "Crtl-F:  unmark text",
                 # "Crtl-A:  unmark text",
                 # "Crtl-P:  unmark text",
                  "Crtl-O:  break line",
                  "Crtl-E:  move cursor to the end of line",
                  "Crtl-A:  move cursor to the beginning of line",
                  "Crtl-B:  move cursor one character to the left",
                  "Crtl-F:  move cursor one character to the right",
                  "Crtl-N:  move cursor to the next line",
                  "Crtl-Pos1:  move cursor to the beginning of text",
                  "Crtl-End:  move cursor to the end of text",
                 # "Crtl-P:  move cursor one line back",
                  "Crtl-I:  insert tab",
                  "Crtl-T:  exchange characters",
                  "Crtl-D:  delete character right of cursor",
                  "Crtl-H:  delete character left of cursor",
                  "Crtl-K:  delete characters from cursor to end of line"     )
    .newl<-tktoplevel();tkwm.geometry(.newl,"+0+15")
    tkpack(tt<-tktext(.newl,height=length(keys)))
    tkwm.title(.newl,paste("shortcuts of relax and text field,","Exit by Return or Escape!"))
    try(tkinsert(tt,"0.0",paste(keys,collapse="\n")))
    abbruch<-function(){tkdestroy(.newl); set.tclvalue("tvscandone",2)}
    tkbind(.newl,"<Escape>", function(){
                 tkdestroy(.newl)
    })
    tkbind(.newl,"<Return>", function(){
                 tkdestroy(.newl)
    })
    tkfocus(.newl)
    melde("ShowShortCuts",2)
  }
  ShowInteractiveIcon<-""


  FindRFns<-function(){
    melde("FindRFns",1)
    frage<-"keyword for function search?"; set.tclvalue("tvinfo",string.sys)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>", function(){
        such<-tclvalue("tvinfo")
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        set.tclvalue("tvmess",paste(tclvalue("tvinfo"),": bitte etwas Geduld!"))
        found<-help.search(such)[[4]][,1:2]
        found<-paste(paste(found[,1],
                           substring("         ",1,pmax(1,10-nchar(found[,1]))),
                           found[,2]),collapse="\n")
        .newl<-tktoplevel();tkwm.geometry(.newl,"+0+15");tkpack(tt<-tktext(.newl))
        tkwm.title(.newl,paste("zu >",such,"< gefundene Funktionen, ",
                               "Beendigung durch Escape oder Return!"))
        try(tkinsert(tt,"0.0",paste(found,collapse="\n")))
        abbruch<-function(){tkdestroy(.newl);set.tclvalue("tvndone",2)}
        tkbind(.newl,"<Return>",abbruch);tkbind(.newl,"<Escape>",abbruch)
        tkfocus(.newl);tkwait.variable("tvndone")
        set.tclvalue("tvmess","relax")
      } # end of function
    )
    melde("FindRFns",2)
  }

  ReloadPlots<-function(){
    melde("ReloadPlots",1)
    if(!no.plots) { exclude.plots(tworkwin); createandshow.all.plots(tworkwin) }

    melde("ReloadPlots",2)
  }

  RemoveSavedPlot<-function(){
    melde("RemoveSavedPlot",1)
    # find cursor line
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    # check whether a plot is defined by cursor line
    if(1==length(grep("p..img src..p2", worktext[line])))   { line <- max(1,line - 2) } else 
      if(0==length(grep("includegraphics",worktext[line]))) { line <- max(1,line - 1) }
    if(0==length(grep("includegraphics",worktext[line]))){ cat("relax warning: no plot for deleting found"); return() }
    out.start <- line
    # extract name of plot
    picname <- sub("^.*includegraphics.*p[2]","p2",worktext[out.start])
    picname <- sub("\\}.*$","",picname)
    if(0==length(grep("^p[2].*[0-9][0-9]",picname))){ cat("relax warning: no plot for deleting found"); return() }
    # delete files 
    f.vec <- list.files(pattern=picname)
    if(0==length(f.vec)) { cat("relax warning: no plot for deleting found"); return() }
    for(pic in f.vec) { cat("file",pic,"deleted"); file.remove(pic) }
    # remove link lines from text
    out.end <- out.start 
    if(1==length(grep("p..img src..p2", worktext[out.start+1]))) out.end <- out.start+1
    if(1==length(grep("p..img src..p2", worktext[out.start+2]))) out.end <- out.start+2  
    worktext<-worktext[-(out.start:out.end)]
    # refresh report window
    if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
    tkdelete(tworkwin,"0.0","end")
    try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
    tksee(tworkwin,"end")
    melde("ak texthervor",1)
    tcl("markclear",tworkwin)
    tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
    tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
    tcl("marklinetypes",tworkwin)
    melde("ak texthervor",2)

    if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


    tkfocus(tworkwin)
    tkmark.set(tworkwin,"insert",paste(out.start-1,".0",sep=""))
    tksee(tworkwin,"insert")
    melde("RemoveSavedPlot",2)
  }

  RemoveALLPLOTS<-function(){
    melde("RemoveALLPLOTS",1)
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    pic.lines.set <- grep("^ *\\\\begin[{]center[}]\\\\includegraphics.*[{]p20.*\\\\end[{]center[}]",
                          worktext) #} 
    if(0==length(pic.lines.set)) return()

    res<-tkmessageBox(message=paste(if(language=="german") 
                                      "Sollen wirklich ALLE Bilder dauerhaft ENTFERNT werden?\n"
                                    else 
                                       "All Pictures will be removed!\nAre you shure?\n",
                                    paste(sub("^.*(p20...*[0-9]+).*$","\\1",worktext[pic.lines.set]),
                                          collapse="\n ")),
                      title="Exit",icon="warning",type="yesnocancel",default="no")
    if("externalptr"==mode(res))  res<-tclvalue(res)
    if(res=="cancel") return(); if(res!="yes") return()

    for(LINE in rev(pic.lines.set)){
      out.start <- LINE
      # extract name of plot
      picname <- sub("^.*includegraphics.*[{]p[2]","p2",worktext[out.start])  
      picname <- sub("\\}.*$","",picname) 
      if(0==length(grep("^p[2].*[0-9][0-9]",picname))){ next }
      # delete files 
      f.vec <- list.files(pattern=picname)
      if(0==length(f.vec)) { next }
      for(pic in f.vec) { cat("file",pic,"deleted"); file.remove(pic) }
      # remove link lines from text
      out.end <- out.start 
      if(1==length(grep("p..img src..p2", worktext[out.start+1]))) out.end <- out.start+1
      if(1==length(grep("p..img src..p2", worktext[out.start+2]))) out.end <- out.start+2  
      worktext<-worktext[-(out.start:out.end)]
    }
    # refresh report window
    if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
    tkdelete(tworkwin,"0.0","end")
    try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
    tksee(tworkwin,"end")
    melde("ak texthervor",1)
    tcl("markclear",tworkwin)
    tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
    tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
    tcl("marklinetypes",tworkwin)
    melde("ak texthervor",2)

    if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


    tkfocus(tworkwin)
    tkmark.set(tworkwin,"insert",paste(out.start-1,".0",sep=""))
    tksee(tworkwin,"insert")
    melde("RemoveALLPLOTS",2)
  }

  ReloadReportWidget<-function(){
    melde("ReloadReportWidget",1)
    revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys); fworkwin<-get("fworkwin",envir=revive.sys)
    del<-0.02; tkpack.forget(tworkwin); Sys.sleep(del)
    tkpack(tworkwin,expand="yes",fill="both"); Sys.sleep(del)
    .Tcl(paste("place ",fworkwin$ID," -relheight 0.35")); Sys.sleep(del)
    .Tcl(paste("place ",fworkwin$ID," -relheight 0.50"))
    melde("ReloadReportWidget",2)
  }

  RefreshChunkNumbers<-function(){
    melde("RefreshChunkNumbers",1)
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    code.start<-grep("^<<(.*)>>=",worktext)
    try(if(0<length(code.start)){ 
           worktext[code.start]<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext[code.start])
           worktext[code.start]<-paste(worktext[code.start]," (",1:length(code.start),")",sep="")
    })
    if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
    tkdelete(tworkwin,"0.0","end")
    try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
    tksee(tworkwin,"end")
    melde("ak texthervor",1)
    tcl("markclear",tworkwin)
    tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
    tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
    tcl("marklinetypes",tworkwin)
    melde("ak texthervor",2)

    if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


    tkmark.set(tworkwin, "insert", paste(line,"0",sep="."))
    tksee(tworkwin,paste(line,"0",sep="."))
    tkfocus(tworkwin)

    melde("RefreshChunkNumbers",1)
  } 

  ListUsedPlots<-function(){
    melde("ListUsedPlots",1)
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    files.used<-NULL
    include.lines<-grep("\\includegraphics",worktext,value=TRUE)
    if(length(include.lines)>0){
      include.lines<-sub("..includegraphics","",include.lines)
      include.lines<-sub("[}].*$","",include.lines)
      include.lines<-sub(".*[{]","",include.lines)
      inclinclude.lines<-sub("[}].*$","",include.lines)
      # cat("files of includegraphics macros of report"); print(include.lines)
      files.used<-c(files.used,include.lines)
    }
    include.lines<-grep("% <p><img src[=]",worktext,value=TRUE)
    if(length(include.lines)>0){
      include.lines<-sub(".*[=].","",include.lines)
      include.lines<-sub(".>.*$","",include.lines)
      # cat("files of html img tags of report:"); print(include.lines)
      include.lines<-sub("....$","",include.lines)
      files.used<-unique(c(files.used,include.lines))
    }
    if(length(files.used)>0){
      cat("graphics files of includegraphics macros or html img tags") 
      print(files.used)
    } else {
      cat("no uses of graphics files of includegraphics macros or html img tags found")
      return()
    }
    files.found<-list.files(pattern="^p20[0-1]+[0-9]+")
    if(length(files.found)==0) {
      cat("no graphics files of directory found")
      return()
    }
    files.found<-sub("\\....$","",files.found); files.found<-sub("\\...$","",files.found)
    h<-files.found[!files.found %in% files.used]
    if(length(h)>0){
      cat("graphics files of the directory not referenced here:")
      print(h)
    }
    h<-files.used[!files.used %in% files.found]
    if(length(h)>0){
      cat("graphics files referenced here but not found in the directory:")
      print(h)
    }
    melde("ListUsedPlots",2)
  }

  FormatTeXLines<-function(){
    melde("FormatTeXLines",1)
    # ermittle neue Zeilen:
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    text.start<-h<-grep("^@",c(worktext,"@"))
    texbegin<-rev(text.start[text.start<line])[1]; if(is.na(texbegin)) return()
    texend<-text.start[text.start>line][1]
    code.start<-grep("^<<(.*)>>=",worktext)
    c.start<-rev(code.start[code.start<line])[1]
    if(!is.na(c.start) &&  c.start>texbegin ) return()
    c.end<-code.start[code.start>line][1];  if(!is.na(c.end))texend<-min(texend,c.end)
    ttext<-worktext[(texbegin+1):(texend-1)]
    # formatiere Zeilen
    ttext<-c("\\documentclass{article}\\pagestyle{empty}\\parindent0mm",
             "\\begin{document}\\LARGE",ttext,"\\end{document}")
    cat(ttext,file="tmptmp.tex",sep="\n")
    # erstelle jpg-bild
    h<-gsub(" ","",gsub(":","",date())); h<-substring(h,4:nchar(h),4:nchar(h))
    jpgname<-paste(c("t",h[6:11],h[4:5],h[1:2],h[14:15],".jpg"),collapse="")
    ppmname<-paste(c("t",h[6:11],h[4:5],h[1:2],h[14:15],".ppm"),collapse="")
    if(substring(version$os,1,5)=="linux"){ #121114
      system("echo q | latex tmptmp.tex; dvips -E tmptmp.dvi")
      system(paste("convert tmptmp.ps ",jpgname))
      system(paste("convert tmptmp.ps tmptmp.ppm; pnmcrop tmptmp.ppm > ",ppmname))
      system("rm tmptmp.tex tmptmp.aux tmptmp.log tmptmp.ppm tmptmp.dvi tmptmp.ps")
    }
    if((.Platform$OS.type=="windows")){ #121114
      if(!(exists("ghostscript")&&2<=nchar(ghostscript)&&0<length(grep("[A-Za-z]",ghostscript)))) 
         return()
         if(!file.exists(paste(ghostscript,".exe",sep=""))){
            cat("ERROR: sorry, ghostscript not found!!!\n"); return()
         }
         ft.path<-file.path(relax.path,"lib")
         ft.path<-ft.path[file.exists(ft.path)][1]
         try(shell("echo q | latex tmptmp.tex"))
         try(shell("dvips -E tmptmp.dvi"))
         try(shell(paste(ghostscript,
    " -dBATCH -sDEVICE=ppmraw -quit -sOutputFile=tmptmp.ppm -dNOPAUSE tmptmp.ps", sep="")))
         try(shell(paste(ghostscript,
    " -dBATCH -sDEVICE=jpg -quit -sOutputFile=tmptmp.jpg -dNOPAUSE tmptmp.ps", sep="")))
         try(shell(file.path(ft.path,"\\pnmcrop   tmptmp.ppm   > tmptmp.crp")))
         try(shell(paste(ft.path,"\\ppmtojpeg tmptmp.crp   > ",jpgname,  sep="")))
         try(shell(paste("cp tmptmp.ppm   > ",ppmname,  sep="")))
         try(shell("del tmptmp.tex tmptmp.aux tmptmp.log tmptmp.ppm tmptmp.crp tmptmp.dvi tmptmp.ps"))
    }
    # passe workttext an, integriere Bild
    insertline<-texend+3
    news<-paste("\n%<!--latex-end-->\n\n% <p><img src=\"",jpgname,"\">\n",sep="")
    line<-paste(texend,".0",sep="")
    if(!exists("tworkwin"))
      tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    try(tkinsert(tworkwin,line,paste(news,collapse="\n")))
    melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))


    line<-paste(texbegin+1,".0",sep=""); news<-"%<!--latex-begin--\n"
    if(!exists("tworkwin"))
      tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    try(tkinsert(tworkwin,line,paste(news,collapse="\n")))
    melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))


    ##zeige Bild im Textfenster an##
    no.plots<-get("no.plots",envir=revive.sys) #081125
    if(!no.plots) createandshow.single.plot(tworkwin,insertline,ppmname,type="ppm") #121114
    ## old: show.single.plot(tworkwin,line,jpgname) 
    melde("FormatTeXLines",2)
  }

  DumpCodeChunk<-function(){
    melde("DumpCodeChunk",1)
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    code.start<-grep("^<<(.*)>>=",worktext)
    if(0==length(code.start)) return()
    code.start<-code.start[code.start<=line]
    if(0<length(code.start)){
      if((code.start<-max(code.start)+1)>length(worktext)) return()
      code.end  <-c(grep("^@",worktext),1+length(worktext))
      code.end  <-min(code.end[code.end>code.start])-1
      code<-worktext[code.start:code.end]
      code<-code[code!=""]
      if(length(weg.ab<-grep("^<<(.*)>>=",code))>0) code<-code[-(weg.ab:length(code))]
      if(length(code)==0 || code[1]=="@") code<-" "
      melde("code:",3,code,"\n")

    }
    if(0<length(code)){
      melde("vor tangle - Code:\n",3,code)
      if(length(grep("<<(.*)>>",code))>0 || length(grep(">#",code))>0){
        code.a    <- grep("^<<(.*)>>=",worktext)
        code.z    <- grep("^@",worktext)
        code.z    <- unlist(sapply(code.a ,function(x,y) y[y>x][1], code.z))
        if(any(h<-is.na(code.z))) code.z<-code.z[!h]
        ################
        # code.n    <- length(worktext)
        # change    <- rep(0,code.n); change[c(code.a ,code.z)]<-1
        # code.ch   <- worktext[1==(cumsum(change)%%2)]
        ## 080311
        ind<-rep(0,length(worktext)); copy<-0
        for(i in seq(ind)){
          if(i %in% code.z) copy<-0
          if(i %in% code.a) copy<-1
          ind[i]<-copy
        }
        code.ch   <- worktext[ind==1]
        ## cat("input-text, Anfaenge, Code.chunks"); print(worktext); print(code.a); print(code.ch)

        code.n    <- length(code.ch)

        code.ch<-gsub("@>>","DoSpCloseKl-esc",gsub("@<<","DoSpOpenKl-esc",code.ch))

        code.ch<-gsub("(.*)<<(.*)>>=(.*)","cOdEdEf\\2",code.ch)
        repeat{
          if(0==length(cand<-grep("<<(.*)>>",code.ch))) break
          code.ch<-unlist(strsplit(gsub("(.*)<<(.*)>>(.*)",
                     "\\1bReAkuSeChUnK\\2bReAk\\3",code.ch),"bReAk"))
        }
        code.ch<-code.ch[code.ch!=""]
        code.n<-length(code.ch)
        melde("code.ch:",3,code.ch,"\n")

        line.typ  <-rep("C",code.n)
        code.a    <-grep("cOdEdEf",code.ch)
        code.ch[code.a]<-substring(code.ch[code.a],8)
        line.typ[code.a]<-"D"
        code.use    <-grep("uSeChUnK",code.ch)
        code.ch[code.use]<-substring(code.ch[code.use],9)
        line.typ[code.use]<-"U"
        code.ext  <-grep("#<file",code.ch)
        line.typ[code.ext]<-"E"
        melde("code.ch:",3,code.ch,"\n")

        code.out<-"##act:##"

        def.names<-code.ch[code.a]
        use.names<- if(length(code.use)>0) code.ch[code.use] else NULL
        code.z<-c(if(length(code.a)>1) code.a[-1]-1, code.n)
        code.ch<-paste(line.typ,code.ch,sep="")
        melde("code.ch:",3,code.ch,"\n")

        melde("vor expand - Code:\n",3,code)
        melde("bearbeite aktuellen Chunk\n",3)
        ###<bestimme Cursorzeile [[line]] von [[tworkwin]]> not a good idea###
        ch.no<-length(grep("^<<(.*)>>=",worktext[1:line]))

        rows      <-c((code.a[ch.no]+1),code.z[ch.no])
        if(all(!is.na(rows))&&rows[1]<=rows[2]){
          rows<-rows[1]:rows[2]
          code.stack<-code.ch[rows]
          max.depth.refinements<-500; i<-1
          repeat{
             if((i<-i+1)>max.depth.refinements){ 
                 cat("ERROR: maximal number of expandations (",max.depth.refinements,
                     ") exceeded\n --- perhaps a unintended recursion ???")
                 return()
             }
             if(0==length(code.stack))break
             typ<-substring(code.stack[1],1,1)
             if("C"==typ||"E"==typ){
               n.lines<-sum(cumprod("C"==substring(code.stack,1,1)))
               code.out<-c(code.out, substring(code.stack[1:n.lines],2))
               code.stack<-code.stack[-(1:n.lines)]
             }
             if(length(code.stack)>0 && "U"==substring(code.stack[1],1,1)){
               if(any(found<-def.names==substring(code.stack[1],2))){
                 found<-seq(along=def.names)[found]; rows<-NULL
                 for(no in found){
                   if((code.a[no]+1)<=code.z[no]) rows<-c(rows,(code.a[no]+1):code.z[no])
                 }
                 code.stack<-c(code.ch[rows],code.stack[-1])
                 melde("found",0,found)
               }else{code.stack<-code.stack[-1]}
             }

          }
        }
        if(length(code.ext)>0){
          code.out<-code.out[code.out!=""]
          code.ext<-rev(grep(">#",code.out))
          found<-TRUE
          repeat{
            if(length(code.ext)==0) break

            if(!found){
              code.out[code.ext[1]]<-paste("# ??",code.out[code.ext[1]])
              cat("ERROR: External Chunk",code.out[code.ext[1]],"not found!!!\n")
              code.ext<-code.ext[-1]
            }

            found<-TRUE
            ext.name <- rev(unlist(strsplit(code.out[code.ext[1]],"#<file:")))[1]
            ext.name <- unlist(strsplit(unlist(strsplit(ext.name,">#"))[1],":"))
            ext.chunk<-ext.name[2]; ext.name <-ext.name[1]
            ext.name.n<-nchar(ext.name)
            if(ext.name.n >4 && ".rev"==substring(ext.name,ext.name.n-3,ext.name.n)){
              ext.name<-substring(ext.name,1,ext.name.n-4)
            }

            if(is.na(as.numeric(ext.chunk))){
              # tld untersuchen
              filename<-paste(ext.name[1],".rev",sep="")
              if(!file.exists(filename)){
                cat("ERROR: file",filename,"for expansion of code chunk not found!!!\n")
                ext.file<-"Error"
              }else{
                ext.file<-try(myscan(file=filename,what="",sep="\n"))
              }
              if("Error"==substring(unlist(ext.file)[1],1,5)){
                found <-FALSE; next
              }
              ext.file <-ext.file[grep("^<<(.*)>>=",ext.file)]
              if(!is.null(ext.file)){ found<-FALSE; next }
              ext.chunk<-grep(ext.chunk,ext.file)
            }

            filename<-paste(ext.name[1],".R",sep="")
            if(!file.exists(filename)){
              cat("Warning: file",filename,"not found!!!\n")
              cat("         file",filename,"is now generated!!!\n")
              try(tangleR(ext.name[1]))
            }
            if(!file.exists(filename)){
              ext.file<-"Error"
            }else{
              ext.file<-try(myscan(file=filename,what="",sep="\n"))
            }
            if("Error"==substring(unlist(ext.file)[1],1,5)){
              found <-FALSE; next
            }

            ext.chunk<-as.numeric(ext.chunk)
            a        <-grep(paste("#", ext.chunk,":",sep=""),ext.file)[1]
            z        <-grep(paste("#:",ext.chunk,    sep=""),ext.file)[1]
            if(is.na(a)){
              found <- FALSE; next
            }
            if(a<=z) ext.file <-ext.file[a:z]

            code.out <-c(code.out[1:(code.ext[1]-1)], ext.file,
                         if(length(code.out)>(code.ext[1]+1))
                           code.out[(code.ext[1]+1):length(code.out)]
                       )

            code.ext<-code.ext[-1]

          }

        }
        code.out<-c(code.out,"##:act##")

        code.out<-gsub("DoSpCloseKl-esc",">>",gsub("DoSpOpenKl-esc","<<",code.out))

        melde("Ende Rtangle-last\n",3)
        code<-code.out[code.out!=""]
        melde("nach expand\n",3,code)
      }

    } else { cat("ERROR: no code found!!!\n"); return() }
    filename<-tkgetSaveFile(filetypes="{{Paper Files} {.R}}",
                            title="name of file to save code chunk?", initialdir=getwd(),
                            defaultextension=".R", initialfile="RCode-act.R")
    if(!is.character(filename)) filename<-tclvalue(filename)
    if(filename=="") return("cancel")

    try.res <- try(cat(code,file=filename,sep="\n"))
  }

  CopyToEnd<-function(){
    melde("CopyToEnd",1)
    news<-tclvalue(tkget(toutwin,"0.0","end"))
    if(1<nchar(news)){
      news<-paste("\n@","output-start",news,"output-end\n",sep="\n")
      news<-gsub("\n+","\n",news)
      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      if(!exists("tworkwin"))
        tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

      pos.to.insert<-"end"
      if(0<length(grep("output-start",news))){
        tail<-rev(strsplit(tclvalue(tkget(tworkwin,"end - 3 lines","end")),"\n")[[1]])
        ltail<-length(tail)
        if( (0==length(grep("<<[*]>>=",tail[1:ltail]))) &&
           any(h<-("output-end"==substring(tail[1:ltail],1,11)))){
           news<-sub(".*output-start\n","",news)
           news<-sub("output-end","",news)
           h<-seq(along=h)[h][1]
           pos.to.insert<-paste("end -",h,"lines")
        }
      }
      try(tkinsert(tworkwin,pos.to.insert,paste(news,collapse="\n")))
      tksee(tworkwin,"end - 0 lines")
      melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      worktext<-""
      if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
      tkdelete(toutwin,"0.0","end")
      try(tkinsert(toutwin,"0.0",paste(worktext,collapse="\n")))
    }
    melde("CopyToEnd",2)
  }

  SetWorkPath<-function(){
    melde("SetWorkPath",1)
    path<-tkchooseDirectory(title="Reportverzeichniswahl",initialdir=getwd())
    if(!is.character(path)) path<-tclvalue(path)
    if(path=="") return() else try.res <- try(setwd(path))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){
      cat("new directory:",getwd(),"\n")
    }else{
      cat("ERROR: change of directory failed!!!\n")
    }
    melde("SetWorkPath",2)
  }

  DumpEnvironment<-function(){
    melde("DumpEnvironment",1)
    filename<-tkgetSaveFile(filetypes="{{Paper Files} {.relax-env}}",
                            title="name of ENVIRONMENT DUMP file?", initialdir=getwd(),
                            defaultextension=".relax-env", 
                            initialfile=sub("rev$","relax-env",workname.sys))
    if(!is.character(filename)) filename<-tclvalue(filename)
    if(filename=="") return("cancel")

    noenv<-TRUE # ignore environments 
    objlist<-ls(envir=revive.env) 
    if(noenv) objlist<-objlist[sapply(objlist,function(x)!is.environment(get(x,envir=revive.env)))]
    objlist<-objlist[substring(objlist,1,6)!="revive"]
    if(length(objlist)==0){
      cat("WARNING: no objects found for saving!!!\n")
      return()
    }
    cat("list of objects:\n"); print(objlist); Sys.sleep(.5)
      

    assign("list.of.env.obj",objlist,envir=revive.env)
    try.res<-try(eval(parse(text=
                    paste(sep="","dump(list.of.env.obj,file=\"",filename,"\")")),envir=revive.env))
    try(eval(parse(text="remove(\"list.of.env.obj\")"),envir=revive.env))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){
      cat(paste("Objects of environment saved in file:",filename,"!\n"))
    } else {
      cat("ERROR: saving of environment objects failed!!!\n")
    }
    melde("DumpEnvironment",2)
  }
  SaveEnvironment<-function(){
    melde("SaveEnvironment",1)
    filename<-tkgetSaveFile(filetypes="{{Paper Files} {.bin-relax-env}}",
                            title="name of BINARY SAVE file?", initialdir=getwd(),
                            defaultextension=".bin-relax-env", 
                            initialfile=sub("rev$","bin-relax-env",workname.sys))
    if(!is.character(filename)) filename<-tclvalue(filename)
    if(filename=="") return("cancel")

    noenv<-FALSE # do not ignore environments 
    objlist<-ls(envir=revive.env) 
    if(noenv) objlist<-objlist[sapply(objlist,function(x)!is.environment(get(x,envir=revive.env)))]
    objlist<-objlist[substring(objlist,1,6)!="revive"]
    if(length(objlist)==0){
      cat("WARNING: no objects found for saving!!!\n")
      return()
    }
    cat("list of objects:\n"); print(objlist); Sys.sleep(.5)
      

    assign("list.of.env.obj",objlist,envir=revive.env)
    try.res<-try(eval(parse(text=
                    paste(sep="","save(list=list.of.env.obj,file=\"",filename,"\")")),envir=revive.env))
    try(eval(parse(text="remove(\"list.of.env.obj\")"),envir=revive.env))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){
      cat(paste("Objects of environment saved in file:",filename,"!\n"))
    } else {
      cat("ERROR: saving of environment objects failed!!!\n")
    }
    melde("SaveEnvironment",2)
  }
  LoadEnvironment<-function(){
    melde("LoadEnvironment",1)
    initialfile<-sub("rev$","bin-relax-env",workname.sys)
    if(!file.exists(initialfile)) initialfile<-sub("rev$","relax-env",workname.sys)
    filename<-tkgetOpenFile(filetypes="{{Paper Files} {.relax-env .bin-relax-env}}",
                            title="Select *relax-env file to be loaded!",
                            defaultextension=".relax-env",
                            initialfile=initialfile,
                            initialdir=getwd())
    if(!is.character(filename)) filename<-tclvalue(filename)
    if(filename=="") return("cancel")

    if(0<length(grep("bin-relax-env$",filename))){
      try.res<-try(load(filename,envir=revive.env))
    } else {
      obj<-scan(filename,""); ind<-grep("<-",obj)-1   
      obj[ind]<-paste(";",obj[ind]); obj<-sub("^;","",paste(obj,collapse=""))
      try.res<-try(eval(parse(text=obj),envir=revive.env))
    }
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){
      cat(paste("Objects of file:",filename,"\nloaded into environment!\n"))
      noenv<-FALSE # do not ignore environments 
      objlist<-ls(envir=revive.env) 
      if(noenv) objlist<-objlist[sapply(objlist,function(x)!is.environment(get(x,envir=revive.env)))]
      objlist<-objlist[substring(objlist,1,6)!="revive"]
      if(length(objlist)==0){
        cat("WARNING: no objects found for saving!!!\n")
        return()
      }
      cat("list of objects:\n"); print(objlist); Sys.sleep(.5)
        

    } else {
      cat("ERROR: saving of environment objects failed!!!\n")
    }
    melde("LoadEnvironment",2)
  }
  CleanEnvironment<-function(){
    melde("CleanEnvironment",1)
    noenv<-FALSE # do not ignore environments 
    objlist<-ls(envir=revive.env) 
    if(noenv) objlist<-objlist[sapply(objlist,function(x)!is.environment(get(x,envir=revive.env)))]
    objlist<-objlist[substring(objlist,1,6)!="revive"]
    if(length(objlist)==0){
      cat("WARNING: no objects found for saving!!!\n")
      return()
    }
    cat("list of objects:\n"); print(objlist); Sys.sleep(.5)
      

    res<-tkmessageBox(message=
        if(language=="german") 
          "Sollen die Objekte der Umgebung wirklich entfernt werden"
        else 
          "Do you want to delete all the objects of the environment?",
        title="Clean Environment?",icon="warning",type="yesnocancel",default="no")
    res<-tclvalue(res)
    if(res=="yes"){    
      objlist<-c(objlist,"objlist"); assign("objlist",objlist,envir=revive.env)
      try.res<-try(eval(parse(text="remove(list=objlist)"),envir=revive.env))
      cat("Objects of environment have been deleted!\n")
      # print(try.res)
    }
    melde("CleanEnvironment",2)
  }

  SaveReport<-function(){
    melde("SaveReport",1)
    filename<-tkgetSaveFile(filetypes="{{Paper Files} {.rev}}",
                            title="name of report file?", initialdir=getwd(),
                            defaultextension=".rev", initialfile=workname.sys)
    if(!is.character(filename)) filename<-tclvalue(filename)
    if(filename=="") return("cancel")

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
    worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
    worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
    if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
       worktext<-c(worktext,"@\n\\end{document}")
    worktext<-TcltoWin.write(worktext)
    try.res <- try(cat(worktext,file=filename,sep="\n"))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){ cat(paste("file",filename,"saved\n")) }

    if(ok){
      workname.sys<-sub(paste(".*",.Platform$file.sep,sep=""),"",filename)
      lworkname.sys<-get("lworkname.sys",envir=revive.sys)
      tkconfigure(lworkname.sys,text=paste(workname.sys,""))
      assign("workname.sys",workname.sys,envir=revive.sys)
      # tkwm.title(TopW,paste("RELAX:",workname.sys))
      tkwm.title(TopW,paste(if(but.Wizardry=="simple") "redit:" else "RELAX:",workname.sys))                      

      olddir<-getwd(); newdir<-sub("(.*)/(.*)","\\1",filename)
      if(olddir!=newdir){
        rev<-grep("includegraphics",worktext,value=TRUE)
        if(0<length(rev)){
          # pics<-sub("(.*)\\{(p.*\\.ps)\\}(.*)","\\2",rev) #081121
          pics<-sub("(.*)\\{(p.*)\\}(.*)","\\2.ps",rev); pics<-sub("\\.ps\\.ps$",".ps",pics)
          pics<-c(pics,sub("ps$","jpg",pics),sub("ps$","ppm",pics)) # 121114
          mess<-NULL
          for(pic in pics) {
            if(file.exists(pic)){
              a<-file.copy(pic,paste(newdir,pic,sep="/"))
              mess<-c(mess,if(a) paste("picture",pic,"copied to",newdir) else
                                 paste("warning: picture", pic,"not copied"))
            }
          }
          cat(mess,sep="\n")
          "process of copying file(s) finished"
        }
      }


      h <- strsplit(filename,.Platform$file.sep)[[1]]
      workname.sys<-rev(h)[1]
      workpath.sys<-paste(h[-length(h)],collapse=.Platform$file.sep)
      setwd(workpath.sys)

      melde(paste("report file",filename,"saved\n"),0)
      melde(paste("w",workname.sys),"cmd.msg")
    } else {
      cat("ERROR: write operation failed!!!\n")
    }
    melde("SaveReport",2)
    return("ok")
  }
  SaveHtml<-function(){
    melde("SaveHtml",1)
    filename<-tkgetSaveFile(filetypes="{{Paper Files} {.rev}}",
                            title="name of report file?", initialdir=getwd(),
                            defaultextension=".rev", initialfile=workname.sys)
    if(!is.character(filename)) filename<-tclvalue(filename)
    if(filename=="") return("cancel")

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
    worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
    worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
    if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
       worktext<-c(worktext,"@\n\\end{document}")
    worktext<-TcltoWin.write(worktext)
    try.res <- try(cat(worktext,file=filename,sep="\n"))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){ cat(paste("file",filename,"saved\n")) }

    if(ok){
      workname.sys<-sub(paste(".*",.Platform$file.sep,sep=""),"",filename)
      lworkname.sys<-get("lworkname.sys",envir=revive.sys)
      tkconfigure(lworkname.sys,text=paste(workname.sys,""))
      assign("workname.sys",workname.sys,envir=revive.sys)
      # tkwm.title(TopW,paste("RELAX:",workname.sys))
      tkwm.title(TopW,paste(if(but.Wizardry=="simple") "redit:" else "RELAX:",workname.sys))                      

      olddir<-getwd(); newdir<-sub("(.*)/(.*)","\\1",filename)
      if(olddir!=newdir){
        rev<-grep("includegraphics",worktext,value=TRUE)
        if(0<length(rev)){
          # pics<-sub("(.*)\\{(p.*\\.ps)\\}(.*)","\\2",rev) #081121
          pics<-sub("(.*)\\{(p.*)\\}(.*)","\\2.ps",rev); pics<-sub("\\.ps\\.ps$",".ps",pics)
          pics<-c(pics,sub("ps$","jpg",pics),sub("ps$","ppm",pics)) # 121114
          mess<-NULL
          for(pic in pics) {
            if(file.exists(pic)){
              a<-file.copy(pic,paste(newdir,pic,sep="/"))
              mess<-c(mess,if(a) paste("picture",pic,"copied to",newdir) else
                                 paste("warning: picture", pic,"not copied"))
            }
          }
          cat(mess,sep="\n")
          "process of copying file(s) finished"
        }
      }


      h <- strsplit(filename,.Platform$file.sep)[[1]]
      workname.sys<-rev(h)[1]
      workpath.sys<-paste(h[-length(h)],collapse=.Platform$file.sep)
      setwd(workpath.sys)

      melde(paste("report file",filename,"saved\n"),0)
      melde(paste("w",workname.sys),"cmd.msg")
    } else {
      cat("ERROR: write operation failed!!!\n")
    }
    {
      melde("weaveRhtml",1)
      weaveRhtml(workname.sys,
                         replace.umlaute=replace.umlaute.sys)
      melde("weaveRhtml",2)
    }
    melde("SaveHtml",2)
    return("ok")
  }
  SaveAsConsoleStyleFile<-function(){
    melde("SaveAsConsoleStyleFile",1)
    filename<-tkgetSaveFile(filetypes="{{TEXT-FILE!} {*.*}}",
                            title="name of file?", initialdir=getwd(),
                            defaultextension="", initialfile=sub(".rev$",".rev",workname.sys))
    if(!is.character(filename)) filename<-tclvalue(filename)
    if(filename=="") return("cancel")

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    if(any(h<-(substring(worktext,1,1)==">"))){ # 091026
      worktext[h]<-sub("^.","@SpLiTtE?<<*>>=SpLiTtE?",worktext[h])
      worktext<-substring(unlist(strsplit(paste("",worktext),"SpLiTtE")),2)
      if(is.vector(line)) line<-line+2*sum(which(h)<=line)
    }

    worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
    worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
    worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
    if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
       worktext<-c(worktext,"@\n\\end{document}")
    worktext<-TcltoWin.write(worktext)
    try.res <- try(cat(worktext,file=filename,sep="\n"))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){ cat(paste("file",filename,"saved\n")) }

    filename<-sub("\\.rev$",".txt",filename)      
    worktext<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext)
    #worktext<-sub("^output-start","\\\\begin{verbatim}",worktext)
    #worktext<-sub("^output-end","\\\\end{verbatim}",worktext)
    worktext<-sub("^output-start","",worktext)
    worktext<-sub("^output-end","",worktext)
    #if(0==length(grep("\\\\end\\{document\\}",worktext))) #2.1.0
    #  worktext<-c(worktext,"@\n\\end{document}")
    ## remove noweb items
    code.ch<-c("@",worktext,"@")
    # translation of code chunks in noweb style
    code.chunk.start<-paste("^<","<.*>",">=",sep="")
    code.a<- grep(code.chunk.start,code.ch)
    if(0<length(code.a)){
      code.z<-grep("^@",code.ch)
      code.z<-unlist(sapply(code.a ,function(x,y)min(y[y>x]),code.z))
      code.n<-length(code.ch)
      change<-rep(0,code.n); change[c(code.a ,code.z)]<-1
      ind<-(1==(cumsum(change)%%2))[-length(change)]
      def.lines<-diff(c(FALSE,ind))>0
      code.ch<-code.ch[!def.lines]; ind<-ind[!def.lines]
      code.ch<-sub(paste("^ *<","<",sep=""),"# <<",code.ch) # Auskommentierung von uses
      code.ch[ind]<-paste("+",code.ch[ind])
      code.ch[ind]<-sub("^[+][ \t]*$"," ",code.ch[ind])
      ind<-diff(c(FALSE,ind))>0
      code.ch[ind]<-paste(">",substring(code.ch[ind],3))
    #  code.ch[ind]<-paste("+",code.ch[ind])
    #  ind<-diff(c(FALSE,ind))>0
    #  code.ch[ind]<-sub("^[+] ","> ",code.ch[ind])
    } else { "no code" }
    worktext<-code.ch
    # translation of code chunks in console style
    code.a<-grep("^>",code.ch)
    if(0<length(code.a)){
      code.z<-grep("^@",code.ch)
      code.z<-unlist(sapply(code.a ,function(x,y)min(y[y>x]),code.z))
      code.n<-length(code.ch)
      change<-rep(0,code.n); change[c(code.a ,code.z)]<-1
      ind<-(1==(cumsum(change)%%2))[-length(change)]
      def.lines<-diff(c(FALSE,ind))>0
      ###  code.ch<-code.ch[!def.lines]; ind<-ind[!def.lines]
      code.ch<-sub(paste("^ *<","<",sep=""),"# <<",code.ch) # Auskommentierung von uses
      code.ch[ind]<-ifelse(!(substring(code.ch[ind],1,1)%in%c(">","+")),paste("+",code.ch[ind]),code.ch[ind])
      code.ch[ind]<-sub("^[+][ \t]*$"," ",code.ch[ind])
      ind<-diff(c(FALSE,ind))>0
      ### code.ch[ind]<-paste(">",substring(code.ch[ind],3))
    } else { "no code" }
    worktext<-code.ch
    worktext<-worktext[-grep("^@",worktext)]
    worktext<-TcltoWin.write(worktext)
    try.res <- try(cat(worktext,file=filename,sep="\n")) ### .txt written
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){ cat(paste("file",filename,"saved\n")) }

    # schreibe  script.R-file
    a<-worktext; ind<-grep("^[+>]",a)
    if(0!=length(ind)){
      a<-a[ind]; ind<-grep("^[>]",a); no<-seq(ind)
      a[ind]<-paste("#",no,": >\n ",substring(a[ind],2),sep="")
      ind<-grep("^[+]",a)
      if(0!=length(ind)) a[ind]<-paste("", substring(a[ind],2))
      filename<-sub(".txt$",".R",filename)
      try.res <- try(cat(a,file=filename,sep="\n")) ### .R
      if(is.function(try.res)){
        ok <- "OK"
      } else {
        if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
        if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
        if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
        if(!is.character(ok)) { ok <- "OK" }
      }
      if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
        ok<-FALSE
        error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
        cat(error.msg,"\n")
        if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
           cat("A warning message stopped the evaluation!",
                 "If you want to\nevaluate the code anyway",
                 "evaluate code by:\n>WarnEval<")
        # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
      } else { ok<-TRUE }


      if(ok){ cat(paste("file",filename,"saved\n")) }
    }

    # schreibe  tex-file
    # verbatim-Separationszeichen festlegen: SEP
    filename<-sub("\\.R$",".tex",filename); filename<-sub("\\.txt$",".tex",filename) 

    .Tcl("set XYZ [encoding system]")
    UTF<-is.UTF<- 0<length(grep("utf",tclvalue("XYZ")))


    SEP<- if(!is.UTF){ eval(parse(text='"\\267"')) }else{ trennzeichen<-"\140" }
    ind<-grep("^[>+]",worktext)
    worktext[ind]<-paste("\\rule{0mm}{0mm}\\\\\\verb",SEP,worktext[ind],SEP,sep="")
    worktext[!ind]<-gsub("\\[\\[(.*)\\$(.*)\\]\\]",paste("[[\\1\\\\","$\\2]]",sep=""),worktext[!ind])
    try.res <- try(cat(worktext,file=filename,sep="\n")) ### .tex
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){ cat(paste("file",filename,"saved\n")) }

    if(ok){
      filename<-sub("\\.tex$","",filename)
      filename<-sub("\\.txt$","",filename)
      filename<-sub("\\.rev$","",filename)
      filename<-sub("\\.R$","",filename)
      filename<-paste(filename,".rev",sep="")
      workname.sys<-sub(paste(".*",.Platform$file.sep,sep=""),"",filename)
      lworkname.sys<-get("lworkname.sys",envir=revive.sys)
      tkconfigure(lworkname.sys,text=paste(workname.sys,""))
      assign("workname.sys",workname.sys,envir=revive.sys)
      # tkwm.title(TopW,paste("RELAX:",workname.sys))
      tkwm.title(TopW,paste(if(but.Wizardry=="simple") "redit:" else "RELAX:",workname.sys))                      

      olddir<-getwd(); newdir<-sub("(.*)/(.*)","\\1",filename)
      if(olddir!=newdir){
        rev<-grep("includegraphics",worktext,value=TRUE)
        if(0<length(rev)){
          # pics<-sub("(.*)\\{(p.*\\.ps)\\}(.*)","\\2",rev) #081121
          pics<-sub("(.*)\\{(p.*)\\}(.*)","\\2.ps",rev); pics<-sub("\\.ps\\.ps$",".ps",pics)
          pics<-c(pics,sub("ps$","jpg",pics),sub("ps$","ppm",pics)) # 121114
          mess<-NULL
          for(pic in pics) {
            if(file.exists(pic)){
              a<-file.copy(pic,paste(newdir,pic,sep="/"))
              mess<-c(mess,if(a) paste("picture",pic,"copied to",newdir) else
                                 paste("warning: picture", pic,"not copied"))
            }
          }
          cat(mess,sep="\n")
          "process of copying file(s) finished"
        }
      }


      h <- strsplit(filename,.Platform$file.sep)[[1]]
      workname.sys<-rev(h)[1]
      workpath.sys<-paste(h[-length(h)],collapse=.Platform$file.sep)
      setwd(workpath.sys)

    } else {
      cat("ERROR: write operation failed!!!\n")
    }
    melde("SaveAsConsoleStyleFile",2)
    return("ok")
  }

  # ermittle neue Zeilen:
  SaveDiffReport<-function(){
    melde("SaveDiffReport",1)
    filename<-tkgetSaveFile(filetypes="{{Paper Files} {.html}}",
                            title="Diff-File-Name?", initialdir=getwd(),
                            defaultextension=".html", initialfile=
                paste("diff",workname.sys,sep="-"))
    if(!is.character(filename)) filename<-tclvalue(filename)
    if(filename=="") return("cancel")

    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    try.res <-try(cat(worktext,file=get("tmp.file.name",envir=revive.sys)
,sep="\n"))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(!ok) return()
    # stelle Differenzen fest
    system(paste("diff ",workname.sys," ",get("tmp.file.name",envir=revive.sys)
," | grep \"^[0-9]\" > ",get("tmp.file.name",envir=revive.sys)
))
    difflines<-scan(get("tmp.file.name",envir=revive.sys)
,"",blank.lines.skip=FALSE,sep="\n")
    difflines<-sub("^([0-9]*)[a-z]","",difflines)
    difflines<-sub(",",":",difflines)
    difflines<-paste("c(",paste(difflines,collapse=","),")")
    try(difflines<-eval(parse(text=difflines)))
    if(length(difflines)>0 || "ERROR"!=substring(unlist(difflines)[1],1,5)){
      # suche Chunks:
      worktext<-c("@",worktext,"@")
      difflines<-difflines+1
      chunkbegins<-grep("(^@)|(^<<(.*)>=)", worktext)
      chunk.log<-rep(FALSE,length(chunkbegins))
      ch.new<- 0<hist(difflines,plot=FALSE,
                      breaks=c(chunkbegins,length(worktext)+1)-0.5)$counts
      ch.new<- ch.new | c(ch.new[-1],FALSE) | c(FALSE,ch.new[-length(ch.new)])
      ch.new<-seq(ch.new)[ch.new]
      # extrahiere Chunks:
      lines<-paste(chunkbegins[ch.new],":",
               c(chunkbegins[-1],chunkbegins[length(chunkbegins)]+1)[ch.new]-1)
      lines<-paste(lines,collapse=",")
      lines<-try(eval(parse(text=paste("c(",lines,")"))))
      if(ch.new[1]) lines<-lines[-1]
      if(ch.new[length(ch.new)]) lines<-lines[-length(lines)]
      difftext<-worktext[unique(lines)]
      # speichere difftext
      ### SaveAsHtml(difftext,filename)
      try.res <-try(cat(difftext,file=filename,sep="\n"))
    }
    melde("SaveDiffReport",2)
  }

  OpenReport<-function(){
    melde("OpenReport",1)
    filename<-tkgetOpenFile(filetypes="{{Paper Files} {.rev}}",
                            title="Select .rev file to be loaded!",
                            defaultextension=".rev", initialfile=workname.sys,
                            initialdir=getwd())
    if(!is.character(filename)) filename<-tclvalue(filename)
    if(filename=="") return("cancel")

    try.res<-try(myscan(filename,"",sep="\n",blank.lines.skip=FALSE))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){
      h <- strsplit(filename,.Platform$file.sep)[[1]]
      workname.sys<-rev(h)[1]
      workpath.sys<-paste(h[-length(h)],collapse=.Platform$file.sep)
      setwd(workpath.sys)

      workname.sys<-sub(paste(".*",.Platform$file.sep,sep=""),"",filename)
      lworkname.sys<-get("lworkname.sys",envir=revive.sys)
      tkconfigure(lworkname.sys,text=paste(workname.sys,""))
      assign("workname.sys",workname.sys,envir=revive.sys)
      # tkwm.title(TopW,paste("RELAX:",workname.sys))
      tkwm.title(TopW,paste(if(but.Wizardry=="simple") "redit:" else "RELAX:",workname.sys))                      

      try.res<-WinToTcl.read(try.res)
          ## Eintrag mit Entfernung des bisherigen Inhalts:
          ## worktext<-paste(try.res, collapse="\n")
          ## <<schreibe [[worktext]] ins Arbeitsfenster>>
      news<-c("",try.res)
      if(!exists("tworkwin"))
        tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

      pos.to.insert<-"end"
      if(0<length(grep("output-start",news))){
        tail<-rev(strsplit(tclvalue(tkget(tworkwin,"end - 3 lines","end")),"\n")[[1]])
        ltail<-length(tail)
        if( (0==length(grep("<<[*]>>=",tail[1:ltail]))) &&
           any(h<-("output-end"==substring(tail[1:ltail],1,11)))){
           news<-sub(".*output-start\n","",news)
           news<-sub("output-end","",news)
           h<-seq(along=h)[h][1]
           pos.to.insert<-paste("end -",h,"lines")
        }
      }
      try(tkinsert(tworkwin,pos.to.insert,paste(news,collapse="\n")))
      tksee(tworkwin,"end - 0 lines")
      melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      worktext<-sub("^\\\\begin\\{verbatim\\}","output-start",worktext)
      worktext<-sub("^\\\\end\\{verbatim\\}","output-end",worktext)

      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      code.start<-grep("^<<(.*)>>=",worktext)
      try(if(0<length(code.start)){ 
             worktext[code.start]<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext[code.start])
             worktext[code.start]<-paste(worktext[code.start]," (",1:length(code.start),")",sep="")
      })
      if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
      tkdelete(tworkwin,"0.0","end")
      try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
      tksee(tworkwin,"end")
      melde("ak texthervor",1)
      tcl("markclear",tworkwin)
      tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
      tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
      tcl("marklinetypes",tworkwin)
      melde("ak texthervor",2)

      if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


      tkmark.set(tworkwin, "insert", paste(line,"0",sep="."))
      tksee(tworkwin,paste(line,"0",sep="."))
      tkfocus(tworkwin)

      RunStart()
      melde("ak texthervor",1)
      tcl("markclear",tworkwin)
      tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
      tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
      tcl("marklinetypes",tworkwin)
      melde("ak texthervor",2)

      if(!no.plots) { exclude.plots(tworkwin); createandshow.all.plots(tworkwin) }

      melde(paste("r", workname.sys), "cmd.msg")
    } else { cat("ERROR: File",filename,"not found!!!\n") }
    melde("OpenReport",2)
  }

  OpenTextFile<-function(){
    melde("OpenTextFile",1)
    filename<-tkgetOpenFile(filetypes="{{TEXT-FILE!} {*.*}}",
                            title="Select file to be loaded!",
                            initialfile=sub(".rev$",".txt",workname.sys),
                            initialdir=getwd())
    if(!is.character(filename)) filename<-tclvalue(filename)
    if(filename=="") return("cancel")

    try.res<-try(myscan(filename,"",sep="\n",blank.lines.skip=FALSE))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){
      filename<-paste(sub(".txt$","",filename),".rev",sep="")
      h <- strsplit(filename,.Platform$file.sep)[[1]]
      workname.sys<-rev(h)[1]
      workpath.sys<-paste(h[-length(h)],collapse=.Platform$file.sep)
      setwd(workpath.sys)

      workname.sys<-sub(paste(".*",.Platform$file.sep,sep=""),"",filename)
      lworkname.sys<-get("lworkname.sys",envir=revive.sys)
      tkconfigure(lworkname.sys,text=paste(workname.sys,""))
      assign("workname.sys",workname.sys,envir=revive.sys)
      # tkwm.title(TopW,paste("RELAX:",workname.sys))
      tkwm.title(TopW,paste(if(but.Wizardry=="simple") "redit:" else "RELAX:",workname.sys))                      

      try.res<-WinToTcl.read(try.res)
          try.res<-sub("^[ \t]*$","\n@",try.res)
          try.res<-sub("^\\+ ","> ",try.res)
          code.line<-"> "==substring(try.res,1,2)
          try.res[code.line]<-substring(try.res[code.line],3)
          first.code.line<-diff(c(FALSE,code.line))>0
          code.chunk.start<-paste("<","<*>",">=",sep="")
          try.res[first.code.line]<-paste(code.chunk.start,try.res[first.code.line],sep="\n")
          first.text.line<-rev(diff(c(FALSE,rev(code.line)))>0)
          try.res[first.text.line]<-paste(try.res[first.text.line],"@",sep="\n")
          
      ## Eintrag mit Entfernung des bisherigen Inhalts:
      ## worktext<-paste(try.res, collapse="\n"); <<schreibe [[worktext]] ins Arbeitsfenster>>
      news<-c("",try.res)
      if(!exists("tworkwin"))
        tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

      pos.to.insert<-"end"
      if(0<length(grep("output-start",news))){
        tail<-rev(strsplit(tclvalue(tkget(tworkwin,"end - 3 lines","end")),"\n")[[1]])
        ltail<-length(tail)
        if( (0==length(grep("<<[*]>>=",tail[1:ltail]))) &&
           any(h<-("output-end"==substring(tail[1:ltail],1,11)))){
           news<-sub(".*output-start\n","",news)
           news<-sub("output-end","",news)
           h<-seq(along=h)[h][1]
           pos.to.insert<-paste("end -",h,"lines")
        }
      }
      try(tkinsert(tworkwin,pos.to.insert,paste(news,collapse="\n")))
      tksee(tworkwin,"end - 0 lines")
      melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      worktext<-sub("^\\\\begin\\{verbatim\\}","output-start",worktext)
      worktext<-sub("^\\\\end\\{verbatim\\}","output-end",worktext)

      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      code.start<-grep("^<<(.*)>>=",worktext)
      try(if(0<length(code.start)){ 
             worktext[code.start]<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext[code.start])
             worktext[code.start]<-paste(worktext[code.start]," (",1:length(code.start),")",sep="")
      })
      if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
      tkdelete(tworkwin,"0.0","end")
      try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
      tksee(tworkwin,"end")
      melde("ak texthervor",1)
      tcl("markclear",tworkwin)
      tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
      tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
      tcl("marklinetypes",tworkwin)
      melde("ak texthervor",2)

      if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


      tkmark.set(tworkwin, "insert", paste(line,"0",sep="."))
      tksee(tworkwin,paste(line,"0",sep="."))
      tkfocus(tworkwin)

      melde("ak texthervor",1)
      tcl("markclear",tworkwin)
      tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
      tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
      tcl("marklinetypes",tworkwin)
      melde("ak texthervor",2)

      if(!no.plots) { exclude.plots(tworkwin); createandshow.all.plots(tworkwin) }

      melde(paste("r", workname.sys), "cmd.msg")
    } else { cat("ERROR: File",filename,"not found!!!\n") }
    melde("OpenTextFile",2)
  }
  OpenRevbook<-function(){
    melde("OpenRevbook",1)
    h<-file.path(relax.path,"rev")
    filename<-tkgetOpenFile(filetypes="{{Paper Files} {.rev}}",
                            title="Choose paper to be loaded",
                            defaultextension=".rev", initialfile=workname.sys,
                            initialdir=h)
    if(!is.character(filename)) filename<-tclvalue(filename)
    if(filename=="") return("cancel")

    try.res<-try(myscan(filename,"",sep="\n",blank.lines.skip=FALSE))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){
      workname.sys<-sub(paste(".*",.Platform$file.sep,sep=""),"",filename)
      lworkname.sys<-get("lworkname.sys",envir=revive.sys)
      tkconfigure(lworkname.sys,text=paste(workname.sys,""))
      assign("workname.sys",workname.sys,envir=revive.sys)
      # tkwm.title(TopW,paste("RELAX:",workname.sys))
      tkwm.title(TopW,paste(if(but.Wizardry=="simple") "redit:" else "RELAX:",workname.sys))                      

      code.ch<-grep("^<<(.*)>>=",try.res)
      if(0<length(code.ch))
         try.res[code.ch]<-paste(try.res[code.ch], " (",1:length(code.ch),")",sep="")
         try.res<-WinToTcl.read(try.res)
          worktext<-paste(try.res, collapse="\n") # alter Inhalt wird entfernt
          if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
          tkdelete(tworkwin,"0.0","end")
          try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
          tksee(tworkwin,"end")
          melde("ak texthervor",1)
          tcl("markclear",tworkwin)
          tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
          tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
          tcl("marklinetypes",tworkwin)
          melde("ak texthervor",2)

          if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


        RunStart()
      melde("ak texthervor",1)
      tcl("markclear",tworkwin)
      tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
      tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
      tcl("marklinetypes",tworkwin)
      melde("ak texthervor",2)

      if(!no.plots) { exclude.plots(tworkwin); createandshow.all.plots(tworkwin) }

      melde(paste("r", workname.sys), "cmd.msg")
    } else { cat("ERROR: File",filename,"not found!!!\n") }
    melde("OpenRevbook",2)
  }
  LoadRwtools<-function(){
    melde("LoadRwtools",1)
    path<-file.path(relax.path,"rev/robj.R")
    fns<-readLines(path) #2.1.0
    try(eval(parse(text=fns),
                 envir=pos.to.env(which(path.package("relax")==searchpaths())))) # 130325 .path.package defunct
    melde("LoadRwtools",2)
  }

  ViewReport.html<-function(){
    melde("ViewReport.html",1)
    n<-nchar(filename<-workname.sys)
    if(is.null(n)||5>n||substring(filename,n-3,n)!=".rev"){
      cat("ERROR: file name not ok!!!\n"); return()
    }
    if( (.Platform$OS.type=="windows") ){
           fname<-paste(sub("rev$","html",filename),sep="")
           browser.exe<- if(browser.sys=="") "start " else browser.sys
           res<-shell(paste(browser.exe, fname),wait=FALSE)
           if(res!=0){
             cat("ERROR: browser hasn't been started successfully \n")
             cat("please, start and load file on foot\n")
           }
      }else{
          fname<-paste(sub("rev$","html",filename),sep="")
          if(browser.sys!=""){
                try(system(paste(browser.sys,fname),wait=FALSE))
          }else{
            if(file.exists(options()$browser))
              try(system(paste(options()$browser,fname),wait=FALSE))
            else if(1==length(system("which firefox",TRUE))) 
                try(system(paste("firefox",fname),wait=FALSE))
              else if(1==length(system("which mozilla",TRUE))) 
                   try(system(paste("mozilla",fname),wait=FALSE))
                 else if(1==length(system("which epiphany",TRUE))) 
                      try(system(paste("epiphany",fname),wait=FALSE))
          }
    }
    melde("ViewReport.html",2)
  }
  EditReport<-function(){
    melde("EditReport",1)
    if(substring(version$os,1,6)=="darwin"  ){
       cat("Sorry: not implemented yet!\n"); return()
    }    
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    worktext<-TcltoWin.write(worktext)
    try.res <-try(cat(worktext,file=get("tmp.file.name",envir=revive.sys)
,sep="\n"))
    if(is.function(try.res)){
      ok <- "OK"
    } else {
      if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
      if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
      if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
      if(!is.character(ok)) { ok <- "OK" }
    }
    if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
      ok<-FALSE
      error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
      cat(error.msg,"\n")
      if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
         cat("A warning message stopped the evaluation!",
               "If you want to\nevaluate the code anyway",
               "evaluate code by:\n>WarnEval<")
      # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
    } else { ok<-TRUE }


    if(ok){
      ##  try.res <- try(system(paste(editor.sys, ..tmp..)))
      cmd<-paste(editor.sys, get("tmp.file.name",envir=revive.sys)
)
      if((.Platform$OS.type=="windows")){
        try.res<-try(shell(cmd,wait=TRUE))
      }else{
        try.res<-try(system(cmd))
      }
      if(is.function(try.res)){
        ok <- "OK"
      } else {
        if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
        if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
        if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
        if(!is.character(ok)) { ok <- "OK" }
      }
      if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
        ok<-FALSE
        error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
        cat(error.msg,"\n")
        if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
           cat("A warning message stopped the evaluation!",
                 "If you want to\nevaluate the code anyway",
                 "evaluate code by:\n>WarnEval<")
        # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
      } else { ok<-TRUE }


      if(ok){
      try.res<-worktext<-try(myscan(get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE))
        if(is.function(try.res)){
          ok <- "OK"
        } else {
          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
          if(!is.character(ok)) { ok <- "OK" }
        }
        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
          ok<-FALSE
          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
          cat(error.msg,"\n")
          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
             cat("A warning message stopped the evaluation!",
                   "If you want to\nevaluate the code anyway",
                   "evaluate code by:\n>WarnEval<")
          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
        } else { ok<-TRUE }


        if(ok){
          if((.Platform$OS.type=="windows"))worktext<-WinToTcl.read(worktext)
          line <-floor(as.numeric(tkindex(tworkwin,"insert")))

          code.start<-grep("^<<(.*)>>=",worktext)
          try(if(0<length(code.start)){ 
                 worktext[code.start]<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext[code.start])
                 worktext[code.start]<-paste(worktext[code.start]," (",1:length(code.start),")",sep="")
          })
          if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
          tkdelete(tworkwin,"0.0","end")
          try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
          tksee(tworkwin,"end")
          melde("ak texthervor",1)
          tcl("markclear",tworkwin)
          tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
          tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
          tcl("marklinetypes",tworkwin)
          melde("ak texthervor",2)

          if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


          tkmark.set(tworkwin, "insert", paste(line,"0",sep="."))
          tksee(tworkwin,paste(line,"0",sep="."))
          tkfocus(tworkwin)

          worktext <- paste(worktext,collapse="\n")
  #       worktext <- sub("\\\\end\{document\}","%end\{document\}",worktext)
        } else {  cat("ERROR: file open operation not successful!!!\n") }
      } else {  cat("ERROR: editor start not successful!!!\n") }
    } else { cat("ERROR: file write operation not successful!!!\n") }
    melde("EditReport",2)
  }

  RunAll<-function(){
    melde("RunAll",1)
    news<-paste("RunAll:",date(),"\n")
    if(!exists("toutwin"))
      toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
    pos.to.insert<-"end"
    ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
    news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
    try(tkinsert(toutwin,pos.to.insert,news))
    tksee(toutwin,"end - 0 lines")
    melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

    news<-paste("\n@\n<<RunAllStartsAndStars>>=\n<<start>>\n<<*>>\n@")
    if(!exists("tworkwin"))
      tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

    pos.to.insert<-"end"
    if(0<length(grep("output-start",news))){
      tail<-rev(strsplit(tclvalue(tkget(tworkwin,"end - 3 lines","end")),"\n")[[1]])
      ltail<-length(tail)
      if( (0==length(grep("<<[*]>>=",tail[1:ltail]))) &&
         any(h<-("output-end"==substring(tail[1:ltail],1,11)))){
         news<-sub(".*output-start\n","",news)
         news<-sub("output-end","",news)
         h<-seq(along=h)[h][1]
         pos.to.insert<-paste("end -",h,"lines")
      }
    }
    try(tkinsert(tworkwin,pos.to.insert,paste(news,collapse="\n")))
    tksee(tworkwin,"end - 0 lines")
    melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

    tkmark.set(tworkwin,"insert","end")
    fWarnEval()
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    h<-grep("RunAllStartsAndStars",worktext)
    if(0<length(h)) worktext<-worktext[1:(h[1]-2)]
    if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
    tkdelete(tworkwin,"0.0","end")
    try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
    tksee(tworkwin,"end")
    melde("ak texthervor",1)
    tcl("markclear",tworkwin)
    tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
    tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
    tcl("marklinetypes",tworkwin)
    melde("ak texthervor",2)

    if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


    melde("RunAll",2)
  }

  RunAllandIncludeResults<-function(){
    melde("RunAllandIncludeResults",1)
    # get report
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

      # set default values
      defaults <- list()
      set.defaults <- function(echo=NULL,eval=NULL,results=NULL,fig=NULL,height=NULL,width=NULL,prefix.string,...){
           defaults <- defaults 
           if(!is.null(echo))     defaults$ECHO     <- echo
           if(!is.null(eval))     defaults$EVAL     <- eval
           if(!is.null(results))  defaults$RESULTS  <- results
           if(!is.null(fig))      defaults$FIG      <- fig
           if(!is.null(height))   defaults$HEIGHT   <- height
           if(!is.null(width))    defaults$WIDTH    <- width
           if(!missing(prefix.string)) {
             prefix.string<-as.character(substitute(prefix.string)) 
             #  prefix.string <- gsub("\\\\","/",prefix.string)
             defaults$PREFIX.STRING    <- prefix.string
           }
           defaults
      }
      ECHO <- TRUE; EVAL <- TRUE; RESULTS <- "verbatim"; FIG=FALSE; HEIGHT="10cm"; WIDTH=NA; PREFIX.STRING <- ""
      defaults <- set.defaults(echo=ECHO,eval=EVAL,results=RESULTS,fig=FIG,height=HEIGHT,width=WIDTH,
                               prefix.string="")
  
      # find names of code chunks  
      idx <- grep(paste("^<","<.*>",">=",sep=""),worktext)
      if(0==length(idx)) return("relax warning: no chunk found!")
      chunk.set <- worktext[idx] <- sub(">>=.*",">>=",worktext[idx])
      sexpr.lines<-grep("\\Sexpr\\{.*\\}",worktext)
      for(l in seq(along=sexpr.lines)){
        # hole Nummer l der Zeilen, die Sexpr-Expressions enthalten 
        cand<-worktext[sexpr.lines[l]]
        # knacke Kandidaten-Zeile an der Stelle auf, an der \Sexpr gefunden wird
        cand<-unlist(strsplit(cand,"\\\\Sexpr"))
        # cand[1] ist der vor der ersten Expression, 
        # cand[i+1] der mit der i-ten Expression beginnt
        # alle Expressions der Zeile werden nacheinander abgearbeitet
        for(j in seq(cand)[-1]){
          # ncandj zeigt die Laenge von Kandidat j an
          ncandj<-nchar(cand[j])
          # sexpr verwaltet den j-ten Kandidaten zeichenweise
          sexpr<-substring(cand[j],1:ncandj,1:ncandj) 
          # es gilt die beendende Klammer von Sexpr zu finden
          brack<-cumsum((sexpr=="{")-(sexpr=="}")) 
          # n.sexpr zeigt die Stelle der schliessenden-Klammer
          n.sexpr<-which(brack==0)[1]; if(is.na(n.sexpr)) next
          # mit n.sexpr greifen wir den vorderen Teil von sexpr und evaluieren
          code <- paste(collapse="",sexpr[1:n.sexpr])
          melde("werte Sexpr aus",2,code)
          # if(identical(revive.env,"")) ... else -> gegenüber weaveR vereinfacht
          result <- try(eval(parse(text=code),envir=revive.env))
          # wenn nichts rauskommt, ist nichts zu modifizieren
          if(0!=length(result)&&!identical(result,"")) { 
            # 101217 auch leere Ergebnisse ersetzen Sexpr!
            # print("---");print(result);print("---")
            # im Fehlerfall muss es eine Meldung geben
            if(class(result)=="try-error"){ 
              result<-paste("[[\\Sexpr-error:",
                            paste(sexpr[1:n.sexpr],collapse=""),"]]",collaspe="")
            }else{
              # bei nummerischen Ergebnissen werden ungewollte Nachkommastellen entfernt
              if(is.numeric(result)) result<-signif(result,digits=options()$digits)
              # Das Ergebnis wird verpackt
              result<-paste("[[",paste(unlist(result),collapse=" "),"]]",sep="")
            }
          }
          # das Ergebnis des j-ten Ausdrucks wird vorn,
          # also wo das Kommando stand eingetragen
          cand[j]<-paste(result, substring(cand[j],n.sexpr+1),sep="")
        }
        worktext[sexpr.lines[l]]<-paste(cand,collapse="")  
      }

      # initialize some variables
      pic.no <- 0; rep.name <- sub(".rev$","",workname.sys)
      for(i.chunk in seq(along=chunk.set)){     
         # get chunk header
         chunk.header <- chunk.set[i.chunk]; LINE <- match(chunk.set[i.chunk],worktext)
         # get code of chunk -- copied from EvalRCode
         code.start<-LINE
         if(LINE < length(worktext) && "@"!=substring(worktext[LINE+1],1,1)){
           if((code.start<-max(code.start)+1)>length(worktext)) return()
           code.end  <-c(grep("^@",worktext),1+length(worktext))
           code.end  <-min(code.end[code.end>code.start])-1
           code<-worktext[code.start:code.end]
           code<-code[code!=""]
           if(length(weg.ab<-grep("^<<(.*)>>=",code))>0) code<-code[-(weg.ab:length(code))]
           if(length(code)==0 || code[1]=="@") code<-" "
           melde("code:",3,code,"\n")

         } else { code <- NULL; code.end <- LINE }
         CODE.END <- code.end
         melde(paste(220,"code.start",code.start),2,code) 
         # set defaults new
         idx <- grep("\\\\SweaveOpts[{].*[}]",worktext[1:LINE])
         if(0<length(idx)){
           txt <- worktext[rev(idx)[1]]
           txt <- sub("^.*\\\\SweaveOpts[{](.*)[}].*$","set.defaults(\\1)",txt)
           defaults <- eval(parse(text=txt))
           for(i in seq(along=defaults)){
             if(is.character(defaults[[i]])) h<-"'" else h <- NULL
             eval(parse(text=paste(toupper(names(defaults))[i],"<-",h,defaults[[i]],h,sep="")))
           }
           if(!is.na(PREFIX.STRING)) rep.name <- paste(sep="",PREFIX.STRING,sub(".rev$","",workname.sys))
         }
         # set local options to default values
           for(i in seq(along=defaults)){
             if(is.character(defaults[[i]])) h<-"'" else h <- NULL
             eval(parse(text=paste(toupper(names(defaults))[i],"<-",h,defaults[[i]],h,sep="")))
           }
         # set local options 
         c.h <- gsub(" ","",chunk.header)
         if(0<length(grep("echo=TRUE", c.h))) ECHO <- TRUE
         if(0<length(grep("echo=FALSE",c.h))) ECHO <- FALSE
         if(0<length(grep("eval=TRUE", c.h))) EVAL <- TRUE
         if(0<length(grep("eval=FALSE",c.h))) EVAL <- FALSE
         if(0<length(grep("results=",c.h))){
           RESULTS <- sub("^.*results=([a-z]*).*$","\\1",c.h)
           if( !(RESULTS %in% c("verbatim","tex","hide"))) RESULTS <- "verbatim"
         }
         if(0<length(grep("fig=TRUE", c.h))) FIG <- TRUE
         if(0<length(grep("fig=FALSE",c.h))) FIG <- FALSE
         if(0<length(grep("height=",c.h))){
           HEIGHT <- sub("^.*height=([0-9.]+[a-z]*).*$","\\1",c.h)
         }
         if(0<length(grep("width=",c.h))){
           WIDTH  <- sub("^.*width=([0-9.]+[a-z]*).*$", "\\1",c.h)
         }
         if(0<length(grep("prefix.string=",c.h))){
           PREFIX.STRING  <- sub("^.*prefix.string=([a-z0-9.\\-]+).*$", "\\1",c.h)
         }
         # remove local option settings
         c.h <- chunk.header
         repeat{
           h<-sub("[,; ]*[a-z]+ *= *[a-z0-9A-Z]+","",c.h)
           if(h==c.h) break else c.h <- h
         }
         worktext[LINE] <- c.h
         # eval option "TRUE" found
         if(EVAL){
           EVAL.RESULT <- try.res <- NULL
           if(0<length(code)){
             melde("vor tangle - Code:\n",3,code)
             if(length(grep("<<(.*)>>",code))>0 || length(grep(">#",code))>0){
               code.a    <- grep("^<<(.*)>>=",worktext)
               code.z    <- grep("^@",worktext)
               code.z    <- unlist(sapply(code.a ,function(x,y) y[y>x][1], code.z))
               if(any(h<-is.na(code.z))) code.z<-code.z[!h]
               ################
               # code.n    <- length(worktext)
               # change    <- rep(0,code.n); change[c(code.a ,code.z)]<-1
               # code.ch   <- worktext[1==(cumsum(change)%%2)]
               ## 080311
               ind<-rep(0,length(worktext)); copy<-0
               for(i in seq(ind)){
                 if(i %in% code.z) copy<-0
                 if(i %in% code.a) copy<-1
                 ind[i]<-copy
               }
               code.ch   <- worktext[ind==1]
               ## cat("input-text, Anfaenge, Code.chunks"); print(worktext); print(code.a); print(code.ch)

               code.n    <- length(code.ch)

               code.ch<-gsub("@>>","DoSpCloseKl-esc",gsub("@<<","DoSpOpenKl-esc",code.ch))

               code.ch<-gsub("(.*)<<(.*)>>=(.*)","cOdEdEf\\2",code.ch)
               repeat{
                 if(0==length(cand<-grep("<<(.*)>>",code.ch))) break
                 code.ch<-unlist(strsplit(gsub("(.*)<<(.*)>>(.*)",
                            "\\1bReAkuSeChUnK\\2bReAk\\3",code.ch),"bReAk"))
               }
               code.ch<-code.ch[code.ch!=""]
               code.n<-length(code.ch)
               melde("code.ch:",3,code.ch,"\n")

               line.typ  <-rep("C",code.n)
               code.a    <-grep("cOdEdEf",code.ch)
               code.ch[code.a]<-substring(code.ch[code.a],8)
               line.typ[code.a]<-"D"
               code.use    <-grep("uSeChUnK",code.ch)
               code.ch[code.use]<-substring(code.ch[code.use],9)
               line.typ[code.use]<-"U"
               code.ext  <-grep("#<file",code.ch)
               line.typ[code.ext]<-"E"
               melde("code.ch:",3,code.ch,"\n")

               code.out<-"##act:##"

               def.names<-code.ch[code.a]
               use.names<- if(length(code.use)>0) code.ch[code.use] else NULL
               code.z<-c(if(length(code.a)>1) code.a[-1]-1, code.n)
               code.ch<-paste(line.typ,code.ch,sep="")
               melde("code.ch:",3,code.ch,"\n")

               melde("vor expand - Code:\n",3,code)
               melde("bearbeite aktuellen Chunk\n",3)
               ###<bestimme Cursorzeile [[line]] von [[tworkwin]]> not a good idea###
               ch.no<-length(grep("^<<(.*)>>=",worktext[1:line]))

               rows      <-c((code.a[ch.no]+1),code.z[ch.no])
               if(all(!is.na(rows))&&rows[1]<=rows[2]){
                 rows<-rows[1]:rows[2]
                 code.stack<-code.ch[rows]
                 max.depth.refinements<-500; i<-1
                 repeat{
                    if((i<-i+1)>max.depth.refinements){ 
                        cat("ERROR: maximal number of expandations (",max.depth.refinements,
                            ") exceeded\n --- perhaps a unintended recursion ???")
                        return()
                    }
                    if(0==length(code.stack))break
                    typ<-substring(code.stack[1],1,1)
                    if("C"==typ||"E"==typ){
                      n.lines<-sum(cumprod("C"==substring(code.stack,1,1)))
                      code.out<-c(code.out, substring(code.stack[1:n.lines],2))
                      code.stack<-code.stack[-(1:n.lines)]
                    }
                    if(length(code.stack)>0 && "U"==substring(code.stack[1],1,1)){
                      if(any(found<-def.names==substring(code.stack[1],2))){
                        found<-seq(along=def.names)[found]; rows<-NULL
                        for(no in found){
                          if((code.a[no]+1)<=code.z[no]) rows<-c(rows,(code.a[no]+1):code.z[no])
                        }
                        code.stack<-c(code.ch[rows],code.stack[-1])
                        melde("found",0,found)
                      }else{code.stack<-code.stack[-1]}
                    }

                 }
               }
               if(length(code.ext)>0){
                 code.out<-code.out[code.out!=""]
                 code.ext<-rev(grep(">#",code.out))
                 found<-TRUE
                 repeat{
                   if(length(code.ext)==0) break

                   if(!found){
                     code.out[code.ext[1]]<-paste("# ??",code.out[code.ext[1]])
                     cat("ERROR: External Chunk",code.out[code.ext[1]],"not found!!!\n")
                     code.ext<-code.ext[-1]
                   }

                   found<-TRUE
                   ext.name <- rev(unlist(strsplit(code.out[code.ext[1]],"#<file:")))[1]
                   ext.name <- unlist(strsplit(unlist(strsplit(ext.name,">#"))[1],":"))
                   ext.chunk<-ext.name[2]; ext.name <-ext.name[1]
                   ext.name.n<-nchar(ext.name)
                   if(ext.name.n >4 && ".rev"==substring(ext.name,ext.name.n-3,ext.name.n)){
                     ext.name<-substring(ext.name,1,ext.name.n-4)
                   }

                   if(is.na(as.numeric(ext.chunk))){
                     # tld untersuchen
                     filename<-paste(ext.name[1],".rev",sep="")
                     if(!file.exists(filename)){
                       cat("ERROR: file",filename,"for expansion of code chunk not found!!!\n")
                       ext.file<-"Error"
                     }else{
                       ext.file<-try(myscan(file=filename,what="",sep="\n"))
                     }
                     if("Error"==substring(unlist(ext.file)[1],1,5)){
                       found <-FALSE; next
                     }
                     ext.file <-ext.file[grep("^<<(.*)>>=",ext.file)]
                     if(!is.null(ext.file)){ found<-FALSE; next }
                     ext.chunk<-grep(ext.chunk,ext.file)
                   }

                   filename<-paste(ext.name[1],".R",sep="")
                   if(!file.exists(filename)){
                     cat("Warning: file",filename,"not found!!!\n")
                     cat("         file",filename,"is now generated!!!\n")
                     try(tangleR(ext.name[1]))
                   }
                   if(!file.exists(filename)){
                     ext.file<-"Error"
                   }else{
                     ext.file<-try(myscan(file=filename,what="",sep="\n"))
                   }
                   if("Error"==substring(unlist(ext.file)[1],1,5)){
                     found <-FALSE; next
                   }

                   ext.chunk<-as.numeric(ext.chunk)
                   a        <-grep(paste("#", ext.chunk,":",sep=""),ext.file)[1]
                   z        <-grep(paste("#:",ext.chunk,    sep=""),ext.file)[1]
                   if(is.na(a)){
                     found <- FALSE; next
                   }
                   if(a<=z) ext.file <-ext.file[a:z]

                   code.out <-c(code.out[1:(code.ext[1]-1)], ext.file,
                                if(length(code.out)>(code.ext[1]+1))
                                  code.out[(code.ext[1]+1):length(code.out)]
                              )

                   code.ext<-code.ext[-1]

                 }

               }
               code.out<-c(code.out,"##:act##")

               code.out<-gsub("DoSpCloseKl-esc",">>",gsub("DoSpOpenKl-esc","<<",code.out))

               melde("Ende Rtangle-last\n",3)
               code<-code.out[code.out!=""]
               melde("nach expand\n",3,code)
             }

             code <- sub("^ *print","base::print",code)
             code <- sub("^ *cat","base::cat",code)
             if(0<length(code)){
               tld <- paste("<","<*>",">=",sep="")
               rh <- c( get("relax.history",envir=revive.sys), list(c("@",tld,code)) )
               assign("relax.history",rh,envir=revive.sys)
             }

             code<-c("options(warn=2)",code)
             try.res <- try(eval(parse(text=code),envir=revive.env))
             options(warn=1)

             if(is.function(try.res)){
               ok <- "OK"
             } else {
               if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
               if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
               if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
               if(!is.character(ok)) { ok <- "OK" }
             }
             if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
               ok<-FALSE
               error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
               cat(error.msg,"\n")
               if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                  cat("A warning message stopped the evaluation!",
                        "If you want to\nevaluate the code anyway",
                        "evaluate code by:\n>WarnEval<")
               # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
             } else { ok<-TRUE }


             if(ok){
               EVAL.RESULT <- try.res
             } else { EVAL.RESULT <- "Error"; cat("sorry, evaluation not successful!!!\n") }
           } ## else { cat("no code found!!!\n") }
         } 

         if(FIG && EVAL){
         melde("225 ...Include picture",2) #<unused: Trash picture of report># 
           pic.no <- pic.no + 1
           bildname<-paste(rep.name,"-out-pic-",pic.no,".ps",sep="")
           if(!is.null(bildname)&&nchar(bildname)>0){
           # check name of picture
             n<-nchar(bildname<-gsub(" ","",bildname))
             bildname<-sub(".ps$","",bildname)
           # postscript:
             psname <-paste(bildname,".ps", sep="")
             try.res<-try({dev.copy(postscript,psname,horizontal=pshorizontal.sys,
                                    width=psdesignwidth.sys,height=psdesignheight.sys);dev.off()})
             if(is.function(try.res)){
               ok <- "OK"
             } else {
               if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
               if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
               if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
               if(!is.character(ok)) { ok <- "OK" }
             }
             if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
               ok<-FALSE
               error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
               cat(error.msg,"\n")
               if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                  cat("A warning message stopped the evaluation!",
                        "If you want to\nevaluate the code anyway",
                        "evaluate code by:\n>WarnEval<")
               # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
             } else { ok<-TRUE }


             if(!ok){ cat("Error: *ps file not generated by dev.copy!!!\n"); return() }
             news<-paste("@\n \\begin{center}","\\includegraphics[",
                         "height=",psheight.sys,"]{",bildname,"}\\end{center}\n",sep="") #081121
           # jpeg:
             jpgname<-paste(bildname,".jpg",sep="")
             if((.Platform$OS.type=="windows")){ # width=width in pixel, 72 dpi
               try.res<-try({dev.copy(jpeg,jpgname,width=jpgdesignsize.sys*72,
                                      height=jpgdesignsize.sys*72,quality=100,pointsize=7);dev.off()})
             }else{
               try.res<-try({dev.copy(bitmap,type="jpeg",jpgname,
                    width=jpgdesignsize.sys,height=jpgdesignsize.sys);dev.off()})
             }
             if(is.function(try.res)){
               ok <- "OK"
             } else {
               if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
               if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
               if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
               if(!is.character(ok)) { ok <- "OK" }
             }
             if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
               ok<-FALSE
               error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
               cat(error.msg,"\n")
               if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                  cat("A warning message stopped the evaluation!",
                        "If you want to\nevaluate the code anyway",
                        "evaluate code by:\n>WarnEval<")
               # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
             } else { ok<-TRUE }


             if(!ok) cat("Error: *jpg file not generated by dev.copy!!!\n")
             news<-paste(news,'\n% <p><img src="',jpgname,'">\n@\n', sep="" )
           # ppm
             ppmname<-paste(bildname,".ppm",sep="")
             # if( <das OS ist Windows>  && 2<=nchar(ghostscript) && !Img.package.found && !no.plots){
             if((.Platform$OS.type=="windows") && 2<=nchar(ghostscript) && !no.plots){
               try.res<-try({dev.copy(bitmap,type="ppmraw",ppmname,res=ppmresolution.sys);dev.off()})
             }
             # if( <das OS ist Linux> && 2<=nchar(ghostscript) && !Img.package.found && !no.plots){ }
             if(substring(version$os,1,5)=="linux" && 2<=nchar(ghostscript) && !no.plots){ # 121113
               try.res<-try({dev.copy(bitmap,type="ppmraw",ppmname,res=ppmresolution.sys);dev.off()})
             }
           # gif+ppm:
             if(substring(version$os,1,6)=="darwin"  && !no.plots){
               gifname<-paste(bildname,".gif",sep="")
               try({dev2bitmap(type="ppmraw",ppmname,res=ppmresolution.sys); dev.off()}) #121218
               try.res<-try({system(paste("convert",jpgname,gifname))})
               if(is.function(try.res)){
                 ok <- "OK"
               } else {
                 if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
                 if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
                 if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
                 if(!is.character(ok)) { ok <- "OK" }
               }
               if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
                 ok<-FALSE
                 error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
                 cat(error.msg,"\n")
                 if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                    cat("A warning message stopped the evaluation!",
                          "If you want to\nevaluate the code anyway",
                          "evaluate code by:\n>WarnEval<")
                 # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
               } else { ok<-TRUE }


               if(!ok) cat("Error: gif file not generated by dev.copy!!!\n")
             }
           }

           # include links
           n <- length(worktext)
           textstart <- c(grep("^@",worktext),n+1)
           textstart <- c(textstart[LINE <= textstart],n)[1]
           news <- unlist(strsplit(news,"\n")) 
           news <- news[grep("includegraphics",news)]
           if(!is.na(WIDTH)) news <- sub("cm",paste("cm,width=",WIDTH,sep=""),news)
           news <- sub("[0-9.]+cm",HEIGHT,news)
           worktext <-  c(worktext[  1:(textstart-1) ],"@",news, 
                          worktext[-(1:(textstart-1))])
         melde("226 ...Include picture",2,LINE,news,textstart)
         }

         melde("R-Output",2,RESULTS,EVAL,ECHO)
         if(EVAL && RESULTS != "hide"){ #<unused: Trash R Output>#
           if(0 < length(try.res <- EVAL.RESULT)){        
             if(!is.null(try.res)&&0<length(try.res)){
               sink(get("tmp.file.name",envir=revive.sys)
);get("print",pos="package:base")(try.res);sink()
               news<-paste(myscan(get("tmp.file.name",envir=revive.sys)
,"",sep="\n"),collapse="\n")
             } else news <- NULL
             melde("231...write result",2,news)
             if(1<=nchar(news)){
               if(length(grep("egin[{]table[}]",news))>0 &&
                  length(grep("generated.*xtable.*package",news))>0){
                 news<-sub("(\n%.*begin[{]table[}])","\noutput-end\n\\1",news)       
                 news<-sub("(.end[{]table[}])","\\1\noutput-start",news)       
                 news<-paste("\n@",paste("output-start",news,sep=""),"output-end\n", sep="\n")
                 news<-sub("output-start\n+output-end","",news)
               } else {    
                 if( RESULTS == "verbatim" ) 
                 # news<-paste("@","\\begin{quote}\\begin{verbatim}",news,"\\end{verbatim}\\end{quote}",sep="\n")
                   news<-paste("@","\\begin{verbatim}",news,"\\end{verbatim}",sep="\n")
                 if( RESULTS == "tex" ) news<-paste("@",news,sep="\n")
               }
               news <- gsub("\n+","\n",news)
               n <- length(worktext)
               textstart <- c(grep("^@",worktext),n+1)
               textstart <- c(textstart[LINE <= textstart],n)[1]
               worktext <- c(worktext[1:(textstart-1)],news, worktext[-(1:(textstart))])
             }
           }
         }

         # worktext[LINE] <- paste("<","<*>",">=",sep="")
         code.style <- "simple" # or "console" or "normal"
         if( LINE <= CODE.END ){
           if(!ECHO){
             worktext <- worktext[-(LINE:CODE.END)]
           } else {
             if(code.style == "console"){
               # remove code chunk name line and ">" followed by "+" in front of the code lines
               worktext[LINE] <- "\\begin{verbatim}"
               if(LINE   < CODE.END) worktext[LINE+1] <- paste(">",worktext[LINE+1],sep=" ")
               if(LINE+1 < CODE.END) {
                 idx <- (LINE+2):CODE.END
                 h <- paste("+",worktext[idx],sep=" "); h <- sub("^[+] *$"," ",h)
                 worktext[idx] <- h
               }
               worktext[CODE.END] <- paste(worktext[CODE.END],"\\end{verbatim}",sep="\n")
             }
             # no code chunk name and empty code chunk
             if( code.style == "simple" && 
                 0 < (h <- length(grep(paste("<","<[ \t]*>",">",sep=""),worktext[LINE]))) && LINE == CODE.END )
               code.style <- "simple console"
             # no code chunk name and no used code chunks within code
             if( code.style == "simple" && 0 < h && LINE < CODE.END &&
                 0 == length(grep(paste("<","<.*>",">",sep=""),worktext[(LINE+1):CODE.END])) ) 
               code.style <- "simple console"
             if(code.style == "simple console"){ 
               # remove code chunk name line and only one ">"
               worktext[LINE] <- "\\begin{verbatim}"
               if(LINE   < CODE.END) worktext[LINE+1] <- paste(">",worktext[LINE+1],sep=" ")
               if(LINE+1 < CODE.END) {
                 idx <- (LINE+2):CODE.END; worktext[idx] <- paste(" ",worktext[idx],sep=" ")
               }
               worktext[CODE.END] <- paste(worktext[CODE.END],"\\end{verbatim}",sep="\n")
             }
           }
         }
         # worktext <- worktext[-LINE]

      }
      idx <- grep("\\SweaveOpts",worktext)
      if(0 < length(idx)) {    
        worktext[idx] <- paste("\\newcommand{\\SweaveOpts}[1]{\\relax}",worktext[idx])
        if(1 < length(idx)) worktext[idx[-1]] <- sub("newcommand","renewcommand",worktext[idx[-1]])
      }
      filename <- paste(rep.name,"-out",".rev",sep="")
        #filename <- paste(rep.name,"-out",".tex",sep="")
        #<speichere [[worktext]] als Textfile ab, Indikator: [[ok]]>#
      base::cat(c(worktext,"@","\\end{document}"),file=filename,sep="\n")
      try(weaveR(filename)) # ,eval_Sexpr=FALSE)
      filename <- paste(rep.name,"-out",".tex",sep="")
      cat("-> '",filename,"' generated",sep="")
      Sys.sleep(0.5)
      ## LatexReport()
      LatexReport(filename=filename)


    melde("RunAllandIncludeResults",2)
  }

  RunBeginEndEnvandIncludeResults<-function(){
    melde("RunBeginEndEnvandIncludeResults",1)
    # from: ProcessBeginEndEnv
    filename <- "local-chunk.rev"
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      worktext <- sub("%.*$","",worktext)
      doc.start   <- grep("^.begin.document.",worktext)
      if(0 == length(doc.start)) { cat("relax warning: no begin document found"); return()}
      env.start <- grep("\\\\begin\\{",worktext) ; env.start <- env.start[env.start != doc.start] # } for symmetry
      env.end   <- grep("\\\\end\\{",worktext)  # } for symmetry
      idx.start <- idx.end <- rep(0,length(worktext))
      idx.start[env.start] <- 1; idx.end  [env.end  ] <- 1
      counter <- 0*idx.start; depth <- 0
      for( i in 1:length(counter)){
        counter[i] <- depth <- max(0, depth - idx.end[i] + idx.start[i])
      }
      #  counter <- cumsum(idx.start) - cumsum(idx.end) # fail if there are too much ends
      extract.start <- which(counter == 0 & seq(along=counter) < line)
      extract.start <- if(0 < length(extract.start)) max(extract.start) + 1 else doc.start + 1
      extract.end   <- which(counter == 0 & seq(along=counter) >= line)
      extract.end   <- if(0 < length(extract.start)) min(extract.end)       else length(worktext)  
      worktext <- c(worktext[1:doc.start],worktext[extract.start:extract.end],"\\end{document}")

    # from: RunAllandIncludeResults
      # set default values
      defaults <- list()
      set.defaults <- function(echo=NULL,eval=NULL,results=NULL,fig=NULL,height=NULL,width=NULL,prefix.string,...){
           defaults <- defaults 
           if(!is.null(echo))     defaults$ECHO     <- echo
           if(!is.null(eval))     defaults$EVAL     <- eval
           if(!is.null(results))  defaults$RESULTS  <- results
           if(!is.null(fig))      defaults$FIG      <- fig
           if(!is.null(height))   defaults$HEIGHT   <- height
           if(!is.null(width))    defaults$WIDTH    <- width
           if(!missing(prefix.string)) {
             prefix.string<-as.character(substitute(prefix.string)) 
             #  prefix.string <- gsub("\\\\","/",prefix.string)
             defaults$PREFIX.STRING    <- prefix.string
           }
           defaults
      }
      ECHO <- TRUE; EVAL <- TRUE; RESULTS <- "verbatim"; FIG=FALSE; HEIGHT="10cm"; WIDTH=NA; PREFIX.STRING <- ""
      defaults <- set.defaults(echo=ECHO,eval=EVAL,results=RESULTS,fig=FIG,height=HEIGHT,width=WIDTH,
                               prefix.string="")
  
      # find names of code chunks  
      idx <- grep(paste("^<","<.*>",">=",sep=""),worktext)
      if(0==length(idx)) return("relax warning: no chunk found!")
      chunk.set <- worktext[idx] <- sub(">>=.*",">>=",worktext[idx])
      sexpr.lines<-grep("\\Sexpr\\{.*\\}",worktext)
      for(l in seq(along=sexpr.lines)){
        # hole Nummer l der Zeilen, die Sexpr-Expressions enthalten 
        cand<-worktext[sexpr.lines[l]]
        # knacke Kandidaten-Zeile an der Stelle auf, an der \Sexpr gefunden wird
        cand<-unlist(strsplit(cand,"\\\\Sexpr"))
        # cand[1] ist der vor der ersten Expression, 
        # cand[i+1] der mit der i-ten Expression beginnt
        # alle Expressions der Zeile werden nacheinander abgearbeitet
        for(j in seq(cand)[-1]){
          # ncandj zeigt die Laenge von Kandidat j an
          ncandj<-nchar(cand[j])
          # sexpr verwaltet den j-ten Kandidaten zeichenweise
          sexpr<-substring(cand[j],1:ncandj,1:ncandj) 
          # es gilt die beendende Klammer von Sexpr zu finden
          brack<-cumsum((sexpr=="{")-(sexpr=="}")) 
          # n.sexpr zeigt die Stelle der schliessenden-Klammer
          n.sexpr<-which(brack==0)[1]; if(is.na(n.sexpr)) next
          # mit n.sexpr greifen wir den vorderen Teil von sexpr und evaluieren
          code <- paste(collapse="",sexpr[1:n.sexpr])
          melde("werte Sexpr aus",2,code)
          # if(identical(revive.env,"")) ... else -> gegenüber weaveR vereinfacht
          result <- try(eval(parse(text=code),envir=revive.env))
          # wenn nichts rauskommt, ist nichts zu modifizieren
          if(0!=length(result)&&!identical(result,"")) { 
            # 101217 auch leere Ergebnisse ersetzen Sexpr!
            # print("---");print(result);print("---")
            # im Fehlerfall muss es eine Meldung geben
            if(class(result)=="try-error"){ 
              result<-paste("[[\\Sexpr-error:",
                            paste(sexpr[1:n.sexpr],collapse=""),"]]",collaspe="")
            }else{
              # bei nummerischen Ergebnissen werden ungewollte Nachkommastellen entfernt
              if(is.numeric(result)) result<-signif(result,digits=options()$digits)
              # Das Ergebnis wird verpackt
              result<-paste("[[",paste(unlist(result),collapse=" "),"]]",sep="")
            }
          }
          # das Ergebnis des j-ten Ausdrucks wird vorn,
          # also wo das Kommando stand eingetragen
          cand[j]<-paste(result, substring(cand[j],n.sexpr+1),sep="")
        }
        worktext[sexpr.lines[l]]<-paste(cand,collapse="")  
      }

      # initialize some variables
      pic.no <- 0; rep.name <- sub(".rev$","",workname.sys)
      for(i.chunk in seq(along=chunk.set)){     
         # get chunk header
         chunk.header <- chunk.set[i.chunk]; LINE <- match(chunk.set[i.chunk],worktext)
         # get code of chunk -- copied from EvalRCode
         code.start<-LINE
         if(LINE < length(worktext) && "@"!=substring(worktext[LINE+1],1,1)){
           if((code.start<-max(code.start)+1)>length(worktext)) return()
           code.end  <-c(grep("^@",worktext),1+length(worktext))
           code.end  <-min(code.end[code.end>code.start])-1
           code<-worktext[code.start:code.end]
           code<-code[code!=""]
           if(length(weg.ab<-grep("^<<(.*)>>=",code))>0) code<-code[-(weg.ab:length(code))]
           if(length(code)==0 || code[1]=="@") code<-" "
           melde("code:",3,code,"\n")

         } else { code <- NULL; code.end <- LINE }
         CODE.END <- code.end
         melde(paste(220,"code.start",code.start),2,code) 
         # set defaults new
         idx <- grep("\\\\SweaveOpts[{].*[}]",worktext[1:LINE])
         if(0<length(idx)){
           txt <- worktext[rev(idx)[1]]
           txt <- sub("^.*\\\\SweaveOpts[{](.*)[}].*$","set.defaults(\\1)",txt)
           defaults <- eval(parse(text=txt))
           for(i in seq(along=defaults)){
             if(is.character(defaults[[i]])) h<-"'" else h <- NULL
             eval(parse(text=paste(toupper(names(defaults))[i],"<-",h,defaults[[i]],h,sep="")))
           }
           if(!is.na(PREFIX.STRING)) rep.name <- paste(sep="",PREFIX.STRING,sub(".rev$","",workname.sys))
         }
         # set local options to default values
           for(i in seq(along=defaults)){
             if(is.character(defaults[[i]])) h<-"'" else h <- NULL
             eval(parse(text=paste(toupper(names(defaults))[i],"<-",h,defaults[[i]],h,sep="")))
           }
         # set local options 
         c.h <- gsub(" ","",chunk.header)
         if(0<length(grep("echo=TRUE", c.h))) ECHO <- TRUE
         if(0<length(grep("echo=FALSE",c.h))) ECHO <- FALSE
         if(0<length(grep("eval=TRUE", c.h))) EVAL <- TRUE
         if(0<length(grep("eval=FALSE",c.h))) EVAL <- FALSE
         if(0<length(grep("results=",c.h))){
           RESULTS <- sub("^.*results=([a-z]*).*$","\\1",c.h)
           if( !(RESULTS %in% c("verbatim","tex","hide"))) RESULTS <- "verbatim"
         }
         if(0<length(grep("fig=TRUE", c.h))) FIG <- TRUE
         if(0<length(grep("fig=FALSE",c.h))) FIG <- FALSE
         if(0<length(grep("height=",c.h))){
           HEIGHT <- sub("^.*height=([0-9.]+[a-z]*).*$","\\1",c.h)
         }
         if(0<length(grep("width=",c.h))){
           WIDTH  <- sub("^.*width=([0-9.]+[a-z]*).*$", "\\1",c.h)
         }
         if(0<length(grep("prefix.string=",c.h))){
           PREFIX.STRING  <- sub("^.*prefix.string=([a-z0-9.\\-]+).*$", "\\1",c.h)
         }
         # remove local option settings
         c.h <- chunk.header
         repeat{
           h<-sub("[,; ]*[a-z]+ *= *[a-z0-9A-Z]+","",c.h)
           if(h==c.h) break else c.h <- h
         }
         worktext[LINE] <- c.h
         # eval option "TRUE" found
         if(EVAL){
           EVAL.RESULT <- try.res <- NULL
           if(0<length(code)){
             melde("vor tangle - Code:\n",3,code)
             if(length(grep("<<(.*)>>",code))>0 || length(grep(">#",code))>0){
               code.a    <- grep("^<<(.*)>>=",worktext)
               code.z    <- grep("^@",worktext)
               code.z    <- unlist(sapply(code.a ,function(x,y) y[y>x][1], code.z))
               if(any(h<-is.na(code.z))) code.z<-code.z[!h]
               ################
               # code.n    <- length(worktext)
               # change    <- rep(0,code.n); change[c(code.a ,code.z)]<-1
               # code.ch   <- worktext[1==(cumsum(change)%%2)]
               ## 080311
               ind<-rep(0,length(worktext)); copy<-0
               for(i in seq(ind)){
                 if(i %in% code.z) copy<-0
                 if(i %in% code.a) copy<-1
                 ind[i]<-copy
               }
               code.ch   <- worktext[ind==1]
               ## cat("input-text, Anfaenge, Code.chunks"); print(worktext); print(code.a); print(code.ch)

               code.n    <- length(code.ch)

               code.ch<-gsub("@>>","DoSpCloseKl-esc",gsub("@<<","DoSpOpenKl-esc",code.ch))

               code.ch<-gsub("(.*)<<(.*)>>=(.*)","cOdEdEf\\2",code.ch)
               repeat{
                 if(0==length(cand<-grep("<<(.*)>>",code.ch))) break
                 code.ch<-unlist(strsplit(gsub("(.*)<<(.*)>>(.*)",
                            "\\1bReAkuSeChUnK\\2bReAk\\3",code.ch),"bReAk"))
               }
               code.ch<-code.ch[code.ch!=""]
               code.n<-length(code.ch)
               melde("code.ch:",3,code.ch,"\n")

               line.typ  <-rep("C",code.n)
               code.a    <-grep("cOdEdEf",code.ch)
               code.ch[code.a]<-substring(code.ch[code.a],8)
               line.typ[code.a]<-"D"
               code.use    <-grep("uSeChUnK",code.ch)
               code.ch[code.use]<-substring(code.ch[code.use],9)
               line.typ[code.use]<-"U"
               code.ext  <-grep("#<file",code.ch)
               line.typ[code.ext]<-"E"
               melde("code.ch:",3,code.ch,"\n")

               code.out<-"##act:##"

               def.names<-code.ch[code.a]
               use.names<- if(length(code.use)>0) code.ch[code.use] else NULL
               code.z<-c(if(length(code.a)>1) code.a[-1]-1, code.n)
               code.ch<-paste(line.typ,code.ch,sep="")
               melde("code.ch:",3,code.ch,"\n")

               melde("vor expand - Code:\n",3,code)
               melde("bearbeite aktuellen Chunk\n",3)
               ###<bestimme Cursorzeile [[line]] von [[tworkwin]]> not a good idea###
               ch.no<-length(grep("^<<(.*)>>=",worktext[1:line]))

               rows      <-c((code.a[ch.no]+1),code.z[ch.no])
               if(all(!is.na(rows))&&rows[1]<=rows[2]){
                 rows<-rows[1]:rows[2]
                 code.stack<-code.ch[rows]
                 max.depth.refinements<-500; i<-1
                 repeat{
                    if((i<-i+1)>max.depth.refinements){ 
                        cat("ERROR: maximal number of expandations (",max.depth.refinements,
                            ") exceeded\n --- perhaps a unintended recursion ???")
                        return()
                    }
                    if(0==length(code.stack))break
                    typ<-substring(code.stack[1],1,1)
                    if("C"==typ||"E"==typ){
                      n.lines<-sum(cumprod("C"==substring(code.stack,1,1)))
                      code.out<-c(code.out, substring(code.stack[1:n.lines],2))
                      code.stack<-code.stack[-(1:n.lines)]
                    }
                    if(length(code.stack)>0 && "U"==substring(code.stack[1],1,1)){
                      if(any(found<-def.names==substring(code.stack[1],2))){
                        found<-seq(along=def.names)[found]; rows<-NULL
                        for(no in found){
                          if((code.a[no]+1)<=code.z[no]) rows<-c(rows,(code.a[no]+1):code.z[no])
                        }
                        code.stack<-c(code.ch[rows],code.stack[-1])
                        melde("found",0,found)
                      }else{code.stack<-code.stack[-1]}
                    }

                 }
               }
               if(length(code.ext)>0){
                 code.out<-code.out[code.out!=""]
                 code.ext<-rev(grep(">#",code.out))
                 found<-TRUE
                 repeat{
                   if(length(code.ext)==0) break

                   if(!found){
                     code.out[code.ext[1]]<-paste("# ??",code.out[code.ext[1]])
                     cat("ERROR: External Chunk",code.out[code.ext[1]],"not found!!!\n")
                     code.ext<-code.ext[-1]
                   }

                   found<-TRUE
                   ext.name <- rev(unlist(strsplit(code.out[code.ext[1]],"#<file:")))[1]
                   ext.name <- unlist(strsplit(unlist(strsplit(ext.name,">#"))[1],":"))
                   ext.chunk<-ext.name[2]; ext.name <-ext.name[1]
                   ext.name.n<-nchar(ext.name)
                   if(ext.name.n >4 && ".rev"==substring(ext.name,ext.name.n-3,ext.name.n)){
                     ext.name<-substring(ext.name,1,ext.name.n-4)
                   }

                   if(is.na(as.numeric(ext.chunk))){
                     # tld untersuchen
                     filename<-paste(ext.name[1],".rev",sep="")
                     if(!file.exists(filename)){
                       cat("ERROR: file",filename,"for expansion of code chunk not found!!!\n")
                       ext.file<-"Error"
                     }else{
                       ext.file<-try(myscan(file=filename,what="",sep="\n"))
                     }
                     if("Error"==substring(unlist(ext.file)[1],1,5)){
                       found <-FALSE; next
                     }
                     ext.file <-ext.file[grep("^<<(.*)>>=",ext.file)]
                     if(!is.null(ext.file)){ found<-FALSE; next }
                     ext.chunk<-grep(ext.chunk,ext.file)
                   }

                   filename<-paste(ext.name[1],".R",sep="")
                   if(!file.exists(filename)){
                     cat("Warning: file",filename,"not found!!!\n")
                     cat("         file",filename,"is now generated!!!\n")
                     try(tangleR(ext.name[1]))
                   }
                   if(!file.exists(filename)){
                     ext.file<-"Error"
                   }else{
                     ext.file<-try(myscan(file=filename,what="",sep="\n"))
                   }
                   if("Error"==substring(unlist(ext.file)[1],1,5)){
                     found <-FALSE; next
                   }

                   ext.chunk<-as.numeric(ext.chunk)
                   a        <-grep(paste("#", ext.chunk,":",sep=""),ext.file)[1]
                   z        <-grep(paste("#:",ext.chunk,    sep=""),ext.file)[1]
                   if(is.na(a)){
                     found <- FALSE; next
                   }
                   if(a<=z) ext.file <-ext.file[a:z]

                   code.out <-c(code.out[1:(code.ext[1]-1)], ext.file,
                                if(length(code.out)>(code.ext[1]+1))
                                  code.out[(code.ext[1]+1):length(code.out)]
                              )

                   code.ext<-code.ext[-1]

                 }

               }
               code.out<-c(code.out,"##:act##")

               code.out<-gsub("DoSpCloseKl-esc",">>",gsub("DoSpOpenKl-esc","<<",code.out))

               melde("Ende Rtangle-last\n",3)
               code<-code.out[code.out!=""]
               melde("nach expand\n",3,code)
             }

             code <- sub("^ *print","base::print",code)
             code <- sub("^ *cat","base::cat",code)
             if(0<length(code)){
               tld <- paste("<","<*>",">=",sep="")
               rh <- c( get("relax.history",envir=revive.sys), list(c("@",tld,code)) )
               assign("relax.history",rh,envir=revive.sys)
             }

             code<-c("options(warn=2)",code)
             try.res <- try(eval(parse(text=code),envir=revive.env))
             options(warn=1)

             if(is.function(try.res)){
               ok <- "OK"
             } else {
               if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
               if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
               if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
               if(!is.character(ok)) { ok <- "OK" }
             }
             if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
               ok<-FALSE
               error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
               cat(error.msg,"\n")
               if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                  cat("A warning message stopped the evaluation!",
                        "If you want to\nevaluate the code anyway",
                        "evaluate code by:\n>WarnEval<")
               # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
             } else { ok<-TRUE }


             if(ok){
               EVAL.RESULT <- try.res
             } else { EVAL.RESULT <- "Error"; cat("sorry, evaluation not successful!!!\n") }
           } ## else { cat("no code found!!!\n") }
         } 

         if(FIG && EVAL){
         melde("225 ...Include picture",2) #<unused: Trash picture of report># 
           pic.no <- pic.no + 1
           bildname<-paste(rep.name,"-out-pic-",pic.no,".ps",sep="")
           if(!is.null(bildname)&&nchar(bildname)>0){
           # check name of picture
             n<-nchar(bildname<-gsub(" ","",bildname))
             bildname<-sub(".ps$","",bildname)
           # postscript:
             psname <-paste(bildname,".ps", sep="")
             try.res<-try({dev.copy(postscript,psname,horizontal=pshorizontal.sys,
                                    width=psdesignwidth.sys,height=psdesignheight.sys);dev.off()})
             if(is.function(try.res)){
               ok <- "OK"
             } else {
               if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
               if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
               if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
               if(!is.character(ok)) { ok <- "OK" }
             }
             if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
               ok<-FALSE
               error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
               cat(error.msg,"\n")
               if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                  cat("A warning message stopped the evaluation!",
                        "If you want to\nevaluate the code anyway",
                        "evaluate code by:\n>WarnEval<")
               # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
             } else { ok<-TRUE }


             if(!ok){ cat("Error: *ps file not generated by dev.copy!!!\n"); return() }
             news<-paste("@\n \\begin{center}","\\includegraphics[",
                         "height=",psheight.sys,"]{",bildname,"}\\end{center}\n",sep="") #081121
           # jpeg:
             jpgname<-paste(bildname,".jpg",sep="")
             if((.Platform$OS.type=="windows")){ # width=width in pixel, 72 dpi
               try.res<-try({dev.copy(jpeg,jpgname,width=jpgdesignsize.sys*72,
                                      height=jpgdesignsize.sys*72,quality=100,pointsize=7);dev.off()})
             }else{
               try.res<-try({dev.copy(bitmap,type="jpeg",jpgname,
                    width=jpgdesignsize.sys,height=jpgdesignsize.sys);dev.off()})
             }
             if(is.function(try.res)){
               ok <- "OK"
             } else {
               if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
               if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
               if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
               if(!is.character(ok)) { ok <- "OK" }
             }
             if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
               ok<-FALSE
               error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
               cat(error.msg,"\n")
               if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                  cat("A warning message stopped the evaluation!",
                        "If you want to\nevaluate the code anyway",
                        "evaluate code by:\n>WarnEval<")
               # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
             } else { ok<-TRUE }


             if(!ok) cat("Error: *jpg file not generated by dev.copy!!!\n")
             news<-paste(news,'\n% <p><img src="',jpgname,'">\n@\n', sep="" )
           # ppm
             ppmname<-paste(bildname,".ppm",sep="")
             # if( <das OS ist Windows>  && 2<=nchar(ghostscript) && !Img.package.found && !no.plots){
             if((.Platform$OS.type=="windows") && 2<=nchar(ghostscript) && !no.plots){
               try.res<-try({dev.copy(bitmap,type="ppmraw",ppmname,res=ppmresolution.sys);dev.off()})
             }
             # if( <das OS ist Linux> && 2<=nchar(ghostscript) && !Img.package.found && !no.plots){ }
             if(substring(version$os,1,5)=="linux" && 2<=nchar(ghostscript) && !no.plots){ # 121113
               try.res<-try({dev.copy(bitmap,type="ppmraw",ppmname,res=ppmresolution.sys);dev.off()})
             }
           # gif+ppm:
             if(substring(version$os,1,6)=="darwin"  && !no.plots){
               gifname<-paste(bildname,".gif",sep="")
               try({dev2bitmap(type="ppmraw",ppmname,res=ppmresolution.sys); dev.off()}) #121218
               try.res<-try({system(paste("convert",jpgname,gifname))})
               if(is.function(try.res)){
                 ok <- "OK"
               } else {
                 if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
                 if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
                 if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
                 if(!is.character(ok)) { ok <- "OK" }
               }
               if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
                 ok<-FALSE
                 error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
                 cat(error.msg,"\n")
                 if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                    cat("A warning message stopped the evaluation!",
                          "If you want to\nevaluate the code anyway",
                          "evaluate code by:\n>WarnEval<")
                 # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
               } else { ok<-TRUE }


               if(!ok) cat("Error: gif file not generated by dev.copy!!!\n")
             }
           }

           # include links
           n <- length(worktext)
           textstart <- c(grep("^@",worktext),n+1)
           textstart <- c(textstart[LINE <= textstart],n)[1]
           news <- unlist(strsplit(news,"\n")) 
           news <- news[grep("includegraphics",news)]
           if(!is.na(WIDTH)) news <- sub("cm",paste("cm,width=",WIDTH,sep=""),news)
           news <- sub("[0-9.]+cm",HEIGHT,news)
           worktext <-  c(worktext[  1:(textstart-1) ],"@",news, 
                          worktext[-(1:(textstart-1))])
         melde("226 ...Include picture",2,LINE,news,textstart)
         }

         melde("R-Output",2,RESULTS,EVAL,ECHO)
         if(EVAL && RESULTS != "hide"){ #<unused: Trash R Output>#
           if(0 < length(try.res <- EVAL.RESULT)){        
             if(!is.null(try.res)&&0<length(try.res)){
               sink(get("tmp.file.name",envir=revive.sys)
);get("print",pos="package:base")(try.res);sink()
               news<-paste(myscan(get("tmp.file.name",envir=revive.sys)
,"",sep="\n"),collapse="\n")
             } else news <- NULL
             melde("231...write result",2,news)
             if(1<=nchar(news)){
               if(length(grep("egin[{]table[}]",news))>0 &&
                  length(grep("generated.*xtable.*package",news))>0){
                 news<-sub("(\n%.*begin[{]table[}])","\noutput-end\n\\1",news)       
                 news<-sub("(.end[{]table[}])","\\1\noutput-start",news)       
                 news<-paste("\n@",paste("output-start",news,sep=""),"output-end\n", sep="\n")
                 news<-sub("output-start\n+output-end","",news)
               } else {    
                 if( RESULTS == "verbatim" ) 
                 # news<-paste("@","\\begin{quote}\\begin{verbatim}",news,"\\end{verbatim}\\end{quote}",sep="\n")
                   news<-paste("@","\\begin{verbatim}",news,"\\end{verbatim}",sep="\n")
                 if( RESULTS == "tex" ) news<-paste("@",news,sep="\n")
               }
               news <- gsub("\n+","\n",news)
               n <- length(worktext)
               textstart <- c(grep("^@",worktext),n+1)
               textstart <- c(textstart[LINE <= textstart],n)[1]
               worktext <- c(worktext[1:(textstart-1)],news, worktext[-(1:(textstart))])
             }
           }
         }

         # worktext[LINE] <- paste("<","<*>",">=",sep="")
         code.style <- "simple" # or "console" or "normal"
         if( LINE <= CODE.END ){
           if(!ECHO){
             worktext <- worktext[-(LINE:CODE.END)]
           } else {
             if(code.style == "console"){
               # remove code chunk name line and ">" followed by "+" in front of the code lines
               worktext[LINE] <- "\\begin{verbatim}"
               if(LINE   < CODE.END) worktext[LINE+1] <- paste(">",worktext[LINE+1],sep=" ")
               if(LINE+1 < CODE.END) {
                 idx <- (LINE+2):CODE.END
                 h <- paste("+",worktext[idx],sep=" "); h <- sub("^[+] *$"," ",h)
                 worktext[idx] <- h
               }
               worktext[CODE.END] <- paste(worktext[CODE.END],"\\end{verbatim}",sep="\n")
             }
             # no code chunk name and empty code chunk
             if( code.style == "simple" && 
                 0 < (h <- length(grep(paste("<","<[ \t]*>",">",sep=""),worktext[LINE]))) && LINE == CODE.END )
               code.style <- "simple console"
             # no code chunk name and no used code chunks within code
             if( code.style == "simple" && 0 < h && LINE < CODE.END &&
                 0 == length(grep(paste("<","<.*>",">",sep=""),worktext[(LINE+1):CODE.END])) ) 
               code.style <- "simple console"
             if(code.style == "simple console"){ 
               # remove code chunk name line and only one ">"
               worktext[LINE] <- "\\begin{verbatim}"
               if(LINE   < CODE.END) worktext[LINE+1] <- paste(">",worktext[LINE+1],sep=" ")
               if(LINE+1 < CODE.END) {
                 idx <- (LINE+2):CODE.END; worktext[idx] <- paste(" ",worktext[idx],sep=" ")
               }
               worktext[CODE.END] <- paste(worktext[CODE.END],"\\end{verbatim}",sep="\n")
             }
           }
         }
         # worktext <- worktext[-LINE]

      }
      idx <- grep("\\SweaveOpts",worktext)
      if(0 < length(idx)) {    
        worktext[idx] <- paste("\\newcommand{\\SweaveOpts}[1]{\\relax}",worktext[idx])
        if(1 < length(idx)) worktext[idx[-1]] <- sub("newcommand","renewcommand",worktext[idx[-1]])
      }
      filename <- paste(rep.name,"-out",".rev",sep="")
        #filename <- paste(rep.name,"-out",".tex",sep="")
        #<speichere [[worktext]] als Textfile ab, Indikator: [[ok]]>#
      base::cat(c(worktext,"@","\\end{document}"),file=filename,sep="\n")
      try(weaveR(filename)) # ,eval_Sexpr=FALSE)
      filename <- paste(rep.name,"-out",".tex",sep="")
      cat("-> '",filename,"' generated",sep="")
      Sys.sleep(0.5)
      ## LatexReport()
      LatexReport(filename=filename)


    melde("RunBeginEndEnvandIncludeResults",2)
  }

  RunStart<-function(){
    melde("RunStart",1)
    news<-paste("RunStart:",date(),"\n")
    if(!exists("toutwin"))
      toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
    pos.to.insert<-"end"
    ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
    news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
    try(tkinsert(toutwin,pos.to.insert,news))
    tksee(toutwin,"end - 0 lines")
    melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

    news<-paste("\n@\n<<RunStarts>>=\n<<start>>\n@")
    if(!exists("tworkwin"))
      tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

    pos.to.insert<-"end"
    if(0<length(grep("output-start",news))){
      tail<-rev(strsplit(tclvalue(tkget(tworkwin,"end - 3 lines","end")),"\n")[[1]])
      ltail<-length(tail)
      if( (0==length(grep("<<[*]>>=",tail[1:ltail]))) &&
         any(h<-("output-end"==substring(tail[1:ltail],1,11)))){
         news<-sub(".*output-start\n","",news)
         news<-sub("output-end","",news)
         h<-seq(along=h)[h][1]
         pos.to.insert<-paste("end -",h,"lines")
      }
    }
    try(tkinsert(tworkwin,pos.to.insert,paste(news,collapse="\n")))
    tksee(tworkwin,"end - 0 lines")
    melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

    tkmark.set(tworkwin,"insert","end")
    fWarnEval()
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    h<-grep("RunStart",worktext)
    if(0<length(h)) worktext<-worktext[1:(rev(h)[1]-2)] ## 060310
    if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
    tkdelete(tworkwin,"0.0","end")
    try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
    tksee(tworkwin,"end")
    melde("ak texthervor",1)
    tcl("markclear",tworkwin)
    tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
    tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
    tcl("marklinetypes",tworkwin)
    melde("ak texthervor",2)

    if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


    melde("RunStart",2)
  }

  PLAYGROUND<-function() playground(revive.env)

  DeleteAll<-function(){
    melde("DeleteAll",1)
    if(!exists("tworkwin"))
      tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

    worktext<-TcltoWin.write(tclvalue(tkget(tworkwin,"0.0","end")))
    get("cat","package:base")(worktext,file="report-UnDo-bak.rev")

    worktext<-paste("% New Report:",date(),"\n")
    if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
    tkdelete(tworkwin,"0.0","end")
    try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
    tksee(tworkwin,"end")
    melde("ak texthervor",1)
    tcl("markclear",tworkwin)
    tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
    tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
    tcl("marklinetypes",tworkwin)
    melde("ak texthervor",2)

    if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


    melde("DeleteAll",2)
  }

  UnDo<-function(){
    melde("UnDo",1)
    if(file.exists("report-UnDo-bak.rev")){
      worktext<-myscan(file="report-UnDo-bak.rev",
                       what="",sep="\n",blank.lines.skip=FALSE)
      if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
      tkdelete(tworkwin,"0.0","end")
      try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
      tksee(tworkwin,"end")
      melde("ak texthervor",1)
      tcl("markclear",tworkwin)
      tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
      tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
      tcl("marklinetypes",tworkwin)
      melde("ak texthervor",2)

      if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


    }else{
      cat("ERROR: file report-UnDo-bak.rev not found!!!\n")
    }
    melde("UnDo",2)
  }

  ConvertEncodingToLocal<-function(){# Error: statt Umlaut ue U-tilde-Viertel
    melde("ConvertEncodingToLocal",1)
    frage<-"Input the old Encoding to be converted to locale encoding (e.g. macintosh):"; 
    set.tclvalue("tvinfo","macintosh")
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>", function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        set.tclvalue("tvmess","relax")
        from <- tclvalue("tvinfo")
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        ICONV <- function(x,from, to){ # iconv zeilenweise 
          local.encoding <- Encoding(iconv("\344","latin1",""))
          if(from=="") from <- local.encoding; if(to=="") to <- local.encoding
          org <- x; x <- sapply(x, function(x) iconv(x, from, to))
          x[is.na(x)] <- org[is.na(x)]; Encoding(x) <- to
          return(x)
        }

        worktext<-ICONV(worktext,from,"") # LATIN1
        line <-floor(as.numeric(tkindex(tworkwin,"insert")))

        code.start<-grep("^<<(.*)>>=",worktext)
        try(if(0<length(code.start)){ 
               worktext[code.start]<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext[code.start])
               worktext[code.start]<-paste(worktext[code.start]," (",1:length(code.start),")",sep="")
        })
        if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
        tkdelete(tworkwin,"0.0","end")
        try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
        tksee(tworkwin,"end")
        melde("ak texthervor",1)
        tcl("markclear",tworkwin)
        tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
        tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
        tcl("marklinetypes",tworkwin)
        melde("ak texthervor",2)

        if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


        tkmark.set(tworkwin, "insert", paste(line,"0",sep="."))
        tksee(tworkwin,paste(line,"0",sep="."))
        tkfocus(tworkwin)
        
    })
    melde("ConvertEncodingToLocal",2)
  }
  ConvertEncodingFromLocal<-function(){
    melde("ConvertEncodingFromLocal",1)
    frage<-"Input the new encoding (e.g. macintosh) to be applied:" 
    set.tclvalue("tvinfo","macintosh")
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>", function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        set.tclvalue("tvmess","relax")
        newencoding <- tclvalue("tvinfo")
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        ICONV <- function(x,from, to){ # iconv zeilenweise 
          local.encoding <- Encoding(iconv("\344","latin1",""))
          if(from=="") from <- local.encoding; if(to=="") to <- local.encoding
          org <- x; x <- sapply(x, function(x) iconv(x, from, to))
          x[is.na(x)] <- org[is.na(x)]; Encoding(x) <- to
          return(x)
        }

        try(worktext <- ICONV(worktext,"",newencoding))
        line <-floor(as.numeric(tkindex(tworkwin,"insert")))

        code.start<-grep("^<<(.*)>>=",worktext)
        try(if(0<length(code.start)){ 
               worktext[code.start]<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext[code.start])
               worktext[code.start]<-paste(worktext[code.start]," (",1:length(code.start),")",sep="")
        })
        if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
        tkdelete(tworkwin,"0.0","end")
        try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
        tksee(tworkwin,"end")
        melde("ak texthervor",1)
        tcl("markclear",tworkwin)
        tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
        tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
        tcl("marklinetypes",tworkwin)
        melde("ak texthervor",2)

        if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


        tkmark.set(tworkwin, "insert", paste(line,"0",sep="."))
        tksee(tworkwin,paste(line,"0",sep="."))
        tkfocus(tworkwin)
        
    })
    melde("ConvertEncodingFromLocal",2)
  }
  ConvertGermanUmlauteToLocalEncoding<-function(){
    melde("ConvertGermanUmlauteToLocalEncoding",1)
    ascii <- function(x) { strtoi(charToRaw(x), 16L) }
    ##lese Arbeitsfenster auf [[worktext]] ein## ohne splitting
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys); worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    .Tcl("set XYZ [encoding system]")
    UTF<-is.UTF<- 0<length(grep("utf",tclvalue("XYZ")))

    changes <- FALSE
    if(is.UTF){ melde("locale encoding is utf8:",3)
      idx <- grep("<[0-9a-f][0-9a-f]>",worktext)
      if( 0 < length(idx) ){ 
        melde("not-utf8-Umlaute found ",3)    
        x <- unlist(strsplit(worktext,"\n"))
        idx <- c(grep("<[0-9a-f][0-9a-f]>",x),grep("\\\\x[8-9ae]",x)) 
        xx <- x[idx] #; print(idx)
        uml <- iconv(c("\xe4","\xf6","\xfc","\xc4","\xd6","\xdc","\xdf"),"latin1","")
        ae<-uml[1];oe<-uml[2];ue<-uml[3];Ae<-uml[4];Oe<-uml[5];Ue<-uml[6];sz<-uml[7]
        # Konflikt: <9a>: dos-Ue==mac-oe => mac-oe wird zuerst ersetzt
        #              DOS:   MAC:   WIN:  hp-ux ( cp850, macintosh, latin1 )
        ae.codes <- c("<84>","<8a>","<e4>","<cc>"); for(code in ae.codes)xx <- gsub(code,ae,xx)
        oe.codes <- c("<94>","<9a>","<f6>","<ce>"); for(code in oe.codes)xx <- gsub(code,oe,xx)
        ue.codes <- c("<81>","<9f>","<fc>","<cf>"); for(code in ue.codes)xx <- gsub(code,ue,xx)
        Ae.codes <- c("<8e>","<80>","<c4>","<d8>"); for(code in Ae.codes)xx <- gsub(code,Ae,xx)
        Oe.codes <- c("<99>","<85>","<d6>","<da>"); for(code in Oe.codes)xx <- gsub(code,Oe,xx)
        # DOS Ue(x9a) will not be found because of MAC oe (x9a)
        Ue.codes <- c("<9a>","<86>","<dc>","<db>"); for(code in Ue.codes)xx <- gsub(code,Ue,xx)
        sz.codes <- c("<e1>","<a7>","<df>","<de>"); for(code in sz.codes)xx <- gsub(code,sz,xx)
        x[idx] <- xx; changes <- TRUE
        worktext <- paste(x,collapse="\n") #; print(worktext)      
      } else { # 121114
        if(substring(version$os,1,6)=="darwin"  && paste(R.version$major,R.version$minor,sep="")<211){
          # find Mac-Umlaute:
          umlloc <- iconv("[\xe4\xf6\xfc\xc4\xd6\xdc\xdf]","latin1","")
          # are there mac-Umlaute
          idx <- grep(umlloc,iconv(worktext,"macintosh",""))
          if(0 < length(idx)){
               x <- unlist(strsplit(worktext,""))
               # where are mac-Umlaute inclusive first charater of multichar representation
               # idx <- grep(umlloc,iconv(x,"macintosh","")); # alt
               idx <- grep(paste("\xc2\xac",umlloc,sep=""),iconv(x,"macintosh",""));
               # extract mac-Umlaute
               xx <- x[idx]
               # convert and repair mac-Umlaute
               xx <- iconv(xx,"macintosh",""); xx <- gsub("\xc2\xac","",xx)
               # repair errors
               if(any(is.na(xx))) xx[is.na(xx)] <- x[idx][is.na(xx)] 
               # save changes
               x[idx] <- xx; worktext <- x; worktext <- paste(worktext,collapse="") 
               changes <- TRUE
          } # end of replacements
        } # end of Mac-Mini
      } # end of not <..> case
    } # end of is.UTF
    is.LATIN1 <-l10n_info()[["Latin-1"]] 
    if(!is.UTF && is.LATIN1){ # 
      melde("locale encoding is Latin-1:",3)
      umlutf8 <- iconv("[\xe4\xf6\xfc\xc4\xd6\xdc\xdf]","latin1","utf8")
      idx <- grep(umlutf8,worktext)
      if( 0 < length(idx) ){ 
        melde("utf8-Umlaute found",3)    
        x <- unlist(strsplit(worktext,"\n")) 
        idx <- grep(umlutf8,x); xx <- x[idx] #; print(idx)      
        uml <-     iconv(c("\xe4","\xf6","\xfc","\xc4","\xd6","\xdc","\xdf"),"latin1","")
        umlutf8 <- iconv(c("\xe4","\xf6","\xfc","\xc4","\xd6","\xdc","\xdf"),"latin1","utf8")
        for(i in 1:7) xx <- gsub(umlutf8[i],uml[i],xx) #; cat("xx-utf-ok"); print(xx)      
        x[idx] <- xx; changes <- TRUE
        worktext <- paste(x,collapse="\n") #; print(worktext)      
      }
      # transform NOT-utf8 Umlaute: --------------------------------------------------
      umldosmac <- paste("[",iconv("\xe4\xf6\xfc\xc4\xd6\xdc\xdf","latin1","cp850"),
                             iconv("\xe4\xf6\xfc\xc4\xd6\xdc\xdf","latin1","macintosh"),"]",sep="")
      idx <- grep(umldosmac,worktext)
      if( 0 < length(idx) ){ 
        melde("dos- or mac-Umlaute found",3)    
        x <- unlist(strsplit(worktext,"\n"))
        idx <- grep(umldosmac,x); xx <- x[idx] #; print(idx)
        uml    <- iconv(c("\xe4","\xf6","\xfc","\xc4","\xd6","\xdc","\xdf"),"latin1","")
        umlmac <- iconv(c("\xe4","\xf6","\xfc","\xc4","\xd6","\xdc","\xdf"),"latin1","macintosh")
        umldos <- iconv(c("\xe4","\xf6","\xfc","\xc4","\xd6","\xdc","\xdf"),"latin1","cp850")
        # Konflikt: <9a>: dos-Ue==mac-oe => mac-umlaute werden zuerst ersetzt
        for(i in 1:7) xx <- gsub(umlmac[i],uml[i],xx) #; cat("xx-mac-ok"); print(xx)      
        for(i in 1:7) xx <- gsub(umldos[i],uml[i],xx) #; cat("xx-dos-ok"); print(xx)      
        x[idx] <- xx; changes <- TRUE
        worktext <- paste(x,collapse="\n") #; print(worktext)      
      }
    }
    if(changes){
        if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
        tkdelete(tworkwin,"0.0","end")
        try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
        tksee(tworkwin,"end")
        melde("ak texthervor",1)
        tcl("markclear",tworkwin)
        tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
        tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
        tcl("marklinetypes",tworkwin)
        melde("ak texthervor",2)

        if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }

    
    }
    melde("ConvertGermanUmlauteToLocalEncoding",2)
  } 

  DefineFKeys<-function(){
    melde("DefineFKeys",1)
    frage<-"number of function key?"; set.tclvalue("tvinfo",1)
    tkconfigure(linfo.tmp,text=frage)
    # tkpack("forget",linfo.name,linfo); Sys.sleep(0.01)
    tkpack("forget",linfo); Sys.sleep(0.01)
    tkpack(linfo.tmp,einfo.tmp,side="left"); Sys.sleep(0.01)
    tkfocus(einfo.tmp)
    tkselection.range(einfo.tmp,"0","end") ## 051219
    tkbind(TopW,"<Escape>",function(){
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


      }
    )


    tkbind(TopW,"<Return>", function(){
        F.no<-tclvalue("tvinfo")
        tkbind(TopW,"<Return>","")
        tkpack("forget",einfo.tmp,linfo.tmp); Sys.sleep(0.01)
        # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
        tkpack(linfo,side="left",fill="x",expand="yes")


        if(!any(F.no==as.character(1:8))) return() else F.no<-as.numeric(F.no)
        .newl<-tktoplevel();tkwm.geometry(.newl,"+0+15");tkpack(tt<-tktext(.newl))
        twin<-tt; f.sonderzeichen<-function(zeichen){
                    function(){
                      tkinsert(twin,"insert",zeichen)
                      if((.Platform$OS.type=="windows"))tkdelete(twin,"insert-1chars")
                      if(substring(version$os,1,6)=="darwin" )tkdelete(twin,"insert-1chars")
                    }
                  }
                  f.umlaut<-function(zeichen){
                    function(){
                      return()
                      #char337<-eval(parse(text='"\\337"'))
                      #if(zeichen==char337 & tclvalue(tkget(twin,"insert-1chars","insert"))=="\\") return()
                      #tkinsert(twin,"insert",zeichen); tkdelete(twin,"insert-2chars")
                    }
                  }
                  tkbind(twin,"<<LKeckig>>", f.sonderzeichen("["))
                  tkbind(twin,"<<RKeckig>>", f.sonderzeichen("]"))
                  tkbind(twin,"<<Tilde>>",   f.sonderzeichen("~"))
                  tkbind(twin,"<<LKgeschw>>",f.sonderzeichen("{"))
                  tkbind(twin,"<<RKgeschw>>",f.sonderzeichen("}"))
                  tkbind(twin,"<<Klammera>>",f.sonderzeichen("@"))
                  tkbind(twin,"<<Pipe>>",    f.sonderzeichen("|"))
                  tkbind(twin,"<<Backsl>>",  f.sonderzeichen("\\"))
                  renewhighlighting<-function(){
                    tworkwin<-get("tworkwin",envir=revive.sys)
                    melde("ak texthervor",1)
                    tcl("markclear",tworkwin)
                    tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
                    tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
                    tcl("marklinetypes",tworkwin)
                    melde("ak texthervor",2)

                  }
                  # tkbind(twin,"<<Klammeraffe>>",renewhighlighting)
                  tkbind(twin,"<Return>",renewhighlighting)

        tkwm.title(.newl,paste("contents of function key","Definition by Escape!"))
        F.buffers<-get("F.buffers",envir=revive.sys)
        try(tkinsert(tt,"0.0",paste(F.buffers[F.no],collapse="\n")))
        abbruch<-function(){tkdestroy(.newl); set.tclvalue("tvscandone",2)}
        tkbind(.newl,"<Escape>", function(){
                   F.buffers[F.no]<-tclvalue(tkget(tt,"0.0","end"))
                   assign("F.buffers",F.buffers,envir=revive.sys)
                   tkdestroy(.newl); set.tclvalue("tvscandone",2)
                 })
        tkfocus(.newl);tkwait.variable("tvndone")
        set.tclvalue("tvmess","relax")
      } # end of function
    )
    melde("DefineFKeys",2)
  }
    
          ##setze Umgebung für Knopf\-funktionen##
          ##setze Umgebung für Testknopf\-funktion## 
          ##setze Umgebung für Zeilen-Funktionen##
  environment(myhead.menu)<-revive.sys

  loadwwwdata<-function(){
    GetWWWData<-function(){
      melde("GetWWWData",1)
      URL<-paste("http://www.wiwi.uni-bielefeld.de/fileadmin/stat/wolf/data")
      download.file(paste(URL,"00Contents",sep="/"),"r.tmp")
      choices<-scan(file="r.tmp","",sep="\n")
      listboxmenu<-function(choices,title="items",addexit=TRUE){
          if(addexit) choices<-c(choices,"EXIT")
          lbmTop<-tktoplevel();tkwm.geometry(lbmTop,"+0+15")
          tkwm.title(lbmTop,"Selection by Click and Return, Quit by Esc")
          ltit<-tklabel(lbmTop,text=title); lbframe<-tkframe(lbmTop)
          lb<-tklistbox(lbframe,height=8,width=60); sb<-tkscrollbar(lbframe)
          tkconfigure(lb,yscrollcommand=function(...) tkset(sb,...))
          tkconfigure(sb ,command=function(...) tkyview(lb,...))
          tkpack(ltit,lbframe)
          tkpack(sb,side="right",fill="y"); tkpack(lb)
          for(i in seq(choices)) tkinsert(lb,"end",choices[i])
          lbmdone      <- tclVar()
          tkbind(lbmTop,"<Return>",function()set.tclvalue(lbmdone,"1"))
          tkbind(lbmTop,"<Escape>",function()set.tclvalue(lbmdone,"0"))
          tkfocus(lb); tkwait.variable(lbmdone)
          choice<-tclvalue(tkget(lb,"active")); tkdestroy(lbmTop)
          if(tclvalue(lbmdone)=="0") return(0)
          ind <- match(choice, choices)
          choice <- if(addexit && ind==length(choices)) 0 else ind
          return(choice)
      }

      choice<-listboxmenu(choices,"gefundene Datensaetze")
      if(0==choice) return() else choice<-sub(" .*$",".R",choices[choice])
      download.file(paste(URL,choice,sep="/"),"www.data.R")
      cmds<-scan(file="www.data.R","",sep="\n")
      ok<-try(eval(parse(text=cmds),envir=revive.env))
      cat("\"",sub(".R$","",choice),"\" geladen\n",sep="")
      melde("GetWWWData",2)
    }
    GetWWWData()
  }

  generate.1.data<-function(){
   (function(){
    DS<-c("sample(99)","sample(99,replace=TRUE)","rnorm(99,mean=0,sd=1)","rexp(99)")
    set.tclvalue<-function(name,value)tclvalue(name)<-as.character(value) # 121217
    listboxmenu<-function(choices,title="items",addexit=TRUE){
        if(addexit) choices<-c(choices,"EXIT")
        lbmTop<-tktoplevel();tkwm.geometry(lbmTop,"+0+15")
        tkwm.title(lbmTop,"Selection by Click and Return, Quit by Esc")
        ltit<-tklabel(lbmTop,text=title); lbframe<-tkframe(lbmTop)
        lb<-tklistbox(lbframe,height=8,width=60); sb<-tkscrollbar(lbframe)
        tkconfigure(lb,yscrollcommand=function(...) tkset(sb,...))
        tkconfigure(sb ,command=function(...) tkyview(lb,...))
        tkpack(ltit,lbframe)
        tkpack(sb,side="right",fill="y"); tkpack(lb)
        for(i in seq(choices)) tkinsert(lb,"end",choices[i])
        lbmdone      <- tclVar()
        tkbind(lbmTop,"<Return>",function()set.tclvalue(lbmdone,"1"))
        tkbind(lbmTop,"<Escape>",function()set.tclvalue(lbmdone,"0"))
        tkfocus(lb); tkwait.variable(lbmdone)
        choice<-tclvalue(tkget(lb,"active")); tkdestroy(lbmTop)
        if(tclvalue(lbmdone)=="0") return(0)
        ind <- match(choice, choices)
        choice <- if(addexit && ind==length(choices)) 0 else ind
        return(choice)
    }

    choice<-listboxmenu(DS,"Items:")
    if(0==length(choice)||is.na(choice)||0==choice) return()
    choice<-DS[choice]
    if(!exists("myscan")) myscan<-get("myscan",envir=revive.sys)



    cmd<-paste("x <-",choice)
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    melde <- get("melde", envir=revive.sys) # 121217
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    last.text<-c(0,grep("^@",worktext)); last.text<-rev(last.text[last.text<=line])[1]
    last.code<-grep("^<<(.*)>>=",worktext)
    if(0<length(last.code)){
      last.code<-rev(last.code[last.code<=line])[1]
      if(is.na(last.code)||last.code<last.text)
                {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="") }
      else
                      {delta<-1; news<-paste(cmd,"\n",sep="")}
    }else{
                {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="")}
    }
    if(!exists("tworkwin"))
      tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

        anzrows<-length(news)
        try(tkinsert(tworkwin,paste(line+1,"0",sep="."),paste(news,collapse="\n")))
        tkmark.set(tworkwin, "insert", paste(line+anzrows,"0",sep="."))
        tksee(tworkwin,paste(line+anzrows,"0",sep="."))
    tkfocus(tworkwin)
    melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

    tkmark.set(tworkwin, "insert", paste(line+delta,"0",sep="."))

   })()
  }

  get.dim.1.data<-function(){
   (function(){
    ds<-function(pos=1,type,mode="numeric",struc=is.vector){
      if(pos>0){
       obj<-ls(pos=pos)
       obj<-obj[unlist(lapply(obj ,function(o,pos)
         exists(o,where=pos,mode="numeric") &&
         eval(parse(text=paste("struc(get(\"",o,"\",pos=",pos,"))",sep=""))),pos))]
      }else{
       obj<-ls(envir=revive.env)
       obj<-obj[unlist(lapply(obj ,function(o)
         exists(o,where=revive.env,mode="numeric") &&
         eval(parse(text=paste("struc(get(\"",o,"\",envir=revive.env))",sep=""))) ))]
      }
      if(0==length(obj)) return(NULL) else return(obj)
    }
    ds.of.R<-function(type="vector"){
      dat<-ls(pos=grep("datasets",search()))
      dat.type<-unlist(lapply(dat,function(x) {       
             num<-mode(x<-eval(parse(text=x)))
             num<-ifelse(is.array(x),"array",num)
             num<-ifelse(is.list(x),"list",num)
             num<-ifelse(is.matrix(x),"matrix",num)
             num<-ifelse(is.data.frame(x),"matrix",num)
             num<-ifelse(num=="numeric","vector",num)
             num }))
      return(dat[dat.type==type])
    }



    DS<-c( ds.of.R("vector"),ds(-1), ds(1),
               ds(which(paste("package","relax"
,sep=":")==search())
),
               "-5:5 # integers from -5 to 5", "rep(7,5) # vector: (7,7,7,7,7)")

    set.tclvalue<-function(name,value)tclvalue(name)<-as.character(value) # 121217
    listboxmenu<-function(choices,title="items",addexit=TRUE){
        if(addexit) choices<-c(choices,"EXIT")
        lbmTop<-tktoplevel();tkwm.geometry(lbmTop,"+0+15")
        tkwm.title(lbmTop,"Selection by Click and Return, Quit by Esc")
        ltit<-tklabel(lbmTop,text=title); lbframe<-tkframe(lbmTop)
        lb<-tklistbox(lbframe,height=8,width=60); sb<-tkscrollbar(lbframe)
        tkconfigure(lb,yscrollcommand=function(...) tkset(sb,...))
        tkconfigure(sb ,command=function(...) tkyview(lb,...))
        tkpack(ltit,lbframe)
        tkpack(sb,side="right",fill="y"); tkpack(lb)
        for(i in seq(choices)) tkinsert(lb,"end",choices[i])
        lbmdone      <- tclVar()
        tkbind(lbmTop,"<Return>",function()set.tclvalue(lbmdone,"1"))
        tkbind(lbmTop,"<Escape>",function()set.tclvalue(lbmdone,"0"))
        tkfocus(lb); tkwait.variable(lbmdone)
        choice<-tclvalue(tkget(lb,"active")); tkdestroy(lbmTop)
        if(tclvalue(lbmdone)=="0") return(0)
        ind <- match(choice, choices)
        choice <- if(addexit && ind==length(choices)) 0 else ind
        return(choice)
    }

    choice<-listboxmenu(DS,"Items:")
    if(0==length(choice)||is.na(choice)||0==choice) return()
    choice<-DS[choice]
    if(!exists("myscan")) myscan<-get("myscan",envir=revive.sys)



    cmd<-paste("x <-",choice)
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    melde <- get("melde", envir=revive.sys) # 121217
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    last.text<-c(0,grep("^@",worktext)); last.text<-rev(last.text[last.text<=line])[1]
    last.code<-grep("^<<(.*)>>=",worktext)
    if(0<length(last.code)){
      last.code<-rev(last.code[last.code<=line])[1]
      if(is.na(last.code)||last.code<last.text)
                {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="") }
      else
                      {delta<-1; news<-paste(cmd,"\n",sep="")}
    }else{
                {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="")}
    }
    if(!exists("tworkwin"))
      tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

        anzrows<-length(news)
        try(tkinsert(tworkwin,paste(line+1,"0",sep="."),paste(news,collapse="\n")))
        tkmark.set(tworkwin, "insert", paste(line+anzrows,"0",sep="."))
        tksee(tworkwin,paste(line+anzrows,"0",sep="."))
    tkfocus(tworkwin)
    melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

    tkmark.set(tworkwin, "insert", paste(line+delta,"0",sep="."))

   })()
  }

  get.dim.2.data<-function(){
   (function(){
    ds<-function(pos=1,type,mode="numeric",struc=is.vector){
      if(pos>0){
       obj<-ls(pos=pos)
       obj<-obj[unlist(lapply(obj ,function(o,pos)
         exists(o,where=pos,mode="numeric") &&
         eval(parse(text=paste("struc(get(\"",o,"\",pos=",pos,"))",sep=""))),pos))]
      }else{
       obj<-ls(envir=revive.env)
       obj<-obj[unlist(lapply(obj ,function(o)
         exists(o,where=revive.env,mode="numeric") &&
         eval(parse(text=paste("struc(get(\"",o,"\",envir=revive.env))",sep=""))) ))]
      }
      if(0==length(obj)) return(NULL) else return(obj)
    }
    ds.of.R<-function(type="vector"){
      dat<-ls(pos=grep("datasets",search()))
      dat.type<-unlist(lapply(dat,function(x) {       
             num<-mode(x<-eval(parse(text=x)))
             num<-ifelse(is.array(x),"array",num)
             num<-ifelse(is.list(x),"list",num)
             num<-ifelse(is.matrix(x),"matrix",num)
             num<-ifelse(is.data.frame(x),"matrix",num)
             num<-ifelse(num=="numeric","vector",num)
             num }))
      return(dat[dat.type==type])
    }



    DS<-c(ds.of.R("matrix"),ds(struc=is.matrix))
    DS<-c(DS,"matrix(sample(100),50,2)","cbind(1:100,rnorm(100,mean=0,sd=1))")
    DS<-c(DS, ds(-1,struc=is.matrix), ds(1,struc=is.matrix) ,
               ds(which(paste("package","relax"
,sep=":")==search())
 ,struc=is.matrix) )

    set.tclvalue<-function(name,value)tclvalue(name)<-as.character(value) # 121217
    listboxmenu<-function(choices,title="items",addexit=TRUE){
        if(addexit) choices<-c(choices,"EXIT")
        lbmTop<-tktoplevel();tkwm.geometry(lbmTop,"+0+15")
        tkwm.title(lbmTop,"Selection by Click and Return, Quit by Esc")
        ltit<-tklabel(lbmTop,text=title); lbframe<-tkframe(lbmTop)
        lb<-tklistbox(lbframe,height=8,width=60); sb<-tkscrollbar(lbframe)
        tkconfigure(lb,yscrollcommand=function(...) tkset(sb,...))
        tkconfigure(sb ,command=function(...) tkyview(lb,...))
        tkpack(ltit,lbframe)
        tkpack(sb,side="right",fill="y"); tkpack(lb)
        for(i in seq(choices)) tkinsert(lb,"end",choices[i])
        lbmdone      <- tclVar()
        tkbind(lbmTop,"<Return>",function()set.tclvalue(lbmdone,"1"))
        tkbind(lbmTop,"<Escape>",function()set.tclvalue(lbmdone,"0"))
        tkfocus(lb); tkwait.variable(lbmdone)
        choice<-tclvalue(tkget(lb,"active")); tkdestroy(lbmTop)
        if(tclvalue(lbmdone)=="0") return(0)
        ind <- match(choice, choices)
        choice <- if(addexit && ind==length(choices)) 0 else ind
        return(choice)
    }

    choice<-listboxmenu(DS,"Items:")
    if(0==length(choice)||is.na(choice)||0==choice) return()
    choice<-DS[choice]
    if(!exists("myscan")) myscan<-get("myscan",envir=revive.sys)



    cmd<-paste("xy <-",choice)
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    melde <- get("melde", envir=revive.sys) # 121217
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    last.text<-c(0,grep("^@",worktext)); last.text<-rev(last.text[last.text<=line])[1]
    last.code<-grep("^<<(.*)>>=",worktext)
    if(0<length(last.code)){
      last.code<-rev(last.code[last.code<=line])[1]
      if(is.na(last.code)||last.code<last.text)
                {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="") }
      else
                      {delta<-1; news<-paste(cmd,"\n",sep="")}
    }else{
                {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="")}
    }
    if(!exists("tworkwin"))
      tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

        anzrows<-length(news)
        try(tkinsert(tworkwin,paste(line+1,"0",sep="."),paste(news,collapse="\n")))
        tkmark.set(tworkwin, "insert", paste(line+anzrows,"0",sep="."))
        tksee(tworkwin,paste(line+anzrows,"0",sep="."))
    tkfocus(tworkwin)
    melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

    tkmark.set(tworkwin, "insert", paste(line+delta,"0",sep="."))

   })()
  }

  removenafromvecx<-function(){
    (function(){
      cmd<-"x<-x[!is.na(x)]"
      if(!exists("myscan")) myscan<-get("myscan",envir=revive.sys)


      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      melde <- get("melde", envir=revive.sys) # 121217
      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      last.text<-c(0,grep("^@",worktext)); last.text<-rev(last.text[last.text<=line])[1]
      last.code<-grep("^<<(.*)>>=",worktext)
      if(0<length(last.code)){
        last.code<-rev(last.code[last.code<=line])[1]
        if(is.na(last.code)||last.code<last.text)
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="") }
        else
                        {delta<-1; news<-paste(cmd,"\n",sep="")}
      }else{
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="")}
      }
      if(!exists("tworkwin"))
        tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

          anzrows<-length(news)
          try(tkinsert(tworkwin,paste(line+1,"0",sep="."),paste(news,collapse="\n")))
          tkmark.set(tworkwin, "insert", paste(line+anzrows,"0",sep="."))
          tksee(tworkwin,paste(line+anzrows,"0",sep="."))
      tkfocus(tworkwin)
      melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      tkmark.set(tworkwin, "insert", paste(line+delta,"0",sep="."))

     })()
  }
  removenafrommatxy<-function(){
    (function(){
      cmd<-"xy<-xy[!apply(is.na(xy),1,any),,drop=FALSE]"
      if(!exists("myscan")) myscan<-get("myscan",envir=revive.sys)


      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      melde <- get("melde", envir=revive.sys) # 121217
      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      last.text<-c(0,grep("^@",worktext)); last.text<-rev(last.text[last.text<=line])[1]
      last.code<-grep("^<<(.*)>>=",worktext)
      if(0<length(last.code)){
        last.code<-rev(last.code[last.code<=line])[1]
        if(is.na(last.code)||last.code<last.text)
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="") }
        else
                        {delta<-1; news<-paste(cmd,"\n",sep="")}
      }else{
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="")}
      }
      if(!exists("tworkwin"))
        tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

          anzrows<-length(news)
          try(tkinsert(tworkwin,paste(line+1,"0",sep="."),paste(news,collapse="\n")))
          tkmark.set(tworkwin, "insert", paste(line+anzrows,"0",sep="."))
          tksee(tworkwin,paste(line+anzrows,"0",sep="."))
      tkfocus(tworkwin)
      melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      tkmark.set(tworkwin, "insert", paste(line+delta,"0",sep="."))

     })()
  }

  choosecol<-function(){
    (function(){
      ds<-function(pos=1,type,mode="numeric",struc=is.vector){
        if(pos>0){
         obj<-ls(pos=pos)
         obj<-obj[unlist(lapply(obj ,function(o,pos)
           exists(o,where=pos,mode="numeric") &&
           eval(parse(text=paste("struc(get(\"",o,"\",pos=",pos,"))",sep=""))),pos))]
        }else{
         obj<-ls(envir=revive.env)
         obj<-obj[unlist(lapply(obj ,function(o)
           exists(o,where=revive.env,mode="numeric") &&
           eval(parse(text=paste("struc(get(\"",o,"\",envir=revive.env))",sep=""))) ))]
        }
        if(0==length(obj)) return(NULL) else return(obj)
      }
      ds.of.R<-function(type="vector"){
        dat<-ls(pos=grep("datasets",search()))
        dat.type<-unlist(lapply(dat,function(x) {       
               num<-mode(x<-eval(parse(text=x)))
               num<-ifelse(is.array(x),"array",num)
               num<-ifelse(is.list(x),"list",num)
               num<-ifelse(is.matrix(x),"matrix",num)
               num<-ifelse(is.data.frame(x),"matrix",num)
               num<-ifelse(num=="numeric","vector",num)
               num }))
        return(dat[dat.type==type])
      }



      DS<-c(ds.of.R("matrix"),ds(struc=is.matrix))
      DS<-c(DS,"matrix(sample(100),50,2)","cbind(1:100,rnorm(100,mean=0,sd=1))")
      DS<-c(DS, ds(-1,struc=is.matrix), ds(1,struc=is.matrix) ,
                 ds(which(paste("package","relax"
,sep=":")==search())
 ,struc=is.matrix) )

      set.tclvalue<-function(name,value)tclvalue(name)<-as.character(value) # 121217
      listboxmenu<-function(choices,title="items",addexit=TRUE){
          if(addexit) choices<-c(choices,"EXIT")
          lbmTop<-tktoplevel();tkwm.geometry(lbmTop,"+0+15")
          tkwm.title(lbmTop,"Selection by Click and Return, Quit by Esc")
          ltit<-tklabel(lbmTop,text=title); lbframe<-tkframe(lbmTop)
          lb<-tklistbox(lbframe,height=8,width=60); sb<-tkscrollbar(lbframe)
          tkconfigure(lb,yscrollcommand=function(...) tkset(sb,...))
          tkconfigure(sb ,command=function(...) tkyview(lb,...))
          tkpack(ltit,lbframe)
          tkpack(sb,side="right",fill="y"); tkpack(lb)
          for(i in seq(choices)) tkinsert(lb,"end",choices[i])
          lbmdone      <- tclVar()
          tkbind(lbmTop,"<Return>",function()set.tclvalue(lbmdone,"1"))
          tkbind(lbmTop,"<Escape>",function()set.tclvalue(lbmdone,"0"))
          tkfocus(lb); tkwait.variable(lbmdone)
          choice<-tclvalue(tkget(lb,"active")); tkdestroy(lbmTop)
          if(tclvalue(lbmdone)=="0") return(0)
          ind <- match(choice, choices)
          choice <- if(addexit && ind==length(choices)) 0 else ind
          return(choice)
      }

      choice<-listboxmenu(DS,"Items:")
      if(0==length(choice)||is.na(choice)||0==choice) return()
      choice<-DS[choice]
      if(!exists("myscan")) myscan<-get("myscan",envir=revive.sys)



      choose.col<-function(obj,title="Which variable?"){
        if(is.character(obj)) try(obj<-eval(parse(text=obj)))
        obj <- as.matrix(obj)
        if(!is.matrix(obj)) return(0)
        if(is.null(choices<-dimnames(obj)[[2]])) choices<-1:ncol(obj)
        listboxmenu<-function(choices,title="items",addexit=TRUE){
            if(addexit) choices<-c(choices,"EXIT")
            lbmTop<-tktoplevel();tkwm.geometry(lbmTop,"+0+15")
            tkwm.title(lbmTop,"Selection by Click and Return, Quit by Esc")
            ltit<-tklabel(lbmTop,text=title); lbframe<-tkframe(lbmTop)
            lb<-tklistbox(lbframe,height=8,width=60); sb<-tkscrollbar(lbframe)
            tkconfigure(lb,yscrollcommand=function(...) tkset(sb,...))
            tkconfigure(sb ,command=function(...) tkyview(lb,...))
            tkpack(ltit,lbframe)
            tkpack(sb,side="right",fill="y"); tkpack(lb)
            for(i in seq(choices)) tkinsert(lb,"end",choices[i])
            lbmdone      <- tclVar()
            tkbind(lbmTop,"<Return>",function()set.tclvalue(lbmdone,"1"))
            tkbind(lbmTop,"<Escape>",function()set.tclvalue(lbmdone,"0"))
            tkfocus(lb); tkwait.variable(lbmdone)
            choice<-tclvalue(tkget(lb,"active")); tkdestroy(lbmTop)
            if(tclvalue(lbmdone)=="0") return(0)
            ind <- match(choice, choices)
            choice <- if(addexit && ind==length(choices)) 0 else ind
            return(choice)
        }

        choice<-listboxmenu(choices,title)
        choice
      }

   cat("choice:",choice)
      if(0==(col.no<-choose.col(choice))) return()
      cmd<-paste("x<-",choice,"[,",col.no,"]",sep="")
      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      melde <- get("melde", envir=revive.sys) # 121217
      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      last.text<-c(0,grep("^@",worktext)); last.text<-rev(last.text[last.text<=line])[1]
      last.code<-grep("^<<(.*)>>=",worktext)
      if(0<length(last.code)){
        last.code<-rev(last.code[last.code<=line])[1]
        if(is.na(last.code)||last.code<last.text)
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="") }
        else
                        {delta<-1; news<-paste(cmd,"\n",sep="")}
      }else{
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="")}
      }
      if(!exists("tworkwin"))
        tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

          anzrows<-length(news)
          try(tkinsert(tworkwin,paste(line+1,"0",sep="."),paste(news,collapse="\n")))
          tkmark.set(tworkwin, "insert", paste(line+anzrows,"0",sep="."))
          tksee(tworkwin,paste(line+anzrows,"0",sep="."))
      tkfocus(tworkwin)
      melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      tkmark.set(tworkwin, "insert", paste(line+delta,"0",sep="."))

     })()
  }

  get.dim.1.stats<-function(){
    (function(){
      methoden  <-c("median of x"                     ="median(x,na.rm=TRUE)"
                    ,"mean of x"                              ="mean(x,na.rm=TRUE)"
                    ,"standard deviation of x"         ="sd(x,rm.na=TRUE)"
                    ,"variance of x"                          ="var(x,na.rm=TRUE)"
                    ,"maximum of x"                   ="max(x,rm.na=TRUE)"
                    ,"minimum of x"                           ="min(x,na.rm=TRUE)"
                    ,"quantiles of x, prob: p"        ="quantile(x,p)"
                    ,"range of x"                             ="range(x,na.rm=TRUE)"
                    ,"inter quartile range of x"              ="IQR(x,na.rm=TRUE)"
                    ,"5 number summary of x"          ="fivenum(x,na.rm=TRUE)"
                    ,"summary statistics of x"        ="summary(x,na.rm=TRUE)"
                    ,"sorted data of x"               ="sort(x)"
                    ,"ranks of x"                                     ="rank(x)"
                    ,"sum of x"               ="sum(x)"
                    ,"cumulative sum of x"                    ="cumsum(x)"
                    ,"length of x"                    ="length(x)"
                    ,"frequency table of x"                   ="table(x)"
                    ,"relative frequencies of x"                      ="table(x)/length(x)"
                    )
      DS<-names(methoden)

      set.tclvalue<-function(name,value)tclvalue(name)<-as.character(value) # 121217
      listboxmenu<-function(choices,title="items",addexit=TRUE){
          if(addexit) choices<-c(choices,"EXIT")
          lbmTop<-tktoplevel();tkwm.geometry(lbmTop,"+0+15")
          tkwm.title(lbmTop,"Selection by Click and Return, Quit by Esc")
          ltit<-tklabel(lbmTop,text=title); lbframe<-tkframe(lbmTop)
          lb<-tklistbox(lbframe,height=8,width=60); sb<-tkscrollbar(lbframe)
          tkconfigure(lb,yscrollcommand=function(...) tkset(sb,...))
          tkconfigure(sb ,command=function(...) tkyview(lb,...))
          tkpack(ltit,lbframe)
          tkpack(sb,side="right",fill="y"); tkpack(lb)
          for(i in seq(choices)) tkinsert(lb,"end",choices[i])
          lbmdone      <- tclVar()
          tkbind(lbmTop,"<Return>",function()set.tclvalue(lbmdone,"1"))
          tkbind(lbmTop,"<Escape>",function()set.tclvalue(lbmdone,"0"))
          tkfocus(lb); tkwait.variable(lbmdone)
          choice<-tclvalue(tkget(lb,"active")); tkdestroy(lbmTop)
          if(tclvalue(lbmdone)=="0") return(0)
          ind <- match(choice, choices)
          choice <- if(addexit && ind==length(choices)) 0 else ind
          return(choice)
      }

      choice<-listboxmenu(DS,"Items:")
      if(0==length(choice)||is.na(choice)||0==choice) return()
      choice<-DS[choice]
      if(!exists("myscan")) myscan<-get("myscan",envir=revive.sys)



      cmd<-methoden[choice]
      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      melde <- get("melde", envir=revive.sys) # 121217
      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      last.text<-c(0,grep("^@",worktext)); last.text<-rev(last.text[last.text<=line])[1]
      last.code<-grep("^<<(.*)>>=",worktext)
      if(0<length(last.code)){
        last.code<-rev(last.code[last.code<=line])[1]
        if(is.na(last.code)||last.code<last.text)
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="") }
        else
                        {delta<-1; news<-paste(cmd,"\n",sep="")}
      }else{
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="")}
      }
      if(!exists("tworkwin"))
        tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

          anzrows<-length(news)
          try(tkinsert(tworkwin,paste(line+1,"0",sep="."),paste(news,collapse="\n")))
          tkmark.set(tworkwin, "insert", paste(line+anzrows,"0",sep="."))
          tksee(tworkwin,paste(line+anzrows,"0",sep="."))
      tkfocus(tworkwin)
      melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      tkmark.set(tworkwin, "insert", paste(line+delta,"0",sep="."))

     })()
  }
  get.dim.1.plots<-function(){
    (function(){
      methoden  <-c("1 dim plot of x"="plot(x)"
       ,"boxplot of x"              ="boxplot(x)"
       ,"jitterplot of x"              ="plot(jitter(x))"
       ,"histogram of x"              ="hist(x,prob=TRUE)"
       ,"histogram of x, breaks at br"           ="hist(x,breaks=br,prob=TRUE)"
       ,"plot of density trace of x, width w"="plot(density(x,width=w),type=\"l\")"
       ,"barplot of x"              ="barplot(x)"
       ,"barplot of relative frequencies of x"              ="plot(table(x)/length(x))"
       ,"barplot of frequencies h.i at x.i"      ="plot(x.i, h.i, type=\"h\")"
       ,"empirical distribution function  of x"              ="plot(ecdf(x))"
       ,"stem and leaf display of x"="if(exists(\"stem.leaf\")) stem.leaf(x) else stem(x)"
      )
      DS<-names(methoden)

      set.tclvalue<-function(name,value)tclvalue(name)<-as.character(value) # 121217
      listboxmenu<-function(choices,title="items",addexit=TRUE){
          if(addexit) choices<-c(choices,"EXIT")
          lbmTop<-tktoplevel();tkwm.geometry(lbmTop,"+0+15")
          tkwm.title(lbmTop,"Selection by Click and Return, Quit by Esc")
          ltit<-tklabel(lbmTop,text=title); lbframe<-tkframe(lbmTop)
          lb<-tklistbox(lbframe,height=8,width=60); sb<-tkscrollbar(lbframe)
          tkconfigure(lb,yscrollcommand=function(...) tkset(sb,...))
          tkconfigure(sb ,command=function(...) tkyview(lb,...))
          tkpack(ltit,lbframe)
          tkpack(sb,side="right",fill="y"); tkpack(lb)
          for(i in seq(choices)) tkinsert(lb,"end",choices[i])
          lbmdone      <- tclVar()
          tkbind(lbmTop,"<Return>",function()set.tclvalue(lbmdone,"1"))
          tkbind(lbmTop,"<Escape>",function()set.tclvalue(lbmdone,"0"))
          tkfocus(lb); tkwait.variable(lbmdone)
          choice<-tclvalue(tkget(lb,"active")); tkdestroy(lbmTop)
          if(tclvalue(lbmdone)=="0") return(0)
          ind <- match(choice, choices)
          choice <- if(addexit && ind==length(choices)) 0 else ind
          return(choice)
      }

      choice<-listboxmenu(DS,"Items:")
      if(0==length(choice)||is.na(choice)||0==choice) return()
      choice<-DS[choice]
      if(!exists("myscan")) myscan<-get("myscan",envir=revive.sys)



      cmd<-methoden[choice]
      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      melde <- get("melde", envir=revive.sys) # 121217
      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      last.text<-c(0,grep("^@",worktext)); last.text<-rev(last.text[last.text<=line])[1]
      last.code<-grep("^<<(.*)>>=",worktext)
      if(0<length(last.code)){
        last.code<-rev(last.code[last.code<=line])[1]
        if(is.na(last.code)||last.code<last.text)
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="") }
        else
                        {delta<-1; news<-paste(cmd,"\n",sep="")}
      }else{
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="")}
      }
      if(!exists("tworkwin"))
        tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

          anzrows<-length(news)
          try(tkinsert(tworkwin,paste(line+1,"0",sep="."),paste(news,collapse="\n")))
          tkmark.set(tworkwin, "insert", paste(line+anzrows,"0",sep="."))
          tksee(tworkwin,paste(line+anzrows,"0",sep="."))
      tkfocus(tworkwin)
      melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      tkmark.set(tworkwin, "insert", paste(line+delta,"0",sep="."))

     })()
  }

  get.dim.2.methods<-function(){
    (function(){
      methoden  <-c("scatter plot of xy[,1:2]"="plot(xy[,1:2])",
                    "correlation of xy" ="cor(xy)",
                    "regression line of xy"  ="plot(xy[,1:2])\nabline(lsfit(xy[,1],xy[,2]))",
                    "mean values of cols of xy" ="apply(xy,2,mean)"
                  )
      DS<-names(methoden)

      set.tclvalue<-function(name,value)tclvalue(name)<-as.character(value) # 121217
      listboxmenu<-function(choices,title="items",addexit=TRUE){
          if(addexit) choices<-c(choices,"EXIT")
          lbmTop<-tktoplevel();tkwm.geometry(lbmTop,"+0+15")
          tkwm.title(lbmTop,"Selection by Click and Return, Quit by Esc")
          ltit<-tklabel(lbmTop,text=title); lbframe<-tkframe(lbmTop)
          lb<-tklistbox(lbframe,height=8,width=60); sb<-tkscrollbar(lbframe)
          tkconfigure(lb,yscrollcommand=function(...) tkset(sb,...))
          tkconfigure(sb ,command=function(...) tkyview(lb,...))
          tkpack(ltit,lbframe)
          tkpack(sb,side="right",fill="y"); tkpack(lb)
          for(i in seq(choices)) tkinsert(lb,"end",choices[i])
          lbmdone      <- tclVar()
          tkbind(lbmTop,"<Return>",function()set.tclvalue(lbmdone,"1"))
          tkbind(lbmTop,"<Escape>",function()set.tclvalue(lbmdone,"0"))
          tkfocus(lb); tkwait.variable(lbmdone)
          choice<-tclvalue(tkget(lb,"active")); tkdestroy(lbmTop)
          if(tclvalue(lbmdone)=="0") return(0)
          ind <- match(choice, choices)
          choice <- if(addexit && ind==length(choices)) 0 else ind
          return(choice)
      }

      choice<-listboxmenu(DS,"Items:")
      if(0==length(choice)||is.na(choice)||0==choice) return()
      choice<-DS[choice]
      if(!exists("myscan")) myscan<-get("myscan",envir=revive.sys)



      cmd<-methoden[choice]
      if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
      tworkwin<-get("tworkwin",envir=revive.sys)
      worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
      if(nchar(worktext)<10000){
        worktext<-strsplit(worktext,"\n")[[1]]
      }else{
        base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
        worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
      }

      melde <- get("melde", envir=revive.sys) # 121217
      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      last.text<-c(0,grep("^@",worktext)); last.text<-rev(last.text[last.text<=line])[1]
      last.code<-grep("^<<(.*)>>=",worktext)
      if(0<length(last.code)){
        last.code<-rev(last.code[last.code<=line])[1]
        if(is.na(last.code)||last.code<last.text)
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="") }
        else
                        {delta<-1; news<-paste(cmd,"\n",sep="")}
      }else{
                  {delta<-4; news<-paste("\n@\n<<*>>=\n",cmd,"\n@\n",sep="")}
      }
      if(!exists("tworkwin"))
        tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

          anzrows<-length(news)
          try(tkinsert(tworkwin,paste(line+1,"0",sep="."),paste(news,collapse="\n")))
          tkmark.set(tworkwin, "insert", paste(line+anzrows,"0",sep="."))
          tksee(tworkwin,paste(line+anzrows,"0",sep="."))
      tkfocus(tworkwin)
      melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      tkmark.set(tworkwin, "insert", paste(line+delta,"0",sep="."))

    })()
  }

  data.fns.menu<-function(){
    myhead.menu(rm.menu=TRUE, menu.no=0)
    # myhead.menu(item="load data via internet",code=loadwwwdata, menu.no=0)   # ok
    myhead.menu(item="1-dim data sets",code=get.dim.1.data,title="Data",menu.no=0)    # ok
    myhead.menu(item="1-dim random numbers",code=generate.1.data,menu.no=0) # ok
    myhead.menu(item="2-dim data sets",code=get.dim.2.data,menu.no=0)    # ok
    myhead.menu(item="delete NAs from vector x",code=removenafromvecx,menu.no=0)
    myhead.menu(item="delete NAs from matrix xy",code=removenafrommatxy,menu.no=0)
    myhead.menu(item="save col of matrix xy as x",code=choosecol,menu.no=0)
    myhead.menu(item="1-dim statistics",code=get.dim.1.stats,title="Methods",menu.no=1)   # ok
    myhead.menu(item="1-dim plots",code=get.dim.1.plots,title="Methods",menu.no=1)   # ok
    myhead.menu(item="2-dim methods",code=get.dim.2.methods,menu.no=1)   # ok
  }

  LAST.WARNING<-"no warnings"
  assign("LAST.WARNING",LAST.WARNING,envir=revive.sys)
  Rversion<-as.numeric(R.version$major)*100+as.numeric(R.version$minor)

  REVFILE            <- "REVFILE"    # eingelesener RevFile
  RCHFILE            <- "RCHFILE"    # eingelesener Chunk-File
  fr.paper.sys       <- "forget"     #
  relax.version.sys<- "relax 1.3.15 - 140310"

  tvexit       <- tclVar("0")
  tvchoice     <- tclVar("0")
  tvndone      <- tclVar("0")
  tvscandone   <- tclVar("0")
  tvinfo       <- tclVar("")      # Variable des Kopf-Entry-Widget
  tvreadline   <- tclVar("0")
  tvmess       <- tclVar("relax")


  revpath.sys  <- getwd()         # Pfad zum Revwebpaperverzeichnis
  secno.sys    <- "0"             # aktuelle SecNo. in der Revdatei
  string.sys   <- ""              # Default-Suchstring
  sweave.args.sys<-""  # Default-Setzung weiterer Sweave-Argumente

  sizes<-c("8-80","10-100","12-120","14-140","18-180","24-240","*-2000")
  tfont.sys    <- "-Adobe-helvetica-Medium-R-Normal--14-140-*" # Text-Schrift
  tfont.sys<-sub("al--.+-.+-\\*",paste("al--",sizes[initial.font.size],"-*",sep=""),tfont.sys)
  outfont.sys  <- "-Adobe-courier-Medium-R-Normal--14-140-*"   # Output-Schrift
  outfont.sys<-sub("al--.+-.+-\\*",paste("al--",sizes[initial.font.size],"-*",sep=""),outfont.sys)
  workname.sys<-"out.rev"

  if((.Platform$OS.type=="windows") | substring(version$os,1,5)=="linux"){
    tkevent.add("<<Down>>",        "<Alt_L><d>")
    tkevent.add("<<Up>>",          "<Alt_L><u>")
    tkevent.add("<<PlanRCode>>",   "<Alt_L><p>")
    tkevent.add("<<EvalRCode>>",   "<Alt_L><e>")
    tkevent.add("<<WarnEval>>",    "<Alt_L><w>")
    tkevent.add("<<Insert>>",      "<Alt_L><i>")
    tkevent.add("<<SavePlot>>",    "<Alt_L><s>")
    tkevent.add("<<RemoveOut>>",   "<Alt_L><r>")
    tkevent.add("<<TrashROutput>>","<Alt_L><t>")
    tkevent.add("<<Help.R>>",      "<Alt_L><h>")
    tkevent.add("<<FindText>>",    "<Alt_L><f>") ## <F3> without virtual event
    tkevent.add("<<FindReportText>>",   "<Control_L><f>")
    tkevent.add("<<GoToLine>>",         "<Control_L><g>")
    tkevent.add("<<SaveReport>>",       "<Control_L><s>")
    tkevent.add("<<ProcessReport>>",    "<Control_L><p>")
    tkevent.add("<<Next>>",        "<Alt_L><n>")
    tkevent.add("<<Back>>",        "<Alt_L><b>")
    ## important for special characters
    tkevent.add("<<RKgeschw>>",    "<Alt_L><0>")
    tkevent.add("<<LKgeschw>>",    "<Alt_L><7>")
    tkevent.add("<<RKeckig>>",     "<Alt_L><9>")
    tkevent.add("<<LKeckig>>",     "<Alt_L><8>")
    tkevent.add("<<Tilde>>",       "<Alt_L><plus>")
    tkevent.add("<<Klammera>>",    "<Alt_L><q>")
    tkevent.add("<<Pipe>>",        "<Alt_L><less>")
    tkevent.add("<<aeumlaut>>",      "<adiaeresis><KeyRelease>")
    tkevent.add("<<oeumlaut>>",      "<odiaeresis><KeyRelease>")
    tkevent.add("<<ueumlaut>>",      "<udiaeresis><KeyRelease>")
    tkevent.add("<<Aeumlaut>>",      "<Adiaeresis><KeyRelease>")
    tkevent.add("<<Oeumlaut>>",      "<Odiaeresis><KeyRelease>")
    tkevent.add("<<Ueumlaut>>",      "<Udiaeresis><KeyRelease>")
    tkevent.add("<<szumlaut>>",      "<ssharp><KeyRelease>")
    tkevent.add("<<Backsl>>",      "<Alt_L><ssharp>")
    #tkevent.add("<<Klammeraffe>>", "<ISO_Level3_Shift><KeyPress><KeyRelease><KeyRelease>") unused
  } 
  if(substring(version$os,1,6)=="darwin" ){
    tkevent.add("<<Down>>",        "<Meta_L><d>") # don't work: deletes characters after Cursor 
    tkevent.add("<<Up>>",          "<Meta_L><u>")
    tkevent.add("<<PlanRCode>>",   "<Meta_L><p>")
    tkevent.add("<<EvalRCode>>",   "<Meta_L><e>")
    tkevent.add("<<Insert>>",      "<Meta_L><i>")
    tkevent.add("<<SavePlot>>",    "<Meta_L><s>")
    tkevent.add("<<RemoveOut>>",   "<Meta_L><r>")
    tkevent.add("<<TrashROutput>>","<Meta_L><t>")
    tkevent.add("<<Help.R>>",      "<Meta_L><h>") # don't work: hide relax window
    tkevent.add("<<FindText>>",    "<Meta_L><f>")
    tkevent.add("<<FindReportText>>",  "<Control_L><f>")
    tkevent.add("<<GoToLine>>",        "<Control_L><g>")
    tkevent.add("<<SaveReport>>",      "<Control_L><s>")
    tkevent.add("<<ProcessReport>>",   "<Control_L><p>")
    tkevent.add("<<Next>>",        "<Meta_L><n>") # Meta_L n opens new window
    tkevent.add("<<Back>>",        "<Meta_L><b>") 
    # tkevent.add("<<WarnEval>>",    "<Meta_L><w>") # willl destroy relax window
    ## important for special characters, aeoe...will work 
    tkevent.add("<<RKgeschw>>",  "<Meta_L><0>")
    tkevent.add("<<LKgeschw>>",  "<Meta_L><7>")
    tkevent.add("<<RKeckig>>",   "<Meta_L><9>")
    tkevent.add("<<LKeckig>>",   "<Meta_L><8>")
    tkevent.add("<<Klammera>>",  "<Meta_L><q>")
    tkevent.add("<<Pipe>>",      "<Meta_L><less>")
    #tkevent.add("<<aeumlaut>>",      "<adiaeresis><KeyRelease>")
    #tkevent.add("<<oeumlaut>>",      "<odiaeresis><KeyRelease>")
    #tkevent.add("<<ueumlaut>>",      "<udiaeresis><KeyRelease>")
    #tkevent.add("<<Aeumlaut>>",      "<Adiaeresis><KeyRelease>")
    #tkevent.add("<<Oeumlaut>>",      "<Odiaeresis><KeyRelease>")
    #tkevent.add("<<Ueumlaut>>",      "<Udiaeresis><KeyRelease>")
    #tkevent.add("<<szumlaut>>",      "<ssharp><KeyRelease>")
    #tkevent.add("<<Backsl>>",        "<Meta_L><ssharp>")
    #tkevent.add("<<Tilde>>",    "<Meta_L><plus>")
    #tkevent.add("<<Klammeraffe>>", "<ISO_Level3_Shift><KeyPress><KeyRelease><KeyRelease>") unused
  } 
         ##definiere Testknopf-Ereignis>>
  implement.but<-
  function(but,frame,mess=" ",side="right",relief="raised",short.cut=TRUE,bwf=1,job=""){ ## 071115
    if(is.character(frame)) frame<-eval(parse(text=frame))
    b<-tkbutton(frame, text=but, relief=relief, pady="-3"
,
                width=floor(10
*bwf), font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"
)
    if(short.cut) tkconfigure(b,underline=0)
    tkbind(b,"<Enter>",function()set.tclvalue("tvmess",mess))
    tkbind(b,"<Leave>",function(){ set.tclvalue("tvmess","relax") })
    tkbind(b,"<Button-1><ButtonRelease>",job) ## 071115
    tkpack(b, side=side)
    assign(but, b, envir=revive.sys)
  }
  melde("Implement.but defined",3)

  TopW<-tktoplevel(); tkwm.geometry(TopW,"+0+15")
  tkwm.title(TopW,paste(
               if(but.Wizardry=="simple") "redit -- simple Report EDITor for statistical analysis:" else
                                          "relax -- Report Editor for Literate Analysis and lateX:",
                        "relax 1.3.15 - 140310"))
  tkwm.protocol(TopW,"WM_DELETE_WINDOW",function(){
                     if(!exists("tworkwin"))
                       tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

                     worktext<-TcltoWin.write(tclvalue(tkget(tworkwin,"0.0","end")))
                     get("cat","package:base")(worktext,file="report-UnDo-bak.rev")

                     cat("... not a nice way to exit relax !!!\n")
                     cat("backup file: report-UnDo-bak.rev\n")
                     if( (.Platform$OS.type=="windows") ){
                        cat("--> to proceed the R session try: Quit R, then don't quit!\n")
                     }
                     if(substring(version$os,1,5)=="linux") cat("-- try: <Ctrl-C> to proceed the R session!\n")
                     set.tclvalue("tvexit","fertig"); tkdestroy(TopW)
             }
  )


  fhead<- tkframe(TopW, relief="raised", bd="1")
  #### finfo<- tkframe(TopW)  ###100208 shift to foutcmds
  # finout   <- tkframe(TopW,height="700",width="700") 
  finout<- tkframe(TopW,height=relaxwindow.height.sys,width=relaxwindow.width.sys) 
  tkpack(fhead, side="top",fill="x")
  #### tkpack(finfo,   side="top") ###100208 shift to foutcmds
  tkpack(finout, side="top",fill="both",expand="1")

  fworkwin<-tkframe(finout)
  fout        <-tkframe(finout)
  fhandle   <-tkframe(finout,bd="2",relief="raised", bg="#BEE3D2",   # #87ff9D",
                                     cursor="sb_v_double_arrow")
  tkplace(fworkwin, relwidth="1",rely="0",height="-1",anchor="nw")
  tkplace(fout, relwidth="1",rely="1",height="-1",anchor="sw")
  tkplace(fhandle, relx="1",width="1500",height="6",anchor="e")

  proc<-paste(
    paste("bind ",TopW$ID," <Configure> {" ),
        paste("set H [winfo height ",finout$ID, "].0" ),
        paste("set Y0  [winfo rooty ",finout$ID, "]" ),
    "}",
    paste("bind",fhandle$ID," <B1-Motion> {"),
      "set fract [expr (%Y -$Y0)/$H]",
      "if { $fract < 0.1 } {",
      "  set fract 0.1",
      "}",
      "if { $fract > 0.95 } {",
      "  set fract 0.95",
      "}",
      paste("place ",fworkwin$ID," -relheight $fract"),
      paste("place ",fhandle$ID," -rely $fract"),
      paste("place ",fout," -relheight [expr 1.0 - $fract]"),
    "}",
    sep="\n")
  .Tcl(proc)

  fract<- 0.5
  tkplace(fworkwin, relheight=fract)
  tkplace(fhandle, rely=fract)
  tkplace(fout, relheight=1-fract)

  fworkcmds<-tkframe(fout, relief="raised", bd="0")
  foutcmds<-tkframe(fout, relief="raised", bd="0") # ,background="#37D70F")
  foutwin   <- tkframe(fout)

  tkpack(fworkcmds,foutcmds,fill="x")
  tkpack(foutwin,fill="both",expand="yes")

  finfo<- tkframe(foutcmds); 
  tkpack(finfo, side="left",fill="both",expand="yes")  ###100208
  ### tkpack(tklabel(foutcmds, text=" Result(s):",pady=#<negativer y-Platz>#,
  ###     font=#<Font für Knöpfe>#), side="left")  #ipady=7, 

  mbFile<-tkmenubutton(fhead, text="File", font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"
,
                       relief="flat", width=10
) # alternativ groove
  mbEdit<-tkmenubutton(fhead, text="Edit", font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"
,
                       relief="flat", width=10
)
  mbOptions<-tkmenubutton(fhead, text="Options", font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"
,
                          relief="flat", width=10
)
  if(but.Wizardry == "all"){
    mbRevweb<-tkmenubutton(fhead, text="Wizardry", font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"
,
                          relief="flat", width=10
)
  }

  f<-function(mess) function()set.tclvalue("tvmess",mess)
  tkbind(mbFile, "<Enter>", f("file operations and exit"))
  tkbind(mbEdit, "<Enter>", f("searching and other operations"))
  tkbind(mbOptions,  "<Enter>", f("change settings"))
  if(but.Wizardry=="all")  tkbind(mbRevweb,  "<Enter>", f("process document and specials"))
  tkbind(mbFile, "<Leave>", function(){ set.tclvalue("tvmess","relax") })
  tkbind(mbEdit, "<Leave>", function(){ set.tclvalue("tvmess","relax") })
  tkbind(mbOptions, "<Leave>", function(){ set.tclvalue("tvmess","relax") })
  if(but.Wizardry=="all")tkbind(mbRevweb, "<Leave>", function(){ set.tclvalue("tvmess","relax") })
  tkbind(finfo, "<Enter>", f("entry / message field"))
  tkbind(finfo, "<Leave>", function(){ set.tclvalue("tvmess","relax") })

  tkpack(mbFile,mbEdit,mbOptions,side="left")
  if(but.Wizardry=="all")tkpack(mbRevweb,side="left")
  implement.but("Help.R", "fhead", "show online documentation of R object",job=fHelp.R)
  implement.but("Examples", "fhead", "show arguments and examples of objects",
                job=fExamples)

  mbFile.menu<-tkmenu(mbFile,font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"
)
  tkconfigure(mbFile, menu=mbFile.menu)

  mbEdit.menu<-tkmenu(mbEdit,font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"
)
  tkconfigure(mbEdit, menu=mbEdit.menu)
  tkadd(mbEdit.menu, "command", command=ShowAboutRelax,
        label="ShowAboutRelax:   what's relax?")
  if(but.Wizardry!="simple") {
    tkadd(mbEdit.menu, "command", command=ShowShortCuts,
          label="ShowShortCuts:   show short cuts for text field")
  }
  tkadd(mbEdit.menu, "command", command=function(){
                                          if(FALSE){
                                            res<-tkmessageBox(message=if(language=="german") "Soll das interaktive Icon aktiert werden?"
                                                                      else "Do you want to activate interactive icon?",
                                                              title="RELAX-ICON",icon="warning",type="yesnocancel",default="no")
                                            if("externalptr"==mode(res))  res<-tclvalue(res)
                                            if(res=="cancel") return()
                                            if(res=="no") return()
                                          }
                                        chair<-function(xcenter=1.2,xspread=1,ycenter=1.2,yspread=1,lwd=4){
                                          # cat("chair: ", xcenter,xspread,ycenter,yspread)
                                          x<-seq(0,1.6,length=(n<-40))
                                          xshift<-.5;ystretch<-1.2; height<-.4/yspread
                                          y<-0.5*cosh(ystretch*x-xshift);
                                          y<-scale(y-y[1]-.06,-ycenter,1/yspread)
                                          x<-scale(x/1.6,-xcenter,1/xspread)
                                          segments(x[1],y[1],x[1]*.1+x[n]*.9,y[1]-height,lwd=lwd,xpd=NA,col="blue")
                                          segments(x[1],y[1]-height,x[n],y[n],lwd=lwd,xpd=NA,col="blue")
                                          lines(x,y,lwd=lwd,xpd=NA,col="blue")
                                        }
                                        tree<-function(xcenter=.51,xspread=1.5,ycenter=0,yspread=2,lwd=3,n=5,n.l=10){
                                          # cat(xcenter,xspread,ycenter,yspread)
                                          yspread<-yspread/2.5; ycenter<-ycenter+.15
                                          leaf<-function(von,start,delta,len,n.l=10,krum.par){
                                           xy.mat<-matrix(von,2,n.l)
                                           h<-start*2*pi/360
                                           xy.delta<-c(cos(h),sin(h))/n.l
                                           delta<-delta*2*pi/360*len/n.l
                                           for(i in 2:n.l){
                                            w<-delta
                                            xy.delta<-matrix(c(cos(w),-sin(w),sin(w),cos(w)),2,2)%*%xy.delta
                                            xy.mat[,i]<-xy.mat[,i-1]+xy.delta
                                           }
                                           delta.y<-qbeta((0:n.l)/n.l,2,3)
                                           delta.y<-.04*cumsum(sin((0:(n.l-1))/(n.l)*2*pi))
                                           xy<-t(xy.mat) # +von-xy.mat[,1])
                                           xy[,1]<-scale(xy[,1]*xspread/2,-xcenter,1)#2/xspread)
                                           xy[,2]<-(xy[,2]*yspread+ycenter)
                                           delta.y<-(delta.y*yspread)
                                           #lines(xy[,1],xy[,2],lwd=3) #lines(xy[,1],xy[,2]+delta.y,lwd=3)
                                           polygon(c(xy[,1],rev(xy[,1])),c(xy[,2],rev(xy[,2]+delta.y)),col="green",border=NA,xpd=NA)
                                          }
                                          x<-seq(-50,90,length=floor(n/2))
                                          for(i in x){
                                            leaf(von=c(0,0),start=i,delta=5,len=10,n.l=n.l)
                                            x<-seq(-50,90,length=ceiling(n/2))
                                          }
                                          for(i in x) leaf(von=c(0,0),start=190-i,delta=-5,len=10,n.l=n.l)
                                          segments(xcenter-0.02*xspread,ycenter,
                                                   xcenter-0.05*xspread,ycenter-1.5*yspread,lwd=1)
                                          segments(xcenter+0.02*xspread,ycenter,
                                                   xcenter+0.05*xspread,ycenter-1.5*yspread,lwd=1)
                                        } # end of tree
                                        redo<-function(...){
                                          ytrsc<-xtrsc<-slider(no=1); xchsc<-slider(no=2)
                                          ychsc<-slider(no=3);hori<-slider(no=4);n.leafs<-slider(no=5)
                                          leaf.sty=slider(no=6)
                                          h<-slider(obj.name="tree.center"); xtrcenter<-h[1]; ytrcenter<-h[2]
                                          h<-slider(obj.name="chair.center"); xchcenter<-h[1]; ychcenter<-h[2]
                                        #  a<-par(family="mono") # 121217
                                          plot(-2:2,-2:2,type="n",axes=FALSE,xlab="",ylab=""); # text(.6,1.8,"relax",cex=4)
                                          #polygon(2*c(-2,-2,2,2),c(hori,-3,-3,hori),col="#fcf1a2",border=NA,xpd=NA) # strand
                                          #polygon(2*c(-2,-2,2,2),c(hori,3,3,hori),col="#c8e6f2",border=NA,xpd=NA)   # air
                                          n <- 1000; xa<-rep(-4,n); xb <- xa + 8
                                          y <- seq(hori,-3.5,length=n); segments(xa,y,xb,y,col=rainbow(n,start=0.11,end=0.18),xpd=NA)
                                          y <- seq(hori, 3.5,length=n); segments(xa,y,xb,y,col=rainbow(n,start=0.53,end=0.70),xpd=NA)
                                          b<-as.numeric(unlist(strsplit(sub(".*([0-9][0-9]:[0-9][0-9]:).*","\\1",date()),":")))
                                          b[2]<-b[2]/60; bb<-b; b[1] <- b[1]+6; b <- (sum(b)%%12)/12
                                          cen <- c(-2*cos(pi*b),sin(pi*b)*(2-hori)+hori+.5)
                                          points(cen[1],cen[2],cex=15,col=heat.colors(18)[6+floor(abs(20*min(b,1-b)))],pch=16,xpd=NA)
                                          h <- 2*pi*(sum(bb)/12); h <- cen + c(sin(h),cos(h))*0.3
                                          arrows(cen[1],cen[2],h[1],h[2],.08,col="blue",lwd=4,xpd=NA)
                                          h <- 2*pi*(bb[2]); h <- cen + c(sin(h),cos(h))*0.4
                                          arrows(cen[1],cen[2],h[1],h[2],.08,col="blue",lwd=4,xpd=NA)
                                          tree(xcenter=xtrcenter,xspread=xtrsc,n=n.leafs,n.l=leaf.sty,
                                               ycenter=ytrcenter,yspread=ytrsc,lwd=3)
                                          chair(xcenter=xchcenter,xspread=xchsc,
                                                ycenter=ychcenter,yspread=ychsc,lwd=3)
                                          # text(1.6,-2,"relax",xpd=NA) # 121217 letters: RELAX:
                                          f <- .1; lwd <- f*30; fr <- 1.5; dx <- -2.5; dy <- 2.5; xdel <- .15
                                          segments(dx-0.5*xdel+fr*f*c(0,0,0,0,1),dy+fr*f*c(0,1,1,2,1),
                                                   dx-0.5*xdel+fr*f*c(0,1,1,1,1),dy+fr*f*c(2,0,1,2,2),lwd=lwd,xpd=NA)#R
                                          dx <- dx + xdel
                                          segments(dx+f*c(0,0,0,0),dy+f*c(0,0,1,2),dx+f*c(0,1,1,1),dy+f*c(2,0,1,2),lwd=lwd,xpd=NA)#E
                                          dx <- dx + xdel
                                          segments(dx+f*c(0,0),dy+f*c(2,0),dx+f*c(0,1),dy+f*c(0,0), lwd=lwd,xpd=NA)#L
                                          dx <- dx + xdel
                                          segments(dx+f*c(0,0.5,.25),dy+f*c(0,2,.8),dx+f*c(0.5,1,0.75),dy+f*c(2,0,.8),
                                                   lwd=lwd,xpd=NA)#A
                                          dx <- dx + xdel
                                          segments(dx+f*c(0,0), dy+f*c(0,2),dx+f*c(1,1),dy+f*c(2,0), lwd=lwd,xpd=NA)#X
                                        #   par(family=a)
                                        }
                                          slider(obj.name="tree.center",obj.value=c(-1.1,.8))
                                          slider(obj.name="chair.center",obj.value=c(.3,-1.2))
                                          set.tree<-function(...){
                                            # cat("click on desired position in graphics device!!!\n")
                                            # xy<-unlist(locator(n=1)); xy<-c(xy$x,xy$y)
                                            xy<-c(runif(1,-1,1),runif(1,-1.5,1.0))
                                            slider(obj.name="tree.center",obj.value=xy)
                                            redo()
                                          }
                                          set.chair<-function(...){ #   xy<-locator(n=1); xy<-c(xy$x,xy$y)
                                            xy<-c(runif(1,-1,1),runif(1,-1.5,0))
                                            slider(obj.name="chair.center",obj.value=xy)
                                            redo()
                                          }
                                          if("tkrplot" %in% list.files(.libPaths())) sl <- gslider else sl <- slider
                                          sl(redo,
                                            c("tree scale factor", "x scale chair",
                                              "design of chair", "size of beach",
                                              "number of leafs", "style of leafs"),
                                            c(.1,.1,.1,-2,7,5),c(5,2.4,2.4,2,17,30),
                                            c(.1,.1,.1,.1,1,1),c(3.5,1.2,1.2,-0.7,9,12),
                                            list(set.tree,set.chair,redo),
                                            but.names=c("tree location","new chair location","update clock")
                                          )
                                          if("tkrplot" %in% list.files(.libPaths())) "relax" else redo()
                                          "relax"
                                        }
,
        label="ShowInteractiveIcon:   show an interactive icon")
  if(but.Wizardry!="simple") {
    tkadd(mbEdit.menu, "separator")
    tkadd(mbEdit.menu, "command", command=ReloadPlots,
          label="ReloadPlots:   reload jpeg-plot of text field")
    tkadd(mbEdit.menu, "command", command=RefreshChunkNumbers,
          label="RefreshChunkNumbers:   refresh the chunks numbers after chunk names")
    tkadd(mbEdit.menu, "command", command=ReloadReportWidget,
          label="ReloadReportWidget:   reconstruct REPORT WINDOW in case of strange appearances")
    tkadd(mbEdit.menu, "separator")
  }
  tkadd(mbEdit.menu, "command", command=GoToLine,
          label="GoToLine:   go to line ... ")
  if(but.Wizardry!="simple") {
    tkadd(mbEdit.menu, "command", command=Replace, 
          label="SearchReplace:  search and replace text strings") #100917
  }
  tkadd(mbEdit.menu, "command", command=FindRFns,
        label="FindRFns:   search R function by keyword")
  tkadd(mbEdit.menu, "command", command=FindReportText,
        label="FindReportText:   search text string in text field  (Crtl+F)")
  if(but.Wizardry!="simple") {
    tkadd(mbEdit.menu, "command", command=FindReportChunk,
          label="FindReportChunk:   search code chunk in text field")
    tkadd(mbEdit.menu, "command", command=FindLaTeXSection,
          label="FindLaTeXSection:   search for \\section, \\subsection and \\subsubsection")
    tkadd(mbEdit.menu, "command", command=ListUsedPlots,
          label="ListUsedPlots:   list graphics files of directory / referenced in the report")
    tkadd(mbEdit.menu, "separator")
    tkadd(mbEdit.menu, "command", command=InsertLaTeXEnv,
          label="InsertLaTeXEnv:   insert LaTeX environment")
  }

  mbOptions.menu<-tkmenu(mbOptions, font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"
)
  tkconfigure(mbOptions, menu=mbOptions.menu)
  tkadd(mbOptions.menu,"command", command=SetOutputLength,
        label="SetOutputLength: define maximal lines of output")
  if(but.Wizardry!="simple") 
    tkadd(mbOptions.menu, "separator")
  tkadd(mbOptions.menu,"command", command=SetFontType,
        label="SetFontType:   define font type")
  tkadd(mbOptions.menu,"command", command=SetFontSize,
        label="SetFontSize:   define font size")
  if(but.Wizardry!="simple") {
    tkadd(mbOptions.menu, "separator")
    tkadd(mbOptions.menu,"command", command=SetRelaxWinSize,
          label="SetRelaxWinSize:  redefine width x height of relax window")
    tkadd(mbOptions.menu, "separator")
    tkadd(mbOptions.menu,"command", command=SetPlotHeight,
          label="SetPlotHeight:  define height of plot (-> latex)")
    tkadd(mbOptions.menu,"command", command=SetPSDesignWidth,
          label="SetPSWidth:   define width of ps-graphics")
    tkadd(mbOptions.menu,"command", command=SetPSDesignHeight,
          label="SetPSHeight:  define height of ps-graphics")
    tkadd(mbOptions.menu,"command", command=SetPSRotation,
          label="SetPSRotation:  define non normal PS rotation")
    tkadd(mbOptions.menu,"command", command=SetJPGSize,
          label="SetJPGSize:  define size of jpeg-graphics")
    tkadd(mbOptions.menu,"command", command=SetPPMResolution,
          label="SetPPMResolution:  define resolution of ppm-graphics")
    tkadd(mbOptions.menu, "separator")
    tkadd(mbOptions.menu,"command", command=ConvertEncodingToLocal,
          label="ConvertEncodingToLocal:  convert encoding to local one")
    tkadd(mbOptions.menu,"command", command=ConvertEncodingFromLocal,
          label="ConvertEncodingFromLocal:  convert encoding to a new one")
    tkadd(mbOptions.menu,"command", command=ConvertGermanUmlauteToLocalEncoding,
          label="ConvertGermanUmlauteToLocalEncoding:  try to convert German-Umlaute to locale encoding")
    tkadd(mbOptions.menu, "separator")
  }
  tkadd(mbOptions.menu,"command", command=ConfigRelax,
        label="Configure Relax: view or change parameters of relax")

  if(but.Wizardry=="all"){
   mbRevweb.menu<-tkmenu(mbRevweb,font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"
)
   tkconfigure(mbRevweb, menu=mbRevweb.menu)
   tkadd(mbRevweb.menu,"command", command=ProcessReport,
        label="ProcessReport:  SaveReport, WeaveReport, and LatexReport")
   tkadd(mbRevweb.menu,"command", command=ViewReport,
        label="ViewReport:   show formated report")
   tkadd(mbRevweb.menu,"command", command=ProcessChunk,
        label="ProcessChunk:  SaveChunk as 'local-chunk', WeaveChunk, and LatexChunk")  ## experimental
   tkadd(mbRevweb.menu,"command", command=ProcessBeginEndEnv,
        label="ProcessBeginEndEnv:  Save Latex Environment as 'local-chunk', Weave, and Latex")  ## experimental
   tkadd(mbRevweb.menu, "separator")
   tkadd(mbRevweb.menu,"command", command=LaTeX.head,
        label="LaTeX.head:   include simple LaTeX head")
   tkadd(mbRevweb.menu,"command", command=LatexReport,
        label="LatexReport:   format Latex file")
   tkadd(mbRevweb.menu,"command", command=ShowLogFile,
        label="ShowLogFile:   open Latex log file by editor")
   tkadd(mbRevweb.menu,"command", command=DvipdfReport,
        label="DvipdfReport:   translate dvi- in pdf-file")
   tkadd(mbRevweb.menu, "separator")
   tkadd(mbRevweb.menu,"command", command=ProcessWithSexpr,
        label="ProcessWithSexpr:   save, weave source file, eval \\Sexpr{...}, latex")
   tkadd(mbRevweb.menu,"command", command=WebReport,
        label="WebReport:   save, weave and tangle source file")
   tkadd(mbRevweb.menu,"command", command=WeaveReport,
        label="WeaveReport:   save, weave source file")
   tkadd(mbRevweb.menu,"command", command=WeaveReportNoCode,
        label="WeaveReportNoCode:   save, weave source file, hide code")
   tkadd(mbRevweb.menu,"command", command=WeaveReportNoText,
        label="WeaveReportNoText:   save, weave source file, hide text")
   tkadd(mbRevweb.menu,"command", command=WeaveReportEchoCode,
        label="WeaveReportEchoCode:   save, weave source file, hide code if code name contains echo=FALSE")
   tkadd(mbRevweb.menu,"command", command=TangleReport,
        label="TangleReport:   save, tangle source file")
   tkadd(mbRevweb.menu,"command", command=TangleReportChunk,
        label="TangleReportChunk:   save, tangle source file ask for root chunk")
   tkadd(mbRevweb.menu,"command", command=TangleReportNoComments,
        label="TangleReportNoComments:   save, tangle source file without added comment lines")
   tkadd(mbRevweb.menu, "separator")
   tkadd(mbRevweb.menu, "command", command=RunAllandIncludeResults,
        label="RunAllandIncludeResults:   simple Sweave: <filename>-out.tex and LaTeX file")
   tkadd(mbRevweb.menu, "command", command=RunBeginEndEnvandIncludeResults,
        label="RunBeginEndEnvandIncludeResults:   simple Sweave of environment: <filename>-out.tex and LaTeX file")
   tkadd(mbRevweb.menu,"command", command=SWEAVE,
        label="SWEAVE:   save, Sweave and LaTeX file")
   tkadd(mbRevweb.menu,"command", command=SWEAVEB,
        label="SWEAVE:   save, Sweave (new process) and LaTeX file")
   ##tkadd(mbRevweb.menu, "separator")
   FTL.ok<-FALSE
   if( (.Platform$OS.type=="windows") ){
         libpath<-file.path(relax.path,"lib")
         if(file.exists(file.path(libpath,"pnmcrop.exe")) &&
            file.exists(file.path(libpath,"ppmtojpeg.exe")) &&
            file.exists(file.path(libpath,"netpbm.dll")) &&
            file.exists(file.path(libpath,"libjpeg.dll"))) FTL.ok<-TRUE
   }
   if(substring(version$os,1,5)=="linux"){
      if(0<length(system("which convert",TRUE,TRUE)) &&
         0<length(system("which dvips",TRUE,TRUE))) FTL.ok<-TRUE
     }
   if(FTL.ok){
     tkadd(mbRevweb.menu, "separator")
     tkadd(mbRevweb.menu, "command", command=FormatTeXLines,
           label="FormatTeXLines:   format text chunk by LaTeX and insert  jpeg-file")
   }

   ##erstelle Menüeintrag für Funktionstasten-Puffer##
   ## tkadd(mbRevweb.menu, "separator")
   ## tkadd(mbRevweb.menu,"command", command=LoadRwtools,
       ## label="LoadR-wtools:   load some functions: stem.leaf, ... and some data sets")
  }

  tkadd(mbFile.menu,"command", command=SetWorkPath,
        label="SetWorkPath:   change working path")
  tkadd(mbFile.menu, "separator")
  tkadd(mbFile.menu,"command", command=OpenReport,
        label="OpenReport:   >>APPEND<< file to text field")
  if(but.Wizardry!="simple") 
    tkadd(mbFile.menu,"command", command=SaveReport,
          label="SaveReport:   save text field to rev file")
  if(but.Wizardry!="simple") 
    tkadd(mbFile.menu, "separator")
  tkadd(mbFile.menu,"command", command=SaveHtml,
        label="SaveHtml:   save text field to a rev file and a HTML file")

  tkadd(mbFile.menu,"command", command=ViewReport.html,
        label="ViewReport.html:   view html representation of report")
  if(but.Wizardry!="simple") {
    tkadd(mbFile.menu,"separator")
    tkadd(mbFile.menu,"command", command=OpenTextFile,
          label="OpenTextFile:  open not-rev-file, translate console style and >>APPEND<< it to text field")
    tkadd(mbFile.menu,"command", command=SaveAsConsoleStyleFile,
          label="SaveAsConsoleStyleFile:   dump report as txt-, tex-file in console style, and as a rev- , R-file")
  #tkadd(mbFile.menu,"command", command=SaveAsPlainTeXFile,
  #      label="SaveAsPlainTeXFile:   dump report as TeX-file in console style, as a rev-file and as an 
  #      R-file")
  }
  tkadd(mbFile.menu,"separator")
  tkadd(mbFile.menu, "command", command=PLAYGROUND,
        label="playground:   open R-window for testing R-code")
  tkadd(mbFile.menu,"command", command=ShowHistory,
        label="ShowHistory:   show history of evaluation in a separate window")
  tkadd(mbFile.menu,"command", command=SaveHistory,
        label="SaveHistory:   save history of evaluation")
  if(but.Wizardry!="simple") {
    tkadd(mbFile.menu,"command", command=StartCodeChunkPlayer,
        label="StartCodeChunkPlayer:     call code chunk player (without saving report field)")
    tkadd(mbFile.menu,"command", command=ConstructDemoFunction,
        label="ConstructDemoFunction:   save, construct demo showing code chunks")
    tkadd(mbFile.menu, "separator")

    tkadd(mbFile.menu,"command", command=LoadEnvironment,
          label="LoadEnvironment:   load objects from dump file into relax environment")
    tkadd(mbFile.menu,"command", command=DumpEnvironment,
          label="DumpEnvironment:   save objects of relax environment as dump file")
    tkadd(mbFile.menu,"command", command=SaveEnvironment,
          label="SaveEnvironment:   save objects of relax environment in binary file")
    tkadd(mbFile.menu,"command", command=CleanEnvironment,
          label="CleanEnvironment:   delete objects of environment")
  }
  tkadd(mbFile.menu,"separator")

  if(but.Wizardry=="all"){
    tkadd(mbFile.menu,"command", command=OpenRevbook,
        label="OpenCompbook:   load compbook from library relax/rev")
    tkadd(mbFile.menu, "separator")
  }
  tkadd(mbFile.menu, "command", label="Exit:   quit RELAX",
        command=Exit)


  tkadd(mbEdit.menu, "separator")
  if(but.Wizardry!="simple") {
    tkadd(mbEdit.menu, "command", command=EditReport,
          label=paste("EditReport:   use editor",editor.sys,"for editing"))
  }
  tkadd(mbEdit.menu, "command", command=DumpCodeChunk,
        label=paste("DumpCodeChunk:   save code chunk in a file"))
  if(but.Wizardry!="simple") {
    tkadd(mbEdit.menu, "command", command=CopyToEnd,
          label=paste("CopyToEnd:   copy output to end of text"))
    tkadd(mbEdit.menu, "separator")
    tkadd(mbEdit.menu, "command", command=RunAll,
          label="RunAll:   run all start- and *-code chunks")
    tkadd(mbEdit.menu, "command", command=RunStart,
          label="RunStart:   run all start-code chunks")
    tkadd(mbEdit.menu, "separator")
    tkadd(mbEdit.menu, "command", command=RemoveSavedPlot,
          label="RemoveSavedPlot:   remove saved plots defined by cursor and plot links from REPORT FIELD")
    tkadd(mbEdit.menu, "command", command=RemoveALLPLOTS,
          label="RemoveALLPLOTS:   remove ALL saved plots of the report from file system and report")
    tkadd(mbEdit.menu, "command", command=DeleteAll,
          label="DeleteAll:   clear text field")
    tkadd(mbEdit.menu, "separator")
  }
  tkadd(mbEdit.menu, "command", command=UnDo,
        label="UnDo:   load report-UnDo-bak.rev (report before last evaluation)")

  melde("frame head defined",3)
  implement.but("TrashROutput", "fworkcmds", "trash output (after line of cursor)",
                job=fTrashROutput) ## 071115
  implement.but("WarnEval",   "fworkcmds", "evaluate code (even if warnings occur)",
                job=fWarnEval)
  implement.but("EvalRCode",   "fworkcmds", "evaluate code (stop in case of warnings)",job=fEvalRCode)
  implement.but("PlanRCode",    "fworkcmds", "insert new empty code chunk",job=fPlanRCode)

  tkpack(tklabel(fworkcmds, text="Line:",
         pady="-3"
,font="-Adobe-helvetica-Medium-R-Normal--12-140-*" # ,foreground="#124800"
),
         # ipady=7,
         side="left")
  llineno.sys<-tklabel(fworkcmds, text=" ",width="4") # 
  tkpack(llineno.sys, side="left")                    ###100208
  lworkname.sys<-tklabel(fworkcmds, text=paste(workname.sys,""),
                      pady="-3"
 ) ### (, relief="ridge")
  tkpack(lworkname.sys, side="left")

  implement.but("Down",   "fworkcmds", "jump downwards",side="left",bwf=0.40,job=fDown)
  implement.but("Up",   "fworkcmds", "jump backwards",side="left",bwf=0.40,job=fUp)

  implement.but("RemoveOut",foutcmds,"clear output field",job=fRemoveOut) ## 071115
  implement.but("FindText",foutcmds,"find text string in report field",job=fFindText)
  implement.but("SavePlot",foutcmds,
            "save plot as ppm, jpg and postscript file",job=fSavePlot)
  implement.but("Insert",foutcmds,"insert output into report",job=fInsert)

  ### tkpack(tklabel(foutcmds, text=" Result(s):",pady=#<negativer y-Platz>#,
  ###     font=#<Font für Knöpfe>#), side="left")  #ipady=7,  ###100208


           ##implementiere Testknopf##
  tworkwin<-tktext(fworkwin, background="#f7fffF", font=tfont.sys) #cff#d0f0ff#c81e6ff23#f0fffF
  try(tkconfigure(tworkwin, undo=1)) # 050704
  workbar<-tkscrollbar(fworkwin)
  tkconfigure(tworkwin,yscrollcommand=function(...) tkset(workbar,...))
  tkconfigure(workbar, command=function(...) tkyview(tworkwin,...))
  tkpack(workbar ,side="right",fill="y")
  tkpack(tworkwin,fill="both",expand="yes")
  tkinsert(tworkwin,"0.0",paste("% New Report:",date(),"\n"))
  ### CLIPBOARD-Funktion PASTE a la windows / clipboard
    ### tkevent.add("<<Paste>>", "<Control_L><v>")
  if(substring(version$os,1,6)=="darwin"  ){
    mac.paste<-function(...){  ## Ctrl-V  # mac-PASTE
       try({.Tcl("clipboard append hello"); .Tcl("clipboard clear")
            news<-base::scan(file=pipe("pbpaste","r"),what="",sep="\n",blank.lines.skip=FALSE)
            tkinsert(tworkwin,"insert",paste(news,collapse="\n"))})
            tksee(tworkwin,"insert - 7 lines");  tksee(tworkwin,"insert + 7 lines") #090706 
    }  #  for shifting view:  tclvalue(tkyview(tworkwin,"scroll","-1","pages"))
    tkbind(tworkwin,"<Control_L><v>",mac.paste)
    mac.copy<-function(...){  ## Ctrl-C  # mac-COPY
       news<-""
       try(news<-tclvalue(.Tcl("if {[catch {clipboard get}]} {set aa empty} {set aa full}")))
       if(news=="empty") return()
       try({news<-tclvalue(.Tcl("set aaa [selection get -selection CLIPBOARD]"))
            base::cat(news,file=get("tmp.file.name",envir=revive.sys)
); system(paste("pbcopy < ",get("tmp.file.name",envir=revive.sys)
))
            .Tcl("clipboard append hello"); .Tcl("clipboard clear")})
    }
    tkbind(tworkwin,"<Control_L><c>",mac.copy)
    tkevent.add("<<extract>>",  "<Control_L><c><KeyRelease>")   # mac-extract
    tkbind(tworkwin,"<<extract>>",mac.copy)
  }else{
    tkevent.add("<<Paste>>",   "<Control_L><v>")
    tkbind(tworkwin,"<<Paste>> { catch {%W insert insert [selection get -selection CLIPBOARD] } }")
  }
  toutwin<-tktext(foutwin,height=8,background="#ffffee", font=outfont.sys) #fff080 #fc8f16a23#ffffcc
  outbar<-tkscrollbar(foutwin)

  tkconfigure(toutwin,yscrollcommand=function(...) tkset(outbar,...))
  tkconfigure(outbar,        command=function(...) tkyview(toutwin,...))
  tkpack(outbar ,side="right",fill="y")
  tkpack(toutwin,fill="both",expand="yes")

  if((.Platform$OS.type=="windows")) #{}
  { # 110505
  # fout<-get("fout",envir=revive.sys); toutwin<-get("toutwin",envir=revive.sys)
   tkbind(fout,"<Configure>","")
   config.fns<-function(...){
    if("1"==tclvalue(tkwinfo("ismapped",toutwin))){
        tkpack("forget",toutwin) #; tkbind(toutwin,"<Configure>","")
        Sys.sleep(.01)
        if("0"==tclvalue(tkwinfo("ismapped",toutwin))){
          tkpack(toutwin,fill="both",expand="yes")
        }
    }
   }
   tkbind(fout,"<Configure>",config.fns)
  }
  if(substring(version$os,1,6)=="darwin"  ){  #121109
    tkbind(toutwin,"<Control_L><v>",mac.paste)
    tkbind(toutwin,"<Control_L><c>",mac.copy) 
    tkbind(toutwin,"<<extract>>",mac.copy) 
  }

  ( # def. namenslose Funktion:
  function(tworkwin){
   tktag.configure(tworkwin,"tld",   foreground="#C8162315126C", relief="raised",
                   borderwidth="2") # alternativ: #aaa222111
   tktag.configure(tworkwin,"tex", background="#e9ffff", 
                   relief="raised", borderwidth="3")
   tktag.configure(tworkwin,"code", background="#ffffff", 
                   foreground="#d21", font=outfont.sys) 
    # ddd222222
   tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
   tktag.configure(tworkwin,"emph", background="#999999999")
   tktag.configure(tworkwin,"jpeg", background="#ffffee",borderwidth="2",relief="raised")
   tktag.configure(tworkwin,"atsign", borderwidth="1",
                #   font="-Adobe-helvetica-Medium-R-Normal--10-100-*",
                   relief="raised"  # ,spacing1="2",spacing3=3,lmargin1=2
   )
   proc<-paste(
    "proc emphline {w mustera musterz} {",
      "scan [$w index end] %d anzzeilen", "set emphline 0\n",
      "for {set i 1} {$i < $anzzeilen} {incr i} {",
        "set actline [$w get $i.0 $i.end]",
        "if {[regexp $mustera $actline]} {",
          "if {[regexp $musterz $actline]} {",
            "uplevel $w tag add emph $i.0 $i.end",
          "}",
        "}",
      "}",
    "}", # type-semantics: 0=text, 1=tld, 2=code, 3=output
    "proc marklinetypes {w} {",
      "scan [$w index end] %d anzzeilen",
      "set type 0\n",
      "for {set i 1} {$i < $anzzeilen} {incr i} {",
       "set zeile [$w get $i.0 $i.end]",
       "set iv [ expr $i-1 ]","set in [ expr $i+1 ]",
       "if {$type==1}                           {\n set type 2\n}",
       "if {[regexp \"^@\"             $zeile]} {\n set type 0\n}",
       "if {[regexp \"^<<.*>>=\"     $zeile]} {\n set type 1\n}",
       "if {[regexp \"^>\"     $zeile]} {\n set type 11\n}", # 091029
       "if {[regexp \"^\\\\\\\\s(u|e)(b|c)\" $zeile]} {\n set type 4\n}",#121127  
       "if {[regexp \"^......img src\"  $zeile]} {\n set type 5\n}",
       "if {[regexp \"^output-start\"  $zeile]}  {\n set type 3\n}",
       "if {[regexp \"^\\\\\\\\begin\\{verbatim\\}\" $zeile]} {\n set type 3\n}",
   "if {$type==0} {\n uplevel $w tag add atsign $i.0 $in.0 \n}",
   "if {$type==4} {\n uplevel $w tag add atsign $i.0 $in.0 \n}",
       # "if {$type==2} {\n uplevel $w tag add code    $i.0    $i.end \n}",
       # "if {$type==2} {\n uplevel $w tag add code    $iv.end $i.end \n}",
       "if {$type==2} {\n uplevel $w tag add code    $iv.end $in.0 \n}",
       # "if {$type==1} {\n uplevel $w tag add tld     $i.0 $i.end \n}",
       "if {$type==1} {\n uplevel $w tag add tld     $i.0 $i.end ",
                                                   "; set type 2 \n}",
       "if {$type==11} {\n uplevel $w tag add tld     $i.0 $i.1 ", # 091029
                       ";  uplevel $w tag add code    $i.1 $i.end ",
                                                   "; set type 2 \n}",
       ### "if {$type==1} {\n uplevel $w tag add tld     $i.0 $i.end \n ; ",
       ###               " \n uplevel $w tag add code $i.end $i.end \n ; ",
       ###                                            "; set type 2 \n } ",
       "if {$type==3} {\n uplevel $w tag add output  $i.0 $i.end \n}",
       "if {$type==4} {\n uplevel $w tag add tex     $i.0 $i.end","; set type 0 \n}",
       #"if {$type==4} {\n uplevel $w tag add tex     $i.0 $in.0", "; set type 0 \n}",
       "if {$type==5} {\n uplevel $w tag add jpeg    $i.0 $i.end","; set type 0 \n}",
       "if {[regexp \"^output-end\"  $zeile]} {\n set type 0\n}",
       "if {[regexp \"^\\\\\\\\end\\{verbatim\\}\" $zeile]} {\n set type 0\n}",
       "$w tag raise sel",
      "}",
    "}",
    "proc markclear w {",
      "$w tag remove tld    1.0 end",
      "$w tag remove tex    1.0 end",
      "$w tag remove code   1.0 end",
      "$w tag remove output 1.0 end",
      "$w tag remove atsign 1.0 end",
    "}", sep="\n")
    .Tcl(proc)
   tktag.bind(tworkwin,"jpeg","<Enter>", #081121
              function() set.tclvalue("tvmess","press RETURN to display jpeg by browser!"))
   tktag.bind(tworkwin,"jpeg","<Leave>",function(){ set.tclvalue("tvmess","relax") })
   ## tktag.bind(tworkwin,"jpeg","<ButtonRelease>",
   tktag.bind(tworkwin,"jpeg","<Return>",
     function(){
       line <-floor(as.numeric(tkindex(tworkwin,"insert")))

       res<-tkmessageBox(message=
               if(language=="german") 
                  paste("jpeg-Datei im Browser angezeigen?")
               else paste("show jpeg graphics by browser"),
               title="Display JPEG?",icon="warning",type="yesnocancel",default="yes")
       res<-tclvalue(res); if(res=="cancel"||res=="no") return()
       fname<-tclvalue(tkget(tworkwin,paste(line,"0",sep="."),paste(line,"40",sep=".")))
       fname<-sub("^.*img src..","",fname)
       fname<-sub("(.jpg).*$","\\1",fname)
       browser.sys<-get("browser.sys",envir=revive.sys)
       cat(fname, "will be displayed by browser in a few seconds!\n")
       if( (.Platform$OS.type=="windows") ){
         browser.exe<- if(browser.sys=="") "start " else browser.sys
             res<-try(shell(paste(browser.exe, fname),wait=FALSE))  
         if(res!=0){ cat("ERROR: browser hasn't been started successfully \n") }
       } else {
         if(browser.sys!=""){
                try(system(paste(browser.sys,fname),wait=FALSE))
         }
       }
       "relax"
    }
  )
  }  # und Aufruf der namenslosen Funktion:
  )(tworkwin)

  .Tcl("         
    # http://wiki.tcl.tk/15612
    proc searchrep {t {replace 1}} {
      set w .sr
      if ![winfo exists $w] {
          toplevel $w
          wm title $w \"Search\"
          grid [label $w.1 -text Find:] [entry $w.f -textvar Find] [
             button $w.bn -text Next -command [list searchrepnext $t]
                   ] -sticky ew
          bind $w.f <Return> [list $w.bn invoke]
          if $replace {
             grid [label $w.2 -text Replace:] [entry $w.r -textvar Replace] [
                button $w.br -text Replace -command [list searchreprep1 $t] 
                   ] -sticky ew
             bind $w.r <Return> [list $w.br invoke]
             grid [label $w.3 -text \"Please:\"] [
                   label $w.4 -text \" use <Next> first!\"] [
                button $w.ba -text \"Replace all\" -command [list searchrepall $t]
                   ] -sticky ew
          }
          grid x [checkbutton $w.i -text \"Ignore case\" -variable IgnoreCase] [
                button $w.c -text Cancel -command \"destroy $w\"
                   ] -sticky ew
          grid $w.i -sticky w
          grid columnconfigure $w 1 -weight 1
             $t tag config hilite -background yellow
      } else {raise $w}
   }
  #-- Find the next instance
   proc searchrepnext w {
       foreach {from to} [$w tag ranges hilite] {
           $w tag remove hilite $from $to
       }
       set cmd [list $w search -count n -- $::Find insert+2c]
       if $::IgnoreCase {set cmd [linsert $cmd 2 -nocase]}
       set pos [eval $cmd]
       if {$pos ne \"\"} {
           $w mark set insert $pos
           $w see insert
           $w tag add hilite $pos $pos+${n}c
       }
   }
  #-- Replace the current instance, and find the next
   proc searchreprep1 w {
       if {[$w tag ranges hilite] ne \"\"} {
           $w delete insert insert+[string length $::Find]c
           $w insert insert $::Replace
           searchrepnext $w
           return 1
       } else {return 0}
   }
  #-- Replace all
   proc searchrepall w {
       set go 1
       while {$go} {set go [searchreprep1 $w]}
   }
  ")

  ## linfo.name <- tklabel(finfo, text=" ", font=#<Font für Knöpfe>#) ##
  linfo      <- tklabel(finfo, text=" ", textvariable="tvmess", relief="ridge")
    ###                      width=#<Infofeldbreite>#,
  # tkpack(linfo.name,linfo,side="left",fill="x",expand="yes")
  tkpack(linfo,side="left",fill="x",expand="yes")

  #tkbind(linfo,"<ButtonPress>",#<generiere interaktives icon>#)
  einfo.tmp <- tkentry(finfo,textvariable="tvinfo",background="#7ff",
                       width=floor(as.numeric(60 ###100208 vorher 80
)*3/5))
  linfo.tmp <- tklabel(finfo)

  melde("frame info defined",3)
  #<entferne [[relax.fns]] nach Beendigung von [[relax]] aus Suchpfad># # 121214

  tkbind(TopW, "<<FindReportText>>", FindReportText)
  tkbind(TopW,"<F3>", function(){  # 050628
        set.tclvalue("tvmess","relax")
        if(!exists("string.sys",envir=revive.sys)) return()
        such.orig<-such<-get("string.sys",envir=revive.sys)
        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
        tworkwin<-get("tworkwin",envir=revive.sys)
        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
        if(nchar(worktext)<10000){
          worktext<-strsplit(worktext,"\n")[[1]]
        }else{
          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
        }

        repl.pat<-gsub("(.)","\\\\\\1","^!$%&/()=?{}}+*#,.-;:\\_[") ## 070830
        repl.pat<-paste("([",repl.pat,"])",collapse="")
        such<-gsub(repl.pat,"\\\\\\1",such)
        assign("string.sys",string.sys,envir=revive.sys)
        if(nchar(such)==0) {
          tcl("findclear",tworkwin); tkfocus(tworkwin); return()
        }

        if(length(found<-grep(such,worktext))>0){
          cline.pos <- as.character(tkindex(tworkwin,"insert"))
          cline <- floor(as.numeric(cline.pos))
          cline.pos<-as.numeric(sub(".*\\.","",cline.pos))

          # if nothing found in the cursor line or later go to the beginning of the document
          line <- if(any(h<-(cline<found))) found[h][1] else found[1]
          # look always for the first occurence, index origin: 1
          h<-worktext[line]; hh<-1:nchar(h)
          line.pos<- -1+which(such.orig==substring(h,hh,hh-1+nchar(such.orig)))[1]
          # if cursor line and match line are equal you have to look at the positions
          if(any(found==cline)){
            h<-worktext[cline]; hh<-1:nchar(h)
            found.pos<- -1+which(such.orig==substring(h,hh,hh-1+nchar(such.orig)))
            if(any(h<-(cline.pos<found.pos))) {line<-cline; line.pos<-found.pos[h][1]}
          }

          tksee(tworkwin,paste(line,".1",sep="")); h<-paste(line,".",line.pos,sep="")
          tkmark.set(tworkwin, "insert", h); tkfocus(tworkwin)
        } else set.tclvalue("tvmess",paste("Warning: search string >",
                                           such.orig,"< not found!!!"))
      } # end of function
  )
  tkbind(TopW, "<<GoToLine>>", GoToLine)

  tkbind(TopW, "<<ProcessReport>>", ProcessReport)

  # tkconfigure(Help.R,command=fHelp.R)
  if( ! substring(version$os,1,6)=="darwin"  ) tkbind(TopW, "<<Help.R>>", fHelp.R) #120111
  # tkconfigure(Up,command=fUp)
  tkbind(TopW, "<<Up>>", fUp)
  # tkconfigure(Down,command=fDown)
  if( ! substring(version$os,1,6)=="darwin"  ) tkbind(TopW, "<<Down>>", fDown) #120111

  # tkconfigure(PlanRCode,command=fPlanRCode)
  tkbind(TopW, "<<PlanRCode>>", fPlanRCode)
  # tkconfigure(EvalRCode,command=fEvalRCode) ## 071115
  tkbind(TopW, "<<EvalRCode>>", fEvalRCode)
  # tkconfigure(WarnEval,command=fWarnEval)
  tkbind(TopW, "<<WarnEval>>", fWarnEval)
  # tkconfigure(TrashROutput,command=fTrashROutput)
  tkbind(TopW, "<<TrashROutput>>", fTrashROutput)

  # tkconfigure(Insert,command=fInsert)
  tkbind(TopW, "<<Insert>>", fInsert)
  # tkconfigure(SavePlot,command=fSavePlot)
  tkbind(TopW, "<<SavePlot>>", fSavePlot)
  # tkconfigure(FindText,command=fFindText)
  tkbind(TopW, "<<FindText>>", fFindText)
  # tkconfigure(RemoveOut,command=fRemoveOut) ## 071115
  tkbind(TopW, "<<RemoveOut>>", fRemoveOut)

  #tkconfigure(CopyToEnd,command=CopyToEnd)
  #tkbind(TopW, "<<CopyToEnd>>", CopyToEnd)

  tkbind(TopW, "<<SaveReport>>", SaveReport)
   ##implementiere Eigenschaften vom Testknopf##
  twin<-tworkwin; f.sonderzeichen<-function(zeichen){
                    function(){
                      tkinsert(twin,"insert",zeichen)
                      if((.Platform$OS.type=="windows"))tkdelete(twin,"insert-1chars")
                      if(substring(version$os,1,6)=="darwin" )tkdelete(twin,"insert-1chars")
                    }
                  }
                  f.umlaut<-function(zeichen){
                    function(){
                      return()
                      #char337<-eval(parse(text='"\\337"'))
                      #if(zeichen==char337 & tclvalue(tkget(twin,"insert-1chars","insert"))=="\\") return()
                      #tkinsert(twin,"insert",zeichen); tkdelete(twin,"insert-2chars")
                    }
                  }
                  tkbind(twin,"<<LKeckig>>", f.sonderzeichen("["))
                  tkbind(twin,"<<RKeckig>>", f.sonderzeichen("]"))
                  tkbind(twin,"<<Tilde>>",   f.sonderzeichen("~"))
                  tkbind(twin,"<<LKgeschw>>",f.sonderzeichen("{"))
                  tkbind(twin,"<<RKgeschw>>",f.sonderzeichen("}"))
                  tkbind(twin,"<<Klammera>>",f.sonderzeichen("@"))
                  tkbind(twin,"<<Pipe>>",    f.sonderzeichen("|"))
                  tkbind(twin,"<<Backsl>>",  f.sonderzeichen("\\"))
                  renewhighlighting<-function(){
                    tworkwin<-get("tworkwin",envir=revive.sys)
                    melde("ak texthervor",1)
                    tcl("markclear",tworkwin)
                    tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
                    tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
                    tcl("marklinetypes",tworkwin)
                    melde("ak texthervor",2)

                  }
                  # tkbind(twin,"<<Klammeraffe>>",renewhighlighting)
                  tkbind(twin,"<Return>",renewhighlighting)

  ##definiere Wirkung von [[Strg Pagedown]] im Reportfenster## ab 1.02 abgeschaltet
  proc<-c(  # 060704
     "proc showbracket {w orta ortb} {",
     "$w tag remove secondbracket    1.0 end",
     "  scan [$w index end] %d numLines",
     "     $w mark set first \"$orta\"",
     "     $w mark set last \"$orta  + 1 chars\"",
     "     uplevel [$w tag add secondbracket first last]",
     "     $w mark set first \"$ortb\"",
     "     $w mark set last \"$ortb  + 1 chars\"",
     "     uplevel [$w tag add secondbracket first last]",
     "}",
      "proc showbracketsclear w {",
        "$w tag remove secondbracket    1.0 end",
      "}"
    )
  .Tcl(paste(proc,collapse="\n"))
  tktag.configure(tworkwin,"secondbracket",background="#000dddfff",relief="raised")
  MarkSecondBracket<-function(){
    tworkwin<-get("tworkwin",envir=revive.sys)
    left.of.cursor<-tclvalue(tkget(tworkwin,"insert - 1 chars"))
    open<-0
    if(left.of.cursor %in% c("{","[","(")){ ## print("Klammer auf")
      open<-1
      a<-tclvalue(tkget(tworkwin,"insert - 1 chars","end"))
      n<-min(nchar(a),1500); aa<-substring(a,1,n); aaa<-substring(aa,1:n,1:n)
      auf<-aaa %in% c("{","[","("); zu<-aaa %in% c("}","]",")")
      count<-cumsum(auf)-cumsum(zu)
      second<-which(count==0)[1]-2
      first<-paste(tclvalue(tkindex(tworkwin,"insert")),"- 1 chars")
      if(!is.na(second)){
        tcl("showbracket",tworkwin,first,paste("insert +",second,"chars"))
      }
    }
    if(left.of.cursor %in% c("}","]",")")){ ## print("Klammer zu"); 
      open<- -1
      a<-tclvalue(tkget(tworkwin,"0.0","insert"))
      n<-nchar(a); aa<-substring(a,max(1,n-1500),n)
      aaa<-substring(aa,nchar(aa):1,nchar(aa):1) 
      auf<-aaa %in% c("{","[","("); zu<-aaa %in% c("}","]",")")
      count<-cumsum(zu)-cumsum(auf)
      second<-which(count==0)[1]
      first<-paste(tclvalue(tkindex(tworkwin,"insert")),"- 1 chars")
      if(!is.na(second)){
        tcl("showbracket",tworkwin,first,paste("insert -",second,"chars"))
      }
    }
    ## if(exists("second")) print(second)  ##; print(open)
    if(open==0) tcl("showbracketsclear",tworkwin)
    ### Zeilennummeraktualisieren: 
    line <-floor(as.numeric(tkindex(tworkwin,"insert")))

    tkconfigure(llineno.sys,text=paste(line))

  }
  tkbind(tworkwin,"<KeyRelease>",MarkSecondBracket)
  ### Zeilennummeraktualisieren falls Mouse-Op: 
  tkbind(tworkwin,"<ButtonRelease>",function(){
      tworkwin<-get("tworkwin",revive.sys)
      line <-floor(as.numeric(tkindex(tworkwin,"insert")))

      tkconfigure(llineno.sys,text=paste(line))
    }
  )

   tkbind(tworkwin,"<Tab><KeyRelease>", function() {
    #  tworkwin<-get("tworkwin",envir=revive.sys)
   ## check, ob Leerzeichen oder nicht vorm Cursor steht
    if(1==length(grep("[a-zA-Z._\\(]", #)
                      tclvalue(tkget(tworkwin,"insert - 2 char")))))
       try(tkdelete(tworkwin,"insert - 1 char")) else return()
    if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
    tworkwin<-get("tworkwin",envir=revive.sys)
    worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
    if(nchar(worktext)<10000){
      worktext<-strsplit(worktext,"\n")[[1]]
    }else{
      base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
      worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
    }

    line<-line.orig<-(tclvalue(tkindex(tworkwin,"insert")))
    pos<-as.numeric(sub("(.*)\\.","",line)); line<-as.numeric(line) 
    w<-worktext[floor(line)]
    anfang<-paste(rev(substring(w,1:pos,1:pos)),collapse="") # aktuelle Zeile 
    anfang<-sub("(^[\\(]{0,1}[a-zA-Z0-9._]*)(.*)","\\1",anfang) # aktuelles Wort
    h<-nchar(anfang); klammer<-"("==substring(anfang,1,1)
    if(klammer && h<2) return()
    anfang<-if(klammer) paste(rev(substring(anfang,2:h,2:h)),collapse="") else
                        paste(rev(substring(anfang,1:h,1:h)),collapse="")
    suchstring<-paste("^",anfang,sep="")

    objekt<-grep(suchstring,ls(envir=revive.env),ignore.case=FALSE,value=TRUE)
    objekt<-c(objekt,apropos(suchstring,ignore.case=FALSE))
    if(length(objekt)==0) return();  if(length(objekt)>40) return()

    if(1==length(objekt)||klammer) { 
      if(klammer) objekt<-substring(suchstring,2) # ohne ^
      objtail<-sub(suchstring,"",objekt)
      objtail<-substring(objekt,nchar(suchstring)) # because of "^", 1+(nchar(.)-1)
      if(0==length(find(objekt))) h<-get(objekt,envir=revive.env) else h<-get(objekt)
      # cat("h"); print(h); cat("objtail"); print(objtail) # der angefuegt werden muss
      if(is.function(h)){
        h<-paste(deparse(args(h)),collapse=" ")[1]
        h<-sub("^function[^\\(]*\\(","(",h);h<-sub(" *NULL$","",h)
        news<-paste("\n\nsyntax of function ",objekt,":\n",h,sep="")
        if(!exists("toutwin"))
          toutwin<-get("toutwin",envir=get("revive.sys",envir=revive.env))
        pos.to.insert<-"end"
        ## news<-paste(gsub("\n+","\n",news),collapse="\n") # 111121
        news<-paste(gsub("\n[ \n\t]*\n","\n",news),collapse="\n")
        try(tkinsert(toutwin,pos.to.insert,news))
        tksee(toutwin,"end - 0 lines")
        melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

      }
      try(tkinsert(tworkwin,line.orig,objtail));    return()
    }
    ooo<-sapply(objekt,function(x) substring(x,1:20,1:20))
    n.equal<-which.min(apply(ooo[,1]==ooo,1,all))-1
    if(n.equal>nchar(anfang)){
      objekt<-paste(ooo[1:n.equal,1],collapse="")
      objtail<-sub(suchstring,"",objekt)
      try(tkinsert(tworkwin,line.orig,objtail));    return()
    } else { print(objekt); return() }
   })
   if( name.complete.sys!=TRUE ) {
     tkbind(tworkwin,"<Tab><KeyRelease>", function() { "relax" } ) 
   }

  ##definiere Funktionstasten-Puffer##
  Execute.cmds<-function(){
    melde("Execute.cmds",1)
    res<-ls(pattern="^cmds$",envir=revive.env)
    if(length(res)>0){
      cmds<-try(eval(parse(text="cmds"),envir=revive.env))
      if(length(cmds)>0){
        eval(parse(text="cmds<-NULL"),envir=revive.env)
        #repeat{
          melde("beginn repeat",1)
          if(length(cmds) == 0) break
          cmd <- substring(cmds[1],1,1); choice<-substring(cmds[1],2)
          cmds<-cmds[-1]; assign("cmds",cmds,envir=revive.env)
          switch(cmd
            ,"s" = {
                     if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
                     tworkwin<-get("tworkwin",envir=revive.sys)
                     worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
                     if(nchar(worktext)<10000){
                       worktext<-strsplit(worktext,"\n")[[1]]
                     }else{
                       base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
                       worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
                     }

                     line<-grep("^<<(.*)>>=",worktext)
                     if(class(try(no<-line[as.numeric(choice)]))!="try-error" && !is.na(no)){
                       line<-paste(no[1],"0",sep=".")
                       tkmark.set(tworkwin,"insert",line);fWarnEval()
                       tksee(tworkwin,"end")
                     }
                    }
            ,"q" = { Exit() }
            ,">" = {
                    news<-paste("\n@\n<<*>>=\n",choice,"",sep="")
                    if(!exists("tworkwin"))
                      tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

                    pos.to.insert<-"end"
                    if(0<length(grep("output-start",news))){
                      tail<-rev(strsplit(tclvalue(tkget(tworkwin,"end - 3 lines","end")),"\n")[[1]])
                      ltail<-length(tail)
                      if( (0==length(grep("<<[*]>>=",tail[1:ltail]))) &&
                         any(h<-("output-end"==substring(tail[1:ltail],1,11)))){
                         news<-sub(".*output-start\n","",news)
                         news<-sub("output-end","",news)
                         h<-seq(along=h)[h][1]
                         pos.to.insert<-paste("end -",h,"lines")
                      }
                    }
                    try(tkinsert(tworkwin,pos.to.insert,paste(news,collapse="\n")))
                    tksee(tworkwin,"end - 0 lines")
                    melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

                     tkmark.set(tworkwin,"insert","end");fWarnEval()
                     tksee(tworkwin,"end")
                   }
            ,"p" = {
                    bildname <- choice
                    f<-function() { 
                      if(!is.null(bildname)&&nchar(bildname)>0){
                      # check name of picture
                        n<-nchar(bildname<-gsub(" ","",bildname))
                        bildname<-sub(".ps$","",bildname)
                      # postscript:
                        psname <-paste(bildname,".ps", sep="")
                        try.res<-try({dev.copy(postscript,psname,horizontal=pshorizontal.sys,
                                               width=psdesignwidth.sys,height=psdesignheight.sys);dev.off()})
                        if(is.function(try.res)){
                          ok <- "OK"
                        } else {
                          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
                          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
                          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
                          if(!is.character(ok)) { ok <- "OK" }
                        }
                        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
                          ok<-FALSE
                          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
                          cat(error.msg,"\n")
                          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                             cat("A warning message stopped the evaluation!",
                                   "If you want to\nevaluate the code anyway",
                                   "evaluate code by:\n>WarnEval<")
                          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
                        } else { ok<-TRUE }


                        if(!ok){ cat("Error: *ps file not generated by dev.copy!!!\n"); return() }
                        news<-paste("@\n \\begin{center}","\\includegraphics[",
                                    "height=",psheight.sys,"]{",bildname,"}\\end{center}\n",sep="") #081121
                      # jpeg:
                        jpgname<-paste(bildname,".jpg",sep="")
                        if((.Platform$OS.type=="windows")){ # width=width in pixel, 72 dpi
                          try.res<-try({dev.copy(jpeg,jpgname,width=jpgdesignsize.sys*72,
                                                 height=jpgdesignsize.sys*72,quality=100,pointsize=7);dev.off()})
                        }else{
                          try.res<-try({dev.copy(bitmap,type="jpeg",jpgname,
                               width=jpgdesignsize.sys,height=jpgdesignsize.sys);dev.off()})
                        }
                        if(is.function(try.res)){
                          ok <- "OK"
                        } else {
                          if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
                          if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
                          if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
                          if(!is.character(ok)) { ok <- "OK" }
                        }
                        if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
                          ok<-FALSE
                          error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
                          cat(error.msg,"\n")
                          if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                             cat("A warning message stopped the evaluation!",
                                   "If you want to\nevaluate the code anyway",
                                   "evaluate code by:\n>WarnEval<")
                          # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
                        } else { ok<-TRUE }


                        if(!ok) cat("Error: *jpg file not generated by dev.copy!!!\n")
                        news<-paste(news,'\n% <p><img src="',jpgname,'">\n@\n', sep="" )
                      # ppm
                        ppmname<-paste(bildname,".ppm",sep="")
                        # if( <das OS ist Windows>  && 2<=nchar(ghostscript) && !Img.package.found && !no.plots){
                        if((.Platform$OS.type=="windows") && 2<=nchar(ghostscript) && !no.plots){
                          try.res<-try({dev.copy(bitmap,type="ppmraw",ppmname,res=ppmresolution.sys);dev.off()})
                        }
                        # if( <das OS ist Linux> && 2<=nchar(ghostscript) && !Img.package.found && !no.plots){ }
                        if(substring(version$os,1,5)=="linux" && 2<=nchar(ghostscript) && !no.plots){ # 121113
                          try.res<-try({dev.copy(bitmap,type="ppmraw",ppmname,res=ppmresolution.sys);dev.off()})
                        }
                      # gif+ppm:
                        if(substring(version$os,1,6)=="darwin"  && !no.plots){
                          gifname<-paste(bildname,".gif",sep="")
                          try({dev2bitmap(type="ppmraw",ppmname,res=ppmresolution.sys); dev.off()}) #121218
                          try.res<-try({system(paste("convert",jpgname,gifname))})
                          if(is.function(try.res)){
                            ok <- "OK"
                          } else {
                            if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
                            if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
                            if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
                            if(!is.character(ok)) { ok <- "OK" }
                          }
                          if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
                            ok<-FALSE
                            error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
                            cat(error.msg,"\n")
                            if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                               cat("A warning message stopped the evaluation!",
                                     "If you want to\nevaluate the code anyway",
                                     "evaluate code by:\n>WarnEval<")
                            # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
                          } else { ok<-TRUE }


                          if(!ok) cat("Error: gif file not generated by dev.copy!!!\n")
                        }
                      }
 
                      # include links
                      if(!is.null(bildname)&&nchar(bildname)>0){
                        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
                        tworkwin<-get("tworkwin",envir=revive.sys)
                        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
                        if(nchar(worktext)<10000){
                          worktext<-strsplit(worktext,"\n")[[1]]
                        }else{
                          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
                          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
                        }

                        ##hole ggf. [[tworkwin]]>>
                        line <-floor(as.numeric(tkindex(tworkwin,"insert")))

                        ##lese Arbeitsfenster auf [[worktext]] ein>>
                        textstart<-grep("^@",worktext)-1; textstart<-textstart[textstart>=line][1]
                        codestart<-grep("^<<(.*)>>=",worktext)-1; codestart<-codestart[codestart>=line][1]
                        if(is.na(codestart))codestart<-Inf; if(is.na(textstart))textstart<-Inf
                        insertline<-if(codestart==textstart) NA else min(codestart,textstart)
                        anzrows<-length(unlist(strsplit(news,"\n")))
                        if(is.na(insertline)){
                            insertline<-"end"
                            try(tkinsert(tworkwin,"end","\n"))
                            try(tkinsert(tworkwin,"end",paste(news,collapse="\n")))
                            tkmark.set(tworkwin, "insert","end - 2 lines")
                            tksee(tworkwin,"end")  # paste(insertline+anzrows,"0",sep="."))
                            insertline<-length(worktext)
                        }else{
                          # in einem Text-Chunks muss ein Kl-Affe eingebaut werden.
                            if(length(grep("<<\\*>>=",news[1]))>0 && codestart < textstart) news<-c(news,"@\n")
                            try(tkinsert(tworkwin,paste(insertline+1,"0",sep="."),paste(news,collapse="\n")))
                            tkmark.set(tworkwin, "insert", paste(insertline+anzrows,"0",sep="."))
                            tksee(tworkwin,paste(insertline+anzrows,"0",sep="."))
                        }
                        ## melde(insertline)
                        melde("ak texthervor",1)
                        tcl("markclear",tworkwin)
                        tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
                        tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
                        tcl("marklinetypes",tworkwin)
                        melde("ak texthervor",2)

                        ##zeige Bilder im Textfenster an##
                        tkfocus(tworkwin)
                        melde("inserted characters: \n",3,substring(news[1:min(7,length(news))],1,80))

                        melde(paste("p", psname), "cmd.msg")
                      }
                    }; f()
                   }
            ,"r" = {
                    choice<-gsub(" ","",choice)
                    if(!file.exists(choice)){
                      cat("ERROR:",choice,"not found!!!\n")
                      ok<-FALSE
                    }else{
                      filename<-choice
                      try.res<-try(myscan(filename,"",sep="\n",blank.lines.skip=FALSE))
                      if(is.function(try.res)){
                        ok <- "OK"
                      } else {
                        if(mode(try.res)=="externalptr"||mode(try.res)=="environment") try.res<-"ok"
                        if(mode(try.res)=="S4") try.res <- "ok" else ok<-try.res[1]  ## 130304
                        if(is.null(ok) ||is.na(ok)|| is.name(ok) || is.list(ok) || is.numeric(ok)) ok <- "OK"
                        if(!is.character(ok)) { ok <- "OK" }
                      }
                      if(0!=length(ok)&&("Error"==substring(ok,1,5) | "Fehler"==substring(ok,1,6))){
                        ok<-FALSE
                        error.msg<-unclass(try.res); error.msg<-sub("options.warn.2.","",error.msg)
                        cat(error.msg,"\n")
                        if(0<length(grep("Warnung",error.msg))||0<length(grep("warning",error.msg)))
                           cat("A warning message stopped the evaluation!",
                                 "If you want to\nevaluate the code anyway",
                                 "evaluate code by:\n>WarnEval<")
                        # cat("sorry, operation failed in:",as.character(sys.call()),"!!!\n") # due to R-2.15.1
                      } else { ok<-TRUE }


                    }
                    if(ok){
                        workname.sys<-sub(paste(".*",.Platform$file.sep,sep=""),"",filename)
                        lworkname.sys<-get("lworkname.sys",envir=revive.sys)
                        tkconfigure(lworkname.sys,text=paste(workname.sys,""))
                        assign("workname.sys",workname.sys,envir=revive.sys)
                        # tkwm.title(TopW,paste("RELAX:",workname.sys))
                        tkwm.title(TopW,paste(if(but.Wizardry=="simple") "redit:" else "RELAX:",workname.sys))                      

                        try.res<-WinToTcl.read(try.res)
                          ## Eintrag mit Entfernung des bisherigen Inhalts:
                          ## worktext<-paste(try.res, collapse="\n")
                          ## <<schreibe [[worktext]] ins Arbeitsfenster>>
                        news<-c("",try.res)
                        if(!exists("tworkwin"))
                          tworkwin<-get("tworkwin",envir=get("revive.sys",envir=revive.env))

                        pos.to.insert<-"end"
                        if(0<length(grep("output-start",news))){
                          tail<-rev(strsplit(tclvalue(tkget(tworkwin,"end - 3 lines","end")),"\n")[[1]])
                          ltail<-length(tail)
                          if( (0==length(grep("<<[*]>>=",tail[1:ltail]))) &&
                             any(h<-("output-end"==substring(tail[1:ltail],1,11)))){
                             news<-sub(".*output-start\n","",news)
                             news<-sub("output-end","",news)
                             h<-seq(along=h)[h][1]
                             pos.to.insert<-paste("end -",h,"lines")
                          }
                        }
                        try(tkinsert(tworkwin,pos.to.insert,paste(news,collapse="\n")))
                        tksee(tworkwin,"end - 0 lines")
                        melde("appended characters: \n",3,substring(news[1:min(7,length(news))],1,80))

                        if(!exists("revive.sys")) revive.sys<-get("revive.sys",envir=revive.env)
                        tworkwin<-get("tworkwin",envir=revive.sys)
                        worktext<-tclvalue(tkget(tworkwin,"0.0","end"))
                        if(nchar(worktext)<10000){
                          worktext<-strsplit(worktext,"\n")[[1]]
                        }else{
                          base::cat(worktext,file=get("tmp.file.name",envir=revive.sys)
)
                          worktext<-myscan(file=get("tmp.file.name",envir=revive.sys)
,"",sep="\n",blank.lines.skip=FALSE)
                        }

                        line <-floor(as.numeric(tkindex(tworkwin,"insert")))

                        code.start<-grep("^<<(.*)>>=",worktext)
                        try(if(0<length(code.start)){ 
                               worktext[code.start]<-sub("^<<(.*)>>=(.*)","<<\\1>>=",worktext[code.start])
                               worktext[code.start]<-paste(worktext[code.start]," (",1:length(code.start),")",sep="")
                        })
                        if(length(worktext)>1) worktext<-paste(worktext,collapse="\n")
                        tkdelete(tworkwin,"0.0","end")
                        try(tkinsert(tworkwin,"0.0",paste(worktext,collapse="\n")))
                        tksee(tworkwin,"end")
                        melde("ak texthervor",1)
                        tcl("markclear",tworkwin)
                        tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
                        tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
                        tcl("marklinetypes",tworkwin)
                        melde("ak texthervor",2)

                        if(!no.plots) { exclude.plots(tworkwin); show.plots.again(tworkwin) }


                        tkmark.set(tworkwin, "insert", paste(line,"0",sep="."))
                        tksee(tworkwin,paste(line,"0",sep="."))
                        tkfocus(tworkwin)

                        RunStart()
                        melde("ak texthervor",1)
                        tcl("markclear",tworkwin)
                        tktag.configure(tworkwin,"output",foreground="#111222999", font=outfont.sys)
                        tktag.configure(tworkwin,"code",  foreground="#ddd222222", font=outfont.sys)
                        tcl("marklinetypes",tworkwin)
                        melde("ak texthervor",2)

                        if(!no.plots) { exclude.plots(tworkwin); createandshow.all.plots(tworkwin) }

                        melde(paste("r", filename), "cmd.msg")
                        Execute.cmds()
                    } else { cat("ERROR: file not found!!!\n") }
                   }
          )

        #}
      } else remove(list="cmds",envir=revive.env)
    }
    melde("Execute.cmds",2)
  }
  tkbind(TopW, "<<Acticmds>>", Execute.cmds)

  ##definiere Bindung zur Warnungsabarbeitung##
  melde("initialization of RELAX finished",3)
  if(!missing(file.name)){
     file.name<-as.character(substitute(file.name)) ## 071115
     file.name<-gsub("\\\\","/",file.name)
     if(0<length(grep("/",file.name))){
       path<-sub("^(.*)/(.*)","\\1",file.name)
       print(path)
       if(nchar(path)>0) try(setwd(path))
       file.name<-sub("^(.*)/","",file.name)
     }
     if(0==length(grep(".rev$",file.name))) file.name<-paste(file.name,"rev",sep=".")
     print(file.name)
     workname.sys<-file.name
     assign("cmds",paste("r",file.name),envir=revive.env)
     Execute.cmds()
  }
  if(0<length(cmds) && cmds[1]!="") {
    assign("cmds",cmds,envir=revive.env)
    Execute.cmds()
  }

  ##definiere Logik zum Eintrag der Zeilennummer##
  data.fns.menu()
  ReloadReportWidget() # to repair defect report widget
  cat( "relax 1.3.15 - 140310" ,"\n")
  if(language=="german"){
    cat("relax Initialisierung abgeschlossen!\nR-Editor wird erneut durch  relax()  gestartet!\n")
  }else{
    cat("initialisation of relax completed!\nrestart relax by: relax()\n enjoy and relax it!!\n")
  }
#  tkwait.variable("tvexit")  # version 1.082
  return()
}
red<-redit<-function(file.name){
  if( missing(file.name)) eval(parse(text="relax(but.Wizardry='simple')"),envir=.GlobalEnv) else {
    file.name<-as.character(substitute(file.name))
    eval(parse(text=paste("relax('",file.name,"',but.Wizardry='simple')",sep="")),envir=.GlobalEnv)
    "redit: starts relax with reduced otions"
  }
}

args.player <- function( fn.call, width.of.label = 20,
                               width.of.entry = 25, main ){
    act.call <- deparse(substitute(fn.call))
    act.call <- sub('^["]{0,1} *',"",act.call)
    act.call <- sub(' *["]{0,1}$',"",act.call)
    args.found <- TRUE
    if( 0 == length(grep("[()]", act.call)) ) { # case only name of a function
      args.found <- FALSE
      # use default methods if exists
      fn.name <- if(exists(h <- paste(act.call,"default",sep="."))) h else act.call
      if(!exists(fn.name)) return(paste("function",fn.name,"not found"))
      # get arguments
      act.call <- deparse( args(fn.name) )
      # remove NULL entry
      if("NULL"==act.call[ h <- length(act.call) ]) act.call <- act.call[-h]
      # substitute "function" by fn.name
      act.call[1] <- sub("^function", fn.name, act.call[1])
      # check ... argument
      # if(0<length(grep("[.]{3}",act.call))) act.call<-sub("[.]{3}"," ",act.call)
    }
    # combine multiline calls
    act.call <- paste(act.call,collapse=" ")
    # get function name
    fn.name <- sub("[(].*","",act.call) #)
    # get argument list
    fn.args <- sub("[^(]*[(](.*)[)].*","\\1",act.call)
    ## print(fn.args)

    ### escape inner ","
    # split string to characters
    n <- nchar(fn.args)
    fn.args <- substring(fn.args,1:n,1:n)
    # ( opens a new level, ) closes a level, "," of level 0 are separators
    state <- 0
    for( i in 1:n ) {
      if( fn.args[i] == "(" ){ state <- state + 1; next }
      if( fn.args[i] == ")" ){ state <- state - 1; next }
      if( 0 < state && fn.args[i] == "," ){
        fn.args[i] <- "cOmMa FoUnD"; next
      }
    } 
    # now a "," is a separator and we can split the string
    fn.args <- paste(fn.args,collapse="")
    fn.args <- strsplit(fn.args,",")
    # and reconstruct escaped commas
    fn.args <- gsub("cOmMa FoUnD",",",unlist(fn.args))
    ### remove blanks in top and end position
    fn.args <- gsub("^ *([^ ].*)","\\1",fn.args)
    fn.args <- gsub("(.*[^ ]) *$","\\1",fn.args)
    # arguments without defaults get a "="
    if(!args.found){
      idx <- grep("[=]", invert=TRUE, fn.args)
      fn.args[idx] <- paste(fn.args[idx],"=")
    }
    ## print(fn.args)

  # Entry Felder schreiben
  n.args <- length(fn.args)
  fn.arg.values <- fn.arg.names <- rep("  ",n.args)
  for( i in 1:n.args ){
    arg <- fn.args[i]
    if( 0 < length(grep("[=]", arg))) fn.arg.names[i] <- sub("[=].*","=",arg) 
    fn.arg.values[i] <- sub(".*[=]","",arg)
  }
  fn.arg.reset <- fn.arg.values <- sub("^ *([^ ].*[^ ]) *$","\\1",fn.arg.values)
  fn.arg.names <- gsub("[= ]","",fn.arg.names)
  ## print(fn.arg.values); print(fn.arg.names) #; return()

  # Fenster einrichten
  Top <- tktoplevel()
  tkwm.title(Top,if(missing(main)) paste("function caller:",fn.name) else main)
  tkwm.geometry(Top,paste("500x",25+25*n.args,"+10+15",sep=""))
  tkpack(arg.frame <- tkframe(Top),fill="x",expand=1) # args
  tkpack(org.frame <- tkframe(Top),fill="x",expand=1) # eval, reset, exit
  tkpack(a <- tkbutton(org.frame, text = "Exit", width=10, height=1,
                       command = function()tkdestroy(Top)),
                       side = "right")
  tkpack(reset.but <- tkbutton(org.frame, text = "Reset", width=10, height=1), side = "right")
  tkpack(eval.but  <- tkbutton(org.frame, text = "Eval", width=10, height=1),  side = "left")

  set.tclvalue <- function(name,value) tclvalue(name)<-as.character(value)
  for( i in 1:n.args ){
    # cat("erstelle arg frame i",i)
    tkpack(a <- tkframe(arg.frame),fill="x")
    tkpack(tklabel(a,text = paste(fn.arg.names[i],"="), 
                   width= width.of.label), side="left")
    vari <- paste("ARG",i,collapse="")
    tkpack(ent <- tkentry(a, textvariable = vari, width= width.of.entry), side="right")
    if(i == 1) tkfocus(ent)
    set.tclvalue(vari,as.character(fn.arg.values[i]))
  }

  tkconfigure(reset.but, command = function(){
    for( i in 1:n.args ){
      vari <- paste("ARG",i,collapse="")
      set.tclvalue(vari,as.character(fn.arg.reset[i]))
    }
  })

  eval.envir <- parent.frame()
  tkconfigure(eval.but, command = function(){
    act.args <- rep(" ",n.args)
    for( i in 1:n.args ){ 
      act.args[i] <- tclvalue( paste("ARG",i,collapse="") )
    }
    act.args <- gsub("^ *([^ ].*)","\\1",act.args)
    act.args <- gsub("(.*[^ ]) *$","\\1",act.args)
    # no args => remove argument
    idx <- 0 == nchar(act.args)
    if( any(idx) ){
      fn.arg.names <- fn.arg.names[!idx]
      act.args <- act.args[!idx]
    }
    # argument name == "..." => remove name
    idx <- "..." == fn.arg.names
    if( any(idx) ) fn.arg.names[idx] <- ""
    fn.arg.names <- ifelse(fn.arg.names == "", "", paste(fn.arg.names,"="))
    # combine call
    fn.call <- paste( fn.name,
                      "(",
                      paste(paste(fn.arg.names,act.args),collapse=", "),
                      ")" )
    cat("evaluation of:\n  ", fn.call,"\n")
    print(try(eval(parse(text=fn.call),envir=eval.envir)))
  })

  invisible(NULL)
}

