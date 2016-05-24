weaveR<-function(in.file,out.file,show.code=TRUE,show.text=TRUE,
                 replace.umlaute=TRUE,eval_Sexpr=FALSE){
  # german documentation of the code:
  # look for file webR.pdf, P. Wolf 050204, 060517, 070307, 070830, 110919
  # require(tcltk) # 140306
  pat.use.chunk    <- paste("<","<(.*)>",">",sep="")
  pat.use.chunk.line <- paste("(.*)",pat.use.chunk,"(.*)",sep="")
  pat.chunk.header <- paste("^<","<(.*)>",">=",sep="")
  pat.verbatim.begin<-"\\\\begin\\{verbatim\\}"
  pat.verbatim.end<-"\\\\end\\{verbatim\\}"
  pat.leerzeile<-"^(\\ )*$"
  pat.KlAffeGG <- paste("@",">",">",sep="")
  pat.KlAffeKK <- paste("@","<","<",sep="")
  pat.Sp.open <- paste("DoSp","OpenKl-esc",sep="")
  pat.Sp.close <- paste("DoSp","CloseKl-esc",sep="")
  pat.Eck.open <- paste("DoEck","OpenKl-esc",sep="")
  pat.Eck.close <- paste("DoEck","CloseKl-esc",sep="")

  .Tcl("set xyz [encoding system]"); UTF<-tclvalue("xyz")
  UTF<-0<length(grep("utf",UTF))
  if(exists("DEBUG")){
    if(UTF) cat("character set: UTF\n") else cat("character set: not utf\n")
  }
  if(!UTF){ 
    char267<-eval(parse(text='"\\267"'))
  }

  if(!file.exists(in.file)) in.file<-paste(in.file,"rev",sep=".")
  if(!file.exists(in.file)){
    cat(paste("ERROR:",in.file,"not found!!??\n"))
    return("Error in weave: file not found")
  }
  # input<-scan(in.file,what="",sep="\n",blank.lines.skip = FALSE)
  input<-readLines(in.file) # 2.1.0
  try(if(replace.umlaute&&UTF && any(idx<-is.na(iconv(input,"","LATIN1")))){  
      # LATIN1-Dok : input<-iconv(input,"LATIN1","") # old
      input[idx] <- iconv( input[idx], "Latin1", "") # 121120
  })
  length.input<-length(input)

  input<-gsub(pat.KlAffeGG,"DoSpCloseKl-esc",gsub(pat.KlAffeKK,"DoSpOpenKl-esc",input))
  input<-gsub("@\\]\\]","DoEckCloseKl-esc",gsub("@\\[\\[","DoEckOpenKl-esc",input))

  empty.index<-grep(pat.leerzeile,input)
  text.start.index<-which("@"==substring(input,1,1))

  a<-rep(0,length.input)
  a[grep(pat.verbatim.begin,input)]<-1
  a[grep(pat.verbatim.end,input)]<- -1
  a<-cumsum(a)
  verb.index<-which(a>0)

  code.start.index<-grep(pat.chunk.header,input)
  use.index<-grep(pat.use.chunk,input)
  use.index<-use.index[is.na(match(use.index,code.start.index))]

  a<-rep(0,length.input)
  a[text.start.index]<- -1; a[code.start.index]<-2
  a<-cbind(c(text.start.index,code.start.index),
    c(rep(-1,length(text.start.index)),rep(1,length(code.start.index))))
  a<-a[order(a[,1]),,drop=FALSE]
  b<-a[a[,2]!=c(-1,a[-length(a[,1]),2]),,drop=FALSE]
  a<-rep(0,length.input); a[b[,1]]<-b[,2]
  a<-cumsum(a); a[code.start.index]<-0 
  ## a[empty.index]<-0 ?? this was not a good idea 070709
  code.index<-which(a>0)
  code.index<-code.index[is.na(match(code.index,use.index))]

  line.typ<-rep("TEXT" ,length.input)
  line.typ[empty.index]<-"EMPTY"
  line.typ[text.start.index]<-"TEXT-START"
  line.typ[verb.index]<-"VERBATIM"
  line.typ[use.index]<-"USE"
  line.typ[code.start.index]<-"HEADER"
  line.typ[code.index]<-"CODE"

  is.code.line<-text.start.indicator<-rep(0,length.input)
  text.start.indicator[1]<-1; text.start.indicator[text.start.index]<-1
  text.start.indicator<-cumsum(text.start.indicator)
  is.code.line[code.start.index]<-0-text.start.indicator[code.start.index]
  is.code.line<-cummin(is.code.line)
  is.code.line<-(text.start.indicator+is.code.line) < 1
  is.code.line[code.start.index]<-FALSE
  ## TSI<<-text.start.index; CSI<<-code.start.index; UI<<-use.index # ->debugging


  code.chunk.names<-code.start.lines<-sub(pat.chunk.header,"\\1",
                                          input[code.start.index])
  use.lines<-input[use.index]
  code.lines<-input[code.index]

  if(eval_Sexpr){ if(!exists("revive.env")) revive.env <- ""
                  # finde alle Textzeilen
                  is.text.line<-line.typ=="TEXT"
                  # extrahiere alle Textzeilen
                  text.lines<-input[is.text.line]
                  # suche Textzeilen mit Sexpr-Expressions
                  sexpr.lines<-grep("\\Sexpr\\{.*\\}",text.lines)
                  # falls es keine Sexpr-Expressions gibt, relax
                  if(0<length(sexpr.lines)){
                    # arbeite Zeilen mit Sexpr-Expressions nacheinander ab
                    for(l in seq(along=sexpr.lines)){
                      # hole Nummer l der Zeilen, die Sexpr-Expressions enthalten 
                      cand<-text.lines[sexpr.lines[l]]
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
                        if(exists("DEBUG")){ print(code) }
                        if(identical(revive.env,"")) 
                           result <- try(eval(parse(text=code),envir=revive.env))  
                        else 
                           result <- try(eval(parse(text=code)))
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
                      text.lines[sexpr.lines[l]]<-paste(cand,collapse="")  
                    }
                    # aktualisiere Textzeilen 
                    input[is.text.line]<-text.lines
                  }

                  # finde alle Textzeilen
                  is.text.line<-line.typ=="TEXT"
                  # extrahiere alle Textzeilen
                  text.lines<-input[is.text.line]
                  # suche Textzeilen mit Splot-Expressions
                  splot.lines<-grep("\\Splot\\{.*\\}",text.lines)
                  # falls es keine Splot-Anweisungen gibt, relax
                  if(0<length(splot.lines)){
                    if(identical(revive.env,"")) revive.sys <- get("revive.sys",revive.env)
                    copy.plot<-function(psname,latexheight,height,width,horizontal,center=TRUE,...){
                      psname<-as.character(substitute(psname))
                      if(exists("revive.sys")){
                        if(missing(latexheight)) latexheight <- get("psheight.sys",revive.sys)
                        if(missing(height))      height <-      get("psdesignheight.sys",revive.sys)
                        if(missing(width))       width <-       get("psdesignwidth.sys",revive.sys)
                        if(missing(horizontal))  horizontal <-  get("pshorizontal.sys",revive.sys)
                      } else {
                        if(missing(latexheight))     latexheight <- "10cm"
                        if(missing(height))          height <-      6
                        if(missing(width))           width <-       6
                        if(missing(horizontal))      horizontal <-  FALSE
                      }
                      if(exists("DEBUG")){ cat(psname,latexheight,height,width,horizontal) }
                      if(0==length(grep("\\.ps$",psname))) psname<-paste(psname,".ps",sep="") 
                      news<-paste(if(center)"\n\\begin{center}","\\includegraphics[",
                                  "height=",latexheight,"]{",psname,"}",
                                  if(center)"\\end{center}\n",sep="") 
                      result<-try({dev.copy(postscript,psname,horizontal=horizontal,
                                  width=width,height=height,...);dev.off()}
                      )
                      if(class(result)=="try-error"){ news <- result } 
                      news
                    }
                    # arbeite Zeilen mit Splot-Expressions nacheinander ab
                    for(l in seq(along=splot.lines)){
                      # hole Nummer l der Zeilen, die Splot-Expressions enthalten 
                      cand<-text.lines[splot.lines[l]]
                      # knacke Kandidaten-Zeile an der Stelle auf, an der \Splot gefunden wird
                      cand<-unlist(strsplit(cand,"\\\\Splot"))
                      # cand[1] ist der vor der ersten Expression, 
                      # cand[i+1] der mit der i-ten Expression beginnt
                      # alle Expressions der Zeile werden nacheinander abgearbeitet
                      for(j in seq(cand)[-1]){
                        # ncandj zeigt die Laenge von Kandidat j an
                        ncandj<-nchar(cand[j])
                        # splot verwaltet den j-ten Kandidaten zeichenweise
                        splot<-substring(cand[j],1:ncandj,1:ncandj) 
                        # es gilt die beendende Klammer von Splot zu finden
                        brack<-cumsum((splot=="{")-(splot=="}")) 
                        # n.splot zeigt die Stelle der schliessenden-Klammer
                        n.splot<-which(brack==0)[1]; if(is.na(n.splot)) next
                        # mit n.splot greifen wir den vorderen Teil von splot und evaluieren
                  ###
                        code <- paste("copy.plot(",paste(collapse="",splot[2:(n.splot-1)]),")")
                        if(exists("DEBUG")){ print(code) }
                  ###
                        result <- try(eval(parse(text=code)))      
                        # wenn nichts rauskommt, ist nichts zu modifizieren
                        if(!identical(result,"")) { 
                          # print("---");print(result);print("---")
                          # im Fehlerfall muss es eine Meldung geben
                          if(class(result)=="try-error"){ 
                            result<-paste("[[\\Splot-error:",
                                          paste(splot[1:n.splot],collapse=""),"]]",collaspe="")
                          }
                        }
                        # das Ergebnis des j-ten Ausdrucks wird vorn, 
                        # also wo das Kommando stand eingetragen
                        cand[j]<-paste(result, substring(cand[j],n.splot+1),sep="")
                      }
                      text.lines[splot.lines[l]]<-paste(cand,collapse="")  
                    }
                    # aktualisiere Textzeilen 
                    input[is.text.line]<-text.lines
                  }
 }
  pos.obj.idx<-grep("^@index.of.objects",input)

  pos.chunk.list<-grep("^@list.of.chunks",input)
    
  ref.infos.found<-FALSE
  # extract lines containing calls of other code chunks
  lines.use<-which(line.typ=="USE"&is.code.line)
  include.use.infos<-0==length(grep("^@no.used.in.infos",input))
  if(include.use.infos&&length(lines.use)>0){
   # find header lines
   names.header<-input[code.start.index]
   # extract set header names: remove "<", ">" and characters not belonging to the name
   names.header.uniq<-sub(paste(pat.chunk.header,".*",sep=""),"\\1",
                          unique(sort(names.header)))

   # lines with uses of code chunks
   names.use.cand<-input[lines.use]; l.u<-n.u<-NULL
   for(ii in seq(along=lines.use)){ # aufknacken mehrerer Chunks in einer Code-Zeile
      h<-names.use.cand[ii] 
      repeat{  # old: last<-sub(paste("^.*<","<(.*)>",">.*",sep=""),"\\1",h) 
        # if(!exists("max.wd")) max.wd<-10; max.wd<-max.wd-1; if(max.wd<1) break
        last<-sub(pat.use.chunk.line,"\\2",h) # extract name
        if(last!=h){ # something found during substitution => chunk use found
          l.u<-c(l.u,lines.use[ii]); n.u<-c(n.u, last)
          # old: h<-sub(paste("^(.*)<","<.*>",">.*",sep=""),"\\1",h) # rm 
          h<-sub(pat.use.chunk.line,"\\1",h) # rm identified chunk use 
          if(nchar(h)==0) break
        } else break # no more chunk uses in line ii
      }
   }
   names.use<-n.u; lines.use<-l.u 

   # chunk uses found: names (names.use) and lines (lines.use)
   # find headers that have been used, their lines and compute used-in-info
   # remove brackets... old: sub(paste("^.*<","<(.*)>",">.*",sep=""),"\\1",names.header)
   names.header<-sub(pat.use.chunk.line,"\\2",names.header) 
   ind<-!is.na(match(names.header,names.use))
   names.header.used<-names.header[ind]; lines.header.used<-code.start.index[ind] # 
   if((anz<-length(names.header.used))>0){ 
     ref.infos.found<-TRUE; used.in.message<-rep("",anz)
     lines.used.in.message<-rep(0,anz)
     # find for each header of names.header.use the numbers of section of their uses
     for(i in 1:anz){
          idx.found<-which(names.header.used[i]==names.use)
          l<-lines.use[idx.found]
          # find number of chunks calling names.header.used[i]
          used.in.no<-unique(unlist(lapply(l,function(x) sum(code.start.index<x)))) #
          # construct message and save line number of input that has to be changed
          used.in.message[i]<-paste("{\\quad$\\subset$ ",
                                    paste(used.in.no,collapse=", "),"}")
          lines.used.in.message[i]<-lines.header.used[i]
     }

 
   }
  }
  obj.set<-obj.index<-NULL

  ind<-0<length(pos.obj.idx) && 0==length(grep("^@index.reduced",input))
  if(ind) {
    # Kandidatensuche
    a<-unlist(strsplit(input[code.index],";"))
    a<-sub("<-.*","",a[grep("<-",a)]) # Zuweisungen suchen
    a<-gsub(" ","",a[!is.na(a)])      # Leerzeichen entfernen
    # Indizes und Unternamen entfernen:
    a<-sub("\\[.*","",a); a<-sub("\\$.*","",a) 
    a<-a[grep("^[a-zA-Z.]",a)]        # Beginn mit Nicht-Ziffer
    a<-a[grep("^[a-zA-Z.0-9_]*$",a)]  # erlaubte Zeichen
    a<-a[nchar(a)>1]                  # nur echte Strings merken
    obj.set<-sort(unique(a))          # Zusammenfassung
  }
  # explizite angegebene Namen
  ind<-0<length(pos.obj.idx) && 0<length( a<-grep("^@index[^o]",input,value=TRUE) )
  if(ind){
    a<-sub("^@index *","",a); a<-gsub(","," ",a); 
    a<-unlist(strsplit(gsub(" +"," ",a)," "))
    obj.set<-sort(unique(c(a,obj.set))) # set of object names
  }
  if(length(obj.set)>0){
    if(0<(anz<-length(code.start.index))){
      obj.used.in<-matrix(0,2,0) 
      for(no in 1:anz){
        # extract code chunk no
        c.start<-code.start.index[no]+1
        c.end<-text.start.index[which(c.start<text.start.index)[1]]-1
        if(is.na(c.end)) c.end<-length(input)  
        if(c.end<c.start) next
        a<-paste("",input[c.start:c.end],"") # code von code chunk no
        # check occurance of all candidats
        h<-sapply(obj.set,function(x) 
           0<length(grep(paste("[^a-zA-Z.0-9]",x,"[^a-zA-Z.0-9]",sep=""),a)))
        if(any(h)) obj.used.in<-cbind(obj.used.in,rbind(which(h),no))
      }
      # obj.used[2,] shows chunk numbers, obj.used[1,] that candidates
      a<-split(obj.used.in[2,],obj.used.in[1,])
      if(length(a)==1&&length(a[[1]])==1) names(a[[1]])<-colnames(obj.used.in)[1]
      # list element i stores the numbers of chunks where object i has been found
      a<-lapply(a,function(x){ 
             x<-paste(names(x)[1],"\\quad$\\in$",paste(x,collapse=", ")) 
                  })
      obj.index<-paste(unlist(a),collapse="\\\\\n"); names(obj.index)<-NULL
    } 

  }

  ## input[text.start.index]<-""
  input[text.start.index]<-paste(
                        "\\ifodd\\value{IsInCodeChunk}",
                        "\\setcounter{IsInCodeChunk}{0}",
                        "\\vspace{-\\parskip}\\par\\hspace*{-\\parindent}",
                        "\\textchunkcommands\\fi",
                        sep="")

  # find section numbers
  no<-seq(along=code.start.index)
  # find first occurences of code chunks
  def.ref.no<-match(gsub("\\ ","",code.start.lines),
                    gsub("\\ ","",code.start.lines))
  # construct modified header lines
  code.start.lines<-paste("\\makemarginno ", 
                          "$\\langle${\\it ",code.start.lines,"}\\ $",def.ref.no,
                          "\\rangle",ifelse(no!=def.ref.no,"+",""),"\\equiv$",sep="")
  # save modified header lines
  input[code.start.index]<-code.start.lines

  # get use lines 
  use.lines<-input[use.index]; is.use.lines.within.code<-is.code.line[use.index]
  # remove and save leeding blanks
  leerzeichen.vor.use<-paste("\\verb|",
                             sub("[^ ](.*)$"," ",use.lines),
                             "|",sep="") ## plus 1 Leerzeichen 
  use.lines<-substring(use.lines,nchar(leerzeichen.vor.use)-7) 
  # loop along use lines
  for(i in seq(along=use.lines)){
    # get single line
    uli<-use.lines[i]
    # split chunk names and other strings
    repeat{
      if(0==length(cand<-grep(pat.use.chunk,uli))) break
      uli.h<-gsub(pat.use.chunk.line,"\\1bReAkuSeChUnK\\2bReAk\\3",uli)
      uli<-unlist(strsplit(uli.h,"bReAk"))
    }
    # find chunk names
    cand<-grep("uSeChUnK",uli); uli<-sub("uSeChUnK","",uli)
    # find chunk numbers of (first) definition
    ref.no<-match(uli[cand],code.chunk.names)
    # include number of definition chunk
    uli[cand]<-paste("$\\langle${\\it ",uli[cand],"} ",ref.no,"$\\rangle$",sep="")
    #  formating code within use references, in code chunk a little different
    if(length(uli)!=length(cand)){
      if(is.use.lines.within.code[i]){ 
        # within code chunks: code (but no the chunk names) has to be escaped
         uli[-cand]<-gsub("\\?","?\\\\verb!?!\\\\verb?",uli[-cand]) #110518
         uli[-cand]<-paste("\\verb?",uli[-cand],"?",sep="") #110518
      }
    }

    use.lines[i]<-paste(uli,collapse="")
  }
  # store modified use lines
  input[use.index]<-ifelse(is.use.lines.within.code,
        paste("\\rule{0mm}{0mm}\\newline",leerzeichen.vor.use,use.lines,"%",sep=""), 
        paste(gsub("[^ ]","",leerzeichen.vor.use),use.lines,sep=""))

  #110518:
  code.lines<-gsub("\\?","?\\\\verb!?!\\\\verb?",code.lines) #110518
  input[code.index]<-paste("\\rule{0mm}{0mm}\\newline\\verb?",
                             " ",code.lines," ","?%",sep="")    
  typ<-"TEXT"
  index<-which(line.typ==typ)
  code.im.text.index<-index[grep("\\[\\[(.*)\\]\\]",input[index])]

  if(0<length(code.im.text.index)){
    lines.to.check<-input[code.im.text.index]
    lines.to.check<-strsplit(lines.to.check," ") # Zerlegung in Worte
    lines.to.check<-unlist(lapply(lines.to.check,function(x){
        ind.cand<-grep("^\\[\\[(.*)\\]\\]$",x)
        if(0<length(ind.cand)){
          cand<-gsub("^\\[\\[(.*)\\]\\]$","\\1",x[ind.cand])
          cand<-gsub("\\[\\[","DoEckOpenKl-esc",cand)
          cand<-gsub("\\]\\]","DoEckCloseKl-esc",cand)
          cand<-gsub("\\\\","BaCkSlAsH",cand)
          cand<-gsub("([#$&_%{}])","\\\\\\1",cand) #2.1.0
          cand<-gsub("BaCkSlAsH","{\\\\char'134}",cand)
          cand<-gsub("\\~","{\\\\char'176}",cand)
          cand<-gsub("\\^","{\\\\char'136}",cand)
          cand<-gsub(pat.Sp.open,"\\\\verb|<<|",cand) # 050612
          cand<-gsub(pat.Sp.close,"\\\\verb|>>|",cand) # 050612
          x[ind.cand]<-paste("{\\tt ",cand,"}",sep="")
        }
        x<-paste(x,collapse=" ")}
    )) # end of unlist(apply(..))

    ind.cand<-grep("\\[\\[(.*)\\]\\]",lines.to.check)
    if(0<length(ind.cand)) {
      # zerlege Zeile in token der Form [[,  ]] und sonstige
      zsplit<-lapply(strsplit(lines.to.check[ind.cand],"\\[\\["),function(x){
         zs<-strsplit(rbind("[[",paste(x[],"aAzsplitAa",sep=""))[-1],"\\]\\]")
         zs<-unlist(lapply(zs,function(y){ res<-rbind("]]",y[])[-1]; res }))
         gsub("aAzsplitAa","",zs)
      })
      # suche von vorn beginnend zusammenpassende [[-]]-Paare
      z<-unlist(lapply(zsplit,function(x){
        repeat{
          cand.sum<-cumsum((x=="[[")-(x=="]]"))
          if(is.na(br.open<-which(cand.sum==1)[1])) break
          br.close<-which(cand.sum==0)
          if(is.na(br.close<-br.close[br.open<br.close][1])) break
          if((br.open+1)<=(br.close-1)){
            h<-x[(br.open+1):(br.close-1)]; 
            h<-gsub("\\\\","BaCkSlAsH",h);
            h<-gsub("([#$&_%{}])","\\\\\\1",h)
            h<-gsub("BaCkSlAsH","{\\\\char'134}",h);
            h<-gsub("\\~","{\\\\char'176}",h) #2.1.0
            h<-gsub(" ","\\\\ ",h) # Leerzeichen nicht vergessen! 060116
            h<-gsub(pat.Sp.open,"\\\\verb|<<|",h) # 050612
            h<-gsub(pat.Sp.close,"\\\\verb|>>|",h) # 050612
          x[(br.open+1):(br.close-1)]<-gsub("\\^","{\\\\char'136}",h)
          }
          x[br.open]<-"{\\tt "; x[br.close]<-"}"
          x<-c(paste(x[1:br.close],collapse=""), x[-(1:br.close)])
        }
        paste(x,collapse="")
      }))
      lines.to.check[ind.cand]<-z
    }

    input[code.im.text.index]<-lines.to.check
  }

  typ<-"HEADER"
  index<-which(line.typ==typ)
  code.im.text.index<-index[grep("\\[\\[(.*)\\]\\]",input[index])]

  if(0<length(code.im.text.index)){
    lines.to.check<-input[code.im.text.index]
    lines.to.check<-strsplit(lines.to.check," ") # Zerlegung in Worte
    lines.to.check<-unlist(lapply(lines.to.check,function(x){
        ind.cand<-grep("^\\[\\[(.*)\\]\\]$",x)
        if(0<length(ind.cand)){
          cand<-gsub("^\\[\\[(.*)\\]\\]$","\\1",x[ind.cand])
          cand<-gsub("\\[\\[","DoEckOpenKl-esc",cand)
          cand<-gsub("\\]\\]","DoEckCloseKl-esc",cand)
          cand<-gsub("\\\\","BaCkSlAsH",cand)
          cand<-gsub("([#$&_%{}])","\\\\\\1",cand) #2.1.0
          cand<-gsub("BaCkSlAsH","{\\\\char'134}",cand)
          cand<-gsub("\\~","{\\\\char'176}",cand)
          cand<-gsub("\\^","{\\\\char'136}",cand)
          cand<-gsub(pat.Sp.open,"\\\\verb|<<|",cand) # 050612
          cand<-gsub(pat.Sp.close,"\\\\verb|>>|",cand) # 050612
          x[ind.cand]<-paste("{\\tt ",cand,"}",sep="")
        }
        x<-paste(x,collapse=" ")}
    )) # end of unlist(apply(..))

    ind.cand<-grep("\\[\\[(.*)\\]\\]",lines.to.check)
    if(0<length(ind.cand)) {
      # zerlege Zeile in token der Form [[,  ]] und sonstige
      zsplit<-lapply(strsplit(lines.to.check[ind.cand],"\\[\\["),function(x){
         zs<-strsplit(rbind("[[",paste(x[],"aAzsplitAa",sep=""))[-1],"\\]\\]")
         zs<-unlist(lapply(zs,function(y){ res<-rbind("]]",y[])[-1]; res }))
         gsub("aAzsplitAa","",zs)
      })
      # suche von vorn beginnend zusammenpassende [[-]]-Paare
      z<-unlist(lapply(zsplit,function(x){
        repeat{
          cand.sum<-cumsum((x=="[[")-(x=="]]"))
          if(is.na(br.open<-which(cand.sum==1)[1])) break
          br.close<-which(cand.sum==0)
          if(is.na(br.close<-br.close[br.open<br.close][1])) break
          if((br.open+1)<=(br.close-1)){
            h<-x[(br.open+1):(br.close-1)]; 
            h<-gsub("\\\\","BaCkSlAsH",h);
            h<-gsub("([#$&_%{}])","\\\\\\1",h)
            h<-gsub("BaCkSlAsH","{\\\\char'134}",h);
            h<-gsub("\\~","{\\\\char'176}",h) #2.1.0
            h<-gsub(" ","\\\\ ",h) # Leerzeichen nicht vergessen! 060116
            h<-gsub(pat.Sp.open,"\\\\verb|<<|",h) # 050612
            h<-gsub(pat.Sp.close,"\\\\verb|>>|",h) # 050612
          x[(br.open+1):(br.close-1)]<-gsub("\\^","{\\\\char'136}",h)
          }
          x[br.open]<-"{\\tt "; x[br.close]<-"}"
          x<-c(paste(x[1:br.close],collapse=""), x[-(1:br.close)])
        }
        paste(x,collapse="")
      }))
      lines.to.check[ind.cand]<-z
    }

    input[code.im.text.index]<-lines.to.check
  }

  typ<-"USE"
  index<-which(line.typ==typ)
  code.im.text.index<-index[grep("\\[\\[(.*)\\]\\]",input[index])]

  if(0<length(code.im.text.index)){
    lines.to.check<-input[code.im.text.index]
    lines.to.check<-strsplit(lines.to.check," ") # Zerlegung in Worte
    lines.to.check<-unlist(lapply(lines.to.check,function(x){
        ind.cand<-grep("^\\[\\[(.*)\\]\\]$",x)
        if(0<length(ind.cand)){
          cand<-gsub("^\\[\\[(.*)\\]\\]$","\\1",x[ind.cand])
          cand<-gsub("\\[\\[","DoEckOpenKl-esc",cand)
          cand<-gsub("\\]\\]","DoEckCloseKl-esc",cand)
          cand<-gsub("\\\\","BaCkSlAsH",cand)
          cand<-gsub("([#$&_%{}])","\\\\\\1",cand) #2.1.0
          cand<-gsub("BaCkSlAsH","{\\\\char'134}",cand)
          cand<-gsub("\\~","{\\\\char'176}",cand)
          cand<-gsub("\\^","{\\\\char'136}",cand)
          cand<-gsub(pat.Sp.open,"\\\\verb|<<|",cand) # 050612
          cand<-gsub(pat.Sp.close,"\\\\verb|>>|",cand) # 050612
          x[ind.cand]<-paste("{\\tt ",cand,"}",sep="")
        }
        x<-paste(x,collapse=" ")}
    )) # end of unlist(apply(..))

    ind.cand<-grep("\\[\\[(.*)\\]\\]",lines.to.check)
    if(0<length(ind.cand)) {
      # zerlege Zeile in token der Form [[,  ]] und sonstige
      zsplit<-lapply(strsplit(lines.to.check[ind.cand],"\\[\\["),function(x){
         zs<-strsplit(rbind("[[",paste(x[],"aAzsplitAa",sep=""))[-1],"\\]\\]")
         zs<-unlist(lapply(zs,function(y){ res<-rbind("]]",y[])[-1]; res }))
         gsub("aAzsplitAa","",zs)
      })
      # suche von vorn beginnend zusammenpassende [[-]]-Paare
      z<-unlist(lapply(zsplit,function(x){
        repeat{
          cand.sum<-cumsum((x=="[[")-(x=="]]"))
          if(is.na(br.open<-which(cand.sum==1)[1])) break
          br.close<-which(cand.sum==0)
          if(is.na(br.close<-br.close[br.open<br.close][1])) break
          if((br.open+1)<=(br.close-1)){
            h<-x[(br.open+1):(br.close-1)]; 
            h<-gsub("\\\\","BaCkSlAsH",h);
            h<-gsub("([#$&_%{}])","\\\\\\1",h)
            h<-gsub("BaCkSlAsH","{\\\\char'134}",h);
            h<-gsub("\\~","{\\\\char'176}",h) #2.1.0
            h<-gsub(" ","\\\\ ",h) # Leerzeichen nicht vergessen! 060116
            h<-gsub(pat.Sp.open,"\\\\verb|<<|",h) # 050612
            h<-gsub(pat.Sp.close,"\\\\verb|>>|",h) # 050612
          x[(br.open+1):(br.close-1)]<-gsub("\\^","{\\\\char'136}",h)
          }
          x[br.open]<-"{\\tt "; x[br.close]<-"}"
          x<-c(paste(x[1:br.close],collapse=""), x[-(1:br.close)])
        }
        paste(x,collapse="")
      }))
      lines.to.check[ind.cand]<-z
    }

    input[code.im.text.index]<-lines.to.check
  }


  # include Referenzinformationen
  if(include.use.infos&&ref.infos.found)
     input[lines.used.in.message]<-paste(input[lines.used.in.message],
                                         used.in.message)

  # merke Header fuer Chunk Index
  chunk.index<-NULL 
  ind<-0<length(pos.chunk.list)
  if(ind){ 
    # Randnummer entfernen
    chunk.index<-sub("\\\\makemarginno.","",input[code.start.index]) 
    # + Zeichen entfernen
    chunk.index<-sub("rangle[+]*.equiv","rangle",chunk.index)    
    first<-match(chunk.index,chunk.index)       # jeweils erste Chunks finden 
    no.ext<-split(seq(along=chunk.index),first) # Nummern gleicher Chunks suchen
    no.ext<-unlist(lapply(no.ext,function(x){
           if(length(x)==1) " " else paste("\\cup",paste(x[-1],collapse="\\cup"))
                                 } ))   # Erweiterungsnummern als String ablegen
    chunk.index<-unique(chunk.index); first<-unique(first)   # gleiche entfernen
    if(0<length(first)){
      # Erweiterungs-Infos einbauen
      chunk.index<-paste(sub(".rangle.*","",chunk.index),no.ext,
                         sub(".*rangle","\\\\rangle",chunk.index),sep="") 
      pageref<-function(no){
        find.label<-function(no){
            label<-paste("CodeChunkLabel",
                paste(LETTERS[1+as.numeric(substring((no+10000),2:5,2:5))],
                              collapse=""),sep="")
        }
          
        label<-sapply(no,find.label)
        paste("p\\pageref{",label,"}",sep="")
      }
      chunk.index<-paste(chunk.index,"\\dotfill",pageref(first))

      chunk.index<-sort(chunk.index)              # sortieren
    }
  }

  if(0<length(chunk.index)){
    chunk.index<-paste(chunk.index,collapse="\\\\") # newline ?
    chunk.index<-paste("{\\paragraph{Code Chunk Index}\\small",
                       "\\rule{0mm}{0mm}\\\\[1.5ex]",
                       chunk.index,"\\\\\\rule{0mm}{0mm}}")
  }

  obj.index<-obj.index[obj.index!=""]
  if(0<length(obj.index)){ 
    obj.index<-gsub("_","\\\\_",obj.index)
    obj.index<-paste(obj.index,collapse="\\\\") # newline ?
    obj.index<-paste(
          "{\\paragraph{Object Index}\\footnotesize\\rule{0mm}{0mm}\\\\[1.5ex]",
          obj.index,"\\\\\\rule{0mm}{0mm}\\par}")
  }

  if(0<length(pos.obj.idx)   &&length(obj.index)  >0) input[pos.obj.idx]<-obj.index
  if(0<length(pos.chunk.list)&&length(chunk.index)>0) input[pos.chunk.list]<-chunk.index

  label<-function(no){
    find.label<-function(no){
        label<-paste("CodeChunkLabel",
            paste(LETTERS[1+as.numeric(substring((no+10000),2:5,2:5))],
                          collapse=""),sep="")
    }
      
    label<-sapply(no,find.label)
    paste("\\label{",label,"}",sep="")
  }
  input[code.start.index]<-paste(input[code.start.index],
                                 label(seq(along=code.start.index)))


  if(replace.umlaute){
   if(!UTF){
   # im Tcl/Tk-Textfenster eingegeben -> iso-8859-1 
   # (see: $ man iso-8859-1 / Latin1 / unicode)
      pc<-eval(parse(text='"\\303"'))  # UTF-8-pre-char
      uml.utf.8 <-eval(parse(text='"\\244\\266\\274\\204\\226\\234\\237"'))
      uml.latin1<-eval(parse(text='"\\344\\366\\374\\304\\326\\334\\337"'))
      input<-chartr(uml.utf.8,uml.latin1,gsub(pc,"",input)) # utfToLatin1
      input<-gsub(substring(uml.latin1,7,7),paste("{\\\\s","s}",sep=""),input) # sz
      uml.pattern<-eval(parse(text='"(\\344|\\366|\\374|\\304|\\326|\\334)"'))
      input<-gsub(uml.pattern,"\\\\\"\\1",input)  # replace Umlaute ae->&aeuml; 
  # replace Umlaute &aeuml;->&auml;
      input<-chartr(substring(uml.latin1,1,6),"aouAOU",input)   
   }else{
    input<-gsub("\303\237",paste("{\\\\s","s}",sep=""),input)
    input<-gsub("(\303\244|\303\266|\303\274|\303\204|\303\226|\303\234)",
                              "\\\\\"\\1",input)
    input<-chartr("\303\244\303\266\303\274\303\204\303\226\303\234", 
                                "aouAOU", input)
   }
   if(exists("DEBUG")){
    cat("german Umlaute replaced\n")
   }
  }
  input<-gsub(pat.Sp.close,">>",gsub(pat.Sp.open,"<<",input))
  input<-gsub(pat.Eck.close,"]]",gsub(pat.Eck.open,"[[",input))

  input[1]<-paste(
       "\\newcounter{Rchunkno}",
       "\\newcounter{IsInCodeChunk}\\setcounter{IsInCodeChunk}{1}",
       "\\newcommand{\\codechunkcommands}{\\relax}",
       "\\newcommand{\\textchunkcommands}{\\relax}",
       "\\newcommand{\\Routputcommands}{\\relax}",
       "\\newcommand{\\makemarginno}",
            "{\\par\\vspace{-0.5\\parskip}\\codechunkcommands",
            "\\stepcounter{Rchunkno}",
            "\\setcounter{IsInCodeChunk}{1}",
            "\\noindent\\hspace*{-3em}",
            "\\makebox[0mm]{\\arabic{Rchunkno}}\\hspace*{3em}}",
       input[1],sep="")

  input<-sub("^\\\\begin\\{verbatim\\}",
             "\\\\par\\\\Routputcommands\\\\begin{verbatim}",input)
  input<-sub("^\\\\end\\{verbatim\\}",
             "\\\\end{verbatim}\\\\textchunkcommands",input)

  if(show.code == "echo"){
    # checke code chunk names
    idx <- grep("echo *= *FALSE",input[code.start.index])
    if(0 < length(idx)){
      remove.chunk.set <- code.start.index[idx]
      for(i in seq(along=remove.chunk.set)){
        # remove code chunks with names containing "echo=FALSE"
        st <- en <- remove.chunk.set
        repeat{
          en <- en + 1; if( length.input < en || line.typ[en] == "TEXT-START" ) break
        }
        en <- en - 1  
        input[st:en] <- " "
      }
    }
  }
  if(show.code==FALSE){
     input[code.index] <-"."
     input[use.index] <-":"
     an<-grep("\\\\begin(.*)\\{document\\}",input)[1]
     if(length(tit<-grep("\\\\maketitle",input))>0) an<-tit
     input[an]<-paste(input[an],
                      "${}^*$ --- only the TEXT of the paper ---\\par")
  }
  if(show.text==FALSE){
     input<-sub("^%.*","%",input)
     an<-grep("\\\\begin(.*)\\{document\\}",input)[1]
     en<-grep("\\\\end(.*)\\{document\\}",input)[1]
     text.index<-which(line.typ=="TEXT")
     text.index<-text.index[an<text.index&text.index<en]
     input[c(text.index, verb.index)] <-"."
     if(length(tit<-grep("\\\\maketitle",input))>0) an<-tit
     input[an]<-paste(input[an],
                      "${}^*$ --- only the CODE of the paper ---\\par")
  }

  if(replace.umlaute && 0<length(ind<-grep(".newline.verb",input))){
     ind2<-grep("langle(.*)rangle",input[ind])
     if(0<length(ind2)) ind<-ind[-ind2]
     if(0<length(ind)){
      inp<-input[ind];  
      if(!UTF){
      # im Tcl/Tk-Textfenster eingegeben -> iso-8859-1 
      # (see: $ man iso-8859-1 / Latin1 / unicode)
        # \"a -> ae, ... oe, ue, Ae, Oe, Ue, ß
        u<-uml.latin1<-unlist(strsplit(eval(parse(
             text='"\\344\\366\\374\\304\\326\\334\\337"')),""))
        inp<-gsub(paste('\\\\"',"a",sep=""),u[1],inp)
        inp<-gsub(paste('\\\\"',"o",sep=""),u[2],inp) 
        inp<-gsub(paste('\\\\"',"u",sep=""),u[3],inp)
        inp<-gsub(paste('\\\\"',"A",sep=""),u[4],inp)
        inp<-gsub(paste('\\\\"',"O",sep=""),u[5],inp)
        inp<-gsub(paste('\\\\"',"U",sep=""),u[6],inp) ##{
        inp<-gsub(paste(".\\\\s","s}",sep=""),u[7],inp)
      }else{
        # pc<-eval(parse(text='"\\303"'))  # UTF-8-pre-char
        uml.utf.8 <-eval(parse(text=
             '"\\303\\244\\303\\266\\303\\274\\303\\204\\303\\226\\303\\234\\303\\237"'))
        u<-substring(uml.utf.8,1:7,1:7)
        inp<-gsub(paste('\\\\"',"a",sep=""),u[1],inp)
        inp<-gsub(paste('\\\\"',"o",sep=""),u[2],inp)
        inp<-gsub(paste('\\\\"',"u",sep=""),u[3],inp)
        inp<-gsub(paste('\\\\"',"A",sep=""),u[4],inp)
        inp<-gsub(paste('\\\\"',"O",sep=""),u[5],inp)
        inp<-gsub(paste('\\\\"',"U",sep=""),u[6],inp) ##{
        inp<-gsub(paste(".\\\\s","s}",sep=""),u[7],inp)
      }
      input[ind]<-inp
     }
  }
  if(exists("DEBUG")){
    cat("german Umlaute in code lines inserted\n")
  }


  if(missing(out.file)||in.file==out.file){
    out.file<-sub("\\.([A-Za-z])*$","",in.file)
  }
  if(0==length(grep("\\.tex$",out.file)))
    out.file<-paste(out.file,".tex",sep="")
  base::cat(input,sep="\n",file=out.file)
  base::cat("weave process finished\n")

}

