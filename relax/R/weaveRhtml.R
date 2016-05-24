weaveRhtml<-function(in.file,out.file,replace.umlaute=TRUE){
  # german documentation of the code:
  # look for file webR.pdf, P. Wolf 060920 / 070309 / 070830 / 071016
  # require(tcltk) # 140306
  pat.use.chunk<-paste("<","<(.*)>",">",sep="")
  pat.use.chunk.line <- paste("(.*)",pat.use.chunk,"(.*)",sep="")
  pat.chunk.header<-paste("^<","<(.*)>",">=",sep="")

  pat.verbatim.begin<-"\\\\begin\\{verbatim\\}"
  pat.verbatim.end<-"\\\\end\\{verbatim\\}"
  pat.leerzeile<-"^(\\ )*$"

  pat.KlAffeGG <- paste("@",">",">",sep="")
  pat.KlAffeKK <- paste("@","<","<",sep="")
  pat.Sp.open <- paste("DoSp","OpenKl-ESC",sep="")
  pat.Sp.close <- paste("DoSp","CloseKl-ESC",sep="")
  pat.Eck.open <- paste("DoEck","OpenKl-ESC",sep="")
  pat.Eck.close <- paste("DoEck","CloseKl-ESC",sep="")

  .Tcl("set xyz [encoding system]"); UTF<-tclvalue("xyz")
  UTF<-0<length(grep("utf",UTF))
  if(exists("DEBUG")){ 
     if(UTF) cat("character set: UTF\n") else cat("character set: not utf\n") 
  }
  get.argument<-function(command,txt,default="",kla="{",kle="}",
      dist=TRUE,not.found.info="no parameter"){
    ## print("get.argument")
    command<-paste("\\\\",command,sep="")
    if(0==length(grep(command,txt))) return(not.found.info)
    txt<-unlist(strsplit(paste(txt,collapse="\n"),command))[-1]
    arg<-lapply(txt,function(x){ 
         n<-nchar(x); if(n<3) return(x)
         x<-substring(x,1:n,1:n)
         h<-which(x==kla)[1]; if(is.na(h)) h<-1
         if(dist)x<-x[h:length(x)]
         k<-which(cumsum((x==kla)-(x==kle))==0)[1]
         ifelse(k<=2,default,paste(x[2:(k-1)],collapse=""))
    })
    arg
  }
  get.head.argument<-function(command,txt,default="",kla="{",kle="}",dist=TRUE){
  ## print("get.head.argument")
    command<-paste("\\\\",command,sep="")
    txt<-unlist(strsplit(paste(txt,collapse="\n"),command))[-1]
    arg<-lapply(txt,function(x){ 
         n<-nchar(x); x<-substring(x,1:n,1:n)
         if(dist)x<-x[which(x==kla)[1]:length(x)]
         k<-which(cumsum((x==kla)-(x==kle))==0)[1]
         paste(x[2:(k-1)],collapse="")
    })
    unlist(arg)
  }
  transform.command<-function(command,txt,atag="<i>",etag="</i>",
                         kla="{",kle="}"){
  ## print("transform.command")
    command<-paste("\\\\",command,sep="")
    ##  if(0==length(grep(command,txt))){print("hallo"); return(txt)}
    txt<-unlist(strsplit(paste(txt,collapse="\n"),command))
    tx<-unlist(lapply(txt[-1],function(x){ 
         n<-nchar(x); if(n<4) return(x)
         x<-substring(x,1:n,1:n)
         an<-which(x==kla)[1]
         en<-which(cumsum((x==kla)-(x==kle))==0)[1]
        if(!is.na(an)&&!is.na(en)) 
         paste(atag,paste(x[(an+1):(en-1)],collapse=""),etag,
              paste(x[-(1:en)],collapse="")) else x
    }))
    unlist(strsplit(c(txt[1],tx),"\n"))
  }
  transform.command.line<-function(command,txt,atag="<i>",etag="</i>",
                         kla="{",kle="}"){
    command<-paste("\\\\",command,sep="")
    if(0==length(ind<-grep(command,txt))){return(txt)}
    txt.lines<-txt[ind]
    txt.lines<-strsplit(txt.lines,command)
    txt.lines<-lapply(txt.lines,function(xxx){ 
      for(i in 2:length(xxx)){
        m<-nchar(xxx[i])
        if(is.na(m)) break
        x.ch<-substring(xxx[i],1:m,1:m); x.info<-rep(0,m)
        x.info<-cumsum((x.ch=="{") - (x.ch=="}"))
        h<-which(x.info==0)[1]
        if(!is.na(h)) {x.ch[1]<-atag; x.ch[h]<- etag }
        xxx[i]<-paste(x.ch,collapse="")
      }
      paste(xxx,collapse="")
    })
    txt[ind]<-unlist(txt.lines)
    txt
  }
  transform.structure.command<-function(command,txt,atag="<i>",etag="</i>",
                         kla="{",kle="}"){
  ## print("transform.structure.command")
    command<-paste("\\\\",command,sep="")
    if(0==length(grep(command,txt))){return(txt)}
    txt<-unlist(strsplit(paste(txt,collapse="\n"),command))
    tx<-unlist(lapply(txt[-1],function(x){ 
         n<-nchar(x); if(n<3) return(x)
         x<-substring(x,1:n,1:n); an<-which(x==kla)[1]
         en<-which(cumsum((x==kla)-(x==kle))==0)[1]
        if(is.na(an)||is.na(en)||(an+1)>(en-1)) x<-paste(x,collapse="") else
           paste(paste(if(an==1)""else x[1:(an-1)],collapse=""),
                          atag,paste(x[(an+1):(en-1)],collapse=""),etag, 
                          paste(if(en==n)""else x[-(1:en)],collapse=""),sep="")
    }))
    unlist(strsplit(c(txt[1],tx),"\n"))
  }

  if(!file.exists(in.file)) in.file<-paste(in.file,"rev",sep=".")
  if(!file.exists(in.file)){
    cat(paste("ERROR:",in.file,"not found!!??\n"))
    return("Error in weaveRhtml: file not found")
  }
  input<-readLines(in.file) 
  try(if(replace.umlaute&&UTF && any(is.na(iconv(input,"","LATIN1")))){  # LATIN1-Dok 
            input<-iconv(input,"LATIN1","")
  })
  input<-gsub("\t","      ",input)
  input<-c(input,"@")
  length.input<-length(input)

  input<-sub('^% (<p><img src="t)',"\\1",input)
  input<-sub('^%(<!--latex-)',"\\1",input)
  h<-grep("^[ ]*%",input)
  if(0<length(h)) input<-input[-h]

  input<-gsub(pat.KlAffeGG,"DoSpCloseKl-ESC",gsub(pat.KlAffeKK,"DoSpOpenKl-ESC",input))
  input<-gsub("@\\]\\]","DoEckCloseKl-ESC",gsub("@\\[\\[","DoEckOpenKl-ESC",input))

  empty.index<-grep(pat.leerzeile,input)
  text.start.index<-which("@"==substring(input,1,1))

  a<-rep(0,length(input))
  an<-grep(pat.verbatim.begin,input)
  if(0<length(an)) {
    a[an]<- 1
    en<-grep(pat.verbatim.end,input); a[en]<- -1
    input[a==1]<-"<pre><FONT COLOR=\"#0000FF\">"
    input[a==-1]<-"</font></pre>"
    a<-cumsum(a)
  }
  verb.index<-which(a>0)
  ##input[verb.index]<-paste(input[verb.index],"<br>") # not used because of pre

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
  a<-cumsum(a); a[code.start.index]<-0; a[empty.index]<-0
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


  input<-sub("\\\\begin\\{center}","<center>",input)
  input<-sub("\\\\end\\{center}","</center>",input)
  input<-sub("\\\\begin\\{quote}","<ul>",input)
  input<-sub("\\\\end\\{quote}","</ul>",input)
  input<-sub("\\\\begin\\{itemize}","<ul>",input)
  input<-sub("\\\\end\\{itemize}","</ul>",input)
  input<-sub("\\\\begin\\{enumerate}","<ul>",input)
  input<-sub("\\\\end\\{enumerate}","</ul>",input)
  input<-sub("\\\\item\\[([^]]*)]","<br>\\1&nbsp;",input)
  input<-sub("\\\\item","</li><li>",input)

  input[text.start.index]<-"<p>"    # vorher: @
  lz<-grep("^[ ]*$",input)
  if(0<length(lz)) input[lz]<-"<br>"

  plz.ind<-grep("\\\\includegraphics",input)
  if(0<length(plz.ind)){
    plz<-input[plz.ind]
    h<-unlist(get.argument("includegraphics",plz))
    #hh<-unlist(get.argument("includegraphics",plz,
    #                     default="height=xxcm",kla="[",kle="]"))
    #l.unit<-match(sub("^.*[0-9]([a-z][a-z]).*$","\\1",hh),c("cm","mm","in"))
    #l.unit<-ifelse(is.na(l.unit),1,l.unit)
    #l.unit<-c(1,.1,2.54)[ifelse(is.na(l.unit),1,l.unit)]
    #hh<-sub("^.*[=]([0-9.]+).*$","\\1",hh)
    #hh<-floor(as.numeric(hh)*35*l.unit)
    #hh<-ifelse(is.na(hh),"",paste(" height=",hh,sep=""))
    #h<-paste("<img SRC=\"",sub(".ps$",".jpg",h),"\"",hh,">",sep="")
    h<-paste("<img SRC=\"",sub("\\.[e]{,1}ps$","",h),".jpg","\">",sep="") #090227
    input[plz.ind]<-h
  }

  code.chunk.names<-code.start.lines<-sub(pat.chunk.header,"\\1",input[code.start.index])
  use.lines<-input[use.index]
  code.lines<-input[code.index]
  ## print(input[code.start.index])

  no<-seq(along=code.start.index)
  def.ref.no<-match(gsub("\\ ","",code.start.lines), gsub("\\ ","",code.start.lines))
  code.start.lines<-paste(
        "<a name=\"codechunk",no,"\"></a>",
        "<a href=\"#codechunk",if(length(no)>0)1+(no%%max(no)),"\">",
        "<br>Chunk:",no," <i>&lt;",code.start.lines,def.ref.no,
        "&gt;",ifelse(no!=def.ref.no,"+",""),"=</i></a>",sep="") 
  input[code.start.index]<-code.start.lines

  use.lines<-input[use.index]; is.use.lines.within.code<-is.code.line[use.index]
  leerzeichen.vor.use<-sub("[^ ](.*)$","",use.lines)
  use.lines<-substring(use.lines,nchar(leerzeichen.vor.use))
  h<-gsub("\\ ","&nbsp;",leerzeichen.vor.use)
  leerzeichen.vor.use<-ifelse(is.use.lines.within.code,h,leerzeichen.vor.use)
  such<-pat.use.chunk.line 
  for(i in seq(along=use.lines)){
    uli<-use.lines[i]
    repeat{
      if(0==length(cand<-grep(pat.use.chunk,uli))) break
      uli.h<-gsub(such,"\\1BrEaKuSeCHUNK\\2BrEaK\\3",uli)
      uli<-unlist(strsplit(uli.h,"BrEaK"))
    }
    cand<-grep("uSeCHUNK",uli); uli<-sub("uSeCHUNK","",uli)
    ref.no<-match(uli[cand],code.chunk.names)
    if(is.use.lines.within.code[i]){
      uli[cand]<-paste("</code>&lt;",uli[cand]," ",ref.no,"&gt;<code>",sep="")
    }else{
      uli[cand]<-paste(" &lt;",uli[cand]," ",ref.no,"&gt; ",sep="") 
    }
   #  formating code within use references weaveRhtml
    if(length(uli)!=length(cand)){
      if(!UTF){ 
        uli[-cand]<-paste("",uli[-cand],"",sep="") #050612
      }else{
        uli[-cand]<-paste("",uli[-cand],"",sep="") #060516
      }
    }
    if(is.use.lines.within.code[i]){
      use.lines[i]<-paste("<code>",paste(uli,collapse=""),"</code>")
    }else{ 
      use.lines[i]<-paste(" ",paste(uli,collapse="")," ")
    }
  }
  input[use.index]<-paste("<br>",leerzeichen.vor.use,use.lines)
  leerzeichen.vor.c<-gsub("\t","      ",code.lines)
  leerzeichen.vor.c<-sub("[^ ](.*)$","",leerzeichen.vor.c)
  leerzeichen.vor.c<-gsub("\\ ","&nbsp;",leerzeichen.vor.c)
  # special case "<letter" has to be handled
  code.lines<-gsub("<(.)","&lt;\\1",code.lines) 
  code.lines<-gsub("\\ ","&nbsp;",code.lines) ## multiple blanks in code lines
  if(!UTF){
    input[code.index]<-paste("<br>",leerzeichen.vor.c,"<code>",code.lines,"</code>")
  }else{
    input[code.index]<-paste("<br>",leerzeichen.vor.c,"<code>",code.lines,"</code>") 
  }

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
          cand<-gsub("\\[\\[","DoEckOpenKl-ESC",cand)
          cand<-gsub("\\]\\]","DoEckCloseKl-ESC",cand)
          cand<-gsub("DoSpOpenKl-ESC","<<",cand) # 050612
          cand<-gsub("DoSpCloseKl-ESC",">>",cand) # 050612
          x[ind.cand]<-paste("<code>",cand,"</code>",sep="")
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
            h<-x[(br.open+1):(br.close-1)]
            h<-gsub(" ","&nbsp;",h) # Leerzeichen nicht vergessen! 060116
            h<-gsub("DoSpOpenKl-ESC","<<",h)
            h<-gsub("DoSpCloseKl-ESC",">>",h)
            x[(br.open+1):(br.close-1)]<-h
          }
          x[br.open]<-"<code>"; x[br.close]<-"</code>"
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
          cand<-gsub("\\[\\[","DoEckOpenKl-ESC",cand)
          cand<-gsub("\\]\\]","DoEckCloseKl-ESC",cand)
          cand<-gsub("DoSpOpenKl-ESC","<<",cand) # 050612
          cand<-gsub("DoSpCloseKl-ESC",">>",cand) # 050612
          x[ind.cand]<-paste("<code>",cand,"</code>",sep="")
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
            h<-x[(br.open+1):(br.close-1)]
            h<-gsub(" ","&nbsp;",h) # Leerzeichen nicht vergessen! 060116
            h<-gsub("DoSpOpenKl-ESC","<<",h)
            h<-gsub("DoSpCloseKl-ESC",">>",h)
            x[(br.open+1):(br.close-1)]<-h
          }
          x[br.open]<-"<code>"; x[br.close]<-"</code>"
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
          cand<-gsub("\\[\\[","DoEckOpenKl-ESC",cand)
          cand<-gsub("\\]\\]","DoEckCloseKl-ESC",cand)
          cand<-gsub("DoSpOpenKl-ESC","<<",cand) # 050612
          cand<-gsub("DoSpCloseKl-ESC",">>",cand) # 050612
          x[ind.cand]<-paste("<code>",cand,"</code>",sep="")
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
            h<-x[(br.open+1):(br.close-1)]
            h<-gsub(" ","&nbsp;",h) # Leerzeichen nicht vergessen! 060116
            h<-gsub("DoSpOpenKl-ESC","<<",h)
            h<-gsub("DoSpCloseKl-ESC",">>",h)
            x[(br.open+1):(br.close-1)]<-h
          }
          x[br.open]<-"<code>"; x[br.close]<-"</code>"
          x<-c(paste(x[1:br.close],collapse=""), x[-(1:br.close)])
        }
        paste(x,collapse="")
      }))
      lines.to.check[ind.cand]<-z
    }

    input[code.im.text.index]<-lines.to.check
  }




  if(replace.umlaute){
   if(!UTF){
    # im Tcl/Tk-Textfenster eingegeben -> iso-8859-1 (man iso-8859-1 / Latin1 / unicode
      pc<-eval(parse(text='"\\303"'))  # UTF-8-pre-char
      uml.utf.8 <-eval(parse(text='"\\244\\266\\274\\204\\226\\234\\237"'))
      uml.latin1<-eval(parse(text='"\\344\\366\\374\\304\\326\\334\\337"'))
      input<-chartr(uml.utf.8,uml.latin1,gsub(pc,"",input)) # utfToLatin1
      input<-gsub(substring(uml.latin1,7,7),"&szlig;",input) # replace sz
      uml.pattern<-eval(parse(text='"(\\344|\\366|\\374|\\304|\\326|\\334)"'))
      input<-gsub(uml.pattern,"&\\1uml;",input)  # replace Umlaute ae->&aeuml; 
      # replace Umlaute &aeuml;->&auml;  
      input<-chartr(substring(uml.latin1,1,6),"aouAOU",input)  
   }else{
      input<-gsub("\303\237","&szlig;",input)
      input<-gsub("(\303\244|\303\266|\303\274|\303\204|\303\226|\303\234)",
                                "&\\1uml;",input)
      input<-chartr("\303\244\303\266\303\274\303\204\303\226\303\234", 
                                "aouAOU", input)
   }
  }
  if(exists("DEBUG")) cat("german Umlaute replaced\n")

  #input<-gsub("DoSpClose Kl-esc",">>",gsub("DoSpOpen Kl-esc","<<",input))
  input<-gsub(pat.Sp.close,"&gt;&gt;",gsub(pat.Sp.open,"&lt;&lt;",input))
  input<-gsub(pat.Eck.close,"]]",gsub(pat.Eck.open,"[[",input))

  # find sections, subsections, subsubsections, paragraphs
  atag<-"<h2>"; etag<-"</h2>"; command<-"section"
  command.n<-nchar(command)+2; command.links<-NULL  
  kla<-"{"; kle<-"}"
  ## print("STRUKTUR")
  if(0<length(com.lines<-grep(paste("^\\\\",command,sep=""),input))){
    sec<-NULL
    for(i in seq(along=com.lines)){
      txt<-input[com.lines[i]+0:2]
      txt<-paste(txt,collapse="\n"); n<-nchar(txt)   
      x<-sub("^ *","",substring(txt,command.n))
      if((n.x<-nchar(x))<3){x<-paste(x,"no title?"); n.x<-nchar(x)}
      x<-substring(x,1:n.x,1:n.x)
      if(x[1]!=kla) {x[2] <-paste(x[1],x[2],sep=""); x[1]<-"{"} #}
      en<-  which(cumsum((x==kla)-(x==kle))==0)[1] ## 090212??
      if(is.na(en)) {en<-length(x); x[en]<-paste(x[en],etag)} else x[en]<-etag
      x[1]<-paste("<a name=\"",command,i,"\">",atag,sep="")
      txt<-paste(x,collapse=""); sec<-c(sec,paste(x[(1+1):(en-1)],collapse=""))
      input[com.lines[i]+0:2]<-unlist(strsplit(txt,"\n"))
    }
    command.links<-paste("<a href=\"#",command,seq(along=com.lines),
                                        "\">",sec,"</a>\n",sep="") 
  }
  #   if(is.na(an)) {an<-1; x[an]<-paste(x[an],etag)} else x[an]<-h
  #   en<-which(cumsum((x=="{")-(x=="}"))==0)[1]
      
  sec.links<-command.links; sec.no<-com.lines
  atag<-"<h3>"; etag<-"</h3>"; command<-"subsection"
  command.n<-nchar(command)+2; command.links<-NULL  
  kla<-"{"; kle<-"}"
  ## print("STRUKTUR")
  if(0<length(com.lines<-grep(paste("^\\\\",command,sep=""),input))){
    sec<-NULL
    for(i in seq(along=com.lines)){
      txt<-input[com.lines[i]+0:2]
      txt<-paste(txt,collapse="\n"); n<-nchar(txt)   
      x<-sub("^ *","",substring(txt,command.n))
      if((n.x<-nchar(x))<3){x<-paste(x,"no title?"); n.x<-nchar(x)}
      x<-substring(x,1:n.x,1:n.x)
      if(x[1]!=kla) {x[2] <-paste(x[1],x[2],sep=""); x[1]<-"{"} #}
      en<-  which(cumsum((x==kla)-(x==kle))==0)[1] ## 090212??
      if(is.na(en)) {en<-length(x); x[en]<-paste(x[en],etag)} else x[en]<-etag
      x[1]<-paste("<a name=\"",command,i,"\">",atag,sep="")
      txt<-paste(x,collapse=""); sec<-c(sec,paste(x[(1+1):(en-1)],collapse=""))
      input[com.lines[i]+0:2]<-unlist(strsplit(txt,"\n"))
    }
    command.links<-paste("<a href=\"#",command,seq(along=com.lines),
                                        "\">",sec,"</a>\n",sep="") 
  }
  #   if(is.na(an)) {an<-1; x[an]<-paste(x[an],etag)} else x[an]<-h
  #   en<-which(cumsum((x=="{")-(x=="}"))==0)[1]
      
  subsec.links<-command.links; subsec.no<-com.lines
  atag<-"<h4>"; etag<-"</h4>"; command<-"subsubsection"
  command.n<-nchar(command)+2; command.links<-NULL  
  kla<-"{"; kle<-"}"
  ## print("STRUKTUR")
  if(0<length(com.lines<-grep(paste("^\\\\",command,sep=""),input))){
    sec<-NULL
    for(i in seq(along=com.lines)){
      txt<-input[com.lines[i]+0:2]
      txt<-paste(txt,collapse="\n"); n<-nchar(txt)   
      x<-sub("^ *","",substring(txt,command.n))
      if((n.x<-nchar(x))<3){x<-paste(x,"no title?"); n.x<-nchar(x)}
      x<-substring(x,1:n.x,1:n.x)
      if(x[1]!=kla) {x[2] <-paste(x[1],x[2],sep=""); x[1]<-"{"} #}
      en<-  which(cumsum((x==kla)-(x==kle))==0)[1] ## 090212??
      if(is.na(en)) {en<-length(x); x[en]<-paste(x[en],etag)} else x[en]<-etag
      x[1]<-paste("<a name=\"",command,i,"\">",atag,sep="")
      txt<-paste(x,collapse=""); sec<-c(sec,paste(x[(1+1):(en-1)],collapse=""))
      input[com.lines[i]+0:2]<-unlist(strsplit(txt,"\n"))
    }
    command.links<-paste("<a href=\"#",command,seq(along=com.lines),
                                        "\">",sec,"</a>\n",sep="") 
  }
  #   if(is.na(an)) {an<-1; x[an]<-paste(x[an],etag)} else x[an]<-h
  #   en<-which(cumsum((x=="{")-(x=="}"))==0)[1]
      
  subsubsec.links<-command.links; subsubsec.no<-com.lines
  atag<-"<br><b>"; etag<-"</b>"; command<-"paragraph"
  command.n<-nchar(command)+2; command.links<-NULL  
  kla<-"{"; kle<-"}"
  ## print("STRUKTUR")
  if(0<length(com.lines<-grep(paste("^\\\\",command,sep=""),input))){
    sec<-NULL
    for(i in seq(along=com.lines)){
      txt<-input[com.lines[i]+0:2]
      txt<-paste(txt,collapse="\n"); n<-nchar(txt)   
      x<-sub("^ *","",substring(txt,command.n))
      if((n.x<-nchar(x))<3){x<-paste(x,"no title?"); n.x<-nchar(x)}
      x<-substring(x,1:n.x,1:n.x)
      if(x[1]!=kla) {x[2] <-paste(x[1],x[2],sep=""); x[1]<-"{"} #}
      en<-  which(cumsum((x==kla)-(x==kle))==0)[1] ## 090212??
      if(is.na(en)) {en<-length(x); x[en]<-paste(x[en],etag)} else x[en]<-etag
      x[1]<-paste("<a name=\"",command,i,"\">",atag,sep="")
      txt<-paste(x,collapse=""); sec<-c(sec,paste(x[(1+1):(en-1)],collapse=""))
      input[com.lines[i]+0:2]<-unlist(strsplit(txt,"\n"))
    }
    command.links<-paste("<a href=\"#",command,seq(along=com.lines),
                                        "\">",sec,"</a>\n",sep="") 
  }
  #   if(is.na(an)) {an<-1; x[an]<-paste(x[an],etag)} else x[an]<-h
  #   en<-which(cumsum((x=="{")-(x=="}"))==0)[1]
      
  parsec.links<-command.links; parsec.no<-com.lines
  sec.typ<-rbind(cbind(c(0,sec.no),1),cbind(c(0,subsec.no),2),
                          cbind(c(0,subsubsec.no),3),cbind(c(0,parsec.no),4))
  sec.typ<-sec.typ[sec.typ[,1]!=0,,drop=FALSE]; contents<-" "
  if(length(sec.typ>2)){
    ind<-order(sec.typ[,1]); sec.typ<-sec.typ[ind,,drop=FALSE]
    links<-c(sec.links,subsec.links,subsubsec.links,parsec.links)[ind]
    # append a column with *section numbers
    sec.typ<-cbind(sec.typ,"0")
    sec.counter<-subsec.counter<-subsubsec.counter<-par.counter<-0
    for(i in 1:nrow(sec.typ)){
     if(sec.typ[i,2]=="1"){
        sec.counter<-sec.counter+1
        subsec.counter<-subsubsec.counter<-par.counter<-0
        sec.typ[i,3]<-sec.counter
     }
     if(sec.typ[i,2]=="2"){
        subsec.counter<-subsec.counter+1
        subsubsec.counter<-par.counter<-0
        sec.typ[i,3]<-paste(sec.counter,".",subsec.counter,sep="")
     }
     if(sec.typ[i,2]=="3"){
        subsubsec.counter<-subsubsec.counter+1
        par.counter<-0
        sec.typ[i,3]<-paste(sec.counter,".",subsec.counter,".",subsubsec.counter,sep="")
     }   
     if(sec.typ[i,2]=="4"){
        par.counter<-par.counter+1
        sec.typ[i,3]<-paste(sec.counter,".",subsec.counter,".",subsubsec.counter,".",
                                      par.counter,sep="")
     }   
    }
    # construct table of contents with links
    contents<-paste(sec.typ[,3],links)
  }

  ## if(  DEBUG-Flag gesetzt  ) print("head")
  head<-grep("^\\\\title|^\\\\author|^\\\\date",input)
  if(0<length(head)){
    h<-min(max(head)+5,length(input))
    head<-input[1:h]
    titel<-get.head.argument("title",head)[1]
    titel<-sub("Report: \\\\rule\\{(.*)\\}","Report: ................",titel)
    autor<-get.head.argument("author",head)[1]
    autor<-sub("File: \\\\jobname.rev",paste("File:",sub(".*/","",in.file)),autor)
    datum<-get.head.argument("date",head)[1]
    if(is.null(datum)) datum<-date()
    ## print(datum)
  } else {
    head<-""; titel<-paste("File:",in.file); autor<-"processed by weaveRhtml"
    datum<-date()
  }
  if(0<length(h<-grep("\\\\begin\\{document\\}",input)))
    input<-input[-(1:h[1])]
    titel.titel<-gsub("\n","--",paste(titel,collapse="--"))
    titel.titel<-gsub("<br>","--",titel.titel)
    titel.titel<-gsub("\\\\","--",titel.titel)
    input[1]<-paste(collapse="\n",
    "<!--  generated by weaveRhtml --><html><head>",
    "<meta content=\"text/html; charset=ISO-8859-1\">",
    "<title>",titel.titel,"</title></head>",
    "<body bgcolor=\"#FFFFFF\">",
    "<h1>",if(!is.null(titel))titel,"</h1>",
    "<h2>",if(!is.null(autor))autor,"</h2>",
    "<h3>",if(!is.null(datum))datum,"</h3>",
    "<h4>",paste(contents,collapse="<br>"),"</h4>"
  )


  if(0<length(h<-grep("\\\\myemph",input))){
    input<-transform.command.line("myemph",input,"<i>","</i>")
  }
  if(0<length(h<-grep("\\\\texttt",input))){
    input<-transform.command.line("texttt",input,"<code>","</code>")
  }
  if(0<length(h<-grep("\\\\emph",input))){
    input<-transform.structure.command("emph",input,
                       atag="<i>",etag="</i>",kla="{",kle="}")
  }
  if(0<length(h<-grep("\\\\textbf",input))){
    input<-transform.structure.command("textbf",input,
                       atag="<b>",etag="</b>",kla="{",kle="}")
  }

  input<-gsub("\\\\newpage","",input)
  input<-gsub("\\\\tableofcontents","",input)
  input<-gsub("\\\\raggedright","",input)
  #  input<-gsub("\\\\\\\\","<br>",input) # brute force
  cand.line <- grep("\\\\",input)
  if(0<length(cand.line)){ 
    cand <- input[cand.line]
    no.text <- c(grep("\\(.*\".*\\.*\"",cand),grep("\\(.*'.*\\.*'",cand)) # ))
    idx <- seq(along=cand); idx <- idx[ !(idx %in% no.text)]
    if(0<length(idx)) cand[idx] <- gsub("\\\\\\\\","<br>",cand[idx])
    input[cand.line] <- cand
  }
  h<-grep("\\\\maketitle|\\\\author|\\\\date|\\\\title|\\\\end\\{document\\}",input)
  if(0<length(h)) input<-input[-h]
    txt<-input
    ind.Rweb<-grep("^<a name=\"codechunk.*<i>&lt;Rweb",txt) ; txt[ind.Rweb]
    ind.p   <-grep(paste("^<","p>",sep=""),txt) ; txt[ind.p]
    if(length(ind.p)>0){
      ind.Rweb.codes<-lapply(ind.Rweb,function(x) (x+1):(ind.p[ind.p>x][1]-1))
      ## if( DEBUG-Flag gesetzt ) print(ind.Rweb.codes)
      Rwebbegin<-paste(c(  # Rweb start 
        '<hr SIZE=3><form onSubmit =" return checkData(this)"',
        ## '  action="http://rweb.stat.umn.edu/cgi-bin/Rweb/Rweb.cgi" ', ## out of order
        '  action="http://pbil.univ-lyon1.fr/cgi-bin/Rweb/Rweb.cgi" ',
        '  enctype="multipart/form-data" method="post" target="_blank"><p>',
        '<label for="filedata" size=20>data frame (load LOCAL file): </label>',
        '<input type="file" name="FileData" id="filedata" size=30 >',
        '<br><br>',
        '<label for="urldata" size=20>data frame (via www -> URL): </label>',
        '<input type="TEXT" name="URLData" size=40 id="urldata" YyYyY>',
        '<br><br>',
        '<label for="Rcode">R-Code:</label>',
        '<textarea name="Rcode" rows="XxXxX" cols="80" id="Rcode">\n'
        ),collapse=" ")
      Rwebend<-  paste(c( # Rweb end
        '</textarea><p>',
        '<input type="submit" value="compute via Rweb">',
        '<input type="reset" value="reset">',
        '</form><hr SIZE=3>',
        ' '),collapse=" ")
      for(i in seq(along=ind.Rweb.codes)){
        h<-sub("<code>","",txt[ind.Rweb.codes[[i]]])
        h<-sub("<br>","",h)
        h<-sub("</code>","",h)
        h<-gsub("[&]nbsp[;]"," ",h)
        data.http<-""    # ; data.loc<-""
        ind<-length(grep("read.table",h))
        if(0<length(ind)) {
          h[ind]<-paste("#",h[ind],"# choose data by input field!!")
          ind<-ind[1]
          if(0<length(grep("http:",h[ind]))) {
             data.http<-sub(".*(http)","\"\\1",h[ind])
             data.http<-paste("value=",sub(".*(\".*\").*$","\\1",data.http),sep="")
          }
        }
        rb<-sub("YyYyY",data.http,Rwebbegin)
        h[1]<-paste(sub("XxXxX",as.character(length(h)),rb),h[1],sep="")
        h[length(h)]<-paste(h[length(h)],Rwebend)
        txt[ind.Rweb.codes[[i]]]<-h
      }
      input<-txt
    }


  input<-gsub("<br>","<br>\n",input)

  if(missing(out.file)||in.file==out.file){
    out.file<-sub("\\.([A-Za-z])*$","",in.file)
  }
  if(0==length(grep("\\.html$",out.file)))
    out.file<-paste(out.file,".html",sep="")
  ## out.file<-"/home/wiwi/pwolf/tmp/out.html"
  get("cat","package:base")(input,sep="\n",file=out.file)
  cat("weaveRhtml process finished\n")

}



