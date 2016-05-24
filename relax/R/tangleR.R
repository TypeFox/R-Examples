tangleR<-
function(in.file,out.file,expand.roots=NULL,expand.root.start=TRUE,
         insert.comments=TRUE,add.extension=TRUE){
  # german documentation of the code:
  # look for file webR.pdf, P. Wolf 050204
  AtSignLineBegin <- "^@"
  pat.chunk.header <- paste("^<","<(.*)>",">=",sep="")
  pat.chunk.header.line <- paste("(.*)",pat.chunk.header,"(.*)",sep="")
  pat.use.chunk      <- paste("<","<(.*)>",">",sep="")
  pat.use.chunk.line <- paste("(.*)",pat.use.chunk,"(.*)",sep="")
  ReplaceStringChHeader <- paste("cOdE","dEf",sep="")
  ReplaceStringChUse <- paste("uSe","ChUnK",sep="")
  ReplaceStringEscLGB <- paste("DoSpOpenKl","-esc",sep="")
  ReplaceStringEscRGB <- paste("DoSpCloseKl","-esc",sep="")
  BREAK <- paste("bRe","Ak",sep="")
  EscLeftLowerBrackets <- paste("@<","<",sep="")
  EscRightGreaterBrackets <- paste("@>",">",sep="")

  if(!file.exists(in.file)) in.file<-paste(in.file,"rev",sep=".")
  if(!file.exists(in.file)){
    cat(paste("ERROR:",in.file,"not found!!??\n"))
    return("Error in tangle: file not found")
  }
  # code.ch<-scan(in.file,sep="\n",what=" ",blank.lines.skip =FALSE)
  code.ch<-readLines(in.file) # 2.1.0

  # append empty text chunk
  code.ch<-c(code.ch,"@")
  # check for code chunk starts
  code.a <- grep(pat.chunk.header,code.ch)
  if(0==length(code.a)){return("Warning: no code found!!!!")}
  # integrate "@" lines 090220
  aa<-as.list(code.ch) 
  aa[code.a]<-lapply(aa[code.a],function(x) c("@",x)) 
  code.ch<-unlist(aa)
  code.n <-length(code.ch)
  # find code chunk starts again
  code.a<- grep(pat.chunk.header,code.ch)
  # find code chunk endings
  #  text chunk starts
  code.z    <-grep(AtSignLineBegin,code.ch)
  #  text chunk starts that follow code chunks
  code.z    <-unlist(sapply(code.a, function(x,y) min(y[y>x]),code.z))
  # find positions of change from text to code chunk or from code chunk to text chunk
  change    <-rep(0,code.n); change[c(code.a ,code.z)] <- 1
  # extract code chunks
  code.ch   <-code.ch[1==(cumsum(change)%%2)]
  # save number of lines
  code.n    <-length(code.ch)

  code.ch<-gsub(EscRightGreaterBrackets,
                ReplaceStringEscRGB,
                gsub(EscLeftLowerBrackets,ReplaceStringEscLGB,code.ch))

  # remove disturbing characters in a code chunk definition line
  code.ch<-gsub(pat.chunk.header.line,
                paste(ReplaceStringChHeader,"\\2",sep=""),
                code.ch)
  ReplaceStringUseLine <- paste("\\1",BREAK,ReplaceStringChUse,"\\2",BREAK,"\\3",sep="")
  repeat{
    if(0==length(cand<-grep(pat.use.chunk,code.ch))) break
    code.ch<-unlist(strsplit(
                gsub(pat.use.chunk.line,ReplaceStringUseLine,code.ch),
                BREAK 
             ))
  }
  code.ch<-code.ch[code.ch!=""]
  code.n<-length(code.ch)
  if(exists("DEBUG")) print(code.ch)

  line.typ  <-rep("C",code.n)
  code.a    <-grep(ReplaceStringChHeader,code.ch)
  code.ch[code.a]<-substring(code.ch[code.a],8)
  line.typ[code.a]<-"D"
  code.use  <-grep(ReplaceStringChUse,code.ch)
  code.ch[code.use]<-substring(code.ch[code.use],9)
  line.typ[code.use]<-"U"

  code.out<-NULL

  def.names<-code.ch[code.a]
  use.names<- code.ch[code.use]
  code.z<-c(if(length(code.a)>1) code.a[-1]-1, code.n)
  code.ch<-paste(line.typ,code.ch,sep="")
  if(exists("DEBUG")) print(code.ch)
  if(expand.root.start){
    if(exists("DEBUG")) cat("bearbeite start\n")
    if(insert.comments) code.out<-c(code.out,"#0:","##start:##")
    if(any(ch.no <-def.names=="start")){
       ch.no     <-seq(along=def.names)[ch.no]; rows<-NULL
       for(i in ch.no)
          if((code.a[i]+1)<=code.z[i]) rows<-c(rows, (code.a[i]+1):code.z[i])
       code.stack<-code.ch[rows]
       repeat{
         if(0==length(code.stack))break
         if("C"==substring(code.stack[1],1,1)){
           n.lines<-sum(cumprod("C"==substring(code.stack,1,1)))
           code.out<-c(code.out, substring(code.stack[1:n.lines],2))
           code.stack<-code.stack[-(1:n.lines)]
         }else{
           if(any(found<-def.names==substring(code.stack[1],2))){
             found<-seq(along=def.names)[found]; rows<-NULL
             for(no in found){
               row.no<-c((code.a[no]+1),code.z[no])
               if(row.no[1]<=row.no[2]) rows<-c(rows,row.no[1]:row.no[2])
             }
             code.stack<-c(code.ch[rows],code.stack[-1])
             cat(found,", ",sep="")
           } else code.stack <-code.stack[-1] # ignore not defined chunks!
           # 051219 
         }

       }
    }
    if(insert.comments) code.out<-c(code.out,"##:start##","#:0")

  }
  root.no<-is.na(match(def.names,use.names))&def.names!="start"
  root.no<-seq(along=root.no)[root.no]
  roots  <-def.names[root.no]
  #if(!is.null(expand.roots)){ # verwende nur wirkliche Wurzeln
  #  h<-!is.na(match(roots,expand.roots))
  #  roots<-roots[h]; root.no<-root.no[h]
  #}
  if(!is.null(expand.roots)){ # verwende gegebene(n) Knoten
    h<-!is.na(match(def.names,expand.roots))
    root.no<-seq(def.names)
    roots<-def.names[h]; root.no<-root.no[h]
  }
  if(0==length(roots)) cat("Sorry, no chunk for expanding found!\n")
  if(exists("DEBUG")) print(roots)

  if(exists("DEBUG")) cat("bearbeite Sektion-Nr./Name\n")
  roots<-unique(roots)
  for(r in seq(along=roots)){
    # if(exists("DEBUG")) cat(root.no[r],":",roots[r],", ",sep="")
    if(exists("DEBUG")) cat(which(def.names==roots[r]),":",roots[r],", ",sep="")
    if(any(ch.no <-def.names==roots[r])){
      ch.no     <-seq(along=def.names)[ch.no]; rows<-NULL
      if(insert.comments) code.out<-c(code.out,
            # # # paste("#",root.no[r],":",sep=""), new 071114
                          paste("##",roots[r],":##",sep=""))
      for(i in ch.no){
         if((code.a[i]+1)<=code.z[i]){ 
           # # # rows<-c(rows, (code.a[i]+1):code.z[i]) # new:
              h<-code.a[i]+1
              rows<-c(rows, h:code.z[i])
              if(insert.comments) code.ch[h]<-paste("C#",i,":NeWlInE",code.ch[h],sep="")
              h<-code.z[i]
              if(insert.comments) code.ch[h]<-paste(code.ch[h] ,"NeWlInEC#:",i,sep="")
         }
      }
      code.stack<-code.ch[rows]
      code.stack<-unlist(strsplit(code.stack,"NeWlInE")) # new   
      repeat{
       if(0==length(code.stack))break
       if("C"==substring(code.stack[1],1,1)){
         n.lines<-sum(cumprod("C"==substring(code.stack,1,1)))
         code.out<-c(code.out, substring(code.stack[1:n.lines],2))
         code.stack<-code.stack[-(1:n.lines)]
       }else{
         def.line<-substring(code.stack[1],2)
         if(any(found<-def.names==def.line)){
           code.stack<-code.stack[-1]
           found<-rev(seq(along=def.names)[found])
           for(no in found){
             row.no<-c((code.a[no]+1),code.z[no])
             if(row.no[1]<=row.no[2]){
               if(insert.comments){ 
                 code.stack<-c(paste("C#"  ,no,":"  ,sep=""      ),
                               paste("C##" ,def.line,":##",sep=""), 
                               code.ch[row.no[1]:row.no[2]],
                               paste("C##:",def.line, "##",sep=""),
                               paste("C#:" ,no      ,sep=""      ),
                               code.stack)
               }else{
                 code.stack<-c(code.ch[row.no[1]:row.no[2]], code.stack)
               }
             }
           } # end of for
         } else code.stack <-code.stack[-1] # ignore not defined chunks!
         # 051219 
       }

      }
      if(insert.comments) code.out<-c(code.out,paste("##:",roots[r],"##",sep="")
                 # # #    ,paste("#:",root.no[r],sep="")
                          )
    }
  }

  code.out<-gsub(ReplaceStringEscRGB,">>",gsub(ReplaceStringEscLGB,"<<",code.out))

  if(missing(out.file)||in.file==out.file){
    out.file<-sub("\\.([A-Za-z])*$","",in.file)
  }
  if(add.extension&&0==length(grep("\\.R$",out.file)))
    out.file<-paste(out.file,".R",sep="")
  get("cat","package:base")(code.out,sep="\n",file=out.file)
  cat("tangle process finished\n")

}

