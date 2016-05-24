spss_export <-
function(table,file.dict=NULL,file.save=NULL,var.keep="ALL",file.runsyntax="C:/Archivos de programa/SPSS/runsyntx.exe",file.data=NULL,run.spss=TRUE,dec="."){


  if (class(table)!="data.frame") table<-as.data.frame(table)  

  names<-names(table) 
  names<-gsub(" ","_",names)
  names<-gsub("\\(","_",names)
  names<-gsub("\\)","_",names)
  names<-gsub("\\.","_",names)
  names(table)<-names
    
  ## If the variable don't have format.SPSS it runs format_corrector in those variables.
  formats<-unlist(lapply(table,function(x){res<-attr(x,"format.SPSS"); if (is.null(res)) return("FORMAT.NULL") else return(res)})) 
    
  no.format<-grep("FORMAT.NULL",formats)
  if (length(no.format)>0) table[,no.format]<-format_corrector(table[,no.format])

  formats<-unlist(lapply(table,function(x) attr(x,"format.SPSS")))    

  table.spss<-table
  for (i in grep("DATE",formats)){
    table.spss[,i]<-as.character(table.spss[,i])
    table.spss[,i]<-gsub("\\)","",table.spss[,i])
    table.spss[,i]<-gsub("\\(","",table.spss[,i])      
  }

  f.tab<-function(x){
    if (any(is.character(x) | is.factor(x))){
      old.attr<-attributes(x)
      x<-gsub("\t"," ",as.character(x)) 
      old.attr$class<-"character"
      attributes(x)<-old.attr
      return(x)
    } else return(x)
  }
  table.spss<-as.data.frame(lapply(table.spss, f.tab))
      
  
  if (is.null(file.data)) file.data<-paste(tempdir(),"dades.txt",sep="\\")
  write.table(table.spss, file = file.data, append = FALSE, quote = FALSE, sep = "\t",
	       	  eol = "\n", na = "", dec = dec, row.names = FALSE, col.names = TRUE)
  
 
  formats.table<-cbind(names,formats)
  aux=""
  for (k in 1:nrow(formats.table))
    aux<-paste(aux,"\n\t",paste(formats.table[k,],collapse="\t"),sep="")

  syntax<-
  paste("
GET DATA  /TYPE = TXT
 /FILE = '",file.data,"'
 /DELCASE = LINE
 /QUALIFIER = ''
 /DELIMITERS = \"\\t\"
 /ARRANGEMENT = DELIMITED
 /FIRSTCASE = 2
 /IMPORTCASE = ALL
 /VARIABLES =",
aux,"
.
CACHE.
EXECUTE.","\n\n",sep="")

  write(syntax,file.syntax<-file.path(tempdir(),"syntax.sps"))

  if (!is.null(file.dict)){
    for (i in 1:length(file.dict)){
      file.append(file.syntax,file.dict[i])
      write("\n\n",file.xxx<-file.path(tempdir(),"xxx.sps"))
      file.append(file.syntax,file.xxx)
    }
  } else{
    internal_dictionary(table,file.xxx<-file.path(tempdir(),"xxx.sps"))
    file.append(file.syntax,file.xxx)
  }
 
  if (!is.null(file.save) && length(file.save)>0){
    syntax.save<-
    paste(
"SAVE OUTFILE='",file.save,"'
 /KEEP ",var.keep,"
 /COMPRESSED.\n\n",sep="")

    write(syntax.save,tempfile.save<-file.path(tempdir(),"syntax.save.sps"))
    file.append(file.syntax,tempfile.save)

    write(paste("GET FILE='",file.save,"'.",sep=""),file.xxx<-file.path(tempdir(),"xxx.sps"))
    file.append(file.syntax,file.xxx)

  }

  syntax<-scan(file.syntax,what="character",sep="@",quote=NULL)
  index.miss<-grep("MISSING VALUES ",syntax)
  if (length(index.miss)>0){
    index.miss.del<-!lapply(strsplit(syntax[index.miss]," "),function(x) x[3])%in%names
    index.miss.del<-index.miss[index.miss.del]
    syntax<-syntax[-index.miss.del]
    write(syntax,file.syntax)
  }
  

  file.show(file.syntax)  
  
  
  if (run.spss){


    shell.exec(file.syntax)
  

    cmd=paste("\"",file.runsyntax,"\" ",file.syntax,sep="")
    shell(cmd,mustWork=FALSE)

  }

  
}
