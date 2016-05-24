"dbb4"<-function()
  {
    OnOk()
    txt4<-0;txt4<<-txt4

    locate.dir<-if (interactive()) choose.dir(getwd(), "Select folder containing oligonucleotide databases files")
    myfiles<-paste(locate.dir,list.files(locate.dir),sep="\\")
    myfiles<<-myfiles
    dimtemp<- vector("list", c(length(myfiles)+3))
    print(paste(length(myfiles)," databases have been found",sep=""))
    print("")    
    for (j in 1:length(myfiles)) { 
      txttemp<-0
      txttemp<-scan(file=myfiles[j],what="character",fill=TRUE,quiet=TRUE)
      a<-which(txttemp==separator)
      if(length(a)==0) txttemp<-unlist(strsplit(txttemp,""))
    if(length(a)==0) a<-which(txttemp==separator)
    if(length(a)==0) tkmessageBox(message="The specified separator is not found in your oligonucleotide database.")
    if(length(a)==0) stop()
    b<-vector(length=length(a))
    for (i in 1:c(length(a)-1)){
    if(i/1000==round(i/1000)) waitGUI(i,c(length(a)-1))
    b[i]<-paste(txttemp[c(a[i]+1):c(a[i+1]-1)],collapse="")}
    b[length(a)]<-paste(txttemp[c(a[length(a)]+1):length(txttemp)],collapse="")
    if(length(a)>=1000) dev.off()
      print(paste("Database ",j," contains ",length(b)," sequences",sep=""))
      if (length(levels(factor(b)))!=length(b)) print("but has sequences in several copies") 
    if (length(levels(factor(b)))!=length(b)) print(paste("Database ",j," contains ",length(levels(factor(b)))," unique sequences",sep=""))
       
   print("")
      dimtemp[[j]]<-levels(factor(b))}
     dimtemp<<-dimtemp
     txt4<-unlist(dimtemp)
      txt4<<-txt4
     val<-vector(length=4);val[]<-0
     val[1]<-length(myfiles)
     if (length(txt1)!=1) val[2]<-1
     if (length(txt2)!=1) val[3]<-1
     if (length(txt3)!=1) val[4]<-1
     val<-sum(val)
     val<<-val
     
   ONDB4name<<-paste("Loaded from ",locate.dir,sep="")
  
   print("Your oligonucleotide databases have been successfully imported")
   
  }
  
