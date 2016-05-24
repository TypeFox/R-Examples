"dbb3"<-function(){
    OnOk()
    txt3<-0;txt3<<-txt3
    fil3=if (interactive()) choose.files(filters = Filters["All",],caption="Select the 3rd oligonucleotide data base")
    ONDB3name<<-fil3
    print("Please wait while loading")
    txt3<-scan(file=fil3,what="character",fill=TRUE,quiet=TRUE)
    a<-which(txt3==separator)
    if(length(a)==0) txt3<-unlist(strsplit(txt3,""))
    if(length(a)==0) a<-which(txt3==separator)
    if(length(a)==0) tkmessageBox(message="The specified separator is not found in your oligonucleotide database.")
    if(length(a)==0) stop()
    b<-vector(length=length(a))
    for (i in 1:c(length(a)-1)){
    if(i/1000==round(i/1000)) waitGUI(i,c(length(a)-1))
    b[i]<-paste(txt3[c(a[i]+1):c(a[i+1]-1)],collapse="")}
    b[length(a)]<-paste(txt3[c(a[length(a)]+1):length(txt3)],collapse="")
        if(length(a)>=1000) dev.off()
         print(paste("Your oligonucleotide database has",length(b),"sequences"))
    c<-levels(factor(b))
    txt3<<-c
    if (length(c)!=length(b)) print("but has sequences in several copies") 
    if (length(c)!=length(b))print(paste("Your oligonucleotide database has",c(length(c)),"unique sequences"))
         print("DataBase 3 successfully imported")
  }