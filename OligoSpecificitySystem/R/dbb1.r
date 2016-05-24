"dbb1"<-function()
  {
    OnOk()
    txt1<-0;txt1<<-txt1
    fil1=if (interactive()) choose.files(filters = Filters["All",],caption="Select the 1st oligonucleotide data base")
    print("Please wait while loading")
    ONDB1name<<-fil1
    txt1<-scan(file=fil1,what="character",fill=TRUE,quiet=TRUE)
    a<-which(txt1==separator)
    if(length(a)==0) txt1<-unlist(strsplit(txt1,""))
    if(length(a)==0) a<-which(txt1==separator)
    if(length(a)==0) tkmessageBox(message="The specified separator is not found in your oligonucleotide database.")
    if(length(a)==0) stop()
    b<-vector(length=length(a))
    for (i in 1:c(length(a)-1)){
    if(i/1000==round(i/1000)) waitGUI(i,c(length(a)-1))
    b[i]<-paste(txt1[c(a[i]+1):c(a[i+1]-1)],collapse="")}
    b[length(a)]<-paste(txt1[c(a[length(a)]+1):length(txt1)],collapse="")
    if(length(a)>=1000) dev.off()
    print(paste("Your oligonucleotide database has",length(b),"sequences"))
    c<-levels(factor(b))
    txt1<<-c
    if (length(c)!=length(b)) print("but has sequences in several copies") 
    if (length(c)!=length(b))print(paste("Your oligonucleotide database has",c(length(c)),"unique sequences"))
    print("DataBase 1 successfully imported")
  }
  
  