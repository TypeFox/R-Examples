"dbb2"<-function()
  { OnOk()
    txt2<-0;txt2<<-txt2
    fil2=if (interactive()) choose.files(filters = Filters["All",],caption="Select the 2nd oligonucleotide data base")
    ONDB2name<<-fil2
    print("Please wait while loading")
    txt2<-scan(file=fil2,what="character",fill=TRUE,quiet=TRUE)
    a<-which(txt2==separator)
    if(length(a)==0) txt2<-unlist(strsplit(txt2,""))
    if(length(a)==0) a<-which(txt2==separator)
    if(length(a)==0) tkmessageBox(message="The specified separator is not found in your oligonucleotide database.")
    if(length(a)==0) stop()
    b<-vector(length=length(a))
    for (i in 1:c(length(a)-1)){
    if(i/1000==round(i/1000)) waitGUI(i,c(length(a)-1))
    b[i]<-paste(txt2[c(a[i]+1):c(a[i+1]-1)],collapse="")}
    b[length(a)]<-paste(txt2[c(a[length(a)]+1):length(txt2)],collapse="")
    if(length(a)>=1000) dev.off()
        print(paste("Your oligonucleotide database has",length(b),"sequences"))
    c<-levels(factor(b))
    txt2<<-c
    if (length(c)!=length(b)) print("but has sequences in several copies") 
    if (length(c)!=length(b))print(paste("Your oligonucleotide database has",c(length(c)),"unique sequences"))
        print("DataBase 2 successfully imported")
  }