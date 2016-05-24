"pcrcal4"<-function(){
val<-vector(length=4);val[]<-0
     val[1]<-length(myfiles)
     if (length(txt1)!=1) val[2]<-1
     if (length(txt2)!=1) val[3]<-1
     if (length(txt3)!=1) val[4]<-1
     val<-sum(val)
     val<<-val
a=summary(factor(c(txt1,txt2,txt3,txt4)),maxsum=length(levels(factor(c(txt1,txt2,txt3,txt4)))))

if(length(which(c(txt1,txt2,txt3,txt4)==0))==val) ff<- length(which(a==val))-1
if(length(which(c(txt1,txt2,txt3,txt4)==0))!=val) ff<- length(which(a==val))

if (nbsequence=="optional") tkmessageBox(message=paste("Oligonucleotides set matchs ",ff," common sequences"))

if (nbsequence!="optional") tkmessageBox(message=paste("Oligonucleotides set matchs ",c(round((ff/as.numeric(nbsequence))*100))," % of common sequences"))
}