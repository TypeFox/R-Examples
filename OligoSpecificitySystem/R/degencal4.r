"degencal4"<-function(){
rr<-1
  if (length(txt1)!=1 & length(txt2)!=1 & length(txt3)!=1 & length(txt4)!=1) rr<- 0
  
  
  if (nbsequence=="optional")tkmessageBox(message=paste("Degenerate oligonucleotide matchs ",c(c(length(levels(factor(c(txt1,txt2,txt3,txt4)))))-rr)," unique sequences"))
  if (nbsequence!="optional")tkmessageBox(message=paste("Degenerate oligonucleotide matchs ",c(round((c(c(c(length(levels(factor(c(txt1,txt2,txt3,txt4)))))-rr)/as.numeric(nbsequence))*100)))," % of unique sequences"))
}