"pcrcal3"<-function(){
temp<-function(){
calc3()
symbols(c(42,58,49,42,58), c(42,42,66,42,42), circles = c(25,25,25,25,25),lwd=c(4,4,4,4,4),bg=c("white","white","white",NA,NA),fg=c("darkgreen","darkred","blue","darkgreen","darkred"), inches = FALSE,xlim=c(1,100),ylim=c(1,100),xaxt="n",yaxt="n",ylab=NA,xlab=NA)
text(20,10,"Oligonucleotide 2",col="darkgreen",cex=1)
text(82,10,"Oligonucleotide 3",col="darkred",cex=1)
text(50,98,"Oligonucleotide 1",col="blue",cex=1)
text(49,54,c("x7"),cex=2.5)
if(nbsequence=="optional") text(50,2,paste("Your oligonucleotides set matchs ", zz1[7], " common sequences (x7)"))
if(nbsequence!="optional") text(50,2,paste("Your oligonucleotides set matchs ", c(round((zz1[7]/as.numeric(nbsequence))*100)), " % of common sequences (x7)"))
}
tt11 <- tktoplevel()
tkwm.title(tt11,"Oligonucleotides system")
img <- tkrplot(tt11,fun=temp,hscale=1.8,vscale=1.8)
CopyToClip <- function(){tkrreplot(img)}
copy.but <- tkbutton(tt11,text="Copy to Clipboard",command=CopyToClip)
tkgrid(img)
tkgrid(copy.but)
tkgrid(tklabel(tt11,text="  ")) 
}