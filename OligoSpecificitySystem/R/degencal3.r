"degencal3"<-function(){
temp<-function(){
calc3()
symbols(c(42,58,49,42,58,49), c(42,42,66,42,42,66), circles = c(25,25,25,24,24,24),lwd=c(4,4,4,4,4,4),bg=c("white","white","white","white","white","white"),fg=c("darkgreen","darkred","blue",NA,NA,NA), inches = FALSE,xlim=c(1,100),ylim=c(1,100),xaxt="n",yaxt="n",ylab=NA,xlab=NA)
text(50,98,"Oligonucleotide 1",col="blue",cex=1)
text(20,10,"Oligonucleotide 2",col="darkgreen",cex=1)
text(82,10,"Oligonucleotide 3",col="darkred",cex=1)
text(50,50,c("Degenerate oligonucleotide matchs"))
if(nbsequence=="optional") text(50,44,paste(sum(zz1)," unique sequences"))
if(nbsequence!="optional") text(50,44,paste(round(((sum(zz1)/as.numeric(nbsequence))*100))," % of unique sequences"))
}
tt11 <- tktoplevel()
tkwm.title(tt11,"Degenrate oligonucleotide with mismatches")
img <- tkrplot(tt11,fun=temp,hscale=1.5,vscale=1.5)
CopyToClip <- function(){tkrreplot(img)}
copy.but <- tkbutton(tt11,text="Copy to Clipboard",command=CopyToClip)
tkgrid(img)
tkgrid(copy.but)
tkgrid(tklabel(tt11,text="  ")) 
}