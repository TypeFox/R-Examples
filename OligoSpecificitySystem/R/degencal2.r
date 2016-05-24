"degencal2"<-function(){
temp<-function(){
calc2()

symbols(c(35,65,35,65), c(50,50,50,50), circles = c(30,30,29,29),lwd=c(4,4,4,4),bg=c("white","white","white","white"),fg=c("blue","darkgreen",NA,NA), inches = FALSE,xlim=c(1,100),ylim=c(1,100),xaxt="n",yaxt="n",ylab=NA,xlab=NA)
text(14,90,"Oligonucleotide 1",col="blue",cex=1)
text(80,90,"Oligonucleotide 2",col="darkgreen",cex=1)

text(48,56,c("Degenerate oligonucleotide matchs"))
if(nbsequence=="optional") text(48,50,paste(sum(zz1)," unique sequences"))
if(nbsequence!="optional") text(48,50,paste(round((sum(zz1)/as.numeric(nbsequence))*100)," % of unique sequences"))

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