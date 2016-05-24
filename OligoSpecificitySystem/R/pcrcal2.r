"pcrcal2"<-function(){
temp<-function(){
calc2()
 ff<-zz1[3]
 ff<<-ff
symbols(c(40,60,40), c(50,50,50), circles = c(30,30,30),lwd=c(4,4,4),bg=c("white","white",NA),fg=c("blue","darkgreen","blue"), inches = FALSE,xlim=c(1,100),ylim=c(1,100),xaxt="n",yaxt="n",ylab=NA,xlab=NA)
text(19,90,"Oligonucleotide 1",col="blue",cex=1)
text(75,90,"Oligonucleotide 2",col="darkgreen",cex=1)

normal<-function(){text(50,61,c("Oligonucleotides set"))
text(50,55,c("matchs"))
text(50,49,ff)
text(50,43,c("common sequences"))
text(50,37,"(intersection)")}

normal_prct<-function(){text(50,61,c("Oligonucleotides set"))
text(50,55,c("matchs"))
text(50,49,paste(c(round((ff/as.numeric(nbsequence))*100))," % of",sep=""))
text(50,43,c("common sequences"))
text(50,37,"(intersection)")}

if(nbsequence=="optional") normal()
if(nbsequence!="optional") normal_prct()

}
tt11 <- tktoplevel()
tkwm.title(tt11,"Oligonucleotides system")
img <- tkrplot(tt11,fun=temp,hscale=1.5,vscale=1.5)
CopyToClip <- function(){tkrreplot(img)}
copy.but <- tkbutton(tt11,text="Copy to Clipboard",command=CopyToClip)
tkgrid(img)
tkgrid(copy.but)
tkgrid(tklabel(tt11,text="  ")) 

}