"impr2d"<-function(){
temp<-function(){
calc2()
symbols(c(40,60,40), c(60,60,60), circles = c(30,30,30),lwd=c(4,4,4),bg=c("white","white",NA),fg=c("blue","darkgreen","blue"), inches = FALSE,xlim=c(1,100),ylim=c(1,100),xaxt="n",yaxt="n",ylab=NA,xlab=NA)
text(19,100,"Oligonucleotide 1",col="blue",cex=1)
text(75,100,"Oligonucleotide 2",col="darkgreen",cex=1)
text(20,59,"x1",cex=2)
text(80,59,c("x2"),cex=2)
if(nbsequence=="optional") text(50,16,c(paste("The improvement of 1st degeneracy (x1) = ",zz1[1])))
if(nbsequence=="optional") text(50,10,c(paste("The improvement of 2nd degeneracy (x2) = ",zz1[3])))

if(nbsequence!="optional") text(50,16,c(paste("The improvement of 1st degeneracy (x1) = ",c(round((zz1[1]/as.numeric(nbsequence))*100))," %")))
if(nbsequence!="optional") text(50,10,c(paste("The improvement of 2nd degeneracy (x2) = ",c(round((zz1[3]/as.numeric(nbsequence))*100))," %")))
}
tt11 <- tktoplevel()
tkwm.title(tt11,"Improvements of degeneracies on oligonucleotide")
img <- tkrplot(tt11,fun=temp,hscale=1.5,vscale=1.5)
CopyToClip <- function(){tkrreplot(img)}
copy.but <- tkbutton(tt11,text="Copy to Clipboard",command=CopyToClip)
tkgrid(img)
tkgrid(copy.but)
tkgrid(tklabel(tt11,text="  ")) 
}