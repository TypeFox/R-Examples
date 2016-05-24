"glob3d"<-function(){
temp<-function(){
calc3()
symbols(c(62,78,69,62,78), c(42,42,66,42,42), circles = c(25,25,25,25,25),lwd=c(4,4,4,4,4),bg=c("white","white","white",NA,NA),fg=c("darkgreen","darkred","blue","darkgreen","darkred"), inches = FALSE,xlim=c(1,100),ylim=c(1,100),xaxt="n",yaxt="n",ylab=NA,xlab=NA)
text(43,10,"Oligonucleotide 2",col="darkgreen",cex=1)
text(87,10,"Oligonucleotide 3",col="darkred",cex=1)
text(70,98,"Oligonucleotide 1",col="blue",cex=1)
text(69,80,"x1",cex=2.5)
text(45,40,"x2",cex=2.5)
text(95,40,"x3",cex=2.5)
text(52,60,"x4",cex=2.5)
text(88,60,"x5",cex=2.5)
text(69,28,"x6",cex=2.5)
text(69,54,"x7",cex=2.5)

normal<-function(){
text(20,93,c(paste("Effectif of x1 = ",zz1[1])))
text(20,87,c(paste("Effectif of x2 = ",zz1[2])))
text(20,81,c(paste("Effectif of x3 = ",zz1[3])))
text(20,75,c(paste("Effectif of x4 = ",zz1[4])))
text(20,69,c(paste("Effectif of x5 = ",zz1[5])))
text(20,63,c(paste("Effectif of x6 = ",zz1[6])))
text(20,57,c(paste("Effectif of x7 = ",zz1[7])))}

normal_prct<-function(){
text(20,93,c(paste("Effectif of x1 = ",c(round((zz1[1]/as.numeric(nbsequence))*100)))))
text(20,87,c(paste("Effectif of x2 = ",c(round((zz1[2]/as.numeric(nbsequence))*100)))))
text(20,81,c(paste("Effectif of x3 = ",c(round((zz1[3]/as.numeric(nbsequence))*100)))))
text(20,75,c(paste("Effectif of x4 = ",c(round((zz1[4]/as.numeric(nbsequence))*100)))))
text(20,69,c(paste("Effectif of x5 = ",c(round((zz1[5]/as.numeric(nbsequence))*100)))))
text(20,63,c(paste("Effectif of x6 = ",c(round((zz1[6]/as.numeric(nbsequence))*100)))))
text(20,57,c(paste("Effectif of x7 = ",c(round((zz1[7]/as.numeric(nbsequence))*100)))))
text(20,51,"   ")
text(20,45,"in % of sequences")
}

if(nbsequence=="optional") normal()
if(nbsequence!="optional") normal_prct()
}
tt11 <- tktoplevel()
tkwm.title(tt11,"All information on oligonucleotides data base")
img <- tkrplot(tt11,fun=temp,hscale=1.5,vscale=1.5)
CopyToClip <- function(){tkrreplot(img)}
copy.but <- tkbutton(tt11,text="Copy to Clipboard",command=CopyToClip)
tkgrid(img)
tkgrid(copy.but)
tkgrid(tklabel(tt11,text="  ")) 
}