Plot.Data.Events <-
function (yy,paciente,inicio,dias,censored,
                            especiales,colevent = "red", colcensor = "blue") 
{
p<-ncol(yy)
N<-nrow(yy)
nn<-length(paciente)
n<-nn
dev.new()
par(bg = "white")
plot(inicio,paciente,xlim=c(-1,(max(dias+inicio)+1)),ylim=c(-0.5,n),
     xlab="Time",ylab="Unit",pch=19,cex=0.4,col = "dark blue", 
     sub = R.version.string)
title(main=list("Graphical Representation of Recurrent Event Data",
      cex = 0.8, font = 2.3,col = "dark blue"))
mtext("Research Group: AVANCE USE R!",
       cex = 0.7, font = 2,col = "dark blue",line=1)
mtext("Software made by: Carlos Martinez",
       cex = 0.6, font = 2,col = "dark red",line=0)

x1<--0.5
y<-0.5
x2<-max(dias+inicio)/5-1
x3<-2*max(dias+inicio)/5-2
legend(x1,y,c("Start"),bty="n",cex=0.6,pch=19,col = "dark blue")
legend(x2,y,c("Event"),bty="n",cex=0.6,pch=4,col = colevent)
legend(x3,y,c("Censored"),bty="n",cex=0.6,pch=0,col = colcensor)
for(i in 1:n) {
temp1<-censored[i]
if (temp1==0) temp1<-0 else temp1<-4
segments(inicio[i],paciente[i],dias[i],paciente[i],col="black", lty = "dotted")
}
a<-1
if(a==1){m<-nrow(especiales)
for(j in 1:m){temp2<-especiales[j,3]
if(temp2==1)temp2<-4 else {if(temp2==0)temp2<-0 else temp2<-0}
temp3<-0
for(i in 1:n){if (especiales[j,1]<=n) temp3[especiales[j,1]==i]<-inicio[i]}
if (temp2==4) {colorx1<-colevent;us<-"\u58"} else {colorx1<-colcensor;us<-"\u4F"}
if (especiales[j,1]<=n)
points((temp3+especiales[j,2]),especiales[j,1],pch=us,col=colorx1,cex=0.5)
}}}
