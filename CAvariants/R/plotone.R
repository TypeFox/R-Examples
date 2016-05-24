plotone <-
function(a1,a2,plottype,things,nthings,nvars,Thingcoord,Varcoord,
                       inertiapc,
                       thinggroup,thinggrlab,vargroup,vargrlab,thinglabels,varlabels,picsize,cex=.8,type , catype ,pos=2) {

plot.new()
plot(Thingcoord[ ,a1], Thingcoord[ ,a2], xlim=picsize, ylim=picsize, 
     xlab=paste("Axis ", a1, "    ", inertiapc[a1], "%", sep=""), 
     ylab=paste("Axis ", a2, "    ", inertiapc[a2], "%", sep=""),  
          asp=1, pch=thinggrlab[[3]][as.integer(thinggroup[[2]])], col=thinggrlab[[4]][as.integer(thinggroup[[2]])],cex=cex,type=type )
abline(h=0,v=0)
if (as.integer(max(thinggroup[[2]]))>1) legend("topleft",thinggrlab[[2]],pch=thinggrlab[[3]])

if ((plottype=="Biplot")|(plottype=="biplot")|(plottype=="bip")|(plottype=="b")) { 
  title(paste(things,"Biplot"))
  for (i in 1:nthings) { 
if ((catype=="NSCA")|(catype=="CA"))
{
lines(c(0,Thingcoord[i,a1]), c(0,Thingcoord[i,a2]), col=thinggrlab[[4]][as.integer(thinggroup[[2]][i])])}
#points(Thingcoord[i,a1],Thingcoord[i,a2],pch=thinggrlab[[3]][as.integer(thinggroup[[2]])],col=thinggrlab[[4]][as.integer(thinggroup[[2]])]) 
} 
 grat <- cbind(Thingcoord[,a1]/picsize[1],Thingcoord[,a1]/picsize[2],Thingcoord[,a2]/picsize[1],Thingcoord[,a2]/picsize[2],0.95)/0.95
 cl <- 1.05/apply(grat,1,max)
  text(cl*Thingcoord[ ,a1], cl*Thingcoord[ ,a2], labels=thinglabels, col=thinggrlab[[4]][as.integer(thinggroup[[2]])], pos=pos, cex=0.75 )

grat2 <- cbind(Varcoord[,a1]/picsize[1],Varcoord[,a1]/picsize[2],Varcoord[,a2]/picsize[1],Varcoord[,a2]/picsize[2],0.95)/0.95
  cl2 <- 1.05/apply(grat2,1,max)
points(Varcoord[ ,a1], Varcoord[ ,a2], asp=1,pch=vargrlab[[3]][as.integer(vargroup[[2]])], col=vargrlab[[4]][as.integer(vargroup[[2]])] )
text(cl2*Varcoord[ ,a1], cl2*Varcoord[ ,a2], labels=varlabels, col=vargrlab[[4]][as.integer(vargroup[[2]])], pos=pos, cex=0.75 )

} else { # french
  title(paste(things,"Plot"))
  points(Varcoord[ ,a1], Varcoord[ ,a2], asp=1, 
         pch=vargrlab[[3]][as.integer(vargroup[[2]])], col=vargrlab[[4]][as.integer(vargroup[[2]])] )
text(Varcoord[ ,a1], Varcoord[ ,a2], labels=varlabels, col=vargrlab[[4]][as.integer(vargroup[[2]])], pos=pos, cex=cex )
text(Thingcoord[ ,a1], Thingcoord[ ,a2], labels=thinglabels, col=thinggrlab[[4]][as.integer(thinggroup[[2]])], pos=pos, cex=cex )
if (as.integer(max(vargroup[[2]]))>1) legend("topright",vargrlab[[2]],pch=vargrlab[[3]]) 

}


}
