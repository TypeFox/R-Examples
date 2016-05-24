persistence.plot <-
function(fit, preds, PI, title)
{
 if (class(fit)!="logforest"  & class(fit)!="LBoost") 
    stop("fit not of class logforest or LBoost")
 tot.PI<-sum(fit$PI.frequency)
 persist.info<-persistence.prep(fit=fit, preds=preds, PI=PI)
 pert.PI<-sum(persist.info$all.freqs)
 percent.pert<-round((pert.PI/tot.PI)*100,digits=1)
 leg<-paste(percent.pert, "% of all", sep="")
 main.sz<-persist.info$main.freq
 m.sz<-max(as.vector(persist.info$all.freqs))
 sizes<-sort(persist.info$sizes)
 rad.mat<-persist.info$rad.mat
 prime.info<-persist.info$primary.info
 primexy<-persist.info$primary.info
 pos.list<-persist.info$position.list
 par(mai= c(0.1,0.1,0.75,0.1))
 grp<-radial.plot(rad.mat[,1],rad.mat[,2],main=title, labels="",
                  line.col=1,lwd=1, radial.lim=c(0,5),show.radial.grid=FALSE,show.grid=FALSE)
 draw.circle(0,0,main.sz*1/m.sz,border="red", col="red")
 text(0,0,PI, cex=0.9, col="black") 
 for (i in 1:length(sizes))
   {
   draw.circle(0,0,sizes[i],border="grey")
   }
 if (length(pos.list)>1)
   {
   for (j in length(pos.list):2)
     {
     pos.mat<-pos.list[[j]]
     if (length(pos.list[[j]])>0)
       {
       for (k in 1:nrow(pos.list[[j]]))
         {
         x1<-as.numeric(pos.mat[k,1]); y1<-as.numeric(pos.mat[k,2])
         x2<-as.numeric(pos.mat[k,3]); y2<-as.numeric(pos.mat[k,4])
         segments(x1,y1,x2,y2,col=1, lwd=1)
         }
       }
     }
   for (j in length(pos.list):2)
     { 
     pos.mat<-pos.list[[j]]
     nms<-rownames(pos.mat)
     occ.tab<-table(nms)
     if (length(pos.list[[j]])>0)
       {
       for (k in 1:nrow(pos.list[[j]]))
         {
         nm<-nms[k]
         x<-as.numeric(pos.mat[k,3]); y<-as.numeric(pos.mat[k,4])
         freq<-as.numeric(pos.mat[k,5])
         circ.sz<-freq*0.5/m.sz
         draw.circle(x,y,circ.sz,col="lightblue1",border="lightblue3")
   if (x>=0) {srt.val<-deg(atan2(y,x))}
   if (x<0 & y>=0) {srt.val<-deg(atan2(y,x))-180}
   if (x<0 & y<0) {srt.val<-180+deg(atan2(y,x))}
         text(x,y,pos.mat[k,6],cex=0.65, col="black", srt=srt.val)
         }
       }
     }
   }
 for(i in 1:nrow(primexy))
   {
   x<-as.numeric(primexy[i,1]); y<-as.numeric(primexy[i,2])
   circ.sz<-as.numeric(primexy[i,3])*0.5/m.sz
   draw.circle(x,y,circ.sz, col="lightblue1", border="lightblue3")
   if (x>=0) {srt.val<-deg(atan2(y,x))}
   if (x<0 & y>=0) {srt.val<-deg(atan2(y,x))-180}
   if (x<0 & y<0) {srt.val<-180+deg(atan2(y,x))}
   text(x,y,primexy[i,4],cex=0.65,col="black", srt=srt.val)
   }
 legend(x=3.25, y=-4.9, c(leg, "interactions"), cex=0.75)
}
