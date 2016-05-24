plotvecs<-function(vecs,
depths=NULL,segme=T,lift=NULL,
modetest=NULL,alpha=NULL,
axes=TRUE,xlim=NULL,ylim=NULL,xaxt=xaxt,col="black",col.axis="black",
modecolors=NULL,modethickness=1,
leafcolors=NULL,leaflift=0,leafsymbo=20,
modelabels=NULL,ptext=0,
yaxt="s",log="",cex.axis=1)
{
#Plots vectors in vec
#
#vecs is nodenum*4-matrix
#vecs[i,1] x-coordi alulle
#vecs[i,2] y-coordi alulle = vecs[i,4] y-coordi lopulle
#vecs[i,3] x-coordi lopulle
#
#plot(c(1,2),c(3,3))  
#segments(1,3,2,3)
#
#plot(c(1,2,3,4),c(3,3,2,2))  
#segments(1,3,2,3)
#segments(3,2,4,2)
#
#vecs<-matrix(0,3,4)    
#vecs[1,]<-c(1,1,4,1)
#vecs[2,]<-c(5,1,6,1)
#vecs[3,]<-c(2,2,3,2)
#
#plot(c(1,4,5,6,2,3),c(1,1,1,1,2,2))
#segments(1,1,4,1)
#segments(5,1,6,1)
#segments(2,2,3,2)

nodenum<-length(vecs[,1])
xcoor<-matrix(0,2*nodenum,1)
ycoor<-matrix(0,2*nodenum,1)

for (i in 1:nodenum){
 xcoor[2*i-1]<-vecs[i,1]
 xcoor[2*i]<-vecs[i,3]
 ycoor[2*i-1]<-vecs[i,2]
 ycoor[2*i]<-vecs[i,4]
}

#ylim<-c(0,max(ycoor)+ptext)
plot(xcoor,ycoor,xlab="",ylab="",axes=axes,xlim=xlim,ylim=ylim,xaxt=xaxt,
col=col,col.axis=col.axis,yaxt=yaxt,log=log,cex.axis=cex.axis)

if (!is.null(leafcolors)){
   xpoint<-matrix(0,nodenum,1)
   ypoint<-matrix(0,nodenum,1)
   leafcol<-matrix("",nodenum,1)
   zahl<-0
   for (no in 1:nodenum){
      if (leafcolors[no]!="black"){
          zahl<-zahl+1
          xpoint[zahl]<-vecs[no,1]+(vecs[no,3]-vecs[no,1])/2
         
          lif<-(depths[no]-1)*lift
          yc<-ycoor[2*no-1]+lif
          ypoint[zahl]<-yc+leaflift
         
          leafcol[zahl]<-leafcolors[no]
      }
   }
   xpoint<-xpoint[1:zahl]
   ypoint<-ypoint[1:zahl]
   leafcol<-leafcol[1:zahl]
   points(xpoint,ypoint,pch=leafsymbo,col=leafcol)
}

if (!is.null(modelabels)){
   for (no in 1:nodenum){
      if (modelabels[no]!=""){
          xpoint<-vecs[no,1]+(vecs[no,3]-vecs[no,1])/2
         
          lif<-(depths[no]-1)*lift
          yc<-ycoor[2*no-1]+lif
          ypoint<-yc+ptext
         
          label<-modelabels[no]
          
          text(xpoint,ypoint,label)

      }
   }
}

if (segme==T){
 
   thick<-1
   lif<-0
   col<-"black"   

   for (i in 1:nodenum){

        if (!is.null(depths))  lif<-(depths[i]-1)*lift
        if (!is.null(modecolors)){
                if (modecolors[i]!="black") thick=modethickness 
                                            #thick<-2.2^(depths[i]-1)  
                col<-modecolors[i]   
        }
        if (!is.null(modetest)){
             col<-4
             if (modetest[i]>0){
                if (modetest[i]>alpha)  col<-2   
                     #red, hyvaksytaan 0-hypoteesi=ei moodia
                     #nodes with red are not a real feature
                else col<-4   #blue
             } 
       }
           #testcrit<-modetest[i]*qnorm(1-alpha/2)
           #if (excmassa>testcrit)

        yc<-ycoor[2*i-1]+lif
        segments(xcoor[2*i-1],yc,xcoor[2*i],yc,col=col,lwd=thick)
        
        #lines(c(xcoor[2*i-1],xcoor[2*i]),c(ycoor[2*i-1],ycoor[2*i]),col=2) 
  
   }
}

#return(t(tc),t(em))
}











