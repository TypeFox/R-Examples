anm.ExpDesign<-function(method= "all", titles =  TRUE, cex.text = 1, mp.col = NULL, interval = 0.5, iter = 30){
old.par <- par(no.readonly = TRUE)
    if(any(method == "all")) method = c("CRD","factorial2by2","factorial2by2by2","nested","RCBD","RIBD","split","split.split", "SPRB","strip","split.block","strip.split","latin","pairs") 
    m<-matrix(nrow=iter,ncol=length(method),rep(method,iter),byrow=TRUE)
    for(i in 1:nrow(m)){
	    dev.hold()
      ExpDesign(method=m[i,],titles=titles, cex.text = cex.text, mp.col = mp.col)
	    dev.flush()
	    Sys.sleep(interval)
    }
on.exit(par(old.par))
invisible()	
}

ExpDesign<-function(method= "all",titles=TRUE, cex.text = 1, mp.col = NULL,...){
if(any(method == "all")) method = c("CRD","factorial2by2","factorial2by2by2","nested","RCBD","RIBD","split","split.split", "SPRB","strip","split.block","strip.split","latin","pairs") 
L<-length(method)
if(L==2) {par(mfrow=c(2,1),mar=c(0.1,1.5,2,1.5))} else
if(L==3) {par(mfrow=c(3,1),mar=c(0.1,1.5,2,1.5))} else
if(L==4) {par(mfrow=c(2,2),mar=c(0.1,1.5,2,1.5))} else
if(L==5|L==6) {par(mfrow=c(3,2),mar=c(0.1,1.5,2,1.5))} else
if(L==7|L==8|L==9) {par(mfrow=c(3,3),mar=c(0.1,1.0,2,1.0))}else
if(L==10|L==11|L==12) {par(mfrow=c(4,3),mar=c(0.1,1.0,2,1.0))}else
if(L==13|L==14) {par(mfrow=c(5,3),mar=c(0.1,1.0,1.5,1.0))}

#CRD
if(any(method=="CRD")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"CRD",""),...)
text(c(3,5.5,8),rep(9,3),c("EU1","EU2","EU3"),cex=cex.text)
segments(5.5,0,5.5,7.5,col="gray")
segments(0,7.5,11,7.5,col="gray")
text(c(3,8),c(4,4),cex=1.2*cex.text,c(expression(A[1]),expression(A[2])))
x.s<-sample(c(3,5.5,8),3,replace=FALSE)
x.e<-sample(c(3,8),3,replace=TRUE)
x.e[3]<-ifelse(x.e[1]==x.e[2]&x.e[1]==3,8,x.e[3])
x.e[3]<-ifelse(x.e[1]==x.e[2]&x.e[1]==8,3,x.e[3])
arrows(x.s,rep(8.5,3),x.e,rep(4.5,3),length=.1)
  }

#Factorial 2 by 2
if(any(method=="factorial2by2")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"2x2 Factorial",""),...)
segments(5.5,0,5.5,10.5,col="gray")
segments(0,5.5,10.5,5.5,col="gray")
text(x=c(1.45,6.4,1.45,6.4),y=c(9.75,4.75,4.75,9.75),c("EU1","EU2","EU3","EU4"),cex=cex.text)
m<-matrix(nrow=4,ncol=2,c(2.75,2.75,7.75,2.75,7.75,7.75,2.75,7.75),byrow=TRUE)
s<-sample(seq(1,4),4,replace=FALSE)
m.s<-m[s,]
text(x=m.s[,1],y=m.s[,2],c(expression(paste(A[1],B[1])),expression(paste(A[2],B[1])),expression(paste(A[1],B[2])),expression(paste(A[2],B[2]))),cex=1.1*cex.text)
  }

#Factorial 2 by 2 by 2
if(any(method=="factorial2by2by2")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"2x2x2 Factorial",""),...)
segments(4,0,4,11,col="gray")
segments(7.25,0,7.25,11,col="gray")
segments(0,4,11,4,col="gray")
segments(0,7.25,11,7.25,col="gray")
text(x=c(rep(c(1.32,4.65,7.92),2),c(1.32,4.65)),y=c(rep(9.9,3),rep(6.8,3),rep(3.55,2)), c("EU1","EU2","EU3","EU4","EU5","EU6","EU7","EU8"),cex=cex.text)
t<-c(expression(paste(A[1],B[1],C[1])),expression(paste(A[2],B[1],C[1])),expression(paste(A[1],B[2],C[1])),expression(paste(A[2],B[2],C[1])),
expression(paste(A[1],B[1],C[2])),expression(paste(A[2],B[1],C[2])),expression(paste(A[1],B[2],C[2])),expression(paste(A[2],B[2],C[2])))
s<-sample(seq(1,8),8,replace=FALSE)
x<-c(rep(c(2.25,5.75,8.75),2),c(2.25,5.75))
y<-c(8.75,8.75,8.75,5.75,5.75,5.75,2.25,2.25)
text(x,y,t[s],cex=cex.text)
  }

#Nested
if(any(method=="nested")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"Nested",""),...)
x0<-c(3,3,8,8,1.75,1.75,4.25,4.25,6.75,6.75,9.25,9.25)
x1<-c(1.75,4.25,6.75,9.25,1.25,2.25,3.75,4.75,6.25,7.25,8.75,9.75)
y0<-c(7.75,7.75,7.75,7.75,5.5,5.5,5.5,5.5,5.5,5.5,5.5,5.5)
y1<-c(5.5,5.5,5.5,5.5,2.75,2.75,2.75,2.75,2.75,2.75,2.75,2.75)
segments(x0,y0,x1,y1)
symbols(x=c(5.5,5.5,5.5),y=c(7.75,5.5,2.75),rectangles=matrix(ncol=2,nrow=3,data=c(9,9,9,.75,.75,.75)),add=TRUE,
inches=FALSE,fg="white",bg="white")
segments(5.5,0,5.5,10.5,col="gray")
a<-c(expression(A[1]),expression(A[2]))
b<-c(expression(B[1]),expression(B[2]))
c<-c(expression(C[1]),expression(C[2]))
s<-sample(seq(1,2),2,replace=FALSE);s1<-sample(seq(1,2),2,replace=FALSE);s2<-sample(seq(1,2),2,replace=FALSE)
s3<-sample(seq(1,2),2,replace=FALSE);s4<-sample(seq(1,2),2,replace=FALSE);s5<-sample(seq(1,2),2,replace=FALSE)
s6<-sample(seq(1,2),2,replace=FALSE)
A<-a[s];B1<-b[s1];B2<-b[s2];C1<-c[s3];C2<-c[s4];C3<-c[s5];C4<-c[s6]
text(x=c(3,8,1.75,4.25,6.75,9.25,1.25,2.25,3.75,4.75,6.25,7.25,8.75,9.75),y=c(7.75,7.75,5.5,5.5,5.5,5.5,2.75,2.75,2.75,2.75,2.75,2.75,2.75,2.75),c(A[1],
A[2],B1[1],B1[2],B2[1],B2[2],C1[1],C1[2],C2[1],C2[2],C3[1],C3[2],C4[1],C4[2]),cex=cex.text)
  }

#RCBD
if(any(method=="RCBD")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"RCBD",""),...)
text(c(2,4.3,6.7,9.2),rep(9,3),c("EU1","EU2","EU3","EU4"),cex=cex.text)
segments(c(5.5,0,3,8),c(0,7.5,0,0),c(5.5,11,3,8),c(11,7.5,7.5,7.5),col="gray")
text(c(3,8),c(9.9,9.9),c("Block 1","Block 2"),cex=.9*cex.text,col="gray")
text(c(2,4.3,6.7,9.2),c(4,4,4,4),cex=1.1*cex.text,c(expression(A[1]),expression(A[2]),expression(A[1]),expression(A[2])))
x.s1<-sample(c(2,4.3),2,replace=FALSE)
x.s2<-sample(c(6.7,9.2),2,replace=FALSE)
x.e1<-sample(c(2,4.3),2,replace=FALSE)
x.e2<-sample(c(6.7,9.2),2,replace=FALSE)
arrows(x.s1,rep(8.5,2),x.e1,rep(4.5,2),length=.1)
arrows(x.s2,rep(8.5,2),x.e2,rep(4.5,2),length=.1)
  }

#RIBD
if(any(method=="RIBD")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"RIBD",""),...)
segments(4,0,4,11,col="gray")
segments(7.25,0,7.25,11,col="gray")
segments(0,7.5,11,7.5,col="gray")
segments(c(2.3,5.7,8.8),c(0,0,0),c(2.3,5.7,8.8),c(7.5,7.5,7.5),col="gray")
text(c(2.3,5.7,8.8),c(9.9,9.9,9.9),c("Block 1","Block 2","Block 3"),cex=.9*cex.text,col="gray")
text(c(1.4,3.15,4.85,6.475,8.025,9.55),c(rep(4,6)),c(expression(A[1]),expression(A[3]),expression(A[2]),expression(A[3]),expression(A[1]),expression(A[2])),cex=cex.text)
text(c(1.4,3.15,4.85,6.475,8.025,9.55),c(rep(9,6)),c("EU1","EU2","EU3","EU4","EU5","EU6"),cex=cex.text)
x.s1<-sample(c(1.4,3.15),2,replace=FALSE)
x.s2<-sample(c(4.85,6.475),2,replace=FALSE)
x.s3<-sample(c(8.025,9.55),2,replace=FALSE)
x.e1<-sample(c(1.4,3.15),2,replace=FALSE)
x.e2<-sample(c(4.85,6.475),2,replace=FALSE)
x.e3<-sample(c(8.025,9.55),2,replace=FALSE)
arrows(x.s1,rep(8.5,2),x.e1,rep(4.5,2),length=.1)
arrows(x.s2,rep(8.5,2),x.e2,rep(4.5,2),length=.1)
arrows(x.s3,rep(8.5,2),x.e3,rep(4.5,2),length=.1)
  }
                                                                  
#Split plot
if(any(method=="split")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"Split Plot",""),...)
segments(5.5,0,5.5,10.5,col="gray")
segments(0,5.5,10.5,5.5,col="gray")
symbols(x=c(3,8),y=c(5.5,5.5),rectangles=matrix(nrow=2,ncol=2,c(1.5,1.5,.75,.75)),add=TRUE,
inches=FALSE,fg="white",bg="white")
a<-c(expression(A[1]),expression(A[2]))
s<-sample(seq(1,2),2,replace=FALSE)
text(x=c(3,8),y=c(5.5,5.5),a[s],cex=1.1*cex.text)
s<-sample(seq(1,2),2,replace=FALSE)
m<-c(expression(B[1]),expression(B[2]))
text(x=c(3,3),y=c(3,8),m[s],cex=1.1*cex.text)
text(x=c(8,8),y=c(3,8),m[s],cex=1.1*cex.text)
  }

#Split split plot
if(any(method=="split.split")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"Split Split Plot",""),...)
sym.y<-c(rep(3,2),rep(5.5,2),rep(8,2))
segments(c(5.5,3,8),c(0,0,0),c(5.5,3,8),c(11,11,11),col=c("black","gray","gray"))
segments(0,5.5,11,5.5,col=1)
symbols(x=rep(c(3,8),3),y=sym.y,rectangles=matrix(nrow=6,ncol=2,rep(1,12)),add=TRUE,inches=FALSE,fg="white",bg="white")
a<-c(expression(A[1]),expression(A[2]))
b<-c(expression(B[1]),expression(B[2]))
c<-c(expression(C[1]),expression(C[2]))
s<-sample(seq(1,2),2,replace=FALSE);s1<-sample(seq(1,2),2,replace=FALSE);s2<-sample(seq(1,2),2,replace=FALSE)
s3<-sample(seq(1,2),2,replace=FALSE);s4<-sample(seq(1,2),2,replace=FALSE);s5<-sample(seq(1,2),2,replace=FALSE)
s6<-sample(seq(1,2),2,replace=FALSE)
A<-a[s];B1<-b[s1];B2<-b[s2];C1<-c[s3];C2<-c[s4];C3<-c[s5];C4<-c[s6]
text(x=c(3,8,3,3,8,8,1.75,4.25,1.75,4.25,6.75,9.25,6.75,9.25),y=c(5.5,5.5,8,3,8,3,8,8,3,3,8,8,3,3),c(A[1],
A[2],B1[1],B1[2],B2[1],B2[2],C1[1],C1[2],C2[1],C2[2],C3[1],C3[2],C4[1],C4[2]),cex=1.1*cex.text)
  }

#Split plot in completely randomized blocks
if(any(method=="SPRB")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"Split Plot Randomized Block",""),...)
segments(c(5.5,3,8),c(0,0,0),c(5.5,3,8),c(11,11,11),col=c("black","gray","gray"),lty=c(2,1,1))
segments(c(0,0,0),c(5.5,3,8),c(11,11,11),c(5.5,3,8),col="gray")
sym.y<-c(rep( 1.85,2),rep(4.3,2),rep(6.7,2),rep(9.2,2))
symbols(x=rep(c(3,8),4),y=sym.y,rectangles=matrix(nrow=8,ncol=2,rep(1,16)),add=TRUE,inches=FALSE,fg="white",bg="white")
split.x<-c(1.85,4.3,6.7,9.2)
split.y<-matrix(ncol=4,nrow=4, data=c(rep( 1.85,4),rep(4.3,4),rep(6.7,4),rep(9.2,4)),byrow=TRUE)
split.txt<-c(expression(B[1]),expression(B[2]))
whole.txt<-c(expression(A[1]),expression(A[2]))
for(i in 1:4){
s1<-sample(seq(1,2),2,replace=FALSE);s2<-sample(seq(1,2),2,replace=FALSE);s3<-sample(seq(1,2),2,replace=FALSE)
text(c(1.9,4.3),y=split.y[i,],split.txt[s1],cex=1.1*cex.text)
text(c(6.7,9.2),y=split.y[i,],split.txt[s2],cex=1.1*cex.text)
text(c(3,8),y=matrix(nrow=4,ncol=2,sym.y,byrow=TRUE)[i,],whole.txt[s3],cex=1.1*cex.text)
  }
mtext(c("Block 1","Block 2","Block 3","Block 4"), side = 2,at=c(1.85,4.3,6.7,9.2),cex=ifelse(L>5,1/(L*0.2*cex.text),0.7*cex.text))
  }

#Strip plot
if(any(method=="strip")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"Strip Plot",""),...)
segments(5.5,0,5.5,10.5,col="gray")
segments(0,5.5,10.5,5.5,col="gray")
symbols(x=c(3,8),y=c(5.5,5.5),rectangles=matrix(nrow=2,ncol=2,c(1.5,1.5,.75,.75)),add=TRUE,
inches=FALSE,fg="white",bg="white")
symbols(x=c(5.5,5.5),y=c(3,8),rectangles=matrix(nrow=2,ncol=2,c(.75,.75,1.5,1.5)),add=TRUE,
inches=FALSE,fg="white",bg="white")
a<-c(expression(A[1]),expression(A[2]))
s<-sample(seq(1,2),2,replace=FALSE)
text(x=c(3,8),y=c(5.5,5.5),a[s],cex=1.1*cex.text)
s<-sample(seq(1,2),2,replace=FALSE)
m<-c(expression(B[1]),expression(B[2]))
text(y=c(3,8),x=c(5.5,5.5),m[s],cex=1.1*cex.text)
  }

#Split block/Strip split plot (2 way) # Littell et al. Mixed Models in SAS
if(any(method=="split.block")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"Split Block",""),...)
segments(c(5.5,3,8),c(0,0,0),c(5.5,3,8),c(11,11,11),col=c("black","gray","gray"))
segments(c(0,0,0),c(5.5,3,8),c(11,11,11),c(5.5,3,8),col=c("black","gray","gray"))
sym.y<-c(rep(1.85,2),rep(4.3,2),rep(6.7,2),rep(9.2,2))
symbols(x=rep(c(3,8),4),y=sym.y,rectangles=matrix(nrow=8,ncol=2,rep(1,16)),add=TRUE,inches=FALSE,fg="white",bg="white")
split.y<-matrix(ncol=4,nrow=2,data=c(rep(1.85,2),rep(4.3,2),rep(6.7,2),rep(9.2,2)))
split.y2<-matrix(ncol=4,nrow=2,data=c(rep(4.3,2),rep(1.85,2),rep(9.2,2),rep(6.7,2)))
split.x<-matrix(ncol=4,nrow=2, data=c(c(1.85,4.3),c(6.7,9.2),c(1.85,4.3),c(6.7,9.2)))
whole.x<-matrix(ncol=2,nrow=4,rep(c(rep(3,2),rep(8,2)),2),byrow=FALSE)
whole.y<-matrix(ncol=2,nrow=4,c(rep(c(1.85,4.3,6.7,9.2),2)),byrow=TRUE)
split.txt<-c(expression(B[1]),expression(B[2]))
whole.txt<-c(expression(A[1]),expression(A[2]))
for(i in 1:4){
s1<-sample(seq(1,2),2,replace=FALSE);s2<-sample(seq(1,2),2,replace=FALSE)
s3<-sample(seq(1,2),2,replace=FALSE)
text(x=split.x[,i],y=split.y[,i],split.txt[s1],cex=1.1*cex.text)
text(x=split.x[,i],y=split.y2[,i],split.txt[s1],cex=1.1*cex.text)
text(whole.x[i,],y=whole.y[i,],whole.txt[s3],cex=1.1*cex.text)
  }
mtext(c("Block 1","Block 2"), side = 2,at=c(3,8),cex=ifelse(L>5,1/(L*0.2*cex.text),0.7*cex.text))
mtext(c("Block 3","Block 4"), side = 4,at=c(3,8),cex=ifelse(L>5,1/(L*0.2*cex.text),0.7*cex.text),line=0)
  }


#Strip split plot (3 way) # Milliken et al. Analysis of Messy Data Vol. 1
if(any(method=="strip.split")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"Strip Split Plot (3 way)",""),...)
sym.y<-c(rep(3,2),rep(5.5,2),rep(8,2))
segments(c(5.5,3,8),c(0,0,0),c(5.5,3,8),c(11,11,11),col=c("black","gray","gray"))
segments(0,5.5,11,5.5,col=1)
symbols(x=rep(c(3,8),3),y=sym.y,rectangles=matrix(nrow=6,ncol=2,rep(1,12)),add=TRUE,inches=FALSE,fg="white",bg="white")
a<-c(expression(A[1]),expression(A[2]))
b<-c(expression(B[1]),expression(B[2]))
c<-c(expression(C[1]),expression(C[2]))
s<-sample(seq(1,2),2,replace=FALSE);s1<-sample(seq(1,2),2,replace=FALSE);s2<-sample(seq(1,2),2,replace=FALSE)
s3<-sample(seq(1,2),2,replace=FALSE);s4<-sample(seq(1,2),2,replace=FALSE);s5<-sample(seq(1,2),2,replace=FALSE)
s6<-sample(seq(1,2),2,replace=FALSE)
A<-a[s];B1<-b[s1];B2<-b[s1];C1<-c[s3];C2<-c[s4];C3<-c[s5];C4<-c[s6]
text(x=c(3,8,3,3,8,8,1.75,4.25,1.75,4.25,6.75,9.25,6.75,9.25),y=c(5.5,5.5,8,3,8,3,8,8,3,3,8,8,3,3),c(A[1],
A[2],B1[1],B1[2],B2[1],B2[2],C1[1],C1[2],C2[1],C2[2],C3[1],C3[2],C4[1],C4[2]),cex=1.1*cex.text)
  }

#Latin squares
if(any(method=="latin")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"Latin Squares",""),...)
segments(4,0,4,11,col="gray")
segments(7.25,0,7.25,11,col="gray")
segments(0,4,11,4,col="gray")
segments(0,7.25,11,7.25,col="gray")
m1<-c("A","B","C","B","C","A","C","A","B")
m2<-c("A","C","B","B","A","C","C","B","A")
m3<-c("B","A","C","C","B","A","A","C","B")
m4<-c("C","B","A","A","C","B","B","A","C")
m<-rbind(m1,m2,m3,m4)
s<-sample(seq(1,4),1)
x<-rep(c(2.25,5.75,8.75),3)
y<-c(8.75,8.75,8.75,5.75,5.75,5.75,2.25,2.25,2.25)
text(x,y,m[s,],cex=cex.text)
  }

#Matched Pairs Design
if(any(method=="pairs")){
plot(seq(1:10),seq(1:10),xaxt="n",yaxt="n",xlab="",ylab="",type="n",main=ifelse(titles==TRUE,"Matched Pairs",""),...)
text(c(3,5.5,8),rep(9,3),c("EU1","EU2","EU3"),col=ifelse(is.null(mp.col),c(1,2,3),mp.col),cex=cex.text)
segments(5.5,0,5.5,7.5,col="gray")
segments(0,7.5,11,7.5,col="gray")
text(c(3,8),c(4,4),cex=1.2*cex.text,c(expression(A[1]),expression(A[2])))
x.s<-sample(c(3,5.5,8),3,replace=FALSE)
x.e<-sample(c(3,8),3,replace=TRUE)
x.e[3]<-ifelse(x.e[1]==x.e[2]&x.e[1]==3,8,x.e[3])
x.e[3]<-ifelse(x.e[1]==x.e[2]&x.e[1]==8,3,x.e[3])
col<-matrix(ncol=1,nrow=3)
for(i in 1:3){
if(x.s[i]==3){col[i]=2} 
if(x.s[i]==5.5){col[i]=3} 
if(x.s[i]==8){col[i]=4}}
if(is.null(mp.col))col=col
else col=mp.col
arrows(x.s,rep(8.5,3),x.e,rep(4.5,3),length=.1,col=col,lty=c(1,2,3))
x.e2<-matrix(ncol=1,nrow=3)
y.s2<-matrix(ncol=1,nrow=3)
for(i in 1:3){
x.e2[i]<-ifelse(x.e[i]==3,8,3) 
if(x.s[i]==3){y.s2[i]=3.2}
if(x.s[i]==5.5){y.s2[i]=2.6}
if(x.s[i]==8){y.s2[i]=2.0}}
arrows(x.e,y.s2,x.e2,y.s2,length=.1,col=col,lty=c(1,2,3))
  }
  }
