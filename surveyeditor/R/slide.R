slide <-
function(id,type,title,content,h,v,size,col,time=NULL,interval=c(-Inf,Inf),isinteger=F,choices=NULL,layout=NULL,n.select=NULL,box.level=NULL,cond=NULL,scale=NULL,grid=NULL,range=NULL,unit=NULL,initial=NULL,increment=NULL,order=c("ascending","descending"),Extra=NULL,Line.draw=NULL){

Annex<-function(A){
rescale<-function(data,data.ends,ends){
new.data<-NULL;for(i in 1:length(data)){
new.data<-c(new.data,(data[i]-data.ends[1])*(ends[2]-ends[1])/(data.ends[2]-data.ends[1])+ends[1])}
return(new.data)}
if(A$type!="t"){
rect(A$loc[1],A$loc[2],A$loc[3],A$loc[4])
if(A$type=="h1"||A$type=="h2"){
Tab<-table(A[[1]]);N<-length(rownames(Tab))
if(A$type=="h2"){Tab<-Tab/sum(Tab)}
barl<-(A$loc[3]-A$loc[1]-(N+1))/N
xL.list<-(0:(N-1))*barl+(1:N)+A$loc[1];xR.list<-(1:N)*(barl+1)+A$loc[1]
yB.list<-rep(A$loc[2],N);yT.list<-rescale(Tab,c(0,max(Tab)),c(A$loc[2],A$loc[4]-2))
for(i in 1:N){rect(xL.list[i],yB.list[i],xR.list[i],yT.list[i],col=A$col[i])}
ticks1<-A$ticks[[1]]
mk<-rescale(ticks1[ticks1>=-0.001 & (max(Tab)+0.001)>=ticks1],c(0,max(Tab)),c(A$loc[2],A$loc[4]-2))
text(A$loc[1]-1,mk,adj=c(1,0.5),cex=1,labels=paste(ticks1[ticks1>=0 & (max(Tab)+0.001)>=ticks1]))
for(i in 1:length(mk)){lines(c(A$loc[1]-1,A$loc[1]),rep(mk[i],2))}
text((xL.list+xR.list)/2,A$loc[2]-2,adj=c(0.5,0.5),cex=1,labels=paste(rownames(Tab)))
}
if(A$type=="p"||A$type=="l"||A$type=="b"){
if(is.numeric(A[[1]])!=T){X<-seq(1,length(A[[1]]),1)}else{X<-A[[1]]}
Y<-A[[2]]
Xs<-rescale(X,c(min(X),max(X)),c(A$loc[1]+2,A$loc[3]-2));Ys<-rescale(Y,c(min(Y),max(Y)),c(A$loc[2]+2,A$loc[4]-2))
if(A$type=="l"||A$type=="b"){lines(Xs,Ys)}
if(A$type=="p"||A$type=="b"){points(Xs,Ys)}
ticks1<-A$ticks[[1]]
if(is.numeric(A[[1]])!=T){ticks1<-round(ticks1,0)}
ticks2<-A$ticks[[2]]
mk1<-rescale(ticks1[ticks1>=(min(X)-0.001) & ticks1<=(max(X)+0.001)],c(min(X),max(X)),c(A$loc[1]+2,A$loc[3]-2))
mk2<-rescale(ticks2[ticks2>=(min(Y)-0.001) & ticks2<=(max(Y)+0.001)],c(min(Y),max(Y)),c(A$loc[2]+2,A$loc[4]-2))
for(i in 1:length(mk1)){lines(rep(mk1[i],2),c(A$loc[2]-1,A$loc[2]))}
for(i in 1:length(mk2)){lines(c(A$loc[1]-1,A$loc[1]),rep(mk2[i],2))}
if(is.numeric(A[[1]])==T){text(mk1,A$loc[2]-2,adj=c(0.5,0.5),cex=1,labels=paste(A[[1]][ticks1[ticks1>=(min(X)-0.001) & ticks1<=(max(X)+0.001)]]))}
if(is.numeric(A[[1]])!=T){text(mk1,A$loc[2]-2,adj=c(0.5,0.5),cex=1,labels=paste(A[[1]][ticks1[ticks1>=(min(X)-0.001) & ticks1<=(max(X)+0.001)]]))}
text(A$loc[1]-2,mk2,adj=c(1,0.5),cex=1,labels=paste(ticks2[ticks2>=(min(Y)-0.001) & ticks2<=(max(Y)+0.001)]))}
if(is.null(A$xlabl)==T){xlabl<-names(A)[1]}else{xlabl<-A$xlabl}
text((A$loc[1]+A$loc[3])/2,A$loc[2]-4,adj=c(0.5,1),cex=1,labels=xlabl)
if(is.null(A$ylabl)==T){ylabl<-names(A)[2]}else{ylabl<-A$ylabl}
text(A$loc[1]-5,(A$loc[2]+A$loc[4])/2,adj=c(1,0.5),cex=1,labels=ylabl)
titlel<-NULL
if(is.null(A$title)==F){titlel<-A$title}
text((A$loc[1]+A$loc[3])/2,A$loc[4]+2,adj=c(0.5,0.5),labels=titlel,font=2)
}
if(A$type=="t"){
rect(A$loc[1],A$loc[2],A$loc[3],A$loc[4])
dx<-(A$loc[3]-A$loc[1])/ncol(A$Table)
dy<-(A$loc[4]-A$loc[2])/nrow(A$Table)
for(i in 1:ncol(A$Table)){
for(j in 1:nrow(A$Table)){
rect(A$loc[1]+(i-1)*dx,A$loc[4]-j*dy,A$loc[1]+i*dx,A$loc[4]-(j-1)*dy)
text(A$loc[1]+i*dx-dx/2,A$loc[4]-j*dy+dy/2,adj=c(0.5,0.5),labels=A$Table[j,i],cex=A$size[3])
}}
for(j in 1:nrow(A$Table)){text(A$loc[1]-1,A$loc[4]-j*dy+dy/2,adj=c(1,0.5),cex=A$size[1],labels=rownames(A$Table)[j])}
for(i in 1:ncol(A$Table)){text(A$loc[1]+i*dx-dx/2,A$loc[4]+1,adj=c(0.5,0),cex=A$size[2],labels=colnames(A$Table)[i])}
}
}


num.display<-function(xL=80,yB=4,dx=5,dy=5,fsp=1.5,conditions=c(-Inf,Inf),integer=F){
keyb<-matrix(NA,nrow=13,ncol=5)
rownames(keyb)<-keyb[,5]<-c(0,".","del",1:9,"enter")
colnames(keyb)<-c("xL","yB","xR","yT","key")
rect(xL,yB,xL+dx*3,yB+dy*4,col="yellow")
for(i in 1:13){
if(i==13){rect(xL,yB-dy,xL+dx*3,yB,col="green")
keyb[i,1:4]<-c(xL,yB-dy,xL+dx*3,yB)
text((2*xL+dx*3)/2,yB-dy+dy/2,adj=c(0.5,0.5),col="black",cex=1.5,labels="enter")}
else{rect(xL+(i+2)%%3*dx,yB+(i-1)%/%3*dy,xL+(i+2)%%3*dx+dx,yB+(i-1)%/%3*dy+dy)
keyb[i,1:4]<-c(xL+(i+2)%%3*dx,yB+(i-1)%/%3*dy,xL+(i+2)%%3*dx+dx,yB+(i-1)%/%3*dy+dy)
text(xL+dx/2+(i+2)%%3*dx,dy/2+yB+(i-1)%/%3*dy,adj=c(0.5,0.5),labels=rownames(keyb)[i],col="darkgreen",cex=1.5)}
}
repeat{
repeat{
answer<-NULL
rect(xL,yB+dy*4+2,xL+dx*3,yB+dy*4+2+dy,col="orange")
sp<-0
repeat{
repeat{
ck<-locator(n=1);if(ck$x>xL & ck$x<(xL+3*dx) & ck$y>(yB-dy) & ck$y<(yB+dy*4)) break}
u<-keyb[(ck$x>as.numeric(keyb[,1]) & ck$x<as.numeric(keyb[,3]) & ck$y>as.numeric(keyb[,2]) & ck$y<as.numeric(keyb[,4])),5]
if(u!="enter"){answer<-c(answer,u)}
if(u=="enter" | u=="del") break
if(u!="enter"){text(xL+dx/4+sp,yB+dy*4+2+dy/2,adj=c(0.5,0.5),col="darkred",cex=1.5,labels=u)}
if(u!="del"){text(xL+dx/4+sp,yB+dy*4+2+dy/2,adj=c(0.5,0.5),col="darkred",cex=1.5,labels=u)}
sp<-sp+fsp}
answer<-suppressWarnings(as.numeric(paste(answer,collapse="")))
if(u=="enter"){E<-T}else{E<-F}
if(is.na(answer)!=T && E==T) break}
int<-T
if(integer==T && (answer!=round(answer,0))){int<-F}else{int<-T}
if(answer>=conditions[1]  && answer<=conditions[2]  && int==T) break}
return(answer)}

result<-NULL
if(type==1){
close.screen(all.screens=TRUE)
split.screen(figs=rbind(c(0,1,0.95, 1), c(0,1,0, 0.95)))
screen(1)
par(mar=c(0,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(50,50,adj=c(0.5,0.5),cex=3,labels=title)
screen(2)
par(mar=c(1,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(h,v,adj=c(0,1),cex=size,col=col,labels=content)
if(is.null(Extra)!=T){for(k in 1:length(Extra)){Annex(A=Extra[[k]])}}
if(is.null(Line.draw)!=T){for(k in 1:length(Line.draw)){eval(parse(text=Line.draw[k]))}}
Sys.sleep(time)
}
if(type==2){
close.screen(all.screens=TRUE)
split.screen(figs=rbind(c(0,1,0.95, 1), c(0,1,0.1, 0.95),c(0,1,0,0.1)))
screen(1)
par(mar=c(0,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(50,50,adj=c(0.5,0.5),cex=3,labels=title)
screen(2)
par(mar=c(0,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(h,v,adj=c(0,1),cex=size,col=col,labels=content)
if(is.null(Extra)!=T){for(k in 1:length(Extra)){Annex(A=Extra[[k]])}}
if(is.null(Line.draw)!=T){for(k in 1:length(Line.draw)){eval(parse(text=Line.draw[k]))}}
screen(3)
par(mar=c(1,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
rect(20,5,80,95,col="green")
text(50,50,adj=c(0.5,0.5),cex=2,labels="Press here to proceed to the next page",col="darkgreen")
repeat{ck<-locator(n=1);
if(ck$x>20 & ck$x<80 & ck$y>5 &ck$y<95) break}
}
if(type==3){
close.screen(all.screens=TRUE)
split.screen(figs=rbind(c(0,1,0.95, 1), c(0,1,0, 0.95)))
screen(1)
par(mar=c(0,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(50,50,adj=c(0.5,0.5),cex=3,labels=title)
screen(2)
par(mar=c(1,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(h,v,adj=c(0,1),cex=size,col=col,labels=content)
if(is.null(Extra)!=T){for(k in 1:length(Extra)){Annex(A=Extra[[k]])}}
if(is.null(Line.draw)!=T){for(k in 1:length(Line.draw)){eval(parse(text=Line.draw[k]))}}
result<-num.display(xL=box.level[1],yB=box.level[2]-27,conditions=interval,integer=isinteger)
}
if(type==4){
close.screen(all.screens=TRUE)
split.screen(figs=rbind(c(0,1,0.95, 1), c(0,1,0, 0.95)))
screen(1)
par(mar=c(0,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(50,50,adj=c(0.5,0.5),cex=3,labels=title)
screen(2)
par(mar=c(1,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(h,v,adj=c(0,1),cex=size,col=col,labels=content)
if(is.null(Extra)!=T){for(k in 1:length(Extra)){Annex(A=Extra[[k]])}}
if(is.null(Line.draw)!=T){for(k in 1:length(Line.draw)){eval(parse(text=Line.draw[k]))}}
dx<-(100-(layout[2]-1))/layout[2];dy<-max(unlist(lapply(strsplit(choices,split="\n"),length)))*5
yT<-box.level
chcmatrix<-matrix(choices,nrow=layout[1],ncol=layout[2],byrow=T);locmatrix<-NULL
repeat{
for(j in 1:layout[1]){
for(i in 1:layout[2]){
locmatrix<-rbind(locmatrix,c(j,i,0+(i-1)*(dx+1),yT-(j-1)*(dy+2),0+(i-1)*(dx+1)+dx,yT-(j-1)*(dy+2)-dy))
rect(0+(i-1)*(dx+1),yT-(j-1)*(dy+2),0+(i-1)*(dx+1)+dx,yT-(j-1)*(dy+2)-dy,col="yellow")
text(0+(i-1)*(dx+1)+dx/2,yT-(j-1)*(dy+2)-dy/2,adj=c(0.5,0.5),cex=2,labels=chcmatrix[j,i])
}}
rect(75,yT-(layout[1]-1)*(dy+2)-9-dy,85,yT-(layout[1]-1)*(dy+2)-2-dy,col="green");text(80,yT-(layout[1]-1)*(dy+2)-5.5-dy,adj=c(0.5,0.5),labels="enter",cex=2)
rect(90,yT-(layout[1]-1)*(dy+2)-9-dy,100,yT-(layout[1]-1)*(dy+2)-2-dy,col="yellow");text(95,yT-(layout[1]-1)*(dy+2)-5.5-dy,adj=c(0.5,0.5),labels="delete",cex=2)
if(n.select>1){text(55,yT-(layout[1]-1)*(dy+2)-5.5-dy,adj=c(0.5,0.5),cex=1.5,labels=paste("You may choose exactly",n.select,"selections"))}
if(n.select==1){text(55,yT-(layout[1]-1)*(dy+2)-5.5-dy,adj=c(0.5,0.5),cex=1.5,labels=paste("Only one selection is permitted"))}
result.set<-NULL
repeat{
ck<-locator(n=1)
for(j in 1:layout[1]){
for(i in 1:layout[2]){
if(ck$x>0+(i-1)*(dx+1) & ck$x<0+(i-1)*(dx+1)+dx & ck$y<yT-(j-1)*(dy+2) & ck$y>yT-(j-1)*(dy+2)-dy)
{rect(0+(i-1)*(dx+1),yT-(j-1)*(dy+2),0+(i-1)*(dx+1)+dx,yT-(j-1)*(dy+2)-dy,col="orange")
text(0+(i-1)*(dx+1)+dx/2,yT-(j-1)*(dy+2)-dy/2,adj=c(0.5,0.5),cex=2,labels=chcmatrix[j,i])
result.set<-c(result.set,(j-1)*layout[2]+i);result.set<-union(result.set,result.set)}
}} #i & j
fB<-yT-(layout[1]-1)*(dy+2)-9-dy;fT<-yT-(layout[1]-1)*(dy+2)-2-dy
if( (ck$x>90 & ck$x<100 & ck$y>fB & ck$y<fT) || (length(result.set)==n.select & (ck$x>75 & ck$x<85 & ck$y>fB & ck$y<fT))) break}
if(ck$x>75 & ck$x<85 & ck$y>fB & ck$y<fT) break}
result<-result.set}
if(type==5){
W<-C<-numeric(0)
for(t in 1:length(cond)){
repeat{
close.screen(all.screens=TRUE)
split.screen(figs=rbind(c(0,1,0.95, 1), c(0,1,0, 0.95)))
screen(1)
par(mar=c(0,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(50,50,adj=c(0.5,0.5),cex=3,labels=title)
screen(2)
par(mar=c(1,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
if(is.null(Extra)!=T){for(k in 1:length(Extra)){Annex(A=Extra[[k]])}}
if(is.null(Line.draw)!=T){for(k in 1:length(Line.draw)){eval(parse(text=Line.draw[k]))}}
text(h[1],v[1],adj=c(0,1),cex=size[1],col=col[1],labels=paste(content[1],cond[t]))
if(length(content)>1){text(h[2:length(content)],v[2:length(content)],adj=c(0,1),cex=size[2:length(content)],col=col[2:length(content)],labels=content[2:length(content)])}
y.top<-box.level
text(0,y.top+10,adj=c(0,1),labels="Single click one of the options",col="green",cex=2)
choices.m<-scale
x.length<-(95-3*(length(scale)-1))/length(scale)
y.length<-10
for(j in 1:length(scale)){
x.left<-3+(x.length+3)*(j-1)
x.right<-3+x.length+(x.length+3)*(j-1)
y.top<-box.level
y.bottom<-y.top-y.length
rect(x.left,y.bottom,x.right,y.top,col="yellow")
text((x.left+x.right)/2,(y.bottom+y.top)/2,adj=c(0.5,0.5),labels=choices.m[j],cex=2)}
select<-locator(n=1);valid.select<-NULL
for(j in 1:length(scale)){
x.left<-3+(x.length+3)*(j-1)
x.right<-3+x.length+(x.length+3)*(j-1)
y.top<-box.level
y.bottom<-y.top-y.length
if(select$x>x.left & select$x<x.right & select$y>y.bottom & select$y<y.top)
{rect(x.left,y.bottom,x.right,y.top,col="orange")
box.selected<-c(x.left,x.right,y.top,y.bottom);selection<-scale[j]
valid.select<-c(valid.select,T)}
else{rect(x.left,y.bottom,x.right,y.top)
valid.select<-c(valid.select,NULL);S<-NULL}
text((x.left+x.right)/2,(y.bottom+y.top)/2,adj=c(0.5,0.5),labels=choices.m[j],cex=2)}
if(is.null(valid.select)!=T){
text(0,y.top+5,adj=c(0,1),labels="Click your selection again to confirm",col="green",cex=2)
confirm<-locator(n=1)
if(confirm$x>box.selected[1] && confirm$x<box.selected[2] && confirm$y<box.selected[3] && confirm$y>box.selected[4]){S<-selection}}
if(is.null(S)!=T) break}
W<-c(W,S)
C<-c(C,cond[t])}
result<-cbind(C,W)}

if(type==6){
repeat{
close.screen(all.screens=TRUE)
split.screen(figs=rbind(c(0,1,0.95, 1), c(0,1,0, 0.95)))
screen(1)
par(mar=c(0,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(50,50,adj=c(0.5,0.5),cex=3,labels=title)
screen(2)
par(mar=c(1,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
if(is.null(Extra)!=T){for(k in 1:length(Extra)){Annex(A=Extra[[k]])}}
if(is.null(Line.draw)!=T){for(k in 1:length(Line.draw)){eval(parse(text=Line.draw[k]))}}
text(h,v,adj=c(0,1),cex=size,col=col,labels=content)
text(0,box.level+7,adj=c(0,1),labels="To select, place the cursor on the below measurement tape and single click the mouse",cex=2,col="green")
r<-sort(range)
for(i in 1:(length(r)-1)){
L<-2+(i-1)*95/(length(r)-1)
R<-2+i*95/(length(r)-1)
rect(L,box.level-3,R,box.level,col="skyblue")
text(L,box.level-8,adj=c(0.5,0),labels=r[i],col="darkblue",cex=2.5)
for(i in 1:as.numeric(grid)){
lines(rep(i*(R-L)/as.numeric(grid)+L,2),c(box.level-3,box.level-2))}}
text(97,box.level-8,adj=c(0.5,0),labels=r[length(r)],col="darkblue",cex=2.5)
if(is.null(unit)!=T){text(50,box.level-13,adj=c(0.5,0.5),labels=paste("in",unit),col="darkblue",cex=2.5)}
rect(5,box.level-25,30,box.level-18,col="green",lwd=2);rect(35,box.level-25,45,box.level-18,col="yellow",lwd=2)
gb<-c(5,box.level-25,30,box.level-18);yb<-c(35,box.level-25,45,box.level-18)
text(17.5,(gb[2]+gb[4])/2,adj=c(0.5,0.5),labels="Submit selection(s)",cex=2)
text(40,(yb[2]+yb[4])/2,adj=c(0.5,0.5),labels="Delete",cex=2)
res<-NULL
repeat{W<-locator(n=1)
if(W$x>=2 & W$x<=97 & W$y>box.level-5 & W$y<box.level){res<-c(res,(W$x-2)*(max(r)-min(r))/95+min(r))
lines(rep(W$x,2),c(box.level+2,box.level-4),col="red",lwd=2)}
if(((W$x>gb[1] && W$x<gb[3] && W$y<gb[4] && W$y>gb[2])&(length(res)==n.select))||(W$x>yb[1] & W$x<yb[3] & W$y>yb[2] & W$y<yb[4])||(length(res)>n.select)) break}
if(W$x>gb[1] && W$x<gb[3] && W$y<gb[4] && W$y>gb[2]) break}
result<-res
}

if(type==7){
repeat{ ##switch
X<-initial
repeat{
repeat{
dmm<-NULL
close.screen(all.screens=TRUE)
split.screen(figs=rbind(c(0,1,0.95, 1), c(0,1,0, 0.95)))
screen(1)
par(mar=c(0,1,0,1)+0)
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
text(50,50,adj=c(0.5,0.5),cex=3,labels=title)
screen(2)
par(mar=c(1,1,0,1))
plot(-10,-10,xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",xaxt="n",yaxt="n")
if(is.null(Extra)!=T){for(k in 1:length(Extra)){Annex(A=Extra[[k]])}}
if(is.null(Line.draw)!=T){for(k in 1:length(Line.draw)){eval(parse(text=Line.draw[k]))}}
text(h[1],v[1],adj=c(0,1),cex=size[1],col=col[1],labels=paste(content[1],X))
if(length(content)>1){text(h[2:length(content)],v[2:length(content)],adj=c(0,1),cex=size[2:length(content)],col=col[2:length(content)],labels=content[2:length(content)])}
text(0,box.level+5,adj=c(0,1),labels="Please click one of the followings",cex=2,col="green")
rect(10,box.level-10,25,box.level,col="yellow",lwd=2);rect(30,box.level-10,45,box.level,col="yellow",lwd=2)
yb<-c(10,box.level-10,25,box.level);nb<-c(30,box.level-10,45,box.level)
text(17.5,box.level-5,adj=c(0.5,0.5),labels="Yes",cex=2.5)
text(37.5,box.level-5,adj=c(0.5,0.5),labels="No",cex=2.5)
repeat{YN<-locator(n=1,type = "n")
if(YN$x>yb[1] & YN$x<yb[3] & YN$y<yb[4] & YN$y>yb[2]) {dmm<-F;rect(yb[1],yb[2],yb[3],yb[4],lwd=2,col="green");text(17.5,box.level-5,adj=c(0.5,0.5),labels="Yes",cex=2.5)}
if(YN$x>nb[1] & YN$x<nb[3] & YN$y<nb[4] & YN$y>nb[2]) {dmm<-T;rect(nb[1],nb[2],nb[3],nb[4],lwd=2,col="green");text(37.5,box.level-5,adj=c(0.5,0.5),labels="No",cex=2.5)}
if(is.null(dmm)!=T) break}
text(0,box.level-12,adj=c(0,1),labels="Please click your selection again to confirm",cex=2,col="green")
CF<-locator(n=1,type="n")
if(CF$x>yb[1] & CF$x<yb[3] & CF$y<yb[4] & CF$y>yb[2]) {dmm<-c(dmm,F)}
if(CF$x>nb[1] & CF$x<nb[3] & CF$y<nb[4] & CF$y>nb[2]) {dmm<-c(dmm,T)}
if (length(dmm)==2 && dmm[1]==dmm[2]) break} #2nd repeat
if(union(dmm,dmm)==T | X<=0) break #1st repeat
if(order=="ascending"){X<-X+increment}
if(order=="descending"){X<-X-increment}
}
if(X<=0){
for(i in 0:10){
rect(0,box.level-20,90,box.level-17,border="white",col="white")
text(0,box.level-17,adj=c(0,1),labels=paste("Error:",content[1],X,"?  Redo starts in",10-i,"second(s)"),cex=1.5,col="red")
Sys.sleep(1)}}
if(X>0) break} ##switch
result<-X
}
if(is.matrix(result)!=T ){result<-matrix(c(rep(id,length(result)),rep(gsub(" ","",title),length(result)),rep(type,length(result)),rep(NA,length(result)),result),ncol=5)}else
{result<-cbind(rep(id,nrow(result)),rep(gsub(" ","",title),nrow(result)),rep(type,nrow(result)),result);result<-as.matrix(result)}
colnames(result)<-c("ID","Question.number","type","Condition.Likert","Response");rownames<-NULL
return(result)
}
