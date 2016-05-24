pave<-function(cordseg,yc,xc,fix.edge=NULL,ydown=TRUE,output="list"){
if (is.data.frame(cordseg)) cordseg<-as.matrix(cordseg)
cordseg<-matrix(cordseg,ncol=2)
d<-dim(cordseg)
if(d[2]!=2 | d[1]!=2) stop("Segment coordinates must be a 2 x 2 matrix")
cordseg<-cordseg[order(cordseg[,1]),]
if (cordseg[1,1]==cordseg[2,1]) cordseg <- cordseg[rev(order(cordseg[, 2])), ]
choice<-c("list", "points", "spdf")
m<-pmatch(output, choice)
if (!is.na(m)) output<-choice[m] else stop("output must be \"list\", \"points\" or \"spdf\"")
flag=TRUE; if (cordseg[1,2]>cordseg[2,2]) flag=FALSE

if(!is.null(fix.edge)) {
    d1<-sqrt((cordseg[1,1]-cordseg[2,1])^2+(cordseg[1,2]-cordseg[2,2])^2)
    d2<-yc*fix.edge
    x3<-cordseg[1,1]+d2*abs(cordseg[1,1]-cordseg[2,1])/d1
    if (flag) {
    y3<-cordseg[1,2]+d2*abs(cordseg[1,2]-cordseg[2,2])/d1
    }
    else y3<-cordseg[1,2]-d2*abs(cordseg[1,2]-cordseg[2,2])/d1
    cordseg[2,]<-cbind(x3,y3)
}

    xd<-seq(cordseg[1,1],cordseg[2,1],l=yc+1)
    yd<-seq(cordseg[1,2],cordseg[2,2],l=yc+1)
    L<-sqrt((xd[1]-xd[2])^2+(yd[1]-yd[2])^2)
    alp<-atan(abs(xd[1]-xd[2])/abs(yd[1]-yd[2]))
        res<-cbind(xd,yd)
        xd2<-xd;yd2<-yd
            for (a in 1:xc-1) {
                if (ydown){
                if(flag) xd2<-xd2+L*cos(alp) else xd2<-xd2-L*cos(alp)
                yd2<-yd2-L*sin(alp)
                res<-rbind(res,cbind(xd2,yd2))
                }
                else
                {
                if (flag) {xd2<-xd2-L*cos(alp);yd2<-yd2+L*sin(alp)}
                else {xd2<-xd2+L*cos(alp);yd2<-yd2+L*sin(alp)}
                res<-rbind(res,cbind(xd2,yd2))
                }
                
            } 
    
    xmat<-matrix(res[,1],nrow=yc+1,ncol=xc+1)
    ymat<-matrix(res[,2],nrow=yc+1,ncol=xc+1)
    
    poList<-rep(list(NA),(yc-1)*(xc-1))
    k<-0
    for (j in 1:xc){
        for (i in 1:yc){
            k<-k+1
            poList[[k]]<-matrix(c(xmat[i,j],ymat[i,j],xmat[i,j+1],ymat[i,j+1],xmat[i+1,j+1],ymat[i+1,j+1],xmat[i+1,j],ymat[i+1,j],xmat[i,j],ymat[i,j]),ncol=2,byrow=TRUE)
        }
    }

    if (output=="spdf") {
        p2<-rep(list(NA),length(poList))
        for (i in 1:length(poList)){
            p1<-Polygon(poList[[i]])
            p2[[i]]<-Polygons(list(p1),ID=i)
        }
        SPs <-SpatialPolygons(p2)
        SPDF <-SpatialPolygonsDataFrame(SPs, data.frame(ID=1:length(poList)))
    }
    else if (output=="points") res else poList 
}
