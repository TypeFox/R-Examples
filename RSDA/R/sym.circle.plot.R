sym.circle.plot <-
function(prin.corre) {
  v<-c("green", "red", "blue", "cyan", "brown", "yellow", 
         "pink","purple","orange","gray" ); 
  msg = paste("Correlation Circle");
  plot(-1.5:1.5, -1.5:1.5, type = "n",xlab = "C1", ylab = "C2",main=msg) 
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  symbols(0,0, circles=1, inches=FALSE, add=TRUE)
  c1=1
  c2=2 
  n<-dim(prin.corre)[1]
  f<-dim(prin.corre)[2]
  CRTI<-matrix(nrow=n, ncol=f)
  CRTI<-prin.corre
  vars <- rownames(prin.corre)
  for(k in 1 : n ) { 
    x1 <- min(CRTI[k,c1],CRTI[k,c2])
    x2 <- max(CRTI[k,c1],CRTI[k,c2])
    y1 <- min(CRTI[k,c2+1],CRTI[k,c2+2])
    y2 <- max(CRTI[k,c2+1],CRTI[k,c2+2])
    if(((x1>0)&&(x2>0)&&(y1>0)&&(y2>0))||((x1<0)&&(x2<0)&&(y1<0)&&(y2<0))) {
      plotX.slice(x1,y2,x2,y1,v,vars,k)
    }
    if(((x1<0)&&(x2<0)&&(y1>0)&&(y2>0))||((x1>0)&&(x2>0)&&(y1<0)&&(y2<0))) {
      plotX.slice(x1,y1,x2,y2,v,vars,k)
    }
    if((y1>0)&&(y2>0)&&(x1<0)&&(x2>0)) {
      plotX.slice(x1,y1,x2,y1,v,vars,k)
    }
    if((y1<0)&&(y2<0)&&(x1<0)&&(x2>0)) {
      plotX.slice(x1,y2,x2,y2,v,vars,k)
    }
    if((x1>0)&&(x2>0)&&(y1<0)&&(y2>0)) {
      plotX.slice(x1,y1,x1,y2,v,vars,k)
    }
    if((x1<0)&&(x2<0)&&(y1<0)&&(y2>0)) {
      plotX.slice(x2,y1,x2,y2,v,vars,k)
    }
    if((x1<0)&&(x2>0)&&(y1<0)&&(y2>0)) {
      plotX.slice(x2,y1,x2,y2,v,vars,k)
    }
  }
}
