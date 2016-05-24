extremexi <-
function(x,y,i1,i2,nt,alpha,xlb,ylb,xnd,ynd){
  #Find extreme of discrete data (xi,yi) with Taylor regression and degree of polynomial nt
  #Choose the desired range i1,i2 of data in order tosearch at interval [x_i1,x_i2] 
  #Choose the level of statistical significance as a number, ie 5 means 5%
  #Plot results with xlb=label for x-axis, ylb=label for y-axis
  #Set number of digits for x-axis -->xnd and for y-axis-->ynd
  #final version
  #Initialize:
  am<-0;lev=(100-alpha)/100;
  xm<-cbind();df<-NULL;
  x1<-x[i1:i2];y1<-y[i1:i2];
  vn<-rep(0,i2-i1);vdn<-rep(0,i2-i1);v0n<-rep(0,i2-i1);vpr<-rep(0,i2-i1);ii<-0;l2<-rep(0,i2-i1);
  aml<-list();amel<-list();yl<-list();
  #Search for all available extreme-points p=x_i
  for (i in (i1:i2)){
    ii<-ii+1;
    pr1<-x1[ii];vpr[ii]<-pr1;
    xm<-cbind();df<-NULL;
    for (i in (1:nt)){xm<-cbind(xm,cbind((x1-pr1)^i))}
    df<-as.data.frame(xm,row.names = NULL, optional = FALSE)
    xnam <- paste0("V", 1:(nt+0))
    fmla <- as.formula(paste("y1 ~ ", paste(xnam, collapse= "+")))
    c1<-lm(fmla,data=df);
    yf1<-predict(c1);yl<-list(yl,yf1);
    ci1<-confint(c1,level = lev);
    am<-matrix(ci1,ncol=2,dimnames=list(c(paste0("a",0:nt)),c(colnames(ci1))));
    aml<-list(aml,am);
    v0n[ii]=abs(c1$coeff[2]);
  }
  mans<-matrix(cbind(v0n,vpr),ncol=2,dimnames=list(c(),c("|a_1|","point")));
  rownames(mans) <- rownames(mans, do.NULL = FALSE, prefix = "#")
  #Find minimum value |a_1| and corresponding extreme x_i=p
  n0<-which(v0n == min(v0n));
  ipc<-mans[which.min(v0n),2];#print(n2);print(x1[n2]);
  ys<-matrix(unlist(yl),ncol=length(x1));
  yf1<-ys[,n0];
  ################
  #Plot results
  par(mfrow=c(1,2))
  ymin<-min(y1);
  ymax<-max(y1);
  dy<-(ymax-ymin)/5;
  plot(x1,y1,ylim=range(y1),col='blue',pch=19,cex=.5,axes=FALSE,ylab=ylb ,xlab=xlb);
  lines(x1,yf1,col='black',lty=1,lwd=2)
  dx<-(x1[length(x1)]-x1[1])/5;
  xticks<-round(c(x1[1],x1[n0],x1[length(x1)]),digits=xnd)
  axis(1,at=xticks,cex.axis=0.7,las=2)
  yticks<-round((seq(ymin,ymax,by=dy)),digits=ynd);
  axis(2,at=yticks,cex.axis=0.7)
  abline(h=yticks, v=xticks, col="gray", lty=3)
  abline(v=x1[n0],lty=2,col="red",lwd=3)
  legend('top',col=c('blue'),pch=c(19),legend=c('data'),bty='n',cex=0.6);
  legend('left',col=c('black','red'),lty=c(1,2),lwd=c(2,3),legend=c(paste0('Taylor fit (',nt,' )'),'extreme'),bty='n',cex=0.7)
  stit<-paste0('Data for [',toString(round(x1[1],digits=2)),',',toString(round(x1[length(x1)],digits=xnd)),'] \n');
  title(paste0(stit,'Taylor Regression n = ',nt,' , a = ',alpha,' % ','\n Find extreme'),cex.main=0.7)
  box()
  plot(vpr,v0n,col='blue',pch=19,cex=.5,axes=FALSE,ylab=expression(paste('|',alpha[1],'|')) ,xlab=expression(rho));
  abline(v=ipc,lty=3,col="blue",lwd=2)
  dp1<-(ipc-vpr[1])/2;
  dp2<-(vpr[length(vpr)]-ipc)/3;
  xticks<-round(c(vpr[1],ipc,vpr[length(vpr)]),digits=xnd)
  axis(1,at=xticks,cex.axis=0.7,las=2)
  dvn0<-(max(v0n)-min(v0n))/5;
  yticks<-c(seq(min(v0n),max(v0n),by=dvn0))
  axis(2,at=yticks,cex.axis=0.7,labels=round(yticks,digits=ynd))
  legend('top',col=c('blue'),pch=c(19),legend=c(expression(paste('|',alpha[1],'|'))),bty='n',cex=0.6);
  title(main=expression(paste('Plot of all available |',alpha[1],'|')),cex.main=0.7)  
  box() 
  n0k<-2*seq(1,nt/2);
  amf<-matrix(matrix(unlist(aml),nrow=nt+1)[,(2*n0-1):(2*n0)],ncol=2);
  par(mfrow=c(1,1))
  ans<-new.env();
  ans$an<-matrix(cbind(amf,0.5*rowSums(amf)),nrow=nt+1,ncol=3,byrow=F,list(c(paste0("a",c(0:nt))),c(colnames(ci1),"an")));
  ans$fextr<-c(n0,vpr[n0]);
  ans;
}
