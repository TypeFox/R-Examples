boxplotlegend <- 
function(X,Y,el,boxinfo,x.shift=4e4,xf=1e4,y.shift=0.2,y.scale=13e4,legend.title="Legend",
         cex.legtit=1, logscale=TRUE, symb=c(1,1,16,3,3),ssize=c(1.5,1.0,0.3,1.0,1.5),
         accentuate=FALSE,cex.scale=0.8)
{
################################################################################
# Boxplot legend (in log10 scale):
#
# X ... X-coordinates
# Y ... Y-coordinates
# el ... variable considered
# boxinfo ... from boxplot(el) or boxplotlog(el)
# x.shift ... shift in x-direction
# xf ... width in x-direction
# y.shift ... shift in y-direction (from title)
# y.scale ... scale in y-direction
# legend.title ... title for legend
# logscale ... if TRUE plot boxplot in log-scale
# cex.legtit ... cex of title for legend
# symb ... symbols to be used (length 5!)
# ssize ... symbol sizes to be used (length 5!)
# cex.scale ... cex for text "log-scale" for scale
################################################################################

q=drop(boxinfo$stats)
out=boxinfo$out
# lower or upper outliers? If no, delete symbols!
low.out=TRUE
up.out=TRUE
if (!any(out>q[5])) {
     up.out=FALSE
     symb=symb[-length(symb)]
     ssize=ssize[-length(ssize)]
}
if (!any(out<q[1])) {
     low.out=FALSE
     symb=symb[-1]
     ssize=ssize[-1]
}


min.el=min(el)
max.el=max(el)

if (logscale){
  out=log10(out)
  qp=log10(c(min.el,q,max.el))
  mtext("log-scale",side=4,at=max(Y)-y.scale*0.4,adj=1,cex=cex.scale)
}
else {
  qp=c(min.el,q,max.el)
}

yc=(qp-qp[1])/(qp[7]-qp[1])-y.shift
out=(out-qp[1])/(qp[7]-qp[1])-y.shift
yct=max(Y)-y.scale*(1-yc)
out=max(Y)-y.scale*(1-out)
xc=max(X)-x.shift
segments(xc-xf,yct[1],xc+xf,yct[1])
segments(xc-xf/2,yct[2],xc+xf/2,yct[2])
segments(xc,yct[2],xc,yct[3])
segments(xc-xf,yct[3],xc+xf,yct[3])
segments(xc-xf,yct[4],xc+xf,yct[4])
segments(xc-xf,yct[5],xc+xf,yct[5])
segments(xc-xf,yct[3],xc-xf,yct[5])
segments(xc+xf,yct[3],xc+xf,yct[5])
segments(xc,yct[5],xc,yct[6])
segments(xc-xf/2,yct[6],xc+xf/2,yct[6])
segments(xc-xf,yct[7],xc+xf,yct[7])
points(rep(xc,length(out)),out,pch=16,cex=0.3)

# text labels
text(xc+1*xf,yct[2],roundpretty(q[1],2),pos=4,cex=cex.scale)
text(xc+1*xf,yct[3],roundpretty(q[2],1),pos=4,cex=cex.scale)
text(xc+1*xf,yct[5],roundpretty(q[4],1),pos=4,cex=cex.scale)
text(xc+1*xf,yct[6],roundpretty(q[5],1),pos=4,cex=cex.scale)

if (up.out) text(xc+1*xf,yct[7],roundpretty(max.el,1),pos=4,cex=cex.scale)
if (low.out) text(xc+1*xf,yct[1],roundpretty(min.el,1),pos=4,cex=cex.scale)

# plot symbols
symbnew <- symb
ssizenew <- ssize
qnew <- q[-3] # without median
if (!accentuate){
  ind=0
  if (low.out) {
    points(xc-3*xf,yct[1]+(yct[2]-yct[1])/2,pch=symb[1],cex=ssize[1])
    qnew <- c(min.el, qnew)
    ind <- 1
  }
  points(xc-3*xf,yct[2]+(yct[3]-yct[2])/2,pch=symb[1+ind],cex=ssize[1+ind])
  points(xc-3*xf,yct[3]+(yct[5]-yct[3])/2,pch=symb[2+ind],cex=ssize[2+ind])
  points(xc-3*xf,yct[5]+(yct[6]-yct[5])/2,pch=symb[3+ind],cex=ssize[3+ind])
  if (up.out) {
    points(xc-3*xf,yct[6]+(yct[7]-yct[6])/2,pch=symb[4+ind],cex=ssize[4+ind])
    qnew <- c(qnew, max.el)
  }
}
else {
  ind=0
  if (low.out) {
    points(xc-3*xf,yct[1]+(yct[2]-yct[1])/2,pch=symb[1],cex=ssize[1])
    qnew <- c(min.el, qnew)
    ind <- 1
  }
  points(xc-3*xf,yct[2]+(yct[3]-yct[2])/2,pch=symb[1+ind],cex=ssize[1+ind])
  points(xc-3*xf,yct[3]+(yct[5]-yct[3])/2,pch=symb[2+ind],cex=ssize[2+ind])
  points(xc-3*xf,yct[5]+(yct[6]-yct[5])/2,pch=symb[3+ind],cex=ssize[3+ind])
  if (up.out) {
    points(xc-3*xf,yct[6],pch=symb[4+ind],cex=ssize[4+ind]/1.3)
    points(xc-3*xf,yct[6]+(yct[7]-yct[6])/2,pch=symb[4+ind],cex=ssize[4+ind])
    points(xc-3*xf,yct[7],pch=symb[4+ind],cex=ssize[4+ind]*1.3)
    symbnew <- c(symbnew,symb[4+ind],symb[4+ind])
    ssizenew <- c(ssizenew[-(4+ind)],ssize[4+ind]/1.3,ssize[4+ind],ssize[4+ind]*1.3)
    qstep <- (max.el-q[4])/3
    qnew <- c(qnew, q[4]+qstep, q[4]+2*qstep, max.el)
  }
}
  
  
text(xc,max(Y),legend.title,cex=cex.legtit)
################################################################################
list(symb=symbnew,ssize=ssizenew,q=qnew)
}

