`plotMvloc` <-
function(est1,est2 = NULL ,est3 = NULL, X=NULL, alim = NULL, color.ell=2:4, color.points=grey(0.5), lty.ell= rep(1,3), pch.ell=rep(16,3), lwd.ell=rep(1,3), cex.ell=rep(1,3),
         pch.points=1,level=0.95,npoints = 100, x.legend, y.legend, cex.legend = 1, pty="s", gap=1, oma.bottom, labels, cex.labels=2, main, ...)
{
opar <-  par(no.readonly = TRUE)
p <- ncol(est1$vcov)

alim <- match.arg(alim,c("both","ellipses"))
    

legend.text<-NULL ; legend.pch<-NULL; legend.col<-NULL; legend.lty<-NULL
    
if (!is.null(est1)) {legend.text<-est1$est.name ; legend.pch<-pch.ell[1]; legend.col<-color.ell[1]; legend.lty<-lty.ell[1]}
if (!is.null(est2)) {legend.text<-c(legend.text,est2$est.name) ; legend.pch<-c(legend.pch,pch.ell[2]); legend.col<-c(legend.col,color.ell[2]); legend.lty<-c(legend.lty,lty.ell[2])}
if (!is.null(est3)) {legend.text<-c(legend.text,est3$est.name) ; legend.pch<-c(legend.pch,pch.ell[3]); legend.col<-c(legend.col,color.ell[3]); legend.lty<-c(legend.lty,lty.ell[3])}
    
l.legend <- length(legend.text)
if (missing(oma.bottom)) oma.bottom <- l.legend*2*cex.legend+2
oma.top <- 4
if (!missing(main)) {oma.top=oma.top+1} else main<-NULL
    
if (missing(labels)) {
        if (!is.null(X) & !is.null(colnames(X)))  {labels <- colnames(X)} 
        else {labels <- paste("var", 1:p)}
    }

par(mfrow=c(p,p),mar=c(gap/2,gap/2,gap/2,gap/2), oma=c(oma.bottom,4,oma.top,4), xpd = NA, pty=pty)

dots <- list(...)
nmdots <- names(dots)

for (j in 1:p)
    {
    for (i in 1:p) 
        {
        if (i==j) { par(usr=c(0,1,0,1))
                    plot(0, 0, xlim=c(0,1),  ylim=c(0,1),type="n", xlab="", ylab="", axes=F)
                    box()
                    text(0.5,0.5,labels[i],cex=cex.labels)} else{
                   ell<-ellipse(est1$vcov, centre=c(est1$location[i],est1$location[j]), t = sqrt(qchisq(level, p)), which = c(i, j), npoints = npoints)     
           
                   if (!is.null(est2)) {ell2<-ellipse(est2$vcov, centre=c(est2$location[i],est2$location[j]), t = sqrt(qchisq(level, p)), which = c(i, j), npoints = npoints) }
                   else{ell2<-NULL}
                   if (!is.null(est3)) {ell3<-ellipse(est3$vcov, centre=c(est3$location[i],est3$location[j]), t = sqrt(qchisq(level, p)), which = c(i, j), npoints = npoints) }
                   else{ell3<-NULL}
                   if (!is.null(X)) {x <- X[,i]; y<-X[,j]}
                   else{x<-NULL; y<-NULL} 
                    switch(alim,
                        "both" = {y.limits <- range(c(y,ell[,2],ell2[,2],ell3[,2]))
                                  x.limits <- range(c(x,ell[,1],ell2[,1],ell3[,1]))
                                  ind<-0},
                         "ellipses"= {y.limits <- range(c(ell[,2],ell2[,2],ell3[,2]))
                                      x.limits <- range(c(ell[,1],ell2[,1],ell3[,1]))
                                      ind.x<- ifelse(x < x.limits[1] | x > x.limits[2], 1,0)
                                      ind.y<- ifelse(y < y.limits[1] | y > y.limits[2], 1,0)
                                      ind<-ind.x+ind.y}
                              )
           
           
           
           plot(ell,type="l",axes=F, ylim=y.limits,xlim=x.limits, col=color.ell[1], lty=lty.ell[1], lwd=lwd.ell[1], xlab="",ylab="")
           box()
           points(est1$location[i],est1$location[j],pch=pch.ell[1], col=color.ell[1], cex=cex.ell[1])
           
           #text(myRes$location[i],myRes$location[j], paste("j=",j,"i=",i), cex=1.5)
           if (j==1) Axis(ell[,1],side=3)
           if (j==p) Axis(ell[,1],side=1)
           if (i==1) Axis(ell[,2],side=2)
           if (i==p) Axis(ell[,2],side=4)
           
           if ((!is.null(x))) {
           points(x[ind==0],y[ind==0],col=color.points, pch=pch.points)
           }
        
           if ((!is.null(est2))) {
           points(ell2,type="l",col=color.ell[2], lty=lty.ell[2], lwd=lwd.ell[2])
           points(est2$location[i],est2$location[j],col=color.ell[2], pch=pch.ell[2], cex=cex.ell[2])
           }
    
           if ((!is.null(est3))) {
           points(ell3,type="l",col=color.ell[3], lty=lty.ell[3], lwd=lwd.ell[3])
           points(est3$location[i],est3$location[j],col=color.ell[3], pch=pch.ell[3], cex=cex.ell[3])
           }
           }
           }
    }
    
    if (!is.null(main)) {
                           font.main <- if ("font.main" %in% nmdots) 
                           dots$font.main
                           else par("font.main")
                           cex.main <- if ("cex.main" %in% nmdots) 
                           dots$cex.main
                           else par("cex.main")
                           mtext(main,3,3,TRUE,0.5, cex = cex.main, font = font.main)}
    invisible(NULL)
    par(xpd=NA)


if (missing(x.legend))
    {
    if (p==2) x.legend <- -0.1
    if (p==3) x.legend <- -0.9
    if (p>=4) x.legend <- -1
    }
if (missing(y.legend))
    {
    if (l.legend==1 & p==2) y.legend <- -0.3
    if (l.legend==2 & p==2) y.legend <- -0.45
    if (l.legend==3 & p==2) y.legend <- -0.6
    if (l.legend==1 & p==3) y.legend <- -0.4
    if (l.legend==2 & p==3) y.legend <- -0.55
    if (l.legend==3 & p==3) y.legend <- -0.65
    if (l.legend==1 & p>=4) y.legend <- -0.45
    if (l.legend==2 & p>=4) y.legend <- -0.60
    if (l.legend==3 & p>=4) y.legend <- -0.85
    }

if (!is.null(x.legend) & !is.null(y.legend))    legend(x.legend,y.legend, legend=legend.text, col=legend.col, pch=legend.pch, lty=legend.lty ,yjust=0, xjust=0.5,bty="n", cex=cex.legend) 
   
par(opar)
}
