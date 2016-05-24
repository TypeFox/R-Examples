#kplot version 3.1
#
#       Version 2: adds optional ordered plotting and improves margin setting
#       Also adds specification of k (for pdf's) as (recycled) vector and 
#       adds marginal PDF plot
#
#       Version 2.1 adds graphical parameter expansion to full length to follow ordering 
#       correctly and also adds control of error bar cols etc. (NB: NA for col.ci suppresses
#       error bar as per arrows() )
#
#       Version 3 adds a generic and default method together with an implementation for the
#       ilab class
#
#       V3.1 is a bug-fix for the ilab method, which was not properly implementing 
#       do.pdf=TRUE
#
# Wish List: 
#       - Implement strata
#       - Add layout as argument instead of pdf.layout

kplot<-function(x, ...) {
        UseMethod("kplot")
}

kplot.ilab<-function(x, ...) {
        pars <- list(...)
        if(is.null(pars$x)) pars$x <- x$data$x
        if(is.null(pars$labels)) pars$labels <- as.character(x$data$org)
        if(is.null(pars$U.lo)) pars$U.lo <- x$data$U.lower
        if(is.null(pars$U.hi)) pars$U.hi <- x$data$U.upper
        if(is.null(pars$main)) pars$main <- x$title[1]
        if(!is.null(pars$do.pdf)) {
                if(pars$do.pdf) {
                        if(any(pars$U.lo != pars$U.hi )) {
                                warning("do.pdf=TRUE does not support asymetry; using mean U")
                                pars$U <- ( x$data$U.lower + x$data$U.upper )/2
                        }
                        if(is.null(pars$U)) pars$U <- x$data$u * x$data$k 
                }
        
        }
        if(!is.na(x$subset)) pars$main <- paste(pars$main, "\nSubset: ", x$subset, sep="")
        
        rv<-do.call(kplot.default, pars)
        invisible(rv)
}

kplot.default<-function(x,U=NULL, labels=names(x),  assigned=NULL, U.assigned=NULL, 
        U.lo=U, U.hi=U, k=2, strata=NULL,
        do.percent=!is.null(assigned) && !do.pdf, 
        ordered=TRUE, order.strata=levels(strata),
        xlim=c(0.5, length(x)+0.5), ylim,
        main=NULL, xlab=NULL, ylab=NULL,
        axis.main=2, axis.pct=4, at=1:length(x), at.main=NULL,
        cex.axis=0.8, las=2, las.pct=1, ylab.line=2.5, ylab.line.pct=2.1,
        ci.width=0.03, col.ci=par("fg"), lty.ci=par("lty"), lwd.ci=par("lwd"),
        pch=21, col=par("fg"), bg="white", add.outliers=FALSE, outlier.offset=0.2,  
        mar=NULL, box=TRUE,
        do.pdf=FALSE, 
                do.individual.pdf=do.pdf, col.pdf=par("fg"), lwd.pdf=1, lty.pdf=1,
                do.total.pdf=TRUE, col.total.pdf=col.pdf[1], lwd.total.pdf=2, lty.total.pdf=1,
                n.pdf=200, pdf.layout=c(4,1), pdf.scale=0.7, pdf.offset=0.05, xlim.pdf, 
                pdf.axis=FALSE, las.pdf=0, mgp.pdf=c(3,0.5,0),
        ...)
{
        oo<-if(ordered) order(x) else 1:length(x)
        
        Lx<-length(x)
        
        #expand graphics pars to full length to simplify ordering
        col <- rep(col, length.out=Lx)
        bg <- rep(bg, length.out=Lx)
        pch <- rep(pch, length.out=Lx)
        col.ci <- rep(col.ci, length.out=Lx)
        lty.ci <- rep(lty.ci, length.out=Lx)
        lwd.ci <- rep(lwd.ci, length.out=Lx)
        
        upper=x+U.hi
        lower=x-U.lo
        if(length(k)<Lx) k <- rep(k, length.out=Lx)
        if(missing(xlim)) xlim <- c(0.5, Lx+0.5)
        if(missing(ylim)) ylim <- range(pretty(na.omit(c(x, upper,lower))))
        if(is.null(at.main)) at.main<-pretty(ylim)
        
                
        if(do.pdf) layout(matrix(c(1,2),ncol=2), widths=pdf.layout)
                else layout(matrix(1))
        
        axis.RHS <- (axis.main==4 | (do.percent & axis.pct==4) )
        mar.adj.do.pdf<-if(do.pdf) c(-0.5,-1.5) else c(0,0)
        if(is.null(mar))
                mar<-if(axis.RHS) c(5,4,4,4+mar.adj.do.pdf[1])+0.1 else c(5,4,4,2+mar.adj.do.pdf[2])+0.1 #if(do.percent)

        
        plot.new()
        par(mar=mar)
        
        plot.window(xlim=xlim, ylim=ylim)

        if(box) box()
        axis(1,at=at, labels=labels[oo], las=las, cex.axis=cex.axis, srt=45)
        if(axis.main) {
                axis(axis.main, at=at.main, las=las, cex.axis=cex.axis)
                mtext(ylab, side=axis.main, line=ylab.line)
        }

        if(do.percent & !is.null(assigned)) {
                at.percent<-pretty( 100*(ylim-assigned)/assigned )
                if(axis.pct) {
                        axis(axis.pct,at=assigned*(1+at.percent/100), labels=paste( at.percent ), las=las.pct, cex.axis=cex.axis)
                        mtext("% Deviation", side=axis.pct, line=ylab.line.pct)
                }
        }
        
        if(!is.null(assigned)) {
                abline(h=assigned, lwd=2)
                if(!is.null(U.assigned)) abline(h=assigned+c(-U.assigned, U.assigned), lwd=1, lty=2)

        }
        
        arrows(at,lower[oo], at, upper[oo], length=ci.width, code=3,angle=90, col=col.ci[oo], lwd=lwd.ci[oo], lty=lty.ci[oo])
        points(at,x[oo], pch=pch[oo], bg=bg[oo], col=col[oo])
        
        if(add.outliers) {
                a.len<-diff(ylim)/6
                high.index<-which( x[oo]>max(ylim) )
                low.index<-which( x[oo]<min(ylim) )
                
                if(length(high.index)>0) {
                        arrows(at[high.index], max(ylim)-a.len, at[high.index], max(ylim), length=ci.width*1.5, angle=20)
                        text(x=at[high.index]-outlier.offset, y=max(ylim)-a.len, labels=paste(x[oo][high.index]), adj=c(0.5,0.5), srt=90, cex=0.7)
                }
                
                if(length(low.index)>0) {
                        arrows(at[low.index], min(ylim)+a.len, at[low.index], min(ylim), length=ci.width*1.5, angle=20)
                        text(x=at[low.index]-outlier.offset, y=min(ylim)+a.len, labels=paste(x[oo][low.index]), adj=c(0.5,0.5), srt=90, cex=0.7)
                }
        }
        
        if(is.null(main))
                main<-paste("kplot -", deparse(substitute(x)))
        title(main=main, xlab=xlab)
        
        if(do.pdf) {
                if(pdf.axis & !axis.RHS) axis(4, at=at.main, las=las.pdf, mgp=mgp.pdf, cex.axis=cex.axis) #adds extra RHS axis
                #Calculate (normal) pdf's
                d.pdf<-matrix(ncol=n.pdf, nrow=Lx)
                x.pdf<-seq(par("usr")[3],par("usr")[4],length.out=n.pdf) #Formerly based on ylim[1:2]

                for(i in 1:Lx) d.pdf[i,]<-dnorm(x.pdf,x[i],U[i]/k[i])/Lx
                d.total<-colSums(d.pdf)

                par(mar=c(mar[1],(pdf.axis & !do.percent)*2,mar[2],1))
                
                plot.new() #Advances 'frame'
                if(missing(xlim.pdf)) {
                        if(do.total.pdf) 
                                xlim.pdf<-c(0,max(pretty(d.total))) 
                        else
                                xlim.pdf<-c(0,max(pretty(d.pdf))) 
                }
                plot.window(xlim=xlim.pdf, ylim=ylim)
                
                if(do.individual.pdf) {
                        col.pdf=rep(col.pdf, length.out=Lx)
                        lwd.pdf=rep(lwd.pdf, length.out=Lx)
                        lty.pdf=rep(lty.pdf, length.out=Lx)
                        for(i in 1:Lx) 
                                lines(d.pdf[i,], x.pdf, col=col.pdf[i], lwd=lwd.pdf[i],lty=lty.pdf[i])
                }
                if(do.total.pdf) 
                        lines(d.total, x.pdf, col=col.total.pdf, lwd=lwd.total.pdf, lty=lty.total.pdf)
                

                if(pdf.axis & !axis.RHS) axis(2, at=at.main, cex.axis=cex.axis, labels=F, line=0.2, ) #adds tick-only axis
        }
        return(invisible(list(order=oo, at=at)))
}



kpoints<-function(x,U=NULL, labels=names(x),  
        U.lo=U, U.hi=U, k=2, strata=NULL,
        ordered=TRUE, order.strata=levels(strata),
        at=1:length(x), 
        ci.width=0.03, col.ci=par("fg"), lty.ci=par("lty"), lwd.ci=par("lwd"),
        pch=21, col=par("fg"), bg="white", add.outliers=FALSE, outlier.offset=0.2,  
        ...)
{
        oo<-if(ordered) order(x) else 1:length(x)
        
        Lx<-length(x)
        
        #expand graphics pars to full length to simplify ordering
        col <- rep(col, length.out=Lx)
        bg <- rep(bg, length.out=Lx)
        pch <- rep(pch, length.out=Lx)
        col.ci <- rep(col.ci, length.out=Lx)
        lty.ci <- rep(lty.ci, length.out=Lx)
        lwd.ci <- rep(lwd.ci, length.out=Lx)
        
        upper=x+U.hi
        lower=x-U.lo
        if(length(k)<Lx) k <- rep(k, length.out=Lx)
        xlim <- par("usr")[1:2]
        ylim <- par("usr")[3:4] 
                
        arrows(at,lower[oo], at, upper[oo], length=ci.width, code=3,angle=90, col=col.ci[oo], lwd=lwd.ci[oo], lty=lty.ci[oo])
        points(at,x[oo], pch=pch[oo], bg=bg[oo], col=col[oo])
        
        if(add.outliers) {
                a.len<-diff(ylim)/6
                high.index<-which( x[oo]>max(ylim) )
                low.index<-which( x[oo]<min(ylim) )
                
                if(length(high.index)>0) {
                        arrows(at[high.index], max(ylim)-a.len, at[high.index], max(ylim), length=ci.width*1.5, angle=20)
                        text(x=at[high.index]-outlier.offset, y=max(ylim)-a.len, labels=paste(x[oo][high.index]), adj=c(0.5,0.5), srt=90, cex=0.7)
                }
                
                if(length(low.index)>0) {
                        arrows(at[low.index], min(ylim)+a.len, at[low.index], min(ylim), length=ci.width*1.5, angle=20)
                        text(x=at[low.index]-outlier.offset, y=min(ylim)+a.len, labels=paste(x[oo][low.index]), adj=c(0.5,0.5), srt=90, cex=0.7)
                }
        }
        
        return(invisible(list(order=oo, at=at)))
}

