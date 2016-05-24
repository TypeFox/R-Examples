duewer.plot<-function(x, ...) UseMethod("duewer.plot")

duewer.plot.default<-function(x,s,mu=median(x),sigma=mad(x), s0=median(s), labels=NA,
        radius=1:3, units=c("z","x"), 
        main, xlab, ylab, 
        xlim, ylim,  at.xax=NULL, at.yax=NULL, aspect, 
        col.contours="lightgrey", lty.contours=par("lty"), lwd.contours=par("lwd"),
        label.contours=T, format.clab="p=%4.3f",
        cex=par("cex"), cex.label=0.7, pos=3, adj=NULL, 
        pos.clab="bottomright", col.clab=col.contours,
        cex.axis=par("cex.axis"), pch=par("pch"), las=par("las"), 
        col=par("col"), bg=par("bg"), ...) {
        #Produces a Duewer concordance/apparent precision plot 
        #of f=s/s0 against z=(x-mu)/sigma, centred on zero, with
        #(semicircular) contours lines at (f^2+s^2)=radius
        
        #units: if units=="z", plot z-score and s/s0; if "x" plot x and s with boundaries scaled to 
        #       radius*sigma
        
        #pos.clab:      Contour labels for prob basis are placed approximately at the location indicated
        #               and adjusted outward appropriately. Options are ‘"top"’, ‘"topright"’, ‘"right"’,
        #               ‘"bottomright"’, ‘"bottom"’, ‘"bottomleft"’, ‘"left"’, ‘"topleft"’.
        #               For basis="prob" the positions are taken as provided. For basis="radius", 
        #               "bottomright" and "bottomleft" are as for "right" and "left" but just below the x-axis,
        #               and "bottom" is replaced with c("bottomright", "bottomleft").
        
        
        probs=1-(2*(1-pnorm(radius)))
        
        pos.clab <- match.arg(pos.clab, choices=c("top", "topright", "right",
                        "bottomright", "bottom", "bottomleft", "left", "topleft"), several.ok=TRUE)
                        
        units <- match.arg(units)
        
        if(units=="x") {
                z<-(x-mu)
                f<-s
                x.scale<-sigma
                if(missing(xlab)) xlab <- expression(x-mu)
                if(missing(ylab)) ylab <- deparse(substitute(s))
        } else {
                z<-(x-mu)/sigma
                f<-s/s0
                x.scale<-1
                if(missing(xlab)) xlab <- expression((x-mu)/sigma)
                if(missing(ylab)) ylab <- expression(s/s[0])
       }       

        if(missing(xlim)) {
                xlim=c(-1,1)*max(c(radius*x.scale,pretty(abs(z))), na.rm=T)
        }

        if(missing(ylim)) {
                ylim=c(0, max(c(radius*x.scale,pretty(abs(f))), na.rm=T) )
        
        }
        
        
        if(missing(aspect)){
                        aspect<-1.0
        }
        
        plot.new()
        
        plot.window(xlim=xlim, ylim=ylim, asp=aspect)
        
        
        #Calculate boundaries
        t<-seq(0,pi, length.out=200)
        max.sx<-0

        xpd <- par(xpd=TRUE)
        for(r in radius) {
                lines(x.x <- r*x.scale*cos(t), s.x <- r*x.scale*sin(t), col=col.contours, lty=lty.contours, lwd=lwd.contours)
                        max.sx <- max(max.sx, max(s.x))
                if(label.contours) {
                        if("bottom" %in% pos.clab) {
                                pos.clab <- unique(c(pos.clab[-match("bottom", pos.clab)], "bottomleft", "bottomright"))
                        }
                        for(i in 1:length(pos.clab)) {
                                l.pos <- switch(pos.clab[i], top=100, topright=50, right=1,
                                                bottomright=1, bottom=1, bottomleft=200, 
                                                left=200, topleft=150)
                                l.adj <- switch(pos.clab[i], top=c(-0.1,-0.3), topright=c(-0.1,-0.3), right=c(-0.1,-1),
                                                bottomright=c(-0.1,1.2), bottom=c(-0.1,-1), bottomleft=c(1.1,1.2), 
                                                left=c(1.1,-1), topleft=c(1,-0.2))
                                text( x.x[l.pos], s.x[l.pos], labels=sprintf(format.clab, probs[match(r, radius)]), 
                                                adj=l.adj, cex=cex.label, col=col.clab)
                        }
                }
        }
        
        par(xpd=xpd)
        
        usr<-par("usr")
        
        segments(c(usr[1],0),c(0,0 ), c(usr[2],0), c(0,usr[4]))
        
        axis(1, pos=0, las=las, cex.axis=cex.axis, at=at.xax)
        
        if(is.null(at.yax)) {
                at.yax <- pretty(c(0, usr[4]), n=4)
                at.yax <- at.yax[at.yax > 0]
        } 
        axis(2, pos=0, las=las, cex.axis=cex.axis, at=at.yax)
        
        points(z,f,pch=pch, col=col, bg=bg, cex=cex)
        
        if(!all(is.na(labels))) text(z, f, as.character(labels), pos=pos, adj=adj, cex=cex.label)
        
        if(!missing(main)) title(main=main)
        if(!missing(ylab)) title(ylab=ylab)
        if(!missing(xlab)) {
                xpd <- par(xpd=TRUE)
                h <- strheight("xlab", units = "user") #uses constant string to keep constant distance
                text(0, -h*5, labels=xlab)
                par(xpd=xpd)
        }
        
        return(invisible(list(mu=mu, sigma=sigma, s0=s0)))

}

