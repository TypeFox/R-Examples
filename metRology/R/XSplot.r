xs.plot<-function(x, ...) UseMethod("xs.plot")

#Provide a camel-friendly alias
XSplot <- xs.plot


xs.plot.default<-function(x,g,s, degfree, labels.arg=NA, mu, sigma, 
        probs=c(0.5, 0.95, 0.99), basis=c("robust","classical"), method=c("chisq","density"),
        main=paste("X-S plot -", basis, "basis"), xlab=deparse(substitute(x)), ylab, 
        contours=TRUE, col.contours="lightgrey", lty.contours=par("lty"), 
        lwd.contours=par("lwd"),
        label.contours=contours, format.clab="p=%3.2f",
        pos.clab="bottomright", col.clab=col.contours, cex.clab=0.7,
        cex.label=0.7, pos=3, adj=NULL, 
        pch=par("pch"), col=par("col"), bg=par("bg"), cex=par("cex"), 
        add=FALSE, ...) {
        #Produces a plot of observed values vs standard deviations
        #as described in ISO 13528:2005
        
        if(missing(g) && missing(s)) stop("One of g and s must be specified")
        if(!missing(g) && !missing(s)) warning("Only one of g and s should be specified: ignoring g")

        basis<-match.arg(basis)
        method<-match.arg(method)
    
        pos.clab <- match.arg(pos.clab, choices=c("top", "topright", "right",
                        "bottomright", "bottom", "bottomleft", "left", "topleft"), several.ok=TRUE)
                        
        if(missing(ylab)) {
                ylab <- if(missing(s))
                                sprintf("s( %s | %s )",deparse(substitute(x)) , deparse(substitute(g)))
                        else
                                deparse(substitute(s)) 
        }

        if(!missing(s)) {
                #s has been provided
                #Check that df is also there and that x and s are the same length:
                if(missing(degfree)) stop("degfree must be specified with s")
                if(length(s) != length(x) ) stop("x and s are different lengths")
                #x must be single values, so set z to x
                z <- x
                n <- degfree+1
                ni <- rep(n, length(s))
        } else {
                #Using g
                z <- tapply(x, g, mean, na.rm=T)
                s <- tapply(x, g, sd, na.rm=T)
                ni <- tapply(x, g, function(x) sum( !is.na(x)) )
                n <- median(ni, na.rm=T) 
                degfree <- n-1
        }
        
        if(basis=="robust") {
                loc.scale <- algA(z, na.rm=TRUE)
        } else {
                loc.scale <- list(mu=mean(na.omit(z)), s=sd(na.omit(z)))
        }
        
        if(missing(mu)) {
                        mu <- loc.scale$mu
        }
        
        if(missing(sigma)) {
                if(basis=="robust") 
                        sigma <- algS(s, degfree, na.rm = FALSE)
                else {
                        real.s <- !is.na(s)
                        dfi <- ni-1
                        sigma <- sqrt(sum(dfi[real.s]*s[real.s]^2)/sum(dfi))
                }
        }        
        

        
        #Calculate contour coordinates
        clist <- list()
        if(contours) {
             if(method=="chisq") {
                for(cli in 1:length(probs)) {
                        p <- probs[cli]
                        #Generate xmin...xmax
                        xmin <- mu-sigma*sqrt(qchisq(p,2)/(n))
                        xmax <- mu+sigma*sqrt(qchisq(p,2)/(n))

                        #Generate sine-spaced sequence for x (produces a smoother contour)
                        b.x <- mu+0.5*(xmax-xmin)*sin(seq(-pi/2, pi/2, length.out=101))

                        #get s
                        b.s.upper <- sigma * exp( (1/sqrt(2*degfree)) * sqrt( pmax(0, qchisq(p, 2) - n * ( ( b.x - mu  ) / sigma )^2 )) )
                        b.s.lower <- sigma * exp( -(1/sqrt(2*degfree)) * sqrt( pmax(0, qchisq(p, 2) - n * ( ( b.x - mu  ) / sigma )^2 )) )

                        clist[[cli]]<- data.frame(b.x=c(b.x, rev(b.x[-c(1, length(b.x))])),
                                                        b.s=c(b.s.lower, rev(b.s.upper[-c(1,length(b.s.upper))]))
                                                  )
                }
            } else {
               #method == "density"
               smax<-sqrt((n-2)/(n))
               dhmax<-.dhelmert(n=n, u=0, s=smax, sigma=1)
               for(cli in 1:length(probs)) {
                        p <- probs[cli]
                        #Generate smin...smax
                        smin<-if(n<=2) 0 else .qdhelmert(n=n, d=dhmax*(1-p), u=0, lower.tail=TRUE)
                        smax<-.qdhelmert(n=n, d=dhmax*(1-p), u=0, lower.tail=FALSE)
                        #Generate sine-spaced sequence (produces a smoother contour)
                        b.s <- 0.5*(smin+smax)+0.5*(smax-smin)*sin(seq(-pi/2, pi/2, length.out=101))
                        #solve for x, omitting the x=0 cases as they cause uniroot problems.
                        b.x <- .qdhelmert(n=n, d=dhmax*(1-p), s=b.s[c(-1,-length(b.s))] )
                        
                        #Complete the ellipsoid
                        b.s<-c(b.s, rev(b.s)[c(-1,-length(b.s))])
                        b.x<-c(0,b.x, 0,-rev(b.x))
                        
                        #Scale for x and s:
                        b.s <- b.s*sigma
                        b.x <- mu + b.x*sigma
                        #Rearrange to give the same plotting order as 'chisq' contour
                        b.s <- c(b.s[149:200], b.s[1:148])
                        b.x <- c(b.x[149:200], b.x[1:148])
                        clist[[cli]]<- data.frame(b.x=b.x, b.s=b.s)
               }
            
            }
            #Get the index for the biggest probability (with the largest contour)
            clist.max <- which.max(probs)

        } else
                clist.max <- NA
        
        #Simple plot setup...
        
        if(!add) {
                if(contours) 
                        plot(c(clist[[clist.max]]$b.x, z), c(clist[[clist.max]]$b.s, s), type="n", 
                                main=main, xlab=xlab, ylab=ylab,  ...)
                else 
                        plot(z, s, type="n", main=main, xlab=xlab, ylab=ylab,  ...)
        }
        
        #Plot contours
        if(contours) {
                col.contours <- rep(col.contours, length.out=length(probs))
                lty.contours <- rep(lty.contours, length.out=length(probs))
                lwd.contours <- rep(lwd.contours, length.out=length(probs))
                col.clab<-rep(col.clab, length.out=length(probs))
                cex.clab<-rep(cex.clab, length.out=length(probs))
                for(cli in 1:length(probs)) {
                        with(clist[[cli]], polygon(b.x,b.s, border=col.contours[cli], lty=lty.contours[cli], lwd=lwd.contours[cli], col=NA))
                        if(label.contours) {
                                lpos.x<-vector()
                                lpos.s<-vector()
                                
                                wmins <- which.min(clist[[cli]]$b.s)
                                wmaxs <- which.max(clist[[cli]]$b.s)
                                
                                lpos.x[c("bottom", "right", "top", "left")] <-
                                         with(clist[[cli]], c(mu, max(b.x), mu,  min(b.x)))
                                lpos.s[c("bottom", "right", "top", "left")] <-
                                        with(clist[[cli]], c(min(b.s), b.s[which.max(b.x)], max(b.s),  b.s[which.max(b.x)]))
                                                #Relies on symmetry to avoid 'first' min problems
                                lpos.s[c("bottomright", "topright")] <-
                                        (lpos.s["right"]+lpos.s[c("bottom", "top")] )/2
                                lpos.x["bottomright"]<-with(clist[[cli]], b.x[wmins:wmaxs][which.min(abs(b.s[wmins:wmaxs] - lpos.s["bottomright"] ))] )
                                lpos.x["topright"]<-with(clist[[cli]], b.x[wmins:wmaxs][which.min(abs(b.s[wmins:wmaxs] - lpos.s["topright"] ))] )
                                
                                lpos.s[c("bottomleft", "topleft")] <- lpos.s[c("bottomright", "topright")]
                                lpos.x[c("bottomleft", "topleft")] <- 2*mu - lpos.x[c("bottomright", "topright")]
                                
                                l.adj <- list(top=c(0.5,-0.5), topright=c(-0.1,-0.2), right=c(-0.1,0.5),
                                                bottomright=c(0,1), bottom=c(0.5,1.3), bottomleft=c(1,1), 
                                                left=c(1.1,0.5), topleft=c(1,-0.2))

                                if(n <= 2 && method=="density") {
                                        #This case has lower bound 0                                        
                                        lpos.s[c("left", "bottomleft", "bottom", "bottomright", "right")] <- 0
                                        lpos.x[c("bottomleft", "bottomright")] <- lpos.x[c("left", "right")]
                                        lpos.x["bottom"] <- NA
                                        if("bottom" %in% pos.clab) pos.clab <- unique(c(pos.clab, c("bottomleft", "bottomright")))
                                        l.adj[["left"]] <- c(1.1,-0.5)
                                        l.adj[["right"]] <- c(-0.1,-0.5)
                                        l.adj[["bottomleft"]] <- c(0.5,1.3)
                                        l.adj[["bottomright"]] <- c(0.5,1.3)
                                }                                       
                                
                                for(pn in pos.clab) {
#                                        text( lpos.x[pn], lpos.s[pn], labels=paste("p=",format(probs[cli], digits=3), sep=""), 
#                                                adj=l.adj[[pn]], cex=cex.clab[cli], col=col.clab[cli])
                                        text( lpos.x[pn], lpos.s[pn], labels=sprintf(format.clab, probs[cli]), 
                                                adj=l.adj[[pn]], cex=cex.clab[cli], col=col.clab[cli])
                                         
                                }
                        }
                }
        }

 

        points(z, s, pch=pch, col=col, bg=bg, cex=cex)
        
        if(!all(is.na(labels.arg))) {
                #labels.arg is not NA; plot labels
                
                if(!is.logical(labels.arg[1])) {
                        labels.arg <- as.character(labels.arg) #Just in case non-character
                } else if(labels.arg[1]) {
                        labels.arg <- if(length(nz<-names(z))) nz 
                                  else paste(1:length(z)) 
                } 
                
                if(is.character(labels.arg[1]) )
                        text(z, s, as.character(labels.arg), pos=pos, adj=adj, cex=cex.label)
        }
        
        return(invisible(list(x=z, y=s, mu=mu, sigma=sigma, clist=clist))) 

}


.dhelmert <- function(n, u=0, s=1, sigma=1) {
        K1<-( n^(n/2) ) / ( 2^((n-2)/2) *sqrt(pi) * gamma((n-1)/2) * sigma^2 )
        K1*(s^(n-2))*exp( (-n/(2*sigma^2))*(s^2+u^2) )
}

.qdhelmert<-function(n,d,u,s,sigma=1, umax=6, smax=6*sigma, lower.tail=TRUE) {
        if(missing(u) && missing(s)) stop("One of u or s must be given")
        find.u <- function(u, n,d,s,sigma) .dhelmert(n=n,u=u,s=s,sigma=sigma )-d
        find.s <- function(s, n,d,u,sigma) .dhelmert(n=n,u=u,s=s,sigma=sigma )-d
        if(missing(u)) {
                #find u
                objective <- find.u
                u <- vector(length=length(s))
                for(i in 1:length(u)) {
                        u[i] <- uniroot(objective, c(0, umax), n=n, s=s[i], d=d, sigma=sigma)$root
                }
                return(u)
        } else if(missing(s)){
                #find s
                s.max<-sigma*sqrt((n-2)/(n))
                #find u
                objective <- find.s
                s <- vector(length=length(u))
                if(lower.tail) interval <- c(0,s.max) 
                        else interval <- c(s.max, smax)
                for(i in 1:length(s)) {
                        s[i] <- uniroot(objective, interval, n=n, u=u[i], d=d, sigma=sigma)$root
                }
                return(s)
        } else {
                stop("Only one of u and s should be given")
        }
}


