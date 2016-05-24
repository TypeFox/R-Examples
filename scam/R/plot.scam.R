
## PLOT FUNCTION based on plot.gam(mgcv)....


plot.scam <- function(x,residuals=FALSE,rug=TRUE,se=TRUE,pages=0,select=NULL,scale=-1,n=100,n2=40,
                     pers=FALSE,theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                     ylim=NULL,xlim=NULL,too.far=0.1,all.terms=FALSE,shade=FALSE,shade.col="gray80",
                shift=0,trans=I,seWithMean=FALSE,unconditional=FALSE,by.resids=FALSE,scheme=0,...)

# Create an appropriate plot for each smooth term of a gam.....
# x is a scam object
# rug determines whether a rug plot should be added to each plot
# se determines whether twice standard error bars are to be added
# pages is the number of pages over which to split output - 0 implies that 
# graphic settings should not be changed for plotting
# scale -1 for same y scale for each plot
#        0 for different y scales for each plot
# n - number of x axis points to use for plotting each term
# n2 is the square root of the number of grid points to use for contouring
# 2-d terms.

{ ######################################
  ## Local function for producing labels
  ######################################

  sub.edf <- function(lab,edf) {
    ## local function to substitute edf into brackets of label
    ## labels are e.g. smooth[[1]]$label
    pos <- regexpr(":",lab)[1]
    if (pos<0) { ## there is no by variable stuff
      pos <- nchar(lab) - 1
      lab <- paste(substr(lab,start=1,stop=pos),",",round(edf,digits=2),")",sep="")
    } else {
      lab1 <- substr(lab,start=1,stop=pos-2)
      lab2 <- substr(lab,start=pos-1,stop=nchar(lab))
      lab <- paste(lab1,",",round(edf,digits=2),lab2,sep="")
    }
    lab
  } ## end of sub.edf


  #########################
  ## start of main function
  #########################

  if (unconditional) {
    if (is.null(x$Vc)) warning("Smoothness uncertainty corrected covariance not available") else 
    x$Vp <- x$Vc ## cov matrix reset to full Bayesian
  }

  w.resid<-NULL
  if (length(residuals)>1) # residuals supplied 
  { if (length(residuals)==length(x$residuals)) 
    w.resid <- residuals else
    warning("residuals argument to plot.gam is wrong length: ignored")
    partial.resids <- TRUE
  } else partial.resids <- residuals # use working residuals or none

  m <- length(x$smooth) ## number of smooth terms

  if (length(scheme)==1) scheme <- rep(scheme,m)
  if (length(scheme)!=m) { 
    warn <- paste("scheme should be a single number, or a vector with",m,"elements")
    warning(warn)
    scheme <- rep(scheme[1],m)
  }

  ## array giving order of each parametric term...
  order <- if (is.list(x$pterms))  unlist(lapply(x$pterms,attr,"order")) else attr(x$pterms,"order")

  if (all.terms) # plot parametric terms as well
  n.para <- sum(order==1) # plotable parametric terms   
  else n.para <- 0 
 
  if (se) ## sort out CI widths for 1 and 2D
  { if (is.numeric(se)) se2.mult <- se1.mult <- se else { se1.mult <- 2;se2.mult <- 1} 
    if (se1.mult<0) se1.mult<-0;if (se2.mult < 0) se2.mult <- 0
  } else se1.mult <- se2.mult <-1
  
  if (se && x$Vp[1,1] < 0) ## check that variances are actually available
  { se <- FALSE
    warning("No variance estimates available")
  }

  if (partial.resids) { ## getting information needed for partial residuals...
    if (is.null(w.resid)) { ## produce working resids if info available
      if (is.null(x$residuals)||is.null(x$weights)) partial.resids <- FALSE else {
        wr <- sqrt(x$weights)
        w.resid <- x$residuals*wr/mean(wr) # weighted working residuals
      }
    }
    if (partial.resids) fv.terms <- predict(x,type="terms") ## get individual smooth effects
  }

  pd <- list(); ## plot data list
  i <- 1 # needs a value if no smooths, but parametric terms ...

  ##################################################
  ## First the loop to get the data for the plots...
  ##################################################

  if (m>0) for (i in 1:m) { ## work through smooth terms
    first <- x$smooth[[i]]$first.para
    last <- x$smooth[[i]]$last.para
    edf <- sum(x$edf[first:last]) ## Effective DoF for this term
    term.lab <- sub.edf(x$smooth[[i]]$label,edf)
    #P <- plot(x$smooth[[i]],P=NULL,data=x$model,n=n,n2=n2,xlab=xlab,ylab=ylab,too.far=too.far,label=term.lab,
    #          se1.mult=se1.mult,se2.mult=se2.mult,xlim=xlim,ylim=ylim,main=main,scheme=scheme[i],...)
    attr(x$smooth[[i]],"coefficients") <- x$coefficients[first:last]   ## relevent coefficients
    P <- plot(x$smooth[[i]],P=NULL,data=x$model,partial.resids=partial.resids,rug=rug,se=se,scale=scale,n=n,n2=n2,
                     pers=pers,theta=theta,phi=phi,jit=jit,xlab=xlab,ylab=ylab,main=main,label=term.lab,
                     ylim=ylim,xlim=xlim,too.far=too.far,shade=shade,shade.col=shade.col,
                     se1.mult=se1.mult,se2.mult=se2.mult,shift=shift,trans=trans,
                     by.resids=by.resids,scheme=scheme[i],...)

    if (is.null(P)) pd[[i]] <- list(plot.me=FALSE) else if (is.null(P$fit)) {
      p <- x$coefficients[first:last]   ## relevent coefficients 
      offset <- attr(P$X,"offset")      ## any term specific offset
      ## get fitted values ....
### addition for univariate shape-constrained smooths ...
     if (inherits(x$smooth[[i]], c("mpi.smooth","mpd.smooth", "cv.smooth", "cx.smooth",      "mdcv.smooth","mdcx.smooth","micv.smooth","micx.smooth"))){
                        q <- ncol(P$X)
                        beta.c <- c(0,exp(p))
                        fit.c <- P$X%*%beta.c # fitted values for the SCOP-splines identifiability constraints
                          # get an intercept that is a difference between the fit with SCOP-spline constraints and with the centering constraint...
                        intercept <- -sum(fit.c)/n 
                        onet <- matrix(rep(1,n),1,n)
                        A <- onet%*%P$X 
                        qrX <- qr(P$X)
                        R <- qr.R(qrX) 
                        qrA <- qr(t(A))
                        R <- R[-1,]
                        RZa <- t(qr.qty(qrA,t(R)))[,2:q] 
                        RZa.inv <- solve(RZa)
                        RZaR <- RZa.inv%*%R
                        beta.a <- RZaR%*%beta.c
                        p <- c(0,beta.a)
                        p <- qr.qy(qrA,p)
                }
### addition for double and single monotonicity...
        if (inherits(x$smooth[[i]], c("tedmi.smooth","tedmd.smooth",                     "tesmi1.smooth","tesmi2.smooth","tesmd1.smooth","tesmd2.smooth")))
                   { p.ident <- x$p.ident[first:last]
                     ii <- p.ident == 1
                     count <- sum(ii) ## the number of exponentiated parameters 
                     iv<-array(0, dim=c(count,1)) # index vector for the exponent. parameters
                     k<-1
                     for (j in 1:length(p))
                        {if (p.ident[j]==1) {iv[k]<-j; k<-k+1}}
                     p[iv] <- exp(p[iv])
                     beta <- x$smooth[[i]]$Zc%*%p
                     fit.c <- P$X%*%beta # fitted values for the SCOP-spline identifiability constraints
                     # get an intercept as a difference between the fit with SCOP-spline constraints and fit with the centering constraint...
                     intercept <- -sum(fit.c)/length(fit.c)
                     onet <- matrix(rep(1,nrow(P$X)),1,nrow(P$X))
                     A <- onet%*%P$X
                     qrX <- qr(P$X)
                     R <- qr.R(qrX)
                     qrA <- qr(t(A))
                     q <- ncol(P$X)
                     if (inherits(x$smooth[[i]], c("tedmi.smooth","tedmd.smooth")))
                     { # get `beta.a' for double monotonicity...
                         R <- R[-1,]
                         RZa <- t(qr.qty(qrA,t(R)))[,2:q] 
                         RZa.inv <- solve(RZa)
                         RZaR <- RZa.inv%*%R
                         beta.a <- RZaR%*%beta
                     }
                     else # get `beta.a' for single monotonicity...
                     { RZa <- t(qr.qty(qrA,t(R)))[,2:q]
                       RZatRZa.inv <- solve(t(RZa)%*%RZa) 
                       Q <- qr.Q(qrX)
                       B1 <- RZatRZa.inv%*%t(RZa)
                       RZaR <- B1%*%R
                       beta.a <- RZaR%*%beta + B1%*%(intercept*colSums(Q))
                     }
                     p <- c(0,beta.a)
                     p <- qr.qy(qrA,p)
                }
## end shape constrained supplement...

      if (is.null(offset)) P$fit <- P$X%*%p else P$fit <- P$X%*%p + offset 
      if (!is.null(P$exclude)) P$fit[P$exclude] <- NA
      if (se && P$se) { ## get standard errors for fit
## shape constrained supplement...
        if (inherits(x$smooth[[i]], c("mpi.smooth","mpd.smooth", "cv.smooth", "cx.smooth",                    "mdcv.smooth","mdcx.smooth","micv.smooth","micx.smooth"))){
                        XZa <- t(qr.qty(qrA,t(P$X)))[,2:q]
                        Ga <- XZa%*%RZaR
                        Vp <- x$Vp.t[c(1,first:last),c(1,first:last)] 
                        Vp.c <- Vp
                        Vp.c[,1] <- rep(0,nrow(Vp))
                        Vp.c[1,] <- rep(0,ncol(Vp))
                        se.fit <- sqrt(rowSums((Ga%*%Vp.c)*Ga))
        } else if (inherits(x$smooth[[i]], c("tedmi.smooth","tedmd.smooth",                     "tesmi1.smooth","tesmi2.smooth","tesmd1.smooth","tesmd2.smooth"))) {
                          XZa <- t(qr.qty(qrA,t(P$X)))[,2:ncol(P$X)]
                          Ga <- XZa%*%RZaR%*%x$smooth[[i]]$Zc
                          Vp <- x$Vp.t[first:last,first:last] 
                          se.fit <- rowSums((Ga%*%Vp)*Ga)^.5
        } else {   ## end scam supplement      
        ## test whether mean variability to be added to variability (only for centred terms)
        if (seWithMean && attr(x$smooth[[i]],"nCons")>0) {
          if (length(x$cmX) < ncol(x$Vp)) x$cmX <- c(x$cmX,rep(0,ncol(x$Vp)-length(x$cmX)))
          X1 <- matrix(x$cmX,nrow(P$X),ncol(x$Vp),byrow=TRUE)
          meanL1 <- x$smooth[[i]]$meanL1
          if (!is.null(meanL1)) X1 <- X1 / meanL1
          X1[,first:last] <- P$X
          se.fit <- sqrt(pmax(0,rowSums((X1%*%x$Vp)*X1)))
        } else se.fit <- ## se in centred (or anyway unconstained) space only
        sqrt(pmax(0,rowSums((P$X%*%x$Vp[first:last,first:last,drop=FALSE])*P$X)))
        if (!is.null(P$exclude)) P$se.fit[P$exclude] <- NA
        } 
      } ## standard errors for fit completed
      if (partial.resids) { P$p.resid <- fv.terms[,length(order)+i] + w.resid }
      if (se && P$se) P$se <- se.fit*P$se.mult  # Note multiplier
      P$X <- NULL
      P$plot.me <- TRUE
      pd[[i]] <- P;rm(P) 
    } else { ## P$fit created directly
      if (partial.resids) { P$p.resid <- fv.terms[,length(order)+i] + w.resid }
      P$plot.me <- TRUE
      pd[[i]] <- P;rm(P)
    }
  } ## end of data setup loop through smooths

  
  ##############################################
  ## sort out number of pages and plots per page 
  ##############################################

  n.plots <- n.para
  if (m>0) for (i in 1:m) n.plots <- n.plots + as.numeric(pd[[i]]$plot.me) 

  if (n.plots==0) stop("No terms to plot - nothing for plot.gam() to do.")

  if (pages>n.plots) pages<-n.plots
  if (pages<0) pages<-0
  if (pages!=0)    # figure out how to display things
  { ppp<-n.plots%/%pages
    if (n.plots%%pages!=0) 
    { ppp<-ppp+1
      while (ppp*(pages-1)>=n.plots) pages<-pages-1
    } 

    # now figure out number of rows and columns
    c <- r <- trunc(sqrt(ppp))
    if (c<1) r <- c <- 1
    if (c*r < ppp) c <- c + 1
    if (c*r < ppp) r <- r + 1  
    oldpar<-par(mfrow=c(r,c))
  
  } else
  { ppp<-1;oldpar<-par()}
  
  if ((pages==0&&prod(par("mfcol"))<n.plots&&dev.interactive())||
       pages>1&&dev.interactive()) ask <- TRUE else ask <- FALSE 
  
  if (!is.null(select)) {
    ask <- FALSE
  }
 
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  #####################################
  ## get a common scale, if required...
  #####################################

  if (scale==-1&&is.null(ylim)) {
    k <- 0
    if (m>0) for (i in 1:m) if (pd[[i]]$plot.me&&pd[[i]]$scale) { ## loop through plot data 
      if (se&&length(pd[[i]]$se)>1) { ## require CIs on plots
        ul<-pd[[i]]$fit+pd[[i]]$se
        ll<-pd[[i]]$fit-pd[[i]]$se
        if (k==0) { 
          ylim <- c(min(ll,na.rm=TRUE),max(ul,na.rm=TRUE));k <- 1
        } else {
          if (min(ll,na.rm=TRUE)<ylim[1]) ylim[1] <- min(ll,na.rm=TRUE)
	  if (max(ul,na.rm=TRUE)>ylim[2]) ylim[2] <- max(ul,na.rm=TRUE)
        }
      } else { ## no standard errors
        if (k==0) {
          ylim <- range(pd[[i]]$fit,na.rm=TRUE);k <- 1
        } else {
          if (min(pd[[i]]$fit,na.rm=TRUE)<ylim[1]) ylim[1] <- min(pd[[i]]$fit,na.rm=TRUE)
          if (max(pd[[i]]$fit,na.rm=TRUE)>ylim[2]) ylim[2] <- max(pd[[i]]$fit,na.rm=TRUE)
        }
      }
      if (partial.resids) { 
        ul <- max(pd[[i]]$p.resid,na.rm=TRUE)
        if (ul > ylim[2]) ylim[2] <- ul
        ll <-  min(pd[[i]]$p.resid,na.rm=TRUE)
        if (ll < ylim[1]) ylim[1] <- ll
      } ## partial resids done
    } ## loop end 
  } ## end of common scale computation
  
  ##############################################################
  ## now plot smooths, by calling plot methods with plot data...
  ##############################################################

  if (m>0) for (i in 1:m) if (pd[[i]]$plot.me&&(is.null(select)||i==select)) {
    plot(x$smooth[[i]],P=pd[[i]],partial.resids=partial.resids,rug=rug,se=se,scale=scale,n=n,n2=n2,
                     pers=pers,theta=theta,phi=phi,jit=jit,xlab=xlab,ylab=ylab,main=main,
                     ylim=ylim,xlim=xlim,too.far=too.far,shade=shade,shade.col=shade.col,
                     shift=shift,trans=trans,by.resids=by.resids,scheme=scheme[i],...)

  } ## end of smooth plotting loop
  
  ####################################################
  ## Finally deal with any parametric term plotting...
  ####################################################

  if (n.para>0) # plot parameteric terms
  { class(x) <- c("scam", "gam","glm","lm") # needed to get termplot to call model.frame.glm 
    if (is.null(select)) {
      attr(x,"para.only") <- TRUE
      termplot(x,se=se,rug=rug,col.se=1,col.term=1,main=attr(x$pterms,"term.labels"),...)
    } else { # figure out which plot is required
      if (select > m) { 
        ## can't figure out how to get this to work with more than first linear predictor
        ## as termplots relies on matching terms to names in original data... 
        select <- select - m # i.e. which parametric term
        term.labels <- attr(x$pterms,"term.labels")
        term.labels <- term.labels[order==1]
        if (select <= length(term.labels)) {
          # if (interactive() && m &&i%%ppp==0) 
          termplot(x,terms=term.labels[select],se=se,rug=rug,col.se=1,col.term=1,...)
        }  
      }
    }
  }
  if (pages>0) par(oldpar)
} ## end plot.scam


#############################################
## below is needed and copied from plots.r of mgcv...
##########################

plot.random.effect <- function(x,P=NULL,data=NULL,label="",se1.mult=1,se2.mult=2,
                     partial.resids=FALSE,rug=TRUE,se=TRUE,scale=-1,n=100,n2=40,
                     pers=FALSE,theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                     ylim=NULL,xlim=NULL,too.far=0.1,shade=FALSE,shade.col="gray80",
                     shift=0,trans=I,by.resids=FALSE,scheme=0,...) {
## plot method for a "random.effect" smooth class
 
  if (is.null(P)) { ## get plotting information...
    if (!x$plot.me) return(NULL) else { ## shouldn't or can't plot 
      raw <- data[x$term][[1]]
      p <- x$last.para - x$first.para + 1
      X <- diag(p)   # prediction matrix for this term
      if (is.null(xlab)) xlabel<- "Gaussian quantiles" else xlabel <- xlab
      if (is.null(ylab)) ylabel <- "effects" else ylabel <- ylab
      if (!is.null(main)) label <- main
      return(list(X=X,scale=FALSE,se=FALSE,raw=raw,xlab=xlabel,ylab=ylabel,
             main=label))

    } ## end of basic plot data production 
  } else { ## produce plot
    qqnorm(P$fit,main=P$main,xlab=P$xlab,ylab=P$ylab,...)
    qqline(P$fit)
  } ## end of plot production
} ## end of plot.random.effect


plot.mgcv.smooth <- function(x,P=NULL,data=NULL,label="",se1.mult=1,se2.mult=2,
                     partial.resids=FALSE,rug=TRUE,se=TRUE,scale=-1,n=100,n2=40,
                     pers=FALSE,theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
                     ylim=NULL,xlim=NULL,too.far=0.1,shade=FALSE,shade.col="gray80",
                     shift=0,trans=I,by.resids=FALSE,scheme=0,...) {
## default plot method for smooth objects `x' inheriting from "mgcv.smooth"
## `x' is a smooth object, usually part of a `gam' fit. It has an attribute
##     'coefficients' containg the coefs for the smooth, but usually these
##     are not needed.
## `P' is a list of plot data. 
##     If `P' is NULL then the routine should compute some of this plot data
##     and return without plotting...  
##     * X the matrix mapping the smooth's coefficients to the values at
##         which the smooth must be computed for plotting.
##     * The values against which to plot.
##     * `exclude' indicates rows of X%*%p to set to NA for plotting -- NULL for none.
##     * se TRUE if plotting of the term can use standard error information.
##     * scale TRUE if the term should be considered by plot.gam if a common
##             y scale is required.
##     * any raw data information.
##     * axis labels and plot titles 
##     As an alternative, P may contain a 'fit' field directly, in which case the 
##     very little processing is done outside the routine, except for partial residual
##     computations.
##     Alternatively return P as NULL if x should not be plotted.
##     If P is not NULL it will contain 
##     * fit - the values for plotting 
##     * se.fit - standard errors of fit (can be NULL)
##     * the values against which to plot
##     * any raw data information
##     * any partial.residuals 
## `data' is a data frame containing the raw data for the smooth, usually the 
##        model.frame of the fitted gam. Can be NULL if P is not NULL.
## `label' is the term label, usually something like e.g. `s(x,12.34)'.
#############################

  sp.contour <- function(x,y,z,zse,xlab="",ylab="",zlab="",titleOnly=FALSE,
               se.plot=TRUE,se.mult=1,trans=I,shift=0,...)   
  ## function for contouring 2-d smooths with 1 s.e. limits
  { gap<-median(zse,na.rm=TRUE)  
    zr<-max(trans(z+zse+shift),na.rm=TRUE)-min(trans(z-zse+shift),na.rm=TRUE) # plotting range  
    n<-10  
    while (n>1 && zr/n<2.5*gap) n<-n-1    
    zrange<-c(min(trans(z-zse+shift),na.rm=TRUE),max(trans(z+zse+shift),na.rm=TRUE))  
    zlev<-pretty(zrange,n)  ## ignore codetools on this one  
    yrange<-range(y);yr<-yrange[2]-yrange[1]  
    xrange<-range(x);xr<-xrange[2]-xrange[1]  
    ypos<-yrange[2]+yr/10
    args <- as.list(substitute(list(...)))[-1]
    args$x <- substitute(x);args$y <- substitute(y)
    args$type="n";args$xlab<-args$ylab<-"";args$axes<-FALSE
    do.call("plot",args)

    cs<-(yr/10)/strheight(zlab);if (cs>1) cs<-1 # text scaling based on height  
  
    tl<-strwidth(zlab);  
    if (tl*cs>3*xr/10) cs<-(3*xr/10)/tl  
    args <- as.list(substitute(list(...)))[-1]
    n.args <- names(args)
    zz <- trans(z+shift) ## ignore codetools for this
    args$x<-substitute(x);args$y<-substitute(y);args$z<-substitute(zz)
    if (!"levels"%in%n.args) args$levels<-substitute(zlev)
    if (!"lwd"%in%n.args) args$lwd<-2
    if (!"labcex"%in%n.args) args$labcex<-cs*.65
    if (!"axes"%in%n.args) args$axes <- FALSE
    if (!"add"%in%n.args) args$add <- TRUE
    do.call("contour",args)
  
    if (is.null(args$cex.main)) cm <- 1 else cm <- args$cex.main
    if (titleOnly)  title(zlab,cex.main=cm) else 
    { xpos<-xrange[1]+3*xr/10  
      xl<-c(xpos,xpos+xr/10); yl<-c(ypos,ypos)   
      lines(xl,yl,xpd=TRUE,lwd=args$lwd)  
      text(xpos+xr/10,ypos,zlab,xpd=TRUE,pos=4,cex=cs*cm,off=0.5*cs*cm)  
    }
    if  (is.null(args$cex.axis)) cma <- 1 else cma <- args$cex.axis
    axis(1,cex.axis=cs*cma);axis(2,cex.axis=cs*cma);box();
    if  (is.null(args$cex.lab)) cma <- 1 else cma <- args$cex.lab  
    mtext(xlab,1,2.5,cex=cs*cma);mtext(ylab,2,2.5,cex=cs*cma)  
    if (!"lwd"%in%n.args) args$lwd<-1
    if (!"lty"%in%n.args) args$lty<-2
    if (!"col"%in%n.args) args$col<-2
    if (!"labcex"%in%n.args) args$labcex<-cs*.5
    zz <- trans(z+zse+shift)
    args$z<-substitute(zz)

    do.call("contour",args)

    if (!titleOnly) {
      xpos<-xrange[1]  
      xl<-c(xpos,xpos+xr/10)#;yl<-c(ypos,ypos)  
      lines(xl,yl,xpd=TRUE,lty=args$lty,col=args$col)  
      text(xpos+xr/10,ypos,paste("-",round(se.mult),"se",sep=""),xpd=TRUE,pos=4,cex=cs*cm,off=0.5*cs*cm)  
    }

    if (!"lty"%in%n.args) args$lty<-3
    if (!"col"%in%n.args) args$col<-3
    zz <- trans(z - zse+shift)
    args$z<-substitute(zz)
    do.call("contour",args)
    
    if (!titleOnly) {
      xpos<-xrange[2]-xr/5  
      xl<-c(xpos,xpos+xr/10);  
      lines(xl,yl,xpd=TRUE,lty=args$lty,col=args$col)  
      text(xpos+xr/10,ypos,paste("+",round(se.mult),"se",sep=""),xpd=TRUE,pos=4,cex=cs*cm,off=0.5*cs*cm)  
    }
  }  ## end of sp.contour

  if (is.null(P)) { ## get plotting information...
    if (!x$plot.me||x$dim>2) return(NULL) ## shouldn't or can't plot
    if (x$dim==1) { ## get basic plotting data for 1D terms 
      raw <- data[x$term][[1]]
      if (is.null(xlim)) xx <- seq(min(raw),max(raw),length=n) else # generate x sequence for prediction
      xx <- seq(xlim[1],xlim[2],length=n)
      if (x$by!="NA")         # deal with any by variables
      { by<-rep(1,n);dat<-data.frame(x=xx,by=by)
        names(dat)<-c(x$term,x$by)
      } else { 
        dat<-data.frame(x=xx);names(dat) <- x$term
      } ## prediction data.frame finished
      X <- PredictMat(x,dat)   # prediction matrix for this term
      if (is.null(xlab)) xlabel<- x$term else xlabel <- xlab
      if (is.null(ylab)) ylabel <- label else ylabel <- ylab
      if (is.null(xlim)) xlim <- range(xx)
      return(list(X=X,x=xx,scale=TRUE,se=TRUE,raw=raw,xlab=xlabel,ylab=ylabel,
             main=main,se.mult=se1.mult,xlim=xlim))
    } else { ## basic plot data for 2D terms
      xterm <- x$term[1]
      if (is.null(xlab)) xlabel <- xterm else xlabel <- xlab
      yterm <- x$term[2]
      if (is.null(ylab)) ylabel <- yterm else ylabel <- ylab
      raw <- data.frame(x=as.numeric(data[xterm][[1]]),
                        y=as.numeric(data[yterm][[1]]))
      n2 <- max(10,n2)
      if (is.null(xlim)) xm <- seq(min(raw$x),max(raw$x),length=n2) else 
        xm <- seq(xlim[1],xlim[2],length=n2)
      if (is.null(ylim)) ym <- seq(min(raw$y),max(raw$y),length=n2) else
        ym <- seq(ylim[1],ylim[2],length=n2)
      xx <- rep(xm,n2)
      yy <- rep(ym,rep(n2,n2))
      if (too.far>0)
      exclude <- exclude.too.far(xx,yy,raw$x,raw$y,dist=too.far) else
      exclude <- rep(FALSE,n2*n2)
      if (x$by!="NA")         # deal with any by variables
      { by <- rep(1,n2^2);dat <- data.frame(x=xx,y=yy,by=by)
        names(dat) <- c(xterm,yterm,x$by)
      } else { 
        dat<-data.frame(x=xx,y=yy);names(dat)<-c(xterm,yterm)
      }  ## prediction data.frame complete
      X <- PredictMat(x,dat)   ## prediction matrix for this term
      if (is.null(main)) { 
        main <- label
      }
      if (is.null(ylim)) ylim <- range(ym) 
      if (is.null(xlim)) xlim <- range(xm) 
      return(list(X=X,x=xm,y=ym,scale=FALSE,se=TRUE,raw=raw,xlab=xlabel,ylab=ylabel,
             main=main,se.mult=se2.mult,ylim=ylim,xlim=xlim,exclude=exclude))
    } ## end of 2D basic plot data production 
  } else { ## produce plot
    if (se) { ## produce CI's
      if (x$dim==1) { 
        if (scheme == 1) shade <- TRUE
        ul <- P$fit + P$se ## upper CL
        ll <- P$fit - P$se ## lower CL
        if (scale==0&&is.null(ylim)) { ## get scale 
          ylimit<-c(min(ll),max(ul))
          if (partial.resids) { 
            max.r <- max(P$p.resid,na.rm=TRUE)
            if ( max.r> ylimit[2]) ylimit[2] <- max.r
            min.r <-  min(P$p.resid,na.rm=TRUE)
            if (min.r < ylimit[1]) ylimit[1] <- min.r
          }
        }
        if (!is.null(ylim)) ylimit <- ylim
         
        ## plot the smooth... 
        if (shade) { 
          plot(P$x,trans(P$fit+shift),type="n",xlab=P$xlab,ylim=trans(ylimit+shift),
                 xlim=P$xlim,ylab=P$ylab,main=P$main,...)
          polygon(c(P$x,P$x[n:1],P$x[1]),
                    trans(c(ul,ll[n:1],ul[1])+shift),col = shade.col,border = NA)
          lines(P$x,trans(P$fit+shift),...)
        } else { ## ordinary plot 
          plot(P$x,trans(P$fit+shift),type="l",xlab=P$xlab,ylim=trans(ylimit+shift),xlim=P$xlim,
                 ylab=P$ylab,main=P$main,...)
          if (is.null(list(...)[["lty"]])) { 
            lines(P$x,trans(ul+shift),lty=2,...)
            lines(P$x,trans(ll+shift),lty=2,...)
          } else { 
            lines(P$x,trans(ul+shift),...)
            lines(P$x,trans(ll+shift),...)
          }
        } ## ... smooth plotted
       
        if (partial.resids&&(by.resids||x$by=="NA")) { ## add any partial residuals
          if (length(P$raw)==length(P$p.resid)) {
            if (is.null(list(...)[["pch"]]))
            points(P$raw,trans(P$p.resid+shift),pch=".",...) else
            points(P$raw,trans(P$p.resid+shift),...) 
          } else {
            warning("Partial residuals do not have a natural x-axis location for linear functional terms")
          }
        } ## partial residuals finished 
	 
        if (rug) { 
          if (jit) rug(jitter(as.numeric(P$raw)),...)
          else rug(as.numeric(P$raw),...)
	} ## rug plot done

      } else if (x$dim==2) { 
        P$fit[P$exclude] <- NA
        if (pers) scheme <- 1
        if (scheme == 1) { ## perspective plot 
          persp(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
                  zlab=P$main,ylim=P$ylim,xlim=P$xlim,theta=theta,phi=phi,...)
        } else if (scheme==2) {
          image(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
                  main=P$main,xlim=P$xlim,ylim=P$ylim,col=heat.colors(50),...)
          contour(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),add=TRUE,col=3,...)
          if (rug) {  
            if (is.null(list(...)[["pch"]])) points(P$raw$x,P$raw$y,pch=".",...) else
            points(P$raw$x,P$raw$y,...)
          }
        } else { ## contour plot with error contours
          sp.contour(P$x,P$y,matrix(P$fit,n2,n2),matrix(P$se,n2,n2),
                     xlab=P$xlab,ylab=P$ylab,zlab=P$main,titleOnly=!is.null(main),
                     se.mult=1,trans=trans,shift=shift,...)
          if (rug) { 
            if (is.null(list(...)[["pch"]]))
            points(P$raw$x,P$raw$y,pch=".",...) else
            points(P$raw$x,P$raw$y,...) 
          }
        } ## counter plot done 
      } else { 
         warning("no automatic plotting for smooths of more than two variables")
      }
    } else { ## no CI's
      if (x$dim==1) { 
        if (scale==0&&is.null(ylim)) { 
          if (partial.resids) ylimit <- range(P$p.resid,na.rm=TRUE) else ylimit <-range(P$fit)
        }
        if (!is.null(ylim)) ylimit <- ylim
        plot(P$x,trans(P$fit+shift),type="l",xlab=P$xlab,
             ylab=P$ylab,ylim=trans(ylimit+shift),xlim=P$xlim,main=P$main,...)
        if (rug) { 
          if (jit) rug(jitter(as.numeric(P$raw)),...)
          else rug(as.numeric(P$raw),...) 
        }
        if (partial.resids&&(by.resids||x$by=="NA")) { 
          if (is.null(list(...)[["pch"]]))
          points(P$raw,trans(P$p.resid+shift),pch=".",...) else
          points(P$raw,trans(P$p.resid+shift),...)
        }
      } else if (x$dim==2) { 
        P$fit[P$exclude] <- NA
        if (!is.null(main)) P$title <- main
        if (pers) scheme <- 1
        if (scheme==1) { 
          persp(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
                          zlab=P$main,theta=theta,phi=phi,xlim=P$xlim,ylim=P$ylim,...)
        } else if (scheme==2) {
          image(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
                  main=P$main,xlim=P$xlim,ylim=P$ylim,col=heat.colors(50),...)
          contour(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),add=TRUE,col=3,...)
          if (rug) {  
            if (is.null(list(...)[["pch"]])) points(P$raw$x,P$raw$y,pch=".",...) else
            points(P$raw$x,P$raw$y,...)
          }
        } else { 
          contour(P$x,P$y,matrix(trans(P$fit+shift),n2,n2),xlab=P$xlab,ylab=P$ylab,
                  main=P$main,xlim=P$xlim,ylim=P$ylim,...)
          if (rug) {  
            if (is.null(list(...)[["pch"]])) points(P$raw$x,P$raw$y,pch=".",...) else
            points(P$raw$x,P$raw$y,...)
          }
        }  
      } else { 
        warning("no automatic plotting for smooths of more than one variable")
      }
    } ## end of no CI code
  } ## end of plot production
}



