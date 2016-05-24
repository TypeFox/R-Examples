plot.flap <- function(x,xlim=NULL,ylim=NULL,lty=NULL,cols=NULL,
                      xlab=NULL,ylab=NULL,main=NULL,main.panel=NULL,
                      groups=NULL,add=FALSE,gloss=FALSE,glind=NULL,
                      extend=0.3,col.gloss=1,cex.gloss=0.8,mfrow=NULL,...) {
#
    qmax <- attr(x,"qmax")
    jmax <- attr(x,"jmax")
    if(is.null(jmax)) jmax <- 1

# If groups is NULL, build it --- with the groups column
# running from 1 to length(x).  Otherwise make sure its entries
# are sensible.
    if(is.null(groups)) {
        groups <- cbind(group=1,as.data.frame(i2qj(1:length(x),qmax,jmax)))
    } else {
        if(!is.data.frame(groups) || ncol(groups) > 3)
            stop("Argument \"groups\" must be a data frame with 1, 2 or 3 columns.\n")
        if(!all(names(groups)%in%c("group","q","j")))
            stop("Names of \"groups\" argument are not right.\n")
        if(ncol(groups)==1) {
            if(names(groups) != "q")
                stop("The column name of a 1-column \"groups\" should be \"q\".\n")
        }
        if(!("j" %in% names(groups))) {
            if(jmax > 1)
                stop(paste("Column \"j\" of \"groups\" must be specified\n",
                           "when \"jmax\" is greater than 1.\n"))
            groups$j <- 1
        }
        if(!("group" %in% names(groups))) groups$group <- 1:nrow(groups)
   }
      
# If jmax is 1 make sure that all j-entries of "groups" are equal to 1.
    if(jmax==1 & !all(groups$j == 1))
        stop(paste("Not all of the entries of \"groups$j\" are\n",
                   "equal to 1, but \"jmax\" is 1.\n"))

# Make sure all entries of "groups" are in the right range.
    if(any(groups$q != round(groups$q)) | any(groups$j != round(groups$j))
        | any(groups$q <= 0) | any(groups$j <= 0))
            stop(paste("Entries of \"groups$q\" and \"groups$j\"\n",
                       "must be postive integers.\n"))
    if(any(groups$q > qmax) | any(groups$j > pmin(jmax,groups$q)))
        stop(paste("Some entries of the \"q\" or \"j\" columns\n",
                   "of \"groups\" are out of range.\n"))
    ng <- length(unique(groups$group))

    if(add & ng > 1) stop(paste("Cannot add to an existing plot when there\n",
                                "is more than one group of traces.\n"))

# Check on whether "gloss" should be done, and if so, should it be made
# (or is it given).
    do.gloss   <- if(is.logical(gloss)) gloss[1] else TRUE
    make.gloss <- do.gloss && is.logical(gloss)

# Make sure "gloss" and "glind", if given, are lists of vectors
# whose lengths match up appropriately with the groups in "groups".
    if(do.gloss) {
        if(is.null(glind)) { # Add gloss for all traces.
    	    glind <- rep(TRUE,nrow(groups))
        } else {
            if(length(glind) != nrow(groups))
                stop("Mismatch in length of \"glind\" and nrow of \"groups\".\n")
            if(!is.logical(glind))
                stop("Argument \"glind\" must be a logical vector.\n")
        }
    }
    if(make.gloss) {
       gloss <- if(jmax==1) {
          paste("q =",groups$q)
       } else {
          paste("q =",groups$q,"j = ",groups$j)
       }
    } else if(do.gloss) {
        if(length(gloss) != nrow(groups))
            stop("Mismatch in length of \"gloss\" and nrow of \"groups\".\n")
        ok <- sapply(1:ng,function(k,groups,gloss){
                              length(gloss[[k]]) == length(groups[[k]]$q)},
                          groups=groups,gloss=gloss)
    }

# Blank out those entries of "gloss" where the values
# of "glind" are FALSE.
    if(do.gloss) gloss[!glind] <- ""

# Set up multiway array of plots.
    adjMfrow <- !add & (is.null(mfrow) || !is.na(mfrow))
    if(adjMfrow) {
        if(is.null(mfrow)) {
            if(ng==1) mfrow <- c(1,1)
            else if(2 <= ng & ng <= 4) mfrow <- c(2,2)
            else mfrow <- c(3,2)
        }
        np <- prod(mfrow)
        oma <- if(ng > 1 & !is.null(main)) c(0,0,2,0) else rep(0,4)
        opar <- par(mfrow=mfrow, oma=oma)
        on.exit(par(opar))
    } else np <- 1

# Set up xlab, ylab, main, main.panel and ylim, 
    if(!add) {
        if(is.null(xlab)) xlab <- ""
        if(is.null(ylab)) ylab <- ""
        if(is.null(main)) main <- ""
        if(ng > 1) {
            if(is.null(main.panel)) {
                main.panel <- paste("group",1:ng)
            } else if(length(main.panel) == 1) {
                main.panel <- rep(main.panel,ng)
            } else if(length(main.panel) != ng) {
	        stop("Mismatch in lengths of \"main.panel\" and \"groups\".\n")
            }
        }
        if(is.null(ylim)) ylim <- attr(x,'ylim')
        if(is.null(ylim))
            stop(paste("Argument ylim is missing and x has no",
                                   "ylim attribute.\n"))
    }

# Set up xlim.  Needed even, if add is FALSE, since we're using
# plot.function/plot.stepfun.
        if(is.null(xlim)) xlim <- attr(x,'tlim')
        if(is.null(xlim))
            stop(paste("Argument xlim is missing and x has no",
                           "\"tlim\" attribute.\n"))
        xlime <- xlim
        if(do.gloss) {
            if(extend < 0 | extend > 1)
                stop(paste("Crazy value",extend,"for \"extend\".\n"))
            xlime[2] <- xlime[2] + extend*diff(xlime)
        }

# Set up line types and colours.
    if(is.null(lty))  lty  <- 1
    if(is.null(cols)) cols <- 1
    lty  <- rep(lty,length.out=nrow(groups))
    cols <- rep(cols,length.out=nrow(groups))

# A couple of auxiliary constructs ...
    startPlot <- function(xlim,xlime,ylim,xlab,ylab,main) {
        plot(0,0,type="n",xlim=xlime,ylim=ylim,
             xlab=xlab,ylab=ylab,main=main,axes=FALSE)
        axis(side=2)
        axis(side=1,at=pretty(xlim))
    }
    if(do.gloss) x0 <- xlim[2] + 0.05*diff(xlim)
    stride <- inherits(x,"pwc.flap")

# Are you ready boots? Start plotting!
    if(ng==1 & !add) {
            startPlot(xlim=xlim,xlime=xlime,ylim=ylim,
                      xlab=xlab,ylab=ylab,main=main)
    }
    for(kg in 1:ng) {
        if(ng > 1) {
            startPlot(xlim=xlim,xlime=xlime,ylim=ylim,
                      xlab=xlab,ylab=ylab,main=main.panel[kg])
        }
        ikg <- groups$group==kg
        qkg <- groups$q[ikg]
        jkg <- groups$j[ikg]
        lkg <- lty[ikg]
        ckg <- cols[ikg]
        K  <- length(qkg)
        for(k in 1:K) {
            i <- qj2i(qkg[k],jkg[k],qmax)
            if(stride) {
                plot(x[[i]], xlim=xlim,lty=lkg[k], col=ckg[k], add=TRUE,
                     do.points=FALSE, ...)
            } else {
                plot(x[[i]], xlim=xlim,lty=lkg[k], col=ckg[k], add=TRUE,...)
            }
# Marginal gloss.
            if(do.gloss) {
                lbl <- gloss[ikg][k]
                xi  <- x[[i]](xlim[2])
                text(x0,xi,labels=lbl,adj=0,cex=cex.gloss,col=col.gloss)
            }
        }
        if(ng > 1) {
            if(kg %% np == 0 | kg == ng) {
                mtext(outer=TRUE,side=3,line=0,text=main,cex=1.2,font=2)
            }
            if(dev.interactive() & kg < ng & kg%%np == 0) readline('Go? ')
        }
    }
    invisible()
}
