

#' Plot outbreaker's results
#'
#' These are the main functions used for generating graphics from the raw
#' output of \code{outbreaker} and \code{outbreaker.parallel}.
#'
#' \itemize{ \item \code{plotChains} is used for plotting MCMCs
#'
#' \item \code{transGraph} plots a graph of inferred ancestries
#'
#' \item \code{plotOutbreak} attempts to synthetize the reconstruction of small
#' outbreaks }
#'
#'
#' @aliases plotChains transGraph plotOutbreak
#'
#' @rdname plot
#'
#' @export
#'
#' @param x the output of \code{outbreaker} or \code{outbreaker.parallel}.
#' @param what a character chains giving the name of the item to be plotted.
#' See \code{names(x$chains)} for possible values. By default, log-posterior
#' values are plotted
#' @param type a character indicating if the chains should be plotted as time
#' series ("series"), or as density ("density").
#' @param burnin an integer indicating the number of MCMC steps to discard
#' before plotting chains.
#' @param dens.all a logical indicating if, in the case of multiple runs, the
#' overall density of the different chains should be plotted in addition to
#' individual densities.
#' @param col a vector of colors to be used to plot different chains.
#' @param lty a vector of integers specifying line types for the different
#' chains.
#' @param lwd same as \code{lty}, but for line width.
#' @param main the title to be added to the plot.
#' @param labels the labels to be used to name the nodes of the graph (cases).
#' @param threshold the minimum support for ancestries to be plotted; 'support'
#' is defined as the frequency of a given ancestor in the posterior
#' distribution; defaults to 0.2.
#' @param thres.hide a threshold of posterior support for displaying
#' ancestries; ancestries with less than this frequency in the posterior are
#' hidden.
#' @param col.pal,edge.col.pal the color palette to be used for the edges
#' (ancestries).
#' @param curved.edges a logical indicating whether edges should be curved.
#' @param col.edge.by a character string indicating which information should be
#' used to color the edges ('dist': genetic distance; 'prob': support for the
#' ancestry)
#' @param annot a character indicating which information should be used to
#' annotate the edges; this can be the distances between ancestors and
#' descendents ("dist") and the posterior support for ancestries ("support");
#' if both are requested, fields will be concatenated.
#' @param sep a character indicating the separator to be used when
#' concatenating several types of annotation.
#' @param cex.bubble a numeric value indicating the size factor for the bubbles
#' representing the generation time distribution.
#' @param edge.max.dist a number indicating the threshold distance bounding the
#' color palette used for the edges; useful to avoid showing edges
#' corresponding to distances larger than a given number.
#' @param lwd.arrow a numeric value indicating the size factor for the arrows.
#' @param xlim the limits of the X axis; if NULL, determined from the data.
#' @param legend a logical indicating if a legend should be plotted for the
#' different runs.
#' @param posi a character string indicating the position of the legend (see
#' \code{?legend}).
#' @param \dots further arguments to be passed to other functions.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @examples
#'
#' data(fakeOutbreak)
#' attach(fakeOutbreak)
#'
#' ## examine MCMC
#' plotChains(res)
#' plotChains(res,type="dens")
#' plotChains(res,type="dens", what="mu1", burnin=2e4)
#'
#' ## represent posterior ancestries
#' transGraph(res, annot="", main="Posterior ancestries")
#' transGraph(res, annot="", main="Posterior ancestries - support > 0.5",
#'    threshold=0.5)
#' if(require(adegenet)){
#' transGraph(res, annot="", main="Posterior ancestries - support > 0.01",
#'    threshold=0.01, col.pal=spectral)
#' }
#' ## summary plot
#' plotOutbreak(res,cex.bubble=0.5, thres.hide=0.5,
#'    main="Outbreak reconstruction")
#'
#'
#' detach(fakeOutbreak)
#'
#'
plotChains <- function(x, what="post", type=c("series","density"), burnin=0, dens.all=TRUE,
                       col=funky(x$n.runs), lty=1, lwd=1, main=what,
                       legend=TRUE, posi="bottomleft", ...){
    ## HANDLE ARGUMENTS ##
    type <- match.arg(type)
    n.runs <- x$n.runs
    col.ori <- col
    if(!what %in% names(x$chains)) stop(paste(what,"is not a column of x$chains"))
    if(!is.null(col)) col <- rep(col, length = n.runs)
    if(!is.null(lty)) lty <- rep(lty, length = n.runs)
    if(!is.null(lwd)) lwd <- rep(lwd, length = n.runs)
    if(is.null(burnin)){
        burnin <- max(x$burnin, x$find.import.at, x$tune.end)
    }

    ## GET DATA TO PLOT ##
    dat <- cbind(x$chains$step[x$chains$run==1],data.frame(split(x$chains[,what], x$chains$run)))
    names(dat) <- c("step", paste(what, 1:n.runs,sep=""))
    if(!any(dat$step>burnin)) stop("burnin is greater than the number of steps in x")
    dat <- dat[dat$step>burnin,,drop=FALSE]

    ## MAKE PLOT ##
    if(type=="series"){
        matplot(dat$step, dat[,-1,drop=FALSE], type="l", col=col, lty=lty, xlab="MCMC iteration", ylab="value", main=main, ...)
    }

    if(type=="density"){
        ## add general density if needed ##
        temp <- lapply(dat[, -1, drop=FALSE], density)
        if(dens.all){
            temp[[n.runs+1]] <- density(unlist(dat[,-1,drop=FALSE]))
            col <- c(col, "black")
            lty <- c(lty, 1)
            lwd <- c(lwd, 3)
            n.runs <- n.runs+1
        }
        range.x <- range(sapply(temp, function(e) e$x))
        range.y <- range(sapply(temp, function(e) e$y))
        plot(1,type="n", xlim=range.x, ylim=range.y, xlab="value", ylab="density", main=main, ...)
        invisible(lapply(1:n.runs, function(i) lines(temp[[i]], col=col[i], lty=lty[i], lwd=lwd[i])))
    }

    ## ADD LEGEND ##
    if(legend){
        legend(posi, fill=col, title="Runs", leg=1:length(col.ori))
    }
    return(invisible())
} # end plotChains






#' @rdname plot
#' @export
transGraph <- function(x, labels=NULL, burnin=x$burnin, threshold=0.2, col.pal=NULL, curved.edges=TRUE,
                       annot=c("dist","support"), sep="/", ...){
    ## CHECKS ##
    ## if(!require(igraph)) stop("igraph is required")
    ## if(!require(adegenet)) stop("adegenet is required")

    ## HANDLE ARGUMENTS ##
    if(burnin> max(x$chains$step)) stop("burnin exceeds the number of chains in the output")
    if(is.null(col.pal)){
        col.pal <- function(n){
            return(grey(seq(1,0,length=n)))
        }
    }

    ## GET ANCESTRY DATA ##
    ances <- x$chains[x$chains$step>=burnin, grep("alpha",names(x$chains)),drop=FALSE]
    tabances <- apply(ances,2,table)
    N <- ncol(ances)


    ## GET DATA FOR GRAPH ##
    ## ancestor, descendents, list of nodes
    to.old <- rep(1:N, sapply(tabances,length))
    from.old <- as.numeric(unlist(lapply(tabances,names)))
    from.old[from.old<1] <- NA
    to.old[to.old<1] <- NA
    isNotNA <- !is.na(from.old) & !is.na(to.old)
    vnames <- sort(unique(c(from.old,to.old)))
    from <- match(from.old,vnames)
    to <- match(to.old,vnames)

    ## support for the ancestries
    support <- unlist(lapply(tabances, function(e) e/sum(e)))[isNotNA]
    edge.col <- num2col(support, col.pal=col.pal, x.min=0, x.max=1)

    ## average dates of infection
    Tinf <- x$chains[x$chains$step>=burnin, grep("Tinf",names(x$chains)),drop=FALSE]
    inf.dates <- apply(Tinf,2,mean)
    names(inf.dates) <- 1:N

    ## remove weakly supported ancestries
    dat <- data.frame(from,to,stringsAsFactors=FALSE)[isNotNA,,drop=FALSE]
    if(sum(support>=threshold)==0) warning("threshold to high - no edge left")
    dat <- dat[support>=threshold, , drop=FALSE]

    ## convert to graph
    out <- graph.data.frame(dat, directed=TRUE, vertices=data.frame(names=vnames, dates=inf.dates[as.character(vnames)]))

    ## get ancestor->descendent mutations ##
    D <- as.matrix(x$D)
    findMut <- function(i){
        if(any(is.na(c(to[i],from[i])))) return(NA)
        return(D[to[i],from[i]])
    }
    nb.mut <- sapply(1:length(to), function(i) findMut(i))
    nb.mut <- nb.mut[isNotNA]


    ## SET PARAMETERS OF THE GRAPH ##
    ## vertices
    if(is.null(labels)){
        V(out)$label <- vnames
    } else {
        V(out)$label <- labels
    }


    ## edges
    E(out)$color <- edge.col[support>=threshold]
    E(out)$support <- support[support>=threshold]
    E(out)$curved <- curved.edges
    E(out)$nb.mut <- nb.mut[support>=threshold]

    ## edge labels
    lab <- ""
    if(!is.null(annot) && length(annot)>0){
        if(any(c("dist","nb.mut","mut") %in% annot)) lab <- E(out)$nb.mut
        if(any(c("support","prob") %in% annot)) lab <- paste(lab, round(E(out)$support,2), sep=sep)
    }
    lab <- sub(paste("^",sep,sep=""),"",lab)
    E(out)$label <- lab


    ## set layout
    attr(out, "layout") <- layout.fruchterman.reingold(out, params=list(minx=V(out)$onset, maxx=V(out)$onset))


    ## MAKE THE PLOT ##
    plot(out, layout=attr(out, "layout"), ...)


    ## RETURN OBJECT ##
    return(invisible(out))

} # end transGraph






################
## plotOutbreak
################
.entropy <- function(p){
    p <- p/sum(p, na.rm=TRUE)
    return(-sum(p*log(p), na.rm=TRUE))
}


#' @rdname plot
#' @export
plotOutbreak <- function(x, burnin=x$burnin, thres.hide=0.2, col=NULL,
                         col.pal=colorRampPalette(c("blue","lightgrey")), edge.col.pal=NULL,
                         col.edge.by="prob", annot=c("dist","prob"), sep="/", cex.bubble=1,
                         edge.max.dist=10,lwd.arrow=2, xlim=NULL, ...){
    ## CHECKS ##
    ## if(!require(adegenet)) stop("adegenet is not installed")
    col.edge.by <- match.arg(col.edge.by, c("dist","prob"))

    ## GET TREE ##
    tre <- get.tTree(x,burnin=burnin)
    N <- length(tre$idx)

    ## GET NUMBER OF MUTATIONS BETWEEN SEQUENCES
    M <- as.matrix(x$D)

    ## GET ALPHA_I ##
    alpha <- x$chains[x$chains$step>burnin, grep("alpha",names(x$chains))]
    f1 <- function(vec){
        sapply(1:N, function(i) mean(vec==i,na.rm=TRUE))
    }

    ## support for ancestries
    alphadat <- apply(alpha,2,f1)

    ## define colors for the individuals
    ## default: based on entropy of ancestries support
    if(is.null(col)){
        entropy <-apply(alphadat, 2, .entropy)
        col <- num2col(entropy, col.pal=col.pal)
    } else {
        entropy <- NULL
    }


    ## PLOT INFECTION DATES ##
    toKeep <- grep("Tinf",names(x$chains))
    Tinf <- x$chains[x$chains$step>burnin, toKeep]
    colnames(Tinf) <- colnames(as.matrix(x$D))

    ## find max date
    if(is.null(xlim)){
        min.date <- min(x$chains[,grep("Tinf",names(x$chains))])
        max.date <- max(tre$inf.curves[[1]][,1])
        xlim <- c(min.date-.5,  max.date+.5)
    }

    ## basic boxplot
    boxplot(Tinf, col=col, horizontal=TRUE, las=1, ylim=xlim, ...)

    ## add infectious periods
    lapply(1:N, function(i) points(tre$inf.curves[[i]][,1], rep(i, nrow(tre$inf.curves[[i]])),
                                   cex=sqrt(tre$inf.curves[[i]][,2])*10*cex.bubble,
                                   col=transp(col)[i], pch=19) )

    ## plot collection dates
    points(x$collec.dates, 1:N, col="black", pch=20)


    ## ADD ANCESTRIES
    if(is.null(edge.col.pal)){
        edge.col.pal <- function(n){
            return(grey(seq(1,0,length=n)))
        }
    }

    ## ## FUNCTION TO FIND THE NUMBER OF MUTATIONS DISPLAYED ##
    ## ## (required to define the color scale for arrows ##
    ## get.nb.mut <- function(from, to){
    ##      support <- alphadat[from,to]
    ##      nb.mut <- M[from,to]
    ##      if(support<thres.hide) nb.mut <- 0
    ##      return(nb.mut)
    ##  }

    ## FUNCTION TO DRAW ARROWS ##
    drawArrow <- function(from, to){
        if(is.na(from)||from<1) return(invisible())
        ## get stuff for arrows ##
        infdays <- apply(Tinf,2,mean)
        ##x.to <- x.from <- infdays[to] # for arrows on day of infection
        x.from <- infdays[from]
        x.to <- infdays[to]
        y.from <- from
        y.to <- to
        support <- alphadat[from,to]
        ## col <- col[from]
        if(col.edge.by=="prob"){
            edge.col <- num2col(support, x.min=0, x.max=1, col.pal=edge.col.pal)
        }
        if(col.edge.by=="dist"){
            edge.col <- num2col(as.integer(M[from,to]), x.min=0, x.max=edge.max.dist, col.pal=edge.col.pal)
        }
        edge.col[support<thres.hide] <- "transparent"
        ## lwd <- round(support*arrow.max)
        ## col.back <- rep("transparent",N)
        ## col.back[lwd>=2] <- "black"

        ## draw arrows ##
        arrows(x.from, y.from, x.to, y.to, col=edge.col, length=0.1, angle=20, lwd=lwd.arrow)

        ## get stuff for annotations ##
        if(support>=thres.hide){
            x.ann <- (x.from + x.to)/2
            y.ann <- 0.15+(y.from + y.to)/2

            lab <- ""
            if(!is.null(annot) && length(annot)>0){
                if(any(c("dist","nb.mut","mut") %in% annot)) lab <- M[from,to]
                if(any(c("support","prob") %in% annot)) lab <- paste(lab, round(support,2), sep=sep)
            }
            lab <- sub(paste("^",sep,sep=""),"",lab)
            text(x.ann,y.ann,lab)
        }
        return(invisible())
    }

    ## DRAW ALL ARROWS ##
    ances <- apply(alpha,2,function(e) table(e)/sum(e))
    names(ances) <- 1:N

    ##  ## need to determine max nb of mutations if col.edge.by=="dist" ##
    ## if(col.edge.by=="dist"){
    ##     allMut <- unlist(lapply(1:N, function(i) sapply(as.integer(names(ances[[i]])), function(from) get.nb.mut(from, i))))
    ##     max.nb.mut <- max(allMut,na.rm=TRUE)
    ## }

    lapply(1:N, function(i) sapply(as.integer(names(ances[[i]])), function(from) drawArrow(from, i)))


    ## BUILT RESULT AND RETURN ##
    res <- list(col=col, col.pal=col.pal, entropy=entropy, edge.col.pal=edge.col.pal)
    return(invisible(res))
} # end plotOutbreak











## #############
## ## epicurves
## #############
## epicurves <- function (x, col=NULL, bg="lightgrey", line.col="white", coef=1, max.arr=5,...) {
##     if(!require(adegenet)) stop("adegenet is required")
##     if(!inherits(x,"tTree")) stop("x is not a tTree object")

##     ## GET USEFUL INFO ##
##     N <- length(x$idx)
##     timeSpan <- range(x$inf.curves[[1]][,1])
##     if(is.null(col)){
##         colPal <- colorRampPalette(c("grey30","blue","green3","orange","red","brown4","purple"))
##         col <- colPal(N)
##     }


##     ## MAKE EMPTY PLOT ##
##     plot(0,0,type="n",xlim=timeSpan+c(-1,1), ylim=c(0,N+1), xlab="Dates",ylab="Individual index",...)
##     rect(timeSpan[1]-2,-2,timeSpan[2]+2,N+2, col=bg)
##     abline(h=1:N, col="white",lwd=3)
##     abline(h=1:N, col=transp(col),lwd=2)
##     abline(v=pretty(timeSpan,4), col=line.col)

##     ## DRAW INFECTIOUSNESS CURVES ##
##     for(i in 1:N){
##         temp <- x$inf.curves[[i]][x$inf.curves[[i]][,2]> 1e-12,,drop=FALSE]
##         x.coord <- c(temp[,1], rev(temp[,1]))
##         y.coord <- c(i+temp[,2]*coef, rep(i,nrow(temp)))
##         polygon(x.coord, y.coord, col=transp(col[i]),border=col[i])
##         points(temp[,1], i+(temp[,2]*coef), type="o", pch=20,cex=0.5, col=col[i])
##     }


##     ## ADD COLLECTION DATES ##
##     points(x$collec.dates, 1:N, pch=20, cex=2, col=col)


##     ## ADD INFECTIONS DATES ##
##     points(x$inf.dates, 1:N, cex=2, col=col)


##     ## ADD INFECTIONS ##
##     arr.w <- (x$p.ances- 1/(N-1)) * max.arr
##     arr.w[arr.w<0.5] <- 0.5
##     arrows(x$inf.date,x$ances, x$inf.dates, x$idx, angle=15, col=col[x$ances], lwd=arr.w)

## } # end epicurves
