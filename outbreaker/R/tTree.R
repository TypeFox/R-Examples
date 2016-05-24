
#' Simple transmission tree from outreaber's output
#'
#' The S3 class \code{tTree} is used for storing simplified transmission trees,
#' obtained from outbreaker's ouptput (functions \code{outbreaker} and
#' \code{outbreaker.parallel}) using \code{get.tTree}. Some additional features
#' are available for \code{tTree} objects, including plotting (\code{plot}),
#' conversion to \code{igraph} graphs (\code{as.igraph}), and identification of
#' mutations on the branches of the tree (\code{findMutations}).
#'
#'
#' @aliases tTree get.tTree plot.tTree as.igraph.tTree findMutations.tTree
#'
#' @export
#'
#' @rdname tTree
#'
#' @param x for \code{get.tTree}, the output of \code{outbreaker} or
#' \code{outbreaker.parallel}. For other functions, a \code{tTree} object.
#' @param burnin an integer indicating the number of steps of the MCMC to be
#' discarded as burnin period. Defaults to 20,000.
#' @param best a character string matching "ancestries" or "tree", indicating
#' which criterion is used to define the consensus tree; "ancestries" retains,
#' for each case, the most supported ancestor; "tree" retains the most
#' supported tree; note that the latter may exist only in the case of very
#' small epidemics.
#' @param y unused - there for compatibility with the generic of \code{plot}.
#' @param edge.col the color used for the edges; overriden if
#' \code{col.edge.by} is provided.
#' @param col.edge.by a character string indicating how edges should be
#' colored. Can be "dist" (by number of mutations), "n.gen" (by number of
#' generations), or "prob" (by posterior support for the ancestries).
#' @param col.pal the palette of colors to be used for edges; if NULL, a grey
#' palette is used, with larger values in darker shades.
#' @param annot same as \code{col.edge.by}, but specifies the information used
#' to annotated the edges; several values can be provided, in which case
#' different fields will be concatenated to generate the annotation.
#' @param sep a character indicating the separator for different field (see
#' \code{annot}).
#' @param dna a DNAbin object containing the aligned sequences of the isolates
#' in the tree.
#' @param \dots further arguments to be passed to other functions.
#' @return \code{tTree} objects are lists with the following components:
#' \itemize{ \item idx: integer, the index of the cases
#'
#' \item collec.dates: the collection dates of the isolates
#'
#' \item idx.dna: the index of the cases to which each DNA sequence corresponds
#'
#' \item ances: the index of the inferred ancestor, for each case
#'
#' \item inf.dates: the inferred infection date, for each case
#'
#' \item p.ances: the posterior probability of the inferred ancestor (i.e.,
#' proportion in the posterior distribution of ancestors)
#'
#' \item nb.mut: the number of mutations between isolates and their inferred
#' ancestor, for each isolate
#'
#' \item n.gen: the number of generations between isolates and their inferred
#' ancestor, for each isolate
#'
#' \item p.gen: the posterior probability of the inferred number of generations
#' between each case and its inferred ancestor
#'
#' \item inf.curves: the infectivity curves for each case }
#'
#' The plot function invisibly returns the conversion of the \code{tTree}
#' object into a \code{igraph} graph.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @keywords classes
#'
#' @examples
#'
#' data(fakeOutbreak)
#' attach(fakeOutbreak)
#'
#' ## represent posterior ancestries
#' if(require(adegenet)){
#' transGraph(res, annot="", main="Posterior ancestries - support > 0.01",
#'    threshold=0.01, col.pal=spectral)
#' }
#' ## get consensus ancestries
#' tre <- get.tTree(res)
#' plot(tre, annot="", main="Consensus ancestries")
#'
#' ## show match data/consensus ancestries
#' col <- rep("lightgrey", 30)
#' col[which(dat$ances != tre$ances)] <- "pink"
#' plot(tre, annot="", vertex.color=col, main="Consensus ancestries")
#' mtext(side=3, text="cases with erroneous ancestries in pink")
#'
#'
#' detach(fakeOutbreak)
#'
#'
get.tTree <- function(x, burnin=2e4, best=c("ancestries","tree")){

    ## HANDLE ARGUMENTS ##
    if(all(x$chains$step<=burnin)) stop("requested burn-in exeeds the number of chains")
    best <- match.arg(best)


    ## CREATE OUTPUT LIST ##
    res <- list()
    res$idx <- 1:length(x$collec.dates)
    res$collec.dates <- x$collec.dates
    res$idx.dna <- x$idx.dna

    ## PROCESS CHAINS ##
    chains <- x$chains[x$chains$step>burnin, ] # discard burn-in

    ## get ancestors ##
    temp <- chains[,grep("alpha",names(chains))]
    if(best=="ancestries"){
        res$ances <- apply(temp,2, function(e) names(table(e))[which.max(as.numeric(table(e)))])
    }
    if(best=="tree"){
        ## be careful that trees order is not scrambled by table(...)
        allTrees <- apply(temp,1,paste,collapse="-")
        allTrees <- factor(allTrees, levels=unique(allTrees))
        best.tree <- which.max(table(allTrees))
        res$ances <- temp[best.tree,,drop=TRUE]
    }
    res$ances <- as.integer(res$ances)
    res$ances[res$ances<1] <- NA

    ## get infection dates ##
    temp <- chains[,grep("Tinf",names(chains))]
    res$inf.dates <- apply(temp,2,median)

    ## get probabilities of ancestries ##
    temp <- chains[,grep("alpha",names(chains))]
    res$p.ances <- apply(temp,2, function(e) max(as.numeric(table(e)))/sum(table(e)))

    ## get ancestor->descendent mutations ##
    D <- as.matrix(x$D)
    findMut <- function(i){
        if(any(is.na(c(res$idx[i],res$ances[i])))) return(NA)
        if(!all(c(res$idx[i],res$ances[i]) %in% res$idx.dna)) return(NA)
        return(D[res$idx[i],res$ances[i]])
    }
    res$nb.mut <- sapply(1:length(res$idx), function(i) findMut(i))
    ##res$nb.mut <- sapply(1:length(res$idx), function(i) D[res$idx[i],res$ances[i]])

    ## get kappa ##
    temp <- chains[,grep("kappa",names(chains))]
    res$n.gen <- apply(temp,2, function(e) names(table(e))[which.max(as.numeric(table(e)))])
    res$n.gen <- as.integer(res$n.gen)
    res$p.gen <- apply(temp,2, function(e) max(as.numeric(table(e)))/sum(table(e)))

    ## get infectiousness curves ##
    timeSpan <- c(min(res$inf.dates),max(res$inf.dates)+length(x$w))
    f1 <- function(infDate, w){
        dens <- rep(0,diff(timeSpan)+1)
        idxStart <- infDate-timeSpan[1]+1
        dens[idxStart:(idxStart+length(w)-1)] <- w
        dates <- seq(timeSpan[1],timeSpan[2])
        res <- data.frame(dates,dens)
        return(as.matrix(res))
    }

    res$inf.curves <- lapply(res$inf.dates, f1, x$w)


    ## SET CLASS AND RETURN ##
    class(res) <- "tTree"
    return(res)
} # end get.tTree






#' @rdname tTree
#' @export
#'
as.igraph.tTree <- function(x, edge.col="black", col.edge.by="prob",
                              col.pal=NULL, annot=c("dist","n.gen","prob"), sep="/", ...){
    ## if(!require(igraph)) stop("package igraph is required for this operation")
    ## if(!require(adegenet)) stop("adegenet is required")
    if(!inherits(x,"tTree")) stop("x is not a tTree object")
    if(!col.edge.by %in% c("dist","n.gen","prob")) stop("unknown col.edge.by specified")

    ## GET DAG ##
    from.old <- x$ances
    to.old <- x$idx
    isNotNA <- !is.na(from.old) & !is.na(to.old)
    vnames <- sort(unique(c(from.old,to.old)))
    from <- match(from.old,vnames)
    to <- match(to.old,vnames)
    dat <- data.frame(from,to,stringsAsFactors=FALSE)[isNotNA,,drop=FALSE]
    out <- graph.data.frame(dat, directed=TRUE, vertices=data.frame(names=vnames, dates=x$inf.dates[vnames]))

    ## SET VARIOUS INFO ##
    E(out)$dist <- x$nb.mut[isNotNA]
    E(out)$prob <- x$p.ances[isNotNA]
    E(out)$n.gen <- x$n.gen[isNotNA]
    E(out)$p.kappa <- x$p.gen[isNotNA]

   ## SET EDGE COLORS ##
    if(is.null(col.pal)){
        col.pal <- function(n){
            return(grey(seq(0.75,0,length=n)))
        }
    }
    if(col.edge.by=="prob") edge.col <- num2col(E(out)$prob, col.pal=col.pal, x.min=0, x.max=1)
    if(col.edge.by=="dist") edge.col <- num2col(E(out)$dist, col.pal=col.pal, x.min=0, x.max=1)
    if(col.edge.by=="n.gen") edge.col <- num2col(E(out)$n.gen, col.pal=col.pal, x.min=0, x.max=1)

    E(out)$color <- edge.col

    ## SET EDGE LABELS ##
    n.annot <- sum(annot %in% c("dist","n.gen","prob"))
    lab <- ""
    if(!is.null(annot) && n.annot>0){
        if("dist" %in% annot) lab <- E(out)$dist
        if("n.gen" %in% annot) lab <- paste(lab, E(out)$n.gen, sep=sep)
        if("prob" %in% annot) lab <- paste(lab, round(E(out)$prob,2), sep=sep)
    }
    lab <- sub(paste("^",sep,sep=""),"",lab)
    E(out)$label <- lab


    ## SET LAYOUT ##
    attr(out, "layout") <- layout.fruchterman.reingold(out, params=list(minx=x$inf.dates, maxx=x$inf.dates))

    return(out)
} # end as.igraph.tTree







## ##################
## ## [.tTree
## ##################
## "[.tTree" <- function(x,i,j,drop=FALSE){
##     res <- x
##     res$idx <- res$idx[i]
##     res$ances <- res$ances[i]
##     res$p.ances <- res$p.ances[i]
##     res$n.gen <- res$n.gen[i]
##     res$p.gen <- res$p.gen[i]
##     res$inf.curves <- res$inf.curves[i]
##     res$collec.dates <- res$collec.dates[i]
##     res$inf.dates <- res$inf.dates[i]
##     res$n <- length(res$ances)
##     res$nb.mut <- res$nb.mut[i]

##     return(res)
## }






#' @rdname tTree
#' @export
#'
findMutations.tTree <- function(x, dna, ...){
    ## CHECKS ##
    ## if(!require(ape)) stop("the ape package is needed")
    if(!inherits(x,"tTree")) stop("x is not a tTree object")

    ## ## function to pull out mutations from sequence a to b ##
    ## f1 <- function(a,b){
    ##     seqa <- as.character(dna[a,])
    ##     seqb <- as.character(dna[b,])
    ##     temp <- which(seqa != seqb)
    ##     ori <- seqa[temp]
    ##     mut <- seqb[temp]
    ##     names(ori) <- names(mut) <- temp
    ##     toRemove <- !ori %in% c('a','t','g','c') | !mut %in% c('a','t','g','c')
    ##     ori <- ori[!toRemove]
    ##     mut <- mut[!toRemove]
    ##     res <- data.frame(ori,mut)
    ##     names(res) <- rownames(dna)[c(a,b)]
    ##     res$short <- paste(row.names(res),":",res[,1],"->",res[,2],sep="")
    ##     return(res)
    ## }

    ## ## get mutations for each ancestry
    ## isNotNA <- which(!(is.na(x$ances) | is.na(x$idx)))
    ## out <- lapply(isNotNA, function(i) f1(x$ances[i], x$idx[i]))

    ## GET PAIRS TO COMPARE ##
    pairs <- cbind(x$ances,x$idx)
    isNotNA <- which(!(is.na(x$ances) | is.na(x$idx)))
    if(length(isNotNA)==0) return()
    pairs <- pairs[isNotNA,,drop=FALSE]

    ## CALL DNABIN METHOD ##
    out <- findMutations(dna, pairs)

    return(out)

} # end findMutations








#' @rdname tTree
#' @export
#'
plot.tTree <- function(x, y=NULL, edge.col="black", col.edge.by="prob",
                              col.pal=NULL, annot=c("dist","n.gen","prob"), sep="/", ...){
    ## if(!require(igraph)) stop("igraph is required")
    ## if(!require(adegenet)) stop("adegenet is required")
    if(!inherits(x,"tTree")) stop("x is not a tTree object")
    if(!col.edge.by %in% c("dist","n.gen","prob")) stop("unknown col.edge.by specified")

    ## get graph ##
    g <- as.igraph(x, edge.col=edge.col, col.edge.by=col.edge.by, col.pal=col.pal, annot=annot, sep=sep)

     ## make plot ##
    plot(g, layout=attr(g,"layout"), ...)

    ## return graph invisibly ##
    return(invisible(g))

} # end plot.tTree





## ###################
## ## refine.mutrate
## ###################
## refine.mutrates <- function(x, min.support=0.5){
##     ## CHECKS ##
##     if(!inherits(x,"tTree")) stop("x is not a tTree object")

##     ## AUXILIARY FUNCTIONS ##

##     ## REFINE THE MUTATION RATE ##
##     ## get only retained ancestries
##     kept <- which((x$p.ances >= min.support)[x$idx.dna])
##     if(length(kept)<1) {
##         warning("No information retained - min.support may be too high")
##         return(NULL)
##     }


##     ## get mutation rate - per genome / per generation
##     mu <- mean(x$nb.mut[kept]/x$n.gen[kept],na.rm=TRUE)


## } # end refine.mutrates
