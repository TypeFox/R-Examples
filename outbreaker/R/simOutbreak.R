############
## disperse
############
## random movement with a square area
## with 'bounce' effect
## xy are the current xy coordinates
## area.size is the max xy coord (length of an edge of the square)
disperse <- function(xy, disp=.1, area.size=10){
    out <- xy + matrix(rnorm(2*nrow(xy), mean=0, sd=disp),ncol=2)

    ## bounce back if too far right/up
    out[out>area.size] <- area.size - (out[out>area.size] - area.size)

    ## bounce back if too far left/down
    out[out<0] <- -out[out<0]

    ## bouncing could theoretically still go wrong
    ## if so, use 'hard edge'
    out[out<0] <- 0
    out[out>area.size] <- area.size

    ## return
    return(out)
} # end disperse


## fun test:
## library(adegenet)
## xy <- matrix(runif(40, min=0, max=10), ncol=2)
## for(i in 1:2000) plot(xy <- disperse(xy), col=transp(funky(20),.8), cex=10, pch=20, main=i, xlim=c(0,10), ylim=c(0,10))



################
## .kernel.expo
################
## compute an exponential kernel for a set of points,
## computing pairwise distances first
## k(L_i,L_{\alpha_i}) = \beta e^{-\theta * h(L_i,L_{\alpha_i})}
.kernel.expo <- function(xy, mean){

    ## compute output
    D <- as.matrix(dist(xy))
    out <- dexp(D, rate=1/mean)
    diag(out) <- 0

    return(out)
} # end .kernel.expo




#####################
## .plot.kernel.expo
#####################
.plot.kernel.expo <- function(mean=1, ...){
    plot(function(x)  return(dexp(x, rate=1/mean)), xlim=c(0,10), ylab="Density")
    return(invisible(NULL))
}




#' Simulation of pathogen genotypes during disease outbreaks
#'
#' The function \code{simOutbreak} implements simulations of disease outbreaks.
#' The infectivity of cases is defined by a generation time distribution. The
#' function \code{as.igraph} allows to convert simulated transmission trees
#' into \code{igraph} objects.
#'
#' @export
#'
#' @rdname simOutbreak
#'
#' @aliases simOutbreak print.simOutbreak [.simOutbreak labels.simOutbreak
#' simOutbreak-class as.igraph.simOutbreak plot.simOutbreak disperse
#'
#' @param R0 the basic reproduction number; to use several groups, provide a
#' vector with several values.
#' @param infec.curve a \code{numeric} vector describing the individual
#' infectiousness at time t=0, 1, \dots{}
#' @param n.hosts the number of susceptible hosts at the begining of the
#' outbreak
#' @param duration the number of time steps for which simulation is run
#' @param seq.length an integer indicating the length of the simulated
#' haplotypes, in number of nucleotides.
#' @param mu.transi the rate of transitions, in number of mutation per site and
#' per time unit.
#' @param mu.transv the rate of transversions, in number of mutation per site
#' and per time unit.
#' @param rate.import.case the rate at which cases are imported at each time
#' step.
#' @param diverg.import the number of time steps to the MRCA of all imported
#' cases.
#' @param spatial a logical indicating if a spatial model should be used.
#' @param disp the magnitude of dispersal (standard deviation of a normal
#' distribution).
#' @param area.size the size of the square area to be used for spatial
#' simulations.
#' @param reach the mean of the exponential kernel used to determine new
#' infections.
#' @param plot a logical indicating whether an animated plot of the outbreak
#' should be displayed; only available with the spatial model.
#' @param group.freq the frequency of the different groups; to use several
#' groups, provide a vector with several values.
#' @param x,object \code{simOutbreak} objects.
#' @param i,j,drop \code{i} is a vector used for subsetting the object. For
#' instance, \code{i=1:3} will retain only the first three haplotypes of the
#' outbreak. \code{j} and \code{drop} are only provided for compatibility, but
#' not used.
#' @param y present for compatibility with the generic 'plot' method. Currently
#' not used.
#' @param col the color of the vertices of the plotted graph.
#' @param edge.col the color of the edges of the plotted graph; overridden by
#' \code{col.edge.by}.
#' @param col.edge.by a character indicating the type of information to be used
#' to color the edges; currently, the only valid value is "dist" (distances, in
#' number of mutations). Other values are ignored.
#' @param vertex.col the colors to be used for the vertices (i.e., cases).
#' @param edge.col.pal the color palette to be used for the edges; if NULL, a
#' grey scale is used, with darker shades representing larger values.
#' @param annot a character indicating the information to be used to annotate
#' the edges; currently accepted values are "dist" (genetic distances, in
#' number of mutations), and "n.gen" (number of generations between cases).
#' @param sep a character used to separate fields used to annotate the edges,
#' whenever more than one type of information is used for annotation.
#' @param \dots further arguments to be passed to other methods
#' @return === simOutbreak class ===\cr \code{simOutbreak} objects are lists
#' containing the following slots:\cr \itemize{ \item n: the number of cases in
#' the outbreak\cr \item dna: DNA sequences in the DNAbin matrix format\cr
#' \item dates: infection dates\cr \item dynam: a data.frame containing, for
#' each time step (row), the number of susceptible, infected, or recovered in
#' the population. \cr \item id: a vector of integers identifying the cases\cr
#' \item ances: a vector of integers identifying infectors ('ancestor')\cr
#' \item nmut: the number of mutations corresponding to each ancestry\cr \item
#' ngen: the number of generations corresponding to each ancestry\cr \item
#' call: the matched call }
#' @author Implementation by Thibaut Jombart \email{t.jombart@@imperial.ac.uk}.
#'
#' Epidemiological model designed by Anne Cori and Thibaut Jombart.
#' @examples
#'
#' \dontrun{
#' dat <- list(n=0)
#'
#' ## simulate data with at least 30 cases
#' while(dat$n < 30){
#'    dat <- simOutbreak(R0 = 2, infec.curve = c(0, 1, 1, 1), n.hosts = 100)
#' }
#' dat
#'
#' ## plot first 30 cases
#' N <- dat$n
#' plot(dat[1:(min(N,30))], main="First 30 cases")
#' mtext(side=3, text="nb mutations / nb generations")
#'
#' ## plot a random subset (n=10) of the first cases
#' x <- dat[sample(1:min(N,30), 10, replace=FALSE)]
#' plot(x, main="Random sample of 10 of the first 30 cases")
#' mtext(side=3, text="nb mutations / nb generations")
#'
#' ## plot population dynamics
#' head(dat$dynam,15)
#' matplot(dat$dynam[1:max(dat$onset),],xlab="time",
#'    ylab="nb of individuals", pch=c("S","I","R"), type="b")
#'
#'
#' ## spatial model
#' w <-  exp(-sqrt((1:40)))
#' x <- simOutbreak(2, w, spatial=TRUE,
#'                  duration=500, disp=0.1, reach=.2)
#'
#' ## spatial model, no dispersal
#' x <- simOutbreak(.5, w, spatial=TRUE,
#'                  duration=500, disp=0, reach=5)
#' }
#'
simOutbreak <- function(R0, infec.curve, n.hosts=200, duration=50,
                        seq.length=1e4, mu.transi=1e-4, mu.transv=mu.transi/2,
                        rate.import.case=0.01, diverg.import=10, group.freq=1,
                        spatial=FALSE, disp=0.1, area.size=10, reach=1,
                        plot=spatial){

    ## HANDLE ARGUMENTS ##
    ## handle group sizes
    if(any(group.freq<0)) stop("negative group frequencies provided")
    group.freq <- group.freq/sum(group.freq)
    K <- length(group.freq)
    ## host.group <- sample(1:K, size=n.hosts, prob=group.freq, replace=TRUE)
    R0 <- rep(R0, length=K) # recycle R0

    ## normalize gen.time
    infec.curve <- infec.curve/sum(infec.curve)
    infec.curve <- c(infec.curve, rep(0, duration)) # make sure dates go all the way
    t.clear <- which(diff(infec.curve<1e-10)==1) # time at which infection is cleared

    ## GENETIC FUNCTIONS ##
    NUCL <- as.DNAbin(c("a","t","c","g"))
    TRANSISET <- list('a'=as.DNAbin('g'), 'g'=as.DNAbin('a'), 'c'=as.DNAbin('t'), 't'=as.DNAbin('c'))
    TRANSVSET <- list('a'=as.DNAbin(c('c','t')), 'g'=as.DNAbin(c('c','t')), 'c'=as.DNAbin(c('a','g')), 't'=as.DNAbin(c('a','g')))


    ## AUXILIARY FUNCTIONS ##
    ## generate sequence from scratch
    seq.gen <- function(){
        ##res <- list(sample(NUCL, size=seq.length, replace=TRUE)) # DNAbin are no longer lists by default
        res <- sample(NUCL, size=seq.length, replace=TRUE)
        class(res) <- "DNAbin"
        return(res)
    }

    ## create substitutions for defined SNPs - no longer used
    substi <- function(snp){
        res <- sapply(1:length(snp), function(i) sample(setdiff(NUCL,snp[i]),1)) # ! sapply does not work on DNAbin vectors directly
        class(res) <- "DNAbin"
        return(res)
    }

    ## create transitions for defined SNPs
    transi <- function(snp){
        res <- unlist(TRANSISET[as.character(snp)])
        class(res) <- "DNAbin"
        return(res)
    }

    ## create transversions for defined SNPs
    transv <- function(snp){
        res <- sapply(TRANSVSET[as.character(snp)],sample,1)
        class(res) <- "DNAbin"
        return(res)
    }

    ## duplicate a sequence (including possible mutations)
    seq.dupli <- function(seq, T){ # T is the number of time units between ancestor and decendent
        ## transitions ##
        n.transi <- rbinom(n=1, size=seq.length*T, prob=mu.transi) # total number of transitions
        if(n.transi>0) {
            idx <- sample(1:seq.length, size=n.transi, replace=FALSE)
            seq[idx] <- transi(seq[idx])
        }

        ## transversions ##
        n.transv <- rbinom(n=1, size=seq.length*T, prob=mu.transv) # total number of transitions
        if(n.transv>0) {
            idx <- sample(1:seq.length, size=n.transv, replace=FALSE)
            seq[idx] <- transv(seq[idx])
        }
        return(seq)
    }

    ## define the group of 'n' hosts
    choose.group <- function(n){
        out <- sample(1:K, size=n, prob=group.freq, replace=TRUE)
        return(out)
    }

    ## handle 'plot' argument ##
    if(plot && !spatial) warning("Plot only available with spatial model")


    ## MAIN FUNCTION ##
    ## initialize results ##
    dynam <- data.frame(nsus=integer(duration+1), ninf=integer(duration+1), nrec=integer(duration+1))
    rownames(dynam) <- 0:duration
    res <- list(n=1, dna=NULL, onset=NULL, id=NULL, ances=NULL, dynam=dynam)
    res$dynam$nsus[1] <- n.hosts-1
    res$dynam$ninf[1] <- 1
    res$onset[1] <- 0
    res$id <- 1 # id of infected individuals
    res$ances <- NA
    res$group <- choose.group(1)
    EVE <- seq.gen()
    res$dna <- matrix(seq.dupli(EVE, diverg.import),nrow=1)
    class(res$dna) <- "DNAbin"
    if(spatial) {
        ## current coordinates
        res$xy <- matrix(runif(n.hosts*2, min=0, max=area.size), ncol=2)

        ## location when infected
        res$inf.xy <- res$xy[1,,drop=FALSE]
    }
    res$status <- c("I", rep("S", n.hosts-1)) # will be I, S, or R


    ## run outbreak ##
    for(t in 1:duration){
        ## DETERMINE NEW INTERNAL INFECTIONS ##
        ## individual force of infection - purely based on symptom onset
        indivForce <- infec.curve[t-res$onset+1]

        ## individual force of infection - spatial case
        if(spatial){
            ## disperse
            res$xy <- disperse(res$xy, disp=disp, area.size=area.size)
            if(plot) {
                myCol <- rep("black", nrow(res$xy))
                myCol[res$status=="I"] <- "red"
                myCol[res$status=="R"] <- "royalblue"
                plot(res$xy, pch=20, cex=6, col=transp(myCol),
                     main=paste("time:", t), xlab="", ylab="")
                ##      main=paste("time:", t), xlab="", ylab="", fg=transp(myCol))
                ## symbols(res$xy, circ=rep(reach,nrow(res$xy)), inches=FALSE, bg=transp(myCol),
                ##      main=paste("time:", t), xlab="", ylab="", fg=transp(myCol))
            }

            ## compute kernels
            k.spa <- .kernel.expo(res$xy, mean=reach)


            ## here, using (force f_i):
            ## f_i = w(t-t_i) x k(i,1) x R0
            ##      + w(t-t_i) x k(i,2) x R0
            ##      + ...
            ##      + w(t-t_i) x k(i,N) x R0
            ##     = w(t-t_i) x \sum_j k(i,j) x R0

            ## keep only contacts with susceptibles
            spa.force <- apply(k.spa[, res$status=="S", drop=FALSE], 1, sum)

            ## keep only values for infected indiv
            spa.force <- spa.force[res$id]
            if(length(indivForce)!=length(spa.force)) warning("temporal and spatial forces of infection have different length") # sanity check
            indivForce <- indivForce *  spa.force

            ## prevent over-shooting R0
            indivForce[indivForce>1] <- 1

            ## note: without scaling, we may risk overshooting R0
            ## as apply(K,1,sum) can be > 1
        }

        ## temporal (spatial) force of infection * R0
        indivForce <- indivForce * R0[res$group]

        ## global force of infection (R0 \sum_j I_t^j / N)
        N <- res$dynam$nrec[t] + res$dynam$ninf[t] + res$dynam$nsus[t] # this may change because of imports
        globForce <- sum(indivForce)/N

        ## stop if no ongoing infection in the population
        if(globForce < 1e-12) break;

        ## compute proba of infection for each susceptible
        p <- 1-exp(-globForce)

        ## number of new infections
        nbNewInf <- rbinom(1, size=res$dynam$nsus[t], prob=p)


        ## HANDLE NEW INTERNAL INFECTIONS ##
        if(nbNewInf>0){
            ## dates of new infections ##
            res$onset <- c(res$onset, rep(t,nbNewInf))

            ## identify the infectors of the new cases ##
            newAnces <- sample(res$id, size=nbNewInf, replace=TRUE, prob=indivForce)
            res$ances <- c(res$ances,newAnces)

            ## find the groups of the new cases ##
            newGroup <- choose.group(nbNewInf)
            res$group <- c(res$group,newGroup)

            ## id of the new cases ##
            if(!spatial){ # non-spatial case - ID doesn't matter
                areSus <- which(res$status=="S") # IDs of susceptibles
                newId <- sample(areSus, size=nbNewInf, replace=FALSE)
                res$id <- c(res$id, newId)
                res$status[newId] <- "I"
            } else {
                for(i in 1:nbNewInf){ # for each new infection
                    areSus <- which(res$status=="S") # IDs of susceptibles
                    newId <- sample(areSus, 1, prob=k.spa[newAnces[i],areSus]) # prob depend on location
                    res$id <- c(res$id, newId)
                    res$status[newId] <- "I"
                    res$inf.xy <- rbind(res$inf.xy, res$xy[newId]) # set coords at infection
                }
            }

            ## dna sequences of the new cases ##
            ## molecular clock / generation
            ## newSeq <- t(sapply(match(newAnces, res$id), function(i) seq.dupli(res$dna[i,], 1)))
            ## molecular clock / time unit
            newSeq <- t(sapply(match(newAnces, res$id), function(i) seq.dupli(res$dna[i,], t-res$onset[match(newAnces, res$id)])))
            res$dna <- rbind(res$dna, newSeq)
        }


        ## IMPORTED CASES ##
        ## number of imported cases
        nbImpCases <- rpois(1, rate.import.case)
        if(nbImpCases>0){
            ## dates of imported cases
            res$onset <- c(res$onset, rep(t, nbImpCases))

            ## ancestries of the imported cases
            res$ances <- c(res$ances, rep(NA, nbImpCases))

            ## id of the imported cases
            newId <- seq(N+1, by=1, length=nbImpCases)
            res$id <- c(res$id, newId)

            ## status of new cases
            res$status[newId] <- "I"

            ## spatial coord of the new id
            if(spatial){
                newXy <- matrix(runif(nbImpCases*2, min=0, max=area.size), ncol=2)
                res$xy <- rbind(res$xy, newXy)
                res$inf.xy <- rbind(res$inf.xy, newXy) # set coords at infection

            }

            ## group of the imported cases
            res$group <- c(res$group, choose.group(nbImpCases))

            ## dna sequences of the new infections
            newSeq <- t(sapply(1:nbImpCases, function(i) seq.dupli(EVE, diverg.import)))
            res$dna <- rbind(res$dna, newSeq)
        }


        ## set recovered status ##
        res$status[res$id[(t-res$onset) >= t.clear]] <- "R"

        ## update nb of infected, recovered, etc.
        res$dynam$nrec[t+1] <- sum(res$status=="R")
        res$dynam$ninf[t+1] <- sum(res$status=="I")
        res$dynam$nsus[t+1] <- sum(res$status=="S")
    } # end for


    ## SHAPE AND RETURN OUTPUT ##
    ## data need to be reordered so that res$id is 1:res$n
    res$n <- nrow(res$dna)
    res$ances <- match(res$ances, res$id)
    res$id <- 1:res$n
    res$xy <- res$inf.xy # don't keep entire distribution, not right order anymore anyway
    res$inf.xy <- NULL # simpler to just call coords 'xy'
    res$status <- NULL # we don't need this

    findNmut <- function(i){
        if(!is.na(res$ances[i]) && res$ances[i]>0){
            out <- dist.dna(res$dna[c(res$id[i],res$ances[i]),], model="raw")*ncol(res$dna)
        } else {
            out <- NA
        }
        return(out)
    }

    ##res$nmut <- sapply(1:res$n, function(i) dist.dna(res$dna[c(res$id[i],res$ances[i]),], model="raw"))*ncol(res$dna)
    res$nmut <- sapply(1:res$n, function(i) findNmut(i))
    res$ngen <- rep(1, length(res$ances)) # number of generations
    res$call <- match.call()
    ## if(tree){
    ##     res$tree <- fastme.ols(dist.dna(res$dna, model="TN93"))
    ##     res$tree <- root(res$tree,"1")
    ## }


    class(res) <- "simOutbreak"
    return(res)

} # end simOutbreak







#' @rdname simOutbreak
#' @export
#'
print.simOutbreak <- function(x, ...){

    cat("\t\n=========================")
    cat("\t\n=   simulated outbreak  =")
    cat("\t\n=  (simOutbreak object) =")
    cat("\t\n=========================\n")

    cat("\nSize :", x$n,"cases (out of", x$dynam$nsus[1],"susceptible hosts)")
    cat("\nGenome length :", ncol(x$dna),"nucleotids")
    cat("\nDate range :", min(x$onset),"-",max(x$onset))
    cat("\nGroup distribution:")
    print(table(x$group))

    cat("\nContent:\n")
    print(names(x))

    return(NULL)
} # end print.simOutbreak







#' @rdname simOutbreak
#' @export
#'
"[.simOutbreak" <- function(x,i,j,drop=FALSE){
    res <- x
    ## trivial subsetting ##
    res$dna <- res$dna[i,,drop=FALSE]
    res$id <- res$id[i]
    res$onset <- res$onset[i]
    res$group <- res$group[i]
    res$n <- nrow(res$dna)
    res$nmut <- x$nmut[i]
    res$ngen <- x$ngen[i]
    res$status <- x$status[i]
    res$xy <- x$xy[i,,drop=FALSE]

    ## non-trivial subsetting ##
    res$ances <- res$ances[i]
    toFind <- !res$ances %in% res$id

    ## function to find Most Recent Ancestor (MRA)
    findMRA <- function(id){
        out <- list(ances=x$ances[id], nmut=x$nmut[id], ngen=x$ngen[id])
        if(is.na(out$ances)) return(out)
        while(!is.na(out$ances) && !out$ances %in% res$id){
            out$nmut <- out$nmut + x$nmut[out$ances]
            out$ngen <- out$ngen + x$ngen[out$ances]
            out$ances <- x$ances[out$ances]
        }
        if(is.na(out$ances)) return(list(ances=NA, nmut=NA, ngen=1))
        return(out)
    }

    ## browse the tree to find indirect ancestries
    if(any(toFind)){
        temp <- sapply(res$id[toFind], findMRA)
        res$ances[toFind] <- unlist(temp["ances",])
        res$nmut[toFind] <- unlist(temp["nmut",])
        res$ngen[toFind] <- unlist(temp["ngen",])
    }

    return(res)
} # end subsetting method





##################
## labels.simOutbreak
##################
#' @rdname simOutbreak
#' @export
#'
labels.simOutbreak <- function(object, ...){
    return(object$id)
}





#' @rdname simOutbreak
#' @export
#'
as.igraph.simOutbreak <- function(x, edge.col="black", col.edge.by="dist", vertex.col="gold",
                                  edge.col.pal=NULL, annot=c("dist","n.gen"), sep="/", ...){
    ## if(!require(igraph)) stop("package igraph is required for this operation")
    ## if(!require(ape)) stop("package ape is required for this operation")
    if(!inherits(x,"simOutbreak")) stop("x is not a tTree object")
    ## if(!require(adegenet)) stop("adegenet is required")
    if(!col.edge.by %in% c("dist","n.gen","prob")) stop("unknown col.edge.by specified")

    ## GET DAG ##
    from <- as.character(x$ances)
    to <- as.character(x$id)
    isNotNA <- !is.na(from) & !is.na(to)
    vnames <- unique(c(from,to))
    vnames <- vnames[!is.na(vnames)]
    dat <- data.frame(from,to,stringsAsFactors=FALSE)[isNotNA,,drop=FALSE]
    out <- graph.data.frame(dat, directed=TRUE, vertices=data.frame(names=vnames))

    ## from <- as.character(x$ances)
    ## to <- as.character(x$id)
    ## dat <- data.frame(from,to,stringsAsFactors=FALSE)[!is.na(x$ances),]
    ## out <- graph.data.frame(dat, directed=TRUE)

    ## SET VERTICE INFO ##
    ## labels
    V(out)$label <- V(out)$name

    ## dates
    names(x$onset) <- x$id
    V(out)$date <- x$onset[V(out)$name]

    ## ## groups
    ## names(x$group) <- x$id
    ## V(out)$group <- x$group[V(out)$name]

    ## colors
    V(out)$color <- vertex.col
    ## V(out)$color <- fac2col(factor(V(out)$group), col.pal=vertex.col.pal)


    ## SET EDGE INFO ##
    ## genetic distances to ancestor
    E(out)$dist <- x$nmut[!is.na(x$ances)]

    ## number of generations to ancestor
    E(out)$ngen <- x$ngen[!is.na(x$ances)]

    ## colors
    if(is.null(edge.col.pal)){
        edge.col.pal <- function(n){
            return(grey(seq(0.75,0,length=n)))
        }
    }
    if(col.edge.by=="dist") edge.col <- num2col(E(out)$dist, col.pal=edge.col.pal, x.min=0, x.max=1)

    ## labels
    n.annot <- sum(annot %in% c("dist","n.gen"))
    lab <- ""
    if(!is.null(annot) && n.annot>0){
        if("dist" %in% annot) lab <- E(out)$dist
        if("n.gen" %in% annot) lab <- paste(lab, x$ngen[!is.na(x$ances)], sep=sep)
    }
    lab <- sub(paste("^",sep,sep=""),"",lab)
    E(out)$label <- lab

    ## SET LAYOUT ##
    attr(out, "layout") <- layout.fruchterman.reingold(out, params=list(minx=V(out)$date, maxx=V(out)$date))

    return(out)
} # end as.igraph.simOutbreak







#' @rdname simOutbreak
#' @export
#'
plot.simOutbreak <- function(x, y=NULL, edge.col="black", col.edge.by="dist", vertex.col="gold",
                              edge.col.pal=NULL, annot=c("dist","n.gen"), sep="/", ...){
    ## if(!require(igraph)) stop("igraph is required")
    ## if(!require(adegenet)) stop("adegenet is required")
    if(!inherits(x,"simOutbreak")) stop("x is not a simOutbreak object")
    if(!col.edge.by %in% c("dist","n.gen")) stop("unknown col.edge.by specified")

    ## get graph ##
    g <- as.igraph(x, edge.col=edge.col, col.edge.by=col.edge.by, vertex.col=vertex.col,
                   edge.col.pal=edge.col.pal, annot=annot, sep=sep)

     ## make plot ##
    plot(g, layout=attr(g,"layout"), ...)

    ## return graph invisibly ##
    return(invisible(g))

} # end plot.simOutbreak






## #####################
## ## seqTrack.simOutbreak
## #####################
## seqTrack.simOutbreak <- function(x, best=c("min","max"), prox.mat=NULL, ...){
##     myX <- dist.dna(x$dna, model="raw")
##     x.names <- labels(x)
##     x.dates <- as.POSIXct(x)
##     seq.length <- ncol(x$dna)
##     myX <- myX * seq.length
##     myX <- as.matrix(myX)
##     prevCall <- as.list(x$call)
##     if(is.null(prevCall$mu)){
##         mu0 <- 0.0001
##     } else {
##         mu0 <- eval(prevCall$mu)
##     }
##     res <- seqTrack(myX, x.names=x.names, x.dates=x.dates, best=best, prox.mat=prox.mat,...)
##     return(res)
## }






## ########################
## ## as.seqTrack.simOutbreak
## ########################
## as.seqTrack.simOutbreak <- function(x){
##     ## x.ori <- x
##     ## x <- na.omit(x)
##     toSetToNA <- x$onset==min(x$onset)
##     res <- list()
##     res$id <- labels(x)
##     res <- as.data.frame(res)
##     res$ances <- x$ances
##     res$ances[toSetToNA] <- NA
##     res$weight <- 1 # ??? have to recompute that...
##     res$weight[toSetToNA] <- NA
##     res$date <- as.POSIXct(x)[labels(x)]
##     res$ances.date <- as.POSIXct(x)[x$ances]

##     ## set results as indices rather than labels
##     res$ances <- match(res$ances, res$id)
##     res$id <- 1:length(res$id)

##     ## SET CLASS
##     class(res) <- c("seqTrack", "data.frame")

##     return(res)
## }








## ###################
## ## sample.simOutbreak
## ###################
## sample.simOutbreak <- function(x, n){
## ##sample.simOutbreak <- function(x, n, rDate=.rTimeSeq, arg.rDate=NULL){
##     ## EXTRACT THE SAMPLE ##
##     res <- x[sample(1:nrow(x$dna), n, replace=FALSE)]


##     ## RETRIEVE SOME PARAMETERS FROM HAPLOSIM CALL
##     prevCall <- as.list(x$call)
##     if(!is.null(prevCall$mu)){
##         mu0 <- eval(prevCall$mu)
##     } else {
##         mu0 <- 1e-04
##     }

##     if(!is.null(prevCall$dna.length)){
##         L <- eval(prevCall$dna.length)
##     } else {
##         L <- 1000
##     }

##     ## truedates <- res$onset
##     ## daterange <- diff(range(res$onset,na.rm=TRUE))

##     ## if(identical(rDate,.rTimeSeq)){
##     ##     sampdates <- .rTimeSeq(n=length(truedates), mu=mu0, L=L, maxNbDays=daterange/2)
##     ## } else{
##     ##     arg.rDate$n <- n
##     ##     sampdates <- do.call(rDate, arg.rDate)
##     ## }
##     ## sampdates <- truedates + abs(sampdates)

##     return(res)
## } # end sample.simOutbreak

