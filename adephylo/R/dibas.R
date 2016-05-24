#########
## dibas ('distance-based group assignment')
#########

dibas <- function (x, ...) UseMethod("dibas")





################
## dibas.matrix
################
dibas.matrix <- function(x, grp, method=c("default","leaveOneOut"), ...){
    method <- match.arg(method)
    ## DECLARE SOME VARIABLES, HANDLE ARGUMENTS ##
    grp <- factor(grp)
    K <- length(LEV <- levels(grp))
    N <- nrow(x)


    ## AUXILIARY FUNCTIONS ##
    ## COMPUTE LOG AND AVOIDS -INF
    logprob <- function(prob){
        res <- log(prob)
        res[res< -1e20] <- -1e20
        return(res)
    }

    ## FUNCTION TO GET A VECTOR OF PAIRWISE DISTANCE WITHIN ONE GROUP
    ## M: matrix of distances
    ## fac: factor
    ## val: level of the chosen group
    getdist.within.grp <- function(M, fac, val){ # val is one level of fac
        temp <- M[fac==val,fac==val]
        return(temp[lower.tri(temp)])
    }


    ## FUNCTION TO GET LIST OF VECTORS OF PAIRWISE DISTANCES WITHIN GROUP, FOR EVERY GROUP
    getdist.within.allgrp <- function(M, fac){
        res <- lapply(LEV, function(e) getdist.within.grp(M, fac, e))
        names(res) <- LEV
        return(res)
    }


    ## FUNCTION TO GET THE DISTANCES OF AN INDIVIDUAL TO THE GROUPS
    getdist.indiv <- function(i){
        return(split(x[i,-i],grp[-i]))
    }

    ## FUNCTION TO COMPUTE MEMBERSHIP PROBA FOR ONE INDIV
    ## i: index of an individual
    ## distrib.param[1,]: vector of the means of with-grp distance distributions
    ## distrib.param[2,]: vector of the sds of with-grp distance distributions
    getproba.ind <- function(i, distrib.param){
        temp <- getdist.indiv(i)
        out <- sapply(1:K, function(k) mean(logprob(dnorm(temp[[k]], distrib.param[1,k],distrib.param[1,k]))))
        return(exp(out)/sum(exp(out)))
    }



    ## CORE COMPUTATIONS ##
    ## DEFAULT: DENSITY BASED ON ENTIRE SAMPLE ##
    if(method=="default"){
        ## get distance data for each group
        temp <- getdist.within.allgrp(x, grp)

        ## parameter of the group distributions
        ## matrix of within-group dist: col=grp, row1=mean,row2=sd
        distrib.param <- sapply(temp, function(e) return(c(mean(e,na.rm=TRUE),sd(e,na.rm=TRUE))))

        ## result: row=indiv, col=groups, values=membership proba
        prob <- t(sapply(1:N, getproba.ind, distrib.param))
    }


    ## LEAVEONEOUT: DENSITY EXCLUDES THE INDIV FOR WHICH PROBA IS SEEKED ##
    if(method=="leaveOneOut"){
        ## get within-group distance data for each individual
        temp <- lapply(1:N, function(i) getdist.within.allgrp(x[-i,-i], grp[-i])) # grp density data for each individual

        ## parameter of the group distributions
        ## list of matrices, one per individual
        ## matrix of within-group dist: col=grp, row1=mean,row2=sd
        distrib.param <- lapply(1:N, function(i) sapply(temp[[i]], function(e) return(c(mean(e,na.rm=TRUE),sd(e,na.rm=TRUE)))))

        ## result: row=indiv, col=groups, values=membership proba
        prob <- t(sapply(1:N, function(i) getproba.ind(i, distrib.param[[i]])))
    }


    ## SHAPE MEMBERSHIP PROBABILITIES MATRIX ##
    colnames(prob) <- LEV
    rownames(prob) <- rownames(x)


    ## FIND GROUP ASSIGNMENTS ##
    temp <- factor(colnames(prob)[apply(prob,1, which.max)])
    annot <- rep(" ", N)
    annot[as.character(grp)!=as.character(temp)] <- "!"
    groups <- data.frame(observed=grp, inferred=temp, annot=annot)
    ##rownames(groups) <- rownames(prob)


    ## BUILD / RETURN RESULT ##
    ## get proportion of correct assignment
    propcorrect <- mean(annot==" ")
    propcorrect.bygroup <- tapply(annot==" ", grp, mean)

    ## get summary of assignments
    grp.tab <- table(observed=groups[,1], assigned=groups[,2])

    ## get assignability
    ## i.e. how many times better than at random is assignment?
    ## 0 = grp very unlikely
    ## 1 = assignment no better than at random
    ## >1 = better than random (e.g. 2 = twice as better as at random)
    temp <- table(grp)/N
    probActualGrp <- sapply(1:N, function(i) prob[i, as.character(grp[i])])
    assign.idx <-     probActualGrp / as.numeric(temp[as.character(grp)])
    assignStat <- list(assign.idx=assign.idx, mean=mean(assign.idx), grp.mean=tapply(assign.idx,grp,mean))


    ##res <- list(prob=prob,groups=groups, mean.correct=propcorrect, prop.correct=propcorrect.bygroup)
    res <- list(prob=prob, groups=groups, mean.ok=propcorrect, grp.tab=grp.tab, assignStat=assignStat)

    return(res)
} # end dibas.matrix








################
## dibas.vector
################
##
## in this one, one distance to a reference point
## is used to defined group membership probabilities
##
dibas.vector <- function(x, grp, method=c("default","leaveOneOut"), n.items=NULL, ...){
    method <- match.arg(method)

    ## DECLARE SOME VARIABLES, HANDLE ARGUMENTS ##
    grp <- factor(grp)
    K <- length(LEV <- levels(grp))
    N <- length(x)
    if(!is.null(n.items)){
        n.items <- round(n.items)
        if(length(n.items)!=N) stop("n.items has a wrong length")
        if(any(n.items<1)) stop("values in n.items cannot be less than 1")
        x <- rep(x, n.items)
        grp <- rep(grp, n.items)
    }


    ## AUXILIARY FUNCTIONS ##
    ## COMPUTE LOG AND AVOIDS -INF
    logprob <- function(prob){
        res <- log(prob)
        res[res< -1e20] <- -1e20
        return(res)
    }


    ## FUNCTION TO COMPUTE MEMBERSHIP PROBA FOR ONE INDIV
    ## i: index of an individual
    ## distrib.mu: vector of the means of with-grp distance distributions
    ## distrib.sigma: vector of the sds of with-grp distance distributions
    getproba.ind <- function(i, leaveOneOut){
        if(leaveOneOut){
             distrib.mu  <- tapply(x[-i], grp[-i], mean, na.rm=TRUE)
             distrib.sigma <- tapply(x[-i], grp[-i], sd, na.rm=TRUE)
         } else {
             distrib.mu  <- tapply(x, grp, mean, na.rm=TRUE)
             distrib.sigma <- tapply(x, grp, sd, na.rm=TRUE)
         }
        out <- sapply(1:K, function(k) logprob(dnorm(x[i], distrib.mu[k],distrib.sigma[k])))
        return(exp(out)/sum(exp(out)))
    }



    ## CORE COMPUTATIONS ##
    prob <- t(sapply(1:length(x), getproba.ind, leaveOneOut=method=="leaveOneOut"))

    ## SHAPE MEMBERSHIP PROBABILITIES MATRIX ##
    colnames(prob) <- LEV
    rownames(prob) <- rownames(x)


    ## FIND GROUP ASSIGNMENTS ##
    temp <- factor(colnames(prob)[apply(prob,1, which.max)])
    annot <- rep(" ", N)
    annot[as.character(grp)!=as.character(temp)] <- "!"
    groups <- data.frame(observed=grp, inferred=temp, annot=annot)
    ##rownames(groups) <- rownames(prob)


    ## BUILD / RETURN RESULT ##
    ## get proportion of correct assignment
    propcorrect <- mean(annot==" ")
    propcorrect.bygroup <- tapply(annot==" ", grp, mean)

    ## get summary of assignments
    grp.tab <- table(observed=groups[,1], assigned=groups[,2])

    ## get assignability
    ## i.e. how many times better than at random is assignment?
    ## 0 = grp very unlikely
    ## 1 = assignment no better than at random
    ## >1 = better than random (e.g. 2 = twice as better as at random)
    temp <- table(grp)/N
    probActualGrp <- sapply(1:N, function(i) prob[i, as.character(grp[i])])
    assign.idx <-     probActualGrp / as.numeric(temp[as.character(grp)])
    assignStat <- list(assign.idx=assign.idx, mean=mean(assign.idx), grp.mean=tapply(assign.idx,grp,mean))


    ##res <- list(prob=prob,groups=groups, mean.correct=propcorrect, prop.correct=propcorrect.bygroup)
    res <- list(prob=prob, groups=groups, mean.ok=propcorrect, grp.tab=grp.tab, assignStat=assignStat)

    return(res)
} # end dibas.vector






###############
## dibas.phylo
###############
dibas.phylo <- function(x, grp, method=c("default","leaveOneOut"), fromRoot=FALSE, metric=c("Abouheif", "nNodes", "patristic", "sumDD"),
                        n.items=NULL, ...){
    ## if(!require(ape)) stop("ape package is required")
    if(!inherits(x,"phylo")) stop("x is not a phylo object")

    metric <- match.arg(metric)

    if(fromRoot){
        res <- dibas.vector(distRoot(x, method=metric), grp=grp, method=method, n.items=n.items)
    } else {
        res <- dibas(distTips(x, method=metric), grp=grp, method=method)
    }

    return(res)
} # end dibas.phylo






##############
## dibas.dist
##############
dibas.dist <- function(x, grp, method=c("default","leaveOneOut"), ...){

    res <- dibas.matrix(as.matrix(x), grp, method)

    return(res)
} # end dibas.phylo







##############################
## simulate data with groups
##############################
simDatGroups <- function(k=2, p=1000, n=10, mu=0, sigma=1){
    ## RECYCLE ARGUMENTS ##
    n <- rep(n, length=k)
    mu <- rep(mu, length=k)
    sigma <- rep(sigma, length=k)


    ## GENERATE DATA ##
    dat <- list()
    for(i in 1:k){
        dat[[i]] <- replicate(p, rnorm(n[i], mu[i], sigma[i]))
    }

    dat <- Reduce(rbind,dat)
    rownames(dat) <- paste("ind", 1:nrow(dat))

    ## SHAPE OUTPUT ##
    grp <- factor(paste("grp", rep(1:k, n)))
    res <- list(dat=dat, grp=grp)
    return(res)
} # end simDatGroups
















########## OLD CODE, USING A DIFFERENT APPROACH ###########
## THIS WAS USING A KERNEL APPROX OF THE DISTRIBUTION OF
## WITHIN GROUP DISTANCES. NOT WORKING BECAUSE THERE COULD BE
## MORE THAN ONE MODE, SO GROUPS COULD BE PRETTY SPLIT ACROSS
## THE PHYLOGENY
##
## ##############
## ## dibas
## ##############
## dibas <- function(x, grp, method=c("default","leaveOneOut","bootstrap"), n.dens=4096, plot=TRUE,
##                       warn.lab=FALSE, dat=NULL, FUN=NULL, n.boot=10, ...){
##     if(!require(ape)) stop("ape package is required")
##     if(!inherits(x,"phylo")) stop("x is not a phylo object")
##     method <- match.arg(method)

##     if(method=="bootstrap" && (is.null(dat) || is.null(FUN))) stop("dat and FUN must be provided for the bootstrap procedure")
##     if(warn.lab && !is.null(dat) && !identical(x$tip.label,rownames(dat))) warning("Tip labels in x and rownames of dat differ \nplease make sure the same order is used in x, grp, and dat")

##     ## DECLARE SOME VARIABLES, HANDLE ARGUMENTS ##
##     grp <- factor(grp)
##     K <- length(LEV <- levels(grp))
##     N <- length(x$tip.label)
##     D <- cophenetic.phylo(x)
##     THRES <- 1e-320 # densities < THRES will be set to THRES to avoid log(x)=-Inf


##     ## RE-ORDER GRP AND DATA MATRIX AFTER TIP LABELS ##
##     if(!is.null(dat)){
##         if(is.null(rownames(dat))) rownames(dat) <- x$tip.label
##         if(!all(x$tip.label %in% rownames(dat))) stop("some tips do not have data matching their label")
##         grp <- grp[match(x$tip.label, rownames(dat))] # grp is assumed to be in the same order as 'dat'
##         dat <- dat[x$tip.label,,drop=FALSE]
##     }

##     #### AUXILIARY FUNCTIONS ####
##     ## FUNCTION TO ESTIMATE A DENSITY AT A SERIES OF POINTS ##
##     compute.dens <- function(dens, values){
##         pred.y <- double(n <- length(values))
##         return(.C("predict_density", dens$x, dens$y, length(dens$x), as.double(values), pred.y, n, PACKAGE="adephylo")[[5]])
##     }


##     ## FUNCTION TO GET A VECTOR OF PAIRWISE DISTANCE WITHIN ONE GROUP ##
##     getdist.within.grp <- function(M, fac, val){ # val is one level of fac
##         temp <- M[fac==val,fac==val]
##         return(temp[lower.tri(temp)])
##     }


##     ## FUNCTION TO GET A VECTOR OF PAIRWISE DISTANCES WITHIN GROUP, FOR ALL GROUPS ##
##     getdist.within.allgrp <- function(M, fac){
##         res <- lapply(LEV, function(e) getdist.within.grp(M, fac, e))
##         names(res) <- LEV
##         return(res)
##     }


##     ## FUNCTION TO GET PROBA FOR ONE INDIV / ONE GROUP ##
##     getprob.ind <- function(i, g, dens.per.grp){ # i: idx of indiv; g: idx of a group
##         temp <- 1:ncol(D)
##         dens <- compute.dens(dens.per.grp[[g]], D[i,grp==LEV[g] & temp!=i])
##         dens[dens < THRES] <- THRES
##         res <- exp(mean(log(dens)))
##         return(res)
##     }


##     ## FUNCTION TO GET PROBA FOR ALL INDIV / ONE GROUP ##
##     if(method=="leaveOneOut"){
##         getprob.grp <- function(g, dens.per.ind.grp){ # g: idx of a group; dens.per.ind.grp is a list giving grp density for each indiv
##             return(sapply(1:N, function(i) getprob.ind(i,g,dens.per.ind.grp[[i]])))
##         }
##     } else {
##         getprob.grp <- function(g, dens.per.grp){ # g: idx of a group
##             return(sapply(1:N, function(i) getprob.ind(i,g,dens.per.grp)))
##         }
##     }


##     ## FUNCTION TO GET A BOOTSTRAPPED TREE AND MATCHING GRP ##
##     getboot.tree.grp <- function(){
##         samp <- sample(1:N,replace=TRUE)
##         tre <- FUN(dat[samp,,drop=FALSE])
##         newgrp <- factor(grp[samp])
##         return(list(tre=tre, grp=newgrp))
##     }


##     #### CORE COMPUTATIONS ####
##     ## DEFAULT: DENSITY BASED ON SAMPLE ##
##     if(method=="default"){
##         dens.dat <- getdist.within.allgrp(D, grp) # density data for each group
##         list.dens <- lapply(dens.dat, density, n=n.dens, ...) # densities for each group
##     }


##     ## LEAVEONEOUT: GRP DENSITY EXCLUDES THE INDIV FOR WHICH PROBA IS SEEKED ##
##     if(method=="leaveOneOut"){
##         dens.dat <- lapply(1:N, function(i) getdist.within.allgrp(D[-i,-i], grp[-i])) # grp density data for each individual
##         list.dens <- lapply(1:N, function(i) lapply(dens.dat[[i]], density, n=n.dens, ...)) # densities for each group
##     }


##     ## BOOTSTRAP: DENSITY BASED ON BOOTSTRAPPIN INDIVIDUALS ##
##     if(method=="bootstrap"){
##         ## GET BOOTSTRAPPED TREES ##
##         list.trees.grp <- lapply(1:n.boot, function(i) getboot.tree.grp())


##         ## GET WITHIN-GROUP DISTANCES FOR EACH BOOTSTRAP SAMPLE ##
##         list.D <- lapply(list.trees.grp, function(e) cophenetic.phylo(e$tre))
##         temp <- lapply(1:n.boot, function(i) getdist.within.allgrp(list.D[[i]], list.trees.grp[[i]]$grp)) # for each replicate, list of distances within for each grp


##         ## GET DENSITIES FOR EACH GROUP ##
##         dens.dat <- lapply(LEV, function(onegroup) unlist(lapply(temp, function(e) e[[onegroup]]))) # density data for each group
##         list.dens <- lapply(dens.dat, density, n=n.dens, ...) # densities for each group
##     }


##     ## PLOT DENSITIES ##
##     if(method != "leaveOneOut" && plot){
##         find.mfrow <- function(i) {
##             nrow <- ceiling(sqrt(i))
##             ncol <- ceiling(i/ceiling(sqrt(i)))
##             return(c(nrow,ncol))
##         }
##         par(mfrow = find.mfrow(K))
##         for(i in 1:K){
##             plot(list.dens[[i]], main=paste("Group:",LEV[i]),xlab="Within-group pairwise distance",ylab="Density", col="blue")
##             points(dens.dat[[i]], rep(0,length(dens.dat[[i]])), pch="|", col="blue")
##         }
##     }


##     ## COMPUTE MEMBERSHIP PROBABILITIES ##
##     prob <- matrix(unlist(lapply(1:K, getprob.grp, list.dens)), ncol=K)
##     prob <- prop.table(prob,1)
##     colnames(prob) <- LEV
##     rownames(prob) <- x$tip.label


##     ## FIND GROUP ASSIGNMENTS ##
##     temp <- factor(colnames(prob)[apply(prob,1, which.max)])
##     annot <- rep(" ", N)
##     annot[as.character(grp)!=as.character(temp)] <- "!"
##     groups <- data.frame(observed=grp, inferred=temp, annot=annot)
##     ##rownames(groups) <- rownames(prob)


##     ## BUILD / RETURN RESULT ##
##     propcorrect <- mean(annot==" ")
##     ## propcorrect.bygroup <- tapply(annot==" ", grp, mean)
##     assignability <- mean((apply(prob,1,max)-.5)/.5)
##     ##res <- list(prob=prob,groups=groups, mean.correct=propcorrect, prop.correct=propcorrect.bygroup)
##     res <- list(prob=prob, groups=groups, assigndex=assignability, mean.correct=propcorrect)

##     return(res)
## } # end dibas
