
######################################################################################################################################
######################################################################################################################################
### Lewis MKv model for phylogenetic inference based on morphological characters that are all variable
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu

lewisMkv <- function(phy, data, include.gamma=FALSE, ngammacats=4, include.beta=FALSE, exclude.sites=NULL, max.tol=.Machine$double.eps^0.25, ncores=NULL) {
    
    if(include.beta == TRUE){
        stop("Beta distribution for allowing unequal stationary frequencies is not yet implemented.")
    }
    
    # Sort according to the order of the tip labels:
    data.sorted <- data[phy$tip.label,]
    if(!is.null(exclude.sites)){
        data.sorted <- data.sorted[!sequence(dim(data.sorted)[2]) %in% exclude.sites]
    }
    cat("Getting starting branch lengths using a Rogers-Swoffordish approximation....", "\n")
    # Starting values for branch lengths:
    new.phy <- RogersSwoffordish(phy, data.sorted)
    new.phy$edge.length[new.phy$edge.length==0] = 1e-8
    
    cat("Finished. Optimizing model parameters...", "\n")
    opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = "100000", "ftol_rel" = max.tol)
    obj <- NULL
    if(include.gamma == TRUE){
        ip <- c(1, new.phy$edge.length)
        upper = c(1, rep(log(5), length(new.phy$edge.length)))
        lower = rep(-21, length(ip))
        search.results <- nloptr(x0=log(ip), eval_f=GetTotalLikelihood, ub=upper, lb=lower, opts=opts, phy=new.phy, data=data.sorted, include.gamma=include.gamma, ngammacats=ngammacats, include.beta=include.beta, ncores=ncores)
        loglik <- -search.results$objective
        shape.est <- exp(search.results$solution[1])
        new.phy$edge.length <- exp(search.results$solution[2:length(search.results$solution)])
        np = 1 + length(new.phy$edge.length)
    }else{
        ip <- c(new.phy$edge.length)
        upper = c(rep(log(5), length(new.phy$edge.length)))
        lower = rep(-21, length(ip))
        search.results <- nloptr(x0=log(ip), eval_f=GetTotalLikelihood, ub=upper, lb=lower, opts=opts, phy=new.phy, data=data.sorted, include.gamma=include.gamma, ngammacats=ngammacats, include.beta=include.beta, ncores=ncores)
        loglik <- -search.results$objective
        shape.est <- NULL
        new.phy$edge.length <- exp(search.results$solution[1:length(search.results$solution)])
        np = 1 + length(new.phy$edge.length)
    }
    cat("Finished.", "\n")
    obj = list(np=np, loglik = loglik, AIC = -2*loglik+2*np, gamma.shape=shape.est, phy=new.phy, data=data.sorted, opts=opts)
    class(obj) = "lewis.mkv"
    return(obj)
}


######################################################################################################################################
######################################################################################################################################
### Print function for the lewis.mkv class:
######################################################################################################################################
######################################################################################################################################

print.lewis.mkv <- function(x,...){
    if(!is.null(x$gamma.shape)){
        output<-data.frame(x$loglik,x$AIC,dim(x$data)[1],dim(x$data)[2], x$gamma.shape, row.names="")
        names(output)<-c("-lnL","AIC", "ntax", "nsites", "gamma")
        cat("\nFit\n")
        print(output)
        cat("\n")
    }else{
        output<-data.frame(x$loglik,x$AIC,dim(x$data)[1],dim(x$data)[2], row.names="")
        names(output)<-c("-lnL","AIC", "ntax", "nsites")
        cat("\nFit\n")
        print(output)
        cat("\n")
    }
    print(x$phy)
    cat("\n")
}


######################################################################################################################################
######################################################################################################################################
### Functions for getting the likelihood
######################################################################################################################################
######################################################################################################################################

GetTotalLikelihood <- function(x, phy, data, include.gamma=FALSE, ngammacats=4, include.beta=FALSE, ncores=NULL) {
    x = exp(x)
    if(include.gamma == TRUE){
        shape = x[1]
        x = x[-1]
    }
    phy$edge.length <- x
    nsites <- dim(data)[2]
    if(include.gamma==TRUE){
        rates.k <- DiscreteGamma(shape, ngammacats)
        final.likelihood.mat <- matrix(0, nrow=ngammacats, ncol=nsites)
        for(k in sequence(ngammacats)){
            final.likelihood.mat[k,] <- -GetSiteLikelihoods(phy=phy, data=data, gamma.rate=rates.k[k], beta.freqs=NULL, ncores=ncores)
        }
        likelihood <- -sum(log(colMeans(exp(final.likelihood.mat))))
    }else{
        final.likelihood <- GetSiteLikelihoods(phy=phy, data=data, gamma.rate=1, beta.freqs=NULL, ncores=ncores)
        likelihood <- sum(final.likelihood)
    }
    return(likelihood)
}


GetSiteLikelihoods <- function(phy, data, gamma.rate=1, beta.freqs=NULL, ncores=NULL) {
    nsites <- dim(data)[2]
    if(is.null(ncores)){
        final.likelihood.vector <- rep(NA, nsites)
        for (site.index in sequence(nsites)) {
            cache <- FormatForLikelihood(phy=phy, data=data, charnum=site.index)
            final.likelihood.vector[site.index] <- LikelihoodCalculation(phy=phy, liks=cache$liks, Q=cache$Q, gamma.rate=gamma.rate)
        }
        return(final.likelihood.vector)
    }else{
        MultiCoreLikelihood <- function(site.index){
            cache <- FormatForLikelihood(phy=phy, data=data, charnum=site.index)
            final.likelihood.vector.tmp <- LikelihoodCalculation(phy=phy, liks=cache$liks, Q=cache$Q, gamma.rate=gamma.rate)
            return(final.likelihood.vector.tmp)
        }
        final.likelihood.vector <- unlist(mclapply(1:nsites, MultiCoreLikelihood, mc.cores=ncores))
        return(final.likelihood.vector)
    }
}


######################################################################################################################################
######################################################################################################################################
### Likelihood calculator
######################################################################################################################################
######################################################################################################################################

LikelihoodCalculation <- function(phy, liks, Q, gamma.rate){
    
    # Normalize the matrix because we are also estimating branch lengths:
    Q <- t(Q) * rep(1/dim(Q)[2], dim(Q)[2])
    scale.factor <- -sum(diag(Q) * rep(1/dim(Q)[2], dim(Q)[2]))
    Q.scaled <- Q * (1/scale.factor)
    Q.scaled <- Q.scaled * gamma.rate
    #####################################################################
    
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    TIPS <- 1:nb.tip
    liks.dummy <- liks
    liks.dummy[TIPS,1] = 1
    liks.dummy[TIPS,2:dim(Q.scaled)[1]] = 0
    comp <- numeric(nb.tip + nb.node)
    comp.dummy <- comp
    phy <- reorder(phy, "pruningwise")
    # Obtain an object of all the unique ancestors
    anc <- unique(phy$edge[,1])
    for (i  in seq(from = 1, length.out = nb.node)) {
        # The ancestral node at row i is called focal
        focal <- anc[i]
        # Get descendant information of focal
        desRows<-which(phy$edge[,1]==focal)
        desNodes<-phy$edge[desRows,2]
        v <- 1
        v.dummy <- 1
        for (desIndex in sequence(length(desRows))){
            v <- v*expm(Q.scaled * phy$edge.length[desRows[desIndex]]) %*% liks[desNodes[desIndex],]
            v.dummy <- v.dummy * expm(Q.scaled * phy$edge.length[desRows[desIndex]]) %*% liks.dummy[desNodes[desIndex],]
        }
        comp[focal] <- sum(v)
        comp.dummy[focal] <- sum(v.dummy)
        liks[focal, ] <- v/comp[focal]
        liks.dummy[focal, ] <- v.dummy/comp.dummy[focal]
    }
    # Specifies the root:
    root <- nb.tip + 1L
    # If any of the logs have NAs restart search:
    if(is.nan(sum(log(comp[-TIPS]))) || is.na(sum(log(comp[-TIPS])))){
        return(10000000000)
    }
    else{
        flat.root <- 1/dim(Q.scaled)[1]
        loglik.num <- (sum(log(comp[-TIPS])) + log(sum(exp(log(flat.root)+log(liks[root,])))))
        loglik.denom <- (sum(log(comp.dummy[-TIPS])) + log(sum(exp(log(flat.root)+log(liks.dummy[root,])))))
        loglik <- loglik.num  - log(1 - exp(loglik.denom))
        if(is.infinite(loglik)){
            return(10000000000)
        }
    }
    return(-loglik)
}



######################################################################################################################################
######################################################################################################################################
### Utility functions
######################################################################################################################################
######################################################################################################################################

FormatForLikelihood <- function(phy, data, charnum){
    k <- 1
    factored <- FactorDataLewis(data, charnum=charnum)
    nl <- ncol(factored)
    obj <- NULL
    nb.tip<-length(phy$tip.label)
    nb.node <- phy$Nnode
    Q <- matrix(1, nrow=nl, ncol=nl)
    diag(Q) <- 0
    diag(Q) <- -rowSums(Q)
    stateTable <- NULL # will hold 0s and 1s for likelihoods of each state at tip and deals with ambiguous characters
    for(column in 1:nl){
        stateTable <- cbind(stateTable,factored[,column])
    }
    colnames(stateTable) <- colnames(factored)
    ancestral <- matrix(0,nb.node,nl) # all likelihoods at ancestral nodes will be 0
    liks <- rbind(stateTable,ancestral) # combine tip likelihoods & ancestral likelihoods
    rownames(liks) <- NULL
    obj$liks <- liks
    obj$Q <- Q
    return(obj)
}


DiscreteGamma <- function (shape, ngammacats){
    quantiles <- qgamma((1:(ngammacats - 1))/ngammacats, shape = shape, rate = shape)
    return(diff(c(0, pgamma(quantiles * shape, shape + 1), 1)) * ngammacats)
}


# A function to find positions of ampersands for separating different states.
# Will allow character state to be greater than one character long.
FindAmpsLewis <- function(string, charnum){
    if(!is.character(string)) return(NULL)
    locs <- NULL # Will hold location values
    for(charnum in 1:nchar(as.character(string))){
        if(substr(string,charnum,charnum) == "&"){
            locs <- c(locs,charnum)
        }
    }
    return(locs)
}


# Function to make factored matrix as levels are discovered.
FactorDataLewis <- function(data,charnum){
    charcol <- charnum
    factored <- NULL # will become the matrix.  Starts with no data.
    lvls <- NULL
    numrows <- length(data[,charcol])
    missing <- NULL
    
    for(row in 1:numrows){
        currlvl <- NULL
        levelstring <- as.character(data[row,charcol])
        ampLocs <- FindAmpsLewis(levelstring, charnum)
        if(length(ampLocs) == 0){ #No ampersands, character is monomorphic
            currlvl <- levelstring
            if(currlvl == "?" || currlvl == "-"){ # Check for missing data
                missing <- c(missing,row) # add to list of taxa with missing values, will fill in entire row later
            }
            else { # Not missing data
                if(length(which(lvls == currlvl)) == 0){# encountered a level not seen yet
                    if(length(factored) == 0){ # Matrix is empty, need to create it
                        factored <- matrix(0,numrows,1)
                        colnames(factored) <- currlvl
                        rownames(factored) <- rownames(data)
                    } else { # matrix already exists, but need to add a column for the new level
                        zerocolumn <- rep(0,numrows)
                        factored <- cbind(factored, zerocolumn)
                        colnames(factored)[length(factored[1,])] <- currlvl
                    }
                    lvls <- c(lvls,currlvl) # add that level to the list
                } # already found this level in another state.  Set the value to one
                whichlvl <- which(lvls == currlvl) # this index number should correspond to the column number of the state
                factored[row,whichlvl] <- 1
            }
        } else { #At least one ampersand found, polymorphic character
            start <- 1
            numlvls <- length(ampLocs)+1
            for(part in 1:numlvls){
                # Pull out level from levelstring
                if(part <= length(ampLocs)){ # Havent reached the last state
                    currlvl <- substr(levelstring,start,(ampLocs[part]-1)) # pull out value between start and the location-1 of the next ampersand
                } else { # Final state in list
                    currlvl <- substr(levelstring,start,nchar(levelstring)) # pull out value between start and the last character of the string
                }
                if(currlvl == "?" || currlvl == "-"){ # Missing data, but polymorphic?
                    missing <- c(missing,row) # add to list of taxa with missing values, will fill in entire row later
                }
                else { # Not missing data
                    if(length(which(lvls == currlvl)) == 0){# encountered a level not seen yet
                        if(length(factored) == 0){ # Matrix is empty, need to create it
                            factored <- matrix(0,numrows,1)
                            colnames(factored) <- currlvl
                            rownames(factored) <- rownames(data)
                        } else { # matrix already exists, but need to add a column for the new level
                            zerocolumn <- rep(0,numrows)
                            factored <- cbind(factored, zerocolumn)
                            colnames(factored)[length(factored[1,])] <- currlvl
                        }
                        lvls <- c(lvls,currlvl) # add that level to the list
                    } # already found this level in another state.  Set the value to one
                    whichlvl <- which(lvls == currlvl) # this index number should correspond to the column number of the state
                    factored[row,whichlvl] <- 1
                    start <- ampLocs[part] + 1
                }
            }
        }
    }
    #Need to deal with any rows with missing data; fill in NA for all columns for that row
    for(missingrows in 1:length(missing)){
        for(column in 1:length(factored[1,])){
            factored[missing[missingrows],column] <- 1 # All states equally likely
        }
    }
    factored <- factored[,order(colnames(factored))]
    return(factored)
}


#This is a function for generating starting branch lengths -- based on Rogers and Swofford (1998), except here parsimony branch lengths are based on ACCTRAN as opposed to the MPR. Note that any sort of ambiguity is treated as complete ambiguity across all states:
RogersSwoffordish <- function(phy, data){
    pre.phyDat <- c()
    for(site.index in sequence(dim(data)[2])){
        tmp <- FactorDataLewis(data, charnum=site.index)
        tmp.vec <- apply(tmp, 1, which.max)
        tmp.vec[which(rowSums(tmp)>1)] = "?"
        pre.phyDat <- cbind(pre.phyDat, tmp.vec)
    }
    dat <- as.matrix(pre.phyDat)
    dat <- phyDat(dat,type="USER", levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"), ambiguity=c("?"))
    mpr.tre <- acctran(phy, dat)
    mpr.tre$edge.length <- mpr.tre$edge.length/dim(data)[2]
    return(mpr.tre)
}


# Based on a post by Liam Revell and Emmanual Paradis on R-sig-phylo
AddTaxonEverywhere <- function(tree,tip.name) {
    tree <- unroot(tree)
    tree$edge.length <- rep(1, nrow(tree$edge))
    new.tip <- compute.brlen(stree(1, tip.label = tip.name), 1)
    trees <- list()
    class(trees) <- "multiPhylo"
    for(i in 1:nrow(tree$edge)){
		trees[[i]] <- bind.tree(tree,new.tip, where=tree$edge[i,2], position=0.5)
		trees[[i]]$edge.length <- NULL
    }
    return(trees)
}


######################################################################################################################################
######################################################################################################################################
### Example run
######################################################################################################################################
######################################################################################################################################

# morph.data <- ReadNexusMorph("AsteralesPollen_MPfixed.nex")
# print(class(morph.data))
# phy <- read.nexus("Asterales.tre")
# tree.set <- AddTaxonEverywhere(phy, "Tubulifloridites_lilliei")
# pp <- lewisMkv(tree.set[[5]], morph.data, exclude.sites=c(17,18), ncores=4, include.gamma=TRUE)


######################################################################################################################################
######################################################################################################################################
### To do
######################################################################################################################################
######################################################################################################################################

#1. Starting values for the tree -- DONE. Uses a modified Rogers-Swofford approach.
#2. Tree search capabilities -- I forsee a TBR algorithm and a multistage optimization for the tree and model parameters
