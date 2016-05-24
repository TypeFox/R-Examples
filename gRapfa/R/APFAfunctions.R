
#library(igraph)
# get XY coordinates for an igraph graph using horizontal dot-like layout
getXY <- function(G) {
    XY  <- layout.sugiyama(G)$layout
    XY1 <- cbind(V(G)$level, -XY[, 1])
    dup <- duplicated(XY1)
    R   <- max(XY1[, 2]) - min(XY1[, 2])
    XY1[dup, 2] <- XY1[dup, 2] + R/5
    return(XY1)
}
# Sample tree function
st <- function(iD) {
     E <- data.frame(lapply(iD, addNA, ifany = TRUE))
    E <- data.frame(sapply(E, unclass))
    nr <- nrow(E)
    nc <- ncol(E)
    max.symb <- max(E)+1
    tablelist <- vector(mode="list", length=nc)
    from <- rep(0, nr)
    for (i in 1:nc) {
      key <- max.symb*from + E[,i]
      t <- table(key)
      from <- match(key, as.integer(names(t)))
      tablelist[[i]] <- t
    }
    no_edges_per_level <- sapply(tablelist, length)
    cs <- cumsum(no_edges_per_level)
    cs <- c(0, cs[-nc])
    ftsc <- matrix(NA, nrow=sum(no_edges_per_level), ncol=4)
    from <- 0
    for (i in 1:nc) {
      cnt <- tablelist[[i]]
      ain <- as.integer(names(cnt))
      symbol <- ain %% max.symb
      from <- max(from) + (ain - symbol)/max.symb
      to <- max(from) + 1:no_edges_per_level[i]
      ftsc[cs[i]+1:no_edges_per_level[i],]<-cbind(from,to,symbol,cnt)
    }
    ftsc <- data.frame(ftsc)
    names(ftsc) <- c("from","to","symbol","count")
    G <- graph.data.frame(ftsc)
    E(G)$symbol <- factor(E(G)$symbol)
    V(G)$size <- 2
    V(G)$level <- c(0, rep(1:nc, times=no_edges_per_level))
    cls <- c("red", "blue", "green", "black", "purple", "yellow",
        "orange", "magenta", colors()[c(73, 142, 91, 450, 547,
            100, 300, 400, 450, 500)])
    V(G)$count <- c(nr, E(G)$count)
    E(G)$color <- cls[as.numeric(unclass(E(G)$symbol))]
    E(G)$arrow.size <- 0.3
    G$nLevels <- sapply(iD, nlevels)
    G$fLevels <- lapply(iD, levels)
    G$merge.level <- 0
    G$N <- nr
    G$layout <- getXY(G)
    G$p <- length(G$nLevels)
    V(G)$label <- V(G)$name
    return(G)
}
# Adds edge probabilities, log-likelihood and dimension to an APFA object.
add.stats <- function(G) {
    if (is.null(G$dim)) {
        E <- get.edges(G, 1:ecount(G))
        E(G)$prob <- E(G)$count/V(G)[E[, 1]]$count
        d <- igraph::degree(G, mode = "out")
        G$dim <- sum(d[d > 0] - 1)  # of dubious validity: better to rely on degrees of freedom calculated by dIC
        G$loglike <- -sum(E(G)$count * log(E(G)$prob))
    }
    return(G)
}


# Calculates various quantities in connection with merging two nodes.
# Input: NS is a node by symbol array, the 1st half of the columns are target nodes, the 2nd half the edge counts.
#        When the corresponding edge is absent, the target node is set to 0.
#        mnode is a vector of nodes to be merged, specified as vertex ids (rather than names). Required to be of length two.

# Output: mmat is a kx2 matrix of vertex ids, containing the merge list
#         if get.map=TRUE, map is a vector of length vcount(G) giving the vertex ids of the vertices after merging
#         if test=TRUE, devtest contains the deviance and df associated with the merging
#         if doMerge=TRUE, the NS returned is the node by symbol array after merging (used in MergeNodes)

merge2nodes <- function(NS, mnode, test=TRUE, get.map=FALSE, doMerge =FALSE) {
  no_symbols <- ncol(NS)/2
  sym.indices <- 1:no_symbols
  cnt.indices <- (no_symbols+1):ncol(NS)
  if (test) doMerge <- TRUE
  cmat <- mnode
  dim(cmat) <- c(1,2)
  mmat <- cmat
  repeat {
    NS.1 <- NS[cmat[,1], sym.indices]
    NS.2 <- NS[cmat[,2], sym.indices]
    minC <- pmin(NS.1, NS.2)
    maxC <- pmax(NS.1, NS.2)
    wC <- (minC>0) & (minC != maxC)
    if (!any(wC)) break else {
       dd <- cbind(minC[wC], maxC[wC])
       mmat <- rbind(mmat, dd); cmat <- dd
    }
  }
  if (test) {
    all.merged <- c(mmat[,1], mmat[,2])
    rs <- rowSums(NS[all.merged,cnt.indices, drop=FALSE])
    LL0 <- sum(NS[all.merged,cnt.indices]*log(NS[all.merged,cnt.indices]/rs),na.rm=TRUE)
  }
  if (doMerge) {
    NS[mmat[,1], cnt.indices] <- NS[mmat[,1], cnt.indices] + NS[mmat[,2], cnt.indices]
    NS[mmat[,2], cnt.indices] <- 0
    A1 <- NS[mmat[,1], sym.indices]; A2 <- NS[mmat[,2], sym.indices]
    NS[mmat[,1], sym.indices] <- A1 + A2*(A1==0)
    NS[mmat[,2], sym.indices] <- 0
  }
  devtest <- c(NA,NA)
  if (test) {
    rs <- rowSums(NS[mmat[,1],cnt.indices, drop=FALSE])
    LL1 <- sum(NS[mmat[,1],cnt.indices]*log(NS[mmat[,1],cnt.indices]/rs),na.rm=TRUE)
    G2 <- 2*(LL0-LL1)
    A <- NS[mmat[,1],cnt.indices, drop=FALSE]
    df <- sum(A!=0)-nrow(A)
    devtest <- c(df, G2)
  }
  m <- NULL
  if (get.map) {
     prev.m <- map <- 1:nrow(NS); map[mmat[,2]] <- mmat[,1]
     repeat {m <- map[prev.m]; if (identical(m, prev.m)) break; prev.m <- m}
  }
  return(list(mmat=mmat, map=m, devtest=devtest, NS=NS))
}

# Merges two nodes (at the same level) in an APFA, returning the resulting APFA.
# nodeset is a vector of vertex names of length two.


MergeNodes <- function(G, nodeset, NS=NULL, setLayout=TRUE) {
	map=NULL
     if (length(V(G)[V(G)$level == G$p])>1) stop("not an APFA - need to contract last level")
     vidset <- match(nodeset, V(G)$name)
     levels <- V(G)[vidset]$level
     if (max(levels) != min(levels))
        stop("Error: attempt to merge nodes at different levels")
     if (levels[1] < G$merge.level)
        stop("Error: attempt to merge nodes at a lower level than a previous merge")

     if (is.null(NS)) NS <- apfa2NS(G)
     nm <- merge2nodes(NS, vidset, doMerge=TRUE, get.map=TRUE); NS <- nm$NS; map=nm$map

     no.symbols <- ncol(NS)/2
     sym.indices <- 1:no.symbols
     cnt.indices <- (no.symbols+1):ncol(NS)

     ind <- cbind(rep(1:nrow(NS), each=no.symbols), rep(1:no.symbols, times=nrow(NS)))
     w <- NS[ind]!=0
     edf <- data.frame(cbind(ind[w,1], NS[ind[w,]], ind[w,2], NS[,cnt.indices][ind[w,]]))
     names(edf) <- c("from", "to", "symbol", "count")

     edf[,1] <- map[edf[,1]]
     edf[,2] <- map[edf[,2]]

     cls <- c("red", "blue", "green", "black", "purple", "yellow", "orange", "magenta",
        colors()[c(73, 142, 91, 450, 547, 100, 300, 400, 450, 500)])

     edf$color <- cls[edf$symbol]
     edf$arrow.size <- 0.3
     G$merge.level <- levels[1]

     vdf <- data.frame(name = 1:vcount(G), level = V(G)$level, count = rowSums(NS[,cnt.indices]))

     G1  <- graph.data.frame(edf, vertices = vdf)
     V(G1)$name <- V(G)$name
     V(G1)$size <- V(G)$size
     G1$nLevels <- G$nLevels
     G1$fLevels <- G$fLevels
     G1$p <- G$p
     G1$merge.level <- G$merge.level
     G1$N <- G$N
     G1   <- G1 - V(G1)[igraph::degree(G1) == 0]
     V(G1)$label <- V(G1)$name
     if (setLayout) G1$layout <- getXY(G1)
     return(G1)
}

#derives a node by symbol array from an APFA, using matrix indexing. to=0 means no corresponding edge.
apfa2NS <- function(G) {
   if (length(V(G)[V(G)$level == G$p])>1) stop("not an APFA - need to contract last level")
   ftsc <-  data.frame(get.edges(G, 1:ecount(G)), E(G)$symbol, E(G)$count)
   names(ftsc) <- c("from", "to", "symbol", "count")
   no.symbols <- max(ftsc$symbol)
   NS <- array(0, dim=c(vcount(G), 2*no.symbols))
   NS[as.matrix(ftsc[,c(1,3)])] <- ftsc[,2]
   NS[,(1+no.symbols):(2*no.symbols)][as.matrix(ftsc[,c(1,3)])] <- ftsc[,4]
   return(NS)
}


MergeSelect <- function(G, NS=NULL, this.level, crit = "BIC", verbose = FALSE) {
    if (length(V(G)[V(G)$level == G$p])>1) stop("not an APFA - need to contract last level")
    if (!(this.level %in% 1:(G$p - 1)))
        stop("Error: level out of range")
    if (verbose)
        print(this.level)
    vnames <- V(G)[V(G)$level == this.level]$name
    if (verbose)
        print(vnames)
    k <- length(vnames)
    if (is.null(NS)) NS <- apfa2NS(G)
    if (k == 1) return(list(G=G,NS=NS))
    Res <- data.frame(cbind(rep(1:k, each = k), rep(1:k, times = k)))
    Res <- Res[Res[, 1] < Res[, 2], ]
    Res$score <- NA
    names(Res) <- c("i", "j", "score")
    repeat {
        for (i in which(is.na(Res$score))){
            Res$score[i] <- dIC(G, vnames[c(Res[i,1], Res[i, 2])], crit, NS)[1]
        }
        if (verbose)
            print(Res)
        wm <- which.min(Res$score)
        if (Res$score[wm] > 0)
            break
        ri <- Res$i[wm]
        rj <- Res$j[wm]
        G <- MergeNodes(G, nodeset=vnames[c(ri, rj)], NS, setLayout=FALSE)
        NS <- apfa2NS(G)
        if (!(vnames[ri] %in% V(G)$name)) {
            x <- ri
            ri <- rj
            rj <- x
        }
        V(G)$label <- ""
        Res$score[(Res$i == ri) | (Res$j == ri)] <- NA
        Res <- Res[(Res$i != rj) & (Res$j != rj), ]
        if (nrow(Res) == 0)
            break
    }
    return(list(G=G,NS=NS))
}

dIC <- function(G, nodeset, crit = "BIC", NS=NULL) {
   if (length(V(G)[V(G)$level == G$p])>1) stop("not an APFA - need to contract last level")
   vidset <- match(nodeset, V(G)$name)
   if (is.null(NS)) NS <- apfa2NS(G)
   dfdev <- merge2nodes(NS, vidset, test=TRUE)$devtest
   if (is.numeric(crit)) k <- crit else {if (crit == "BIC") k <- log(G$N) else k <- 2}
   delta.IC <- dfdev[2] - k * dfdev[1]
   return(c(dIC = delta.IC, dev = dfdev[2], df = dfdev[1]))
}

select.IC <- function(dat, crit = "BIC", verbose = FALSE) {
    G <- st(dat)
    G <- contract.last.level(G)
    NS <- apfa2NS(G)
    for (i in 1:(G$p - 1)){
        Gns <- MergeSelect(G, NS, i, crit, verbose)
        G <- Gns$G
        NS <- Gns$NS
    }
    G$layout <- getXY(G)
    return(G)
}

contract.last.level <- function(G) {
    no.levels <- G$p
    last.nodes <- V(G)[V(G)$level == no.levels]
    map <- 1:vcount(G)
    map[last.nodes] <- min(last.nodes)
    G1 <- contract.vertices(G, map, vertex.attr.comb = "first")
    G1 <- G1 - V(G1)[igraph::degree(G1) == 0]
    G1$layout <- getXY(G1)
    return(G1)
}

simulateAPFA <- function(g, Nsim = 1000) {
    if (is.null(E(g)$prob))
        g <- add.stats(g)
    EP <- data.frame(get.edges(g, 1:ecount(g)), E(g)$prob, E(g)$symbol)
    names(EP) <- c("from", "to", "prob", "symbol")
    EP$cum <- unsplit(lapply(split(EP, EP$from), function(x) cumsum(x$prob)), EP$from)
    nodedata <- matrix(nrow = Nsim, ncol = g$p + 1)
    symdata <- matrix(nrow = Nsim, ncol = g$p)
    nodedata[, 1] <- 1
    for (i in 1:g$p) {
        vset <- unique(nodedata[, i])
        for (v in vset) {
            num <- sum(nodedata[, i] == v)
            pr <- runif(num)
            cumP <- EP$cum[EP$from == v]
            to <- EP$to[EP$from == v]
            symb <- EP$symbol[EP$from == v]
            kron <- kronecker(pr, cumP, "<=")
            dim(kron) <- c(length(cumP), num)
            cS <- 1 + length(cumP) - colSums(kron)
            nodedata[nodedata[, i] == v, i + 1] <- to[cS]
            symdata[nodedata[, i] == v, i] <- symb[cS]
        }
    }
    symdata <- data.frame(symdata)
    names(symdata) <- paste("X", 1:g$p, sep = "")
    symdata <- data.frame(symdata)
    for (i in 1:ncol(symdata)) symdata[,i] <- factor(g$fLevels[[i]][symdata[,i]], levels=g$fLevels[[i]])
    names(symdata) <- names(g$fLevels)
    return(symdata)
}

select.beagle <- function(A, m=4, b=0.2, dir = '', row.marker = FALSE, col.hap = FALSE){
    A <- data.frame(lapply(A, addNA, ifany=TRUE))
    flevels <- lapply(A, levels)
    nlevels <- sapply(A, nlevels)
    A <- data.frame(lapply(A, unclass))  # added 30/10/13 to avoid problems with factor level labels with blanks
    if (!(row.marker & col.hap)) A <- t(A) 	
    dat <- as.data.frame(cbind('M',1:nrow(A), data.matrix(A))) # add needed columns to data
    write.table(dat, file = 'dat.bgl', col.names = F, row.names = F, quote = F)# save as bgl file
    system( paste('java -Xmx800m -jar ',dir,'beagle.jar data=',dir, 'dat.bgl scale=',m,' shift=',b,' out=D', sep='')) # run in BEAGLE
		
    RD <- scan("D.dat.bgl.dag.gz", skip = 1, # select required column to plot the graph.
      		what = list(Level = 0, Marker = "", Parent = 0, Child = 0, Allele = 0, Count = 0), flush = TRUE)
    Data <- data.frame(Level = RD$Level, Marker = RD$Marker, Parent = RD$Parent,
					Child = RD$Child, Allele = RD$Allele, Count = RD$Count)
    V.name <- function(dfn) {
           data.frame(from = paste(dfn$Level, '.', dfn$Parent, sep = ''),
	       	 to = paste(dfn$Level+1, '.', dfn$Child, sep = ''))
    }
    FT <- as.matrix(V.name(Data))
    #require(igraph)
    G <- graph.edgelist(FT)
    G$scale <- m
    G$shift <- b
    G$fLevels <- flevels
    G$nLevels <- nlevels

    E(G)$symbol<- factor(Data$Allele)

    FT <- data.frame(get.edges(G, 1:ecount(G)))
    names(FT) <- c("from","to")
    V(G)$level <- as.integer(V(G)$name)

    E(G)$label <- ""
    E(G)$count <- Data$Count
    FT <- data.frame(get.edges(G, 1:ecount(G)), E(G)$count)
    names(FT) <- c("from","to","count")
    a <- aggregate(FT$count, by=list(FT$from), sum)
    V(G)$count <- c(a$x, 0)
    E(G)$arrow.size <- 0.3
    V(G)$size  <- 2
    G$p <- max(V(G)$level)
    V(G)$label <- ""
    cls <- c("red", "blue", "green", "black", "purple", "yellow", "orange", "magenta",
        colors()[c(73, 142, 91, 450, 547, 100, 300, 400, 450, 500)])
    E(G)$color <- cls[as.numeric(unclass(E(G)$symbol))]
    G$layout <- getXY(G)
    G <- add.stats(G)
    return(G)
}
# ---------------------------------------------------------------------

       # MODEL COMPARISON

# ---------------------------------------------------------------------

# G is a fitted APFA and dat is a conformable dataset.
# returns the log-likelihood and the per-symbol log-likelihood.
# Uses the edge probabilities from G.
# see Thollard (2000)

LogLike.APFA <- function(G, dat, complete.cases=TRUE) {
    dat <- data.frame(lapply(dat, addNA, ifany=TRUE))
    subseteq <- function(a,b) length(setdiff(b,a))==0
    conformable <- all(mapply(subseteq, G$fLevels, lapply(dat, levels)))
    if (!conformable) stop("non-conformable")
    a <- data.frame(get.edges(G, 1:ecount(G)), E(G)$symbol)
    names(a) <- c("from", "to", "symbol")
    max.symb <- max(a$symbol, na.rm=TRUE)
    a$hash <- max.symb * (a$from - 1) + a$symbol
    from   <- 1
    E      <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat))
    for (i in 1:G$p){
        x <- unclass(dat[,i])
        key <- max.symb * (from - 1) + x
        e   <- match(key, a$hash)
        from<- a$to[e]
        E[ ,i]   <-  e
    }
    G <- add.stats(G)
    p.zero <- sum(!complete.cases(E)) # Number of obs that cannot be generated by G
    if (complete.cases) E <- E[complete.cases(E),]
    LL <- -sum(log(E(G)[E]$prob))  # the log-likelihood
    psLL <- -mean(log(E(G)[E]$prob))  # the per-symbol log-likelihood
    return(list(LL = LL, psLL = psLL, p.zero=p.zero))
}

# G is an APFA and D is a conformable dataset
# fits G to the dataset, ie the edge probabilities are calculated using D.
# Is this definition of conformability correct?

fit.APFA <- function(G, dat) {
  dat <- data.frame(lapply(dat, addNA, ifany=TRUE))
  conformable <- identical(G$fLevels, lapply(dat, levels))
  if (!conformable) stop("non-conformable")
  a <- data.frame(get.edges(G, 1:ecount(G)),E(G)$symbol)
  names(a)<- c('from', 'to', 'symbol')
  max.symb <- max(a$symbol)
  a$hash <- max.symb*(a$from-1) + a$symbol
  from <- 1
  E <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat))
  for (i in 1:G$p) {
    key <- max.symb*(from-1) + as.numeric(dat[,i])
    e <- match(key, a$hash)
    from <- a$to[e]
    E[ ,i] <- e
  }
  E <- E[complete.cases(E),]
  E(G)$count <- table(E)
  G$N <- sum(E(G)$count)
  FT <- data.frame(get.edges(G, 1:ecount(G)), E(G)$count)
  names(FT) <- c("from","to","count")
  a <- aggregate(FT$count, by=list(FT$from), sum)
  V(G)$count <- c(a$x, 0)
  G <- add.stats(G)
  return(G)
}

#Kullback-Leibler divergence A and B are two conformable APFA
KL <- function(A, B) {
    conformable <- identical(A$fLevels, B$fLevels)
    if (!conformable) stop("non-conformable")

    fl <- as.list(A$nLevels)
    for (i in 1:length(fl)) fl[[i]] <- 1:fl[[i]]
    Xa <- expand.grid(fl)  # full product space

    a <- data.frame(get.edges(A, 1:ecount(A)), E(A)$symbol)
    names(a) <- c("from", "to", "symbol")

    max.symb <- max(a$symbol)

    a$hash <- max.symb * (a$from - 1) + a$symbol
    from <- 1
    Ea <- Eb <- matrix(NA, nrow=nrow(Xa), ncol=ncol(Xa))
    for (i in 1:A$p) {
        key <- max.symb * (from - 1) + as.numeric(Xa[, i])
        e <- match(key, a$hash)
        from <- a$to[e]
        Ea[ ,i] <-  e
    }

    ok <- complete.cases(Ea)
    Xa <- Xa[ok, ]  # sample space X(A)
    Ea <- Ea[ok, ]  # corresponding root-to-sink paths in A

    b <- data.frame(get.edges(B, 1:ecount(B)), E(B)$symbol)
    names(b) <- c("from", "to", "symbol")

    max.symb <- max(b$symbol)

    b$hash   <- max.symb * (b$from - 1) + b$symbol
    from     <- 1
    
    for (i in 1:B$p) {
        key <- max.symb * (from - 1) + as.numeric(Xa[, i])
        e   <- match(key, b$hash)
        from<- b$to[e]
        Eb[ ,i]  <- e
    }

    ok <- complete.cases(Eb)
    if (any(!ok)) {
       # This is a kludge or fudge. Use the intersection = X(A) \cap X(B)
        Eb <- Eb[ok, ]
        Ea <- Ea[ok, ]
    }
    A <- add.stats(A)
    B <- add.stats(B)
    Alogprob <- log(E(A)[Ea]$prob)
    dim(Alogprob) <- dim(Ea)
    Blogprob <- log(E(B)[Eb]$prob)
    dim(Blogprob) <- dim(Eb)
    P <- exp(rowSums(Alogprob))
    Q <- exp(rowSums(Blogprob))
    P <- P/sum(P) # Q <- Q/sum(Q) removed
    KLdiv <- sum(P * log(P/Q))
    return(list(KLdiv=KLdiv, welldefined=all(ok)))
}

# crit <- integer values
# beagle = TRUE, if beagle model needs to be compared
# K-fold cross-validation
cross.validate <- function(Data, K=10, crit = NULL, beagle = TRUE, dir=''){
        N <- nrow(Data)
        logN <- round(log(N))
        if(is.null(crit)){
            crit <- c('AIC', 'BIC')
        }
        size <- N %/% K
        set.seed(15)
        rdm <- runif(N)
        ranked <- rank(rdm)
        block <- (ranked-1) %/% size+1
        block <- as.factor(block)

        if(beagle){
            mb <- cbind(rep(1:4, each=6), rep(seq(0,1, by=0.2), 4))
            psLL <- pzero <- matrix(0, nrow = K, ncol = length(crit)+nrow(mb))
        }else{
            psLL <- pzero <- matrix(0, nrow = K, ncol = length(crit))
        }
        for (k in 1:K) {
            train <- Data[block!=k,]
            test  <- Data[block==k,]
            A <- LL <- list()
            for(i in 1:length(crit)){
                A[[i]]      <- select.IC(train, crit = crit[i])
                LL[[i]]     <- LogLike.APFA(A[[i]], test)
                psLL[k, i]    <- LL[[i]]$psLL
                pzero[k, i] <- sum(LL[[i]]$p.zero)
            }
            if(beagle){
                for(j in 1:nrow(mb)){
                    A[[i+j]]      <- select.beagle(train, m = mb[j,1], b = mb[j,2], dir='')
                    LL[[i+j]]     <- LogLike.APFA(A[[i+j]], test)
                    psLL[k, i+j]    <- LL[[i+j]]$psLL
                    pzero[k, i+j] <- sum(LL[[i+j]]$p.zero)
                }
            }
            if(beagle){
                names(A) <- names(LL) <-
                colnames(psLL) <- colnames(pzero) <- c(paste('crit:',crit, sep=''),paste('m:',mb[,1],',b:',mb[,2],sep=''))
            }else{
                names(A) <- names(LL) <-
                colnames(psLL) <- colnames(pzero) <- paste('crit:',crit, sep='')
            }
            save(psLL, file='cvpsLL.rdata')
        }
        return(list(mean.psLL = colMeans(psLL), pzero = pzero))
}
