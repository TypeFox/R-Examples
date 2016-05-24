`imbspace` <-
function (obj, data, depth = 1, L1.breaks = NULL, 
plot = TRUE, fixed = NULL, minimal = 1, maximal = 6,
M=250, raw.profile=NULL) 
{
    if (class(obj) != "cem.match") 
	stop("obj must be of class `cem.match'")
    if (!is.null(raw.profile) & class(raw.profile) != "L1profile") 
	stop("raw.profile must be of class `L1profile'")
    
    grouping <- obj$grouping
    if (!is.null(grouping) & !is.null(names(grouping))) {
        gn <- names(grouping)
        n.gn <- length(gn)
        for (g in 1:n.gn) {
            if (!is.null(data)) 
			data[[gn[g]]] <- group.var(data[[gn[g]]], grouping[[g]])
            if (!is.null(obj$breaks)) 
			obj$breaks[gn[g]] <- NULL
        }
    }

	if (is.null(L1.breaks)) {
		L1.breaks <- obj$imbalance$L1$breaks
		if (is.null(L1.breaks)) {
			vars <- colnames(data)
			nv <- length(vars)
			L1.breaks <- vector(nv, mode = "list")
			for (i in 1:nv) {
				L1.breaks[[i]] <- pretty(range(as.numeric(data[[i]]), na.rm = TRUE), 
											 n = nclass.scott(as.numeric(na.omit(data[[i]]))), 1)
				names(L1.breaks) <- vars
			}
		}
	}

    vnames <- obj$vars
    nv <- length(vnames)
    v.num <- 1:nv

    b.seq <- vector(nv, mode = "list")
    names(b.seq) <- vnames
    if (!is.null(fixed)) {
        idx <- match(fixed, vnames)
        if (length(idx) > 0) 
		v.num <- v.num[-idx]
    }
	tmp.min <- 2
	tmp.max <- 7
	if(!is.list(minimal)){
	 tmp.min <- minimal+1
	 minimal <- vector(nv, mode="list")
	 for (i in v.num) 
	  minimal[[i]] <- tmp.min
	}
	if(!is.list(maximal)){
		tmp.max <- maximal+1
		maximal <- vector(nv, mode="list")
		for (i in v.num) 
		maximal[[i]] <- tmp.max
	}
    for (i in v.num) {
		vna <- vnames[i]
		min.br <- tmp.min
		max.br <- tmp.max
		
		nuval <- length(unique(data[[vna]]))

        if( (nuval==2) | !(is.numeric(data[[vna]]) | is.integer(data[[vna]]) | is.logical(data[[vna]])) )
			max.br <- nuval+1
		
		if (!is.null(minimal[[vna]]))  
			min.br <- minimal[[vna]] + 1
		min.br <- max(tmp.min, min.br)
		if (!is.null(maximal[[vna]])) 
			max.br <- maximal[[vna]] + 1
		max.br <- min(tmp.max, max.br)
		b.seq[[i]] <- min.br:max.br
    }
    relax <- NULL
    for (d in 1:depth) {
        relax <- c(relax, combn(v.num, d, simplify = F))
    }
    n.relax <- length(relax)
    n.comb <- 0
    for (r in 1:n.relax) {
        v.idx <- relax[[r]]
        n.comb <- n.comb + dim(expand.grid(b.seq[v.idx]))[1]
    }
    max.k <- n.comb
	imb0 <- NULL
	medianL1 <- NULL
	medianCP <- NULL	

	if(is.null(raw.profile)){
	 cat("\nCalculating L1 profile for the raw data...\n")
	 imb0 <- L1.profile(obj$groups, data[, obj$vars], M=M, plot=FALSE)
     medianL1 <- median(imb0$L1)
	 medianCP <- imb0$CP[[ which(imb0$L1>medianL1)[1] ]]
	} else {
		imb0 <- raw.profile
		medianL1 <- raw.profile$medianL1
		medianCP <- raw.profile$medianCP
	}
	
    cat(sprintf("Executing %d different relaxations\n", n.comb))
    tab <- as.data.frame(matrix(NA, n.comb + 2, 2 * obj$n.groups + 
								2 + 1))
    colnames(tab) <- c(paste("G", obj$g.names, sep = ""), 
					   paste("PercG", obj$g.names, sep = ""), "Relaxed", "ML1", "var")
	n.coars <- dim(tab)[1]
	coars <- vector(n.coars, mode="list")
    tab[1, 1:obj$n.groups] <- obj$tab[2, ]
    tab[1, (obj$n.groups + 1):(2 * obj$n.groups)] <- obj$tab[2, ]/obj$tab[1, ] * 100
    tab$Relaxed <- "<start>"
    tab$var <- "<start>"
    IDX <- which(obj$matched)
	coars[[1]] <- obj$breaks

	tab[1, "ML1"] <- L1.meas(obj$groups, data[, obj$vars], breaks=medianCP, weights=obj$w)$L1

	tab[2, 1:obj$n.groups] <- obj$tab[1, ]
    tab[2, (obj$n.groups + 1):(2 * obj$n.groups)] <- obj$tab[1, ]/obj$tab[1, ] * 100
    tab$Relaxed[2] <- "<raw>"
    tab$var[2] <- "<raw>"
	raw.breaks <- obj$breaks
	for(i in names(raw.breaks))
		raw.breaks[[ i ]] <- range(data[[ i ]], na.rm=TRUE)
	
	coars[[2]] <- raw.breaks
	tab[2, "ML1"] <- medianL1

    
    last <- 0
    k <- 1
    K.tab <- 3
    adv <- 0
	pb <- txtProgressBar(min = 1, max = n.relax, initial = 1, style = 3)

    for (r in 1:n.relax) {
		setTxtProgressBar(pb, r)
        v.idx <- relax[[r]]
        brk <- expand.grid(b.seq[v.idx])
        r1 <- dim(brk)[1]
        c1 <- dim(brk)[2]
        newcut <- obj$breaks
        for (i in 1:r1) {
            X <- obj$X
            for (j in 1:c1) {
                vna <- colnames(brk)[j]
                newcut[[vna]] <- brk[i, j]
                tmp <- reduce.var(as.numeric(data[[vna]]), brk[i, j])
                X[vna] <- tmp$x
                newcut[[vna]] <- tmp$breaks
            }
            tmp.obj <- cem.match(data = X, verbose = 0)
			coars[[K.tab]] <- newcut
            tmp.obj$groups <- obj$groups
            tmp.obj$g.names <- obj$g.names
            tmp.obj$n.groups <- obj$n.groups
            tmp.obj$group.idx <- obj$group.idx
            tmp.obj$group.len <- obj$group.len
            mstrata <- find.strata(tmp.obj)$mstrata
            tmp.obj$mstrata <- mstrata
            tmp.obj$matched <- !is.na(mstrata)
			tmp.obj$n <- length(tmp.obj$matched)
            tmp.obj$treatment <- obj$treatment
			tmp.obj$tab <- cem.summary(tmp.obj)
			tmp.obj$w <- cem.weights(tmp.obj)
            IDX <- which(tmp.obj$matched)

			tab[K.tab, "ML1"] <- L1.meas(obj$groups, data[, obj$vars], breaks=medianCP, weights=tmp.obj$w)$L1
			

            tab[K.tab, "var"] <- vna
            tmp.tab <- cem.summary(obj = tmp.obj, verbose = 0)
            tab[K.tab, 1:(2 * obj$n.groups)] <- as.numeric(c(tmp.tab[2, 
															 ], tmp.tab[2, ]/tmp.tab[1, ] * 100))
            r.str <- NULL
            for (v in 1:c1) {
                r.str <- c(r.str, sprintf("%s(%d)", colnames(brk)[v], 
										  brk[i, v] - 1))
            }
            tab$Relaxed[K.tab] <- paste(r.str, collapse = ", ")
            k <- k + 1
            K.tab <- K.tab + 1

        }
    }
	close(pb)
    idx <- order(tab[, "ML1"])
    tab <- tab[idx, ]
    rownames(tab) <- 1:(dim(tab)[1])
    tab$var <- factor(tab$var)
	out <- list(space = tab)
	out$L1breaks <- L1.breaks
	out$raw.profile <- imb0
	out$tab <- obj$tab
    out$medianCP <- medianCP
	out$medianL1 <- medianL1
	out$coars <- coars[idx]
	out$n.coars <- n.coars
	out$match <- obj
    class(out) <- "imbalance.space"
	
    if (plot) 
	  plot(out,data=data,explore=interactive())
	
    return(invisible(out))
}


