`relax.cem` <-
function (obj, data, depth = 1, verbose = 1, L1.breaks = NULL, 
plot = TRUE, fixed = NULL, shifts = NULL, minimal = NULL, 
use.coarsened = TRUE, eval.imbalance=TRUE, use.weights=FALSE, ...)
{
    if (class(obj) != "cem.match") 
	stop("obj must be of class `cem.match'")
    if(is.null(obj$X))
    stop("run cem with argument keep.all=TRUE")
    L1data <- NULL
    if(eval.imbalance)
	L1data <- data
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
    if(eval.imbalance){
		if (is.null(L1.breaks)) {
			L1.breaks <- obj$imbalance$L1$breaks
			if (is.null(L1.breaks)) {
				vars <- colnames(data)
				nv <- length(vars)
				L1.breaks <- vector(nv, mode = "list")
				for (i in 1:nv) {
                   if(!is.factor(data[[i]]))
					L1.breaks[[i]] <- pretty(range(data[[i]], na.rm = TRUE), 
											 n = nclass.scott(na.omit(data[[i]])), 1)
					
				}
                names(L1.breaks) <- vars
			}
		}
    }
    vnames <- obj$vars
    nv <- length(vnames)
    v.num <- 1:nv
    n.sh <- length(shifts)
    b.seq <- vector(nv, mode = "list")
    if (!is.null(fixed)) {
        idx <- match(fixed, vnames)
        if (length(idx) > 0) 
		v.num <- v.num[-idx]
    }
    for (i in v.num) {
        vna <- vnames[i]
        if (!is.null(obj$breaks[[vna]]) || (is.null(obj$breaks[[vna]]))) {
            if (!is.null(obj$breaks[[vna]])) {
                n.br <- length(obj$breaks[[vna]])
            }
            else {
                if (use.coarsened) 
				n.br <- nclass.FD(obj$X[[vna]])
                else n.br <- nclass.FD(as.numeric(data[[vna]]))
            }
            if (n.br <= 1) 
			n.br <- length(unique(obj$X[[vna]]))
            if (n.br > 1) {
                br.seq <- NULL
                min.br <- 2
                if (!is.null(minimal[[vna]])) 
				min.br <- minimal[[vna]] + 1
                min.br <- max(2, min.br)
                if (n.br > 10) {
					if ((min.br >= 2) & (min.br < 10)) 
                    br.seq <- c(seq(10, n.br - 1, min.br), 9:min.br)
                }
                else {
					br.seq <- n.br:min.br
                }
                br.seq <- sort(unique(br.seq))
                if (use.coarsened) 
				len.uni <- length(unique(obj$X[[vna]]))
                else len.uni <- length(unique(as.numeric(data[[vna]])))
                idx.max.br <- which(br.seq <= len.uni)
                b.seq[[i]] <- br.seq[idx.max.br]
            }
            else {
                b.seq[[i]] <- 1
            }
        }
    }
    names(b.seq) <- vnames
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
    if (n.sh > 0) {
        n.comb <- n.comb * 2 * n.sh
        max.k <- max.k * (1 + 2 * n.sh)
    }
    cat(sprintf("Executing %d different relaxations\n", n.comb))

	pb <- txtProgressBar(min = 1, max = n.relax, initial = 1, style = 3)

    tab <- as.data.frame(matrix(NA, n.comb + 1, 2 * obj$n.groups + 
								2 + 1))
    colnames(tab) <- c(paste("G", obj$g.names, sep = ""), paste("PercG", 
																obj$g.names, sep = ""), "Relaxed", "L1", "var")
    tab[1, 1:obj$n.groups] <- obj$tab[2, ]
    tab[1, (obj$n.groups + 1):(2 * obj$n.groups)] <- obj$tab[2, 
	]/obj$tab[1, ] * 100
    tab$Relaxed <- "<start>"
    tab$var <- "<start>"
    IDX <- which(obj$matched)
    if(eval.imbalance){
        if(use.weights){
            tab[1, "L1"] <- L1.meas(group=obj$groups[IDX],
                  data=L1data[IDX, obj$vars],
                  breaks=L1.breaks, weights=obj$w[IDX])$L1
        } else {
            tab[1, "L1"] <- L1.meas(group=obj$groups[IDX],
             data=L1data[IDX, obj$vars], breaks=L1.breaks)$L1
        }
    }
    last <- 0
    k <- 1
    K.tab <- 2
    adv <- 0
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
                if (use.coarsened) 
				tmp <- reduce.var(obj$X[[vna]], brk[i, j])
                else tmp <- reduce.var(as.numeric(data[[vna]]), 
									   brk[i, j])
                X[vna] <- tmp$x
                newcut[[vna]] <- tmp$breaks
            }
            tmp.obj <- cem.match(data = X, verbose = verbose)
            tmp.obj$groups <- obj$groups
            tmp.obj$g.names <- obj$g.names
            tmp.obj$n.groups <- obj$n.groups
            tmp.obj$group.idx <- obj$group.idx
            tmp.obj$group.len <- obj$group.len
            mstrata <- find.strata(tmp.obj)$mstrata
            tmp.obj$mstrata <- mstrata
            tmp.obj$matched <- !is.na(mstrata)
            if(use.weights){
                tmp.obj$n <- dim(data)[1]
                tmp.obj$treatment <- obj$treatment
                tmp.obj$baseline.group <- obj$baseline.group
                tmp.obj$tab <- cem.summary(obj=tmp.obj, verbose = verbose)
                tmp.w <- cem.weights(tmp.obj)
            }
            IDX <- which(tmp.obj$matched)
            if(eval.imbalance){
                if(use.weights){
                    tab[K.tab, "L1"] <- L1.meas(group=obj$groups[IDX],
                    data=L1data[IDX, tmp.obj$vars], breaks=L1.breaks,
                    weights=tmp.w[IDX])$L1
                } else {
                    tab[K.tab, "L1"] <- L1.meas(group=obj$groups[IDX],
                    data= L1data[IDX, tmp.obj$vars], breaks=L1.breaks)$L1
                }
            }
            tab[K.tab, "var"] <- vna
            tmp.tab <- cem.summary(obj = tmp.obj, verbose = verbose)
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
#            if (verbose == 1) {
#                iter <- as.integer(k/max.k * 10)
#                if (iter %in% c(2, 4, 6, 8, 10)) {
#					if (last != iter) 
#                    cat(sprintf("[%2d%%]", iter * 10))
#                }
#                else cat(".")
#                last <- iter
#            }
            if (n.sh > 0) {
                tmp.obj$breaks <- newcut
                tmp.obj$drop <- obj$drop
                tmp.obj$treatment <- obj$treatment
                tmp.obj$k2k <- FALSE
                class(tmp.obj) <- "cem.match"
                s.tmp.obj <- shift.cem(tmp.obj, data = data, 
									   shifts = shifts, verbose = verbose, plot = FALSE)
                IDX <- which(s.tmp.obj$matched)
                if(eval.imbalance)
				tab[K.tab, "L1"] <- L1.meas(obj$groups[IDX], 
											L1data[IDX, s.tmp.obj$vars], L1.breaks)$L1
                tmp.tab <- s.tmp.obj$tab
                tab[K.tab, 1:(2 * obj$n.groups)] <- as.numeric(c(tmp.tab[2, 
																 ], tmp.tab[2, ]/tmp.tab[1, ] * 100))
                tab$Relaxed[K.tab] <- sprintf("S:%s", r.str)
                K.tab <- K.tab + 1
#                if (verbose == 1) {
#					for (kk in k:(k + 2 * n.sh)) {
#iter <- as.integer(kk/max.k * 10)
#						if (iter %in% c(2, 4, 6, 8, 10)) {
#							if (last != iter) 
#							cat(sprintf("[%2d%%]", iter * 10))
#						}
#						else cat(".")
#						last <- iter
#					}
#                }
                k <- k + 2 * n.sh
            }
            if (verbose > 1) {
                cat(r.str)
                cat("\n")
            }
        }
    }
	close(pb)

    idx <- order(tab[, 1])
    tab <- tab[idx, ]
    rownames(tab) <- 1:(dim(tab)[1])
    tab$var <- factor(tab$var)
    out <- vector(obj$n.groups, mode = "list")
    names(out) <- paste("G", obj$g.names, sep = "")
    out[[1]] <- tab
    for (i in 2:obj$n.groups) out[[i]] <- tab[order(tab[, i]), 
	]
    class(out) <- "relax.cem"
    if(eval.imbalance)
	out$L1breaks <- L1.breaks
    if (plot) 
	plot(out, ...)
    return(invisible(out))
}
