.capwords <-
function (s, strict = FALSE) 
{
    cap <- function(s) paste(toupper(substring(s, 1, 1)), {
        s <- substring(s, 2)
        if (strict) 
            tolower(s)
        else s
    }, sep = "", collapse = " ")
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
.remove.space <-
function (s)
{
	sapply(strsplit(s, split=" ", fixed=TRUE), paste0, collapse="")
}
.check.input.consistency <-
function (cbias, div, warn = FALSE, family.sep) 
{
    err.FUN <- stop
    if (warn) 
        err.FUN <- warning
    seq.cbias <- if (which(colnames(cbias) == "Type") == 1 || !is.null(rownames(cbias))) 
		rownames(cbias) else cbias[, 1]
		
    seq.div <- unique(c(.get.TE.sub(as.character(div[, 1]), species = 1, family.sep=family.sep), 
						.get.TE.sub(as.character(div[, 1]), species = 2, family.sep=family.sep)))
    extra <- seq.cbias[!(seq.cbias %in% seq.div)]
    if (length(extra) > 0) 
        err.FUN("Sequences ", paste(extra, collapse = " "), " are in CB but not in Divergence.")
    extra <- seq.div[!(seq.div %in% seq.cbias)]
    if (length(extra) > 0) 
        err.FUN("Sequences ", paste(extra, collapse = " "), " are in Divergence but not in CB.")
    sp.cbias <- colnames(cbias)[-(1:which(colnames(cbias)=="Type"))]
    sp.div <- unique(c(as.character(div[, 3]), as.character(div[, 4])))
    extra <- sp.cbias[!(sp.cbias %in% sp.div)]
    if (length(extra) > 0) 
        err.FUN("Species ", paste(extra, collapse = " "), " are in CB but not in Divergence.")
    extra <- sp.div[!(sp.div %in% sp.cbias)]
    if (length(extra) > 0) 
        err.FUN("Species ", paste(extra, collapse = " "), " are in Divergence but not in CB.")
}
.check.species <-
function (vhica.obj, user.species = NULL, tree.species = NULL) 
{
    species <- user.species
    data.species <- unique(unlist(strsplit(names(vhica.obj$reg), 
        split = "X")))
    if (is.null(species)) {
        species <- data.species
    }
    if (!all(data.species %in% species)) {
        warning("More species in the data than specified in the species vector.")
        species <- unique(c(species, data.species))
    }
    if (!all(species %in% data.species)) {
        warning("Some of the species are not in the data. Dropped.")
        species <- species[!species %in% data.species]
    }
    if (!is.null(tree.species) && !all(species %in% tree.species)) {
        warning("Labels in the phylogenetic tree do not match the data set. Better not to plot the tree.")
    }
    if (!is.null(tree.species)) {
        species <- species[match(tree.species, species)]
    }
    if (is.null(names(species))) {
        names(species) <- species
    }
    names(species)[is.na(names(species))] <- species[is.na(names(species))]
    return(species)
}
.element.present <-
function (vhica.obj, element, species = NULL, skip.void = FALSE) 
{
    element.present <- NULL
    for (cross in names(vhica.obj$reg)) {
        sp12 <- unlist(strsplit(cross, "X"))
        target.table <- vhica.obj$reg[[cross]][[vhica.obj$target]]
        if (nrow(target.table) > 0) {
            for (index.target in 1:nrow(target.table)) {
                fullname <- rownames(target.table)[index.target]
                if (.get.TE.fam(fullname, family.sep=vhica.obj$family.sep) == element) {
                  if (.get.TE.sub(fullname, sub.only = TRUE, family.sep=vhica.obj$family.sep) == 
                    "") {
                    element.present <- c(element.present, sp12)
                  }
                  else {
                    if (!is.na(target.table[index.target, "resid"])) {
                      element.present <- c(element.present, paste(sp12, 
                        c(.get.TE.sub(fullname, species = 1, sub.only = TRUE, family.sep=vhica.obj$family.sep), 
						  .get.TE.sub(fullname, species = 2, sub.only = TRUE, family.sep=vhica.obj$family.sep)), 
						 sep = vhica.obj$family.sep))
                    }
                  }
                }
            }
        }
    }
    element.present <- unique(element.present)
    if (!skip.void && !is.null(species)) {
        already.there <- unique(sapply(element.present, function(el) {
            strsplit(el, vhica.obj$family.sep, fixed = TRUE)[[1]][1]
        }))
        element.present <- c(element.present, species[!(species %in% 
            already.there)])
    }
    if (!is.null(species)) {
        species.only <- sapply(element.present, function(el) {
            strsplit(el, vhica.obj$family.sep, fixed = TRUE)[[1]][1]
        })
        subspecies.only <- sapply(element.present, function(el) {
            strsplit(el, vhica.obj$family.sep, fixed = TRUE)[[1]][2]
        })
        element.present <- element.present[order(match(species.only, 
            species), subspecies.only)]
    }
    return(element.present)
}
.get.TE.fam <-
function (seqname, family.sep) 
{
    sapply(strsplit(seqname, split = family.sep, fixed=TRUE), function(s) s[1])
}
.get.TE.sub <-
function (seqname, species = 1, sub.only = FALSE, family.sep) 
{
    if (!species %in% c(1, 2)) {
        stop("Species should be 1 or 2")
    }
    .get.TE.sub.intra.subonly <- function(seq) {
        sp <- strsplit(seq, split = family.sep, fixed=TRUE)[[1]]
        if (length(sp) == 1) 
            return("")
        if (length(sp) == 2) 
            return(sp[2])
        if (length(sp) == 3) 
            return(sp[1 + species])
        stop(paste0("Error: sequence name ", seq, " not properly formatted"))
    }
    .get.TE.sub.intra <- function(seq) {
        sp <- strsplit(seq, split = family.sep, fixed=TRUE)[[1]]
        if (length(sp) == 1) 
            return(sp)
        if (length(sp) == 2) 
            return(paste(sp[1], sp[2], sep = family.sep))
        if (length(sp) == 3) 
            return(paste(sp[1], sp[1 + species], sep = family.sep))
        stop(paste0("Error: sequence name ", seq, " not properly formatted"))
    }
    FUN.sub <- if (sub.only) 
        .get.TE.sub.intra.subonly
    else .get.TE.sub.intra
    return(sapply(seqname, FUN.sub))
}
.get.TE.fam.longseq <-
function(seqname, species.sep, family.sep)
{
	sapply(seqname, function(s) {
		fs <- strsplit(s, split=family.sep, fixed=TRUE)[[1]]
		fulln <- strsplit(fs[1], split=species.sep, fixed=TRUE)[[1]][1]
		if (length(fs) == 1) 
			return(fulln)
		else
			return(paste(fulln, fs[length(fs)], sep=family.sep))
	})
}
.make.col.obj <-
function (n = 1000, max.col = "blue", min.col = "red", mid.col = "white", 
    range = c(-5, 5), threshold = c(-1, 1) * abs(log10(0.05)), 
    threshcol = 0.1, colsqueeze=1, extr = c(-1000, 1000)) 
{
    .make.half <- function(nn, thr, ran, ext) {
        thr <- abs(thr)
        ran <- abs(ran)
        ext <- abs(ext)
        nn.t <- max(round(threshcol * nn), 2)
        nn.r <- max(nn - nn.t - 2, 2)
        ans <- 0
        if (thr > 0) 
            ans <- c(ans, seq(0, thr, length.out = nn.t)[-1])
        if (ran > thr) 
            ans <- c(ans, seq(thr, ran, length.out = nn.r)[-1])
        s <- 1/colsqueeze
        fun <- function(x)  ifelse(x < thr, x*x*(s-1)/thr + x*(2-s), (x*x*(s-1) - x*(s*(thr+ran) - 2*thr) + ran*thr*(s-1))/(thr-ran))
        ans <- fun(ans)
        ans <- ans[ans >= 0 & ans <= ran]       
        ans <- c(ans, ext)			
        return(ans)
    }
    ans <- list()
    range <- sort(range)
    if (length(threshold) > 0) 
        threshold <- sort(threshold)
    extr <- sort(extr)
    if (range[1] < extr[1]) 
        extr[1] <- range[1] - 0.01
    if (range[2] > extr[2]) 
        extr[2] <- range[2] + 0.01
    if (length(range) != 2) 
        stop("range should be a vector of size 2")
    if (length(extr) != 2) 
        stop("extr should be a vector of size 2")
    if (length(threshold) > 2) 
        stop("threshold should be of size 0, 1, and 2")
    if (length(threshold) == 0) 
        threshold <- extr
    if (sum(threshold < 0) == 0) 
        threshold <- c(0, threshold)
    if (sum(threshold > 0) == 0) 
        threshold <- c(threshold, 0)
    if (n < 7) 
        stop("At least 7 colors are necessary")
    ans$breaks <- sort(unique(c(-.make.half(round(n/2), ext = extr[1], 
        ran = range[1], thr = threshold[1]), .make.half(n - round(n/2) + 
        1, ext = extr[2], ran = range[2], thr = threshold[2]))))
    ans$col <- c(colorRampPalette(c(min.col, mid.col))(sum(ans$breaks < 
        0) - 1), colorRampPalette(c(mid.col, max.col))(sum(ans$breaks >= 
        0)))
    return(ans)
}
.plot.caption <-
function (col.obj, main = "", p.adjust.method = "none", nslices = 1000, 
    thresh.lines = NA) 
{
    par(mar = c(4, 1, 2, 1))
    if (!requireNamespace("plotrix", quietly=TRUE)) {
        frame()
        return
    }
    compl <- ""
    if (p.adjust.method != "none") {
        compl <- paste(" (", .capwords(p.adjust.method), ")", 
            sep = "")
    }
    ticks <- pretty(floor(col.obj$breaks[2]):ceiling(col.obj$breaks[length(col.obj$breaks) - 
        1]))
    ticks <- ticks[ticks >= round(col.obj$breaks[2]) & ticks <= 
        round(col.obj$breaks[length(col.obj$breaks) - 1])]
    plot(NULL, type = "n", yaxt = "n", xlab = paste("p-value", 
        compl, sep = ""), ylab = "", xlim = range(ticks), ylim = c(0, 
        0.4), bty = "n", xaxt = "n")
    axis(1, at = ticks, labels = parse(text = sapply(ticks, function(t) if (t == 
        0) 
        "1"
    else paste0("10^-", abs(t)))), las = 2)
    .fun.getcol <- function(pp) {
        col.obj$col[which(col.obj$breaks > pp)[1] - 1]
    }
    plotrix::gradient.rect(min(ticks), 0, max(ticks), 0.2, nslices = nslices, 
        col = sapply(seq(min(ticks), max(ticks), length.out = nslices), 
            .fun.getcol))
    if (!is.na(thresh.lines[1])) {
        abline(v = thresh.lines, lty = 3)
    }
    title(main)
}
.plot.matrix <-
function (pmatrix, species, elements, zlim = range(pmatrix, na.rm = TRUE), 
    col.obj = .make.col.obj(n = 1000), na.col = "gray", grid.col = "darkgray", 
    ...) 
{
	if(sum(!is.na(pmatrix)) == 0)
		warning("Nothing to plot. Check the dataset and/or the element name.")
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    realx <- 0.5
    for (sp in species) {
        nn <- length(grep(pattern = sp, x = elements))
        realx <- c(realx, realx[length(realx)] + seq(0, 1, length.out = nn + 
            1)[-1])
    }
    ccol <- col.obj$col
    cbreaks <- col.obj$breaks
    if (!is.na(na.col)) {
        dummy.val <- min(cbreaks) - 1
        ccol <- c(na.col, ccol)
        cbreaks <- c(dummy.val, cbreaks)
        pmatrix[is.na(pmatrix)] <- dummy.val
    }
    image(x = realx, y = max(realx) - rev(realx) + 0.5, z = t(pmatrix[nrow(pmatrix):1, 
        ]), axes = FALSE, col = ccol, breaks = cbreaks, zlim = zlim, 
        ...)
    if (!is.na(grid.col)) {
        abline(h = seq(from = 0.5, by = 1, length.out = length(species) + 
            1), col = grid.col)
        abline(v = seq(from = 0.5, by = 1, length.out = length(species) + 
            1), col = grid.col)
    }
    axis(2, at = 1:length(species), labels = rev(names(species)), 
        las = 2, lwd.ticks = 0, lwd = 0, family = "mono")
    axis(3, at = 1:length(species), labels = names(species), 
        las = 2, lwd.ticks = 0, lwd = 0, family = "mono")
}
.plot.phylo <-
function (tree, species = "", horizontal = FALSE, show.tip.label = FALSE, 
    ...) 
{
    if (!requireNamespace("ape", quitely=TRUE)) 
        stop("Cannot plot trees without the package ape")
    shift <- if (length(tree$tip.label) < 15) 
        22/(length(tree$tip.label)^1.5)
    else 0
    if (horizontal) {
        par(mar = c(shift, 0.1, shift, 3.6))
        plot(ape::rotateConstr(tree, rev(species)), direction = "rightwards", 
            show.tip.label = show.tip.label, ...)
    }
    else {
        par(mar = c(3.6, shift, 0.1, shift))
        plot(tree, direction = "downwards", show.tip.label = show.tip.label, 
            ...)
    }
}
.plot.regression <-
function (reg, xlim = range(c(reg$model[, 2], reg[[length(reg)]][, 
    1]), na.rm=TRUE), ylim = range(c(reg$model[, 1]), reg[[length(reg)]][, 
    2], na.rm=TRUE), xlab = names(reg$model)[2], ylab = names(reg$model)[1], 
    reg.line = TRUE, elements = rownames(reg[[length(reg)]]), 
    pch.gene = 1, pch.element = 2, col.gene = "black", col.element = "black", 
    element.names = TRUE, lty.reg = 2, col.reg = "black", pval = NA, 
    lty.pval = 3, col.pval = "red", unilat = -1, ...) 
{
    plot(NULL, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
        ...)
    points(reg$model[, 2], reg$model[, 1], pch = pch.gene, col = col.gene)
    if (reg.line) {
        abline(reg, lty = lty.reg, col = col.reg)
    }
    if (!is.na(pval)[1]) {
        if (unilat == 0) {
            pval.lines <- qnorm(c(pval/2, 1 - pval/2))
        }
        else {
            pval.lines <- sign(unilat) * qnorm(1 - pval)
        }
        sdr <- sd(reg$residuals)
        for (pl in pval.lines) {
            abline(coef(reg)[1] + sdr * pl, coef(reg)[2], lty = lty.pval, 
                col = col.pval)
        }
    }
    points(reg[[length(reg)]][elements, 1], reg[[length(reg)]][elements, 
        2], pch = pch.element, col = col.element)
    in.elements <- elements %in% rownames(reg[[length(reg)]])

    if (element.names && sum(in.elements > 0)) {
        pos <- ifelse(reg[[length(reg)]][elements, 1] > mean(xlim), 
            2, 4)
    NAs <- is.na(pos)
        text(reg[[length(reg)]][elements, 1][!NAs], reg[[length(reg)]][elements, 
            2][!NAs], pos = pos[in.elements][!NAs], labels = elements[!NAs])
    }
}
.prepare.phylo <-
function (treefile) 
{
    if (is.null(treefile)) {
        warning("No tree file specfied: the phylogeny will not be plotted")
        return(NA)
    }
    if (!requireNamespace("ape", quietly=TRUE)) {
        warning("The ape library is not available: the phylogeny will not be plotted")
        return(NA)
    }
    try(tree <- ape::read.tree(treefile))
    if (inherits(tree, "try-error")) {
        return(NA)
    }
    return(tree)
}
.read.codon.bias <-
function (file, reference = "Gene") 
{
    rawdata <- read.table(file, header = TRUE, row.names = NULL)
    if (ncol(rawdata) < 3) {
        stop("Not enough columns in file ", file)
    }
    type.column <- which(sapply(1:ncol(rawdata), function(x) {
        is.factor(rawdata[, x])
    }))
    if (length(type.column) != 1) {
        stop("Only one column should be a factor in file ", file)
    }
    colnames(rawdata)[type.column] <- "Type"
    if (length(levels(rawdata$Type)) > 2) {
        stop("Only two type levels allowed in file ", file)
    }
    if (!(reference %in% levels(rawdata$Type))) {
        stop("No reference type ", reference, " in file ", file)
    }
    rownames(rawdata) <- rawdata[, 1]
    return(rawdata)
}
.read.divergence <-
function (file, divergence = "dS") 
{
    rawdata <- read.table(file, header = TRUE, row.names = NULL)
    if (ncol(rawdata) != 4) {
        stop("Number of columns different than 4 in file ", file)
    }
    rawdata[, 1] <- as.character(rawdata[, 1])
    return(rawdata)
}
.seq.codon.bias.clean <-
function(ll, gene.sep, species.sep, family.sep)
{
	ans <- list()
	names(ll) <- sapply(strsplit(basename(names(ll)), split=gene.sep, fixed=TRUE), 
		function(s) .remove.space(s[1]))
	for (i in 1:length(ll)) {
		namsp <- .get.TE.fam.longseq(names(ll[[i]]), species.sep=species.sep, family.sep=family.sep)
		subf <- .get.TE.sub(namsp, family.sep=family.sep, sub.only=TRUE)
		sp <- .get.TE.fam(namsp, family.sep=family.sep)
		if (all(subf == "")) {
			xx <- ll[[i]]
			names(xx) <- sp
			genename <- names(ll)[i]
			if (genename %in% names(ans)) 
				stop("Gene ", genename, " duplicated in the data set.")
			ans[[genename]] <- xx
		} else {
			for (ss in unique(subf)) {
				if (ss=="") stop("Inconsistent sub-family naming in gene ", names(ll)[i], ".")
				xx <- ll[[i]][subf==ss]
				names(xx) <- sp[subf==ss]
				genename <- paste(names(ll)[i], ss, sep=family.sep)
				if (genename %in% names(ans)) 
					stop("Gene ", genename, " duplicated in the data set.")				
				ans[[genename]] <- xx
			}
		}
	}
	return(ans)
}
.seq.codon.bias.format <-
function(CUB, name, tag, species.ref) 
{
	stopifnot(
		all(is.na(CUB)) || is.vector(CUB, mode="numeric"),
		length(CUB) > 0,
		!is.null(names(CUB)),
		nchar(name) > 0,
		nchar(tag) > 0,
		length(species.ref) > 0)
	if (any(nchar(names(CUB)) < 1))
		stop("In sequence file ", name, ": species names too short.")
	if (length(unique(names(CUB))) != length(names(CUB))) {
		stop("In sequence file ", name, ": Duplicated species names: ", names(CUB)[duplicated(names(CUB))], ".")
	}
	if (!all(names(CUB) %in% species.ref)) {
		stop("In sequence file ", name, ": Unknown species: ", paste(names(CUB)[!(names(CUB) %in% species.ref)], collapse=","), ".")
	}
	CUB <- CUB[species.ref]
	names(CUB) <- species.ref
	data.frame(Type=tag, as.list(CUB), row.names=name, check.names=FALSE)
}
.seq.codon.bias <-
function(gene.fasta, target.fasta, method="ENC", ref.name="Gene", targ.name="TE", species.sep="_", gene.sep=".", family.sep=".")
{
    if (requireNamespace("parallel", quietly=TRUE)) {
		mymclapply <- parallel::mclapply
	} else {
        mymclapply <- lapply
    }	
    CUBgenes <- mymclapply(gene.fasta, function(genefile) CUB(genefile, method=method))
    names(CUBgenes) <- gene.fasta
    CUBtarget <- mymclapply(target.fasta, function(genefile) CUB(genefile, method=method))
    names(CUBtarget) <- target.fasta
	listgenes <- .seq.codon.bias.clean(CUBgenes, gene.sep=gene.sep, species.sep=species.sep, family.sep=family.sep)
	listtarget <- .seq.codon.bias.clean(CUBtarget, gene.sep=gene.sep, species.sep=species.sep, family.sep=family.sep)
	
	species <- unique(unlist(lapply(c(listgenes, listtarget), names)))
	
	ans.genes <- do.call(rbind, mymclapply(1:length(listgenes), function(i) .seq.codon.bias.format(listgenes[[i]], names(listgenes)[i], tag=ref.name, species.ref=species)))
	ans.targ <- do.call(rbind, mymclapply(1:length(listtarget), function(i) .seq.codon.bias.format(listtarget[[i]], names(listtarget)[i], tag=targ.name, species.ref=species)))
	return(rbind(ans.genes, ans.targ))
}
.seq.divergence <-
function(sequence.fasta, divergence="dS", method="LWL85", pairwise=FALSE, species.sep="_", gene.sep=".", family.sep=".", max.lim=3)
{
    if (requireNamespace("parallel", quietly=TRUE)) {
		mymclapply <- parallel::mclapply
	} else {
        mymclapply <- lapply
    }	
	listseq <- mymclapply(sequence.fasta, function(genefile) {
		ans <- div(genefile, method=method, pairwise=pairwise, max.lim=max.lim)
		names(ans)[which(names(ans)=="div")] <- divergence
		seqn <- rep(.remove.space(strsplit(basename(genefile), split=gene.sep, fixed=TRUE)[[1]][1]), nrow(ans))
		
		fullsp1 <- .get.TE.fam.longseq(as.character(ans$sp1), species.sep=species.sep, family.sep=family.sep)
		fullsp2 <- .get.TE.fam.longseq(as.character(ans$sp2), species.sep=species.sep, family.sep=family.sep)
		ans$sp1 <- .get.TE.fam(fullsp1, family.sep=family.sep)
		ans$sp2 <- .get.TE.fam(fullsp2, family.sep=family.sep)
		subf1 <- .get.TE.sub(fullsp1, sub.only=TRUE, family.sep=family.sep)
		subf2 <- .get.TE.sub(fullsp2, sub.only=TRUE, family.sep=family.sep)
		if (any(xor(subf1 == "", subf2 == ""))) stop("Inconsistent species/sub-family naming")
		seqn <- ifelse(subf1=="", seqn, ifelse(subf1==subf2, paste(seqn, subf1, sep=family.sep), paste(seqn, subf1, subf2, sep=family.sep)))
		data.frame(seq=seqn, ans)
	})
	do.call(rbind, listseq)
}
.reference.regression <-
function (cbias, div, reference = "Gene", divergence = "dS", 
    CB.as.x = TRUE, warn = FALSE, family.sep=".") 
{
    full.list <- .tables2list(cbias, div, warn = warn, reference = reference, 
        divergence = divergence, family.sep=family.sep)
    if (requireNamespace("parallel", quietly=TRUE)) {
		mymclapply <- parallel::mclapply
	} else {
        mymclapply <- lapply
    }
    return(mymclapply(full.list, FUN = function(cross) {
        cross2 <- cross[cross$Type == reference, ]
        cross.TE <- cross[cross$Type != reference, ]
        meanCB <- 0.5 * cross2$CB1 + 0.5 * cross2$CB2
        div <- cross2$div
        names(div) <- names(meanCB) <- rownames(cross2)
        meanCB.TE <- 0.5 * cross.TE$CB1 + 0.5 * cross.TE$CB2
        div.TE <- cross.TE$div
        names(meanCB.TE) <- names(div.TE) <- rownames(cross.TE)
        ans <- NULL
        resid.TE <- NULL
        if (CB.as.x) {
            ans <- lm(div ~ meanCB)
            resid.TE <- div.TE - (meanCB.TE * coef(ans)[2] + 
                coef(ans)[1])
            ans[[levels(cross$Type)[2]]] <- data.frame(meanCB = meanCB.TE, 
                div = div.TE, resid = resid.TE, rel.res = resid.TE/sd(resid(ans)))
        } else {
            ans <- lm(meanCB ~ div)
            resid.TE <- meanCB.TE - (div.TE * coef(ans)[2] + 
                coef(ans)[1])
            ans[[levels(cross$Type)[2]]] <- data.frame(div = div.TE, 
                meanCB = meanCB.TE, resid = resid.TE, rel.res = resid.TE/sd(resid(ans)))
        }
        rownames(ans[[levels(cross$Type)[2]]]) <- rownames(cross.TE)
        return(ans)
    }))
}
.reverse.sub <-
function (seqname, family.sep) 
{
    sp <- strsplit(seqname, split = family.sep, fixed=TRUE)
    return(unlist(lapply(sp, function(ss) {
        if (length(ss) == 3) return(paste0(ss[1], family.sep, ss[3], 
            family.sep, ss[2]))
        if (length(ss) == 2) return(paste0(ss[1], family.sep, ss[2]))
        if (length(ss) == 1) return(ss[1])
        stop("Error: sequence name not properly formatted")
    })))
}
.stat.matrix <-
function (vhica.obj, element, elements, p.adjust.method = "none", 
    H1.test = "bilat") 
{
    ans <- matrix(NA, ncol = length(elements), nrow = length(elements))
    colnames(ans) <- rownames(ans) <- elements
    for (index.TE1 in 1:(length(elements) - 1)) {
        for (index.TE2 in (index.TE1 + 1):length(elements)) {
            TE1 <- elements[index.TE1]
            TE2 <- elements[index.TE2]
            decomp <- strsplit(c(TE1, TE2), split = vhica.obj$family.sep, fixed = TRUE)
            sp.TE1 <- decomp[[1]][1]
            sp.TE2 <- decomp[[2]][1]
            if (!paste(sp.TE1, "X", sp.TE2, sep = "") %in% names(vhica.obj$reg)) {
                tmp <- index.TE1
                new.index.TE1 <- index.TE2
                new.index.TE2 <- index.TE1
                TE1 <- elements[new.index.TE1]
                TE2 <- elements[new.index.TE2]
                decomp <- strsplit(c(TE1, TE2), split = vhica.obj$family.sep, fixed = TRUE)
                sp.TE1 <- decomp[[1]][1]
                sp.TE2 <- decomp[[2]][1]
            }
            else {
                new.index.TE1 <- index.TE1
                new.index.TE2 <- index.TE2
            }
            sub.TE1 <- if (length(decomp[[1]]) == 2) 
                decomp[[1]][2]
            else ""
            sub.TE2 <- if (length(decomp[[2]]) == 2) 
                decomp[[2]][2]
            else ""
            crossname <- paste(sp.TE1, sp.TE2, sep = "X")
            element.table <- vhica.obj$reg[[crossname]][[vhica.obj$target]]
            if (sub.TE1 == "") {
                if (sub.TE2 == "") {
                  linename <- element
                }
                else {
                  linename <- "DOESNOTEXIST"
                }
            }
            else {
                if (sub.TE1 == sub.TE2) {
                  linename <- paste(element, sub.TE1, sep = vhica.obj$family.sep)
                }
                else {
                  if (sub.TE2 == "") {
                    linename <- "DOESNOTEXIST"
                  }
                  else {
                    linename <- paste0(element, vhica.obj$family.sep, sub.TE1, 
                      vhica.obj$family.sep, sub.TE2)
                  }
                }
            }
            if (linename %in% rownames(vhica.obj$reg[[crossname]]$TE)) {
                norm.resid <- element.table[linename, "rel.res"]
                pval <- NA
                if (H1.test == "lower") {
                  p.val <- pnorm(norm.resid)
                }
                else if (H1.test == "bilat") {
                  p.val <- 2 * pnorm(-abs(norm.resid))
                  if (norm.resid > 0) 
				      p.val <- -p.val
                }
                else if (H1.test == "greater") {
                  p.val <- -pnorm(-norm.resid)
                }
                else {
                  stop("H1.test ", H1.test, " incorrect. Should be \"lower\", \"bilat\", or \"greater\".")
                }
                ans[TE1, TE2] <- ans[TE2, TE1] <- p.val
            }
        }
    }
    corrected.p <- log10(p.adjust(abs(ans[upper.tri(ans)]), method = p.adjust.method))
    corrected.p <- ifelse(ans[upper.tri(ans)] > 0, corrected.p, -corrected.p)
    ans[upper.tri(ans)] <- corrected.p
    ans[lower.tri(ans)] <- t(ans)[lower.tri(ans)]
    return(ans)
}
.tables2list <-
function (cbias, div, check = TRUE, keep.absent = FALSE, warn = FALSE, 
    reference = "Gene", divergence = "dS", family.sep) 
{
    if (requireNamespace("parallel", quietly=TRUE)) {
		mymclapply <- parallel::mclapply
	} else {
        mymclapply <- lapply
    }
    .make.unitary.table <- function(nn, sp1, sp2, sub.div) {
        if (nn %in% cbias[, 1]) {
            cc <- data.frame(Type = cbias[nn, "Type"], CB1 = cbias[nn, 
                sp1], CB2 = cbias[nn, sp2], div = sub.div[nn, 
                divergence], name = nn)
        }
        else {
            rev.species <- (sub.div[nn, "sp1"] == sp2)
            te1 <- .get.TE.sub(nn, species = 1, family.sep=family.sep)
            te2 <- .get.TE.sub(nn, species = 2, family.sep=family.sep)
            if (rev.species) {
                cc <- data.frame(Type = cbias[te1, "Type"], CB1 = cbias[te2, 
                  sp1], CB2 = cbias[te1, sp2], div = sub.div[nn, 
                  divergence], name = .reverse.sub(nn, family.sep=family.sep))
            }
            else {
                cc <- data.frame(Type = cbias[te1, "Type"], CB1 = cbias[te1, 
                  sp1], CB2 = cbias[te2, sp2], div = sub.div[nn, 
                  divergence], name = nn)
            }
        }
        return(cc)
    }
    if (check) 
        .check.input.consistency(cbias, div, warn = warn, family.sep=family.sep)
    ans <- list()
    species <- unique(c(colnames(cbias)[-c(1, 2)], as.character(div[, 
        3]), as.character(div[, 4])))
    for (index.sp1 in 1:(length(species) - 1)) {
        for (index.sp2 in (index.sp1 + 1):length(species)) {
            sp1 <- species[index.sp1]
            sp2 <- species[index.sp2]
            cross <- paste(sp1, sp2, sep = "X")
            sub.div <- div[(div[, 3] == sp1 & div[, 4] == sp2) | 
                (div[, 3] == sp2 & div[, 4] == sp1), ]
            compnames <- rownames(sub.div) <- as.character(sub.div[, 1])
            tt <- do.call(rbind, mymclapply(compnames, .make.unitary.table, 
                sp1 = sp1, sp2 = sp2, sub.div = sub.div))
            rownames(tt) <- tt$name
            tt$name <- NULL
            if (!keep.absent) {
                tt <- tt[!is.na(tt[, 4]), ]
            }
            ans[[cross]] <- tt
        }
    }
    return(ans)
}
.checkseq <-
function(seq, gene.name="") {
	# 0 check if the object makes sense
	stopifnot(
		length(seq) > 0,
		is.list(seq),
		all(sapply(seq, function(s) "SeqFastadna" %in% class(s))))
	# 1 check if all sequences have the same size
	ll <- sapply(seq, length)
	if (max(ll) != min(ll)) {
		warning(gene.name, " Sequences do not have the same length. Adding as many n as necessary.")
		seq <- lapply(seq, function(s) {ans <- c(s, rep("n", max(ll)-length(s))); class(ans) <- "SeqFastadna"; ans})
	}
	#2 check if sequence length is a multiple of 3
	mll <- max(ll)
	if (mll %% 3 != 0) {
		warning(gene.name, " Sequence length ", mll, " is not a multiple of 3. Truncating.")
		seq <- lapply(seq, function(s) {ans <- s[1:(3*(mll%/%3))]; class(ans) <- "SeqFastadna"; ans})
	}
	return(seq)
}
.ENC <-
function(seq, numcode=1, Wright.corr=TRUE) 
{
	stopifnot(
		"SeqFastadna" %in% class(seq),
		requireNamespace("seqinr"))
	yy <- seqinr::ucoweight(seq, numcode=numcode)
    yy.filt <- yy[sapply(yy,sum) > 1 & names(yy) != "*"]
    Fc <- sapply(yy.filt, function(x) {n <- sum(x); (n*sum((x/n)^2)-1)/(n-1)}) 
    SF <- sapply(yy.filt, length)
    SF <- SF[Fc != 0]
    Fc <- Fc[Fc != 0]
	ans <- 2 + sum(SF==2)/mean(Fc[SF==2]) + sum(SF==3)/mean(Fc[SF==3]) + sum(SF==4)/mean(Fc[SF==4]) + sum(SF==6)/mean(Fc[SF==6])
	ans[!is.finite(ans)] <- NA
	if (Wright.corr) {
		if (sum(SF==3)==0) # No Ile 
			ans <- ans + 1/(0.5*mean(Fc[SF==2])+0.5*mean(Fc[SF==4]))
		if (sum(SF==2) == 0 | sum(SF==4) == 0 | sum(SF==6) == 0) 
			ans <- NA
		if (is.finite(ans) && ans > 61) ans <- 61
	}	
	return(ans)
}
.LWL85 <-
function(seq, sq1=names(seq)[1], sq2=names(seq)[2], pairwise=TRUE, max.lim=max.lim)
{
	stopifnot(
		all(sapply(seq, function(s) "SeqFastadna" %in% class(s))),
		length(sq1) == length(sq2),
		all(c(sq1,sq2) %in% names(seq)),
		requireNamespace("seqinr"))
	if (!pairwise) {
		ali <- seqinr::as.alignment(nb=length(seq), nam=names(seq), seq=sapply(seq, function(s) paste(s, collapse="")))
		ks <- as.matrix(seqinr::kaks(ali)$ks)
		candidate <- ks[cbind(sq1, sq2)]
		return(ifelse (is.na(candidate) | candidate > max.lim, NA, candidate))
	} else {
		return(sapply(1:length(sq1), function(i) {
				subseq <- seq[c(sq1[i], sq2[i])]
				subali <- seqinr::as.alignment(nb=2, nam=c(sq1[i], sq2[i]), seq=sapply(subseq, function(s) paste(s, collapse="")))
				candidate <- seqinr::kaks(subali)$ks[1]
				return(if(is.na(candidate) || candidate > max.lim) NA else candidate)
			}))
	}
}
