
mstmap <- function(object, ...)
    UseMethod("mstmap")

mstmap.default <- function(object, ...)
      stop("Only \"data.frame\" and \"cross\" methods are available. See specific help files for more details.")

mstmap.cross <- function(object, chr, id = "Genotype", bychr = TRUE, suffix = "numeric",
                   anchor = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
                   p.value=1e-6, noMap.dist = 15.0, noMap.size = 0, miss.thresh = 1.0,
                   mvest.bc = FALSE, detectBadData = FALSE, return.imputed = FALSE, trace = FALSE, ...) {
    if(!(id %in% names(object$pheno)))
         stop("The unique identifier for the genotypes, ", deparse(substitute(id)),
              ", cannot be found in the object")
    if(!inherits(object, c("bc","dh","riself","bcsft")))
        stop("Cross object must inherit from one of the classes \"bc\", \"dh\",\"riself\",\"bcsft\". See ?mstmap.cross for more details.")
    if(!(dist.fun %in% c("haldane","kosambi")))
        stop("Distance function needs to be \"haldane\" or \"kosambi\" (see ?mstmap.cross).")
    if(!(objective.fun %in% c("COUNT","ML")))
        stop("Objective function needs to be \"COUNT\" or \"ML\" (see ?mstmap.cross).")
    if((miss.thresh > 1) | miss.thresh < 0)
        stop("Missing value threshold is required to be between 0 and 1.")
    if(class(object)[1] %in% c("bc","dh","riself")){
        pop.type <- "DH"
        allele <- c("A","B","")
    }
    if(class(object)[1] %in% "bcsft"){
        scheme <- attr(object, "scheme")
        if(scheme[1] > 0)
            stop("\"bcsft\" type is restricted to selfed populations only.")
        err <- lapply(object$geno, function(el){
            if(any(c(4,5) %in% el$data))
                stop("Only selfed populations with fully informative markers allowed. See ?mstmap.cross for more details.")})
        pop.type <- paste("RIL", scheme[2], sep = "")
        allele <- c("A","X","B")
    }
    if(trace) trace <- "MSToutput.txt"
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
        trace <- TRUE
    }
    oldObject <- NULL
    if(!missing(chr)) {
        oldObject <- subset(object, paste("-", chr, sep =""))
        object <- subset(object, chr)
    }
    geno <- as.character(object$pheno[[id]])
    lmat <- lapply(object$geno, function(el, geno, allele) {
                   rownames(el$data) <- geno
                   el$data[el$data == 1] <- allele[1]
                   el$data[el$data == 2] <- allele[2]
                   el$data[el$data == 3] <- allele[3]
                   el$data[is.na(el$data)] <- "U"
                   t(el$data)
               }, geno, allele)
    mn <- markernames(object)
    if(!bychr) {
        names(object$geno) <- ochr <- paste("ALL", 1:length(object$geno), sep = "")
        object$geno$"L"$data <- do.call("cbind", lapply(object$geno, function(el) el$data))
        object$geno$"L"$map <- unlist(lapply(object$geno, function(el) el$map))
        class(object$geno$"L") <- class(object$geno[[1]])
        lmat <- list(do.call("rbind", lmat))
        names(object$geno$"L"$map) <- rownames(lmat[[1]]) <- mn
        object <- subset(object, paste("-", ochr, sep = ""))
    }
    nm <- names(object$geno)
    param <- list("population_type"=pop.type,
                "distance_function"=dist.fun,
                "cut_off_p_value"=as.double(p.value),
                "no_map_dist"=as.double(noMap.dist),
                "no_map_size"=as.integer(noMap.size),
                "missing_threshold"=as.double(miss.thresh),
                "estimation_before_clustering"=as.integer(mvest.bc),
                "detect_bad_data"=as.integer(detectBadData),
                "objective_function"=objective.fun,
                trace=as.integer(trace))
    omit.list <- list()
    for(i in 1:length(lmat)){
        mst <- .Call("mst", param, as.data.frame(lmat[[i]], stringsAsFactors = FALSE))
        if (class(object)[1] == "riself") {
            imf <- switch(dist.fun, haldane = imf.h, kosambi = imf.k)
            mf <- switch(dist.fun, haldane = mf.h, kosambi = mf.k)
            mst <- lapply(mst, function(el, imf, mf){
                nm <- names(el$map)
                rf <- mf(diff(el$map))
                rf <- (rf/2)/(1 - rf)
                el$map <- cumsum(c(0, imf(rf)))
                names(el$map) <- nm
                el
            }, imf, mf)
        }
        chrw <- nm[i]
        nmm <- names(object$geno[[chrw]]$map)
        if(anchor)
            mst <- lapply(mst, function(el, nmm){
                  nmms <- nmm[nmm %in% names(el$map)]
                  lens <- 1:length(nmms)
                  apos <- pmatch(nmms, names(el$map))
                  rpos <- pmatch(nmms, names(rev(el$map)))
                  if(sum(abs(rpos - lens)) < sum(abs(apos - lens))) {
                      el$map <- rev(el$map)
                      el$map <- abs(el$map - el$map[1])
                  }
                  el }, nmm)
        map <- lapply(mst, function(el) el$map)
        ordm <- pmatch(names(unlist(map)), nmm)
        if(any(is.na(ordm)))
            omit.list[[i]] <- object$geno[[chrw]]$data[,is.na(ordm), drop = FALSE]
        ordm <- ordm[!is.na(ordm)]
        object$geno[[chrw]]$data <- object$geno[[chrw]]$data[,ordm, drop = FALSE]
        object$geno[[chrw]]$map <- unlist(map)
        if(length(map) > 1){
            spl <- sapply(map, function(el) names(el)[length(el)])
            spl <- list(spl[1:(length(spl) - 1)])
            names(spl) <- chrw
            cobject <- subset(object, chr = chrw)
            cobject <- breakCross(cobject, split = spl, suffix = suffix)
            object$geno <- object$geno[!(names(object$geno) %in% chrw)]
            chrw <- names(cobject$geno)
            object$geno <- c(object$geno, cobject$geno)
        }
        if(return.imputed){
            imp <- lapply(mst, function(el){
                names(el)[2] <- "data"
                el$data <- t(el$data)
                el$map <- el$map[names(el$map) %in% dimnames(el$data)[[2]]]
                el$data <- el$data[,pmatch(names(el$map), dimnames(el$data)[[2]]),drop = FALSE]
                el})
            names(imp) <- chrw
            if(!is.null(object$imputed.geno))
                object$imputed.geno <- object$imputed.geno[!(names(object$imputed.geno) %in% chrw)]
            object$imputed.geno <- c(object$imputed.geno, imp)
            object$imputed.geno <- object$imputed.geno[mixedorder(names(object$imputed.geno))]
        }
    }
    if(exists("omit.list"))
        object$omit <- do.call("cbind", omit.list[!sapply(omit.list, is.null)])
    object$geno <- c(object$geno, oldObject$geno)
    object$geno <- object$geno[mixedorder(names(object$geno))]
    object
}

breakCross <- function(cross, split = NULL, suffix = "numeric", sep = "."){
    if(is.null(split))
        stop("Split cannot be null.")
    chr <- names(split)
    if(!all(chr %in% names(nmar(cross))))
        stop("Some linkage group names do not exist in linkage map.")
    if(is.character(suffix)){
        if(!all(suffix %in% c("numeric","alpha")))
            stop("Character values for post names must be a vector of either \"numeric\" or \"alpha\".")
        if((length(suffix) == 1) & length(chr) > 1)
            suffix <- rep(suffix, length(chr))
        suffix <- as.list(suffix)
        names(suffix) <- chr
    } else if(is.list(suffix)){
       if(!all(names(suffix) %in% names(split)))
           stop("Some linkage group names in argument suffix do not match names in argument split.")
       if(any(is.na(pm <- pmatch(chr, names(suffix)))))
           warning("Some linkage group names missing from suffix argument .. using \"numeric\".")
           pnam <- split
           pnam[is.na(pm)] <- "numeric"
           pe <- pm[!is.na(pm)]
           pnam[!is.na(pm)] <- suffix[pe]
           suffix <- pnam
   } else stop("Post names are of the wrong type, see ?breakCross.")
   for(i in chr){
        cmap <- names(cross$geno[[i]]$map)
        if(any(is.na(mmatch <- pmatch(split[[i]], cmap))))
            stop("Some marker names do not exist in specified linkage groups.")
        pm <- c(mmatch, length(cmap))
        if(any(c("numeric","alpha") %in% suffix[[i]]))
            p.nam <- paste(i, switch(suffix[[i]], numeric = 1:length(pm),
                                     alpha = LETTERS[1:length(pm)]), sep = sep)
        else p.nam <- suffix[[i]]
        if(length(p.nam) != length(pm))
            stop("Length of linkage group post names does match number of splits.")
        slg <- rep(p.nam, times = c(pm[1], diff(pm)))
        lg <- lapply(p.nam, function(el, gen, slg) {
            res <- list()
            res$data <- gen$data[,slg %in% el, drop = FALSE]
            res$map <- gen$map[slg %in% el]
            res$map <- res$map - res$map[1]
            class(res) <- class(gen)
            res
        }, gen = cross$geno[[i]], slg)
        names(lg) <- p.nam
        cross$geno <- cross$geno[!(names(cross$geno) %in% i)]
        cross$geno <- c(cross$geno, lg)
    }
    cross$geno <- cross$geno[mixedorder(names(cross$geno))]
    cross
}

mergeCross <- function(cross, merge = NULL, gap = 5){
   if(is.null(merge))
       stop("Need list of linkage groups to merge.")
   if(any(sapply(merge, length) < 2))
       stop("Number of linkages groups to merge must be 2 or greater")
   if(any(!(unlist(merge) %in% names(nmar(cross)))))
       stop("Some listed linkage groups do not appear in cross object")
   for(i in 1:length(merge)){
      nm <- nmar(cross)
      if(is.null(lnam <- names(merge[i])))
          lnam <- paste(merge[[i]], sep = "")
      whl <- pmatch(merge[[i]], names(nm))
      sc <- cross$geno[whl]
      sc[[1]]$data <- do.call("cbind", lapply(sc, function(el) el$data))
      mapd <- sc[[1]]$map
      for(j in 2:length(sc))
          mapd <- c(mapd, sc[[j]]$map + max(mapd)[1] + gap)
      sc[[1]]$map <- mapd
      names(sc[[1]]$map) <- dimnames(sc[[1]]$data)[[2]]
      sc <- sc[1]
      names(sc)[1] <- lnam
      cross$geno <- cross$geno[-whl]
      cross$geno <- c(cross$geno, sc)
   }
   cross$geno <- cross$geno[mixedorder(names(cross$geno))]
   cross
}

subsetCross <- function(cross, chr, ind, ...){
    if(!inherits(cross, c("bc","dh","riself","bcsft")))
        stop("Cross object must inherit from one of the classes \"bc\", \"dh\",\"riself\",\"bcsft\". See ?subsetCross for more details.")
    if(!missing(chr))
        cross <- subset(cross, chr = chr, ...)
    if(!missing(ind)){
        n.ind <- nind(cross)
        if(!(is.numeric(ind) | is.logical(ind)))
            stop("Argument ind can only be logical or numeric.")
        cross <- subset(cross, ind = ind, ...)
        if(is.numeric(ind) & all(ind < 0))
            ind <- (1:n.ind)[ind]
        type <- c("co.located","seg.distortion","missing")
        if(any(wh <- type %in% names(cross))){
            type <- type[wh]
            for(i in type){
                cross[[i]]$data <- cross[[i]]$data[ind,,drop = FALSE]
                tcross <- cross
                if(i %in% c("seg.distortion","missing")){
                    tcross$geno[[i]]$data <- cross[[i]]$data
                    tcross$geno[[i]]$map <- 1:ncol(cross[[i]]$data)
                    class(tcross$geno[[i]]) <- "A"
                    tab <- geno.table(tcross, chr = i, scanone.output = TRUE)
                    if(class(cross)[1] == "bcsft") tab <- tab[,1:(ncol(tab) - 2)]
                    cross[[i]]$table[,4:ncol(cross[[i]]$table)] <- tab[,3:ncol(tab)]
                }
            }
        }
    }
    cross
}

mstmap.data.frame <- function(object, pop.type= "DH", dist.fun = "kosambi",
                   objective.fun= "COUNT", p.value=1e-6, noMap.dist=15.0, noMap.size=0,
                   miss.thresh=1.0, mvest.bc=FALSE, detectBadData=FALSE,
                   as.cross = TRUE, return.imputed = FALSE, trace=FALSE, ...) {
    if(trace) trace <- "MSToutput.txt"
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
        trace <- TRUE
    }
    if(!all(sapply(object, is.character)) & !all(sapply(object, is.numeric)))
        stop("Columns of the marker data.frame must be all of type \"character\" or all of \"numeric\".")
    RILN <- paste("RIL",1:20, sep = "")
    allow.pop <- c("BC","DH","ARIL",RILN)
    if(!(pop.type %in% allow.pop))
        stop("Population type needs to be \"BC\",\"DH\",\"ARIL\" or \"RILn\" (see ?mstmap.data.frame).")
    if(pop.type %in% RILN) rnum <- substring(pop.type, 4, nchar(pop.type))
    if(!(dist.fun %in% c("haldane","kosambi")))
        stop("Distance function needs to be \"haldane\" or \"kosambi\" (see ?mstmap.data.frame).")
    if(!(objective.fun %in% c("COUNT","ML")))
        stop("Objective function needs to be \"COUNT\" or \"ML\" (see ?mstmap.data.frame).")
    alleles <- unique(unlist(lapply(object, unique)))
    allow.list <- c("A","a","B","b","-","U")
    ptype <- pop.type
    if(pop.type %in% c("BC","DH","ARIL")){
        if(all(sapply(object, is.numeric))){
            if(any(apply(object, 2, function(el) length(el[is.na(el)]))))
                stop("Numeric input cannot contain missing values")
            if(any(apply(object, 2, function(el) el < 0 | el > 1)))
                stop("Ony values between 0 and 1 are allowed if input is numeric.")
            if(as.cross)
                stop("Cross object cannot be returned for numeric input.")
        } else {
            if(!all(alleles %in% allow.list))
                stop("Non-allowable allele encodings in DH marker set")
            cons <- c("1", "1", "2", "2", NA, NA, NA)
        }
        pop.type <- "DH"
    } else {
        if(is.numeric(object))
            stop("Numeric input is not available for RILn populations.")
        allow.list <- c(allow.list, "X")
        if(!all(alleles %in% allow.list))
            stop("Non-allowable allele encodings in RIL marker set")
        cons <- c("1", "1", "3", "3", NA, NA, "2")
    }
    if(length(grep(" ", rownames(object)))){
        rownames(object) <- gsub(" ", "-", rownames(object))
        warning("Replacing spaces in genotype names with a "-" separator\n")
    }
    if(length(grep(" ", names(object)))){
        rownames(object) <- gsub(" ", "-", names(object))
        warning("Replacing spaces in marker nxames with a "-" separator\n")
    }
    dist.fun <- match.arg(dist.fun)
    objective.fun <- match.arg(objective.fun)
    param <- list("population_type"=pop.type,
                "distance_function"=dist.fun,
                "cut_off_p_value"=as.double(p.value),
                "no_map_dist"=as.double(noMap.dist),
                "no_map_size"=as.integer(noMap.size),
                "missing_threshold"=as.double(miss.thresh),
                "estimation_before_clustering"=as.integer(mvest.bc),
                "detect_bad_data"=as.integer(detectBadData),
                "objective_function"=objective.fun,
                trace=as.integer(trace))
    mst <- .Call("mst", param, object)
    map <- lapply(mst, function(el) el$map)
    marks <- unlist(lapply(map, names))
    ordm <- pmatch(marks, rownames(object))
    if(any(is.na(ordm)))
        omit.list <- t(object[is.na(ordm), drop = FALSE])
    ordm <- ordm[!is.na(ordm)]
    object <- object[ordm,]
    chn <- paste("L", 1:length(map), sep = "")
    if (ptype == "ARIL") {
        imf <- switch(dist.fun, haldane = imf.h, kosambi = imf.k)
        mf <- switch(dist.fun, haldane = mf.h, kosambi = mf.k)
        map <- lapply(map, function(el, imf, mf){
            nams <- names(el)
            rf <- mf(diff(el))
            rf <- (rf/2)/(1 - rf)
            el <- cumsum(c(0, imf(rf)))
            names(el) <- nams
            el
        }, imf, mf)
    }
    if(as.cross){
        co <- list()
        spl <- rep(1:length(map), times = sapply(map, length))
        object <- as.matrix(object)
        al <- allow.list %in% alleles
        cons <- cons[al]
        al <- allow.list[al]
        for(i in 1:length(al))
            object[object == al[i]] <- cons[i]
        dn <- dimnames(object)[[1]]
        object <- apply(object, 2, function(el) as.numeric(el))
        dimnames(object)[[1]] <- dn
        datm <- lapply(split.data.frame(object, spl), t)
        mo <- lapply(1:length(map), function(el, datm, map){
            temp <- list()
            temp$data <- datm[[el]]
            temp$map <- map[[el]]
            class(temp) <- "A"
            temp
        }, datm, map)
        co$geno <- mo
        names(co$geno) <- chn
        if(return.imputed){
             co$imputed.geno <- lapply(mst, function(el){
                names(el)[2] <- "data"
                el$data <- t(el$data)
                el$map <- el$map[names(el$map) %in% dimnames(el$data)[[2]]]
                el$data <- el$data[,pmatch(names(el$map), dimnames(el$data)[[2]]),drop = FALSE]
                el})
             names(co$imputed.geno) <- chn
         }
        co$pheno <- data.frame(Genotype = factor(dimnames(object)[[2]]))
        wp <- (1:23)[allow.pop %in% ptype]
        class(co) <- c(c("bc","dh","riself",rep("f2",20))[wp],"cross")
        if(wp %in% 4:23) co <- convert2bcsft(co, F.gen = wp - 3, estimate.map = FALSE)
        object <- co
    } else {
        do <- list()
        chrv <- rep(chn, times = sapply(map, length))
        dist <- unlist(map)
        do$geno <- cbind.data.frame(markers = marks, chr = chrv, dist = dist, object)
        rownames(do$geno) <- NULL
        if(return.imputed){
            imp <- lapply(mst, function(el){
                nm <- names(el$map)[names(el$map) %in% dimnames(el$imputed_values)[[1]]]
                pm <- pmatch(nm, dimnames(el$imputed_values)[[1]])
                el$imputed_values <- el$imputed_values[pm, , drop = FALSE]
                el$imputed_values })
            chri <- rep(chn, times = sapply(imp, function(el) dim(el)[1]))
            marki <- unlist(lapply(imp, rownames))
            disti <- dist[pmatch(marki, names(dist))]
            do$imputed.geno <- cbind.data.frame(markers = marki, chr = chri, disti = disti, do.call("rbind", imp))
            rownames(do$imputed.geno) <- NULL
        }
        object <- do
    }
    if(exists("omit.list"))
        object$omit <- omit.list
    object
}

heatMap <- function (x, chr, mark, what = c("both", "lod", "rf"),
  lmax = 12, rmin = 0, markDiagonal = FALSE, color = rev(colorRampPalette(brewer.pal(11,"Spectral"))(256)), ...)
{
    opar <- par(no.readonly = TRUE)
    if (!inherits(x, "cross"))
        stop("Input should have class \"cross\".")
    what <- match.arg(what)
    if ("onlylod" %in% names(attributes(x$rf)) && attr(x$rf, "onlylod")) {
        onlylod <- TRUE
        what <- "lod"
    }
    else onlylod <- FALSE
    if (!missing(chr))
        x <- subset(x, chr = chr)
    n.mar <- nmar(x)
    if(!missing(mark)){
        if(!is.list(mark))
            mark <- list(mark)
        if(length(mark) == 1){
            if(length(n.mar) != 1)
                mark <- rep(mark, length(n.mar))
            names(mark) <- names(n.mar)
        }
        if(!all(names(mark) %in% names(n.mar)))
            stop("Names of marker subsets do not match names of linkage groups.")
        for(i in names(mark)){
            if(max(mark[[i]]) > length(x$geno[[i]]$map))
               stop("Range of marker subset greater than the number of markers in linkage group")
            x$geno[[i]]$map <- x$geno[[i]]$map[mark[[i]]]
            x$geno[[i]]$data <- x$geno[[i]]$data[,mark[[i]]]
        }
    }
    if (!("rf" %in% names(x))) {
        warning("Running est.rf.")
        x <- est.rf(x)
    }
    g <- x$rf
    old.xpd <- par("xpd")
    old.las <- par("las")
    par(las = 1)
    dots <- list(...)
    cols <- color
    if(is.null(dots$bigplot))
        dots$bigplot <- c(0.05, 0.85, 0.15, 0.9)
    if(is.null(dots$smallplot))
        dots$smallplot <- c(0.92, 0.94, 0.15, 0.9)
    if (what == "lod" && !onlylod) {
        g[lower.tri(g)] <- t(g)[lower.tri(g)]
        diag(g) <- lmax
        g[!is.na(g) & g > lmax] <- lmax
        lseq <- seq(0, lmax, length.out = 256)
        lmin <- min(g, na.rm = TRUE)
        cols <- cols[lmin - lseq < 0]
        iargs <- list(x = 1:nrow(g), y = 1:ncol(g), z = t(g), ylab = "Markers",
                      xlab = "Markers", col = cols, legend.args = list(side = 4,
                                          text = "LOD", las = 0, line = 2))
        do.call("image.plot", c(iargs, dots))
    }
    else if (what == "rf") {
        g[upper.tri(g)] <- t(g)[upper.tri(g)]
        diag(g) <- rmin
        g[!is.na(g) & g < rmin] <- rmin
        maxg <- max(g, na.rm = TRUE)
        rseq <- seq(rmin, 0.5, length.out = 256)
        cols <- rev(cols)[maxg - rseq > 0]
        clen <- 256*(maxg/0.5 - 1)
        clen <- ifelse(clen < 0, 0, clen)
        cols <- c(cols, colorRampPalette(c(cols[length(cols)], "white"))(clen))
        iargs <- list(x = 1:nrow(g), y = 1:ncol(g), z = t(g), ylab = "Markers",
                      xlab = "", col = cols, legend.args = list(side = 4,
                      text = "Recombination Fractions", las = 0, line = 2, crt = 180))
        do.call("image.plot", c(iargs, dots))
    }
    else {
        dots$bigplot[1] <- 0.15
        diag(g) <- lmax
        g[!is.na(g) & g > lmax] <- lmax
        glt <- g[lower.tri(g)]
        g[lower.tri(g)] <- NA
        lseq <- seq(0, lmax, length.out = 256)
        lmin <- min(g, na.rm = TRUE)
        coll <- cols[lmin - lseq < 0]
        iargs <- list(x = 1:nrow(g), y = 1:ncol(g), z = t(g), ylab = "",
                      xlab = "", col = coll, legend.args = list(side = 4,
                                                    text = "Linkage", las = 0, line = 2))
        do.call("image.plot", c(iargs, dots))
        g[lower.tri(g)] <- glt
        g[upper.tri(g)] <- NA
        diag(g) <- NA
        g[!is.na(g) & g < rmin] <- rmin
        dots$smallplot[1:2] <- c(0.05, 0.07)
        maxg <- max(g, na.rm = TRUE)
        rseq <- seq(rmin, 0.5, length.out = 256)
        colr <- rev(cols)[maxg - rseq > 0]
        clen <- 256*(maxg/0.5 - 1)
        clen <- ifelse(clen < 0, 0, clen)
        colr <- c(colr, colorRampPalette(c(colr[length(colr)], "white"))(clen))
        iargs <- list(x = 1:nrow(g), y = 1:ncol(g), z = t(g), ylab = "",
                      xlab = "Markers", col = colr, add = TRUE, legend.args = list(side = 2,
                      text = "Recombination", las = 0, line = 1.5))
        do.call("image.plot", c(iargs, dots))
 }
    par(plt = dots$bigplot)
    if (markDiagonal) {
        for (i in 1:ncol(g)) segments(i + c(-0.5, -0.5, -0.5,
            +0.5), i + c(-0.5, +0.5, -0.5, -0.5), i + c(-0.5,
            +0.5, +0.5, +0.5), i + c(+0.5, +0.5, -0.5, +0.5))
    }
    n.mar <- nmar(x)
    n.chr <- nchr(x)
    a <- c(0.5, cumsum(n.mar) + 0.5)
    abline(v = a, xpd = TRUE, col = "white", lwd = 0)
    abline(h = a, xpd = FALSE, col = "white", lwd = 0)
    a <- par("usr")
    wh <- cumsum(c(0.5, n.mar))
    chrnam <- names(x$geno)
    chrpos <- (wh[-1] + wh[-length(wh)])/2
    par(plt = dots$bigplot)
    c.args <- c(list(side = 3, at = chrpos, labels = chrnam,
                     tick = FALSE, line = -0.6), dots$axis.args)
    do.call("axis", c.args)
    c.args$side <- 4
    do.call("axis", c.args)
    par(opar)
    if ("main" %in% names(dots))
        title(main = dots$main, line = 2.5)
    else {
        if (what == "lod")
            title(main = "Pairwise LOD scores", line = 2.5)
        else if (what == "rf")
            title(main = "Recombination fractions", line = 2.5)
        else title("Pairwise recombination fractions and LOD scores", line = 2.5)
    }
    invisible()
}

quickEst <- function(object, chr, map.function = "kosambi", ...){
    if (!any(class(object) == "cross"))
        stop("Input should have class \"cross\".")
    if (missing(chr))
        chr <- names(nmar(object))
    imf <- switch(map.function, kosambi = imf.k, haldane = imf.h,
                  morgan = imf.m, cf = imf.cf)
    nm <- nmar(object)
    for(i in chr){
        temp <- subset(object, chr = i)
        est <- est.rf(temp)$rf
        nc <- dim(est)[1]
        er <- est[cbind(2:nc,1:(nc - 1))]
        temp$geno[[i]]$map <- c(0,cumsum(imf(er)))
        names(temp$geno[[i]]$map) <- dimnames(temp$geno[[i]]$data)[[2]]
        tempa <- argmax.geno(temp, step = 0, map.function = map.function, ...)
        tempa$geno[[i]]$data <- tempa$geno[[i]]$argmax
        tempa$geno[[i]] <- tempa$geno[[i]][-3]
        esta <- est.rf(tempa)$rf
        era <- esta[cbind(2:nc,1:(nc - 1))]
        if(class(object)[1] == "riself")
            era <- (era/2)/(1 - era)
        object$geno[[i]]$map <- c(0,cumsum(imf(era)))
        names(object$geno[[i]]$map) <- dimnames(object$geno[[i]]$data)[[2]]
    }
    object
}

genClones <- function(object, chr, tol = 0.9, id = "Genotype"){
    pairsfun <- function(mat1 = mat1, mat2 = mat2){
        ng <- c(unique(mat1),unique(mat2))
        ng2 <- ng[!is.na(ng)]^2
        mat1[is.na(mat1)] <- mat2[is.na(mat2)] <- pi
        apply(mat1*mat2, 1, function(el, ng2) {
             is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
             ma <- sum(rep(1, length(el))[el %in% ng2])
             nm <- sum(rep(1, length(el))[!(el %in% ng2) & is.wholenumber(el)])
             na <- sum(rep(1, length(el))[el == (pi^2)])
             c(ma, nm, na)
       }, ng2 = ng2)
    }
    if(!inherits(object, "cross"))
        stop(deparse(substitute(object)), " must be of class \"cross\"")
    if(!missing(chr))
        object <- subset(object, chr)
    if(!(id %in% names(object$pheno)))
         stop("The unique identifier for the genotypes, ", deparse(substitute(id)),
              ", cannot be found in the object")
    gmat <- pull.geno(object)
    gnam <- as.character(object$pheno[[id]])
    nm <- nmar(object)
    cgm <- cg <- comparegeno(object)
    cg[upper.tri(cg, diag = TRUE)] <- NA
    wh <- which(cg > tol, arr.ind = TRUE)
    if(dim(wh)[1] == 0){
        message("There are no genotype pairs with matching allele proportions greater than ", tol,".")
        return(list(cgm = cgm))
    }
    ind1 <- wh[,1]; ind2 <- wh[,2]
    cgd <- cbind.data.frame(G1 = gnam[ind1], G2 = gnam[ind2])
    cgd$coef <- round(cg[wh], 4)
    cgd <- cbind.data.frame(cgd, t(pairsfun(gmat[ind1,, drop = FALSE],gmat[ind2,, drop = FALSE])))
    names(cgd)[4:6] <- c("match", "diff","na.both")
    cgd$na.one <- ncol(gmat) - apply(cgd[,4:6], 1, sum)
    group <- rep(1, nrow(cgd))
    ## organize groups
    if(nrow(cgd) > 1){
        for(el in 2:nrow(cgd)){
            ind <- rep(1:(el - 1), 2)
            if(any(wg <- as.character(unlist(cgd[1:(el - 1),1:2])) %in% as.character(unlist(cgd[el,1:2])))){
                wg <- ind[wg]
                if(length(ug <- unique(group[wg])) > 1){
                    mug <- min(ug)
                    oth <- ug[ug != mug]
                    group[group %in% oth] <- mug
                    ug <- mug
                }
                group[el] <- ug
            }
            else group[el] <- max(group) + 1
        }
    }
    cgd$group <- group
    cgd <- cgd[order(cgd$group),]
    rownames(cgd) <- NULL
    dimnames(cgm) <- list(gnam, gnam)
    list(cgm = cgm, cgd = cgd)
}

fixClones <- function(object, gc, id = "Genotype", consensus = TRUE){
    if(!inherits(object, "cross"))
        stop(deparse(substitute(object)), " must be of class \"cross\"")
    if(missing(gc))
        stop("Argument \"gc\" cannot be missing.")
    if(!is.data.frame(gc))
        stop("Argument \"gc\" needs to be a data frame.")
    if(!all(c("G1","G2","group") %in% names(gc)))
        stop("Argument \"gc\" must have names \"G1\", \"G2\" and \"group\", see ?fixClones.")
    if(!(id %in% names(object$pheno)))
         stop("The unique identifier for the genotypes, ", deparse(substitute(id)),
              ", cannot be found in the object")
    ph <- as.character(object$pheno[[id]])
    nm <- nmar(object)
    sl <- split(gc, gc$group)
    bm <- pull.geno(object)
    rownames(bm) <- ph
    nd <- unlist(pull.map(object))
    nc <- rep(names(nm), times = nm)
    cl <- sapply(object$geno, function(el) class(el))
    colnames(bm) <- paste(nc, markernames(object), nd, cl, sep = ";")
    gomit <- list()
    for(i in 1:length(sl)){
        gn <- unique(as.character(unlist(sl[[i]][,c("G1","G2")])))
        gn <- gn[mixedorder(gn)]
        if(!all(gn %in% ph))
            stop("Some genotypes in clone set cannot be found in cross object.")
        wh <- (1:length(ph))[ph %in% gn]
        if(consensus){
            td <- apply(bm[wh, ], 2, function(mm){
                mm <- mm[!is.na(mm)]
                if(length(mmu <- unique(mm)) > 1 | !length(mmu))
                    NA else mmu
            })
            kp <- rownames(bm) %in% gn[1]
            bm[kp,] <- td
            gomit[[i]] <- gn[2:length(gn)]
        }
        else {
            len <- apply(bm[wh,], 1, function(el) length(el[!is.na(el)]))
            kp <- gn[len == max(len)][1]
            gomit[[i]] <- gn[!(gn %in% kp)]
        }
        rownames(bm)[kp] <- paste(gn, collapse = "_")
    }
    gomit <- unlist(gomit)
    object$pheno[[id]] <- factor(rownames(bm))
    object$geno <- lapply(split.data.frame(t(bm), nc), function(el){
        temp <- list()
        spl.c <- strsplit(rownames(el), ";")
        temp$data <- as.matrix(t(el))
        temp$map <- as.numeric(sapply(spl.c, "[", 3))
        names(temp$map) <- dimnames(temp$data)[[2]] <- sapply(spl.c, "[", 2)
        class(temp) <- unique(sapply(spl.c, "[", 4))
        temp
    })
    object <- subset(object, ind = !(object$pheno[[id]] %in% gomit))
    names(object$geno) <- names(nm)
    object
}

pp.init <- function(seg.thresh = 0.05, seg.ratio = NULL, miss.thresh = 0.1, max.rf = 0.25, min.lod = 3){
       if(!is.null(seg.thresh)){
           if(!(is.numeric(seg.thresh) | is.character(seg.thresh)))
               stop("seg.thresh argument of wrong type, see ?pullCross pr ?pushCross.")
           if(is.character(seg.thresh)){
               if(seg.thresh != "bonf")
                   stop("seg.thresh argument can only be numeric or equal to \"bonf\", see ?pullCross or ?pushCross")
           }
           if(is.numeric(seg.thresh)){
               if(seg.thresh >= 1 | seg.thresh < 0 )
                   stop("seg.thresh cannot exceed 1 and must be between 0 and 1")
           }
       }
       if(!is.null(seg.ratio)){
           seg.thresh <- NULL
           if(!is.character(seg.ratio))
               stop("seg.ratio argument of wrong type, see ?pullCross or ?pushCross.")
           if(!length(grep(":", seg.ratio)))
               stop("Alelle ratios must be split by \":\".")
       }
       if(!is.numeric(miss.thresh))
           stop("miss.thresh argument of wrong type, see ?pullCross or ?pushCross.")
       if(miss.thresh > 1 | miss.thresh < 0)
           stop("seg.thresh cannot exceed 1 and must be between 0 and 1")
       if(max.rf > 0.5)
           stop("Maximum recombination fraction cannot exceeed 0.5.")
       if(min.lod <= 0)
           stop("Maximum LOD score cannot be less than or equal to zero.")
       list(seg.thresh = seg.thresh, seg.ratio = seg.ratio, miss.thresh = miss.thresh, max.rf = max.rf, min.lod = min.lod)
   }

pullCross <- function(object, chr, type = c("co.located","seg.distortion","missing"), pars = NULL, replace = FALSE, ...){
    if(!is.null(object[[type]])){
        if(dim(object$pheno)[1] != dim(object[[type]]$data)[1])
            stop("Number of genotypes in linkage map does not match external ", type, " data.")
    }
    if(is.null(pars))
        pars <- pp.init()
    if(!is.list(pars))
        stop("argument pars must be a list object one or more named elements matching the results of pp.init()")
    pars <- do.call("pp.init", pars)
    type <- match.arg(type)
    if(!(class(object)[1] %in% c("bc","dh","riself","bcsft")))
        stop("Pulling of markers is not supported for this population type, see ?pullCross.")
    oldObject <- NULL
    if (!missing(chr)) {
        oldObject <- subset(object, paste("-", chr, sep =""))
        object <- subset(object, chr = chr)
    }
    nm <- nmar(object)
    chr <- names(nm)
    chrn <- rep(chr, times = nm)
    geno <- genot <- do.call("cbind", lapply(object$geno, function(el) el$data))
    rownames(genot) <- as.character(object$pheno[[1]])
    markn <- colnames(geno)
    whm <- FALSE
    if(type == "co.located"){
        dm <- findDupMarkers(object, chr, exact.only = FALSE)
        if(!is.null(dm)){
            cm <- unlist(dm)
            whm <- markn %in% cm
            km <- names(dm)
            markl <- chrl <- list()
            for(i in 1:length(dm)){
                markl[[i]] <- c(km[i], dm[[i]])
                chrl[[i]] <- chrn[markn %in% markl[[i]]]
            }
            lens <- sapply(dm, function(x) length(x) + 1)
            bins <- rep(1:length(lens), times = lens)
            tab <- cbind.data.frame(bins = bins, chr = unlist(chrl), mark = unlist(markl))
        }
    }
    if(type %in% c("seg.distortion","missing")){
        tab <- geno.table(object, scanone.output = TRUE)
        nc <- ncol(tab)
        if(class(object)[1] == "bcsft")
            nc <- nc - 2
        if(type == "seg.distortion"){
            if(!is.null(seg.thresh <- pars$seg.thresh)){
                pval <- 10^(-tab$neglog10P)
                if(is.character(seg.thresh)){
                    if(!(seg.thresh %in% "bonf"))
                        stop("Thresholding only exists for \"bonf\", see ?pull.cross")
                        seg.thresh <- 0.05/totmar(object)
                    }
                else if(seg.thresh > 1)
                    stop("Numerical threshold cannot exceed a p.value of 1.")
                whm <- pval < seg.thresh
            }
            else if(!is.null(seg.ratio <- pars$seg.ratio)) {
                ctab <- tab[,5:nc]
                sr <- as.numeric(unlist(strsplit(seg.ratio, ":")))
                if(length(sr) != ncol(ctab))
                    stop("Wrong number of ratio elements for cross type.")
                psr <- sr/sum(sr)
                props <- t(apply(ctab, 1, function(el) el/sum(el)))
                pmx <- apply(props, 1, max)
                pmn <- apply(props, 1, min)
                whm <- (pmx > max(psr)[1]) | (pmn < min(psr)[1])
            } else stop("Need a threshold or a ratio (see ?pullCross).")
        }
        if(type == "missing")
            whm <- tab$missing > pars$miss.thresh
        tab <- cbind.data.frame(mark = markn[whm], tab[whm,1:nc])
    }
    if(any(whm)){
        rownames(tab) <- NULL
        object <- drop.markers(object, markn[whm])
        if(is.null(object[[type]]) | replace)
            object[[type]]  <- list()
        object[[type]]$table <- rbind(object[[type]]$table, tab)
        object[[type]]$data <- cbind(object[[type]]$data, genot[,whm, drop = FALSE])
    } else message("No markers were found with these characteristics.")
    object$geno <- c(object$geno, oldObject$geno)
    object$geno <- object$geno[mixedorder(names(object$geno))]
    object
}

pushCross <- function(object, type = c("co.located","seg.distortion","missing","unlinked"), unlinked.chr = NULL, pars = NULL, ...){
    fixObject <- function(x, type){
        cc <- class(x)
        x <- unclass(x)
        x[[type]] <- NULL
        x <- x[!sapply(x, is.null)]
        class(x) <- cc
        x }
    if(type != "unlinked"){
       if(dim(object$pheno)[1] != dim(object[[type]]$data)[1])
           stop("Number of genotypes in linkage map does not match external ", type, " data.")
    }
    cs <- attr(object, "scheme")
    if(is.null(pars))
         pars <- pp.init()
    if(!is.list(pars))
        stop("Argument pars must be a list object with one or more named elements matching the results of pp.init()")
    pars <- do.call("pp.init", pars)
    type <- match.arg(type)
    if(!(class(object)[1] %in% c("bc","dh","riself","bcsft")))
        stop("Pushing of markers is not supported for this population type, see ?pushCross.")
    if(is.null(object[[type]]) & !(type %in% "unlinked"))
        stop("There are no markers of this type to push back.")
    if(type == "co.located"){
        odat <- object$co.located$data
        mtab <- object$co.located$table
        mnam <- markernames(object)
        bm <- split(mtab, mtab$bins)
        om <- sapply(bm, function(el) as.character(el$mark[1]))
        if(!all(mwh <- om %in% mnam)){
            warning("Some co-located marker anchors do not exist in map .. using matching set only.")
            bm <- bm[mwh]
            om <- om[mwh]
        }
        chrn <- rep(names(nmar(object)), times = nmar(object))
        nb <- sapply(bm, nrow)
        chrn <- chrn[pmatch(om, mnam)]
        uc <- unique(chrn)
        for(i in uc){
            bmc <- bm[chrn %in% i]
            markc <- as.character(unlist(sapply(bmc, function(el) el$mark[2:length(el$mark)])))
            omap <- object$geno[[i]]$map
            for (j in 1:length(bmc)){
                whm <- (1:length(omap))[names(omap) == bmc[[j]]$mark[1]]
                nmap <- c(omap[1:whm], rep(omap[whm], length(bmc[[j]]$mark) -
                                           1))
                if (length(omap) != whm)
                    nmap <- c(nmap, omap[(whm + 1):length(omap)])
                names(nmap)[(whm + 1):(whm + length(bmc[[j]]$mark) -
                                       1)] <- as.character(bmc[[j]]$mark)[-1]
                omap <- nmap
            }
            whc <- dimnames(odat)[[2]] %in% markc
            object$geno[[i]]$data <- cbind(object$geno[[i]]$data,
                                             odat[, whc, drop = FALSE])
            object$geno[[i]]$data <- object$geno[[i]]$data[, names(nmap)]
            object$geno[[i]]$map <- omap
        }
        object <- fixObject(object, type)
        return(object)
    }
    if(type %in% c("seg.distortion","missing")){
        push.object <- oe <- object
        push.object$geno <- list()
        whm <- rep(TRUE, ncol(object[[type]]$data))
        if(type == "seg.distortion"){
            tabs <- object$seg.distortion$table
            if(!is.null(seg.thresh <- pars$seg.thresh)){
                pval <- 10^(-tabs$neglog10P)
                if(length(seg.thresh) == 1)
                    whm <- pval > seg.thresh
                else if(length(seg.thresh) == 2)
                  whm <- (pval > min(seg.thresh)) & (pval < max(seg.thresh))
                else stop("A numerical seg.thresh should contain no more than two elements.")
            }
            else if(!is.null(seg.ratio <- pars$seg.ratio)) {
                sr <- as.numeric(unlist(strsplit(seg.ratio, ":")))
                ptab <- tabs[,6:ncol(tabs)]
                if(length(sr) != ncol(ptab))
                    stop("Wrong number of ratio elements for cross type.")
                psr <- sr/sum(sr)
                ptab <- t(apply(ptab, 1, function(el) el/sum(el)))
                pmx <- apply(ptab, 1, max)
                pmn <- apply(ptab, 1, min)
                whm <- pmx < max(psr)[1] | pmn > min(psr)[1]
            } else stop("Need a threshold or ratio (see ?pushCross).")
        }
        if(type %in% "missing"){
            tabs <- object$missing$table
            whm <- tabs$missing < pars$miss.thresh
        }
        if(any(whm)){
            inds <- (1:ncol(object[[type]]$data))[whm]
            push.object$geno[[type]]$data <- object[[type]]$data[,whm, drop = FALSE]
            push.object$geno[[type]]$map <- 1:ncol(push.object$geno[[type]]$data)
            names(push.object$geno[[type]]$map) <- dimnames(push.object$geno[[type]]$data)[[2]]
            class(push.object$geno[[type]]) <- class(object$geno[[1]])
        } else stop("There are no markers to push back with these characteristics")
        object[[type]]$data <- object[[type]]$data[,!whm, drop = FALSE]
        object[[type]]$table <- object[[type]]$table[!whm,]
        if(dim(object[[type]]$data)[2] == 0)
            object <- fixObject(object, type)
        object$geno <- c(object$geno, push.object$geno)
    }
    if(type %in% "unlinked"){
        if(is.null(unlinked.chr))
            stop("Argument unlinked.chr cannot be NULL.")
        if(!all(unlinked.chr %in% names(nmar(object))))
            stop("Some names in unlinked.chr do not match linkage group names of object.")
        push.object <- subset(object, chr = unlinked.chr)
        oe <- subset(object, paste("-", unlinked.chr, sep =""))
    }
    oe$geno <- lapply(oe$geno, function(el){
        hm <- floor(el$map[length(el$map)]/10)
        if(hm != 0){
            whp <- seq(el$map[1], el$map[length(el$map)], length.out = 2 + hm)
            whp <- whp[2:(length(whp) - 1)]
            pm <- c()
            for(i in 1:length(whp)){
                dm <- abs(el$map - whp[i])
                pm[i] <- (1:length(el$map))[dm == min(dm)][1]
            }
        } else pm <- NULL
        whp <- unique(c(1,pm,length(el$map)))
        el$data <- el$data[,whp, drop = FALSE]
        el$map <- el$map[whp]
        el })
    nme <- totmar(oe)
    mn <- markernames(push.object)
    chrn <- rep(names(nmar(oe)), times = nmar(oe))
    oe$geno <- c(oe$geno, push.object$geno)
    chru <- paste("UL", 1:length(mn), sep = "")
    chro <- c(chrn, chru)
    erf <- est.rf(oe)$rf
    lod <- erf
    lod[lower.tri(erf)] <- t(erf)[lower.tri(erf)]
    erf[upper.tri(erf)] <- t(erf)[upper.tri(erf)]
    diag(erf) <- 1; diag(lod) <- 0;
    tme <- totmar(oe)
    ind <- (nme + 1):tme
    for(i in (nme + 1):tme){
        if(!(chro[i] %in% unique(chrn))){
            wh <- (erf[,i] <= pars$max.rf) & (lod[,i] > pars$min.lod)
            link <- chro %in% chrn
            whl <- wh[link]; whr <- wh[!link]
            if(any(whl) | any(whr)){
                if(any(whl)){
                    chrt <- chro[link][whl]
                    erfs <- erf[,i][link][whl]
                }
                else if(any(whr)) {
                    chrt <- chro[!link][whr]
                    erfs <- erf[,i][!link][whr]
                }
                chi <- chrt[erfs == min(erfs)][1]
                chro[!link][whr] <- chro[chro == chro[i]] <- chi
            }
        }
    }
    for(i in 1:length(mn))
        object <- movemarker(object, mn[i], chro[nme + i])
    nm <- nmar(object)
    lenu <- grep("UL",names(nm))
    names(object$geno)[lenu] <- paste("UL",1:length(lenu),sep= "")
    object$geno <- object$geno[mixedorder(names(object$geno))]
    attr(object, "scheme") <- cs
    object

}

combineMap <- function(..., id = "Genotype", keep.all = TRUE){
    mapl <- list(...)
    marku <- unlist(lapply(mapl, function(el) markernames(el)))
    tm <- table(marku)
    if(any(tm > 1))
        stop("Non-unique markers between linkage maps.")
    if(length(unique(sapply(mapl, function(el) class(el)[1]))) > 1)
        stop("Classes of maps need to be identical.")
    scheme <- lapply(mapl, function(el) attr(el, "scheme"))
    if(any(!sapply(scheme, is.null))){
        sc <- do.call("rbind", scheme)
        if((nrow(sc) != length(mapl)) | !all(duplicated(sc)[2:nrow(sc)]))
            stop("Mismatched cross schemes in linkage maps")
    }
    if(!all(sapply(mapl, function(el) id %in% names(el$pheno))))
        stop("Some linkage maps do not contain column\"", id, "\".")
    mapl <- lapply(mapl, function(el){
        names(el$geno) <- gsub("x","X", names(el$geno))
        el})
    maplb <- lapply(mapl, function(el){
        mapb <- do.call("cbind", lapply(el$geno, function(x) x$data))
        mdist <- unlist(pull.map(el))
        chrs <- rep(names(el$geno), times = nmar(el))
        cl <- sapply(el$geno, function(el) class(el))
        dimnames(mapb)[[2]] <- paste(chrs, markernames(el), as.character(mdist), cl, sep = ";")
        mapb <- as.data.frame(mapb)
        mapb[[id]] <- as.character(el$pheno[[id]])
        mapb
    })
    phelb <- lapply(mapl, function(el) el$pheno)
    mapm <- maplb[[1]]
    phem <- phelb[[1]]
    for(i in 1:(length(maplb) - 1)) {
        mapm <- merge(mapm, maplb[[i + 1]], by = id, all = keep.all)
        phem <- merge(phem, phelb[[i + 1]], by = id, all = keep.all)
    }
    nams <- mapm[[id]]
    mxo <- mixedorder(nams)
    mapm <- mapm[mxo,]
    phem <- phem[mixedorder(phem[[id]]),,drop = FALSE]
    mapm <- mapm[,!(names(mapm) %in% id)]
    spl.m <- strsplit(names(mapm), ";")
    ch <- sapply(spl.m, "[", 1)
    mapm <- t(mapm)
    mapf <- list()
    mapf$geno <- lapply(split.data.frame(mapm, ch), function(el, nams){
        temp <- list()
        spl.c <- strsplit(rownames(el), ";")
        temp$data <- as.matrix(t(el))
        rownames(temp$data) <- nams
        temp$map <- as.numeric(sapply(spl.c, "[", 3))
        names(temp$map) <- dimnames(temp$data)[[2]] <- sapply(spl.c, "[", 2)
        mo <- order(temp$map)
        temp$map <- temp$map[mo]
        temp$data <- temp$data[,mo, drop = FALSE]
        class(temp) <- unique(sapply(spl.c, "[", 4))
        temp
    }, nams[mxo])
    class(mapf) <- class(mapl[[1]])
    attr(mapf, "scheme") <- scheme[[1]]
    mapf$pheno <- phem
    mapf$geno <- mapf$geno[mixedorder(names(mapf$geno))]
    mapf
}

## combineMap <- function(..., id = "Genotype", keep.all = TRUE, merge.by = "genotype"){
##     mapl <- list(...)
##     if(merge.by %in% "genotype"){
##         marku <- unlist(lapply(mapl, function(el) markernames(el)))
##         if(any(table(marku) > 1))
##             stop("Non-unique markers between linkage maps.")
##     } else {
##         genu <- unlist(lapply(mapl, function(el) as.character(el$pheno[[id]])))
##         if(any(table(genu) > 1))
##             stop("Non-unique genotypes between linkage maps.")
##     }
##     if(length(unique(sapply(mapl, function(el) class(el)[1]))) > 1)
##         stop("Classes of maps need to be identical.")
##     scheme <- lapply(mapl, function(el) attr(el, "scheme"))
##     if(any(!sapply(scheme, is.null))){
##         sc <- do.call("rbind", scheme)
##         if((nrow(sc) != length(mapl)) | !all(duplicated(sc)[2:nrow(sc)]))
##             stop("Mismatched cross schemes in linkage maps")
##     }
##     if(!all(sapply(mapl, function(el) id %in% names(el$pheno))))
##         stop("Some linkage maps do not contain column\"", id, "\".")
##     mapl <- lapply(mapl, function(el){
##         names(el$geno) <- gsub("x","X", names(el$geno))
##         el})
##     maplb <- lapply(mapl, function(el, merge.by){
##         mapb <- do.call("cbind", lapply(el$geno, function(x) x$data))
##         mdist <- unlist(pull.map(el))
##         chrs <- rep(names(el$geno), times = nmar(el))
##         cl <- sapply(el$geno, function(el) class(el))
##         dimnames(mapb)[[2]] <- paste(chrs, markernames(el), as.character(mdist), cl, sep = ";")
##         dimnames(mapb)[[1]] <- as.character(el$pheno[[id]])
##         if(merge.by %in% "marker"){
##             mapb <- as.data.frame(t(mapb))
##             mapb[[merge.by]] <- markernames(el)
##             mapb[["dims"]] <- dimnames(mapb)[[1]]
##         } else {
##             mapb <- as.data.frame(mapb)
##             mapb[[merge.by]] <- as.character(el$pheno[[id]])
##         }
##         mapb
##     }, merge.by)
##     phelb <- lapply(mapl, function(el) el$pheno)
##     mapm <- maplb[[1]]
##     phem <- phelb[[1]]
##     for(i in 1:(length(maplb) - 1)) {
##         mapm <- merge(mapm, maplb[[i + 1]], by = merge.by, all = keep.all)
##         phem <- merge(phem, phelb[[i + 1]], by = id, all = keep.all)
##     }
##     if(length(wh <- grep("dims", names(mapm)))){
##         dnam <- apply(mapm[,wh], 1, function(el){
##             el <- el[!is.na(el)]
##             el[1] })
##         omit <- c(names(mapm)[wh], "marker")
##         mapm <- mapm[,!(names(mapm) %in% omit)]
##         rownames(mapm) <- dnam
##     } else {
##         rownames(mapm) <- mapm[["genotype"]]
##         mapm <- t(mapm[,!(names(mapm) %in% "genotype")])
##     }
##     nams <- dimnames(mapm)[[2]]
##     mxo <- mixedorder(nams)
##     mapm <- mapm[,mxo]
##     phem <- phem[mixedorder(phem[[id]]),,drop = FALSE]
##     spl.m <- strsplit(rownames(mapm), ";")
##     ch <- sapply(spl.m, "[", 1)
##     mapf <- list()
##     mapf$geno <- lapply(split.data.frame(mapm, ch), function(el, nams){
##         temp <- list()
##         spl.c <- strsplit(rownames(el), ";")
##         if(nrow(el) == 1) print(as.matrix(t(el)))
##         temp$data <- as.matrix(t(el))
##         rownames(temp$data) <- nams
##         temp$map <- as.numeric(sapply(spl.c, "[", 3))
##         names(temp$map) <- dimnames(temp$data)[[2]] <- sapply(spl.c, "[", 2)
##         mo <- order(temp$map)
##         temp$map <- temp$map[mo]
##         temp$data <- temp$data[,mo, drop = FALSE]
##         class(temp) <- unique(sapply(spl.c, "[", 4))
##         temp
##     }, nams[mxo])
##     class(mapf) <- class(mapl[[1]])
##     attr(mapf, "scheme") <- scheme[[1]]
##     mapf$pheno <- phem
##     mapf$geno <- mapf$geno[mixedorder(names(mapf$geno))]
##     mapf
## }

statMark <- function(cross, chr, stat.type = c("marker","interval"), map.function = "kosambi"){
    if (!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")
    if(!(class(cross)[1] %in% c("bc","dh","riself","bcsft")))
        stop("This function is not suitable for this population type, see ?statMark.")
    if(any(!(stat.type %in% c("marker","interval"))))
        stop("Value for stat.type argument does not match allowable names, see ?statMak.")
    if (!missing(chr))
        cross <- subset(cross, chr = chr)
    nm <- nmar(cross)
    mrk <- list()
    if("marker" %in% stat.type) {
        mrk$marker <- geno.table(cross, scanone.output = TRUE)
        mrk$marker$dxo <- unlist(lapply(cross$geno, function(el) {
            if(ncol(el$data) < 3) sxo <- rep(0, ncol(el$data))
            else {
                sxo <- sapply(3:ncol(el$data), function(i, el) {
                    out <- el$data[,i-2] == el$data[,i]
                    left <- el$data[,i-2] != el$data[,i-1]
                    sm <- out + left
                    sm[sm != 2] <- 0
                    sm[sm == 2] <- 1
                    sum(sm, na.rm = TRUE)
                }, el)
                c(0, sxo, 0)
            }}))
    }
    if("interval" %in% stat.type){
        ic <- sapply(cross$geno, function(el) length(el$map) > 1)
        icross <- subset(cross, chr = names(nmar(cross))[ic])
        nmo <- nmar(icross)
        rf <- lod <- c()
        for(i in 1:length(nmo)){
            sc <- subset(icross, chr = names(nmo)[i])
            erf <- est.rf(sc)$rf
            snm <- nmo[i]
            inds <- cbind(1:(snm - 1), 2:snm)
            rf <- c(rf, erf[inds[,2:1, drop = FALSE]])
            lod <- c(lod, erf[inds])
        }
        mrk$interval <- cbind.data.frame(erf = rf, lod = lod)
        mrk$interval$dist <-  unlist(lapply(icross$geno, function(el) diff(el$map)))
        mf <- switch(map.function, kosambi = mf.k, haldane = mf.h,
                                   morgan = mf.m, cf = mf.cf)
        mrk$interval$mrf <- mf(mrk$interval$dist)
        mrk$interval$recomb <- unlist(lapply(icross$geno, function(el){
            sapply(2:ncol(el$data), function(i, el)
                  sum(abs(el$data[,i] - el$data[,i - 1]) > 0, na.rm = TRUE), el)
        }))
        rownames(mrk$interval) <- unlist(lapply(icross$geno, function(el){
                  len <- length(el$map)
                  paste("(",names(el$map)[1:(len - 1)],",", names(el$map)[2:len],")", sep = "")
              }))
        chr <- rep(names(nmo), times = nmo - 1)
        pos <- unlist(lapply(icross$geno, function(el) el$map[1:(length(el$map)-1)] + diff(el$map)/2), use.names = FALSE)
        mrk$interval <- cbind.data.frame(chr = chr, pos = pos, mrk$interval)
        mrk$interval$chr <- as.character(mrk$interval$chr)
        if(any(nm == 1)){
            sc1 <- names(nm)[nm == 1]
            mark1 <- markernames(subset(cross, chr = sc1))
            for(i in 1:length(sc1)) {
                mrk$interval[nrow(mrk$interval) + 1,1] <- sc1[i]
                mrk$interval$pos[nrow(mrk$interval)] <- 0
                rownames(mrk$interval)[nrow(mrk$interval)] <- paste("(",mark1[i],")", sep = "")
            }
            mrk$interval <- mrk$interval[mixedorder(mrk$interval$chr),]
        }
    }
    mrk
}

profileMark <- function(cross, chr, stat.type = "marker", use.dist = TRUE, map.function = "kosambi",
                        crit.val = NULL, display.markers = FALSE, mark.line = FALSE, ...){
    if (!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")
    if(!(class(cross)[1] %in% c("bc","dh","riself","bcsft")))
        stop("This function is not suitable for this population type, see ?profileMark.")
    if (!missing(chr))
        cross <- subset(cross, chr)
    cross$geno <- cross$geno[mixedorder(names(nmar(cross)))]
    stypes <- c("seg.dist","miss","prop","dxo","erf","lod","dist","mrf","recomb")
    mtypes <- c("marker","interval")
    if(any(!(stat.type %in% c(mtypes, stypes))))
        stop("Value for stat.type argument does not match allowable names, see ?profileMark.")
    rmtypes <- rep(mtypes, times = c(4,5))
    if(!any(stypes %in% stat.type)){
        call.type <- mtypes[mtypes %in% stat.type]
        stat.type <- stypes[rmtypes %in% call.type]
    }
    else call.type <- unique(rmtypes[stypes %in% stat.type])
    if(!is.null(crit.val)){
        if(crit.val != "bonf")
            stop("Argument crit.val can only be \"bonf\".")
    }
    stat.type <- stypes %in% stat.type
    stat <- stato <- statMark(cross, names(nmar(cross)), call.type, map.function)
    nm <- nmar(cross)
    if(use.dist) {
        dist <- lapply(cross$geno, function(el) {
            dist <- diff(el$map)/100
            tel <- c(0.05, dist)
            names(tel)[1] <- names(el$map)[1]
            tel})
        dist <- cumsum(unlist(dist))
    }
    else dist <- 1:sum(nm)
    adist <- dist
    mn <- markernames(cross)
    chrn <- factor(rep(names(nm), times = nm), levels = names(nm))
    plis <- list
    if("marker" %in% names(stat)) {
        nlp <- stat$marker$neglog10P
        mtype <- stat.type[1:4]
        nc <- ncol(stat$marker)
        mtype <- c(mtype[1:2],rep(mtype[3],nc - 5), mtype[4])
        props <- paste("Prop. ",names(summary(cross)$typing.freq), sep = "")
        mnam <- c("Seg. Distortion","Missing",props,"Double Crossovers")[mtype]
        mnam <- paste("Marker: ", mnam, sep = "")
        stat$marker <- cbind.data.frame(val = unlist(stat$marker[,(3:nc)[mtype]]))
        stat$marker$dist <- rep(dist, length(mnam))
        stat$marker$nams <- rep(mnam, each = length(dist))
        stat$marker$lg <- rep(chrn, length(mnam))
        if(!is.null(crit.val)){
            crit.val <- 0.05/totmar(cross)
            critm <- mn
            logmc <- 10^(-nlp) > crit.val
            critm[logmc] <- ""
            stato$marker$crit.val <- logmc
            stat$marker$mark <- rep(critm, length(mnam))
        }
    }
    if("interval" %in% names(stat)){
        ni <- nm - 1
        ni[ni == 0] <- 1
        itype <- stat.type[5:9]
        inam <- c("Est. Recomb. Frac.", "LOD", "Map Dist.", "Map Recomb. Frac.","# Recombinations")[itype]
        inam <- paste("Interval: ", inam, sep = "")
        intn <- rownames(stat$interval)
        lod <- stat$interval$lod
        stat$interval <- cbind.data.frame(val = unlist(stat$interval[,(3:7)[itype]]))
        idist <- unlist(lapply(split(dist, chrn), function(el){
            if(length(el) == 1) el
            else el[1:(length(el)-1)] + diff(el)/2
        }))
        stat$interval$nams <- rep(inam, each = length(idist))
        stat$interval$dist <- rep(idist, length(inam))
        stat$interval$lg <- rep(rep(names(nm), times = ni), length(inam))
        if(!is.null(crit.val)){
            crit.val <- -log10(0.05/(sum(nm - 1)))
            criti <- intn
            logic <- lod > crit.val
            criti[logic] <- ""
            stato$interval$crit.val <- logic
            stat$interval$mark <- rep(criti, length(inam))
        }
    }
    stat <- do.call("rbind.data.frame", stat)
    rownames(stat) <- NULL
    if(!display.markers){
        labs <- rep("",length(adist))
        wh <- c(1, cumsum(nm)[1:(length(nm) - 1)] + 1)
        labs[wh] <- paste(names(nm), "1", sep = ".")
    } else labs <- mn
    print(xyplot(val ~ dist | nams, data = stat, groups = stat$lg,
                 panel = function(x, y, groups = groups, subscripts = subscripts, crit.mark, nams, adist, mark.line, ...) {
                     if(mark.line)
                         panel.abline(v = adist, lty = 2, col = gray(0.8))
                     if(length(grep("Prop.", unique(nams[subscripts]))))
                         panel.abline(h = 0.5, lty = 2, col = gray(0.5))
                     if(!is.null(crit.mark)){
                         panel.text(x, y, labels = crit.mark[subscripts], ...)
                     }
                     panel.superpose(x, y, groups = groups, subscripts = subscripts, ...)
                 }, ylab = "", scales = list(y = list(relation = "free"),
                  x = list(labels = labs, rot = 45, at = adist, cex = 0.6)), xlab = "",
                 mark.line = mark.line, crit.mark = stat$mark, nams = stat$nams, adist = adist, ...))
    invisible(stato)
}

statGen <- function(cross, chr, bychr = TRUE, stat.type = c("xo","dxo","miss"), id = "Genotype"){
    if (!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")
    if(!(class(cross)[1] %in% c("bc","dh","riself","bcsft")))
        stop("This function not suitable for this population type, ?see statGen.")
    if(!(id %in% names(cross$pheno)))
        stop("The unique identifier for the genotypes, ", deparse(substitute(id)),
             ", cannot be found in the object")
    if(any(!(stat.type %in% c("xo","dxo","miss"))))
        stop("Value for stat.type argument does not match allowable names, see ?statGen.")
    if (!missing(chr))
        cross <- subset(cross, chr = chr)
    nch <- nchr(cross)
    nm <- nmar(cross)
    gnam <- as.character(cross$pheno[[id]])
    cnt <- list()
    if("xo" %in% stat.type) {
        cc <- sapply(cross$geno, function(el) length(el$map) > 1)
        ccross <- subset(cross, chr = names(nmar(cross))[cc])
        cnt$xo <- do.call("cbind", lapply(ccross$geno, function(el)
             apply(el$data, 1, function(rows) {
                 ad <- abs(diff(rows[!is.na(rows)]))
                 length(ad[ad != 0])
             })))
    }
    if("dxo" %in% stat.type){
        pc <- sapply(cross$geno, function(el) length(el$map) > 2)
        dcross <- subset(cross, chr = names(nmar(cross))[pc])
        cnt$dxo <- do.call("cbind", lapply(dcross$geno, function(el) {
            sxo <- sapply(3:ncol(el$data), function(i, el) {
                out <- el$data[,i-2] == el$data[,i]
                left <- el$data[,i-2] != el$data[,i-1]
                sm <- out + left
                sm[sm != 2] <- 0
                sm[sm == 2] <- 1
                sm
            }, el)
            apply(sxo, 1, sum, na.rm = TRUE)
        }))
    }
    if("miss" %in% stat.type)
        cnt$miss <- do.call("cbind", lapply(cross$geno, function(el)
            apply(el$data, 1, function(rows) {
                rows[is.na(rows)] <- 0
                length(rows[rows == 0])
            })))
    cnt <- lapply(cnt, function(el, gnam, bychr){
        if(!bychr) {
            el <- apply(el, 1, sum)
            names(el) <- gnam
        } else rownames(el) <- gnam
        el
    }, gnam, bychr)
    cnt
}

profileGen <- function(cross, chr, bychr = TRUE, stat.type = c("xo","dxo","miss"), id = "Genotype", xo.lambda = NULL, ...){
    if (!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")
    if(!(class(cross)[1] %in% c("bc","dh","riself","bcsft")))
        stop("This function not suitable for this population type, see ?profileGen.")
    if(!(id %in% names(cross$pheno)))
         stop("The unique identifier for the genotypes, ", deparse(substitute(id)),
              ", cannot be found in the object")
    if (!missing(chr))
        cross <- subset(cross, chr = chr)
     nm <- nmar(cross)
     cnt <- statGen(cross, names(nm), bychr, stat.type, id)
     pnam <- c("Crossovers","Double Crossovers","Missing")
     stat.type <- pnam[pmatch(names(cnt), c("xo","dxo","miss"))]
     ph <- as.character(cross$pheno[[id]])
     if(("xo" %in% names(cnt)) & !is.null(xo.lambda)) {
         if(bychr) cxo <- apply(cnt$xo, 1, sum) else cxo <- cnt$xo
         xo.lambda <- ppois(cxo, xo.lambda, lower.tail = FALSE) < 0.05/length(ph)
     }
     nch <- nchr(cross)
     inter <- floor(seq(1, length(ph), length.out = 20))
     labs <- ph[inter]
     dots <- list(...)
    if (!is.na(pmatch("cex", names(dots))))
        p.cex <- dots$cex else p.cex <- par("cex")
    if(bychr){
        no.splits <- nch*length(cnt)
        stat.type <- paste(rep(names(nm),length(cnt)), rep(stat.type, each = nch), sep = ": ")
    } else no.splits <- length(cnt)
     stat.type <- rep(stat.type, each = length(ph))
     phv <- rep(ph, no.splits)
     noind <- rep(1:length(ph), no.splits)
     cntl <- unlist(cnt)
     form <- cntl ~ noind | stat.type
     grp <- rep(1, length(noind))
     cntn <- cntl; noindn <- noind
     plot.xo <- FALSE
     if(!is.null(xo.lambda) & any(xo.lambda)){
        inds <- (1:length(phv))[rep(xo.lambda,no.splits)]
        grp[inds] <- 2
        cntn[-inds] <- noindn[-inds] <- phv[-inds] <- NA
        ylims <- lapply(split(cntl, stat.type), function(el) c(0, max(el)[1] + 0.15*max(el)[1]))
        plot.xo <- TRUE
    } else ylims <- lapply(split(cntl, stat.type), function(el) c(0, max(el)[1]))
    xo.panel <- function(x, y, subscripts = subscripts, cntn, noindn, phv, plot.xo, ...){
        panel.xyplot(x, y, ...)
        panel.segments(x, 0, x, y, ...)
        if(plot.xo)
            panel.text(noindn[subscripts], cntn[subscripts], labels = phv[subscripts], adj = c(-0.2,-0.2), srt = 45, ...)
    }
    print(xyplot(form, groups = factor(grp),
                 panel = panel.superpose, panel.groups = xo.panel,
                 cntn = cntn, noindn = noindn, phv = phv, plot.xo = plot.xo,
                 scales = list(x = list(labels = labs, rot = 45, at = inter, cex = p.cex),
                     y = list(relation = "free")), xlab = "Genotypes", ylab = "",ylim = ylims, ...))
    invisible(list(stat = cnt, xo.lambda = xo.lambda))
}


