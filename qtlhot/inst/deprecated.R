#############################################################################################################
## Legacy routines to be purged.
#############################################################################################################
hotspot.scan <- function(cross, highobj, lod.thr = NULL, quant.level = NULL, window = NULL,
                         verbose = FALSE)
{
  hotsize(highobj, lod.thr, window, quant.level)
}
#############################################################################################################
hotspot.plot <- function(hotobj, quant.thr = NULL, maps = NULL, main = "")
{
  plot(hotobj, ..., by.chr = TRUE, quant.thr = quant.thr, maps = maps, title = main)
}
#############################################################################################################
##############################################################################
OGetCisCandReg <- function(cand.reg, scan, lod.thr, drop = 1.5) {
  all.trait.nms <- cand.reg[, 1]
  cand.reg <- cand.reg[cand.reg[, 2] == cand.reg[, 4],]
  n <- nrow(cand.reg)
  trait.nms <- names(scan)[-c(1, 2)]
  is.cis <- rep(FALSE, n)
  peak.pos.lower <- rep(NA, n)
  peak.pos.upper <- rep(NA, n)
  for (i in 1 : n) {
    phys.pos <- cand.reg[i, 3]
    trait.index <- which(trait.nms == cand.reg[i, 1]) + 2
    sscan <- scan[scan[, 1] == cand.reg[i, 4], c(1, 2, trait.index)]
    lod.interval <- lodint(sscan, chr = cand.reg[i, 4], drop = drop)
    peak.pos.lower[i] <- lod.interval[1, 2]
    peak.pos.upper[i] <- lod.interval[nrow(lod.interval), 2]
    lod.phys.pos <- sscan[which.min(abs(sscan[, 2] - cand.reg[i, 3])), 3]
    if (phys.pos >= peak.pos.lower[i] & phys.pos <= peak.pos.upper[i] & 
        lod.phys.pos >= lod.thr) {
      is.cis[i] <- TRUE
    }
  }
  out <- data.frame(cand.reg, peak.pos.lower, peak.pos.upper)
  index <- match(out[is.cis, 1], all.trait.nms)
  list(cis.reg = out[is.cis,], cis.index = index)
}
##############################################################################
OGetCandReg <- function(scan, annot, traits, lod.thr, drop = 1.5,
                        nms = names(scan)[-c(1,2)])
  
{
  ## want to use highlod instead of scan below.
  ## need to decode highlod.
  ## currently this only gets max over genome; want max by chr, yes?
  traits <- unique(traits)
  n <- length(traits)
  out <- data.frame(matrix(NA, n, 6))
  names(out) <- c("gene", "phys.chr", "phys.pos", "peak.chr", "peak.pos",
                  "peak.lod")

  for (i in 1:n) {
    ii <- match(traits[i], annot[, 1])
    out[i, 1:3] <- annot[ii, c(1, 3, 5)]
    trait.chr <- annot[ii, 3]
    trait.pos <- annot[ii, 5]
    if (!is.na(trait.pos)) {
      peak.index <- which.max(scan[, traits[i]])
      peak.lod <- scan[peak.index, traits[i]]
      peak.chr <- scan[peak.index, 1]
      peak.pos <- scan[peak.index, 2]
      if (peak.lod >= lod.thr) {
        out[i, 4] <- peak.chr
        out[i, 5] <- peak.pos
        out[i, 6] <- peak.lod
      }
    }     
    cat(" ", i, "\n")
  }
  subset(out, !is.na(out[, 4]))
}
##############################################################################
OGetCoMappingTraits <- function(traits, scan, lod.thr, drop = 1.5) {
  n <- nrow(traits)
  out <- vector(mode = "list", length = n)
  names(out) <- traits[, 1]
  nms1 <- names(scan)[-c(1, 2)]
  for (i in 1:n) {
    peak.pos <- traits[i, 5]
    sub.scan <- scan[scan[, 1] == traits[i, 4], ]
    peak.lods <- apply(sub.scan[, -c(1, 2)], 2, max)
    aux <- which(peak.lods >= lod.thr)
    aux <- aux[-match(traits[i, 1], names(aux))]
    nms2 <- names(aux)
    is.target <- rep(FALSE, length(nms2))
    for (j in 1:length(nms2)) {
      trait.index <- match(nms2[j], nms1)
      sscan <- sub.scan[, c(1, 2, trait.index + 2)]
      lod.interval <- lodint(sscan, chr = traits[i, 4], drop = drop)
      lb <- lod.interval[1, 2]
      ub <- lod.interval[nrow(lod.interval), 2] ## was 3; error if there are tied LOD scores.
      if(peak.pos >= lb & peak.pos <= ub)
        is.target[j] <- TRUE
    }
    out[[i]] <- nms2[is.target]
    cat(" ", i, "\n")
  }
  out
}
##############################################################################
## Run permuations
lod.quantile.permutations.2 <- function(cross, pheno.col, N, n.perm, lod.thr, 
                                        batch.effect, window, seed = 123456789, verbose = FALSE)
{

  ## Set up pseudo-random number generation seeds.
  set.seed(seed[[1]])
  all.seeds <- sample(c(98765:987654), n.perm, replace=FALSE)

  ## Set up matrices to record values.
  nphe.cross <- nphe(cross)
  N <- N[N <= nphe.cross]
  quants <- 1 - (N - 1)/nphe.cross
  l.N <- length(N)
  max.lod.quant <- matrix(NA, n.perm, l.N)
  dimnames(max.lod.quant)[[2]] <- paste(paste(round(quants,4)*100, "%", sep=""), 
                                        N, sep="_")
  l.lod.thr <- length(lod.thr)
  max.N <- matrix(NA, n.perm, l.lod.thr)
  max.lod.quant <- matrix(NA, n.perm, l.N)
  max.N.window <- matrix(NA, n.perm, l.lod.thr)

  if(!is.null(batch.effect))
    cross <- subset(cross, ind = !is.na(cross$pheno[[batch.effect]]))

  covars <- sexbatch.covar(cross, batch.effect)

  n.ind <- nind(cross)

  for(i in 1:n.perm){
    perm.cross <- cross
    set.seed(all.seeds[i])
    tmp <- sample(c(1:n.ind), n.ind, replace=FALSE)
    perm.cross$pheno <- cross$pheno[tmp,]

    per.scan <- scanone(perm.cross, pheno.col=pheno.col, method="hk", 
                        addcovar=covars$addcovar, intcovar=covars$intcovar)

    ## Elias' quantiles.
    quant <- apply(per.scan[,-(1:2)], 1, quantile, quants)
    max.lod.quant[i,] <- apply(quant,1,max)

    ## Jansen's count.
    max.N[i,] <- apply(count.thr(per.scan, lod.thr, droptwo = TRUE), 1, sum)

    ## Smoothed count.
    maxlod <-  apply(per.scan[,-(1:2)], 2, tapply, per.scan[,1], max)
    chrs <- dimnames(maxlod)[[1]]
    chr <- factor(per.scan[,1], levels = chrs)
    pos <- as.numeric(per.scan[,2])

    maxlod.pos <- maxlod.thr.pos <- vector("list", length(chrs))
    names(maxlod.pos) <- names(maxlod.thr.pos) <- chrs
    for(k in seq(length(chrs))) {
      scan.out.bychr <- per.scan[per.scan[,1] == chrs[k], ]
      maxlod.pos[[k]] <- apply(scan.out.bychr[,-(1:2)], 2, function(a,b)
                               mean(b[!is.nan(a) & a==max(a, na.rm=TRUE)]), scan.out.bychr[,2])
    }

    for(j in seq(along = lod.thr)){
      for(k in chrs)
        maxlod.thr.pos[[k]] <- maxlod.pos[[k]][maxlod[k,] >= lod.thr[j]]
      neqtl.pos <- smoothall(maxlod.thr.pos, thechr = chr, thepos = pos, window = window)
      max.N.window[i,j] <- max(neqtl.pos[,3])
    }
    print(i)
  }
  list(max.lod.quant=max.lod.quant, 
       max.N=max.N,
       max.N.window=max.N.window)
}
##############################################################################
OWW.perm <- function(scanmat, lod.thrs, n.perm, verbose = FALSE)
{
  ## Permute all but first column.
  nlod <- length(lod.thrs)
  max.WW <- matrix(NA, n.perm, nlod)
  dimnames(max.WW) <- list(NULL, as.character(lod.thrs))

  for(i in 1:n.perm) {
    if(verbose)
      cat("\n", i)

    ## Separately permute columns of the scan object separately by trait.
    mycat("sample...", verbose, last = "")
    scans <- apply(scanmat, 2, sample)
    ## Original code had reshuffle of same values. Better to start over each time.
    ##scanmat <- apply(scanmat, 2, sample)

    ## Get the maximum spurious hotspot size (N-method) across genome
    ## for different QTL mapping significance levels.
    mycat("sum...", verbose, last = "")
    max.WW[i,] <- apply(count.thr(scans, lod.thrs, FALSE), 1, max)
  }
  max.WW
}
##############################################################################
ONL.N.perm <- function(cross, Nmax, n.perm, lod.thrs, drop = 1.5,
                      verbose = FALSE, init.seed = 0) 
{
  set.seed(init.seed)
  n.phe <- nphe(cross)
  n.ind <- nind(cross)
  Nmax <- min(Nmax, n.phe)
  
  Ns <- seq(Nmax)
  quants <- 1 - (Ns - 1) / n.phe
  n.lod <- length(lod.thrs)

  max.N <- matrix(NA, n.perm, n.lod)
  dimnames(max.N) <- list(NULL, as.character(lod.thrs))

  max.lod.quant <- matrix(NA, n.perm, Nmax)
  dimnames(max.lod.quant) <- list(NULL, as.character(Ns))

  for(i in 1:n.perm){
    if(verbose)
      cat("\n", i)
    
    ## permute rows of the phenotype data matrix
    perm.cross <- cross
    tmp <- sample(c(1:n.ind), n.ind, replace=FALSE)
    perm.cross$pheno <- cross$pheno[tmp,]

    ## perform mapping analysis in the permuted data
    mycat("scanone...", verbose, last = "")
    ## NB: scanone groups phenos in batches based on missing data patterns.
    scanmat <- scanone(perm.cross, pheno.col = c(1:n.phe), method = "hk")
    chr <- scanmat[, 1]
    scanmat <- as.matrix(scanmat[, -(1:2), drop = FALSE])

    ## apply LOD drop interval to per.scan
    mycat("LOD drop...", verbose, last = "")
    scanmat <- set.to.zero.beyond.drop.int(chr, scanmat, min(lod.thrs), drop)

    ## get lod quantiles at each genomic location
    ## rows indexes the quants associated with the Ns
    ## columns indexes the genomic positions
    mycat("quantile...", verbose, last = "")
    quant <- apply(scanmat, 1, quantile, quants)

    ## get the maximum lod-quantile across the genome
    ## rows indexes the permutations
    ## columns indexes the Ns
    mycat("quant...", verbose, last = "")
    max.lod.quant[i,] <- apply(quant, 1, max)

    ## Get the maximum spurious hotspot size (N-method) across genome
    ## for different QTL mapping significance levels.
    mycat("max...", verbose, last = "")
    max.N[i,] <- apply(count.thr(scanmat, lod.thrs, FALSE), 1, max)
  }
  list(max.lod.quant = max.lod.quant, max.N = max.N)
}
######################################################################

## This function takes a scanone object (scan) and for each trait and each
## chromosome it: (1) finds out the max LOD peak; (2) if the max LOD value
## is smaller than the map threshold (thr) it set to zero all LOD values;
## (3) if the max LOD value is higher than thr, it computes the LOD drop 
## interval around the peak and set to zero all LOD values outside the LOD 
## LOD interval.
## 
set.to.zero.beyond.drop.int <- function(chr, scanmat, thr, drop = 1.5)
{
  mylodint <- function(x, drop = 1.5)
  {
    drops <- x < (max(x) - drop)
    if(any(drops))
      x[drops] <- 0
    x
  }

  mychr <- levels(chr)
  for(i in mychr) {
    pointer <- which(chr == i)
    if(length(pointer)) {
      if(max(scanmat[pointer,]) < thr)
        ## Zero out if no peak above threshold.
        scanmat[pointer,] <- 0
      else
        ## Zero out below drop from peak.
        scanmat[pointer,] <-
          apply(scanmat[pointer,, drop = FALSE], 2, mylodint, drop)
    }
  }
  scanmat
}
######################################################################
qtlhot.scan <- function(...) quant.scanone(...)
quant.plot <- function(...)
{
  x <- quant.lod(...)

  thr.level <- x$thr.level
  lod.thr <- x$lod.thr
  n.probs <- length(probs)

  tmp.plot <- function(x.vals, quant, x.crit, probs, level, is.quantile = FALSE, main = "",
                       add.level = FALSE)
  {
    n.probs <- length(probs)
    thr.level <- min(which(probs >= level))
    quant.thr <- rev(quant[thr.level,])[thr.level]

    xlabs <- "single trait LOD threshold"
    if(is.quantile)
      xlabs <- paste(xlabs, "quantile")
    
    plot(range(x.vals), c(0, max(quant)), type = "n", xlab = "", ylab = "")
    mtext(xlabs, 1, 2)
    mtext("hotspot size", 2, 2)
    abline(v = x.crit, col = "darkgray", lty = 2)
    abline(h = quant.thr, col = "darkgray", lty = 2)
    mtext(ceiling(quant.thr), 2, at = quant.thr, las = 2, cex = 0.5)
    for(i in seq(along = probs)) {
      lines(rev(sort(x.vals)), quant[i,],
            lwd = 1 + 2 * (round(probs[i] - level, 2) == 0))
    }
    text(x.crit, rev(quant[n.probs,])[thr.level] + 5, 1 - max(probs), adj = 0)
    text(x.crit, rev(quant[1,])[thr.level] - 5, 1 - min(probs), adj = 1)
    text(x.vals[n.probs], quant.thr + 5, 1 - level, adj = 1)
    if(add.level)
      main <- paste(main, "\n hotspot size significance level =", 1 - max(probs), "to", 1 - min(probs))
    mtext(main, 3, 0.5)

    quant.thr
  }

  tmpar <- par(mfrow = c(2,2), mar = c(3.1,3.1,2.1,0.1))
  ## Jansen method, smoothing.
  x$quant.thr <- tmp.plot(lod.thrs, x$quant.N.window, lod.thr, probs, level, FALSE,
                            "Jansen method 5cM window")
  
  tmp.plot(probs, x$quant.N.window, level, probs, level, TRUE,
           "Jansen method 5cM window")

  tmp.plot(lod.thrs, x$quant.N, lod.thr, probs, level, FALSE,
           "Jansen method per locus")
  tmp.plot(probs, x$quant.N, level, probs, level, TRUE,
           "Jansen method per locus")
  par(tmpar)
  
  ## Chaibub Neto method.
  plot(c(1,x$max.quant), c(min(x$quant, na.rm = TRUE), max(x$quant, na.rm = TRUE)), type = "n",
       xlab = "significant hotspot size with given threshold",
            ylab = "hotspot LOD score threshold",
            log = "xy")
  abline(h = lod.thr, col = "darkgray", lty = 2)
  abline(v = x$quant.thr, col = "darkgray", lty = 2)
  mtext(ceiling(x$quant.thr), 1, at = x$quant.thr, las = 2)
  for(i in seq(along = probs)) {
    tmp <- (x$quant[i,seq(x$max.quant)] >= lod.thr)
    if(any(tmp))
      lines(seq(x$max.quant)[tmp], x$quant[i,seq(x$max.quant)][tmp],
            lwd = 1 + 2 * (round(probs[i] - 0.95, 2) == 0))
    if(any(!tmp))
      lines(seq(x$max.quant)[!tmp], x$quant[i,seq(x$max.quant)][!tmp],
            lwd = 1 + 2 * (round(probs[i] - 0.95, 2) == 0),
            col = "darkgray")
  }
  
  n.thr2 <- length(lod.thrs) / 2
  text(n.thr2 + 1, x$quant[n.probs, n.thr2], 1 - max(probs), adj = 0)
  text(n.thr2 - 1, x$quant[1,n.thr2], 1 - min(probs), adj = 1)
  title(paste("hotspot LOD threshold by hotspot size\nsignificance level =",
              1 - max(probs), "to", 1 - min(probs)))
  
  invisible(x)
} 
