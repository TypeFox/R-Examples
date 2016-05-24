######################################################################
# diag.R
#
# Brian S Yandell
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: loci.qtlnet, est.qtlnet
# 
######################################################################

loci.qtlnet <- function(qtlnet.object, chr.pos=TRUE, merge.qtl = 10,
                        ...)
{
  ## Get QTL loci for average network.
  
  cross <- qtlnet.object$cross
  ## Make sure cross object has genotype probabilities.
  if (!("prob" %in% names(cross$geno[[1]]))) {
      warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }

  pheno.net.str <- get.averaged.net(qtlnet.object, ...)

  ## Extract needed attributes from qtlnet.object.
  pheno.nms <- attr(qtlnet.object, "pheno.names")
  addcov <- attr(qtlnet.object, "addcov") 
  intcov <- attr(qtlnet.object, "intcov") 
  thr <- attr(qtlnet.object, "threshold")
  method <- attr(qtlnet.object, "method")

  n.pheno <- length(pheno.nms)
  ss <- list()
  for(i in seq(n.pheno)){
    parents <- get.parents(pheno=pheno.nms[i], pheno.net.str=pheno.net.str)
    ss[[i]] <- get.ss(cross, i, parents, addcov[[i]], intcov[[i]], thr[[i]],
                      method)$ss
  }
  if(merge.qtl > 0) {
    ## If QTL within merge.qtl of mean pos for that chr, use mean pos.
    chr <- unlist(sapply(ss, function(x) x$chr))
    pos <- unlist(sapply(ss, function(x) x$pos))
    mpos <- tapply(pos, chr, mean)

    ## Need to re-get markers if !chr.pos.
    if(!chr.pos) {
      name.pos <- mpos
      for(j in names(mpos[!is.na(mpos)])) {
        map <- pull.loci(qtlnet.object$cross, j)
        wh <- names(which.min(abs(map - mpos[j]))[1])
        if(substring(wh, 1, 3) == "loc")
          wh <- paste("c", j, ".", wh, sep = "")
        name.pos[j] <- wh
      }
    }
    for(i in seq(n.pheno)) {
      ss.chr <- as.character(ss[[i]]$chr)
      dif <- mpos[ss.chr] - ss[[i]]$pos
      is.close <- abs(dif) <= merge.qtl
      if(any(is.close)) {
        ss[[i]]$pos <- ss[[i]]$pos + is.close * dif
        if(!chr.pos)
          row.names(ss[[i]])[is.close] <- name.pos[ss.chr][is.close]
      }
    }
  }
    
  QTLnodes <- list()
  for(i in seq(n.pheno)){
    markers <- row.names(ss[[i]])
    le.markers <- length(markers)
    if(le.markers > 0){ 
      if(chr.pos){
        QTLnodes[[i]] <- paste("chr", ss[[i]][,1], "@", round(ss[[i]][,2], 2), sep = "")
      }
      else{
        QTLnodes[[i]] <- markers
      }
    }
    if(le.markers == 0){
      QTLnodes[[i]] <- markers
    }   
  }
  names(QTLnodes) <- pheno.nms
  QTLnodes
}
######################################################################
est.qtlnet <- function(qtlnet.object, ..., verbose = TRUE)
{
  ## Get fit for average network.
  ## This does not incorporate QTL estimates same way as loci.qtlnet.
  
  cross <- qtlnet.object$cross
  ## Make sure cross object has genotype probabilities.
  if (!("prob" %in% names(cross$geno[[1]]))) {
      warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }

  pheno.net.str <- get.averaged.net(qtlnet.object, ...)

  ## Extract needed attributes from qtlnet.object.
  pheno.nms <- attr(qtlnet.object, "pheno.names")
  addcov <- attr(qtlnet.object, "addcov") 
  intcov <- attr(qtlnet.object, "intcov") 
  thr <- attr(qtlnet.object, "threshold")
  method <- attr(qtlnet.object, "method")

  n.pheno <- length(pheno.nms)
  est <- list()
  for(i in seq(n.pheno)){
    if(verbose)
      cat(pheno.nms[i], "\n")
    
    parents <- get.parents(pheno=pheno.nms[i], pheno.net.str=pheno.net.str)
    ss <- get.ss(cross, i, parents, addcov[[i]], intcov[[i]], thr[[i]],
                 method)
    crossi <- ss$cross
    covar <- ss$covar
    intcov.names <- dimnames(ss$intcovar)[[2]]
    cov.names <- dimnames(covar)[[2]]
    ss <- ss$ss
    
    if(nrow(ss)) {
      ## Determine additive and interactive covariates.
      qtl <- makeqtl(crossi, ss$chr, ss$pos, what = "prob")
      qs <- paste("Q", seq(length(ss$chr)), sep = "", collapse = "+")
      if(length(intcov.names))
        qs <- paste("(", paste(intcov.names, collapse = "+"), ")*",
                    "(", qs, ")", sep = "")
      
      form <- formula(paste("y ~", qs,
                            ifelse(is.null(cov.names), "",
                                   paste(  "+", paste(cov.names, collapse = "+")))))
      est[[i]] <- fitqtl(crossi, i, qtl, covar, form, "hk",
                         get.ests = TRUE)$ests$ests
    }
    else { ## No QTL!
      form <- formula(paste("y ~", paste(cov.names, collapse = "+")))
      data <- as.data.frame(covar)
      data$y <- crossi$pheno[[i]]
      est[[i]] <- coef(lm(form, data))
    }
  }
  names(est) <- pheno.nms
  est
}

######################################################################
get.ss <- function(cross, phenoi, parents, addcov, intcov, thr, method)
{
  ## Get SS for model. Includes chr and pos.

  ## Reduce to nonmissing data.
  cov.names <- unique(c(addcov, intcov))
  tmp <- unique(c(phenoi, parents, cov.names))
  pheno.na <- apply(cross$pheno[, tmp, drop = FALSE], 1,
                    function(x) any(is.na(x)))
  crossi <- subset(cross, ind = !pheno.na)

  ## Determine parent covariates.
  pacov.dat <- NULL
  if(!is.null(parents))
    pacov.dat <- as.matrix(crossi$pheno[, parents, drop = FALSE])

  ## Determine additive and interactive covariates.
  addcov.dat <- create.cov.matrix(crossi, cov.names = cov.names)
  addcov.dat <- cbind(pacov.dat, addcov.dat)
  if(!is.null(addcov.dat))
    dimnames(addcov.dat)[[2]] <- make.names(dimnames(addcov.dat)[[2]])
  intcov.dat <- create.cov.matrix(crossi, cov.names = intcov)
  if(!is.null(intcov.dat))
    dimnames(intcov.dat)[[2]] <- make.names(dimnames(intcov.dat)[[2]])
    
  ss <- scanone.summary(crossi, phenoi, addcov.dat, intcov.dat, thr, method)
  list(ss = ss, cross = crossi, covar = addcov.dat, intcovar = intcov.dat)
}
######################################################################
get.parents <- function(pheno, pheno.net.str)
{
  aux1 <- which(pheno.net.str[,2] == pheno)
  pa <- pheno.net.str[pheno.net.str[,2] == pheno, 1]
  if(length(pa) == 0) return(NULL)
  else return(pa)
}
######################################################################
pull.loci <- function(cross, i)
{
  ## Taken directly from R/qtl scanone.
  stp <- attr(cross$geno[[i]]$prob, "step")
  oe <- attr(cross$geno[[i]]$prob, "off.end")
  if ("stepwidth" %in% names(attributes(cross$geno[[i]]$prob))) 
    stpw <- attr(cross$geno[[i]]$prob, "stepwidth")
  else stpw <- "fixed"
  map <- create.map(cross$geno[[i]]$map, stp, oe, stpw)
}
